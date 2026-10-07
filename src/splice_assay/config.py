"""Settings: every threshold and drawing constant in one place.

The defaults are the conventions of the analysis these figures were designed for (docs/methods.md). Change them
with `Settings(min_pairs=5)`, `settings.replace(...)`, or a JSON file of overrides (`Settings.from_json`).
"Case" and "reference" are the two groups compared (by default tumour and normal; see Dataset).
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, fields, replace
from pathlib import Path

import numpy as np

TIME_UNITS = {"days": None, "months": 12.0, "years": 1.0}
RIDGE_SCOPES = ("none", "clinical", "molecular", "all")
DRAWING = {"km_max_years", "km_tick_years", "gene_model_min_frac", "nested_biotypes", "exclude_transcript_types",
           "gtf_flank", "case_label", "reference_label", "formats", "dpi", "cohorts_per_page", "q_mark_below",
           "imprecise_ci_ratio"}


def _split_rule(v, name: str):
    """'median', 'mean' or a number (a numeric string from the command line counts)."""
    if isinstance(v, str) and v.strip().lower() in ("median", "mean"):
        return v.strip().lower()
    try:
        x = float(v)
    except (TypeError, ValueError):
        raise ValueError(f'Settings.{name} must be "median", "mean" or a number, not {v!r}') from None
    if not np.isfinite(x):
        raise ValueError(f"Settings.{name} must be finite")
    return x


@dataclass(frozen=True)
class Settings:
    # ------------------------------------------------------------------ tumour vs normal
    min_pairs: int = 10               # the paired test needs at least this many tumour-normal pairs
    min_group: int = 10               # the unpaired test needs at least this many case and reference patients (one
                                      # value each: the mean of a patient's samples in the group)
    alpha: float = 0.05               # significance level of the hit rules
    min_abs_delta: float = 0.10       # a hit needs round(|delta median PSI|, 12) > this
    robust_01: bool = True            # a hit must survive dropping PSI values of exactly 0 or 1
    gex_min_abs_delta: float = 1.0    # an expression hit needs |delta median| > this (on log2 expression: two-fold)
    hit_min_abs_delta: float = 0.20   # a HIT-index hit needs |delta median| > this (the -1..1 scale is twice PSI's);
                                      # HIT events have no 0/1 robustness check
    fdr_min_family: int = 10          # q values (BH within a gene) only for families of at least this many tests
    # ------------------------------------------------------------------ survival
    coverage_frac: float = 0.5        # PSI observed in at least this fraction of the cohort's survival samples
    min_off_modal: int = 10           # ... and at least this many observed values away from the modal value
    km_split: object = "median"       # the KM split of an event's values (PSI or HIT index): "median", "mean" or a
                                      # number; the high arm is above it (at it, when the median is the highest
                                      # value: stats.survival.km_cut)
    km_split_expression: object = "median"  # the same for host-gene expression (a number is in expression units)
    km_min_group: int = 10            # log-rank: at least this many patients in each arm
    km_min_events: int = 10           # ... and at least this many events
    cox_min_n: int = 30               # Cox: at least this many patients
    cox_min_events: int = 10          # ... and at least this many events
    cox_low_power_events: int = 20    # a Cox fit with fewer events is marked low power (a note, ‡ in the forest);
                                      # the fit stands (0 = never)
    cox_min_events_per_term: float = 0.0  # opt-in: no Cox fit with fewer events per estimated term (e.g. 5; 0 = no
                                          # such gate)
    cox_ridge: str = "none"           # opt-in ridge (L2) penalty on Cox terms: "none", "clinical" (clinical
                                      # covariates), "molecular" (PSI or the HIT index, and host expression) or "all"
    cox_ridge_penalty: float = 1.0    # its strength: lambda/2 * beta^2 per penalized term, beta the log HR per SD (per
                                      # level for categories), like a normal prior with SD 1/sqrt(lambda) on each
    low_psi_variance_sd: float = 0.002  # no survival test when SD(PSI) of the cohort is below this
    nr_refit_step_size: float = 0.5   # a Newton-Raphson failure is refitted once with this step size
    psi_step: float = 0.10            # the Cox coefficient is per +0.10 PSI
    psi_hr_unit: str = "sd"           # the PSI (or HIT index) hazard ratio shown: per "sd" or per "iqr" of the fit
                                      # cohort's values (the tables hold both)
    covariate_min_complete: float = 0.8  # a clinical variable recorded for fewer of a cohort's patients is left out
                                         # of that cohort's model (else it would shrink the cohort)
    level_min_patients: int = 10      # a category of a clinical covariate with fewer patients (or no events) in a
                                      # fit cohort is merged into its neighbour (stage I into II); 0 = never
    cox_events_per_term: float = 10.0  # a Cox fit with fewer events per estimated term is flagged as an overfit
                                       # risk (a note only; 0 = never)
    narrow_psi_below: float = 0.05    # a Cox fit whose PSI spread in the fit cohort is below this is flagged as a
                                      # narrow PSI range: its HR covers a few PSI points (a note only; 0 = never)
    narrow_psi_measure: str | None = None  # the spread checked: "iqr" or "sd" (default: psi_hr_unit, the HR's unit)
    ci_z: float = 1.96                # normal quantile of every 95% interval the package computes
    ph_test: bool = True              # proportional-hazards (Schoenfeld) tests of every Cox term and of the KM split
    ph_note_below: float = 0.05       # a PH test p below this adds a note on the page and in the tables; the result
                                      # stands (0 = never)
    time_unit: str = "days"           # unit of survival.time: days, months or years
    days_per_year: float = 365.25     # days in a year, for converting survival times to years
    # ------------------------------------------------------------------ correlation (opt-in: probe --correlation)
    corr_min_n: int = 20              # Spearman's rho (an event with its host gene's expression, or two events of a
                                      # gene, per cohort) needs at least this many patients with both values
    corr_note_above: float = 0.7      # |rho| at or above this is called strong: a frame in correlation.png, named in
                                      # the probe report (0 = never)
    min_observed: float | None = None  # the probe's gene maps and correlations keep the events observed in at least
                                      # this fraction of a cohort's survival samples, in some cohort (0 = all;
                                      # default: coverage_frac, the share the tests need)
    # ------------------------------------------------------------------ numerics
    round_decimals: int = 12          # PSI and differences are rounded to this many decimals before comparisons
    # ------------------------------------------------------------------ drawing
    km_max_years: float = 10.0        # KM axes end at this many years or at the first tick past the
                                      # last time, whichever is earlier
    km_tick_years: float = 2.0        # a KM axis tick every this many years
    gene_model_min_frac: float = 0.10  # collapsed gene model: keep exons used by >= this share of transcripts
    nested_biotypes: tuple = ("snoRNA", "scaRNA")  # genes of these types inside the host gene get a track of
                                      # their own under its model
    exclude_transcript_types: tuple = ("retained_intron",)  # transcripts of these types are left out of the
                                      # collapsed gene model
    gtf_flank: int = 20000            # GTF records are read within this many nt of the events
    case_label: str | None = None     # display name of the case group (default: as spelled in the data)
    reference_label: str | None = None  # display name of the reference group
    cohorts_per_page: int = 6         # more cohorts on a page are split into balanced pages of at most this many,
                                      # each with the full forest (0 = one page)
    q_mark_below: float = 0.05        # a q below this is marked * on the pages, beside the p-based fill (0 = never)
    imprecise_ci_ratio: float = 8.0   # the probe's overview and gene maps draw a survival cell pale and hatched when
                                      # its 95% CI spans more than this ratio (high / low) or its fit is low power,
                                      # so a large but unreliable HR does not stand out (0 = never)
    formats: tuple = ("svg", "pdf", "png")  # the file types of each figure (any of the three; probe.pdf
                                      # is always written)
    dpi: int = 400                    # resolution of the PNG files

    def __post_init__(self):
        for name in ("min_pairs", "min_group", "min_off_modal", "km_min_group", "km_min_events", "cox_min_n",
                     "cox_min_events", "cox_low_power_events", "round_decimals", "dpi", "gtf_flank", "fdr_min_family",
                     "level_min_patients", "cohorts_per_page", "corr_min_n"):
            if int(getattr(self, name)) != getattr(self, name) or getattr(self, name) < 0:
                raise ValueError(f"Settings.{name} must be a non-negative integer")
        if self.corr_min_n < 3:
            raise ValueError("Settings.corr_min_n must be at least 3")
        for name in ("corr_note_above", "min_observed"):
            if getattr(self, name) is not None and not 0 <= getattr(self, name) <= 1:
                raise ValueError(f"Settings.{name} must lie in [0, 1]")
        for name in ("alpha", "coverage_frac", "gene_model_min_frac", "covariate_min_complete"):
            if not 0 < getattr(self, name) <= 1:
                raise ValueError(f"Settings.{name} must lie in (0, 1]")
        for name in ("min_abs_delta", "low_psi_variance_sd", "gex_min_abs_delta", "cox_events_per_term",
                     "narrow_psi_below", "ph_note_below", "hit_min_abs_delta", "q_mark_below"):
            if getattr(self, name) < 0:
                raise ValueError(f"Settings.{name} must be >= 0")
        if self.imprecise_ci_ratio != 0 and not self.imprecise_ci_ratio > 1:
            raise ValueError("Settings.imprecise_ci_ratio must be 0 (never) or above 1")
        if self.cox_ridge not in RIDGE_SCOPES:
            raise ValueError(f"Settings.cox_ridge must be one of {', '.join(RIDGE_SCOPES)}")
        if self.cox_min_events_per_term < 0:
            raise ValueError("Settings.cox_min_events_per_term must be >= 0")
        for name in ("psi_step", "ci_z", "days_per_year", "km_max_years", "km_tick_years", "nr_refit_step_size",
                     "cox_ridge_penalty"):
            if not getattr(self, name) > 0:
                raise ValueError(f"Settings.{name} must be > 0")
        for name in ("km_split", "km_split_expression"):
            object.__setattr__(self, name, _split_rule(getattr(self, name), name))
        if self.narrow_psi_measure not in (None, "iqr", "sd"):
            raise ValueError('Settings.narrow_psi_measure must be "iqr" or "sd"')
        if self.psi_hr_unit not in ("iqr", "sd"):
            raise ValueError('Settings.psi_hr_unit must be "iqr" or "sd"')
        if self.time_unit not in TIME_UNITS:
            raise ValueError(f"Settings.time_unit must be one of {sorted(TIME_UNITS)}")
        bad = set(self.formats) - {"svg", "pdf", "png"}
        if bad:
            raise ValueError(f"Settings.formats: unsupported {sorted(bad)}")
        for name in ("nested_biotypes", "exclude_transcript_types", "formats"):
            object.__setattr__(self, name, tuple(getattr(self, name)))

    # ------------------------------------------------------------------ helpers
    @property
    def narrow_measure(self) -> str:
        """The spread the narrow-range note checks: narrow_psi_measure, else the HR's unit."""
        return self.narrow_psi_measure or self.psi_hr_unit

    @property
    def observed_frac(self) -> float:
        """The share of a cohort's survival samples an event must be observed in for the gene maps and correlations:
        min_observed, else coverage_frac."""
        return self.coverage_frac if self.min_observed is None else self.min_observed

    def replace(self, **changes) -> "Settings":
        return replace(self, **changes)

    def to_dict(self) -> dict:
        """Plain Python values (JSON-ready: tuples as lists, numpy numbers as Python numbers)."""
        return {k: (list(v) if isinstance(v, tuple) else v.item() if isinstance(v, np.generic) else v)
                for k, v in asdict(self).items()}

    def changed(self) -> dict:
        """The analysis settings (gates and thresholds, not drawing) that differ from the defaults: name -> (value,
        default). Reports and figures list them so a relaxed run is never mistaken for a default one."""
        default = Settings()
        return {f.name: (getattr(self, f.name), getattr(default, f.name)) for f in fields(self)
                if f.name not in DRAWING and getattr(self, f.name) != getattr(default, f.name)}

    def changed_text(self) -> str:
        """'min_pairs 7 (default 10), km_split mean (default median)', or ''."""
        f = (lambda x: f"{x:g}" if isinstance(x, (int, float)) and not isinstance(x, bool) else str(x))
        return ", ".join(f"{k} {f(v)} (default {f(d)})" for k, (v, d) in self.changed().items())

    @classmethod
    def from_dict(cls, d: dict) -> "Settings":
        if not isinstance(d, dict):
            raise ValueError('expected an object of settings, e.g. {"min_pairs": 5}')
        known = {f.name for f in fields(cls)}
        unknown = sorted(set(d) - known)
        if unknown:
            raise ValueError(f"unknown setting(s): {', '.join(unknown)}")
        return cls(**d)

    @classmethod
    def from_json(cls, path) -> "Settings":
        return cls.from_dict(json.loads(Path(path).read_text()))

    def years(self, t) -> np.ndarray:
        """Survival times converted to years."""
        t = np.asarray(t, float)
        if self.time_unit == "days":
            return t / self.days_per_year
        return t / TIME_UNITS[self.time_unit]
