"""Run the statistics for every event x cohort (case vs reference) and event x cohort x endpoint (survival)."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from .config import Settings
from .dataset import Dataset, InputError
from .events import opt_in, quantity
from .stats.survival import CoxModel, expression_cell, survival_cell
from .stats.tissue import compare_groups


# ============================================================================================ designs
@dataclass
class GroupDesign:
    """Sample IDs per cohort: pairs (case and reference, aligned), all case and all reference samples."""

    cohorts: dict

    @classmethod
    def from_dataset(cls, ds: Dataset) -> "GroupDesign":
        out = {}
        for c in ds.cohorts:
            s = ds.samples[ds.samples.cohort.eq(c)]
            p = ds.pairs[ds.pairs.cohort.eq(c)]
            out[c] = dict(pair_c=p.case_sample.to_numpy(), pair_r=p.reference_sample.to_numpy(),
                          cases=s.sample_id[s.role.eq("case")].to_numpy(),
                          refs=s.sample_id[s.role.eq("reference")].to_numpy())
        return cls(out)


@dataclass
class SurvivalDesign:
    """Per cohort: the survival samples (and their clinical rows); per endpoint: positions of patients with a valid
    row, their time and event."""

    cohorts: dict

    @classmethod
    def from_dataset(cls, ds: Dataset, endpoints, clinical_columns=()) -> "SurvivalDesign":
        out = {}
        missing = [c for c in clinical_columns if ds.clinical is None or c not in ds.clinical.columns]
        if missing:
            raise InputError(f"clinical column(s) not found: {', '.join(missing)}"
                             + ("" if ds.clinical is not None else " (no clinical table was given)"))
        for c in ds.cohorts:
            base = ds.samples[ds.samples.cohort.eq(c) & ds.samples.survival_cohort]
            d = dict(samples=base.sample_id.to_numpy(), ep={}, clinical=None)
            if clinical_columns:
                d["clinical"] = ds.clinical.reindex(base.patient_id.to_numpy())[list(clinical_columns)] \
                    .reset_index(drop=True)
            for ep in endpoints:
                sv = ds.survival[ds.survival.endpoint.eq(ep)].set_index("patient_id")
                pos = np.flatnonzero(base.patient_id.isin(sv.index).to_numpy())
                pid = base.patient_id.to_numpy()[pos]
                d["ep"][ep] = dict(pos=pos, time=sv.time.reindex(pid).to_numpy(float),
                                   event=sv.event.reindex(pid).to_numpy(int))
            out[c] = d
        return cls(out)


# ============================================================================================ results
@dataclass
class Results:
    """groups: one row per event x cohort; survival: one row per event x cohort x endpoint; cox_terms: one row per
    Cox model term of every fitted cell."""

    groups: pd.DataFrame
    survival: pd.DataFrame
    cox_terms: pd.DataFrame = field(default_factory=pd.DataFrame)
    settings: Settings = field(default_factory=Settings)
    model: CoxModel = field(default_factory=CoxModel)

    def cells(self) -> pd.DataFrame:
        """Both tables side by side: one row per event x cohort x endpoint."""
        if self.survival.empty:
            return self.groups.copy()
        return self.survival.merge(self.groups, on=["event_id", "gene", "cohort"], how="left")

    def write(self, out_dir) -> dict[str, Path]:
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        paths = dict(groups=out / "group_tests.csv", survival=out / "survival.csv", cox_terms=out / "cox_terms.csv")
        self.groups.to_csv(paths["groups"], index=False, lineterminator="\n")
        self.survival.to_csv(paths["survival"], index=False, lineterminator="\n")
        self.cox_terms.to_csv(paths["cox_terms"], index=False, lineterminator="\n")
        return paths

    @classmethod
    def read(cls, folder) -> "Results":
        f = Path(folder)
        rd = lambda n: pd.read_csv(f / n, low_memory=False, float_precision="round_trip")  # noqa: E731
        terms = rd("cox_terms.csv") if (f / "cox_terms.csv").stat().st_size > 1 else pd.DataFrame()
        return cls(rd("group_tests.csv"), rd("survival.csv"), terms)


GROUP_FIRST = ["event_id", "gene", "cohort", "paired_status", "paired_n_pairs", "paired_delta_median", "paired_hl",
               "paired_p", "paired_q", "paired_state_01", "paired_hit_status", "paired_hit", "unpaired_status",
               "unpaired_n_case", "unpaired_n_reference", "unpaired_delta_median", "unpaired_hl", "unpaired_p",
               "unpaired_q", "unpaired_state_01",
               "unpaired_hit_status", "unpaired_hit", "group_hit", "within_patient_support", "no_within_patient_check",
               "composition_sensitive"]
SV_FIRST = ["event_id", "gene", "cohort", "endpoint", "n_survival", "n_obs", "frac_obs", "off_modal", "eligible",
            "cutoff", "km_split", "km_status", "km_n", "n_low", "n_high", "events_low", "events_high", "logrank_hr",
            "km_p", "km_q", "cox_model", "cox_status", "cox_n", "cox_events", "psi_iqr", "cox_beta", "cox_se", "cox_p",
            "cox_q", "hr_per_step",
            "ci_low_step", "ci_high_step", "hr_per_iqr", "ci_low_iqr", "ci_high_iqr", "cox_events_per_term",
            "psi_narrow", "ph_p", "km_ph_p", "km_notes"]


def event_settings(event_type, s: Settings) -> Settings:
    """The rules of an event type. HIT index events (-1..1, not PSI) have no 0/1 robustness check, because -1, 0 and
    1 carry meaning, and their own effect threshold (hit_min_abs_delta); every other type uses `s` as it is."""
    if str(event_type).upper() == "HIT":
        return s.replace(robust_01=False, min_abs_delta=s.hit_min_abs_delta)
    return s


def default_events(ds: Dataset, include_hit: bool = False) -> list[str]:
    """The events analysed when none are named: every event except HIT-index ones, unless include_hit."""
    return [e for e in ds.psi.index if include_hit or not opt_in(ds.events.at[e, "event_type"])]


def _with_family(df: pd.DataFrame, ds: Dataset) -> pd.DataFrame:
    """A temporary _family column, the quantity of each row's event: HIT-index tests form their own q families, so
    including them never changes the q values of the PSI events."""
    if not len(df):
        return df.assign(_family=pd.Series(dtype=object))
    return df.assign(_family=df.event_id.map(ds.events.event_type).map(quantity))


def gene_fdr(df: pd.DataFrame, test: str, by: list[str], min_family: int) -> pd.DataFrame:
    """Benjamini-Hochberg q of `<test>_p` within each family: the rows sharing `by` (e.g. a gene, or a gene and an
    endpoint) whose `<test>_status` is 'tested'. Adds `<test>_q` and `<test>_q_tests` (the family size); a family
    smaller than `min_family` gets no q."""
    from scipy.stats import false_discovery_control
    q = pd.Series(np.nan, index=df.index)
    n = pd.Series(np.nan, index=df.index)
    if f"{test}_p" in df.columns and len(df):
        ok = df[f"{test}_status"].eq("tested") & df[f"{test}_p"].notna()
        for _, idx in df[ok].groupby(by, sort=False).groups.items():
            n.loc[idx] = len(idx)
            if len(idx) >= min_family:
                q.loc[idx] = false_discovery_control(df.loc[idx, f"{test}_p"].to_numpy(float), method="bh")
    return df.assign(**{f"{test}_q": q, f"{test}_q_tests": n})


def _order(df: pd.DataFrame, first: list[str]) -> pd.DataFrame:
    """The standard columns first, every one present (empty when nothing was tested), then the others."""
    df = df.assign(**{c: np.nan for c in first if c not in df.columns})
    return df[first + [c for c in df.columns if c not in first]]


def _names(x, what, known) -> list[str]:
    x = [str(v) for v in ([x] if isinstance(x, str) else x)]
    bad = [v for v in x if v not in known]
    if bad:
        raise InputError(f"{what} not in the data: {', '.join(bad[:5])}")
    return x


# ============================================================================================ run
def analyse(ds: Dataset, events=None, endpoints=None, cohorts=None, settings: Settings | None = None,
            model: CoxModel | None = None, include_hit: bool = False) -> Results:
    """Case-vs-reference and survival statistics for the chosen events, cohorts and endpoints (default: all; without
    named events, HIT-index events only with include_hit).

    q values: Benjamini-Hochberg within each gene, separately for the paired tests, the all-samples tests, the KM
    tests and the Cox PSI terms (the survival ones per endpoint), over the events and cohorts analysed together;
    HIT-index events form families of their own. For gene-wide q, analyse all of a gene's events (panel and probe
    do)."""
    s = settings or Settings()
    model = model or CoxModel()
    if events is None:
        events = default_events(ds, include_hit)
        if not events and len(ds.psi):
            raise InputError("only HIT-index events in the data: add --include-hit (include_hit=True), or name them")
    else:
        events = _names(events, "event(s)", set(ds.psi.index))
    cohorts = ds.cohorts if cohorts is None else _names(cohorts, "cohort(s)", set(ds.cohorts))
    if endpoints is None:
        endpoints = ds.endpoints
    else:
        if ds.survival is None:
            raise InputError("no survival table was given")
        endpoints = _names(endpoints, "endpoint(s)", set(ds.endpoints))
    gd = GroupDesign.from_dataset(ds)
    sd = SurvivalDesign.from_dataset(ds, endpoints, model.clinical_columns) if endpoints else None
    g_rows, sv_rows, t_rows = [], [], []
    for eid in events:
        gene = ds.events.at[eid, "gene"]
        etype = ds.events.at[eid, "event_type"]
        se = event_settings(etype, s)                   # HIT index events have their own rules
        v = ds.psi.loc[eid]
        host_gene = ds.events.at[eid, "expression_gene"]
        host = None
        if ds.expression is not None and host_gene in ds.expression.index:
            host = ds.expression.loc[host_gene]
        for c in cohorts:
            d = gd.cohorts[c]
            g_rows.append(dict(event_id=eid, gene=gene, cohort=c,
                               **compare_groups(v, d["pair_c"], d["pair_r"], d["cases"], d["refs"], se)))
            if sd is None:
                continue
            base = sd.cohorts[c]
            x_base = v.reindex(base["samples"]).to_numpy(float)
            h_base = None if host is None else host.reindex(base["samples"]).to_numpy(float)
            for ep in endpoints:
                E = base["ep"][ep]
                rec = dict(event_id=eid, gene=gene, cohort=c, endpoint=ep)
                row, terms = survival_cell(x_base, E["pos"], E["time"], E["event"], h_base, se, model,
                                           base["clinical"], quantity=quantity(etype))
                if ds.expression is not None and host is None and model.expression and \
                        row["cox_status"] != "coverage_gate":
                    row = {k: v_ for k, v_ in row.items() if not k.startswith(("cox_", "hr_", "ci_", "ph_", "expr_"))
                           and k not in ("psi_iqr", "psi_sd")}
                    row.update(cox_model=model.describe(True), cox_status="no_host_expression")
                    terms = []
                sv_rows.append(dict(rec, **row))
                t_rows += [dict(rec, **t) for t in terms]
    groups = _with_family(pd.DataFrame(g_rows), ds)
    for test in ("paired", "unpaired"):                    # q within each gene, per test (see docs/methods.md)
        groups = gene_fdr(groups, test, ["gene", "_family"], s.fdr_min_family)
    groups = _order(groups.drop(columns="_family"), GROUP_FIRST)
    sv = _with_family(pd.DataFrame(sv_rows), ds)
    for test in ("km", "cox"):
        sv = gene_fdr(sv, test, ["gene", "endpoint", "_family"], s.fdr_min_family)
    sv = _order(sv.drop(columns="_family"), SV_FIRST) if sv_rows else pd.DataFrame(columns=SV_FIRST)
    if len(sv):
        with np.errstate(invalid="ignore"):
            km_sig = sv.km_status.eq("tested") & (sv.get("km_p", np.nan) < s.alpha)
            cox_sig = sv.cox_status.eq("tested") & (sv.get("cox_p", np.nan) < s.alpha)
        sv["survival_hit"] = (km_sig | cox_sig).to_numpy()
    return Results(groups, sv, pd.DataFrame(t_rows), s, model)


# ============================================================================================ host-gene expression
@dataclass
class ExpressionResults:
    """The same statistics for host-gene expression itself. groups: one row per gene x cohort; survival: one row per
    gene x cohort x endpoint (KM split by Settings.km_split_expression, Cox on expression + the model's clinical
    covariates); cox_terms: every term of every fitted model."""

    groups: pd.DataFrame
    survival: pd.DataFrame
    cox_terms: pd.DataFrame = field(default_factory=pd.DataFrame)
    settings: Settings = field(default_factory=Settings)
    model: CoxModel = field(default_factory=CoxModel)

    def cells(self) -> pd.DataFrame:
        if self.survival.empty:
            return self.groups.copy()
        return self.survival.merge(self.groups, on=["gene", "cohort"], how="left")

    def write(self, out_dir) -> dict[str, Path]:
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        paths = dict(groups=out / "expression_group_tests.csv", survival=out / "expression_survival.csv",
                     cox_terms=out / "expression_cox_terms.csv")
        self.groups.to_csv(paths["groups"], index=False, lineterminator="\n")
        self.survival.to_csv(paths["survival"], index=False, lineterminator="\n")
        self.cox_terms.to_csv(paths["cox_terms"], index=False, lineterminator="\n")
        return paths


GEX_SV_FIRST = ["gene", "cohort", "endpoint", "n_survival", "n_obs", "frac_obs", "off_modal", "eligible", "cutoff",
                "km_split", "km_status", "km_n", "n_low", "n_high", "events_low", "events_high", "logrank_hr", "km_p",
                "cox_model", "cox_status", "cox_n", "cox_events", "expr_sd", "cox_beta", "cox_se", "cox_p", "hr_per_sd",
                "ci_low_sd", "ci_high_sd", "cox_events_per_term", "ph_p", "km_ph_p", "km_notes"]


def analyse_expression(ds: Dataset, genes=None, endpoints=None, cohorts=None, settings: Settings | None = None,
                       model: CoxModel | None = None) -> ExpressionResults:
    """Case-vs-reference, KM and Cox statistics of host-gene expression (default: the expression gene of every
    event). The group tests are those of PSI without the 0/1 check; a hit needs |delta| > gex_min_abs_delta. Cox:
    expression per SD plus the model's clinical covariates and strata."""
    if ds.expression is None:
        raise InputError("no expression table was given")
    s = settings or Settings()
    sg = s.replace(min_abs_delta=s.gex_min_abs_delta, robust_01=False)
    model = model or CoxModel()
    want = list(dict.fromkeys(ds.events.expression_gene)) if genes is None else [str(g) for g in
                                                                                 ([genes] if isinstance(genes, str)
                                                                                  else genes)]
    genes = [g for g in want if g in ds.expression.index]
    if not genes:
        raise InputError(f"no expression for gene(s) {', '.join(want[:5])}")
    cohorts = ds.cohorts if cohorts is None else _names(cohorts, "cohort(s)", set(ds.cohorts))
    endpoints = (ds.endpoints if endpoints is None else _names(endpoints, "endpoint(s)", set(ds.endpoints))) \
        if ds.survival is not None else []
    gd = GroupDesign.from_dataset(ds)
    sd = SurvivalDesign.from_dataset(ds, endpoints, model.clinical_columns) if endpoints else None
    g_rows, sv_rows, t_rows = [], [], []
    for g in genes:
        v = ds.expression.loc[g]
        for c in cohorts:
            d = gd.cohorts[c]
            g_rows.append(dict(gene=g, cohort=c, **compare_groups(v, d["pair_c"], d["pair_r"], d["cases"], d["refs"],
                                                                  sg)))
            if sd is None:
                continue
            base = sd.cohorts[c]
            x_base = v.reindex(base["samples"]).to_numpy(float)
            for ep in endpoints:
                E = base["ep"][ep]
                row, terms = expression_cell(x_base, E["pos"], E["time"], E["event"], s, model, base["clinical"])
                rec = dict(gene=g, cohort=c, endpoint=ep)
                sv_rows.append(dict(rec, **row))
                t_rows += [dict(rec, **t) for t in terms]
    groups = _order(pd.DataFrame(g_rows), ["gene", "cohort"] + GROUP_FIRST[3:])
    sv = _order(pd.DataFrame(sv_rows), GEX_SV_FIRST) if sv_rows else pd.DataFrame(columns=GEX_SV_FIRST)
    return ExpressionResults(groups, sv, pd.DataFrame(t_rows), s, model)

