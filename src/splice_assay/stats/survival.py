"""Survival per cohort and endpoint: a median-split Kaplan-Meier / log-rank test and a Cox model on PSI.

Survival cohort   the cohort's survival samples (one case sample per patient). Its median PSI is the KM cut,
                  shared by every endpoint; the coverage gate is applied to it.
Coverage gate     PSI observed in >= coverage_frac of the survival samples and >= min_off_modal observed values away
                  from the modal value; otherwise neither test runs ('coverage_gate').
Endpoint cohort   survival samples whose patient has a valid row for the endpoint and an observed PSI.
KM                The cut is the median of the survival samples' values (Settings.km_split: median, mean or a
                  number; km_split_expression for expression). PSI <= cut is the low arm. Needs km_min_group
                  patients per arm and km_min_events events. The log-rank HR is (O/E high) / (O/E low).
Cox               h(t) = h0(t) exp(b * PSI/0.10 + g * z(host expression) + covariates), Efron ties, no penalty
                  (lifelines). Host expression enters when an expression table is given (CoxModel.expression);
                  clinical covariates when named (CoxModel). Patients missing any model variable are left out.
                  Needs cox_min_n patients, cox_min_events events and min_off_modal PSI values away from the mode.
                  HR per IQR = exp(b * IQR/0.10), IQR of PSI over the fit cohort (NaN when the IQR is 0); HR per
                  SD the same with the SD. Settings.psi_hr_unit picks the one the pages show.
Low variance      SD(PSI) < low_psi_variance_sd stops the tests ('low_psi_variance').
Rare levels       a category with fewer than level_min_patients patients (or no events) in the fit cohort merges
                  into its neighbour before the fit (stage I into 'I–II'); cox_notes lists each merge.
Flags             notes only, the fit runs: fewer events per estimated term than cox_events_per_term (overfit
                  risk), a PSI spread (IQR or SD) in the fit cohort below narrow_psi_below (narrow PSI range), and
                  a proportional-hazards test below ph_note_below for any model term (cox_notes) or for the KM
                  split (km_notes). The PH tests are Schoenfeld tests on the Kaplan-Meier time scale (lifelines).
Failed fits       any lifelines convergence warning, an exception, a non-finite estimate or se <= 0 marks the fit
                  'failed' (no estimate is reported). A Newton-Raphson failure is refitted once with a smaller step.
Expression        `expression_cell` runs the same KM and Cox on host-gene expression itself (Cox: z(expression) +
                  the clinical covariates; HR per SD), for comparing the splicing and the expression association.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np
import pandas as pd

from ..config import Settings
from .core import coverage, logrank

NAN = float("nan")
NR_FAILURE_TEXT = "Newton-Raphson failed to converge"


@dataclass(frozen=True)
class CoxModel:
    """What the Cox model adjusts PSI for.

    expression   host-gene expression, z-scored in the fit cohort (when an expression table is given)
    covariates   clinical columns: numeric ones per SD (z-scored in the fit cohort; per unit with scale=False),
                 others as categories against a baseline level
    categorical  covariates to treat as categories even when numeric (e.g. a numeric stage code)
    strata       clinical columns used as strata (a separate baseline hazard per level)
    baseline     the reference level of categorical covariates, e.g. {"stage": "I"}; default the most common level
                 of the cohort (also when the named level is absent there)
    """
    expression: bool = True
    covariates: tuple = ()
    categorical: tuple = ()
    strata: tuple = ()
    scale: bool = True
    baseline: tuple = ()

    def __post_init__(self):
        for name in ("covariates", "categorical", "strata"):
            v = getattr(self, name)
            object.__setattr__(self, name, (v,) if isinstance(v, str) else tuple(v))
        b = self.baseline.items() if isinstance(self.baseline, dict) else self.baseline
        object.__setattr__(self, "baseline", tuple(sorted((str(k), str(v)) for k, v in b)))
        both = set(self.covariates) & set(self.strata)
        if both:
            raise ValueError(f"a column cannot be both a covariate and a stratum: {', '.join(sorted(both))}")

    @property
    def clinical_columns(self) -> list[str]:
        return list(dict.fromkeys(list(self.covariates) + list(self.strata)))

    def describe(self, with_expression: bool) -> str:
        terms = ["PSI"] + (["host expression"] if with_expression else []) + list(self.covariates)
        return " + ".join(terms) + (f"; strata: {', '.join(self.strata)}" if self.strata else "")

    def to_dict(self) -> dict:
        return dict(expression=self.expression, covariates=list(self.covariates), categorical=list(self.categorical),
                    strata=list(self.strata), scale=self.scale, baseline=dict(self.baseline))

    def with_clinical(self, covariates=(), categorical=(), strata=(), baseline=None) -> "CoxModel":
        """This model plus more clinical terms."""
        return CoxModel(self.expression, tuple(dict.fromkeys(self.covariates + tuple(covariates))),
                        tuple(dict.fromkeys(self.categorical + tuple(categorical))),
                        tuple(dict.fromkeys(self.strata + tuple(strata))), self.scale,
                        {**dict(self.baseline), **dict(baseline or {})})


# ============================================================================================ Cox fitting
def fit_cox(df: pd.DataFrame, term: str, extra: dict[str, str] | None = None, ph: bool = True,
            refit_step_size: float = 0.5, strata: list[str] | None = None) -> dict:
    """lifelines CoxPHFitter (Efron ties, no penalty). `term` is the tested coefficient; `extra` maps other
    coefficients to output prefixes; `summary` holds every coefficient. A Newton-Raphson failure is refitted once
    with `refit_step_size`."""
    out, nr = _fit_once(df, strata, term, extra, ph, None)
    if out["status"] == "failed" and nr:
        out2, _ = _fit_once(df, strata, term, extra, ph, refit_step_size)
        if out2["status"] == "tested":
            return dict(out2, fit_note=f"refit_step_size_{refit_step_size}")
        return dict(out2, fit_warning=f"{out.get('fit_warning', '')} || refit step_size {refit_step_size}: "
                                      f"{out2.get('fit_warning', '')}"[:400],
                    fit_note=f"refit_step_size_{refit_step_size}_failed")
    return dict(out, fit_note="")


def _fit_failure_warning(w) -> bool:
    """A lifelines ConvergenceWarning, or a numeric RuntimeWarning raised inside lifelines, fails a fit."""
    return "Convergence" in w.category.__name__ or (issubclass(w.category, RuntimeWarning)
                                                     and "lifelines" in str(w.filename))


def _fit_once(df, strata, term, extra, ph, step_size) -> tuple[dict, bool]:
    from lifelines import CoxPHFitter
    from lifelines.statistics import proportional_hazard_test

    cph = CoxPHFitter(penalizer=0.0)
    opts = None if step_size is None else {"step_size": step_size}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            cph.fit(df, duration_col="time", event_col="event", strata=strata, fit_options=opts)
        except Exception as e:  # noqa: BLE001
            return dict(status="failed", fit_warning=f"{type(e).__name__}: {e}"[:200]), \
                type(e).__name__ == "ConvergenceError"
    msgs = [str(w.message)[:120] for w in caught if _fit_failure_warning(w)]
    nr = any(NR_FAILURE_TEXT in str(w.message) for w in caught)
    s = cph.summary
    out = dict(beta=s.loc[term, "coef"], se=s.loc[term, "se(coef)"], p=s.loc[term, "p"], fit_warning=" | ".join(msgs))
    for col, pre in (extra or {}).items():
        out[f"{pre}_hr"], out[f"{pre}_p"] = s.loc[col, "exp(coef)"], s.loc[col, "p"]
    if msgs or not np.all(np.isfinite([out["beta"], out["se"]])) or not out["se"] > 0:
        return _failed(out), nr
    out["summary"] = s[["coef", "se(coef)", "p"]].copy()
    if ph:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                res = proportional_hazard_test(cph, df, time_transform="km")
                out["ph_p"] = float(res.summary.loc[term, "p"])
                out["ph_terms"] = {str(k): float(v) for k, v in res.summary["p"].items()}
        except Exception:  # noqa: BLE001
            out["ph_p"], out["ph_terms"] = NAN, {}
    return dict(out, status="tested"), False


def _failed(out: dict) -> dict:
    """A failed fit keeps only its warning: no estimate of any coefficient is reported."""
    blank = {k: (v if k == "fit_warning" else NAN) for k, v in out.items()}
    return dict(blank, status="failed")


# ============================================================================================ gates
def low_psi_variance(x, s: Settings) -> bool:
    x = np.asarray(x, float)
    return bool(len(x) >= 2 and np.std(x, ddof=1) < s.low_psi_variance_sd)


def km_cut(obs: np.ndarray, rule) -> tuple[float, str]:
    """The KM split of the observed values and how it was chosen: 'median', 'mean', or 'set' (a number given in
    Settings.km_split). The high arm is above the split."""
    if not len(obs):
        return NAN, ""
    if rule == "median":
        return float(np.median(obs)), "median"
    if rule == "mean":
        return float(np.mean(obs)), "mean"
    return float(rule), "set"


def km_gate(high: np.ndarray, e: np.ndarray, s: Settings) -> bool:
    return min(high.sum(), (~high).sum()) >= s.km_min_group and e.sum() >= s.km_min_events


def _is_numeric(v: pd.Series) -> bool:
    if pd.api.types.is_bool_dtype(v):
        return False
    if pd.api.types.is_numeric_dtype(v):
        return True
    return bool(pd.to_numeric(v, errors="coerce").notna().all())


def _complete_columns(C: pd.DataFrame, need, s: Settings, notes: list) -> list[str]:
    """The clinical columns recorded for enough of the fit cohort; the others are left out, with a note."""
    used = list(need)
    for c in need:
        frac = float(C[c].notna().mean()) if len(C) else 0.0
        if frac < s.covariate_min_complete:
            used.remove(c)
            notes.append(f"{c} left out ({frac:.0%} recorded)")
    return used


STAGE_ORDER = {"0": 0, "I": 1, "II": 2, "III": 3, "IV": 4}


def _level_order(levels) -> dict | None:
    """Sort keys when every level is a stage numeral (0, I-IV) or a number; None for unordered categories."""
    keys = {}
    for lev in levels:
        if lev in STAGE_ORDER:
            keys[lev] = STAGE_ORDER[lev]
        else:
            try:
                keys[lev] = float(lev)
            except ValueError:
                return None
    return keys


def merge_rare_levels(v: pd.Series, event: np.ndarray, min_n: int) -> tuple[pd.Series, dict, list[str]]:
    """Merge the levels of a categorical covariate that are too rare to estimate in the fit cohort: fewer than
    `min_n` patients, or no events. Ordered levels (stage numerals, numbers) merge into the adjacent level with fewer
    patients, the rarest first, so 'I' joins 'II' as 'I–II'; unordered levels merge into the most common level
    ('white+asian'). Returns the merged values, the label of each original level, and one note per merge.
    `min_n` = 0 merges nothing."""
    levels = sorted(v.unique())
    label = {lev: lev for lev in levels}
    if min_n <= 0 or len(levels) < 2:
        return v, label, []
    ev = np.asarray(event, int)
    n = {lev: int((v == lev).sum()) for lev in levels}
    d = {lev: int(ev[(v == lev).to_numpy()].sum()) for lev in levels}
    size = (lambda g: sum(n[x] for x in g))
    rare = (lambda g: size(g) < min_n or sum(d[x] for x in g) == 0)
    name = (lambda g: "–".join(dict.fromkeys([g[0], g[-1]])))
    notes = []
    order = _level_order(levels)
    if order is not None:
        groups = [[lev] for lev in sorted(levels, key=order.get)]
        while len(groups) > 1 and any(rare(g) for g in groups):
            i = min((j for j, g in enumerate(groups) if rare(g)), key=lambda j: (size(groups[j]), j))
            j = min((k for k in (i - 1, i + 1) if 0 <= k < len(groups)), key=lambda k: (size(groups[k]), k))
            notes.append(f"{name(groups[i])} merged into {name(groups[j])} ({size(groups[i])} patients, "
                         f"{sum(d[x] for x in groups[i])} events)")
            lo, hi = sorted((i, j))
            groups[lo:hi + 1] = [groups[lo] + groups[hi]]
        for g in groups:
            for lev in g:
                label[lev] = name(g)
    else:
        top = sorted(levels, key=lambda k: (-n[k], k))[0]
        small = [lev for lev in levels if lev != top and rare([lev])]
        if small:
            for lev in [top] + small:
                label[lev] = "+".join([top] + small)
            notes.append(f"{', '.join(small)} merged into {top} ({size(small)} patients, "
                         f"{sum(d[x] for x in small)} events)")
    return v.map(label), label, notes


def _clinical_design(df: pd.DataFrame, meta: list, notes: list, Cf: pd.DataFrame, model: CoxModel, used,
                     s: Settings) -> list:
    """Add the model's clinical covariates (numeric per SD, or categories against a baseline, rare levels merged) and
    strata to a design frame; return the strata columns."""
    for i, c in enumerate(model.covariates):
        if c not in used:
            continue
        v = Cf[c]
        if c not in model.categorical and _is_numeric(v):
            v = pd.to_numeric(v).astype(float).to_numpy()
            sd = float(np.std(v, ddof=1)) if len(v) > 1 else 0.0
            if not sd > 0:
                notes.append(f"{c} constant")
                continue
            df[f"cov{i}"] = (v - v.mean()) / sd if model.scale else v
            meta.append(dict(column=f"cov{i}", term=c, kind="numeric", level="", reference="",
                             unit="SD" if model.scale else "unit"))
        else:
            v = v.astype(str).str.strip()
            v, label, merges = merge_rare_levels(v, df["event"].to_numpy(), s.level_min_patients)
            notes.extend(f"{c} {m}" for m in merges)
            counts = v.value_counts()
            ref = sorted(counts.index, key=lambda k: (-counts[k], k))[0]
            want = dict(model.baseline).get(c)
            if want is not None:
                if label.get(want, want) in counts.index:
                    ref = label.get(want, want)
                else:
                    notes.append(f"{c} baseline {want} absent; {ref} used")
            levels = [k for k in sorted(counts.index) if k != ref]
            if not levels:
                notes.append(f"{c} has one level")
            for j, lev in enumerate(levels):
                df[f"cov{i}_{j}"] = (v == lev).astype(float).to_numpy()
                meta.append(dict(column=f"cov{i}_{j}", term=c, kind="categorical", level=lev, reference=ref,
                                 unit=f"vs {ref}"))
    strata_cols = []
    for k, c in enumerate([c for c in model.strata if c in used]):
        df[f"stratum{k}"] = Cf[c].astype(str).str.strip().to_numpy()
        strata_cols.append(f"stratum{k}")
    return strata_cols


def km_ph_test(t: np.ndarray, e: np.ndarray, high: np.ndarray) -> float:
    """Proportional-hazards test of the KM split: the Schoenfeld test (Kaplan-Meier time scale) of a Cox model on the
    high-arm indicator. A small p means the hazard ratio between the arms changes over follow-up (e.g. curves that
    cross), which the log-rank test averages over. NaN when it cannot be computed."""
    from lifelines import CoxPHFitter
    from lifelines.statistics import proportional_hazard_test
    df = pd.DataFrame({"time": np.asarray(t, float), "event": np.asarray(e, int), "high": np.asarray(high, float)})
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cph = CoxPHFitter(penalizer=0.0).fit(df, duration_col="time", event_col="event")
            return float(proportional_hazard_test(cph, df, time_transform="km").summary.loc["high", "p"])
    except Exception:  # noqa: BLE001
        return NAN


def _ph_note(fit: dict, meta: list, s: Settings) -> list[str]:
    """A note naming the model terms whose proportional-hazards test has p below ph_note_below (the fit stands)."""
    if not s.ph_note_below > 0:
        return []
    ph = fit.get("ph_terms") or {}
    bad = [(m, ph[m["column"]]) for m in meta if m["column"] in ph and ph[m["column"]] < s.ph_note_below]
    if not bad:
        return []
    from ..plot.style import fp
    name = (lambda m: m["term"] + (f" {m['level']}" if m["level"] else ""))
    return [f"non-proportional hazards: {', '.join(f'{name(m)} (p {fp(p)})' for m, p in bad)}"]


def _km_ph(r: dict, t, e, high, s: Settings) -> None:
    """The KM split's proportional-hazards p (km_ph_p) and its note (km_notes) on a tested KM row."""
    r["km_ph_p"] = km_ph_test(t, e, high) if s.ph_test else NAN
    if s.ph_note_below > 0 and r["km_ph_p"] < s.ph_note_below:
        from ..plot.style import fp
        r["km_notes"] = f"non-proportional hazards between the arms (p {fp(r['km_ph_p'])})"


def _fit_flags(r: dict, n_terms: int, s: Settings, spread: float | None = None, quantity: str = "PSI") -> list[str]:
    """Notes on a fitted model, which runs anyway: fewer events per estimated term than cox_events_per_term (an
    overfit risk), and, for PSI, a spread in the fit cohort below narrow_psi_below (the HR covers a few PSI points).
    Sets cox_events_per_term and psi_narrow on the row."""
    out = []
    r["cox_events_per_term"] = epv = r["cox_events"] / n_terms if n_terms else NAN
    if s.cox_events_per_term > 0 and epv < s.cox_events_per_term:
        out.append(f"{epv:.1f} events per term (< {s.cox_events_per_term:g})")
    if spread is not None:
        r["psi_narrow"] = bool(s.narrow_psi_below > 0 and spread < s.narrow_psi_below)
        if r["psi_narrow"]:
            out.append(f"narrow {quantity} range ({s.narrow_measure.upper()} {spread:.3f} < {s.narrow_psi_below:g})")
    return out


def _terms(fit: dict, meta: list, s: Settings) -> tuple[list[dict], list[str]]:
    """Every term of a fitted model (HR and CI per its unit), and the unstable ones (se > 3)."""
    terms, unstable = [], []
    summ = fit["summary"]
    ph = fit.get("ph_terms") or {}
    for m in meta:
        c, sei, p = (float(summ.loc[m["column"], k_]) for k_ in ("coef", "se(coef)", "p"))
        with np.errstate(over="ignore"):
            terms.append(dict({k_: v for k_, v in m.items() if k_ != "column"}, coef=c, se=sei, hr=float(np.exp(c)),
                              ci_low=float(np.exp(c - s.ci_z * sei)), ci_high=float(np.exp(c + s.ci_z * sei)), p=p,
                              ph_p=ph.get(m["column"], NAN)))
        if m["kind"] not in ("psi", "gex") and sei > 3:     # a sparse level or a nearly separated covariate
            unstable.append(f"{m['term']}{' ' + m['level'] if m['level'] else ''}")
    return terms, unstable


# ============================================================================================ one cohort
def survival_cell(x_base: np.ndarray, ep_pos: np.ndarray, time: np.ndarray, event: np.ndarray,
                  host_base: np.ndarray | None, s: Settings, model: CoxModel | None = None,
                  clinical: pd.DataFrame | None = None, quantity: str = "PSI") -> tuple[dict, list[dict]]:
    """KM and Cox statistics of one event in one cohort for one endpoint, and the Cox terms.

    x_base      PSI over the cohort's survival samples (NaN = missing)
    ep_pos      positions (into x_base) of the samples whose patient has a valid row for the endpoint
    time, event aligned with ep_pos
    host_base   host-gene expression over the survival samples, or None
    clinical    the model's clinical columns over the survival samples (rows aligned with x_base), or None
    """
    model = model or CoxModel()
    x_base = np.asarray(x_base, float)
    use_expr = host_base is not None and model.expression
    cv = coverage(x_base, s.round_decimals)
    obs = x_base[np.isfinite(x_base)]
    cut, how = km_cut(obs, s.km_split)
    r = dict(n_survival=len(x_base), n_obs=cv["n_obs"], frac_obs=cv["frac_obs"], off_modal=cv["off_modal"],
             modal_share=cv["modal_share"], cutoff=cut, km_split=how,
             eligible=bool(cv["frac_obs"] >= s.coverage_frac and cv["off_modal"] >= s.min_off_modal),
             n_endpoint=len(ep_pos), cox_model=model.describe(use_expr).replace("PSI", quantity, 1))
    if not r["eligible"]:
        return dict(r, km_status="coverage_gate", cox_status="coverage_gate"), []
    keep = np.isfinite(x_base[ep_pos])
    idx = ep_pos[keep]
    t, e, x = np.asarray(time, float)[keep], np.asarray(event, int)[keep], x_base[idx]
    high = x > r["cutoff"]
    low_var = low_psi_variance(x, s)
    # ---------------------------------------------------------------- KM
    r.update(km_n=len(x), n_low=int((~high).sum()), n_high=int(high.sum()), events_low=int(e[~high].sum()),
             events_high=int(e[high].sum()))
    if low_var:
        r["km_status"] = "low_psi_variance"
    elif not km_gate(high, e, s):
        r["km_status"] = "too_few_patients_or_events"
    else:
        lr = logrank(t, e, high)
        r.update(km_status="tested", o_high=lr["o_high"], e_high=lr["e_high"], o_low=lr["o_low"], e_low=lr["e_low"],
                 km_var=lr["var"], km_chi2=lr["chi2"], km_p=lr["p"], logrank_hr=lr["logrank_hr"])
        _km_ph(r, t, e, high, s)
    # ---------------------------------------------------------------- Cox: the fit cohort
    need = model.clinical_columns
    if need and clinical is None:
        return dict(r, cox_status="no_clinical_table"), []
    ok = np.ones(len(idx), bool)
    g = None
    notes = []
    if use_expr:
        g = np.asarray(host_base, float)[idx]
        ok &= np.isfinite(g)
    C, used = None, list(need)
    if need:
        C = clinical.iloc[idx].reset_index(drop=True)
        used = _complete_columns(C, need, s, notes)
        for c in used:
            ok &= C[c].notna().to_numpy()
    if len(used) < len(need):
        m_used = CoxModel(model.expression, tuple(c for c in model.covariates if c in used), model.categorical,
                          tuple(c for c in model.strata if c in used), model.scale, dict(model.baseline))
        r["cox_model"] = m_used.describe(use_expr).replace("PSI", quantity, 1)
    r["cox_n_dropped"] = int((~ok).sum())
    tt, ee, xx = t[ok], e[ok], x[ok]
    # ---------------------------------------------------------------- Cox: the design
    df, why = pd.DataFrame({"time": tt, "event": ee, "psi10": xx / s.psi_step}), "too_few_patients_or_events"
    unit = "PSI" if quantity == "PSI" else "HIT"           # the value's short name in units ("per IQR (0.04 HIT)")
    meta = [dict(column="psi10", term=quantity, kind="psi", level="", reference="", unit=f"+{s.psi_step:g} {unit}")]
    if use_expr:
        g = g[ok]
        if len(g) >= 3 and np.std(g) > 0:
            df["gex_z"] = (g - g.mean()) / g.std(ddof=1)
            meta.append(dict(column="gex_z", term="host expression", kind="expression", level="", reference="",
                             unit="SD"))
        else:
            df, why = None, ("constant_expression" if len(g) >= 3 else why)
    strata_cols = []
    if df is not None and need:
        strata_cols = _clinical_design(df, meta, notes, C[ok].reset_index(drop=True), model, used, s)
    n_off = coverage(xx, s.round_decimals)["off_modal"]
    r.update(cox_n=int(len(xx)), cox_events=int(ee.sum()),
             psi_iqr=float(np.subtract(*np.percentile(xx, [75, 25]))) if len(xx) else NAN,
             psi_sd=float(np.std(xx, ddof=1)) if len(xx) > 1 else NAN, cox_notes="; ".join(notes))
    low_var_fit = low_psi_variance(xx, s)
    if low_var_fit or df is None or r["cox_n"] < s.cox_min_n or r["cox_events"] < s.cox_min_events or \
            n_off < s.min_off_modal:
        r["cox_status"] = "low_psi_variance" if low_var_fit else (why if df is None else "too_few_patients_or_events")
        return r, []
    fit = fit_cox(df, "psi10", {"gex_z": "expr"} if "gex_z" in df else None, ph=s.ph_test,
                  refit_step_size=s.nr_refit_step_size, strata=strata_cols or None)
    r.update(cox_status=fit["status"], cox_fit_note=fit.get("fit_note", ""), cox_fit_warning=fit.get("fit_warning", ""))
    if fit["status"] != "tested":
        return r, []
    b, se, k = fit["beta"], fit["se"], r["psi_iqr"] / s.psi_step
    iqr = (lambda v: float(np.exp(v * k)) if k > 0 else NAN)         # undefined when the IQR is 0
    k_sd = r["psi_sd"] / s.psi_step
    sd = (lambda v: float(np.exp(v * k_sd)) if k_sd > 0 else NAN)
    r.update(cox_beta=b, cox_se=se, cox_p=fit["p"], hr_per_step=float(np.exp(b)),
             ci_low_step=float(np.exp(b - s.ci_z * se)), ci_high_step=float(np.exp(b + s.ci_z * se)),
             hr_per_iqr=iqr(b), ci_low_iqr=iqr(b - s.ci_z * se), ci_high_iqr=iqr(b + s.ci_z * se),
             hr_per_sd=sd(b), ci_low_sd=sd(b - s.ci_z * se), ci_high_sd=sd(b + s.ci_z * se),
             ph_p=fit.get("ph_p", NAN))
    if "gex_z" in df:
        r.update(expr_hr_per_sd=fit["expr_hr"], expr_p=fit["expr_p"])
    terms = []
    summ = fit["summary"]
    unstable = []
    ph = fit.get("ph_terms") or {}
    for m in meta:
        c, sei, p = (float(summ.loc[m["column"], k_]) for k_ in ("coef", "se(coef)", "p"))
        with np.errstate(over="ignore"):
            terms.append(dict({k_: v for k_, v in m.items() if k_ != "column"}, coef=c, se=sei, hr=float(np.exp(c)),
                              ci_low=float(np.exp(c - s.ci_z * sei)), ci_high=float(np.exp(c + s.ci_z * sei)), p=p,
                              ph_p=ph.get(m["column"], NAN)))
        if m["kind"] == "psi":                              # the tested term per IQR (or per SD: psi_hr_unit)
            u, kk = (("iqr", k) if s.psi_hr_unit == "iqr" else ("sd", k_sd))
            terms.append(dict(term=quantity, kind=f"psi_{u}", level="", reference="",
                              unit=f"{u.upper()} ({r[f'psi_{u}']:.3g} {unit})", coef=c * kk, se=sei * kk,
                              hr=r[f"hr_per_{u}"], ci_low=r[f"ci_low_{u}"], ci_high=r[f"ci_high_{u}"], p=p,
                              ph_p=ph.get(m["column"], NAN)))
        elif sei > 3:                                       # a sparse level or a nearly separated covariate
            unstable.append(f"{m['term']}{' ' + m['level'] if m['level'] else ''}")
    flags = _fit_flags(r, len(meta), s, r["psi_iqr"] if s.narrow_measure == "iqr" else r["psi_sd"], quantity)
    flags += _ph_note(fit, meta, s)
    r["cox_notes"] = "; ".join(x for x in [r.get("cox_notes", "")] + flags
                               + ([f"unstable: {', '.join(unstable)} (se > 3)"] if unstable else []) if x)
    return r, terms


def expression_cell(g_base: np.ndarray, ep_pos: np.ndarray, time: np.ndarray, event: np.ndarray, s: Settings,
                    model: CoxModel | None = None, clinical: pd.DataFrame | None = None) -> tuple[dict, list[dict]]:
    """KM and Cox statistics of host-gene expression itself in one cohort for one endpoint, and the Cox terms.

    The same rules as `survival_cell`, with expression in place of PSI: the coverage gate on the survival samples, a
    median split (expression <= cut is the low arm), and Cox on z(expression) plus the model's clinical covariates and
    strata (`model.expression` is ignored). The HR is per SD of expression in the fit cohort.
    """
    model = model or CoxModel()
    g_base = np.asarray(g_base, float)
    cv = coverage(g_base, s.round_decimals)
    obs = g_base[np.isfinite(g_base)]
    describe = " + ".join(["expression"] + list(model.covariates)) + (
        f"; strata: {', '.join(model.strata)}" if model.strata else "")
    cut, how = km_cut(obs, s.km_split_expression)
    r = dict(n_survival=len(g_base), n_obs=cv["n_obs"], frac_obs=cv["frac_obs"], off_modal=cv["off_modal"],
             modal_share=cv["modal_share"], cutoff=cut, km_split=how,
             eligible=bool(cv["frac_obs"] >= s.coverage_frac and cv["off_modal"] >= s.min_off_modal),
             n_endpoint=len(ep_pos), cox_model=describe)
    if not r["eligible"]:
        return dict(r, km_status="coverage_gate", cox_status="coverage_gate"), []
    keep = np.isfinite(g_base[ep_pos])
    idx = ep_pos[keep]
    t, e, x = np.asarray(time, float)[keep], np.asarray(event, int)[keep], g_base[idx]
    high = x > r["cutoff"]
    flat = len(x) < 2 or not np.ptp(x) > 0
    # ---------------------------------------------------------------- KM
    r.update(km_n=len(x), n_low=int((~high).sum()), n_high=int(high.sum()), events_low=int(e[~high].sum()),
             events_high=int(e[high].sum()))
    if flat:
        r["km_status"] = "constant_expression"
    elif not km_gate(high, e, s):
        r["km_status"] = "too_few_patients_or_events"
    else:
        lr = logrank(t, e, high)
        r.update(km_status="tested", o_high=lr["o_high"], e_high=lr["e_high"], o_low=lr["o_low"], e_low=lr["e_low"],
                 km_var=lr["var"], km_chi2=lr["chi2"], km_p=lr["p"], logrank_hr=lr["logrank_hr"])
        _km_ph(r, t, e, high, s)
    # ---------------------------------------------------------------- Cox
    need = model.clinical_columns
    if need and clinical is None:
        return dict(r, cox_status="no_clinical_table"), []
    ok, notes = np.ones(len(idx), bool), []
    C, used = None, list(need)
    if need:
        C = clinical.iloc[idx].reset_index(drop=True)
        used = _complete_columns(C, need, s, notes)
        for c in used:
            ok &= C[c].notna().to_numpy()
    if len(used) < len(need):
        r["cox_model"] = " + ".join(["expression"] + [c for c in model.covariates if c in used]) + (
            f"; strata: {', '.join(c for c in model.strata if c in used)}" if any(c in used for c in model.strata)
            else "")
    r["cox_n_dropped"] = int((~ok).sum())
    tt, ee, xx = t[ok], e[ok], x[ok]
    sd = float(np.std(xx, ddof=1)) if len(xx) > 1 else 0.0
    r.update(cox_n=int(len(xx)), cox_events=int(ee.sum()), expr_sd=sd,
             expr_iqr=float(np.subtract(*np.percentile(xx, [75, 25]))) if len(xx) else NAN)
    if not sd > 0:
        r.update(cox_status="constant_expression", cox_notes="; ".join(notes))
        return r, []
    df = pd.DataFrame({"time": tt, "event": ee, "gex_z": (xx - xx.mean()) / sd})
    meta = [dict(column="gex_z", term="expression", kind="gex", level="", reference="", unit="SD")]
    strata_cols = _clinical_design(df, meta, notes, C[ok].reset_index(drop=True), model, used, s) if need else []
    r["cox_notes"] = "; ".join(notes)
    if r["cox_n"] < s.cox_min_n or r["cox_events"] < s.cox_min_events:
        r["cox_status"] = "too_few_patients_or_events"
        return r, []
    fit = fit_cox(df, "gex_z", None, ph=s.ph_test, refit_step_size=s.nr_refit_step_size, strata=strata_cols or None)
    r.update(cox_status=fit["status"], cox_fit_note=fit.get("fit_note", ""), cox_fit_warning=fit.get("fit_warning", ""))
    if fit["status"] != "tested":
        return r, []
    b, se = fit["beta"], fit["se"]
    r.update(cox_beta=b, cox_se=se, cox_p=fit["p"], hr_per_sd=float(np.exp(b)), ci_low_sd=float(np.exp(b - s.ci_z * se)),
             ci_high_sd=float(np.exp(b + s.ci_z * se)), ph_p=fit.get("ph_p", NAN))
    terms, unstable = _terms(fit, meta, s)
    flags = _fit_flags(r, len(meta), s) + _ph_note(fit, meta, s)
    r["cox_notes"] = "; ".join(x for x in [r["cox_notes"]] + flags
                               + ([f"unstable: {', '.join(unstable)} (se > 3)"] if unstable else []) if x)
    return r, terms
