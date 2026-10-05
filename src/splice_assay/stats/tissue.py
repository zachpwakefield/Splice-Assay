"""Case vs reference per cohort (e.g. tumour vs normal): paired (within patients) and unpaired (all samples), each
with a hit call.

All effects are case minus reference.

paired     exact signed-rank test of the within-patient differences (zeros dropped); delta = median difference;
           needs `min_pairs` pairs
unpaired   Mann-Whitney test of all case against all reference patients, one value per patient and group (a patient
           with several samples in a group, e.g. replicate aliquots, enters with their mean, so that each patient
           counts once per group); delta = median(case) - median(reference); needs `min_group` patients of each. The
           paired test keeps the pairs' own samples
0/1 check  the test is repeated after dropping PSI values of exactly 0 or 1 (for pairs: pairs where either value is;
           unpaired: sample values, before the per-patient mean); state 'not_applicable' (nothing dropped),
           'untestable' (the reduced set fails the size gate), 'pass' (p < alpha, same sign, |delta| > min_abs_delta)
           or 'fail' (otherwise, also when the reduced set shows no difference at all)
hit        p < alpha, round(|delta|, 12) > min_abs_delta, sign(HL) = sign(delta) and a 0/1 state of pass or
           not_applicable
matched / other   composition diagnostics on the unpaired side: the case values of paired patients against all
           reference values ('matched'), and the remaining case patients against the matched ones ('other')
within-patient support   the paired test was possible and either the paired design or the matched comparison is a hit
no within-patient check  a hit in a cohort where the paired test was not possible
Values are rounded to round_decimals before they are ranked or compared (floating-point noise only).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from ..config import Settings
from .core import ROUND, boundary, exact_signed_rank, hl_one_sample, hl_two_sample, mann_whitney

NAN = float("nan")


def paired_row(d: np.ndarray, min_pairs: int, round_decimals: int) -> dict:
    r = dict(n_pairs=len(d))
    if len(d) < max(min_pairs, 1):
        return dict(r, status="too_few_pairs")
    dr = np.round(d, round_decimals)
    if np.all(dr == 0):
        return dict(r, status="all_differences_zero", delta_median=0.0)
    p, n_nonzero = exact_signed_rank(d, round_decimals)
    return dict(r, status="tested", n_nonzero=n_nonzero, delta_median=float(np.median(d)), hl=hl_one_sample(d),
                mean_delta=float(np.mean(d)), n_up=int((dr > 0).sum()), n_down=int((dr < 0).sum()), p=p)


def unpaired_row(x: np.ndarray, y: np.ndarray, min_group: int, round_decimals: int = ROUND) -> dict:
    """The Mann-Whitney comparison of two sets of values (one per patient). Ties and the constant check use the
    values rounded to `round_decimals`; the medians and HL use the values as given."""
    r = dict(n_case=len(x), n_reference=len(y))
    if len(x) < max(min_group, 1) or len(y) < max(min_group, 1):
        return dict(r, status="too_few_samples")
    xr, yr = np.round(x, round_decimals), np.round(y, round_decimals)
    if np.ptp(np.concatenate([xr, yr])) == 0:
        return dict(r, status="constant")
    u = mann_whitney(xr, yr)
    return dict(r, status="tested", median_case=float(np.median(x)), median_reference=float(np.median(y)),
                delta_median=float(np.median(x) - np.median(y)), hl=hl_two_sample(x, y),
                auc=float(u.statistic / (len(x) * len(y))), p=float(u.pvalue))


def per_patient(x: np.ndarray, patients=None) -> np.ndarray:
    """One value per patient: the mean of each patient's values in `x` (aligned with `patients`), in the order of their
    first value; `x` itself when every patient has one value or `patients` is None."""
    if patients is None or not len(x):
        return x
    codes, uniq = pd.factorize(np.asarray(patients, dtype=object))
    if len(uniq) == len(x):
        return x
    return np.bincount(codes, weights=x, minlength=len(uniq)) / np.bincount(codes, minlength=len(uniq))


def patient_values(values, samples, patients) -> pd.DataFrame:
    """The values the unpaired test compares for a group, with their patients: the observed values of `samples`
    (sample IDs; `values` maps sample_id -> value), one per patient (`per_patient`; `patients` aligned with
    `samples`). Columns patient_id, value and n_samples (the observed samples behind the value), in the order of each
    patient's first observed sample."""
    x = values.reindex(samples).to_numpy(float)
    ok = np.isfinite(x)
    p = np.asarray(patients, dtype=object)[ok]
    codes, uniq = pd.factorize(p)
    return pd.DataFrame({"patient_id": np.asarray(uniq, dtype=object), "value": per_patient(x[ok], p),
                         "n_samples": np.bincount(codes, minlength=len(uniq))})


def state_01(main: dict, r01: dict | None, n_dropped: int, s: Settings) -> str:
    """Robustness to dropping PSI 0/1 of a tested row, given the row of the reduced set."""
    if main.get("status") != "tested":
        return "main_not_tested"
    if n_dropped == 0:
        return "not_applicable"
    if r01 is not None and r01.get("status") in ("all_differences_zero", "constant"):
        return "fail"                                  # enough values remain, and they show no difference
    if r01 is None or r01.get("status") != "tested":
        return "untestable"
    d, d1 = main["delta_median"], r01["delta_median"]
    ok = (r01["p"] < s.alpha) and np.sign(d1) == np.sign(d) and np.round(abs(d1), s.round_decimals) > s.min_abs_delta
    return "pass" if ok else "fail"


def hit_status(r: dict, s: Settings) -> str:
    """'hit' or the first reason it is not one (not_tested > not_significant > small_effect > direction_ambiguous >
    robustness_fail / robustness_untestable)."""
    if r.get("status") != "tested":
        return "not_tested"
    if not r["p"] < s.alpha:
        return "not_significant"
    if not np.round(abs(r["delta_median"]), s.round_decimals) > s.min_abs_delta:
        return "small_effect"
    if np.sign(r["delta_median"]) != np.sign(r["hl"]):
        return "direction_ambiguous"
    if r.get("state_01") == "fail":
        return "robustness_fail"
    if r.get("state_01") == "untestable":
        return "robustness_untestable"
    return "hit"


def _paired(t: np.ndarray, n: np.ndarray, s: Settings) -> dict:
    r = paired_row(t - n, s.min_pairs, s.round_decimals)
    if not s.robust_01:
        return dict(r, n_dropped_0_1=0, state_01=state_01(r, None, 0, s))
    keep = ~boundary(t, s.round_decimals) & ~boundary(n, s.round_decimals)
    r01 = paired_row(t[keep] - n[keep], s.min_pairs, s.round_decimals) if (~keep).any() else None
    nd = int((~keep).sum())
    return dict(r, n_dropped_0_1=nd, state_01=state_01(r, r01, nd, s), delta_median_01=(r01 or {}).get(
        "delta_median", NAN), p_01=(r01 or {}).get("p", NAN))


def _unpaired(x: np.ndarray, y: np.ndarray, s: Settings, px=None, py=None) -> dict:
    """The unpaired design on the observed sample values x and y, whose patients are px and py (None: one patient per
    value): the test compares one value per patient; the 0/1 check drops sample values first, then takes the means."""
    r = unpaired_row(per_patient(x, px), per_patient(y, py), s.min_group, s.round_decimals)
    if not s.robust_01:
        return dict(r, n_dropped_0_1=0, state_01=state_01(r, None, 0, s))
    bx, by = boundary(x, s.round_decimals), boundary(y, s.round_decimals)
    nd = int(bx.sum() + by.sum())
    sub = (lambda p, keep: None if p is None else np.asarray(p, dtype=object)[keep])  # noqa: E731
    r01 = unpaired_row(per_patient(x[~bx], sub(px, ~bx)), per_patient(y[~by], sub(py, ~by)), s.min_group,
                       s.round_decimals) if nd else None
    return dict(r, n_dropped_0_1=nd, state_01=state_01(r, r01, nd, s),
                delta_median_01=(r01 or {}).get("delta_median", NAN), p_01=(r01 or {}).get("p", NAN))


def matched_hit_status(m: dict, unpaired_delta: float, s: Settings) -> str:
    """The matched comparison is a hit in the unpaired hit's direction (not_tested > not_significant > small_effect >
    other_direction_than_unpaired > direction_ambiguous > robustness_fail_or_untestable)."""
    if m.get("status") != "tested":
        return "not_tested"
    dm = np.round(m["delta_median"], s.round_decimals)
    if not m["p"] < s.alpha:
        return "not_significant"
    if not abs(dm) > s.min_abs_delta:
        return "small_effect"
    if np.sign(dm) != np.sign(np.round(unpaired_delta, s.round_decimals)):
        return "other_direction_than_unpaired"
    if np.sign(m["hl"]) != np.sign(dm):
        return "direction_ambiguous"
    if m["state_01"] not in ("pass", "not_applicable"):
        return "robustness_fail_or_untestable"
    return "hit"


def compare_groups(values, pair_c, pair_r, cases, refs, s: Settings, case_patients=None,
                   ref_patients=None) -> dict:
    """All case-vs-reference statistics of one event in one cohort.

    `values` maps sample_id -> PSI (a pandas Series; NaN = missing); the other arguments are arrays of sample IDs:
    the case and reference sample of each pair (aligned), every case and every reference sample of the cohort.
    `case_patients` and `ref_patients` (aligned with `cases` and `refs`) name each sample's patient, so that the
    unpaired design compares one value per patient; None treats every sample as its own patient.
    """
    get = lambda ids: values.reindex(ids).to_numpy(float)  # noqa: E731
    if len(refs) == 0:
        out = dict(paired_status="no_reference_samples", unpaired_status="no_reference_samples")
        out.update(paired_hit_status="not_tested", paired_hit=False, unpaired_hit_status="not_tested",
                   unpaired_hit=False, group_hit=False, matched_hit_status="not_tested",
                   no_within_patient_check=False,
                   within_patient_support=False, composition_sensitive=False)
        return out
    t, n = get(pair_c), get(pair_r)
    ok = np.isfinite(t) & np.isfinite(n)
    pr = _paired(t[ok], n[ok], s)
    xs, ys = get(cases), get(refs)
    okx, oky = np.isfinite(xs), np.isfinite(ys)
    x, y = xs[okx], ys[oky]
    pc = None if case_patients is None else np.asarray(case_patients, dtype=object)
    py = None if ref_patients is None else np.asarray(ref_patients, dtype=object)[oky]
    ur = _unpaired(x, y, s, None if pc is None else pc[okx], py)
    out = {f"paired_{k}": v for k, v in pr.items()}
    out.update({f"unpaired_{k}": v for k, v in ur.items()})
    out.update(unpaired_n_case_samples=len(x), unpaired_n_reference_samples=len(y))
    # composition diagnostics: the case values of paired patients (matched) vs all reference values; the other case
    # patients vs the matched ones
    if pc is None:
        mt = get(pair_c)
        mt, pm = mt[np.isfinite(mt)], None
        other_ids = np.setdiff1d(np.asarray(cases, dtype=object).astype(str),
                                 np.asarray(pair_c, dtype=object).astype(str))
        ot, po = get(other_ids), None
        ot = ot[np.isfinite(ot)]
    else:
        of = dict(zip(np.asarray(cases, dtype=object).astype(str), pc))
        paired = {of.get(str(c)) for c in np.asarray(pair_c, dtype=object)}
        in_pair = np.array([p in paired for p in pc], dtype=bool)
        mt, pm = xs[in_pair & okx], pc[in_pair & okx]
        ot, po = xs[~in_pair & okx], pc[~in_pair & okx]
    m = _unpaired(mt, y, s, pm, py)
    mv, ov = per_patient(mt, pm), per_patient(ot, po)
    o = unpaired_row(ov, mv, s.min_group, s.round_decimals)
    out.update(matched_status=m["status"], matched_n_case=len(mv), matched_n_reference=m["n_reference"],
               matched_delta_median=m.get("delta_median", NAN), matched_hl=m.get("hl", NAN),
               matched_p=m.get("p", NAN), matched_state_01=m["state_01"], other_status=o["status"],
               other_n_case=len(ov), other_n_matched=len(mv), other_delta_median=o.get("delta_median", NAN),
               other_hl=o.get("hl", NAN), other_p=o.get("p", NAN))
    # hits
    out["paired_hit_status"] = hit_status(pr, s)
    out["unpaired_hit_status"] = hit_status(ur, s)
    out["paired_hit"] = out["paired_hit_status"] == "hit"
    out["unpaired_hit"] = out["unpaired_hit_status"] == "hit"
    out["group_hit"] = out["paired_hit"] or out["unpaired_hit"]
    out["matched_hit_status"] = matched_hit_status(m, ur.get("delta_median", NAN), s)
    pt, ut = pr["status"] == "tested", ur["status"] == "tested"
    out["no_within_patient_check"] = bool(out["group_hit"] and not pt)
    out["within_patient_support"] = bool(pt and (out["paired_hit"] or out["matched_hit_status"] == "hit"))
    dp, du = (np.round(r.get("delta_median", NAN), s.round_decimals) for r in (pr, ur))
    out["composition_sensitive"] = bool(pt and ut and (np.round(abs(dp - du), s.round_decimals) > s.min_abs_delta
                                                       or np.sign(dp) != np.sign(du)))
    return out
