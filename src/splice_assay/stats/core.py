"""Core statistics.

The arithmetic follows the locked analysis these figures come from (docs/methods.md) line for line, so its tables are
reproduced exactly; tests check every function against an independent implementation (brute-force enumeration,
scipy or lifelines).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

ROUND = 12  # decimals kept before PSI values or differences are ranked or compared (floating-point noise only)


def boundary(x, round_decimals: int = ROUND) -> np.ndarray:
    """PSI exactly 0 or 1 after rounding."""
    r = np.round(np.asarray(x, dtype=float), round_decimals)
    return (r == 0) | (r == 1)


# ============================================================================================ tumour vs normal
def exact_signed_rank(differences, round_decimals: int = ROUND) -> tuple[float, int]:
    """Two-sided exact conditional Wilcoxon signed-rank p (tie-aware, zeros dropped). Returns (p, n_nonzero).

    Doubled average ranks are integers, so a probability-mass dynamic program enumerates every sign assignment;
    ties are handled exactly, not by a normal approximation.
    """
    d = np.round(np.asarray(differences, dtype=float), round_decimals)
    d = d[d != 0]
    if not len(d):
        return 1.0, 0
    w = np.rint(2 * stats.rankdata(np.abs(d), method="average")).astype(int)
    total, pos = int(w.sum()), int(w[d > 0].sum())
    stat = min(pos, total - pos)
    pmf = np.array([1.0])
    for wi in w:
        new = np.zeros(len(pmf) + wi)
        new[: len(pmf)] += 0.5 * pmf
        new[wi:] += 0.5 * pmf
        pmf = new
    return min(1.0, 2 * float(pmf[: stat + 1].sum())), len(d)


def hl_one_sample(d) -> float:
    """Hodges-Lehmann pseudo-median: median of the Walsh averages (d_i + d_j) / 2, i <= j."""
    d = np.asarray(d, dtype=float)
    i, j = np.triu_indices(len(d))
    return float(np.median((d[i] + d[j]) / 2))


def hl_two_sample(x, y) -> float:
    """Hodges-Lehmann shift: median of all pairwise differences x_i - y_j."""
    return float(np.median(np.subtract.outer(np.asarray(x, float), np.asarray(y, float))))


def mann_whitney(x, y):
    """Two-sided Mann-Whitney U test, normal approximation with continuity correction (tie-corrected)."""
    return stats.mannwhitneyu(x, y, alternative="two-sided", method="asymptotic", use_continuity=True)


# ============================================================================================ survival
def logrank(time, event, high, strata=None) -> dict:
    """Two-group log-rank, optionally stratified (sums of the within-stratum O - E and variance).

    `high` is a boolean array. The log-rank HR is (O/E of the high group) / (O/E of the low group).
    """
    time, event, high = (np.asarray(a) for a in (time, event, high))
    strata = np.zeros(len(time), dtype=int) if strata is None else pd.factorize(np.asarray(strata))[0]
    O = {1: 0.0, 0: 0.0}
    E = {1: 0.0, 0: 0.0}
    V = 0.0
    for s in np.unique(strata):
        m = strata == s
        t, e, g = time[m], event[m].astype(bool), high[m].astype(bool)
        ts = np.sort(t)
        t1 = np.sort(t[g])
        ut, d = np.unique(t[e], return_counts=True)
        if not len(ut):
            continue
        d1 = pd.Series(g[e].astype(int)).groupby(t[e]).sum().reindex(ut).to_numpy()
        n = len(ts) - np.searchsorted(ts, ut, side="left")
        n1 = len(t1) - np.searchsorted(t1, ut, side="left")
        e1 = d * n1 / n
        v = np.where(n > 1, d * (n1 / n) * (1 - n1 / n) * (n - d) / np.maximum(n - 1, 1), 0.0)
        O[1] += d1.sum()
        E[1] += e1.sum()
        O[0] += (d - d1).sum()
        E[0] += (d - e1).sum()
        V += v.sum()
    chi2 = (O[1] - E[1]) ** 2 / V if V > 0 else np.nan
    hr = (O[1] / E[1]) / (O[0] / E[0]) if min(E.values()) > 0 and O[0] > 0 else np.nan
    return dict(o_high=O[1], e_high=E[1], o_low=O[0], e_low=E[0], var=V, chi2=chi2,
                p=float(stats.chi2.sf(chi2, 1)) if np.isfinite(chi2) else np.nan, logrank_hr=hr)


def km_curve(time, event) -> pd.DataFrame:
    """Product-limit estimate at each distinct time, with the number at risk just before it."""
    time, event = np.asarray(time, float), np.asarray(event, int)
    ut = np.unique(time)
    ts = np.sort(time)
    at_risk = len(ts) - np.searchsorted(ts, ut, side="left")
    deaths = pd.Series(event).groupby(time).sum().reindex(ut).to_numpy()
    surv = np.cumprod(1 - deaths / at_risk)
    return pd.DataFrame({"time": ut, "at_risk": at_risk, "events": deaths, "surv": surv})


def km_with_ci(time, event, z: float = 1.96) -> pd.DataFrame:
    """km_curve plus a log-log (exponential Greenwood) confidence band."""
    cv = km_curve(time, event)
    n, d, s = cv.at_risk.to_numpy(float), cv.events.to_numpy(float), cv.surv.to_numpy()
    g = np.cumsum(np.where(n > d, d / (n * np.maximum(n - d, 1)), 0.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        ll, se = np.log(-np.log(s)), np.sqrt(g) / np.abs(np.log(s))
        lo, hi = np.exp(-np.exp(ll + z * se)), np.exp(-np.exp(ll - z * se))
    cv["lo"], cv["hi"] = np.where(np.isfinite(lo), lo, s), np.where(np.isfinite(hi), hi, s)
    return cv


# ============================================================================================ coverage
def coverage(values, round_decimals: int = ROUND) -> dict:
    """Observed fraction, the count of observed values away from the exact modal value, and the spread."""
    values = np.asarray(values, dtype=float)
    v = values[np.isfinite(values)]
    if not len(v):
        return dict(n_obs=0, frac_obs=0.0, off_modal=0, modal_share=np.nan, q10=np.nan, q90=np.nan)
    vals, counts = np.unique(np.round(v, round_decimals), return_counts=True)
    return dict(n_obs=len(v), frac_obs=len(v) / len(values), off_modal=int(len(v) - counts.max()),
                modal_share=float(counts.max() / len(v)), q10=float(np.quantile(v, .1)), q90=float(np.quantile(v, .9)))
