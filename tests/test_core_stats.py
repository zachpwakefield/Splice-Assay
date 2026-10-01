"""Each core statistic against an independent implementation."""
import itertools

import numpy as np
import pytest
from scipy import stats

from splice_assay.stats.core import (coverage, exact_signed_rank, hl_one_sample, hl_two_sample, km_curve, km_with_ci,
                                     logrank, mann_whitney)


def brute_signed_rank(d):
    """Two-sided exact signed-rank p by enumerating all sign assignments (average ranks, zeros dropped)."""
    d = np.round(np.asarray(d, float), 12)
    d = d[d != 0]
    r = stats.rankdata(np.abs(d))
    w_obs = r[d > 0].sum()
    total = r.sum()
    sums = np.array([sum(ri for ri, s in zip(r, signs) if s) for signs in itertools.product([0, 1], repeat=len(d))])
    stat = min(w_obs, total - w_obs)
    return min(1.0, 2 * float(np.mean(sums <= stat + 1e-9)))


@pytest.mark.parametrize("seed", range(12))
def test_signed_rank_matches_enumeration_with_ties_and_zeros(seed):
    rng = np.random.default_rng(seed)
    n = rng.integers(3, 11)
    d = rng.choice([-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.25, 0.3], size=n)   # ties and zeros on purpose
    p, n_nonzero = exact_signed_rank(d)
    assert n_nonzero == int((d != 0).sum())
    if n_nonzero == 0:
        assert p == 1.0
    else:
        assert p == pytest.approx(brute_signed_rank(d), abs=1e-12)


@pytest.mark.parametrize("seed", range(5))
def test_signed_rank_matches_scipy_without_ties(seed):
    d = np.random.default_rng(seed).normal(0.05, 0.1, 15)
    p, _ = exact_signed_rank(d)
    assert p == pytest.approx(stats.wilcoxon(d, method="exact").pvalue, rel=1e-10)


def test_signed_rank_ignores_float_noise():
    d = np.array([0.127 - 0.027, 0.1, -0.2, 0.3, 1e-15])       # 0.1 twice (one with noise); 1e-15 counts as 0
    p, n = exact_signed_rank(d)
    assert n == 4
    assert p == pytest.approx(brute_signed_rank([0.1, 0.1, -0.2, 0.3]), abs=1e-12)


def test_hodges_lehmann():
    rng = np.random.default_rng(3)
    d = rng.normal(size=9)
    walsh = [(d[i] + d[j]) / 2 for i in range(9) for j in range(i, 9)]
    assert hl_one_sample(d) == pytest.approx(np.median(walsh))
    x, y = rng.normal(size=7), rng.normal(size=5)
    assert hl_two_sample(x, y) == pytest.approx(np.median([a - b for a in x for b in y]))


def test_mann_whitney_is_scipy_asymptotic_with_continuity():
    rng = np.random.default_rng(4)
    x, y = np.round(rng.uniform(size=30), 2), np.round(rng.uniform(size=25), 2)
    ref = stats.mannwhitneyu(x, y, alternative="two-sided", method="asymptotic", use_continuity=True)
    assert mann_whitney(x, y).pvalue == ref.pvalue


def _surv(seed, n=80):
    rng = np.random.default_rng(seed)
    t = np.round(rng.exponential(5, n), 1) + 0.1                      # ties on purpose
    e = rng.uniform(size=n) < 0.7
    g = rng.uniform(size=n) < 0.5
    return t, e.astype(int), g


@pytest.mark.parametrize("seed", range(5))
def test_logrank_matches_lifelines(seed):
    from lifelines.statistics import logrank_test
    t, e, g = _surv(seed)
    r = logrank(t, e, g)
    ref = logrank_test(t[g], t[~g], e[g], e[~g])
    assert r["p"] == pytest.approx(ref.p_value, rel=1e-9)
    assert r["chi2"] == pytest.approx(ref.test_statistic, rel=1e-9)
    assert r["o_high"] + r["o_low"] == e.sum()
    assert r["logrank_hr"] == pytest.approx((r["o_high"] / r["e_high"]) / (r["o_low"] / r["e_low"]))


def test_stratified_logrank_sums_strata():
    t, e, g = _surv(9, 120)
    st = np.repeat(["a", "b"], 60)
    a, b = logrank(t[:60], e[:60], g[:60]), logrank(t[60:], e[60:], g[60:])
    r = logrank(t, e, g, strata=st)
    for k in ("o_high", "e_high", "var"):
        assert r[k] == pytest.approx(a[k] + b[k])


@pytest.mark.parametrize("seed", range(3))
def test_km_and_band_match_lifelines(seed):
    from lifelines import KaplanMeierFitter
    t, e, _ = _surv(seed)
    z = stats.norm.ppf(0.975)
    cv = km_with_ci(t, e, z)
    kmf = KaplanMeierFitter().fit(t, e)
    sf = kmf.survival_function_.iloc[:, 0].reindex(cv.time)
    assert np.allclose(cv.surv, sf.to_numpy(), atol=1e-12)
    ci = kmf.confidence_interval_.reindex(cv.time)
    inner = (cv.surv > 0) & (cv.surv < 1)
    lo = ci[[c for c in ci.columns if "lower" in c][0]].to_numpy()
    hi = ci[[c for c in ci.columns if "upper" in c][0]].to_numpy()
    assert np.allclose(cv.lo[inner], lo[inner], atol=1e-10)
    assert np.allclose(cv.hi[inner], hi[inner], atol=1e-10)
    assert (km_curve(t, e).at_risk.iloc[0]) == len(t)


def test_coverage():
    v = np.array([0.5, 0.5, 0.5, 0.2, np.nan, 0.5000000000001])
    c = coverage(v)
    assert c["n_obs"] == 5 and c["frac_obs"] == pytest.approx(5 / 6)
    assert c["off_modal"] == 1 and c["modal_share"] == pytest.approx(0.8)
    assert coverage(np.array([np.nan, np.nan]))["n_obs"] == 0
