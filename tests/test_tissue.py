import numpy as np
import pandas as pd
import pytest

from splice_assay import Settings
from splice_assay.stats.core import exact_signed_rank
from splice_assay.stats.tissue import compare_groups, hit_status, paired_row, state_01, unpaired_row

S = Settings()


def test_paired_row_states():
    assert paired_row(np.zeros(5), 10, 12)["status"] == "too_few_pairs"
    assert paired_row(np.zeros(12), 10, 12)["status"] == "all_differences_zero"
    d = np.linspace(-0.1, 0.4, 12)
    r = paired_row(d, 10, 12)
    assert r["status"] == "tested" and r["p"] == exact_signed_rank(d)[0]
    assert r["delta_median"] == pytest.approx(np.median(d))


def test_unpaired_row_states():
    assert unpaired_row(np.ones(5), np.ones(20), 10)["status"] == "too_few_samples"
    assert unpaired_row(np.ones(20), np.ones(20), 10)["status"] == "constant"
    x, y = np.linspace(0.3, 0.9, 20), np.linspace(0.1, 0.6, 15)
    r = unpaired_row(x, y, 10)
    assert r["delta_median"] == pytest.approx(np.median(x) - np.median(y))


@pytest.mark.parametrize("row, expected", [
    (dict(status="too_few_pairs"), "not_tested"),
    (dict(status="tested", p=0.2, delta_median=0.3, hl=0.3, state_01="pass"), "not_significant"),
    (dict(status="tested", p=0.01, delta_median=0.1, hl=0.1, state_01="pass"), "small_effect"),
    (dict(status="tested", p=0.01, delta_median=0.127 - 0.027, hl=0.1, state_01="pass"), "small_effect"),
    (dict(status="tested", p=0.01, delta_median=0.2, hl=-0.01, state_01="pass"), "direction_ambiguous"),
    (dict(status="tested", p=0.01, delta_median=0.2, hl=0.2, state_01="fail"), "robustness_fail"),
    (dict(status="tested", p=0.01, delta_median=0.2, hl=0.2, state_01="untestable"), "robustness_untestable"),
    (dict(status="tested", p=0.01, delta_median=-0.2, hl=-0.1, state_01="not_applicable"), "hit"),
])
def test_hit_rule(row, expected):
    assert hit_status(row, S) == expected


def test_state_01():
    main = dict(status="tested", delta_median=0.2)
    assert state_01(main, None, 0, S) == "not_applicable"
    assert state_01(main, None, 3, S) == "untestable"
    assert state_01(main, dict(status="tested", delta_median=0.15, p=0.01), 3, S) == "pass"
    assert state_01(main, dict(status="tested", delta_median=-0.15, p=0.01), 3, S) == "fail"
    assert state_01(dict(status="too_few_pairs"), None, 0, S) == "main_not_tested"


def test_zero_one_values_make_the_check_untestable():
    ids_t = [f"t{i}" for i in range(12)]
    ids_n = [f"n{i}" for i in range(12)]
    v = pd.Series([1.0] * 4 + [0.9] * 8 + [0.2] * 12, index=ids_t + ids_n)   # 4 tumours at exactly 1
    v.iloc[:12] = v.iloc[:12] + np.r_[np.zeros(4), np.linspace(0, 0.05, 8)]
    r = compare_groups(v, np.array(ids_t), np.array(ids_n), np.array(ids_t), np.array(ids_n), S)
    assert r["paired_n_dropped_0_1"] == 4 and r["paired_state_01"] == "untestable"
    assert r["paired_hit_status"] == "robustness_untestable"


def test_synthetic_cohorts(results):
    tn = results.groups.set_index(["event_id", "cohort"])
    se = tn.loc[("SYN1:SE:1", "COH1")]
    assert se.paired_hit and se.unpaired_hit and se.within_patient_support
    assert tn.loc[("SYN1:SE:1", "COH2")].paired_status == "too_few_pairs"
    assert tn.loc[("SYN1:SE:1", "COH2")].no_within_patient_check          # unpaired hit, 6 pairs
    assert not tn.loc[("SYN1:SE:1", "COH1")].no_within_patient_check
    assert not tn.loc[("SYN1:SE:1", "COH4")].no_within_patient_check      # no normals, no hit
    assert tn.loc[("SYN1:SE:1", "COH4")].paired_status == "no_reference_samples"
    assert not tn.loc[("SYN1:SE:1", "COH3")].group_hit


def test_results_do_not_depend_on_sample_order(tables):
    import splice_assay as sa
    shuffled = {k: (v.sample(frac=1, random_state=1) if k in ("samples", "psi", "survival") else v)
                for k, v in tables.items()}
    a = sa.analyze(sa.Dataset.from_tables(**tables), events=["SYN1:SE:1"])
    b = sa.analyze(sa.Dataset.from_tables(**shuffled), events=["SYN1:SE:1"])
    cols = ["paired_p", "paired_delta_median", "unpaired_p", "unpaired_delta_median", "paired_hit_status"]
    pd.testing.assert_frame_equal(a.groups[cols], b.groups[cols])
    num = ["cutoff", "km_p", "logrank_hr", "cox_beta", "cox_se", "hr_per_iqr"]
    assert np.allclose(a.survival[num].to_numpy(float), b.survival[num].to_numpy(float), rtol=1e-9, equal_nan=True)


def _means(v, samples, patients):
    """Per-patient means, computed independently of the package (pandas groupby)."""
    return pd.Series(v[list(samples)].to_numpy(), index=list(patients)).groupby(level=0, sort=False).mean().to_numpy()


def test_unpaired_design_compares_one_value_per_patient():
    """A patient's several samples in one group enter the unpaired test as their mean: the gate counts patients, and
    p is Mann-Whitney (scipy) on the per-patient means."""
    from scipy.stats import mannwhitneyu
    rng = np.random.default_rng(3)
    case = [(f"P{i}", f"P{i}-T{k}") for i in range(14) for k in range(1 + i % 3)]      # 1 to 3 samples per patient
    ref = [(f"R{i}", f"R{i}-N{k}") for i in range(12) for k in range(1 + i % 2)]
    vals = np.r_[rng.uniform(0.5, 0.9, len(case)), rng.uniform(0.3, 0.7, len(ref))]
    v = pd.Series(vals, index=[sid for _, sid in case + ref])
    (cp, cs), (rp, rs) = (map(np.array, zip(*case)), map(np.array, zip(*ref)))
    r = compare_groups(v, np.array([]), np.array([]), cs, rs, S, cp, rp)
    x, y = _means(v, cs, cp), _means(v, rs, rp)
    assert (r["unpaired_n_case"], r["unpaired_n_reference"]) == (14, 12)
    assert (r["unpaired_n_case_samples"], r["unpaired_n_reference_samples"]) == (len(cs), len(rs))
    want = mannwhitneyu(x, y, alternative="two-sided", method="asymptotic", use_continuity=True).pvalue
    assert r["unpaired_p"] == pytest.approx(want, rel=1e-12)
    assert r["unpaired_delta_median"] == pytest.approx(np.median(x) - np.median(y))
    # ten samples of each of two patients are two patients: not enough for a test
    t_ids, n_ids = np.array([f"t{i}" for i in range(10)]), np.array([f"n{i}" for i in range(10)])
    two = pd.Series(np.r_[0.8 + 0.001 * np.arange(10), 0.4 + 0.001 * np.arange(10)], index=np.r_[t_ids, n_ids])
    r2 = compare_groups(two, np.array([]), np.array([]), t_ids, n_ids, S, np.array(["A"] * 10), np.array(["B"] * 10))
    assert r2["unpaired_status"] == "too_few_samples" and r2["unpaired_n_case"] == 1 and not r2["group_hit"]
    r3 = compare_groups(two, np.array([]), np.array([]), t_ids, n_ids, S)       # without patients: twenty samples
    assert r3["unpaired_status"] == "tested"


def test_unpaired_0_1_check_drops_sample_values_before_the_means():
    """P0's tumour samples are 1.0 and 0.8: the check drops the 1.0 and keeps P0 at 0.8."""
    from scipy.stats import mannwhitneyu
    cs = np.array(["P0-a", "P0-b"] + [f"P{i}" for i in range(1, 12)])
    cp = np.array(["P0", "P0"] + [f"P{i}" for i in range(1, 12)])
    rs = np.array([f"R{i}" for i in range(12)])
    v = pd.Series(np.r_[1.0, 0.8, np.linspace(0.6, 0.75, 11), np.linspace(0.3, 0.5, 12)], index=np.r_[cs, rs])
    r = compare_groups(v, np.array([]), np.array([]), cs, rs, S, cp, rs)
    assert r["unpaired_n_dropped_0_1"] == 1 and r["unpaired_n_case"] == 12
    x01 = np.r_[0.8, np.linspace(0.6, 0.75, 11)]
    assert r["unpaired_p_01"] == pytest.approx(mannwhitneyu(x01, np.linspace(0.3, 0.5, 12), alternative="two-sided",
                                                            method="asymptotic").pvalue, rel=1e-12)


def test_patient_values_keep_their_patients():
    from splice_assay.stats.tissue import patient_values
    v = pd.Series([0.2, 0.4, 0.9, np.nan, 0.5], index=["a1", "a2", "b1", "c1", "d1"])
    t = patient_values(v, ["a1", "b1", "a2", "c1", "d1"], ["A", "B", "A", "C", "D"])
    assert t.patient_id.tolist() == ["A", "B", "D"] and t.n_samples.tolist() == [2, 1, 1]
    assert np.allclose(t.value, [0.3, 0.9, 0.5])


def test_patients_change_nothing_with_one_sample_each(ds):
    """The synthetic data have one sample per patient and group: naming the patients gives the same statistics."""
    from splice_assay.analysis import GroupDesign
    gd = GroupDesign.from_dataset(ds)
    for e in ("SYN1:SE:1", "SYN2:MXE:1"):
        v = ds.psi.loc[e]
        for c in ds.cohorts:
            d = gd.cohorts[c]
            a = compare_groups(v, d["pair_c"], d["pair_r"], d["cases"], d["refs"], S)
            pd.testing.assert_series_equal(pd.Series(a, dtype=object), pd.Series(gd.compare(v, c, S), dtype=object))


def test_values_equal_after_rounding_are_ties():
    """Floating-point noise neither makes a difference nor counts as one."""
    assert unpaired_row(np.full(12, 0.1 + 0.2), np.full(12, 0.3), 10)["status"] == "constant"
    r = paired_row(np.r_[np.full(5, 0.3 - (0.1 + 0.2)), np.full(6, 0.2)], 10, 12)
    assert (r["n_nonzero"], r["n_up"], r["n_down"]) == (6, 6, 0)


def test_empty_groups_are_not_tested_whatever_the_minimum():
    s0 = S.replace(min_group=0, min_pairs=0)
    v = pd.Series([np.nan, np.nan], index=["t1", "n1"])
    r = compare_groups(v, np.array(["t1"]), np.array(["n1"]), np.array(["t1"]), np.array(["n1"]), s0)
    assert r["unpaired_status"] == "too_few_samples" and r["paired_status"] == "too_few_pairs"


def test_state_01_when_the_reduced_set_shows_no_difference():
    """Enough values remain but none differ: the effect does not survive (fail); too few remain: untestable."""
    main_row = dict(status="tested", delta_median=0.2)
    assert state_01(main_row, dict(status="all_differences_zero", delta_median=0.0), 3, S) == "fail"
    assert state_01(main_row, dict(status="constant"), 3, S) == "fail"
    assert state_01(main_row, dict(status="too_few_samples"), 3, S) == "untestable"
    assert state_01(main_row, dict(status="too_few_pairs"), 3, S) == "untestable"
