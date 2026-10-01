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
    a = sa.analyse(sa.Dataset.from_tables(**tables), events=["SYN1:SE:1"])
    b = sa.analyse(sa.Dataset.from_tables(**shuffled), events=["SYN1:SE:1"])
    cols = ["paired_p", "paired_delta_median", "unpaired_p", "unpaired_delta_median", "paired_hit_status"]
    pd.testing.assert_frame_equal(a.groups[cols], b.groups[cols])
    num = ["cutoff", "km_p", "logrank_hr", "cox_beta", "cox_se", "hr_per_iqr"]
    assert np.allclose(a.survival[num].to_numpy(float), b.survival[num].to_numpy(float), rtol=1e-9, equal_nan=True)
