import warnings

import numpy as np
import pandas as pd
import pytest

from splice_assay import Settings
from splice_assay.stats.survival import CoxModel, fit_cox, survival_cell

S = Settings()


def cohort(seed=0, n=150, missing=0.0):
    rng = np.random.default_rng(seed)
    x = np.round(rng.beta(4, 3, n), 3)
    x[rng.uniform(size=n) < missing] = np.nan
    lp = np.nan_to_num(2.0 * (x - 0.57))
    t = np.maximum(1, np.round(rng.exponential(2000 * np.exp(-lp))))
    c = rng.uniform(300, 4000, n)
    e = (t <= c).astype(int)
    t = np.minimum(t, c)
    return x, np.arange(n), t, e


def test_cox_m0_equals_lifelines():
    from lifelines import CoxPHFitter
    x, pos, t, e = cohort()
    r, _ = survival_cell(x, pos, t, e, None, S)
    assert r["cox_status"] == "tested" and r["cox_model"] == "PSI"
    ref = CoxPHFitter().fit(pd.DataFrame({"time": t, "event": e, "psi10": x / 0.1}), "time", "event").summary
    assert r["cox_beta"] == pytest.approx(ref.loc["psi10", "coef"], rel=1e-9)
    assert r["cox_se"] == pytest.approx(ref.loc["psi10", "se(coef)"], rel=1e-9)
    iqr = np.subtract(*np.percentile(x, [75, 25]))
    assert r["hr_per_iqr"] == pytest.approx(np.exp(r["cox_beta"] * iqr / 0.1))
    assert r["ci_low_iqr"] == pytest.approx(np.exp((r["cox_beta"] - 1.96 * r["cox_se"]) * iqr / 0.1))


def test_cox_m1_adjusts_for_z_scored_expression():
    from lifelines import CoxPHFitter
    x, pos, t, e = cohort(1)
    g = np.random.default_rng(5).normal(5, 1, len(x))
    r, _ = survival_cell(x, pos, t, e, g, S)
    z = (g - g.mean()) / g.std(ddof=1)
    ref = CoxPHFitter().fit(pd.DataFrame({"time": t, "event": e, "psi10": x / 0.1, "gex_z": z}), "time",
                            "event").summary
    assert r["cox_model"] == "PSI + host expression"
    assert r["cox_beta"] == pytest.approx(ref.loc["psi10", "coef"], rel=1e-9)
    assert r["expr_hr_per_sd"] == pytest.approx(ref.loc["gex_z", "exp(coef)"], rel=1e-9)


def test_hr_per_iqr_is_undefined_when_the_iqr_is_zero():
    x, pos, t, e = cohort(6)
    x = np.where(np.arange(len(x)) % 5 == 0, 0.9, 0.3)               # 80% at one value: IQR 0
    r, _ = survival_cell(x, pos, t, e, None, S)
    assert r["cox_status"] == "tested" and r["psi_iqr"] == 0
    assert np.isnan(r["hr_per_iqr"]) and np.isnan(r["ci_low_iqr"]) and np.isfinite(r["hr_per_step"])


def test_km_split_and_cutoff_use_the_whole_survival_cohort():
    x, pos, t, e = cohort(2)
    ep = pos[::2]                                    # only half the patients have this endpoint
    r, _ = survival_cell(x, ep, t[::2], e[::2], None, S)
    assert r["cutoff"] == np.median(x)               # cut over every survival sample, not the endpoint cohort
    high = x[ep] > r["cutoff"]
    assert r["n_high"] == high.sum() and r["events_high"] == e[::2][high].sum()


def test_gates():
    x, pos, t, e = cohort(3, missing=0.6)
    assert survival_cell(x, pos, t, e, None, S)[0]["km_status"] == "coverage_gate"
    x2 = np.full(150, 0.5)
    x2[:30] = 0.501
    r, _ = survival_cell(x2, np.arange(150), t, e, None, S)
    assert r["km_status"] == "low_psi_variance" and r["cox_status"] == "low_psi_variance"
    x3, pos3, t3, e3 = cohort(4, n=40)
    r, _ = survival_cell(x3, pos3, t3, np.zeros(40, int), None, S)
    assert r["km_status"] == "too_few_patients_or_events"
    assert r["cox_status"] == "too_few_patients_or_events"
    r, _ = survival_cell(x3, pos3, t3, e3, np.full(40, 3.0), S)
    assert r["cox_status"] == "constant_expression"


def test_separation_fails_the_fit():
    n = 60
    x = np.r_[np.full(30, 0.1), np.full(30, 0.9)]
    t = np.r_[np.arange(30) + 100.0, np.arange(30) + 1.0]    # every high-PSI patient dies first
    df = pd.DataFrame({"time": t, "event": np.ones(n, int), "psi10": x / 0.1})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = fit_cox(df, "psi10", ph=False)
    assert r["status"] == "failed"
    assert np.isnan(r["beta"]) and r["fit_warning"]


def test_survival_hit_flag(results):
    sv = results.survival.set_index(["event_id", "cohort", "endpoint"])
    assert sv.loc[("SYN1:SE:1", "COH1", "OS")].survival_hit
    assert set(results.survival.cox_model) == {"PSI + host expression"}     # analyse leaves out the HIT index


def test_covariates_and_strata_equal_a_direct_lifelines_fit(ds):
    """PSI + host expression + age (per SD) + stage (vs most common) with strata by sex, refitted by hand."""
    from lifelines import CoxPHFitter
    import splice_assay as sa
    m = CoxModel(covariates=("age", "stage"), strata=("sex",))
    res = sa.analyse(ds, events=["SYN1:SE:1"], cohorts=["COH1"], endpoints=["OS"], model=m)
    row = res.survival.iloc[0]
    assert row.cox_status == "tested" and row.cox_model == "PSI + host expression + age + stage; strata: sex"
    # the same fit by hand
    base = ds.samples[ds.samples.cohort.eq("COH1") & ds.samples.survival_cohort]
    sv = ds.survival[ds.survival.endpoint.eq("OS")].set_index("patient_id")
    b = base[base.patient_id.isin(sv.index)].copy()
    b["psi"] = ds.psi.loc["SYN1:SE:1"].reindex(b.sample_id).to_numpy()
    b["g"] = ds.expression.loc["SYN1"].reindex(b.sample_id).to_numpy()
    b = b.join(ds.clinical[["age", "stage", "sex"]], on="patient_id")
    b = b.dropna(subset=["psi", "g", "age", "stage", "sex"])
    z = lambda v: (v - v.mean()) / v.std(ddof=1)  # noqa: E731
    ref = b.stage.value_counts().idxmax()
    df = pd.DataFrame({"time": sv.time.reindex(b.patient_id).to_numpy(), "event": sv.event.reindex(b.patient_id).to_numpy(),
                       "psi10": b.psi.to_numpy() / 0.1, "gex_z": z(b.g).to_numpy(), "age": z(b.age).to_numpy(),
                       "sex": b.sex.to_numpy()})
    for lev in sorted(set(b.stage) - {ref}):
        df[f"stage_{lev}"] = (b.stage == lev).astype(float).to_numpy()
    fit = CoxPHFitter().fit(df, "time", "event", strata=["sex"]).summary
    assert row.cox_n == len(df)
    assert row.cox_beta == pytest.approx(fit.loc["psi10", "coef"], rel=1e-9)
    assert row.cox_se == pytest.approx(fit.loc["psi10", "se(coef)"], rel=1e-9)
    t = res.cox_terms
    one = lambda term, level: t[t.term.eq(term) & t.level.eq(level)].iloc[0]  # noqa: E731
    assert one("age", "").coef == pytest.approx(fit.loc["age", "coef"], rel=1e-9)
    for lev in sorted(set(b.stage) - {ref}):
        assert one("stage", lev).coef == pytest.approx(fit.loc[f"stage_{lev}", "coef"], rel=1e-9)
        assert one("stage", lev).reference == ref


def test_model_without_expression_and_missing_clinical_table(ds, tables):
    import splice_assay as sa
    r = sa.analyse(ds, events=["SYN1:SE:1"], cohorts=["COH1"], endpoints=["OS"],
                   model=CoxModel(expression=False)).survival.iloc[0]
    assert r.cox_model == "PSI"
    bare = sa.Dataset.from_tables(**{k: v for k, v in tables.items() if k != "clinical"})
    with pytest.raises(sa.InputError, match="no clinical table"):
        sa.analyse(bare, events=["SYN1:SE:1"], model=CoxModel(covariates=("age",)))


def test_rare_levels_merge_into_their_neighbours():
    """Fewer than level_min_patients patients, or no events: stage I joins II and IV joins III, the rarest first."""
    from splice_assay.stats.survival import merge_rare_levels
    v = pd.Series(["I"] * 2 + ["II"] * 50 + ["III"] * 40 + ["IV"] * 3)
    e = np.r_[np.zeros(2), np.ones(93)]
    out, label, notes = merge_rare_levels(v, e, 10)
    assert label == {"I": "I–II", "II": "I–II", "III": "III–IV", "IV": "III–IV"}
    assert notes == ["I merged into II (2 patients, 0 events)", "IV merged into III (3 patients, 3 events)"]
    assert out.value_counts().to_dict() == {"I–II": 52, "III–IV": 43}
    # a rare middle level joins its smaller neighbour; a level without events is rare however large
    assert merge_rare_levels(pd.Series(["I"] * 30 + ["II"] * 4 + ["III"] * 20), np.ones(54), 10)[1]["II"] == "II–III"
    assert merge_rare_levels(pd.Series(["I"] * 15 + ["II"] * 40), np.r_[np.zeros(15), np.ones(40)], 10)[1]["I"] == "I–II"
    # unordered levels join the most common one; 0 merges nothing
    race = pd.Series(["white"] * 100 + ["black"] * 20 + ["asian"] * 5)
    assert merge_rare_levels(race, np.ones(125), 10)[1] == {"asian": "white+asian", "black": "black",
                                                            "white": "white+asian"}
    assert merge_rare_levels(v, e, 0)[2] == []


def test_a_rare_baseline_is_merged_and_the_model_fits(ds):
    """Stage I named as baseline but held by 2 censored patients: it joins II, the reference becomes I–II, and the note
    says so (without the merge the stage I contrast has no estimate)."""
    from dataclasses import replace
    import splice_assay as sa
    sv = ds.survival[ds.survival.endpoint.eq("OS")].set_index("patient_id")
    pts = [p for p in ds.samples[ds.samples.cohort.eq("COH1") & ds.samples.survival_cohort].patient_id if p in sv.index]
    cl = ds.clinical.copy()
    cl["stage"] = cl["stage"].astype(object)
    for k, p in enumerate(pts):
        cl.loc[p, "stage"] = "II" if k % 2 else "III"
    cl.loc[[p for p in pts if sv.at[p, "event"] == 0][:2], "stage"] = "I"
    m = CoxModel(covariates=("stage",), categorical=("stage",), baseline={"stage": "I"})
    res = sa.analyse(replace(ds, clinical=cl), events=["SYN1:SE:1"], cohorts=["COH1"], endpoints=["OS"], model=m)
    row = res.survival.iloc[0]
    assert row.cox_status == "tested"
    assert "stage I merged into II (2 patients, 0 events)" in row.cox_notes
    assert set(res.cox_terms[res.cox_terms.term.eq("stage")].reference) == {"I–II"}


def test_flags_note_few_events_per_term_and_a_narrow_psi_range():
    """Notes only, the fit is unchanged: events per estimated term below cox_events_per_term, and a PSI spread (IQR
    or SD) below narrow_psi_below; 0 turns either off."""
    x, pos, t, e = cohort(3)
    r, _ = survival_cell(x, pos, t, e, None, S)
    assert r["cox_status"] == "tested" and r["cox_events_per_term"] == r["cox_events"]      # one term: PSI
    assert not r["psi_narrow"] and r["cox_notes"] == ""
    s = Settings(cox_events_per_term=1000, narrow_psi_below=0.5)
    r2, _ = survival_cell(x, pos, t, e, None, s)
    assert r2["cox_status"] == "tested" and r2["cox_p"] == r["cox_p"]
    assert f"{r2['cox_events']:.1f} events per term (< 1000)" in r2["cox_notes"]
    assert r2["psi_narrow"] and f"narrow PSI range (IQR {r2['psi_iqr']:.3f} < 0.5)" in r2["cox_notes"]
    r3, _ = survival_cell(x, pos, t, e, None, s.replace(narrow_psi_measure="sd"))
    assert f"narrow PSI range (SD {r3['psi_sd']:.3f} < 0.5)" in r3["cox_notes"]
    r4, _ = survival_cell(x, pos, t, e, None, Settings(cox_events_per_term=0, narrow_psi_below=0))
    assert r4["cox_notes"] == "" and not r4["psi_narrow"]


def crossing(seed=4, n=400):
    """PSI whose effect reverses over follow-up: the high arm has a high hazard for a year, then a low one."""
    rng = np.random.default_rng(seed)
    x = np.round(rng.uniform(0, 1, n), 3)
    u = rng.exponential(size=n)                                     # cumulative hazard at the event
    t = np.where(x > 0.5, np.where(u < 2.0, u / 2.0, 1.0 + (u - 2.0) / 0.1), u / 0.5)
    c = rng.uniform(0.2, 6, n)
    return x, np.arange(n), np.round(np.minimum(t, c) * 365.25) + 1, (t <= c).astype(int)


def test_proportional_hazards_are_noted_not_enforced():
    x, pos, t, e = crossing()
    r, terms = survival_cell(x, pos, t, e, None, S)
    assert r["km_status"] == "tested" and r["cox_status"] == "tested"
    assert r["km_ph_p"] < 0.05 and r["km_notes"].startswith("non-proportional hazards between the arms (p ")
    assert r["ph_p"] < 0.05 and "non-proportional hazards: PSI (p " in r["cox_notes"]
    assert all(np.isfinite(term["ph_p"]) for term in terms)
    quiet, _ = survival_cell(x, pos, t, e, None, Settings(ph_note_below=0))
    assert quiet["km_p"] == r["km_p"] and quiet["cox_p"] == r["cox_p"] and quiet["ph_p"] == r["ph_p"]
    assert "km_notes" not in quiet and "non-proportional" not in quiet["cox_notes"]


def test_stratified_models_get_a_proportional_hazards_test(ds):
    import splice_assay as sa
    m = CoxModel(covariates=("age",), strata=("sex",))
    res = sa.analyse(ds, events=["SYN1:SE:1"], cohorts=["COH1"], endpoints=["OS"], model=m)
    assert np.isfinite(res.survival.iloc[0].ph_p) and res.cox_terms.ph_p.notna().all()


def test_km_split_rules():
    """The KM split is the median by default, else the mean or a value set in the settings; Cox does not change."""
    from splice_assay.stats.survival import expression_cell
    x, pos, t, e = cohort(4)
    med, _ = survival_cell(x, pos, t, e, None, S)
    assert med["cutoff"] == np.median(x) and med["km_split"] == "median"
    mean, _ = survival_cell(x, pos, t, e, None, Settings(km_split="mean"))
    assert np.isclose(mean["cutoff"], x.mean()) and mean["km_split"] == "mean"
    assert mean["n_high"] == (x > x.mean()).sum() and mean["cox_p"] == med["cox_p"]
    fixed, _ = survival_cell(x, pos, t, e, None, Settings(km_split="0.5"))        # a number from the command line
    assert fixed["cutoff"] == 0.5 and fixed["km_split"] == "set" and fixed["n_high"] == (x > 0.5).sum()
    assert survival_cell(x, pos, t, e, None, Settings(km_split=2.0))[0]["km_status"] != "tested"   # one arm empty
    with pytest.raises(ValueError, match='"median", "mean" or a number'):
        Settings(km_split="max")
    g = np.log2(1 + 30 * x)                                    # expression has its own rule
    assert expression_cell(g, pos, t, e, Settings(km_split=0.5))[0]["km_split"] == "median"
    row, _ = expression_cell(g, pos, t, e, Settings(km_split_expression="mean"))
    assert row["km_split"] == "mean" and np.isclose(row["cutoff"], g.mean())


def test_hr_per_sd_setting():
    """The tables hold the HR per IQR and per SD; psi_hr_unit picks the model term shown, and the narrow-range note
    follows it unless narrow_psi_measure is set."""
    x, pos, t, e = cohort(5)
    r, terms = survival_cell(x, pos, t, e, None, S)
    assert np.isclose(r["hr_per_sd"], np.exp(r["cox_beta"] * r["psi_sd"] / 0.10))
    assert np.isclose(r["hr_per_iqr"], np.exp(r["cox_beta"] * r["psi_iqr"] / 0.10))
    assert np.isclose(r["ci_low_sd"], np.exp((r["cox_beta"] - 1.959963984540054 * r["cox_se"]) * r["psi_sd"] / 0.10))
    assert [m["kind"] for m in terms] == ["psi", "psi_iqr"]
    sd = Settings(psi_hr_unit="sd")
    r2, terms2 = survival_cell(x, pos, t, e, None, sd)
    assert [m["kind"] for m in terms2] == ["psi", "psi_sd"] and terms2[1]["unit"].startswith("SD (")
    assert terms2[1]["hr"] == r2["hr_per_sd"] and r2["hr_per_iqr"] == r["hr_per_iqr"] and r2["cox_p"] == r["cox_p"]
    assert S.narrow_measure == "iqr" and sd.narrow_measure == "sd"
    assert Settings(psi_hr_unit="sd", narrow_psi_measure="iqr").narrow_measure == "iqr"
    with pytest.raises(ValueError, match="psi_hr_unit"):
        Settings(psi_hr_unit="range")
