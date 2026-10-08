import numpy as np
import pandas as pd
import pytest

import splice_assay as sa
from splice_assay.clinical import auto_clinical, clean_sex, clean_stage, detect
from splice_assay.cli import main
from splice_assay.probe import combine, pick_cohorts, probe

FAST = sa.Settings(formats=("png",), dpi=60)


def test_detect_and_clean_clinical_columns():
    cols = ["patient_id", "gender", "ajcc_pathologic_stage", "age_at_diagnosis", "clinical_stage"]
    assert detect(cols) == {"age": "age_at_diagnosis", "sex": "gender", "stage": "ajcc_pathologic_stage"}
    assert detect(["stage", "Age"]) == {"age": "Age", "stage": "stage"}           # 'stage' is not an age column
    assert [clean_stage(v) for v in ["Stage IIA", "IIIc", "Stage IV", "I", "Stage 0", "missing", "[Not Available]",
                                     "Stage X", np.nan]][:5] == ["II", "III", "IV", "I", "0"]
    assert all(pd.isna(clean_stage(v)) for v in ["missing", "[Not Available]", "Stage X", np.nan])
    # stage numbers: integers, floats (a numeric column with gaps), text with a substage
    assert [clean_stage(v) for v in [1, 2.0, np.int64(3), "4", "2B", "Stage 3a", "1.0", 0]] == \
        ["I", "II", "III", "IV", "II", "III", "I", "0"]
    assert all(pd.isna(clean_stage(v)) for v in [5, 1.5, "-1", True, "IS"])
    assert [clean_sex(v) for v in ["FEMALE", "m", "Male"]] == ["female", "male", "male"]
    assert pd.isna(clean_sex("not reported"))


def test_auto_model_from_a_tcga_style_table(tables):
    cl = tables["clinical"].rename(columns={"sex": "gender", "stage": "ajcc_pathologic_stage"})
    cl["ajcc_pathologic_stage"] = "Stage " + cl.ajcc_pathologic_stage.astype(str) + "A"
    cl.loc[cl.index[:4], "ajcc_pathologic_stage"] = "[Not Available]"
    ds = sa.Dataset.from_tables(**dict(tables, clinical=cl))
    ds2, model, found = auto_clinical(ds)
    assert found == {"age": "age", "sex": "gender", "stage": "ajcc_pathologic_stage"}
    assert model.covariates == ("age", "sex", "stage") and dict(model.baseline) == {"sex": "female", "stage": "I"}
    assert set(ds2.clinical.stage.dropna()) == {"I", "II", "III"} and ds2.clinical.stage.isna().sum() == 4
    assert ds.clinical is not ds2.clinical and "stage" not in ds.clinical.columns       # the input is not changed


def test_a_mostly_missing_covariate_is_left_out_of_that_cohort(tables):
    cl = tables["clinical"].copy()
    coh1 = set(tables["samples"].patient_id[tables["samples"].cohort.eq("COH1")])
    cl.loc[cl.patient_id.isin(coh1), "stage"] = np.nan                          # no stage at all in COH1
    ds = sa.Dataset.from_tables(**dict(tables, clinical=cl))
    m = sa.CoxModel(covariates=("age", "stage"), baseline={"stage": "I"})
    sv = sa.analyze(ds, events=["SYN1:SE:1"], endpoints=["OS"], model=m).survival.set_index("cohort")
    assert sv.at["COH1", "cox_status"] == "tested" and "stage left out (0% recorded)" in sv.at["COH1", "cox_notes"]
    assert sv.at["COH1", "cox_model"] == "PSI + host expression + age"
    assert sv.at["COH3", "cox_model"] == "PSI + host expression + age + stage"


def test_pick_cohorts_prefers_the_adjusted_p(ds):
    base = sa.analyze(ds, events=["SYN1:SE:1"], endpoints=["OS"])
    adj = sa.analyze(ds, events=["SYN1:SE:1"], endpoints=["OS"], model=sa.CoxModel(covariates=("age",)))
    cells = combine(base, adj, "OS")
    got = pick_cohorts(cells, "SYN1:SE:1", 2)
    want = cells[cells.adj_cox_status.eq("tested")].sort_values("adj_cox_p").cohort.head(2).tolist()
    assert got == want


def test_probe_outputs(ds, gtf_path, tmp_path):
    logged = []
    res = probe(ds, genes=["SYN1"], settings=FAST, gtf=gtf_path, out_dir=tmp_path, top=2, log=logged.append)
    drawn = [x.strip().split(":")[0] for x in logged if x.startswith("  ")]  # the PDF's pages after the overview
    assert drawn == ["gene map", "page 1", "page 2", "page 3", "expression page"]   # the expression page last
    assert res.paths["gene_maps"] == tmp_path / "gene_map_SYN1.png" and res.paths["gene_maps"].exists()
    ev = res.events
    assert list(ev["rank"]) == [1, 2, 3] and set(ev.event_id) == {"SYN1:SE:1", "SYN1:A3SS:1", "SYN1:RI:1"}
    assert ev.adj_cox_p05.is_monotonic_decreasing or ev.iloc[0].adj_cox_p05 >= ev.iloc[1].adj_cox_p05
    for k in ("report", "events", "cells", "overview", "pdf"):
        assert res.paths[k].exists()
    assert len(list((tmp_path / "pages").glob("*.png"))) == 4                # and the expression page
    assert res.paths["expression_pages"] == tmp_path / "pages" / "SYN1_expression.png"
    report = res.paths["report"].read_text()
    assert "chance alone gives about" in report and "Benjamini" in report and "(pages/001_" in report
    cells = pd.read_csv(res.paths["cells"])
    assert {"cox_q", "adj_cox_q", "adj_hr_per_iqr", "group_hit"} <= set(cells.columns)
    assert cells.event_id.nunique() == 3 and cells.cohort.nunique() == 4


def test_probe_every_event_without_adjustment(ds, tmp_path):
    res = probe(ds, settings=FAST, adjusted=None, out_dir=tmp_path, max_pages=1, log=lambda *_: None)
    assert len(res.events) == 8 and "adj_cox_p" not in res.cells.columns          # the HIT index is left out
    assert "2 HIT-index events were left out" in res.paths["report"].read_text()
    names = sorted(p.name for p in (tmp_path / "pages").glob("*.png"))      # one event page, then its gene's expression
    assert len(names) == 2 and names[1] == f"{res.events.gene.iloc[0]}_expression.png"


def test_cli_probe_and_default_panel(tmp_path, capsys):
    out = tmp_path / "ex"
    assert main(["example", str(out)]) == 0
    data = out / "data"
    assert main(["probe", str(data), "--gene", "SYN2", "--out", str(tmp_path / "pr"), "--settings",
                 str(_fast(tmp_path))]) == 0
    text = capsys.readouterr().out
    assert "1. MXE:1 (SYN2)" in text and (tmp_path / "pr" / "report.md").exists()
    # panel with nothing but the event: cohorts picked, OS, model rows found in the clinical table
    assert main(["panel", str(data), "--event", "SYN1:SE:1", "--out", str(tmp_path / "pn"), "--settings",
                 str(_fast(tmp_path))]) == 0
    text = capsys.readouterr().out
    assert "cohorts shown (most promising of 4)" in text and "age + sex + stage" in text
    csv = sorted((tmp_path / "pn").glob("*.csv"))
    assert [p.name.startswith("SYN1_expression_") for p in csv] == [False, True]     # the page, its gene's expression
    fig = pd.read_csv(csv[0], low_memory=False)
    assert fig[fig.panel.eq("cox_detail")].tag.nunique() == 3


def _fast(tmp_path):
    p = tmp_path / "fast.json"
    p.write_text('{"formats": ["png"], "dpi": 60}')
    return p


def test_top_all_shows_every_tested_cohort(ds, tmp_path):
    from splice_assay.probe import ALL, parse_top
    assert parse_top("all") == ALL and parse_top("5") == 5
    with pytest.raises(Exception):
        parse_top("0")
    res = probe(ds, events=["SYN1:SE:1"], settings=FAST, adjusted=None, out_dir=tmp_path, top=ALL,
                log=lambda *_: None)
    tested = set(res.cells[res.cells.cox_status.eq("tested") | res.cells.km_status.eq("tested")].cohort)
    assert set(res.pages[0][2]) == tested and len(tested) == 4
    assert "every cohort with a test" in res.paths["report"].read_text()
    assert sorted(p.name for p in (tmp_path / "pages").glob("*.png")) == ["001_SYN1_SE_1.png",
                                                                          "SYN1_expression.png"]


def test_many_cohorts_get_a_short_title_and_stem(tables):
    s = tables["samples"].copy()
    name = {c: f"Cohort with a long descriptive name {i}" for i, c in enumerate(sorted(s.cohort.unique()))}
    s["cohort"] = s.cohort.map(name)
    first = sorted(name.values())[0]                             # a fifth cohort: half of the first one's patients
    pts = sorted(s.patient_id[s.cohort.eq(first)].unique())
    s.loc[s.patient_id.isin(pts[: len(pts) // 2]), "cohort"] = "Cohort with a long descriptive name 5"
    ds5 = sa.Dataset.from_tables(**dict(tables, samples=s))
    p = sa.event_panel(ds5, "SYN1:SE:1", ds5.cohorts, "OS", settings=FAST)
    assert "_5cohorts_" in p.stem and len(p.stem) < 60
    assert "SE:1 · 5 cohorts · OS" in [t.get_text() for t in p.figure.texts]


def test_probe_when_no_cox_model_can_be_fitted(ds, tmp_path):
    s = sa.Settings(formats=("png",), dpi=60, cox_min_events=10_000)      # e.g. a small subset: too few deaths
    res = probe(ds, genes=["SYN1"], settings=s, out_dir=tmp_path, top=1, max_pages=1, log=lambda *_: None)
    assert res.cells.cox_status.ne("tested").all() and "cox_p" in res.cells.columns
    assert (res.events.cox_p05 == 0).all() and res.paths["report"].exists()
    assert res.events.measurable.any()                                    # KM still ran


def test_relaxed_gates_are_reported(ds, tmp_path):
    s = sa.Settings(formats=("png",), dpi=60, min_pairs=7, min_group=7, cox_min_events=5)
    assert s.changed() == {"min_pairs": (7, 10), "min_group": (7, 10), "cox_min_events": (5, 10)}
    assert sa.Settings(dpi=60, formats=("png",)).changed() == {}                # drawing settings do not count
    res = probe(ds, genes=["SYN1"], settings=s, out_dir=tmp_path, top=1, max_pages=1, log=lambda *_: None)
    report = res.paths["report"].read_text()
    assert "Settings changed from the defaults:** min_pairs 7 (default 10), min_group 7 (default 10), " \
           "cox_min_events 5 (default 10)" in report
    p = sa.event_panel(ds, "SYN1:SE:1", ["COH1"], "OS", settings=s)
    assert "Settings changed from the defaults" in " ".join(t.get_text() for t in p.figure.texts)


def test_flags_are_reported(ds, tmp_path):
    s = sa.Settings(formats=("png",), dpi=60, narrow_psi_below=1.0)            # every PSI range counts as narrow
    res = probe(ds, genes=["SYN1"], settings=s, out_dir=tmp_path, top=1, max_pages=1, log=lambda *_: None)
    t = res.cells[res.cells.cox_status.eq("tested")]
    assert len(t) and t.psi_narrow.astype(bool).all() and "narrow PSI range (SD" in t.cox_notes.iloc[0]
    a = t[t.adj_cox_status.eq("tested")]
    assert len(a) and a.adj_psi_narrow.astype(bool).all() and a.adj_cox_events_per_term.notna().all()
    rep = (tmp_path / "report.md").read_text()
    assert "## Notes on the survival tests" in rep and f"{len(t)} of {len(t)} base-model fits" in rep


def test_page_parts_are_balanced():
    from splice_assay.plot.panel_common import page_parts
    assert [len(p) for p in page_parts(range(27), 6)] == [6, 6, 5, 5, 5]
    assert page_parts(["a", "b"], 6) == [["a", "b"]] and page_parts(list("abcd"), 0) == [list("abcd")]
    assert sum(page_parts(list(range(13)), 4), []) == list(range(13))


def test_many_cohorts_are_split_over_pages(ds, tmp_path):
    s = sa.Settings(formats=("png",), dpi=60, cohorts_per_page=2)
    res = probe(ds, genes=["SYN1"], settings=s, out_dir=tmp_path, top=4, max_pages=1, log=lambda *_: None)
    assert [st[-3:] for _, st, _ in res.pages] == ["_p1", "_p2"] and len(sum([c for *_, c in res.pages], [])) == 4
    assert res.events.page.iloc[0].endswith("_p1.png")
    p = sa.event_panel(ds, "SYN1:SE:1", ["COH1"], "OS", settings=s, part=(1, 2))
    assert any("page 1 of 2" in t.get_text() for t in p.figure.texts) and p.stem.endswith("_p1of2")


def test_probe_several_endpoints_on_the_command_line(tables, tmp_path, capsys):
    d = tmp_path / "data"
    d.mkdir()
    for name in ("samples", "psi", "events", "survival", "clinical"):
        tables[name].to_csv(d / f"{name}.csv", index=False)
    assert main(["probe", str(d), "--gene", "SYN1", "--endpoint", "OS", "--endpoint", "DSS", "--max-pages", "1",
                 "--settings", str(_fast(tmp_path)), "--out", str(tmp_path / "pr")]) == 0
    assert (tmp_path / "pr" / "OS" / "report.md").exists() and (tmp_path / "pr" / "DSS" / "report.md").exists()
    assert "== DSS" in capsys.readouterr().out


def test_the_default_adjusted_model_builds_on_the_base_model(ds, tmp_path):
    """Without host expression in the base model (--no-expression), the adjusted model has none either, and the
    report names both models as fitted."""
    res = probe(ds, genes=["SYN1"], settings=FAST, model=sa.CoxModel(expression=False), out_dir=tmp_path, max_pages=0,
                gex=False, log=lambda *_: None)
    assert set(res.cells.cox_model.dropna()) == {"PSI"}
    adj = set(res.cells.adj_cox_model.dropna())
    assert adj and all(m.startswith("PSI + age") and "host expression" not in m for m in adj)
    report = res.paths["report"].read_text()
    assert "(the forest on every page): Cox PSI." in report
    assert "(the model rows on every page and the ranking): Cox PSI + age + sex + stage" in report


def test_panel_model_rows_build_on_the_page_model(tmp_path, capsys):
    """panel's default model rows: the page's model plus the age, sex and stage found, keeping --no-expression and
    --strata; --detail keeps them too."""
    from splice_assay import example
    example.write(tmp_path / "ex")
    data = str(tmp_path / "ex" / "data")
    fast = tmp_path / "fast.json"
    fast.write_text('{"formats": ["png"], "dpi": 60}')
    base = ["panel", data, "--event", "SYN1:SE:1", "--cohort", "COH1", "--endpoint", "OS", "--no-gex", "--settings",
            str(fast)]
    for extra, want in ((["--no-expression"], "model rows: Cox PSI + age + sex + stage; found"),
                        (["--strata", "stage"], "model rows: Cox PSI + host expression + age + sex; strata: stage"),
                        (["--detail"], "model rows: Cox PSI + host expression + age + sex + stage; found")):
        assert main(base + extra + ["--out", str(tmp_path / extra[0].strip("-"))]) == 0
        assert want in capsys.readouterr().out, extra
    rows = pd.read_csv(next((tmp_path / "no-expression").glob("*.csv")), low_memory=False)
    terms = set(rows[rows.panel.eq("cox_detail")].term.dropna())
    assert {"PSI", "age", "sex", "stage"} <= terms and "host expression" not in terms


def test_a_variable_the_base_model_has_is_not_added_again(tables):
    """TCGA-style names: the base model's gender, age_at_diagnosis or stage stratum is not added again as the cleaned
    sex, age or stage (that would duplicate the variable and fail every fit)."""
    from splice_assay.clinical import adjusted
    cl = tables["clinical"].rename(columns={"sex": "gender", "stage": "ajcc_pathologic_stage",
                                            "age": "age_at_diagnosis"})
    ds2, auto, found = auto_clinical(sa.Dataset.from_tables(**dict(tables, clinical=cl)))
    for base, want in ((sa.CoxModel(covariates=("gender",)), "PSI + host expression + gender + age + stage"),
                       (sa.CoxModel(strata=("ajcc_pathologic_stage",)),
                        "PSI + host expression + age + sex; strata: ajcc_pathologic_stage"),
                       (sa.CoxModel(covariates=("age_at_diagnosis",)),
                        "PSI + host expression + age_at_diagnosis + sex + stage")):
        m = adjusted(base, auto, found)
        assert m.describe(True) == want
        sv = sa.analyze(ds2, events=["SYN1:SE:1"], endpoints=["OS"], model=m).survival
        assert sv.cox_status.eq("tested").any() and not sv.cox_status.eq("failed").any()


def test_numeric_stage_codes_are_kept(tables):
    """A stage column coded 1-3 is cleaned to I-III, so the adjusted model keeps stage."""
    cl = tables["clinical"].copy()
    cl["stage"] = cl.stage.map({"I": 1, "II": 2, "III": 3}).astype(float)
    cl.loc[cl.index[0], "stage"] = np.nan                                    # a gap makes the column float
    ds2, model, found = auto_clinical(sa.Dataset.from_tables(**dict(tables, clinical=cl)))
    assert "stage" in model.covariates and set(ds2.clinical.stage.dropna()) == {"I", "II", "III"}
    want = auto_clinical(sa.Dataset.from_tables(**tables))[0].clinical.stage
    assert ds2.clinical.stage.iloc[1:].equals(want.iloc[1:])


def test_stage_numbers_are_read_only_from_overall_stage_columns(tables):
    """T, N, M and summary stages are not taken for the overall stage, and numbers are read only from a column named
    as an overall stage (a SEER summary stage 1-4 is not stage I-IV)."""
    assert detect(["pathologic_t_stage", "overall_stage"]) == {"stage": "overall_stage"}
    assert "stage" not in detect(["seer_summary_stage", "clinical_n_stage", "t_stage"])
    assert detect(["stage_event_pathologic_stage"]) == {"stage": "stage_event_pathologic_stage"}
    assert clean_stage("02") == "II" and pd.isna(clean_stage("02", numbers=False))
    assert clean_stage("Stage IIB", numbers=False) == "II" and pd.isna(clean_stage(3, numbers=False))
    cl = tables["clinical"].drop(columns="stage").assign(
        stage_event_pathologic_stage=tables["clinical"].stage.map({"I": 1, "II": 2, "III": 3}))
    ds2, model, found = auto_clinical(sa.Dataset.from_tables(**dict(tables, clinical=cl)))
    assert found["stage"] == "stage_event_pathologic_stage" and "stage" not in model.covariates   # numbers: not read


def test_the_base_models_own_columns_are_not_rewritten(tables):
    """A column the base model names (here a raw stage) is left as it is by the automatic adjustment, so the forest
    has the same numbers with or without model rows."""
    cl = tables["clinical"].copy()
    cl["stage"] = "Stage " + cl.stage.astype(str) + "A"
    ds = sa.Dataset.from_tables(**dict(tables, clinical=cl))
    ds2, auto, found = auto_clinical(ds, keep=("stage",))
    assert "stage" not in found and ds2.clinical.stage.equals(ds.clinical.stage) and auto.covariates == ("age", "sex")
    m = sa.CoxModel(covariates=("stage",))
    a, b = (sa.analyze(d, events=["SYN1:SE:1"], endpoints=["OS"], model=m).survival for d in (ds, ds2))
    assert np.array_equal(a.cox_p.to_numpy(float), b.cox_p.to_numpy(float), equal_nan=True)


def test_the_probe_names_the_models_as_fitted_without_expression(tables, tmp_path):
    """No expression for the host gene: the report names both models without it and counts the fits."""
    ex = tables["expression"]
    ds = sa.Dataset.from_tables(**dict(tables, expression=ex[ex.gene.ne("SYN1")]))
    res = probe(ds, genes=["SYN1"], settings=FAST, out_dir=tmp_path, max_pages=0, gex=False, log=lambda *_: None)
    report = res.paths["report"].read_text()
    assert "(the forest on every page): Cox PSI." in report
    assert "(the model rows on every page and the ranking): Cox PSI + age + sex + stage" in report
    assert "**Host expression left out**" in report
    assert set(res.cells.cox_model.dropna()) == {"PSI"}


def test_low_power_reaches_the_ranking(ds, tmp_path):
    s = FAST.replace(cox_low_power_events=10_000)                     # every fit counts as low power
    res = probe(ds, genes=["SYN1"], settings=s, out_dir=tmp_path, max_pages=0, gex=False, log=lambda *_: None)
    cox = res.events[res.events.best_model.isin(["adjusted", "base"])]
    assert len(cox) and cox.best_low_power.all()
    rep = res.paths["report"].read_text()
    assert " ‡ |" in rep and "**Low power**" in rep and "‡ low power (fewer than 10000 Cox events)" in rep
    top = res.events.iloc[0]                                          # every hit is in a low-power fit
    assert top.adj_cox_p05 > 0 and top.adj_cox_p05_low_power == top.adj_cox_p05
    assert f"| {top.adj_cox_p05} ({top.adj_cox_p05_low_power}‡) |" in rep


def test_ranking_counts_low_power_hits_last(ds):
    """Cox hits in fits that are not low power rank the events; low-power hits (‡) only break ties. The best cohort and
    the page cohorts prefer a hit that is not low power."""
    from splice_assay.probe import pick_cohorts, rank_events
    ev = ["SYN1:SE:1", "SYN1:A3SS:1", "SYN1:RI:1", "SYN2:MXE:1"]
    adj = {ev[0]: [0.01, 0.02, 0.5, 0.5], ev[1]: [0.03, 0.04, 0.001, 0.5], ev[2]: [0.04, 0.5, 0.002, 0.003],
           ev[3]: [0.5, 0.5, 0.01, 0.5]}
    low = [False, False, True, True]                                  # cohorts C and D have few events
    cells = pd.DataFrame([dict(event_id=e, cohort=c, cox_status="tested", cox_p=0.5, hr_per_sd=1.2, cox_q=0.9,
                               cox_low_power=lo, adj_cox_status="tested", adj_cox_p=p, adj_hr_per_sd=1.5,
                               adj_cox_low_power=lo, km_status="gate", km_p=np.nan, logrank_hr=np.nan,
                               group_hit=False, within_patient_support=False, survival_hit=False, psi_sd=0.1)
                          for e in ev for c, p, lo in zip("ABCD", adj[e], low)])
    r = rank_events(cells, ds, 0.05).set_index("event_id")
    assert list(r.sort_values("rank").index) == [ev[1], ev[0], ev[2], ev[3]]   # 2+1‡, 2, 1+2‡, 0+1‡
    assert r.loc[ev[1], ["adj_cox_p05", "adj_cox_p05_low_power", "best_cohort", "best_low_power"]].tolist() == \
        [3, 1, "A", False]                                            # A (p 0.03), not C (p 0.001, low power)
    assert r.loc[ev[3], "best_cohort"] == "C" and r.loc[ev[3], "best_low_power"]   # only a low-power hit
    assert pick_cohorts(cells, ev[1], 3) == ["A", "B", "C"] and pick_cohorts(cells, ev[3], 2) == ["C", "A"]


def _graded(cox_p=0.5, hr=1.5, km_p=0.5, lhr=1.5, d=0.2, hit=False, ci=(1.1, 2.0), paired=True, cohort="X"):
    """A probe cell with a base Cox fit, a KM test and a within-patient change (no adjusted model)."""
    return dict(event_id="e", cohort=cohort, endpoint="OS", cox_status="tested", cox_p=cox_p, hr_per_sd=hr,
                ci_low_sd=ci[0], ci_high_sd=ci[1], cox_low_power=False, km_status="tested", km_p=km_p,
                logrank_hr=lhr, paired_status="tested" if paired else "too_few_pairs", paired_delta_median=d,
                paired_hit=hit, unpaired_status="too_few_samples")


def test_evidence_grades():
    """Cox (as the overview shows it), KM and the group change: their directions and significance give the grade."""
    from splice_assay.probe import evidence
    s = sa.Settings()
    g = (lambda **kw: evidence(pd.Series(_graded(**kw)), s))                   # noqa: E731
    assert g(cox_p=0.01, km_p=0.01) == ("A", "Cox ↑  KM ↑  T/N (↑)")
    assert g(cox_p=0.01, km_p=0.01, hit=True) == ("A+", "Cox ↑  KM ↑  T/N ↑")
    assert g(cox_p=0.01, km_p=0.01, hr=0.6, lhr=0.7, d=-0.2, hit=True)[0] == "A+"     # all down: one story too
    assert g(cox_p=0.01)[0] == "B" and g(km_p=0.01)[0] == "C" and g(km_p=0.01, hit=True)[0] == "C+"
    assert g(cox_p=0.01, d=-0.2) == ("D", "Cox ↑  KM (↑)  T/N (↓)")          # the tumour shift goes with lower hazard
    assert g(cox_p=0.01, d=-0.2, hit=True)[0] == "E"                         # significantly so
    assert g(cox_p=0.01, km_p=0.01, lhr=0.7)[0] == "E"                       # KM against Cox
    assert g(cox_p=0.01, ci=(0.3, 9.0)) == ("", "Cox (↑)  KM (↑)  T/N (↑)")   # an imprecise fit does not count
    assert g(cox_p=0.01, paired=False) == ("B", "Cox ↑  KM (↑)  T/N –")      # no group test: it cannot disagree
    assert g(hit=True)[0] == "" and g()[0] == ""                             # no survival signal
    both = dict(_graded(cox_p=0.01, d=0.02), unpaired_status="tested", unpaired_delta_median=-0.15, unpaired_hit=True)
    assert evidence(pd.Series(both), s) == ("E", "Cox ↑  KM (↑)  T/N ↓")     # the all-samples hit decides
    assert g(cox_p=0.01, d=1e-17) == ("B", "Cox ↑  KM (↑)  T/N =")           # rounding noise has no direction
    assert g(cox_p=0.01, d=-0.03) == ("B", "Cox ↑  KM (↑)  T/N =")           # a tiny change neither (under 0.05)
    assert g(cox_p=0.01, lhr=0.95) == ("B", "Cox ↑  KM =  T/N (↑)")          # nor an HR within 1.1-fold of 1
    assert g(km_p=0.01, hr=0.97) == ("C", "Cox =  KM ↑  T/N (↑)")
    assert evidence(pd.Series(_graded(cox_p=0.01, lhr=0.95)), s.replace(evidence_hr_band=1))[0] == "D"
    hit_index = pd.Series(_graded(cox_p=0.01, d=-0.08))                     # the HIT index: twice PSI's scale
    assert evidence(hit_index, s, 2.0)[0] == "B" and evidence(hit_index, s)[0] == "D"
    with pytest.raises(ValueError, match="evidence_hr_band"):
        sa.Settings(evidence_hr_band=0.9)
    with pytest.raises(ValueError, match="evidence_min_delta"):
        sa.Settings(evidence_min_delta=float("nan"))
    for lhr, o, e, mark in [(0.0, 0, 3.2, "↓"), (np.nan, 9, 5.8, "↑")]:     # an arm without events: no HR
        one_arm = pd.Series(dict(_graded(km_p=0.001, lhr=lhr), o_high=o, e_high=e))
        assert evidence(one_arm, s) == ("C" if mark == "↑" else "D", f"Cox (↑)  KM {mark}  T/N (↑)")
    assert g(cox_p=0.01, km_p=0.01, hit="True")[0] == "A+" and g(cox_p=float("nan"), km_p=0.01)[0] == "C"
    adj = dict(_graded(cox_p=0.5, hr=1.5), adj_cox_status="tested", adj_cox_p=0.01, adj_hr_per_sd=0.6,
               adj_ci_low_sd=0.4, adj_ci_high_sd=0.9, adj_cox_low_power=False)
    assert evidence(pd.Series(adj), s) == ("D", "Cox ↓  KM (↑)  T/N (↑)")    # the adjusted fit is the one shown
    iqr = dict(_graded(cox_p=0.01), hr_per_iqr=np.nan, ci_low_iqr=np.nan, ci_high_iqr=np.nan)
    assert evidence(pd.Series(iqr), s.replace(psi_hr_unit="iqr"))[1].startswith("Cox ↑")   # an IQR of 0
    wide = dict(iqr, ci_low_sd=0.1, ci_high_sd=22.0)                         # ... and the per-SD CI is 220-fold
    assert evidence(pd.Series(wide), s.replace(psi_hr_unit="iqr")) == evidence(pd.Series(wide), s) == \
        ("", "Cox (↑)  KM (↑)  T/N (↑)")


def test_page_cohorts_are_the_graded_ones():
    """Every cohort graded A+ to C, best grade first (then by p), at most cohorts_per_page; else the strongest."""
    from splice_assay.probe import add_evidence, page_cohorts
    s = sa.Settings()
    cells = add_evidence(pd.DataFrame([
        _graded(cohort="P", cox_p=0.001),                                     # B, the smallest p
        _graded(cohort="V", cox_p=0.01),                                      # B: after P
        _graded(cohort="Q", cox_p=0.02, km_p=0.01),                           # A
        _graded(cohort="R", km_p=0.03),                                       # C
        _graded(cohort="S", cox_p=0.04, km_p=0.04, hit=True),                 # A+
        _graded(cohort="T", cox_p=0.002, d=-0.3),                             # D: not on the page
        _graded(cohort="U", cox_p=0.3)]), s)                                  # no grade
    assert list(cells.columns[:5]) == ["event_id", "cohort", "endpoint", "evidence", "evidence_lines"]
    assert cells.evidence.tolist() == ["B", "B", "A", "C", "A+", "D", ""]
    assert page_cohorts(cells, "e", s) == (["S", "Q", "P", "V", "R"], "")
    assert page_cohorts(cells, "e", s.replace(cohorts_per_page=2)) == (["S", "Q"], "2 of 5 cohorts graded A–C")
    none = add_evidence(pd.DataFrame([_graded(cohort="P", cox_p=0.3), _graded(cohort="Q", cox_p=0.002, d=-0.3)]), s)
    assert page_cohorts(none, "e", s) == (["Q"], "no cohort graded A–C")    # the strongest by p


def test_probe_pages_show_the_graded_cohorts(ds, tmp_path):
    """By default each page shows its event's cohorts graded A+ to C; --top N keeps the N most promising by p. The
    grades are in cells.csv, counted in events.csv and explained in the report."""
    logged = []
    res = probe(ds, genes=["SYN1"], settings=FAST, out_dir=tmp_path / "a", max_pages=3, log=logged.append)
    c = pd.read_csv(res.paths["cells"], keep_default_na=False)
    assert {"evidence", "evidence_lines"} <= set(c.columns)
    graded = c[c.evidence.isin(["A+", "A", "B+", "B", "C+", "C"])]
    for e in res.events.itertuples():
        mine = graded[graded.event_id.eq(e.event_id)]
        assert e.evidence_a_c == len(mine)
        line = next(x for x in logged if x.startswith("  page") and f": {e.event_id} (" in x)
        shown = line.split("(", 1)[1].rstrip(")").split(", ")
        assert set(shown) == set(mine.cohort) if len(mine) else len(shown) == 1
    report = res.paths["report"].read_text()
    assert "**Evidence grades**" in report and "| Graded A–C |" in report
    assert "every cohort graded A+ to C" in report
    logged.clear()
    probe(ds, genes=["SYN1"], settings=FAST, out_dir=tmp_path / "b", max_pages=1, top=2, log=logged.append)
    assert len(next(x for x in logged if x.startswith("  page")).split("(", 1)[1].split(", ")) == 2
    logged.clear()                                                          # named cohorts: still graded
    named = probe(ds, genes=["SYN1"], cohorts=["COH3", "COH4"], settings=FAST, out_dir=tmp_path / "c", max_pages=3,
                  log=logged.append)
    for e in named.events.itertuples():
        line = next(x for x in logged if x.startswith("  page") and f": {e.event_id} (" in x)
        ok = named.cells[named.cells.event_id.eq(e.event_id) & named.cells.evidence.isin(["A+", "A", "B+", "B", "C+",
                                                                                           "C"])]
        shown = line.split("(", 1)[1].rstrip(")").split(", ")
        assert set(shown) == set(ok.cohort) if len(ok) else len(shown) == 1


def test_panel_rows_show_the_page_model_when_nothing_is_left_to_add(tmp_path, capsys):
    """--covariate age, sex and stage already on the page: the model rows show that model (nothing is found to add)."""
    from splice_assay import example
    example.write(tmp_path / "ex")
    fast = tmp_path / "fast.json"
    fast.write_text('{"formats": ["png"], "dpi": 60}')
    assert main(["panel", str(tmp_path / "ex" / "data"), "--event", "SYN1:SE:1", "--cohort", "COH1", "--endpoint", "OS",
                 "--no-gex", "--covariate", "age", "--covariate", "sex", "--covariate", "stage", "--settings",
                 str(fast),
                 "--out", str(tmp_path / "fig")]) == 0
    rows = pd.read_csv(next((tmp_path / "fig").glob("*.csv")), low_memory=False)
    assert {"PSI", "age", "sex", "stage"} <= set(rows[rows.panel.eq("cox_detail")].term.dropna())
    capsys.readouterr()


def test_the_report_says_what_changed_settings_mean(ds, tmp_path):
    """Relaxed gates are named as such; a ridge or another HR unit is not a relaxed gate. The overview's caption names
    a penalty only when the fits it shows have one."""
    n = iter(range(100))
    run = (lambda s, **kw: probe(ds, genes=["SYN1"], settings=s, out_dir=tmp_path / str(next(n)), max_pages=0,
                                 gex=False, log=lambda *_: None, **kw))
    rep = (lambda s: run(s).paths["report"].read_text())
    assert "relaxed gates" not in rep(FAST.replace(cox_ridge="clinical", psi_hr_unit="iqr"))
    assert "relaxed gates" in rep(FAST.replace(low_psi_variance_sd=0.0))
    assert "relaxed gates" in rep(FAST.replace(min_pairs=5))
    from splice_assay.probe import overview
    caption = (lambda r, s: next(x.get_text() for x in overview(r.cells, r.events, "OS", s).texts
                                 if x.get_text().startswith("colour:")))
    s = FAST.replace(cox_ridge="clinical")
    assert "(adjusted model where fitted; ridge λ 1 on clinical terms)" in caption(run(s), s)
    r = run(FAST)                                        # the counts at the right of each row
    fig = overview(r.cells, r.events, "OS", FAST)
    def rows(f):
        by = {}
        for t in f.axes[0].texts:
            by.setdefault(t.get_position()[1], []).append(t.get_text())
        return by
    by_row = rows(fig)
    assert by_row[-0.25] == ["p < 0.05", "chance", "HR"]
    for i, e in enumerate(r.events[r.events.measurable].itertuples()):        # each row: k/n, chance, directions
        c = r.cells[r.cells.event_id.eq(e.event_id)]
        adj = c.adj_cox_status.eq("tested")
        p = c.adj_cox_p.where(adj, c.cox_p.where(c.cox_status.eq("tested")))
        hr = c.adj_hr_per_sd.where(adj, c.hr_per_sd)[p < 0.05]
        k, m = int((p < 0.05).sum()), int(p.notna().sum())
        assert by_row[i + 0.5] == [f"{k}/{m}", f"{0.05 * m:.1f}", f"{(hr > 1).sum()}↑ {(hr < 1).sum()}↓" if k else ""]
    strict = overview(r.cells, r.events, "OS", FAST.replace(alpha=0.001))
    assert "(0.1% of them)" in " ".join(t.get_text() for t in strict.texts)
    assert {rows(strict)[i + 0.5][1] for i in range(int(r.events.measurable.sum()))} == {"0.004"}   # not "0.0"
    assert "ridge" not in caption(run(s, adjusted=None), s)                  # no clinical term was penalized


def test_the_report_says_which_terms_a_ridge_reached(ds, tmp_path):
    """One line, worded for the scope and the terms the penalty reached; none when it reached no fitted model."""
    n = iter(range(100))

    def line(s, **kw):
        r = probe(ds, genes=["SYN1"], settings=s, out_dir=tmp_path / str(next(n)), max_pages=0, gex=False,
                  log=lambda *_: None, **kw)
        found = [x for x in r.paths["report"].read_text().splitlines() if x.startswith("- **Penalty:**")]
        assert len(found) <= 1
        return found[0] if found else ""
    assert line(FAST.replace(cox_ridge="clinical")).startswith(
        "- **Penalty:** ridge λ 1 on the clinical terms of the adjusted model: they are shrunk jointly toward HR 1, "
        "so they adjust PSI only partly: its HR stays closer to the HR without them")
    assert line(FAST.replace(cox_ridge="molecular")).startswith(
        "- **Penalty:** ridge λ 1 on the PSI and host-expression terms of the base and adjusted models: they are "
        "shrunk jointly toward HR 1, though a single HR (PSI's too) can move away from 1")
    every = line(FAST.replace(cox_ridge="all", cox_ridge_penalty=2.0))
    assert every.startswith("- **Penalty:** ridge λ 2 on all terms of the base and adjusted models: they are shrunk "
                            "jointly") and "; the clinical terms adjust PSI only partly" in every
    # clinical terms of the base model alone count too; PSI alone is one term
    assert "; the clinical terms adjust PSI only partly" in line(FAST.replace(cox_ridge="all"), adjusted=None,
                                                                 model=sa.CoxModel(covariates=("age",)))
    assert line(FAST.replace(cox_ridge="all"), adjusted=None, model=sa.CoxModel(expression=False)).startswith(
        "- **Penalty:** ridge λ 1 on the PSI term of the base model: its HR is shrunk toward 1.")
    assert line(FAST.replace(cox_ridge="clinical"), adjusted=None) == ""      # nothing clinical to penalize
    assert line(FAST.replace(cox_ridge="molecular", cox_min_n=10_000)) == ""  # nothing fitted


def test_the_report_names_the_fits_the_events_per_term_minimum_stopped(ds, tmp_path):
    rep = probe(ds, genes=["SYN1"], settings=FAST.replace(cox_min_events_per_term=1000), out_dir=tmp_path,
                max_pages=0, gex=False, log=lambda *_: None).paths["report"].read_text()
    assert "**Not fitted: too few events per term** (fewer than 1000 events per estimated term" in rep
    assert "(except the fits the events-per-term minimum stopped)" in rep and "**Penalty:**" not in rep


def test_report_notes_name_every_flagged_fit(ds):
    """The low-power note names fits with p < alpha in the base or the adjusted model (the ranking uses the adjusted
    one), and an object column with gaps (as in real data) raises no pandas warning; the host-expression line names
    every cohort with p < alpha, ‡ marking low power."""
    from splice_assay.probe import _flag_lines, _gex_cox_line
    cells = pd.DataFrame({
        "event_id": ["SYN1:SE:1", "SYN1:SE:1", "SYN1:A3SS:1"], "cohort": ["COH1", "COH2", "COH1"],
        "cox_status": ["tested"] * 3, "cox_p": [0.01, 0.2, 0.3], "cox_events": [12, 12, 50],
        "adj_cox_status": ["tested"] * 3, "adj_cox_p": [0.02, 0.01, 0.01], "adj_cox_events": [12, 12, 50],
        "km_status": ["gate"] * 3, "psi_narrow": pd.Series([True, np.nan, False], dtype=object)})
    lines = _flag_lines(cells, ds, FAST, sa.CoxModel(covariates=("age",)))
    low = next(x for x in lines if x.startswith("- **Low power**"))
    assert "with p < 0.05 (base or adjusted model): SE:1 COH1, SE:1 COH2." in low      # COH2: adjusted p only
    assert any(x.startswith("- **Narrow PSI range**") and "1 of 3 base-model fits; with p < 0.05: SE:1 COH1." in x
               for x in lines)
    gex = pd.DataFrame({"gene": "SYN1", "cohort": [f"C{i}" for i in range(8)], "cox_status": "tested",
                        "cox_p": np.linspace(0.001, 0.04, 8), "hr_per_sd": 1.5, "cox_low_power": [False] * 7 + [True]})
    line = _gex_cox_line(gex, FAST)
    assert line.startswith("- **Cox:** 8 cohorts tested; expression has p < 0.05 in 8 (C0 HR 1.50, C1 HR 1.50, ")
    assert line.endswith(", C7 HR 1.50 ‡); ‡ fewer than 20 events.")
