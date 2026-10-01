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
    sv = sa.analyse(ds, events=["SYN1:SE:1"], endpoints=["OS"], model=m).survival.set_index("cohort")
    assert sv.at["COH1", "cox_status"] == "tested" and "stage left out (0% recorded)" in sv.at["COH1", "cox_notes"]
    assert sv.at["COH1", "cox_model"] == "PSI + host expression + age"
    assert sv.at["COH3", "cox_model"] == "PSI + host expression + age + stage"


def test_pick_cohorts_prefers_the_adjusted_p(ds):
    base = sa.analyse(ds, events=["SYN1:SE:1"], endpoints=["OS"])
    adj = sa.analyse(ds, events=["SYN1:SE:1"], endpoints=["OS"], model=sa.CoxModel(covariates=("age",)))
    cells = combine(base, adj, "OS")
    got = pick_cohorts(cells, "SYN1:SE:1", 2)
    want = cells[cells.adj_cox_status.eq("tested")].sort_values("adj_cox_p").cohort.head(2).tolist()
    assert got == want


def test_probe_outputs(ds, gtf_path, tmp_path):
    res = probe(ds, genes=["SYN1"], settings=FAST, gtf=gtf_path, out_dir=tmp_path, top=2, log=lambda *_: None)
    ev = res.events
    assert list(ev["rank"]) == [1, 2, 3] and set(ev.event_id) == {"SYN1:SE:1", "SYN1:A3SS:1", "SYN1:RI:1"}
    assert ev.adj_cox_p05.is_monotonic_decreasing or ev.iloc[0].adj_cox_p05 >= ev.iloc[1].adj_cox_p05
    for k in ("report", "events", "cells", "overview", "pdf"):
        assert res.paths[k].exists()
    assert len(list((tmp_path / "pages").glob("*.png"))) == 3
    report = res.paths["report"].read_text()
    assert "chance alone gives about" in report and "Benjamini" in report and "(pages/001_" in report
    cells = pd.read_csv(res.paths["cells"])
    assert {"cox_q", "adj_cox_q", "adj_hr_per_iqr", "group_hit"} <= set(cells.columns)
    assert cells.event_id.nunique() == 3 and cells.cohort.nunique() == 4


def test_probe_every_event_without_adjustment(ds, tmp_path):
    res = probe(ds, settings=FAST, adjusted=None, out_dir=tmp_path, max_pages=1, log=lambda *_: None)
    assert len(res.events) == 4 and "adj_cox_p" not in res.cells.columns
    assert len(list((tmp_path / "pages").glob("*.png"))) == 1


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
    fig = pd.read_csv(next((tmp_path / "pn").glob("*.csv")), low_memory=False)
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
    page = next((tmp_path / "pages").glob("*.png")).name
    assert page == "001_SYN1_SE_1.png"


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
    s = sa.Settings(formats=("png",), dpi=60, min_pairs=7, min_group=7, cox_min_events=10)
    assert s.changed() == {"min_pairs": (7, 10), "min_group": (7, 10), "cox_min_events": (10, 20)}
    assert sa.Settings(dpi=60, formats=("png",)).changed() == {}                # drawing settings do not count
    res = probe(ds, genes=["SYN1"], settings=s, out_dir=tmp_path, top=1, max_pages=1, log=lambda *_: None)
    report = res.paths["report"].read_text()
    assert "Settings changed from the defaults:** min_pairs 7 (default 10), min_group 7 (default 10), " \
           "cox_min_events 10 (default 20)" in report
    p = sa.event_panel(ds, "SYN1:SE:1", ["COH1"], "OS", settings=s, gex=False)
    assert "Settings changed from the defaults" in " ".join(t.get_text() for t in p.figure.texts)


def test_flags_are_reported(ds, tmp_path):
    s = sa.Settings(formats=("png",), dpi=60, narrow_psi_below=1.0)            # every PSI range counts as narrow
    res = probe(ds, genes=["SYN1"], settings=s, out_dir=tmp_path, top=1, max_pages=1, log=lambda *_: None)
    t = res.cells[res.cells.cox_status.eq("tested")]
    assert len(t) and t.psi_narrow.astype(bool).all() and "narrow PSI range (IQR" in t.cox_notes.iloc[0]
    a = t[t.adj_cox_status.eq("tested")]
    assert len(a) and a.adj_psi_narrow.astype(bool).all() and a.adj_cox_events_per_term.notna().all()
    rep = (tmp_path / "report.md").read_text()
    assert "## Flags on the Cox fits" in rep and f"{len(t)} of {len(t)} base-model fits" in rep


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
    p = sa.event_panel(ds, "SYN1:SE:1", ["COH1"], "OS", settings=s, gex=False, part=(1, 2))
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
