import json

import pandas as pd

from splice_assay.cli import main


def test_example_then_commands(tmp_path, capsys):
    out = tmp_path / "ex"
    assert main(["example", str(out)]) == 0
    data = out / "data"
    assert (out / "results" / "survival.csv").exists() and len(list((out / "figures").glob("*.png"))) == 6
    assert {p.name for p in (out / "data").glob("*.csv")} == {"samples.csv", "psi.csv", "expression.csv"}   # two + 1
    assert main(["validate", str(data)]) == 0
    assert "COH1" in capsys.readouterr().out
    assert main(["analyze", str(data), "--out", str(tmp_path / "res"), "--event", "SYN1:SE:1",
                 "--endpoint", "OS"]) == 0
    sv = pd.read_csv(tmp_path / "res" / "survival.csv")
    assert set(sv.endpoint) == {"OS"} and set(sv.event_id) == {"SYN1:SE:1"}
    assert set(sv.km_split) == {"median"}
    assert main(["analyze", str(data), "--out", str(tmp_path / "res7"), "--event", "SYN1:SE:1", "--endpoint", "OS",
                 "--km-split", "0.7"]) == 0
    sv = pd.read_csv(tmp_path / "res7" / "survival.csv")
    assert set(sv.cutoff) == {0.7} and set(sv.km_split) == {"set"}
    assert main(["analyze", str(data), "--out", str(tmp_path / "x"), "--km-split-expression", "max"]) == 2
    assert 'km_split_expression must be "median", "mean" or a number' in capsys.readouterr().err
    (tmp_path / "s.json").write_text(json.dumps({"formats": ["svg"], "dpi": 100}))
    assert main(["panel", str(data), "--event", "SYN1:RI:1", "--cohort", "COH1", "--endpoint", "OS",
                 "--gtf", str(data / "annotation.gtf"), "--out", str(tmp_path / "fig"),
                 "--settings", str(tmp_path / "s.json")]) == 0
    assert (tmp_path / "fig" / "SYN1_RI_1_COH1_OS.svg").exists()
    assert main(["panel", str(data), "--event", "SYN1:SE:1", "--cohort", "COH1", "--cohort", "COH3", "--endpoint", "OS",
                 "--out", str(tmp_path / "assay"), "--settings", str(tmp_path / "s.json"), "--detail-covariate", "age",
                 "--detail-covariate", "stage", "--baseline", "stage=I"]) == 0
    fig = pd.read_csv(tmp_path / "assay" / "SYN1_SE_1_COH1_COH3_OS.csv", low_memory=False)
    band = fig[fig.panel.eq("cox_detail")]
    assert set(band.tag) == {"SYN1:SE:1|COH1", "SYN1:SE:1|COH3"} and set(band.reference.dropna()) >= {"I"}
    assert fig[fig.panel.eq("forest_axis")].cox_model.iloc[0] == "Cox: PSI + host expression"
    pd.DataFrame([dict(events="SYN1:SE:1;SYN1:A3SS:1", cohorts="COH1", endpoint="DSS", stem="two"),
                  dict(events="SYN2:MXE:1", cohorts="COH3;COH1", endpoint="OS", stem="")]).to_csv(
        tmp_path / "spec.csv", index=False)
    assert main(["panels", str(data), "--spec", str(tmp_path / "spec.csv"), "--gtf", str(data / "annotation.gtf"),
                 "--out", str(tmp_path / "many"), "--settings", str(tmp_path / "s.json")]) == 0
    assert (tmp_path / "many" / "two.svg").exists()
    assert (tmp_path / "many" / "SYN2_MXE_1_COH3_COH1_OS.svg").exists()
    assert main(["cox", str(data), "--event", "SYN1:SE:1", "--cohort", "COH1", "--endpoint", "OS",
                 "--covariate", "age", "--covariate", "stage", "--strata", "sex", "--out", str(tmp_path / "cox"),
                 "--settings", str(tmp_path / "s.json")]) == 0
    out = capsys.readouterr().out
    assert "PSI + host expression + age + stage; strata: sex" in out and "vs I" in out
    assert (tmp_path / "cox" / "SYN1_SE_1_COH1_OS_cox.csv").exists()


def test_column_mapping_on_the_command_line(tmp_path, capsys):
    from splice_assay import example
    t = example.make()
    d = tmp_path / "d"
    d.mkdir()
    t["samples"].rename(columns={"cohort": "cancer", "group": "sample_type"}).assign(
        sample_type=lambda x: x.sample_type.map({"tumour": "Primary Tumor", "normal": "Solid Tissue Normal"})).to_csv(
        d / "samples.csv", index=False)
    for k in ("psi", "events", "survival"):
        t[k].to_csv(d / f"{k}.csv", index=False)
    assert main(["validate", str(d), "--column", "cohort=cancer", "--column", "group=sample_type",
                 "--case", "Primary Tumor", "--reference", "Solid Tissue Normal"]) == 0
    assert "case = Primary Tumor, reference = Solid Tissue Normal" in capsys.readouterr().out


def test_input_errors_exit_2(tmp_path, capsys):
    tmp_path.joinpath("samples.csv").write_text("sample_id,patient_id,cohort\na,b,c\n")
    tmp_path.joinpath("psi.csv").write_text("event_id,sample_id,psi\ne,a,0.5\n")
    tmp_path.joinpath("events.csv").write_text("event_id,gene\ne,G\n")
    assert main(["validate", str(tmp_path)]) == 2
    assert "missing column(s) group" in capsys.readouterr().err


def test_example_as_separate_tables(tmp_path, capsys):
    """The default two-table example and --separate give the same data and the same statistics."""
    import splice_assay as sa
    from splice_assay import example
    assert main(["validate", str(example.write(tmp_path / "two")["samples"].parent)]) == 0
    assert "survival: from the samples table (DSS, OS)" in capsys.readouterr().out
    sep = example.write(tmp_path / "sep", separate=True)
    assert {p.name for p in sep["samples"].parent.glob("*.csv")} == {f"{t}.csv" for t in
                                                                    ("samples", "psi", "events", "survival",
                                                                     "clinical", "expression")}
    a, b = (sa.Dataset.from_dir(tmp_path / f / "data") for f in ("two", "sep"))
    key = ["event_id", "cohort", "endpoint"]
    x, y = (sa.analyze(d, cohorts=["COH1"]).survival.sort_values(key).reset_index(drop=True) for d in (a, b))
    pd.testing.assert_frame_equal(x, y)


def _example_data(tmp_path):
    from splice_assay import example
    example.write(tmp_path / "ex")
    fast = tmp_path / "fast.json"
    fast.write_text('{"formats": ["png"], "dpi": 60}')
    return str(tmp_path / "ex" / "data"), str(fast)


def test_settings_errors_are_input_errors(tmp_path, capsys):
    data, _ = _example_data(tmp_path)
    for text, msg in (('{"min_pair": 5}', "unknown setting(s): min_pair"), ("{", "--settings"),
                      ('{"min_pairs": -1}', "min_pairs must be a non-negative integer"),
                      ('{"min_pairs": Infinity}', "--settings"), ("[1, 2]", "expected an object of settings")):
        (tmp_path / "bad.json").write_text(text)
        assert main(["analyze", data, "--out", str(tmp_path / "r"), "--settings", str(tmp_path / "bad.json")]) == 2
        assert msg in capsys.readouterr().err
        assert main(["validate", data, "--settings", str(tmp_path / "bad.json")]) == 2       # validate checks them
        capsys.readouterr()
    (tmp_path / "ok.json").write_text('{"min_pairs": 5}')
    assert main(["validate", data, "--settings", str(tmp_path / "ok.json")]) == 0
    assert "settings changed from the defaults: min_pairs 5 (default 10)" in capsys.readouterr().out


def test_a_psi_table_given_alone_keeps_the_folder_events_table(tmp_path, capsys):
    """--table psi=FILE with a folder holding events.csv: the events come from the folder, as in Dataset.from_dir."""
    from splice_assay import example
    data = example.write(tmp_path / "sep", separate=True)["samples"].parent
    other = tmp_path / "elsewhere.csv"
    pd.read_csv(data / "psi.csv").to_csv(other, index=False)
    (tmp_path / "fast.json").write_text('{"formats": ["png"], "dpi": 60}')
    assert main(["panel", str(data), "--table", f"psi={other}", "--event", "SYN1:SE:1", "--cohort", "COH1",
                 "--endpoint", "OS", "--no-gex", "--no-detail", "--out", str(tmp_path / "fig"), "--settings",
                 str(tmp_path / "fast.json")]) == 0
    assert main(["probe", str(data), "--table", f"psi={other}", "--gene", "SYN1", "--no-adjust", "--no-gex",
                 "--max-pages", "0", "--out", str(tmp_path / "probe"), "--settings", str(tmp_path / "fast.json")]) == 0
    capsys.readouterr()


def test_results_keep_their_model_and_settings(tmp_path, capsys):
    import pytest
    import splice_assay as sa
    data, _ = _example_data(tmp_path)
    assert main(["analyze", data, "--out", str(tmp_path / "res"), "--event", "SYN1:SE:1", "--covariate", "age",
                 "--no-expression", "--hr-unit", "sd"]) == 0
    assert "analysis: " in capsys.readouterr().out
    back = sa.Results.read(tmp_path / "res")
    assert back.model == sa.CoxModel(expression=False, covariates=("age",)) and back.settings.psi_hr_unit == "sd"
    import numpy as np                                          # numpy numbers in the settings are written as numbers
    res = sa.analyze(sa.Dataset.from_dir(data, event_ids=["SYN1:SE:1"]), settings=sa.Settings(min_pairs=np.int64(8)))
    res.write(tmp_path / "np")
    assert sa.Results.read(tmp_path / "np").settings.min_pairs == 8
    assert set(back.survival.cox_model.dropna()) == {"PSI + age"}
    (tmp_path / "res" / "analysis.json").unlink()                     # as written by earlier versions
    with pytest.warns(UserWarning, match="the defaults are assumed"):
        assert sa.Results.read(tmp_path / "res").model == sa.CoxModel()


def test_validate_notes_patients_with_several_samples(tmp_path, capsys):
    import splice_assay as sa
    from splice_assay import example
    t = example.make()
    s = t["samples"]
    extra = s[s.group.eq("tumour") & s.cohort.eq("COH4")].head(3).assign(sample_id=lambda d: d.sample_id + "-b",
                                                                     survival_cohort=False)   # COH4: no normals
    s2 = pd.concat([s.assign(survival_cohort=s.group.eq("tumour")), extra], ignore_index=True)
    d = tmp_path / "d"
    d.mkdir()
    s2.to_csv(d / "samples.csv", index=False)
    for k in ("psi", "events", "survival", "clinical"):
        t[k].to_csv(d / f"{k}.csv", index=False)
    assert main(["validate", str(d)]) == 0
    assert "note: patients 3 with several Tumour samples in one cohort" in capsys.readouterr().out
    assert sa.Dataset.from_dir(d).samples.patient_id.duplicated().any()


def test_panels_detail_adds_the_clinical_terms_found(tmp_path, capsys):
    """panels --detail draws the model rows as panel does: the page's model plus the age, sex and stage found."""
    data, fast = _example_data(tmp_path)
    pd.DataFrame([dict(events="SYN1:SE:1", cohorts="COH1", endpoint="OS", stem="one")]).to_csv(
        tmp_path / "spec.csv", index=False)
    assert main(["panels", data, "--spec", str(tmp_path / "spec.csv"), "--out", str(tmp_path / "many"), "--settings",
                 fast, "--detail", "--no-gex"]) == 0
    assert "model rows: Cox PSI + host expression + age + sex + stage; found" in capsys.readouterr().out
    rows = pd.read_csv(tmp_path / "many" / "one.csv", low_memory=False)
    assert {"PSI", "age", "sex", "stage"} <= set(rows[rows.panel.eq("cox_detail")].term.dropna())
    assert main(["panels", data, "--spec", str(tmp_path / "spec.csv"), "--out", str(tmp_path / "plain"), "--settings",
                 fast, "--no-gex"]) == 0                         # without --detail: no model rows
    assert "model rows" not in capsys.readouterr().out
    assert not pd.read_csv(tmp_path / "plain" / "one.csv", low_memory=False).panel.eq("cox_detail").any()


def test_categorical_implies_a_covariate_and_a_clash_is_an_input_error(tmp_path, capsys):
    data, _ = _example_data(tmp_path)
    base = ["cox", data, "--event", "SYN1:SE:1", "--cohort", "COH1", "--endpoint", "OS"]
    assert main(base + ["--categorical", "stage"]) == 0                    # stage enters, as categories
    assert "Cox PSI + host expression + stage" in capsys.readouterr().out
    assert main(base + ["--covariate", "sex", "--strata", "sex"]) == 2
    assert "a column cannot be both a covariate and a stratum: sex" in capsys.readouterr().err


def test_the_earlier_spelling_still_works(tmp_path, capsys):
    """analyze is the name; analyse (the earlier spelling) still runs, on the command line and in Python."""
    import splice_assay as sa
    from splice_assay.analysis import analyse_expression, analyze_expression
    assert sa.analyse is sa.analyze and analyse_expression is analyze_expression
    data, _ = _example_data(tmp_path)
    assert main(["analyse", data, "--out", str(tmp_path / "old"), "--event", "SYN1:SE:1", "--endpoint", "OS"]) == 0
    assert (tmp_path / "old" / "survival.csv").exists()
    capsys.readouterr()


def test_panel_forest_is_the_same_with_and_without_model_rows(tmp_path, capsys):
    """--covariate stage on raw values ("Stage IA"): the model rows' cleaning must not change the forest."""
    data, fast = _example_data(tmp_path)
    s = pd.read_csv(f"{data}/samples.csv")
    s["stage"] = "Stage " + s.stage.astype(str) + "A"
    s.to_csv(f"{data}/samples.csv", index=False)
    base = ["panel", data, "--event", "SYN1:SE:1", "--cohort", "COH1", "--endpoint", "OS", "--no-gex", "--covariate",
            "stage", "--settings", fast]
    hr = []
    for extra, out in ((["--no-detail"], "plain"), ([], "rows")):
        assert main(base + extra + ["--out", str(tmp_path / out)]) == 0
        f = pd.read_csv(next((tmp_path / out).glob("*.csv")), low_memory=False)
        hr.append(f[f.panel.eq("forest")].set_index("cohort").hr_per_sd)
    pd.testing.assert_series_equal(hr[0], hr[1])
    assert "found age = age" in capsys.readouterr().out                     # stage is the page's own: not re-added


def test_ridge_and_events_per_term_on_the_command_line(tmp_path, capsys):
    data, fast = _example_data(tmp_path)
    assert main(["validate", data, "--ridge", "clinical", "--ridge-penalty", "2"]) == 0
    assert "cox_ridge clinical (default none), cox_ridge_penalty 2 (default 1)" in capsys.readouterr().out
    assert main(["cox", data, "--event", "SYN1:SE:1", "--cohort", "COH1", "--endpoint", "OS", "--covariate", "age",
                 "--ridge", "clinical"]) == 0
    assert "Cox PSI + host expression + age; ridge λ 1 (clinical terms)" in capsys.readouterr().out
    assert main(["analyze", data, "--out", str(tmp_path / "r"), "--event", "SYN1:SE:1", "--endpoint", "OS",
                 "--min-events-per-term", "1000"]) == 0
    sv = pd.read_csv(tmp_path / "r" / "survival.csv")
    assert set(sv.cox_status) <= {"too_few_events_per_term", "too_few_patients_or_events", "coverage_gate"}
    assert sv.cox_status.eq("too_few_events_per_term").any()
    import pytest
    with pytest.raises(SystemExit):                                          # argparse refuses an unknown scope
        main(["analyze", data, "--out", str(tmp_path / "x"), "--ridge", "lasso"])
    capsys.readouterr()
