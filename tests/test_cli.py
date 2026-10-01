import json

import pandas as pd

from splice_assay.cli import main


def test_example_then_commands(tmp_path, capsys):
    out = tmp_path / "ex"
    assert main(["example", str(out)]) == 0
    data = out / "data"
    assert (out / "results" / "survival.csv").exists() and len(list((out / "figures").glob("*.png"))) == 4
    assert main(["validate", str(data)]) == 0
    assert "COH1" in capsys.readouterr().out
    assert main(["analyse", str(data), "--out", str(tmp_path / "res"), "--event", "SYN1:SE:1",
                 "--endpoint", "OS"]) == 0
    sv = pd.read_csv(tmp_path / "res" / "survival.csv")
    assert set(sv.endpoint) == {"OS"} and set(sv.event_id) == {"SYN1:SE:1"}
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
