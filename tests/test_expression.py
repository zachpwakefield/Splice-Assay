import numpy as np
import pandas as pd

import splice_assay as sa
from splice_assay.analysis import analyse_expression
from splice_assay.cli import main
from splice_assay.probe import probe
from splice_assay.stats.survival import CoxModel, expression_cell

FAST = sa.Settings(formats=("png",), dpi=60)


def test_expression_statistics(ds):
    m = CoxModel().with_clinical(("age", "stage"), baseline={"stage": "I"})
    r = analyse_expression(ds, endpoints=["OS"], model=m)
    assert set(r.groups.gene) == {"SYN1", "SYN2", "SYN3"} and r.groups.cohort.nunique() == 4
    assert (r.groups.paired_n_dropped_0_1.fillna(0) == 0).all()              # no PSI 0/1 check on expression
    sv = r.survival.set_index(["gene", "cohort"])
    one = sv.loc[("SYN1", "COH1")]                                           # simulated: higher expression, worse OS
    assert one.cox_status == "tested" and one.hr_per_sd > 1 and one.cox_p < 0.05
    assert one.cox_model == "expression + age + stage" and one.km_status == "tested"
    t = r.cox_terms[(r.cox_terms.gene == "SYN1") & (r.cox_terms.cohort == "COH1")]
    assert list(t.kind[:2]) == ["gex", "numeric"] and np.isclose(t.hr.iloc[0], one.hr_per_sd)
    assert sv.loc[("SYN1", "COH4")].cox_status == "tested"                  # no normals: survival still tested
    g = r.groups.set_index(["gene", "cohort"]).loc[("SYN1", "COH4")]
    assert g.unpaired_status == "no_reference_samples"


def test_expression_cell_gates():
    s, n = sa.Settings(), 60
    t, e, pos = np.linspace(100, 3000, n), np.r_[np.ones(30, int), np.zeros(30, int)], np.arange(n)
    assert expression_cell(np.full(n, np.nan), pos, t, e, s)[0]["cox_status"] == "coverage_gate"
    assert expression_cell(np.full(n, 5.0), pos, t, e, s)[0]["km_status"] == "coverage_gate"   # nothing off the mode
    rng = np.random.default_rng(3)                                        # higher expression, shorter survival
    x = rng.normal(5, 1, 200)
    tt, ee = rng.exponential(1000 * np.exp(-0.5 * (x - 5))), (rng.uniform(size=200) < 0.7).astype(int)
    row, terms = expression_cell(x, np.arange(200), tt, ee, s)
    assert row["cox_status"] == "tested" and terms[0]["kind"] == "gex" and row["hr_per_sd"] > 1


def test_expression_page(ds, results, tables, tmp_path):
    """The host gene's expression has its own page: a forest of every cohort, then the cohorts given in full; the
    splicing page has no expression rows."""
    m = CoxModel().with_clinical(("age",))
    p = sa.expression_panel(ds, "SYN1", ["COH1", "COH4"], "OS", settings=FAST, model=m, out_dir=tmp_path)
    cells = p.table[p.table.panel.eq("gex_cell")].set_index("cohort")
    assert set(cells.index) == {"COH1", "COH4"} and cells.at["COH1", "hr_per_sd"] > 1
    assert {"gex_pairs", "gex_all_samples", "gex_km_curve", "gex_cox_detail", "gex_forest"} <= set(p.table.panel)
    forest = p.table[p.table.panel.eq("gex_forest")]
    assert set(forest.cohort) == {"COH1", "COH2", "COH3", "COH4"} and forest.shade.sum() == 2   # all cohorts, 2 drawn
    assert "gex_groups_not_tested" in set(p.table.panel)                      # COH4 has no normals
    assert p.paths["png"].name == "SYN1_expression_COH1_COH4_OS.png" and p.provenance["call"]["gene"] == "SYN1"
    texts = [t.get_text() for t in p.figure.texts]
    assert "COH1 · SYN1 expression" in texts and any("expression · COH1 and COH4 · OS" in t for t in texts)
    q = sa.event_panel(ds, "SYN1:SE:1", ["COH1", "COH4"], "OS", settings=FAST, results=results, detail=m)
    assert not q.table.panel.str.startswith("gex_").any()
    assert sa.expression_panel(sa.Dataset.from_tables(**{**tables, "expression": None}), "SYN1", ["COH1"], "OS",
                               settings=FAST) is None


def test_probe_and_cli_expression(ds, tmp_path, capsys):
    res = probe(ds, genes=["SYN1"], settings=FAST, out_dir=tmp_path / "pr", top=1, max_pages=1, log=lambda *_: None)
    x = pd.read_csv(res.paths["expression"])
    assert set(x.gene) == {"SYN1"} and {"hr_per_sd", "km_p", "unpaired_p"} <= set(x.columns)
    assert "## Host-gene expression" in res.paths["report"].read_text()
    out = tmp_path / "ex"
    assert main(["example", str(out)]) == 0
    fast = tmp_path / "fast.json"
    fast.write_text('{"formats": ["png"], "dpi": 60}')
    for flag, want in (([], True), (["--no-gex"], False)):
        o = tmp_path / f"pn{len(flag)}"
        assert main(["panel", str(out / "data"), "--event", "SYN1:SE:1", "--cohort", "COH1", "--no-detail", "--out",
                     str(o), "--settings", str(fast)] + flag) == 0
        assert bool(list(o.glob("SYN1_expression_*.png"))) == want             # its own page, or none
        for f in o.glob("SYN1_SE_1_*.csv"):
            assert not pd.read_csv(f, low_memory=False).panel.astype(str).str.startswith("gex_").any()
    assert res.paths["expression_pages"].name == "SYN1_expression.png"
