import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

import splice_assay as sa
from splice_assay.analysis import gene_fdr

FAST = sa.Settings(formats=("png",), dpi=60)


def test_gene_fdr_families():
    rng = np.random.default_rng(0)
    df = pd.DataFrame(dict(gene=np.repeat(["A", "B"], 12), endpoint=np.tile(np.repeat(["OS", "DSS"], 6), 2),
                           cox_status="tested", cox_p=rng.uniform(size=24)))
    df.loc[3, "cox_status"] = "coverage_gate"                       # not tested: outside every family
    out = gene_fdr(df, "cox", ["gene", "endpoint"], min_family=5)
    for (g, ep), d in out.groupby(["gene", "endpoint"]):
        t = d[d.cox_status.eq("tested")]
        assert np.allclose(t.cox_q, false_discovery_control(t.cox_p.to_numpy(), method="bh"))
        assert (t.cox_q_tests == len(t)).all()
    assert np.isnan(out.loc[3, "cox_q"])
    small = gene_fdr(df, "cox", ["gene", "endpoint"], min_family=6)  # the family with the untested row has 5
    assert small[small.gene.eq("A") & small.endpoint.eq("OS")].cox_q.isna().all()
    assert small[small.gene.eq("A") & small.endpoint.eq("DSS")].cox_q.notna().all()


def test_analyze_gives_q_within_each_gene(ds):
    s = sa.Settings(fdr_min_family=2)
    r = sa.analyze(ds, endpoints=["OS"], settings=s)
    for test, table, by in (("paired", r.groups, ["gene"]), ("unpaired", r.groups, ["gene"]),
                            ("km", r.survival, ["gene", "endpoint"]), ("cox", r.survival, ["gene", "endpoint"])):
        assert f"{test}_q" in table.columns
        for _, d in table[table[f"{test}_status"].eq("tested")].groupby(by):
            if len(d) >= 2:
                assert np.allclose(d[f"{test}_q"], false_discovery_control(d[f"{test}_p"].to_numpy(), method="bh"))
    # the default family minimum (10) leaves the synthetic genes' small families without q (SYN3's HIT-index events,
    # which would form a family of their own, are left out by default)
    g = sa.analyze(ds, endpoints=["OS"]).groups
    assert g.paired_q.isna().all() and not g.event_id.str.contains(":HIT:").any()


def test_hit_index_forms_its_own_families(ds):
    """Including the HIT index changes no PSI q value: HIT-index tests are families of their own."""
    s = sa.Settings(fdr_min_family=2)
    a = sa.analyze(ds, endpoints=["OS"], settings=s)
    b = sa.analyze(ds, endpoints=["OS"], settings=s, include_hit=True)
    for x, y in ((a.groups, b.groups), (a.survival, b.survival)):
        hit = y.event_id.str.contains(":HIT:")
        assert hit.any() and not x.event_id.str.contains(":HIT:").any()
        pd.testing.assert_frame_equal(x.reset_index(drop=True), y[~hit].reset_index(drop=True), check_dtype=False)
    h = b.groups[b.groups.event_id.str.contains(":HIT:") & b.groups.paired_status.eq("tested")]
    assert set(h.paired_q_tests) == {len(h)} and h.paired_q.notna().all()


def test_panel_prints_q(ds):
    s = sa.Settings(formats=("png",), dpi=60, fdr_min_family=2)
    m = sa.CoxModel().with_clinical(("age",))
    p = sa.event_panel(ds, "SYN1:SE:1", ["COH1", "COH3"], "OS", settings=s, detail=m)
    texts = " ".join(t.get_text() for t in p.figure.texts) + " ".join(t.get_text() for a in p.figure.axes
                                                                       for t in a.texts)
    assert "\nq " in texts and "PSI q " in texts and "Benjamini–Hochberg within SYN1 for OS" in texts
    assert "KM 12" in texts or "KM " in texts
    printed = p.table[p.table.panel.eq("printed")]
    assert {"paired q (BH within gene)", "log-rank q (BH within gene)"} <= set(printed.what)
    assert p.table[p.table.panel.eq("forest")].cox_q.notna().any()
