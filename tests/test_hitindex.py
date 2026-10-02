import numpy as np
import pandas as pd
import pytest

import splice_assay as sa
from splice_assay import InputError, example
from splice_assay.cli import main
from splice_assay.hitindex import import_hitindex, parse_id

FAST = sa.Settings(formats=("png",), dpi=60)
COLS = ["event_id", "gene", "gene_id", "chrom", "strand", "event_type", "constant", "variable", "label", "source_id"]


def test_parse_hitindex_ids():
    assert parse_id("ENSG00000204371.13;chr6:31879759-31880123;afe") == ("ENSG00000204371", "chr6", 31879758,
                                                                          31880123, "AFE")
    assert parse_id("SYNG0003.1;chr3:303001-303100;HIT")[4] == "HIT"
    for bad in ("ENSG1;chr6:10-20;se", "ENSG1#chr6:10-20#-", "ENSG1;chr6:20-10;afe"):
        with pytest.raises(InputError):
            parse_id(bad)


def test_import_reproduces_the_example(tmp_path):
    """The example's HITindex matrices import to exactly its SYN3 events and values (strand and symbol from the GTF;
    numbered 5' to 3'; AFE/ALE carry the gene's other first/last exons)."""
    example.write(tmp_path)
    mats = sorted((tmp_path / "hitindex").glob("*.csv"))
    ev, psi = import_hitindex(mats, tmp_path / "data" / "annotation.gtf")
    t = example.make()
    want = t["events"][t["events"].event_type.isin(["AFE", "ALE", "HIT"])].reset_index(drop=True)
    pd.testing.assert_frame_equal(ev[COLS], want[COLS])
    got = psi.set_index(["event_id", "sample_id"]).psi
    ref = t["psi"].dropna(subset=["psi"]).set_index(["event_id", "sample_id"]).psi.loc[got.index]
    assert np.array_equal(got.to_numpy(), ref.to_numpy())
    assert psi.psi[psi.event_id.str.contains(":HIT:")].min() < 0                # the HIT index is signed
    one, _ = import_hitindex(mats[2:], tmp_path / "data" / "annotation.gtf", genes=["SYNG0003"])
    assert list(one.event_id) == ["SYN3:HIT:0001", "SYN3:HIT:0002"]
    with pytest.raises(InputError, match="none of"):
        import_hitindex(mats, tmp_path / "data" / "annotation.gtf", genes=["NOPE"])


@pytest.mark.parametrize("layout", ["events+long", "events+wide", "one_table_wide", "one_table_long"])
def test_import_cli_and_append(tmp_path, capsys, monkeypatch, layout):
    """--append adds HITindex events to rMATS-like tables in whichever layout the folder has; the rows already there
    are kept exactly."""
    monkeypatch.delenv("SPLICE_ASSAY_GTF", raising=False)
    example.write(tmp_path)
    gtf = str(tmp_path / "data" / "annotation.gtf")
    mats = [str(m) for m in sorted((tmp_path / "hitindex").glob("*.csv"))]
    t = example.make()
    rm = t["events"][~t["events"].event_type.isin(["AFE", "ALE", "HIT"])].drop(columns="source_id")
    psi = t["psi"][t["psi"].event_id.isin(rm.event_id)]
    wide = psi.pivot(index="event_id", columns="sample_id", values="psi").reset_index()
    out = tmp_path / "tables"
    out.mkdir()
    if layout.startswith("events+"):
        rm.to_csv(out / "events.csv", index=False)
        (psi if layout == "events+long" else wide).to_csv(out / "psi.csv", index=False)
    elif layout == "one_table_wide":
        rm.merge(wide, on="event_id").to_csv(out / "psi.csv", index=False)
    else:
        psi.merge(rm, on="event_id").to_csv(out / "psi.csv", index=False)
    text = (lambda: pd.read_csv(out / "psi.csv", dtype=str, keep_default_na=False))
    before = text()
    assert main(["import-hitindex", *mats, "--gtf", gtf, "--out", str(out), "--append"]) == 0
    after = text()                                    # the old rows first, exactly as they were; new columns may follow
    pd.testing.assert_frame_equal(after.iloc[:len(before)][list(before.columns)], before)
    for name in ("samples", "survival", "clinical", "expression"):
        t[name].to_csv(out / f"{name}.csv", index=False)
    ds = sa.Dataset.from_dir(out)
    assert len(ds.events) == 10 and set(ds.events.event_type) >= {"AFE", "ALE", "HIT", "SE"}
    assert ds.psi.loc["SYN3:HIT:0001"].min() < 0
    ref = sa.Dataset.from_tables(**t)
    pd.testing.assert_frame_equal(ds.psi.sort_index(), ref.psi.loc[ds.psi.index].sort_index()[ds.psi.columns],
                                  check_names=False)
    capsys.readouterr()
    assert main(["import-hitindex", *mats, "--gtf", gtf, "--out", str(out), "--append"]) == 2   # same IDs again
    assert "already in" in capsys.readouterr().err


def test_import_cli(tmp_path, capsys, monkeypatch):
    monkeypatch.delenv("SPLICE_ASSAY_GTF", raising=False)
    example.write(tmp_path)
    mats = [str(m) for m in sorted((tmp_path / "hitindex").glob("*.csv"))]
    gtf = str(tmp_path / "data" / "annotation.gtf")
    assert main(["import-hitindex", *mats, "--gtf", gtf, "--out", str(tmp_path / "fresh")]) == 0
    out = capsys.readouterr().out
    assert "imported: 2 AFE, 2 ALE, 2 HIT" in out and "next: add samples.csv" in out
    assert {p.name for p in (tmp_path / "fresh").iterdir()} == {"events.csv", "psi.csv"}
    assert main(["import-hitindex", *mats, "--gtf", gtf, "--gene", "SYN3", "--one-table",
                 "--out", str(tmp_path / "one")]) == 0
    assert {p.name for p in (tmp_path / "one").iterdir()} == {"psi.csv"}
    assert main(["import-hitindex", *mats, "--out", str(tmp_path / "x")]) == 2       # no GTF, no strand
    assert "needs --gtf" in capsys.readouterr().err
    (tmp_path / "tsv").mkdir()
    pd.DataFrame(dict(event_id=["a"], sample_id=["s"], psi=[0.5])).to_csv(tmp_path / "tsv" / "psi.tsv", sep="\t")
    assert main(["import-hitindex", *mats, "--gtf", gtf, "--out", str(tmp_path / "tsv"), "--append"]) == 2
    assert "adds to CSV tables" in capsys.readouterr().err


def test_hit_index_has_its_own_rules(tables):
    """HIT values may be negative; no 0/1 robustness check; a hit needs |delta| > hit_min_abs_delta (0.20)."""
    ds = sa.Dataset.from_tables(**tables)
    g = sa.analyse(ds, events=["SYN3:HIT:0001"]).groups.set_index("cohort").loc["COH1"]
    assert g.paired_hit_status == "hit" and g.paired_state_01 == "not_applicable" and g.paired_delta_median > 0.2
    strict = sa.analyse(ds, events=["SYN3:HIT:0001"], settings=sa.Settings(hit_min_abs_delta=0.5))
    assert strict.groups.set_index("cohort").at["COH1", "paired_hit_status"] == "small_effect"
    sv = sa.analyse(ds, events=["SYN3:HIT:0002"], endpoints=["OS"])
    assert (sv.cox_terms.term.eq("HIT index") & sv.cox_terms.kind.eq("psi_iqr")).any()
    assert sv.survival.cox_model.str.startswith("HIT index +").all()
    bad = tables["psi"].copy()
    bad.loc[bad.event_id.eq("SYN3:HIT:0001"), "psi"] = 1.5
    with pytest.raises(InputError, match=r"\[-1, 1\] \(HIT index\)"):
        sa.Dataset.from_tables(**dict(tables, psi=bad))
    bad = tables["psi"].copy()
    bad.loc[bad.event_id.eq("SYN1:SE:1"), "psi"] = -0.2
    with pytest.raises(InputError, match=r"\[0, 1\]"):
        sa.Dataset.from_tables(**dict(tables, psi=bad))


def test_hit_and_afe_pages(tables, tmp_path):
    ds = sa.Dataset.from_tables(**tables)
    p = sa.event_panel(ds, "SYN3:HIT:0002", ["COH1"], "OS", settings=FAST)
    texts = [t.get_text() for t in p.figure.texts] + [t.get_text() for a in p.figure.axes for t in a.texts]
    assert any(t.startswith("split at median HIT index") for t in texts)
    assert any("HIT index from −1" in t for t in texts) and any("Region measured by HIT index" == t for t in texts)
    assert not any(t.startswith("Junction of") for t in texts)                 # no arcs, so no arc legend
    a = sa.event_panel(ds, "SYN3:AFE:0001", ["COH1"], "OS", settings=FAST)
    d = a.table[a.table.panel.eq("event_detail")].text.iloc[0]
    assert d.startswith("AFE:0001: alternative first exon chr3:300,001–300,200 (200 nt) (1 other first exon drawn)")
    with pytest.raises(InputError, match="same quantity"):
        sa.event_panel(ds, ["SYN3:HIT:0001", "SYN3:AFE:0001"], ["COH1"], "OS", settings=FAST)


def test_hit_is_left_out_unless_asked(tables, tmp_path, capsys):
    """analyse and probe leave out HIT-index events unless include_hit; a HIT event named explicitly is analysed, with
    its gene's HIT events as its q family."""
    from splice_assay.probe import probe, select_events
    ds = sa.Dataset.from_tables(**tables)
    hits = [e for e in ds.psi.index if ":HIT:" in e]
    assert select_events(ds, genes=["SYN3"]) == ["SYN3:AFE:0001", "SYN3:AFE:0002", "SYN3:ALE:0001", "SYN3:ALE:0002"]
    assert select_events(ds, genes=["SYN3"], include_hit=True)[-2:] == hits
    assert select_events(ds, events=hits[:1]) == hits[:1]
    only = sa.Dataset.from_tables(**tables, event_ids=hits)
    with pytest.raises(InputError, match="only HIT-index events for SYN3: add --include-hit"):
        select_events(only, genes=["SYN3"])
    with pytest.raises(InputError, match="only HIT-index events in the data"):
        sa.analyse(only)
    assert set(sa.analyse(only, include_hit=True, endpoints=["OS"], cohorts=["COH1"]).groups.event_id) == set(hits)
    res = probe(ds, genes=["SYN3"], settings=FAST, adjusted=None, out_dir=tmp_path / "p", max_pages=1,
                include_hit=True, log=lambda *_: None)
    assert len(res.events) == 6 and "HIT-index events form families of their own" in res.paths["report"].read_text()
    plain = sa.analyse(ds, endpoints=["OS"], cohorts=["COH1"])
    with pytest.raises(InputError, match="include_hit=True"):
        sa.event_panel(ds, hits[1], ["COH1"], "OS", results=plain, settings=FAST)
    texts = (lambda p: " ".join(t.get_text() for t in p.figure.texts))
    assert "over all its HIT-index events × cohorts" in texts(sa.event_panel(ds, hits[1], ["COH1"], "OS",
                                                                               settings=FAST))
    assert "over all its PSI events (HIT index apart) × cohorts" in texts(
        sa.event_panel(ds, "SYN3:AFE:0001", ["COH1"], "OS", settings=FAST))
    assert "over all its events × cohorts" in texts(sa.event_panel(ds, "SYN1:SE:1", ["COH1"], "OS", settings=FAST))
    # the command line
    example.write(tmp_path / "ex")
    data = str(tmp_path / "ex" / "data")
    assert main(["validate", data]) == 0
    assert "note: 2 HIT-index event(s); analyse and probe leave them out unless --include-hit" in capsys.readouterr().out
    for flag, n in (([], 0), (["--include-hit"], 2)):
        assert main(["analyse", data, "--out", str(tmp_path / f"a{n}"), "--cohort", "COH1", "--endpoint", "OS",
                     *flag]) == 0
        sv = pd.read_csv(tmp_path / f"a{n}" / "survival.csv")
        assert sv.event_id.str.contains(":HIT:").sum() == n
