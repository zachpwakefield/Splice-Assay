import shutil

import pandas as pd
import pytest

import splice_assay as sa
from splice_assay import InputError, example
from splice_assay.cli import main
from splice_assay.events import geometry
from splice_assay.probe import probe
from splice_assay.protein import (Form, ProteinCache, build_cache, compatible, forms, labels, protein_change,
                                  protein_changes, tsl_rank)

FAST = sa.Settings(formats=("png",), dpi=60)


@pytest.fixture(scope="module")
def cache_dir(tmp_path_factory):
    d = tmp_path_factory.mktemp("proteins")
    for k, v in example.protein_cache().items():
        v.to_csv(d / f"{k}.csv", index=False)
    return d


@pytest.fixture(scope="module")
def cache(cache_dir):
    return ProteinCache.load(cache_dir)


def _event(e, **change):
    r = dict(example.EVENTS.set_index("event_id").loc[e].to_dict(), **change)
    return r, geometry(e, r)


def _tx(*exons):
    return pd.DataFrame(dict(start=[a for a, _ in exons], end=[b for _, b in exons]))


def test_forms_follow_spliceimpactr():
    f = forms(_event("SYN1:SE:1")[1])                          # 1-based closed
    assert f["INC"] == Form(((102001, 102100), (103001, 103150), (104001, 104300)), ())
    assert f["EXC"] == Form(((102001, 102100), (104001, 104300)), ((103001, 103150),))
    f = forms(_event("SYN1:A3SS:1")[1])                        # long / short exon with the flanking exon
    assert f["INC"].inc == ((100001, 100200), (101971, 102100)) and f["EXC"].exc == ((101971, 102000),)
    assert f["EXC"].inc == ((100001, 100200), (102001, 102100))
    assert forms(_event("SYN1:RI:1")[1])["EXC"].exc == ((103151, 104000),)
    f = forms(_event("SYN2:MXE:1")[1])                         # minus strand: the PSI exon is required in INC
    assert (204001, 204120) in f["INC"].inc and f["INC"].exc == ((203001, 203090),)
    # from the event's own exons: SE and RI have no exclusion form
    assert forms(_event("SYN1:SE:1")[1], flanks=False) == {"INC": Form(((103001, 103150),), ())}
    assert set(forms(_event("SYN2:MXE:1")[1], flanks=False)) == {"INC", "EXC"}


def test_compatibility_reasons_match_spliceimpactr():
    se = Form(((101, 200), (301, 400), (501, 600)), ())
    assert compatible(_tx((101, 200), (301, 400), (501, 600)), se, "SE", "INC", "+").context == "junction_chain"
    assert compatible(_tx((51, 200), (301, 400), (501, 700)), se, "SE", "INC", "+").ok    # outer ends are free
    assert compatible(_tx((101, 200), (301, 400), (450, 460), (501, 600)), se, "SE", "INC", "+").reason == \
        "splice_junction_chain_mismatch"                                                 # an exon in between
    assert compatible(_tx((101, 200), (301, 390), (501, 600)), se, "SE", "INC", "+").reason == \
        "required_interval_not_covered"
    skip = Form(((101, 200), (501, 600)), ((301, 400),))
    assert compatible(_tx((101, 200), (301, 400), (501, 600)), skip, "SE", "EXC", "+").reason == \
        "forbidden_interval_overlap"
    assert compatible(_tx((101, 200), (501, 600)), skip, "SE", "EXC", "+").ok
    one = Form(((301, 400),), ())                                                         # the exon alone
    assert compatible(_tx((101, 200), (301, 400), (501, 600)), one, "SE", "INC", "+").context == "partial_event"
    assert compatible(_tx((301, 400), (501, 600)), one, "SE", "INC", "+").reason == "internal_exon_boundary_mismatch"
    a5 = Form(((101, 230),), ())                       # A5SS, plus strand: ends at the site, a downstream exon
    assert compatible(_tx((101, 230), (501, 600)), a5, "A5SS", "INC", "+").context == "splice_site_only"
    assert compatible(_tx((101, 260), (501, 600)), a5, "A5SS", "INC", "+").reason == \
        "alternative_splice_boundary_mismatch"
    assert compatible(_tx((101, 230),), a5, "A5SS", "INC", "+").reason == "alternative_splice_boundary_mismatch"
    ri = Form(((101, 200), (201, 300), (301, 400)), ())
    assert compatible(_tx((101, 400),), ri, "RI", "INC", "+").context == "retained_intron"
    assert compatible(_tx((101, 200), (301, 400)), ri, "RI", "INC", "+").reason == "required_interval_not_covered"
    assert [tsl_rank(x) for x in ("1", "2 (assigned to previous version 1)", "NA", "", "10", 3)] == [1, 2, 6, 6, 6, 3]


def test_protein_change_per_event_type(cache):
    se = protein_change(cache, "SYN1:SE:1", *_event("SYN1:SE:1"))
    assert se.status == "pair" and (se.inc.transcript, se.exc.transcript) == ("SYN1-201", "SYN1-202")
    assert se.effect["kind"] == "insertion" and se.effect["short"] == "+50 aa, in frame"
    assert len(se.inc.protein) - len(se.exc.protein) == 50
    assert se.inc.event_aa == [(67, 117)] and se.exc.insert_after == 67           # 200 coding nt before the exon
    assert labels(se.only_inc) == ["LIG_SYN_1", "disordered region"] and labels(se.only_exc) == []
    assert se.in_event == ["LIG_SYN_1", "disordered region"]
    assert "at the event: LIG_SYN_1 / disordered region" in se.summary()
    a3 = protein_change(cache, "SYN1:A3SS:1", *_event("SYN1:A3SS:1"))
    assert a3.status == "pair" and a3.inc.transcript == "SYN1-203" and a3.effect["short"] == "+10 aa, in frame"
    ri = protein_change(cache, "SYN1:RI:1", *_event("SYN1:RI:1"))           # the retaining transcript is noncoding
    assert ri.status == "pair" and not ri.inc.coding and ri.inc.context == "retained_intron"
    assert ri.effect["kind"] == "frameshift" and "850 nt after residue 117 of SYN1-201" in ri.effect["text"]
    assert [i.transcript for i in ri.shown] == ["SYN1-201"]
    mx = protein_change(cache, "SYN2:MXE:1", *_event("SYN2:MXE:1"))
    assert mx.effect["kind"] == "exchange" and mx.effect["frame"] == "in frame" and len(mx.shown) == 2
    assert mx.match.pairs.context_similarity.iloc[0] == 1.0


def test_fallbacks_and_no_match(cache):
    # the cassette with other flanking exons: no exact INC transcript, so the exon alone is used
    r, g = _event("SYN1:SE:1", constant="100000-100200;104000-104300")
    pc = protein_change(cache, "x", r, g)
    assert pc.inc.context == "partial_event" and pc.inc.transcript in ("SYN1-201", "SYN1-203")
    # a novel exon: only the exclusion form is annotated; the effect is read from its coding coordinates
    r, g = _event("SYN1:SE:1", variable="103010-103150")
    pc = protein_change(cache, "x", r, g)
    assert pc.status == "exclusion_only" and "inclusion would add 140 nt after residue 67" in pc.effect["text"]
    # nothing annotated there
    r, g = _event("SYN1:SE:1", constant="110000-110100;112000-112100", variable="111000-111050")
    assert protein_change(cache, "x", r, g).status == "no_transcript"
    r, g = _event("SYN1:SE:1", gene="NOPE", gene_id="NOPE", constant="900000-900100;902000-902100",
                  variable="901000-901050")
    pc = protein_change(cache, "x", r, g)
    assert pc.status == "gene_not_in_cache" and not pc.ok and pc.row()["protein_status"] == "gene_not_in_cache"


def test_panel_draws_the_band_only_with_a_match(ds, results, gtf_path, cache, tmp_path):
    p = sa.event_panel(ds, "SYN1:SE:1", ["COH1"], "OS", gtf=gtf_path, settings=FAST, results=results, proteins=cache,
                       out_dir=tmp_path, stem="with")
    t = p.table[p.table.panel.eq("protein")]
    assert set(t.what) >= {"effect", "isoform", "event_residues", "domain", "insert_after"}
    assert set(t[t.what.eq("isoform")].isoform) == {"SYN1-201", "SYN1-202"}
    assert p.provenance["call"]["protein_status"] == {"SYN1:SE:1": "pair"}
    q = sa.event_panel(ds, "SYN1:SE:1", ["COH1"], "OS", gtf=gtf_path, settings=FAST, results=results)
    assert "protein" not in set(q.table.panel) and p.figure.get_figheight() > q.figure.get_figheight()
    none = {"SYN1:SE:1": protein_change(cache, "SYN1:SE:1", *_event("SYN1:SE:1", constant="110000-110100;"
                                                                   "112000-112100", variable="111000-111050"))}
    r = sa.event_panel(ds, "SYN1:SE:1", ["COH1"], "OS", gtf=gtf_path, settings=FAST, results=results, proteins=none)
    assert "protein" not in set(r.table.panel) and r.figure.get_figheight() == q.figure.get_figheight()


def test_probe_and_cli_with_proteins(ds, cache, cache_dir, tmp_path, capsys):
    res = probe(ds, genes=["SYN1"], settings=FAST, out_dir=tmp_path / "pr", top=1, max_pages=1, proteins=cache,
                log=lambda *_: None)
    assert {"protein_status", "protein_change", "protein_summary"} <= set(res.events.columns)
    assert res.paths["proteins"].exists() and "Protein (suggested)" in res.paths["report"].read_text()
    tab, _ = protein_changes(cache, ds, ["SYN1:SE:1", "SYN2:MXE:1"])
    assert list(tab.protein_status) == ["pair", "pair"] and tab.effect_short.iloc[1] == "exon swap, in frame"
    out = tmp_path / "ex"
    assert main(["example", str(out)]) == 0
    capsys.readouterr()
    assert main(["proteins", str(out / "data"), "--proteins", str(out / "proteins"), "--gene", "SYN1",
                 "--out", str(tmp_path / "p.csv")]) == 0
    text = capsys.readouterr().out
    assert "SYN1:SE:1: inclusion adds 50 aa in frame" in text and (tmp_path / "p.csv").exists()
    assert main(["proteins", str(out / "data"), "--event", "SYN1:SE:1"]) == 2          # no cache given
    assert main(["panel", str(out / "data"), "--event", "SYN2:MXE:1", "--cohort", "COH3", "--no-detail",
                 "--proteins", str(cache_dir), "--out", str(tmp_path / "pn"), "--settings", str(_fast(tmp_path))]) == 0
    assert "protein SYN2:MXE:1: drawn" in capsys.readouterr().out


def test_build_cache_needs_r_and_the_rds_files(tmp_path, monkeypatch):
    from splice_assay.protein import cache_prefix
    if shutil.which("Rscript"):
        with pytest.raises(InputError, match="no <name>.gtf.rds with its <name>_sequences.rds"):
            build_cache(tmp_path, tmp_path / "out")
    for name in ("mouse_gencode_vM34", "human_gencode_v45"):
        (tmp_path / f"{name}.gtf.rds").touch()
    (tmp_path / "mouse_gencode_vM34_sequences.rds").touch()
    assert cache_prefix(tmp_path) == "mouse_gencode_vM34"                      # any species or release
    (tmp_path / "human_gencode_v45_sequences.rds").touch()
    with pytest.raises(InputError, match="2 annotations"):
        cache_prefix(tmp_path)
    monkeypatch.setattr("shutil.which", lambda *_: None)
    with pytest.raises(InputError, match="Rscript was not found"):
        build_cache(tmp_path, tmp_path / "out")


def _fast(tmp_path):
    p = tmp_path / "fast.json"
    p.write_text('{"formats": ["png"], "dpi": 60}')
    return p
