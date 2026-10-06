from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

import splice_assay as sa
from splice_assay.plot.genemap import gene_map, map_events, squeeze
from splice_assay.probe import combine, probe, rank_events

FAST = sa.Settings(formats=("png",), dpi=60)


def _texts(fig) -> list[str]:
    return [t.get_text() for t in fig.texts] + [t.get_text() for a in fig.axes for t in a.texts]


@pytest.fixture(scope="module")
def probed(ds):
    base = sa.analyze(ds, events=["SYN1:SE:1", "SYN1:A3SS:1", "SYN1:RI:1"], endpoints=["OS"])
    cells = combine(base, None, "OS")
    return cells, rank_events(cells, ds, 0.05)


def test_squeeze_shortens_long_introns_and_the_longest_exons():
    knots, drawn, note = squeeze([(0, 100), (5100, 5200)], 0, 5200)
    assert note == "long introns shortened, exons to scale" and list(knots) == [0, 100, 5100, 5200]
    assert np.diff(drawn)[0] == 100 and np.diff(drawn)[2] == 100                 # exons to scale
    assert np.diff(drawn)[1] == pytest.approx(200)                                # the intron: the exons' length
    knots, drawn, note = squeeze([(0, 100), (150, 300)], 0, 300)                  # a short intron stays
    assert note == "" and list(drawn) == [0, 100, 150, 300]
    k, d, _ = squeeze([(0, 50), (1000, 1100), (1200, 1210), (9000, 9100)], 0, 9100)
    gaps = np.diff(d)[[1, 3, 5]]
    assert np.all(np.diff(d) > 0) and gaps[2] > gaps[0] > gaps[1]                 # monotone, longer stays longer
    exons = [(i * 1000, i * 1000 + 100) for i in range(8)] + [(8000, 12000)]      # one exon 40x the others
    k, d, note = squeeze(exons, 0, 12000)
    seg = dict(zip(zip(k[:-1], k[1:]), np.diff(d)))
    assert note == "long introns and the longest exons shortened"
    assert seg[(0.0, 100.0)] == pytest.approx(100) and 300 < seg[(8000.0, 12000.0)] < 4000 / 4
    micro = [(i * 1000, i * 1000 + 15) for i in range(5)] + [(5000, 5150)]       # a median of 15 nt
    k, d, _ = squeeze(micro, 0, 5150)
    assert dict(zip(zip(k[:-1], k[1:]), np.diff(d)))[(5000.0, 5150.0)] == pytest.approx(150)   # kept to scale


def test_map_events_genomic_order_and_what_is_left_out(ds, probed):
    _, ranked = probed
    ids, left = map_events(ds, "SYN1", ranked)
    span = {e: sa.geometry(e, ds.events.loc[e]).span for e in ids}
    assert sorted(ids, key=lambda e: span[e]) == ids and not any(left.values())   # plus strand: 5' first
    best2, left = map_events(ds, "SYN1", ranked, max_events=2)
    assert left["beyond"] == 1 and set(best2) == set(ranked.event_id[:2])
    flip = ds.events.assign(strand="-")                                           # the same events, minus strand
    assert map_events(replace(ds, events=flip), "SYN1", ranked)[0] == ids[::-1]
    moved = ds.events.copy()
    moved.loc["SYN1:RI:1", "chrom"] = "chr8"                                      # another chromosome
    moved.loc["SYN1:A3SS:1", "variable"] = "102010-102030"                        # does not extend an exon
    ids, left = map_events(replace(ds, events=moved), "SYN1", ranked)
    assert ids == ["SYN1:SE:1"] and left["elsewhere"] == 1 and left["unreadable"] == 1


def test_gene_map_rows_and_layers(ds, probed, gtf_path):
    cells, ranked = probed
    fig = gene_map(ds, "SYN1", ranked, cells, "OS", FAST, gtf=str(gtf_path))
    texts = _texts(fig)
    assert {"SYN1", "gene map · 3 events · OS", "SE:1", "A3SS:1", "RI:1", "PSI = inclusion", "SNORD900",
            "long introns shortened, exons to scale", "Cox HR per SD · OS"} <= set(texts)
    y = {t.get_text(): t.get_position()[1] for t in fig.texts}
    order = map_events(ds, "SYN1", ranked)[0]
    labels = [ds.events.at[e, "label"] for e in order]
    assert sorted(labels, key=lambda lab: -y[lab]) == labels                      # rows top to bottom, 5' to 3'
    assert "ρ between events" not in texts                                        # no correlations given
    from splice_assay.correlation import correlations
    c = correlations(ds, events=order)
    with_c = _texts(gene_map(ds, "SYN1", ranked, cells, "OS", FAST, gtf=str(gtf_path), corr=c))
    assert {"ρ between events", "Spearman ρ, median of cohorts"} <= set(with_c) and "expr." not in with_c
    pairs = c[c.kind.eq("event") & c.corr_status.eq("tested")]
    from splice_assay.plot.grid import rho_text
    for med in pairs.groupby(["event_id", "partner"]).rho.median():
        assert rho_text(med) in with_c                                            # each pair's median, printed
    bare = _texts(gene_map(ds, "SYN1", ranked, cells, "OS", FAST))                # no GTF, no expression
    assert texts.count("SYN1") == 2 and bare.count("SYN1") == 1                   # title, then the gene's row
    assert any("no gene model" in t for t in bare)


def test_rarely_observed_events_are_left_out(ds, probed):
    cells, ranked = probed
    from splice_assay.probe import observed
    sparse = cells.assign(frac_obs=cells.frac_obs.where(cells.event_id.ne("SYN1:RI:1"), 0.2))
    keep = observed(sparse, FAST)
    assert "SYN1:RI:1" not in keep and "SYN1:SE:1" in keep
    assert "SYN1:RI:1" in observed(sparse, FAST.replace(min_observed=0.2))
    assert observed(sparse, FAST.replace(coverage_frac=0.2)) >= {"SYN1:RI:1"}          # the default follows the gate
    ids, left = map_events(ds, "SYN1", ranked, keep=keep)
    assert "SYN1:RI:1" not in ids and left["rarely observed"] == 1
    fig = gene_map(ds, "SYN1", ranked, sparse, "OS", FAST, keep=keep)
    assert any("1 event not drawn (observed in under 50% of every cohort's survival samples)" in x
               for x in _texts(fig))


def test_a_long_gene_is_drawn_in_a_window(ds, probed, tmp_path):
    """A gene far longer than its events whose exons all lie outside the events' window (no exon in view)."""
    cells, ranked = probed
    from splice_assay import example
    lines = [x for x in example.gtf_text().splitlines() if "SYN1" not in x or "SNORD" in x]
    attrs = 'gene_id "SYNG0001.1"; gene_type "lncRNA"; gene_name "SYN1";'
    lines.append(f"chr7\tS\tgene\t1\t400000\t.\t+\t.\t{attrs}")
    for t in range(3):                                                           # exons far from the events
        ta = attrs + f' transcript_id "T{t}"; transcript_type "lncRNA";'
        lines.append(f"chr7\tS\ttranscript\t1\t400000\t.\t+\t.\t{ta}")
        for a in (1, 300001, 399001):
            lines.append(f"chr7\tS\texon\t{a}\t{a + 299}\t.\t+\t.\t{ta}")
    gtf = tmp_path / "long.gtf"
    gtf.write_text("\n".join(lines) + "\n")
    fig = gene_map(ds, "SYN1", ranked, cells, "OS", FAST, gtf=str(gtf))
    assert any(" of a 399 kb gene" in t for t in _texts(fig))


def test_probe_writes_a_map_per_gene_and_survives_a_bad_event(ds, tables, gtf_path, tmp_path):
    res = probe(ds, genes=["SYN1", "SYN2"], settings=FAST, gtf=gtf_path, out_dir=tmp_path, max_pages=0,
                log=lambda *_: None)
    best = res.events[res.events.measurable]
    assert [p.name for p in res.paths["gene_maps"]] == [f"gene_map_{g}.png" for g in dict.fromkeys(best.gene)]
    assert all(p.exists() for p in res.paths["gene_maps"])
    rep = res.paths["report"].read_text()
    assert "| `gene_map_<GENE>.png` (SYN1, SYN2) | Each gene's model" in rep and "the gene maps, " in rep
    ev = tables["events"].copy()                                                   # an RI that is not the gap
    ev.loc[ev.event_id.eq("SYN1:RI:1"), "variable"] = "103150-103900"
    bad = sa.Dataset.from_tables(**dict(tables, events=ev))
    out = probe(bad, genes=["SYN1"], settings=FAST, out_dir=tmp_path / "bad", max_pages=0, log=lambda *_: None)
    assert out.paths["gene_maps"].exists()                                         # drawn without the bad event
