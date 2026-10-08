import numpy as np
import pandas as pd
import pytest
from scipy.stats import spearmanr

import splice_assay as sa
from splice_assay.cli import main
from splice_assay.correlation import CORR_FIRST, correlation_figure, correlations, strong
from splice_assay.probe import overview, probe

FAST = sa.Settings(formats=("png",), dpi=60)
SYN1 = ["SYN1:SE:1", "SYN1:A3SS:1", "SYN1:RI:1"]


def _texts(fig) -> list[str]:
    return [t.get_text() for t in fig.texts] + [t.get_text() for a in fig.axes for t in a.texts]


def test_spearman_in_each_cohorts_survival_samples(ds):
    """rho over the survival samples (one case sample per patient) with both values; pairs only within a gene."""
    c = correlations(ds, events=SYN1 + ["SYN2:MXE:1"])
    assert list(c.columns[:len(CORR_FIRST)]) == CORR_FIRST and set(c.kind) == {"expression", "event"}
    assert len(c[c.kind.eq("expression")]) == 4 * 4 and len(c[c.kind.eq("event")]) == 3 * 4   # SYN2: one event
    assert set(c[c.kind.eq("event")].gene) == {"SYN1"}
    smp = ds.samples.sample_id[ds.samples.cohort.eq("COH1") & ds.samples.survival_cohort]
    x, y = ds.psi.loc["SYN1:SE:1"].reindex(smp), ds.expression.loc["SYN1"].reindex(smp)
    ok = x.notna() & y.notna()
    rho, p = spearmanr(x[ok], y[ok])
    row = c[c.kind.eq("expression") & c.event_id.eq("SYN1:SE:1") & c.cohort.eq("COH1")].iloc[0]
    assert row.partner == "SYN1" and row.corr_status == "tested" and row.corr_n == ok.sum()
    assert row.rho == pytest.approx(rho, rel=1e-12) and row.corr_p == pytest.approx(p, rel=1e-12)
    pair = c[c.kind.eq("event") & c.event_id.eq("SYN1:SE:1") & c.partner.eq("SYN1:A3SS:1") & c.cohort.eq("COH2")]
    assert len(pair) == 1 and pair.corr_status.iloc[0] == "tested"
    fam = c.groupby(["gene", "kind"]).agg(n=("corr_q_tests", "max"), q=("corr_q", "count"))   # BH per gene and kind
    assert fam.loc[("SYN1", "expression")].tolist() == [12, 12] and fam.loc[("SYN1", "event")].tolist() == [12, 12]
    assert fam.loc[("SYN2", "expression")].tolist() == [4, 0]              # a family under fdr_min_family: no q


def test_gates_and_missing_expression(ds, tables):
    few = correlations(ds, events=SYN1, settings=sa.Settings(corr_min_n=10_000))
    assert set(few.corr_status) == {"too_few_patients"} and few.rho.isna().all()
    no_expr = sa.Dataset.from_tables(**{k: v for k, v in tables.items() if k != "expression"})
    assert set(correlations(no_expr, events=SYN1).kind) == {"event"}
    only = correlations(ds, events=SYN1, within_gene=False)
    assert set(only.kind) == {"expression"}
    ev = tables["events"].copy()
    ev["expression_gene"] = np.where(ev.event_id.eq("SYN1:SE:1"), "NOT_MEASURED", ev.gene)
    other = sa.Dataset.from_tables(**dict(tables, events=ev))
    c = correlations(other, events=["SYN1:SE:1"])
    assert set(c.corr_status) == {"no_expression"} and set(c.partner) == {"NOT_MEASURED"}
    with pytest.raises(ValueError):
        sa.Settings(corr_note_above=1.5)
    assert sa.Settings(corr_min_n=10).changed() == {"corr_min_n": (10, 20)}


def test_strong_pairs_and_the_figure(ds):
    """Two alternative first exons of one gene sum to about 1: rho near -1, strong in every cohort."""
    afe = ["SYN3:AFE:0001", "SYN3:AFE:0002"]
    c = correlations(ds, events=afe)
    st = strong(c, FAST)
    assert len(st) == 4 and (st.kind == "event").all() and (st.rho < -0.9).all()
    assert strong(c, FAST.replace(corr_note_above=0)).empty
    ranked = pd.DataFrame({"rank": [1, 2], "event_id": afe, "label": ["AFE:0001", "AFE:0002"], "gene": "SYN3",
                           "event_type": "AFE"})
    fig = correlation_figure(c, ranked, FAST)
    texts = _texts(fig)
    assert {"SYN3", "with SYN3 expression", "between its events", "AFE:0001", "AFE:0002", "×"} <= set(texts)
    values = [t.get_text() for t in fig.axes[0].texts if t.get_text().startswith("−.9")]
    assert len(values) == 4                                                # the pair in each cohort, printed
    assert {"p < 0.05", "p ≥ 0.05", "|ρ| ≥ 0.7"} <= set(texts)                # the key, not the caption
    assert not any("value black" in t for t in texts)
    assert correlation_figure(c.iloc[:0], ranked, FAST) is None
    two = c[c.cohort.isin(sorted(c.cohort.unique())[:2])]                       # a narrow grid, long labels
    long = ranked.assign(label=["AFE:0001 the first exon of a much longer form", "AFE:0002 the other first exon"])
    narrow = correlation_figure(two, long, FAST)
    assert all(a.get_position().x1 <= 1 + 1e-9 for a in narrow.axes)            # the colour bar within the figure


def test_probe_with_correlation(ds, tmp_path):
    res = probe(ds, genes=["SYN1"], settings=FAST, out_dir=tmp_path, max_pages=0, correlation=True,
                log=lambda *_: None)
    for k in ("correlations", "correlation"):
        assert res.paths[k].exists()
    assert {"expr_rho", "expr_rho_p", "expr_rho_q", "expr_rho_n"} <= set(res.cells.columns)
    assert "best_expr_rho" in pd.read_csv(res.paths["events"]).columns
    c = pd.read_csv(res.paths["correlations"])
    one = res.cells.set_index(["event_id", "cohort"]).loc[("SYN1:SE:1", "COH1")]
    want = c[c.kind.eq("expression") & c.event_id.eq("SYN1:SE:1") & c.cohort.eq("COH1")].rho.iloc[0]
    assert one.expr_rho == pytest.approx(want)
    rep = res.paths["report"].read_text()
    assert "## Correlation (Spearman, within each cohort)" in rep and "| `correlations.csv` |" in rep
    assert "Alternative first (last) exons" not in rep                    # SYN1 has none
    plain = probe(ds, genes=["SYN1"], settings=FAST, out_dir=tmp_path / "plain", max_pages=0, log=lambda *_: None)
    assert "correlations" not in plain.paths and "expr_rho" not in plain.cells and plain.correlations is None
    assert "## Correlation" not in plain.paths["report"].read_text()
    assert not (tmp_path / "plain" / "correlation.png").exists()


def test_the_overview_shows_host_expression_in_its_own_row(ds, tmp_path):
    res = probe(ds, genes=["SYN1"], settings=FAST, out_dir=tmp_path, max_pages=0, log=lambda *_: None)
    gex = pd.read_csv(res.paths["expression"])
    fig = overview(res.cells, res.events, "OS", FAST, gex_cells=gex, gex_model="expression + age")
    labels = [t.get_text() for t in fig.axes[0].get_yticklabels()]
    assert labels[-1] == "SYN1 expression" and labels[:-1] == [f"{r.rank}. {r.label}  {r.gene}"
                                                                 for r in res.events.itertuples()]
    cap = " ".join(fig.texts[i].get_text() for i in range(1, len(fig.texts)))     # the caption lines, in order
    assert "last row: host-gene expression" in cap and "(Cox on expression + age)" in cap
    assert "|Δ median| > 1" in cap
    plain = overview(res.cells, res.events, "OS", FAST)
    assert "SYN1 expression" not in [t.get_text() for t in plain.axes[0].get_yticklabels()]
    assert not any(t.startswith("last row") for t in _texts(plain))
    assert "the host gene's expression in its own row" in res.paths["report"].read_text()
    off = probe(ds, genes=["SYN1"], settings=FAST, out_dir=tmp_path / "off", max_pages=0, gex=False,
                log=lambda *_: None)
    assert "in its own row" not in off.paths["report"].read_text()
    # an event whose host gene is named differently: the hosts map finds its row
    renamed = gex.assign(gene="SYN1_HOST")
    fig = overview(res.cells, res.events, "OS", FAST, gex_cells=renamed,
                   hosts={e: "SYN1_HOST" for e in res.events.event_id})
    assert [t.get_text() for t in fig.axes[0].get_yticklabels()][-1] == "SYN1_HOST expression"


def test_rarely_observed_events_stay_out_of_the_maps_and_correlations(ds, tmp_path):
    from dataclasses import replace
    psi = ds.psi.copy()
    psi.loc["SYN1:RI:1", psi.columns[np.arange(psi.shape[1]) % 5 != 0]] = np.nan    # observed in 1 sample of 5
    res = probe(replace(ds, psi=psi), genes=["SYN1"], settings=FAST, out_dir=tmp_path, max_pages=0,
                correlation=True, log=lambda *_: None)
    assert "SYN1:RI:1" in set(res.events.event_id)                                # still ranked (not measurable)
    assert "SYN1:RI:1" not in set(res.correlations.event_id) | set(res.correlations.partner)
    rep = res.paths["report"].read_text()
    assert "1 is observed in under 50% of every cohort's survival samples and left out of the gene maps and " \
        "correlations (`--min-observed`)." in rep
    psi.loc[:, psi.columns[np.arange(psi.shape[1]) % 5 != 0]] = np.nan              # every event rarely observed
    none = probe(replace(ds, psi=psi), genes=["SYN1"], settings=FAST, out_dir=tmp_path / "none", max_pages=0,
                 correlation=True, log=lambda *_: None)
    assert len(none.correlations) == 0 and "expr_rho" in pd.read_csv(none.paths["cells"]).columns
    assert "every event is observed in under 50% of every cohort's survival samples" in \
        none.paths["report"].read_text()
    from splice_assay.cli import _settings
    import argparse
    assert _settings(argparse.Namespace(min_observed=0.2)).min_observed == 0.2
    with pytest.raises(ValueError, match="min_observed"):
        sa.Settings(min_observed=1.5)


def test_cli_probe_correlation(tmp_path, capsys):
    out = tmp_path / "ex"
    assert main(["example", str(out)]) == 0
    (tmp_path / "fast.json").write_text('{"formats": ["png"], "dpi": 60}')
    assert main(["probe", str(out / "data"), "--gene", "SYN1", "--correlation", "--no-adjust", "--max-pages", "0",
                 "--out", str(tmp_path / "pr"), "--settings", str(tmp_path / "fast.json")]) == 0
    text = capsys.readouterr().out
    assert "correlation: 12 event x cohort cells with host expression, 12 event pair x cohort cells tested" in text
    assert (tmp_path / "pr" / "correlations.csv").exists() and (tmp_path / "pr" / "correlation.png").exists()
    c = pd.read_csv(tmp_path / "pr" / "correlations.csv")
    assert np.isfinite(c.rho).all()


def test_hit_pairs_form_their_own_family_and_events_are_named_once(ds):
    """Pairs with a HIT index get q families of their own, so including them changes no PSI q; an event named twice is
    correlated (and probed) once."""
    psi_only = ["SYN3:AFE:0001", "SYN3:AFE:0002", "SYN3:ALE:0001", "SYN3:ALE:0002"]
    s = FAST.replace(fdr_min_family=2)
    a = correlations(ds, events=psi_only, settings=s)
    b = correlations(ds, events=psi_only + ["SYN3:HIT:0001", "SYN3:HIT:0002"], settings=s)
    key = ["kind", "event_id", "partner", "cohort"]
    both = a.merge(b, on=key, suffixes=("", "_hit"))
    assert len(both) == len(a) and np.allclose(both.corr_q, both.corr_q_hit, equal_nan=True)
    hit = b.event_id.str.contains("HIT") | b.partner.str.contains("HIT")
    sizes = b.assign(hit=hit).groupby(["kind", "hit"]).corr_q_tests.unique().map(list).to_dict()
    assert sizes == {("event", False): [24], ("event", True): [36], ("expression", False): [16],
                     ("expression", True): [8]}                       # 6 PSI pairs, 9 with a HIT index; x 4 cohorts
    twice = correlations(ds, events=["SYN1:SE:1", "SYN1:SE:1"])
    assert len(twice) == 4 and set(twice.kind) == {"expression"}
    res = probe(ds, events=["SYN1:SE:1", "SYN1:SE:1"], settings=FAST, adjusted=None, correlation=True,
                log=lambda *_: None)
    assert len(res.cells) == 4 and len(res.correlations) == 4


def test_the_report_on_alternative_first_and_last_exons(ds, tmp_path):
    """Two AFEs (ALEs) of a gene: ρ near −1 by construction, said so; strong pairs listed per pair, not compared with
    chance; without the expression page the report does not point at it."""
    res = probe(ds, genes=["SYN3"], settings=FAST, out_dir=tmp_path, max_pages=0, correlation=True, gex=False,
                log=lambda *_: None)
    rep = res.paths["report"].read_text()
    assert "Two such exons sum to 1, so their ρ is near −1 by construction" in rep
    assert "AFE:0001–AFE:0002 in 4 of 4 cohorts (ρ −" in rep and "not compared with chance" in rep
    assert "expression page" not in rep.split("## Correlation")[1].split("## Files")[0]


def test_empty_correlation_and_settings(ds, tables, tmp_path):
    no_expr = sa.Dataset.from_tables(**{k: v for k, v in tables.items() if k != "expression"})
    res = probe(no_expr, events=["SYN2:MXE:1"], settings=FAST, adjusted=None, out_dir=tmp_path, max_pages=0,
                correlation=True, log=lambda *_: None)
    assert res.correlations.empty and list(res.correlations.columns) == CORR_FIRST
    assert "correlation" not in res.paths and res.paths["correlations"].exists()
    assert "**Nothing to correlate:**" in res.paths["report"].read_text()
    with pytest.raises(ValueError):
        sa.Settings(corr_min_n=2)
    rep = (lambda s, **kw: probe(ds, genes=["SYN1"], settings=s, out_dir=tmp_path / str(len(kw)), max_pages=0,
                                 gex=False, log=lambda *_: None, **kw).paths["report"].read_text())
    assert "relaxed gates" not in rep(FAST.replace(corr_min_n=10))         # correlation did not run
    assert "relaxed gates" in rep(FAST.replace(corr_min_n=10), correlation=True)


def test_the_overview_layout(ds, tmp_path):
    """One row per event, the gene rows after a half-row gap; the captions wrap to the grid instead of widening the
    figure; event labels in their type's colour."""
    from splice_assay.plot.style import TYPE_COLOR
    res = probe(ds, genes=["SYN1"], settings=FAST, out_dir=tmp_path, max_pages=0, log=lambda *_: None)
    plain = overview(res.cells, res.events, "OS", FAST)
    nr = int(res.events.measurable.sum())
    assert plain.axes[0].get_ylim() == (nr, 0) and plain.get_figwidth() == pytest.approx(4.5)
    gex = pd.read_csv(res.paths["expression"])
    with_row = overview(res.cells, res.events, "OS", FAST, gex_cells=gex)
    assert with_row.axes[0].get_ylim() == (nr + 1.5, 0) and with_row.get_figwidth() == pytest.approx(4.5)
    caption_lines = (lambda f: len([t for t in f.texts if t.get_fontsize() == 5.8]))      # not the key's labels
    assert caption_lines(with_row) > caption_lines(plain) >= 2                  # wrapped, not one long line
    assert with_row.get_figheight() > plain.get_figheight()
    colours = {t.get_text(): t.get_color() for t in plain.axes[0].get_yticklabels()}
    assert colours[next(k for k in colours if "SE:1" in k)] == TYPE_COLOR["SE"]


def test_the_overview_hatches_imprecise_cells(ds, tmp_path):
    """Every shown fit beyond the CI ratio is hatched; none with the setting at 0."""
    from matplotlib.collections import LineCollection

    from splice_assay.plot import grid as G
    res = probe(ds, genes=["SYN1"], settings=FAST, out_dir=tmp_path, max_pages=0, gex=False, log=lambda *_: None)
    hatched = (lambda s: sum(isinstance(c, LineCollection) for c in overview(res.cells, res.events, "OS", s)
                             .axes[0].collections))
    tight = FAST.replace(imprecise_ci_ratio=1.001)                              # every CI is wider than that
    ev = set(res.events[res.events.measurable].event_id)
    shown = [G.shown_fit(x, tight) for _, x in res.cells[res.cells.event_id.isin(ev)].iterrows()]
    assert hatched(tight) == sum(f.imprecise for f in shown) == sum(f.tested for f in shown) > 0
    assert hatched(FAST.replace(imprecise_ci_ratio=0)) == 0


def _pairs(fig) -> set[tuple[str, str]]:
    """The (left, right) labels of a correlation figure's pair rows: the texts level with each "×"."""
    at = {}
    for t in fig.texts:
        at.setdefault(round(t.get_position()[1], 6), []).append(t)
    return {tuple(t.get_text() for t in sorted(row, key=lambda t: t.get_position()[0]) if t.get_text() != "×")
            for row in at.values() if any(t.get_text() == "×" for t in row)}


def test_the_correlation_grid_shares_its_rows(ds):
    """With too few rows for everything: the best-ranked events' expression rows take at most half while pairs wait,
    and each gene's strongest pairs the rest, the genes taking turns; fewer pairs leave the expression rows more."""
    syn3 = ["SYN3:AFE:0001", "SYN3:AFE:0002", "SYN3:ALE:0001", "SYN3:ALE:0002"]
    c = correlations(ds, events=SYN1 + syn3)
    ranked = pd.DataFrame({"rank": range(1, 8), "event_id": syn3 + SYN1, "gene": ["SYN3"] * 4 + ["SYN1"] * 3,
                           "label": [e.split(":", 1)[1] for e in syn3 + SYN1],
                           "event_type": [e.split(":")[1] for e in syn3 + SYN1]})
    rho = c[c.kind.eq("event")].groupby(["event_id", "partner"]).rho.apply(lambda r: r.abs().max())
    lab = dict(zip(ranked.event_id, ranked.label))
    strongest = {g: [tuple(lab[e] for e in ab) for ab in rho[[a.startswith(g) for a, _ in rho.index]]
                     .sort_values(ascending=False).index] for g in ("SYN3", "SYN1")}
    fig = correlation_figure(c, ranked, FAST, max_rows=9)        # 7 expression rows, 9 pairs: 4 and 5 shown
    texts = _texts(fig)
    assert "with SYN3 expression" in texts and "with SYN1 expression" not in texts
    want = strongest["SYN3"][:3] + strongest["SYN1"][:2]                          # the genes take turns
    assert {tuple(sorted(p)) for p in _pairs(fig)} == {tuple(sorted(p)) for p in want}
    assert {("AFE:0001", "AFE:0002"), ("ALE:0001", "ALE:0002")} <= _pairs(fig)   # the two pairs near -1
    assert any("9 of 16 rows: the best-ranked events' expression rows, then each gene's strongest pairs" in t
               for t in texts)
    few = c[c.kind.eq("expression") | c.event_id.isin(SYN1)]                      # 7 expression rows, 3 pairs
    texts = _texts(fig := correlation_figure(few, ranked, FAST, max_rows=8))      # 5 of the 7, and every pair
    assert {"with SYN3 expression", "with SYN1 expression"} <= set(texts) and len(_pairs(fig)) == 3
    assert any("8 of 10 rows: the best-ranked events' expression rows, then every pair" in t for t in texts)
    mixed = ranked.assign(rank=[1, 5, 6, 7, 2, 3, 4])                            # SYN3's events 1 and 5 to 7
    texts = _texts(correlation_figure(few, mixed, FAST, max_rows=8))
    assert {"1", "2", "3", "4", "5", "AFE:0001", "AFE:0002"} <= set(texts)      # ranks 1 to 5
    assert not {"6", "7", "ALE:0001", "ALE:0002"} & set(texts)                   # 6 and 7 left out
    texts = _texts(fig := correlation_figure(c[c.event_id.isin(syn3)], ranked, FAST, max_rows=8))   # 4 and 6
    assert len(_pairs(fig)) == 4
    assert any("8 of 10 rows: every expression row, then each gene's strongest pairs" in t for t in texts)
    texts = _texts(correlation_figure(c[c.kind.eq("event")], ranked, FAST, max_rows=5))   # no expression rows
    assert any("5 of 9 rows: each gene's strongest pairs (all" in t for t in texts)
    assert not any(" rows: " in t for t in _texts(correlation_figure(c, ranked, FAST)))   # all 16 fit
    assert correlation_figure(c, ranked, FAST, max_rows=0) is None                # nothing to draw
    assert correlation_figure(c, ranked.assign(event_id="other"), FAST) is None