import numpy as np
import pytest

import splice_assay as sa
from splice_assay import InputError

FAST = sa.Settings(formats=("svg", "png"), dpi=100)


def test_single_event_panel(ds, results, gtf_path, tmp_path):
    p = sa.event_panel(ds, "SYN1:SE:1", ["COH1", "COH2"], "OS", gtf=gtf_path, out_dir=tmp_path, settings=FAST,
                       results=results)
    assert set(p.paths) == {"svg", "png", "csv", "json"} and all(v.stat().st_size > 0 for v in p.paths.values())
    t = p.table
    assert {"pairs", "all_samples", "km_curve", "km_at_risk", "forest", "schematic_gene_block"} <= set(t.panel)
    # the forest shows exactly the analysis numbers
    f = t[t.panel.eq("forest")].set_index("cohort")
    tn = results.groups.set_index(["event_id", "cohort"])
    sv = results.survival.set_index(["event_id", "cohort", "endpoint"])
    for c in f.index:
        assert f.at[c, "unpaired_delta_median"] == pytest.approx(tn.loc[("SYN1:SE:1", c), "unpaired_delta_median"],
                                                                 nan_ok=True)
        assert f.at[c, "hr_per_iqr"] == pytest.approx(sv.loc[("SYN1:SE:1", c, "OS"), "hr_per_iqr"], nan_ok=True)
    # COH1: pairs and all samples; COH2 has 6 pairs: all samples only
    assert set(t[t.panel.eq("all_samples")].tag) == {"SYN1:SE:1|COH1", "SYN1:SE:1|COH2"}
    assert set(t[t.panel.eq("pairs")].tag) == {"SYN1:SE:1|COH1"}
    assert p.provenance["inputs_sha256"]["psi"] and p.provenance["settings"]["min_pairs"] == 10


def test_two_events_with_highlight(ds, results, gtf_path):
    p = sa.event_panel(ds, ["SYN1:SE:1", "SYN1:RI:1"], "COH1", "DSS", gtf=gtf_path, settings=FAST, results=results,
                       highlight={("SYN1:SE:1", "COH1"): "A", ("SYN1:RI:1", "COH1"): "B"}, highlight_title="Tier")
    f = p.table[p.table.panel.eq("forest")]
    assert set(f[f.cohort.eq("COH1")].label) == {"A", "B"}
    texts = [t.get_text() for t in p.figure.texts]
    assert "SE:1 · COH1 · Tier A" in texts and "Tier cell (letter: tier)" in texts


def test_minus_strand_is_drawn_five_to_three(ds, results, gtf_path):
    p = sa.event_panel(ds, "SYN2:MXE:1", "COH3", "OS", gtf=gtf_path, settings=FAST, results=results)
    ax = p.figure.axes[0]
    lo, hi = ax.get_xlim()
    assert lo > hi
    assert "Other exon" in [t.get_text() for t in p.figure.texts]


def test_no_gtf_and_no_normals(ds, results):
    p = sa.event_panel(ds, "SYN1:A3SS:1", "COH4", "OS", settings=FAST, results=results)
    assert "groups_not_tested" in set(p.table.panel)
    assert "schematic_gene_block" not in set(p.table.panel)


def test_outputs_are_byte_identical(ds, results, gtf_path, tmp_path):
    s = sa.Settings(formats=("svg", "pdf", "png"), dpi=100)
    a = sa.event_panel(ds, "SYN1:SE:1", "COH1", "OS", gtf=gtf_path, out_dir=tmp_path / "a", settings=s,
                       results=results)
    b = sa.event_panel(ds, "SYN1:SE:1", "COH1", "OS", gtf=gtf_path, out_dir=tmp_path / "b", settings=s)
    for k in ("svg", "pdf", "png", "csv", "json"):
        assert a.paths[k].read_bytes() == b.paths[k].read_bytes(), k


def test_panel_errors(ds):
    with pytest.raises(InputError, match="one or two"):
        sa.event_panel(ds, ["SYN1:SE:1", "SYN1:RI:1", "SYN1:A3SS:1"], "COH1", "OS")
    with pytest.raises(InputError, match="share gene"):
        sa.event_panel(ds, ["SYN1:SE:1", "SYN2:MXE:1"], "COH1", "OS")
    with pytest.raises(InputError, match="endpoint"):
        sa.event_panel(ds, "SYN1:SE:1", "COH1", "PFI")
    with pytest.raises(InputError, match="cohort"):
        sa.event_panel(ds, "SYN1:SE:1", "COH9", "OS")


def test_km_panel_values(ds, results):
    p = sa.event_panel(ds, "SYN1:SE:1", "COH1", "OS", settings=FAST, results=results)
    arms = p.table[p.table.panel.eq("km_arms")].iloc[0]
    sv = results.survival.set_index(["event_id", "cohort", "endpoint"]).loc[("SYN1:SE:1", "COH1", "OS")]
    assert arms.n_high == sv.n_high and arms.events_low == sv.events_low
    first = p.table[p.table.panel.eq("km_at_risk") & p.table.time_years.eq(0)]
    assert first.n_at_risk.sum() == sv.km_n
    assert np.all(np.diff(p.table[p.table.panel.eq("km_curve") & p.table.arm.eq("Low PSI")].surv) <= 1e-12)


def test_cox_model_figure(ds, tmp_path):
    p = sa.cox_model_figure(ds, "SYN1:SE:1", "COH1", "OS", model=sa.CoxModel(covariates=("age", "sex", "stage")),
                            settings=FAST, out_dir=tmp_path)
    terms = p.table[p.table.panel.eq("term")]
    assert list(terms.kind[:3]) == ["psi_iqr", "expression", "numeric"] and set(terms.term) >= {"sex", "stage"}
    assert p.paths["png"].exists()
    with pytest.raises(InputError, match="not fitted"):
        sa.cox_model_figure(ds, "SYN1:SE:1", "COH4", "DSS", settings=sa.Settings(cox_min_events=10_000))


def test_data_without_reference_samples(tables, gtf_path):
    s = tables["samples"]
    only = dict(tables, samples=s[s.group.eq("tumour")],
                psi=tables["psi"][tables["psi"].sample_id.isin(s.sample_id[s.group.eq("tumour")])])
    ds = sa.Dataset.from_tables(**only)
    assert not ds.has_reference
    p = sa.event_panel(ds, "SYN1:SE:1", "COH1", "OS", gtf=gtf_path, settings=FAST, gex=False)
    assert "pairs" not in set(p.table.panel) and "all_samples" not in set(p.table.panel)
    assert len(p.figure.axes) == 3                                   # schematic, KM, HR forest


def test_custom_group_labels_reach_the_figure(tables):
    s = tables["samples"].assign(group=tables["samples"].group.map({"tumour": "Metastasis", "normal": "Primary"}))
    ds = sa.Dataset.from_tables(**dict(tables, samples=s), case="Metastasis", reference="Primary")
    p = sa.event_panel(ds, "SYN1:SE:1", "COH1", "OS", settings=FAST)
    texts = {t.get_text() for t in p.figure.texts}
    assert "Higher in metastasis" in texts
    ticks = {t.get_text() for ax in p.figure.axes for t in ax.get_xticklabels()}
    assert {"Primary", "Metastasis"} <= ticks


def test_panel_with_model_band(ds, results, gtf_path, tmp_path):
    """One figure: the event panel plus the full Cox model of each cohort shown."""
    clinical = sa.CoxModel().with_clinical(("age", "sex", "stage"), baseline={"stage": "I", "sex": "female"})
    p = sa.event_panel(ds, "SYN1:SE:1", ["COH1", "COH3"], "OS", gtf=gtf_path, results=results, detail=clinical,
                       settings=FAST, out_dir=tmp_path)
    band = p.table[p.table.panel.eq("cox_detail")]
    assert set(band.tag) == {"SYN1:SE:1|COH1", "SYN1:SE:1|COH3"}
    assert set(band[band.term.eq("stage")].reference) == {"I"} and set(band[band.term.eq("sex")].reference) == {"female"}
    # the band's numbers are the stand-alone model's numbers
    alone = sa.cox_model_figure(ds, "SYN1:SE:1", "COH1", "OS", model=clinical, settings=FAST).table
    a = band[band.tag.eq("SYN1:SE:1|COH1")].set_index("label").hr
    b = alone[alone.panel.eq("term")].set_index("label").hr
    assert np.allclose(a.sort_index().to_numpy(), b.sort_index().to_numpy())
    # the forest keeps the panel's own model (PSI + host expression)
    assert p.table[p.table.panel.eq("forest_axis")].cox_model.iloc[0] == "Cox: PSI + host expression"
    assert p.paths["png"].exists()


def test_model_band_reports_cells_it_cannot_fit(ds, results):
    strict = sa.Settings(formats=("svg",), dpi=100, cox_min_events=10_000)
    p = sa.event_panel(ds, "SYN1:SE:1", ["COH1"], "OS", results=results, detail=True, settings=strict)
    assert "cox_detail_not_fitted" in set(p.table.panel)


def test_stacked_layout_for_three_or_more_cohorts(ds, results):
    """With models for three or more cells, each row carries its own model and the forest moves below."""
    clinical = sa.CoxModel().with_clinical(("age", "stage"), baseline={"stage": "I"})
    p = sa.event_panel(ds, "SYN1:SE:1", ["COH1", "COH2", "COH3"], "OS", results=results, detail=clinical,
                       settings=FAST, gex=False)
    assert set(p.table[p.table.panel.eq("cox_detail")].tag) == {f"SYN1:SE:1|COH{i}" for i in (1, 2, 3)}
    axes = p.figure.axes
    forest_hr = [a for a in axes if a.get_xlabel().startswith("HR per IQR")][0]
    models = [a for a in axes if a.get_xlabel() == "Hazard ratio (95% CI)"]
    assert len(models) == 3
    # the models are aligned (same x position and width) and share one axis; the forest sits below them all
    pos = [a.get_position() for a in models]
    assert len({round(q.x0, 6) for q in pos}) == 1 and len({round(q.width, 6) for q in pos}) == 1
    assert len({a.get_xlim() for a in models}) == 1
    assert forest_hr.get_position().y1 < min(q.y0 for q in pos)
    side = sa.event_panel(ds, "SYN1:SE:1", ["COH1", "COH2", "COH3"], "OS", results=results, detail=clinical,
                          settings=FAST, layout="side")
    assert side.figure.get_size_inches()[0] == p.figure.get_size_inches()[0]
    with pytest.raises(InputError, match="layout"):
        sa.event_panel(ds, "SYN1:SE:1", "COH1", "OS", results=results, layout="grid")


def test_long_model_notes_wrap_and_are_cut():
    from splice_assay.plot.panel import MODEL_LINES, _model_lines
    text = "PSI + host expression + age + sex + stage (" + "; ".join(["stage I merged into II (2 patients, 0 events)"] * 6) + ")"
    lines = _model_lines(text, 2.5, 5.6)
    assert len(lines) == MODEL_LINES and lines[-1].endswith(" …")
    assert _model_lines("PSI + host expression", 2.5, 5.6) == ["PSI + host expression"]
