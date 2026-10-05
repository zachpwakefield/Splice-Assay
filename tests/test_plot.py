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
        assert f.at[c, "hr_per_sd"] == pytest.approx(sv.loc[("SYN1:SE:1", c, "OS"), "hr_per_sd"], nan_ok=True)
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
    assert list(terms.kind[:3]) == ["psi_sd", "expression", "numeric"] and set(terms.term) >= {"sex", "stage"}
    assert p.paths["png"].exists()
    with pytest.raises(InputError, match="not fitted"):
        sa.cox_model_figure(ds, "SYN1:SE:1", "COH4", "DSS", settings=sa.Settings(cox_min_events=10_000))


def test_data_without_reference_samples(tables, gtf_path):
    s = tables["samples"]
    only = dict(tables, samples=s[s.group.eq("tumour")],
                psi=tables["psi"][tables["psi"].sample_id.isin(s.sample_id[s.group.eq("tumour")])])
    ds = sa.Dataset.from_tables(**only)
    assert not ds.has_reference
    p = sa.event_panel(ds, "SYN1:SE:1", "COH1", "OS", gtf=gtf_path, settings=FAST)
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
    # the reason quotes the settings of the fit, not the page's
    from splice_assay.plot.model import model_terms
    res = sa.analyze(ds, events=["SYN1:SE:1"], cohorts=["COH1"], endpoints=["OS"],
                     settings=sa.Settings(cox_min_events_per_term=1000))
    with pytest.raises(InputError, match=r"too few events per term \(\d+ events for 2 terms; needs 1000 per term\)"):
        model_terms(ds, "SYN1:SE:1", "COH1", "OS", results=res)


def test_expression_rows_say_why_a_model_was_not_fitted(ds):
    g = sa.expression_panel(ds, "SYN1", ["COH1"], "OS", settings=FAST.replace(cox_min_events_per_term=1000))
    msg = g.table[g.table.panel.eq("gex_cox_detail_not_fitted")].message.tolist()
    assert len(msg) == 1 and msg[0].startswith("not fitted: too few events per term (") and \
        msg[0].endswith(" events for 1 term; needs 1000 per term)")


def test_the_forest_names_the_penalty_of_the_results_it_shows(ds):
    """results= fitted with a ridge: the forest's model line says so, whatever the page's settings."""
    clinical = sa.CoxModel().with_clinical(("age",))
    res = sa.analyze(ds, events=["SYN1:SE:1"], endpoints=["OS"], settings=FAST.replace(cox_ridge="clinical"),
                     model=clinical)
    p = sa.event_panel(ds, "SYN1:SE:1", "COH1", "OS", results=res, settings=FAST)
    assert p.table[p.table.panel.eq("forest_axis")].cox_model.iloc[0] == \
        "Cox: PSI + host expression + age; ridge λ 1 (clinical terms)"


def test_stacked_layout_for_three_or_more_cohorts(ds, results):
    """With models for three or more cells, each row carries its own model and the forest moves below."""
    clinical = sa.CoxModel().with_clinical(("age", "stage"), baseline={"stage": "I"})
    p = sa.event_panel(ds, "SYN1:SE:1", ["COH1", "COH2", "COH3"], "OS", results=results, detail=clinical,
                       settings=FAST)
    assert set(p.table[p.table.panel.eq("cox_detail")].tag) == {f"SYN1:SE:1|COH{i}" for i in (1, 2, 3)}
    axes = p.figure.axes
    forest_hr = [a for a in axes if a.get_xlabel().startswith("HR per SD")][0]
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


def test_km_header_notes_non_proportional_hazards():
    from matplotlib.figure import Figure
    from splice_assay.plot import km as K
    rng = np.random.default_rng(0)
    t, e, high = rng.exponential(3, 60), (rng.uniform(size=60) < 0.7).astype(int), rng.uniform(size=60) > 0.5
    row = dict(logrank_hr=1.5, km_p=0.01, events_high=int(e[high].sum()), events_low=int(e[~high].sum()), cutoff=0.5)
    for ph, shown in ((0.003, True), (0.4, False)):
        fig = Figure(figsize=(3, 3))
        ax = fig.add_axes([0.2, 0.3, 0.7, 0.5])
        K.draw(fig, ax, t, e, high, dict(row, km_ph_p=ph), "Overall survival", [], "t", sa.Settings())
        assert any("non-proportional hazards (p 0.003)" in x.get_text() for x in ax.texts) == shown


def test_km_end_labels_stay_apart_when_the_curves_end_level():
    """Both arms reach 0 (every patient died): the two labels are still nudged apart, not drawn on each other."""
    from matplotlib.figure import Figure
    from splice_assay.plot import km as K
    t, e, high = np.arange(1.0, 9.0), np.ones(8, int), np.arange(8) % 2 == 1
    row = dict(logrank_hr=1.0, km_p=0.9, events_high=4, events_low=4, cutoff=0.5)
    fig = Figure(figsize=(3, 3))
    ax = fig.add_axes([0.2, 0.3, 0.7, 0.5])
    K.draw(fig, ax, t, e, high, row, "Overall survival", [], "t", sa.Settings())
    y = {x.get_text(): x.get_position()[1] for x in ax.texts if x.get_text() in ("Low PSI", "High PSI")}
    assert len(y) == 2 and y["High PSI"] - y["Low PSI"] >= 0.09 - 1e-12 and min(y.values()) >= 0.04 - 1e-12


def _texts(p):
    return [t.get_text() for t in p.figure.texts] + [t.get_text() for a in p.figure.axes for t in a.texts]


def test_km_header_names_the_split_and_the_title_the_event(ds, results):
    from splice_assay.plot.style import fcut
    p = sa.event_panel(ds, ["SYN1:SE:1", "SYN1:RI:1"], "COH1", "OS", settings=FAST, results=results)
    cut = results.survival.set_index(["event_id", "cohort", "endpoint"]).at[("SYN1:SE:1", "COH1", "OS"), "cutoff"]
    head = [t for t in _texts(p) if t.startswith("split at ")]
    assert len(head) == 2 and head[0].splitlines()[0] == f"split at median PSI {fcut(cut)}"
    detail = list(p.table[p.table.panel.eq("event_detail")].text)
    assert detail[0].startswith("SE:1: cassette exon chr7:103,001–103,150 (150 nt) between exons chr7:102,001–102,100")
    assert any(d.startswith("RI:1: retained intron chr7:103,151–104,000") for d in detail)
    assert all(d in _texts(p) for d in detail)                                 # drawn under the title
    q = sa.event_panel(ds, "SYN1:SE:1", "COH1", "OS", settings=FAST.replace(km_split=0.7))
    assert [t for t in _texts(q) if t.startswith("split at ")][0].splitlines()[0] == "split at PSI 0.7 (set)"
    assert q.table[q.table.panel.eq("km_arms")].n_low.iloc[0] < p.table[p.table.panel.eq("km_arms")].n_low.iloc[0]


def test_proportional_hazards_daggers(ds, tmp_path):
    """A dagger follows each test whose proportional-hazards p is below ph_note_below: the KM log-rank p, the forest
    CI and the p of each model row; the legend says what it means. None at 0."""
    m = sa.CoxModel().with_clinical(("age",))
    every = FAST.replace(ph_note_below=1.0)                      # every finite PH p is below 1
    p = sa.event_panel(ds, "SYN1:SE:1", ["COH1"], "OS", settings=every, detail=m)
    texts = _texts(p)
    assert "non-proportional hazards (p < 1)" in texts
    km = next(t for t in texts if t.startswith("split at ")).splitlines()
    assert "log-rank p" in km[1] and km[1].endswith(" †")
    assert km[-1].startswith("† non-proportional hazards (p ")
    assert "†" in [t.get_text() for a in p.figure.axes for t in a.texts]          # the forest CI
    band = p.table[p.table.panel.eq("cox_detail")]
    assert band.ph_marked.all() and sum(t.endswith(" †") for t in texts) == len(band)
    f = sa.cox_model_figure(ds, "SYN1:SE:1", "COH1", "OS", model=m, settings=every)
    assert f.table[f.table.panel.eq("term")].ph_marked.all()
    g = sa.expression_panel(ds, "SYN1", ["COH1"], "OS", settings=every, model=m)
    assert "non-proportional hazards (p < 1)" in _texts(g)
    for page in (sa.event_panel(ds, "SYN1:SE:1", ["COH1"], "OS", settings=FAST.replace(ph_note_below=0), detail=m),
                 sa.expression_panel(ds, "SYN1", ["COH1"], "OS", settings=FAST.replace(ph_note_below=0), model=m)):
        assert not any("†" in t or "non-proportional" in t for t in _texts(page))


def test_q_marks(ds, tmp_path):
    """* follows a q below q_mark_below wherever q is shown (group view, KM header, forest CI, model header and the
    tested term's p), apart from the p-based fill; the legend says so; q_mark_below = 0 turns it off."""
    from splice_assay.plot.style import Q_MARK
    from splice_assay.probe import overview, probe
    m = sa.CoxModel().with_clinical(("age",))
    s = FAST.replace(fdr_min_family=2)                         # the synthetic families are under 10
    p = sa.event_panel(ds, "SYN3:HIT:0002", ["COH1"], "OS", settings=s, detail=m)
    texts = _texts(p)
    assert Q_MARK == "*" and "q < 0.05 (Benjamini–Hochberg)" in texts
    km = next(t for t in texts if t.startswith("split at ")).splitlines()
    assert any(x.startswith("q ") and x.endswith(" *") for x in km)
    assert any(t.startswith("*") for a in p.figure.axes for t in [x.get_text() for x in a.texts])   # the forest
    assert any(" · HIT index q " in t and t.endswith(" *") for t in texts)                         # model header
    assert p.table[p.table.panel.eq("forest")].q_marked.all()
    band = p.table[p.table.panel.eq("cox_detail") & p.table.q.notna()]
    assert len(band) == 1 and band.q_marked.all() and any(t.endswith(" * †") or t.endswith(" *") for t in texts)
    g = sa.event_panel(ds, "SYN1:SE:1", ["COH1"], "OS", settings=s)          # a group view's q
    assert any("\nq " in t and t.splitlines()[-1].endswith(" *") for t in _texts(g))
    off = sa.event_panel(ds, "SYN3:HIT:0002", ["COH1"], "OS", settings=s.replace(q_mark_below=0), detail=m)
    assert not any(t == "*" or t.endswith(" *") or " * " in t or t.startswith("* ") for t in _texts(off))
    assert sa.Settings(q_mark_below=0.1).changed() == {}                     # a drawing setting
    res = probe(ds, genes=["SYN3"], settings=s, adjusted=None, out_dir=tmp_path, max_pages=1, include_hit=True,
                log=lambda *_: None)
    fig = overview(res.cells, res.events, "OS", s)
    assert any(t.get_text() == "*" for a in fig.axes for t in a.texts)
    assert "*: q < 0.05" in " ".join(t.get_text() for t in fig.texts)


@pytest.mark.parametrize("unit", ["sd", "iqr"])
def test_hr_unit_on_the_pages(ds, tmp_path, capsys, unit):
    """The HR per SD (the default) or per IQR (psi_hr_unit, --hr-unit iqr) in the forest, the model rows, the probe
    and its report."""
    from splice_assay.cli import main
    from splice_assay.probe import probe
    s = FAST if unit == "sd" else FAST.replace(psi_hr_unit="iqr")
    U, other = unit.upper(), ("iqr" if unit == "sd" else "sd")
    p = sa.event_panel(ds, "SYN1:SE:1", ["COH1"], "OS", settings=s, detail=sa.CoxModel().with_clinical(("age",)))
    assert f"HR per {U}\nof PSI (95% CI)" in [a.get_xlabel() for a in p.figure.axes]
    f = p.table[p.table.panel.eq("forest")]
    assert f"hr_per_{unit}" in f and f"hr_per_{other}" not in f and f[f"hr_per_{unit}"].notna().any()
    assert any(str(x).startswith(f"PSI (per {U}, ") for x in p.table[p.table.panel.eq("cox_detail")].label)
    res = probe(ds, genes=["SYN1"], settings=s, adjusted=None, out_dir=tmp_path / "p", max_pages=1,
                log=lambda *_: None)
    assert f"best_hr_per_{unit}" in res.events and f"| HR per {U} | p |" in res.paths["report"].read_text()
    from splice_assay import example
    example.write(tmp_path / "ex")
    assert main(["probe", str(tmp_path / "ex" / "data"), "--gene", "SYN2", "--no-adjust", "--out",
                 str(tmp_path / "cli"), "--max-pages", "1", "--settings", str(_fast_json(tmp_path))]
                + ([] if unit == "sd" else ["--hr-unit", "iqr"])) == 0
    import pandas as pd
    assert f"best_hr_per_{unit}" in pd.read_csv(tmp_path / "cli" / "events.csv").columns
    capsys.readouterr()


def _fast_json(tmp_path):
    p = tmp_path / "fast.json"
    p.write_text('{"formats": ["png"], "dpi": 60}')
    return p


def test_the_cox_figure_names_the_model_as_fitted(tables):
    """The header is the model fitted in that cohort: a covariate left out there is not named, and the HIT index is
    called so."""
    cl = tables["clinical"].copy()
    coh1 = sorted(set(tables["samples"].patient_id[tables["samples"].cohort.eq("COH1")]))
    cl.loc[cl.patient_id.isin(coh1[::2]), "stage"] = np.nan                     # stage recorded for half of COH1
    ds = sa.Dataset.from_tables(**dict(tables, clinical=cl))
    f = sa.cox_model_figure(ds, "SYN1:SE:1", "COH1", "OS", model=sa.CoxModel(covariates=("age", "stage")),
                            settings=FAST)
    assert next(t for t in _texts(f) if t.startswith("Cox: ")).startswith("Cox: PSI + host expression + age · ")
    head = f.table[f.table.panel.eq("header")].iloc[0]                          # the header's numbers in the CSV
    assert head.cox_model == "PSI + host expression + age" and head.cox_n > 0 and not head.cox_low_power
    assert f.provenance["inputs_sha256"]["expression"]                          # the model adjusts for it
    h = sa.cox_model_figure(ds, "SYN3:HIT:0002", "COH1", "OS", model=sa.CoxModel(covariates=("age",)), settings=FAST)
    assert next(t for t in _texts(h) if t.startswith("Cox: ")).startswith("Cox: HIT index + host expression + age · ")


def test_provenance_records_the_forest_model_and_the_q_family(ds, tables):
    """call.model is the forest's model (call.detail_model the model rows'); the PSI of the gene's other events, which
    sets the printed q values, is hashed."""
    m = sa.CoxModel().with_clinical(("age",))
    p = sa.event_panel(ds, "SYN1:SE:1", ["COH1"], "OS", settings=FAST, detail=m)
    assert p.provenance["call"]["model"] == sa.CoxModel().to_dict()
    assert p.provenance["call"]["detail_model"] == m.to_dict()
    assert p.provenance["input_rows"]["psi_q_family"] == 2                      # SYN1:A3SS:1 and SYN1:RI:1
    psi = tables["psi"].copy()
    other = psi.event_id.eq("SYN1:RI:1")
    psi.loc[other, "psi"] = psi.loc[other, "psi"].to_numpy()[::-1]
    q = sa.event_panel(sa.Dataset.from_tables(**dict(tables, psi=psi)), "SYN1:SE:1", ["COH1"], "OS", settings=FAST,
                       detail=m)
    a, b = p.provenance["inputs_sha256"], q.provenance["inputs_sha256"]
    assert a["psi"] == b["psi"] and a["psi_q_family"] != b["psi_q_family"]


def test_the_forest_table_holds_the_proportional_hazards_p(ds, results):
    p = sa.event_panel(ds, "SYN1:SE:1", ["COH1"], "OS", settings=FAST.replace(ph_note_below=1.0), results=results)
    f = p.table[p.table.panel.eq("forest")].set_index("cohort")
    sv = results.survival[results.survival.event_id.eq("SYN1:SE:1") & results.survival.endpoint.eq("OS")]
    for r in sv.itertuples():
        if r.cohort in f.index and r.cox_status == "tested" and np.isfinite(r.hr_per_sd):
            assert f.at[r.cohort, "ph_p"] == r.ph_p and f.at[r.cohort, "ph_marked"]


def test_a_highlight_row_for_an_absent_cohort(ds, results, tables):
    """A cohort not in the data is an input error (a typo), unless a subset left it out: then it is skipped."""
    marks = {("SYN1:SE:1", "COH2"): "B", ("SYN1:SE:1", "COH1"): "A"}
    with pytest.raises(InputError, match="highlight: cohort.s. coh9 not in the data"):
        sa.event_panel(ds, "SYN1:SE:1", ["COH1"], "OS", settings=FAST, results=results,
                       highlight={("SYN1:SE:1", "coh9"): "B"})
    sub = sa.Dataset.from_tables(**tables, where=["cohort=COH1"])
    with pytest.warns(UserWarning, match="COH2 not in the subset"):
        p = sa.event_panel(sub, "SYN1:SE:1", ["COH1"], "OS", settings=FAST, highlight=marks)
    with pytest.raises(InputError, match="cohort.s. coh9 not in the data"):     # a typo is one also under a subset
        sa.event_panel(sub, "SYN1:SE:1", ["COH1"], "OS", settings=FAST, highlight={("SYN1:SE:1", "coh9"): "B"})
    f = p.table[p.table.panel.eq("forest")]
    assert "COH2" not in set(f.cohort) and f.set_index("cohort").at["COH1", "label"] == "A"


def test_the_group_view_table_names_each_patient(ds, results):
    p = sa.event_panel(ds, "SYN1:SE:1", ["COH1"], "OS", settings=FAST, results=results)
    pts = p.table[p.table.panel.eq("all_samples")]
    smp = ds.samples[ds.samples.cohort.eq("COH1") & ds.samples.role.eq("case")]
    case = pts[pts.group.eq("case")]
    assert set(case.patient_id) <= set(smp.patient_id) and case.n_samples.eq(1).all()
    psi = ds.psi.loc["SYN1:SE:1"].reindex(smp.set_index("patient_id").loc[case.patient_id, "sample_id"]).to_numpy()
    assert np.allclose(case.psi.to_numpy(float), psi)


def test_a_model_other_than_the_results_model_is_refused(ds, results):
    """The page would name one model and draw another's numbers."""
    other = sa.CoxModel(expression=False)
    with pytest.raises(InputError, match=r"the model \(PSI\) is not the one results= was computed with \(PSI \+ host "
                                         r"expression; they differ in expression: False vs True\)"):
        sa.event_panel(ds, "SYN1:SE:1", ["COH1"], "OS", settings=FAST, results=results, model=other)
    with pytest.raises(InputError, match="is not the one results= was computed with"):
        sa.cox_model_figure(ds, "SYN1:SE:1", "COH1", "OS", settings=FAST, results=results, model=other)
    one = sa.analyze(ds, events=["SYN1:SE:1"], cohorts=["COH1"], endpoints=["OS"])
    with pytest.raises(InputError, match="results= does not cover cohort.s. COH2 for OS"):
        sa.event_panel(ds, "SYN1:SE:1", ["COH1", "COH2"], "OS", settings=FAST, results=one)
    p = sa.event_panel(ds, "SYN1:SE:1", ["COH1"], "OS", settings=FAST, results=results, model=sa.CoxModel())
    assert p.provenance["call"]["model"] == sa.CoxModel().to_dict()


def test_the_cox_figure_hashes_the_q_family_of_its_results(ds, results):
    f = sa.cox_model_figure(ds, "SYN1:SE:1", "COH1", "OS", settings=FAST, results=results)
    assert f.provenance["input_rows"]["psi_q_family"] == 2          # SYN1:A3SS:1 and SYN1:RI:1 share its q family
    alone = sa.cox_model_figure(ds, "SYN1:SE:1", "COH1", "OS", settings=FAST)
    assert alone.provenance["input_rows"]["psi_q_family"] is None


def test_low_power_fits_are_marked_on_the_page(ds, results):
    """‡ follows the forest CI of a fit with fewer than cox_low_power_events events; the legend says so."""
    from splice_assay.plot.forest import LOW_MARK
    p = sa.event_panel(ds, "SYN1:SE:1", ["COH1"], "OS", settings=FAST.replace(cox_low_power_events=10_000),
                       results=results)
    f = p.table[p.table.panel.eq("forest")]
    drawn = f[f.hr_per_sd.notna()]
    assert len(drawn) and drawn.low_power_marked.all()
    assert any(LOW_MARK in t.get_text() for a in p.figure.axes for t in a.texts)
    assert "fewer than 10000 events: low power" in _texts(p)
    q = sa.event_panel(ds, "SYN1:SE:1", ["COH1"], "OS", settings=FAST.replace(cox_low_power_events=0), results=results)
    assert not q.table[q.table.panel.eq("forest")].low_power_marked.any()
    assert not any(LOW_MARK in t.get_text() for a in q.figure.axes for t in a.texts)


def test_no_host_expression_gives_a_psi_only_fit(tables):
    """Without expression for the host gene the Cox model is PSI alone, with a note, and the page names it so."""
    ex = tables["expression"]
    ds = sa.Dataset.from_tables(**dict(tables, expression=ex[ex.gene.ne("SYN1")]))
    sv = sa.analyze(ds, events=["SYN1:SE:1"], endpoints=["OS"]).survival
    t = sv[sv.cox_status.eq("tested")]
    assert len(t) and t.cox_model.eq("PSI").all()
    assert t.cox_notes.str.contains("host expression left out (no values)", regex=False).all()
    p = sa.event_panel(ds, "SYN1:SE:1", ["COH1"], "OS", settings=FAST)
    assert p.table[p.table.panel.eq("forest_axis")].cox_model.iloc[0] == "Cox: PSI"
    # expression missing for one cohort only: PSI alone there, and the forest note says where
    s = tables["samples"]
    coh3 = set(s.sample_id[s.cohort.eq("COH3")])
    ds3 = sa.Dataset.from_tables(**dict(tables, expression=ex[~ex.sample_id.isin(coh3)]))
    p3 = sa.event_panel(ds3, "SYN1:SE:1", ["COH1"], "OS", settings=FAST)
    note = p3.table[p3.table.panel.eq("forest_axis")].cox_model.iloc[0]
    assert note.startswith("Cox: PSI + host expression; PSI alone in COH3"), note


def test_the_forest_note_names_fits_without_expression(tables):
    """Two events of one page, one of them without expression for its host gene: the note says where PSI is alone."""
    from splice_assay.analysis import describe_model
    ev = tables["events"].copy()
    ev["expression_gene"] = np.where(ev.event_id.eq("SYN1:A3SS:1"), "NOPE", ev.gene)      # others: their own gene
    ds = sa.Dataset.from_tables(**dict(tables, events=ev))
    p = sa.event_panel(ds, ["SYN1:SE:1", "SYN1:A3SS:1"], ["COH1"], "OS", settings=FAST)
    note = p.table[p.table.panel.eq("forest_axis")].cox_model.iloc[0]
    assert note.startswith("Cox: PSI + host expression; PSI alone in "), note
    assert describe_model(sa.CoxModel(), ds, ["SYN1:SE:1", "SYN1:A3SS:1"]) == \
        "PSI + host expression; PSI alone where the host gene has no expression"
    assert describe_model(sa.CoxModel(), ds, ["SYN1:A3SS:1"]) == "PSI"
