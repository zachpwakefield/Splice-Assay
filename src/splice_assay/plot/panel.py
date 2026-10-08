"""The event panel: schematic on top; per cell shown (event x cohort) a case-vs-reference view (matched pairs beside
all samples) and a KM panel; a forest of every cohort on the right; optionally the full Cox model of each cell shown
(the model band); a legend at the bottom. One endpoint per figure.

Every drawn number comes from `analysis.analyze` on the same Dataset; the drawing recomputes pair sets, arms and
events from the data and stops if they disagree with the analysis (that would be a bug).
"""
from __future__ import annotations

import hashlib
import warnings

import numpy as np
import pandas as pd

from .. import provenance as P
from ..analysis import Results, analyze, analyze_expression, event_settings
from ..annotation import GeneModel, gene_model, read_gtf
from ..config import Settings
from ..dataset import Dataset, InputError, read_table
from .. import events as EV
from ..events import geometry
from ..stats.survival import CoxModel, km_high, ridge_text, ties_high
from ..stats.tissue import patient_values
from . import forest as F
from . import km as K
from . import protein as PR
from . import schematic as SCH
from . import style as S
from . import tissue as T
from .model import SCALE_NOTE, _range, display_rows, draw_terms, model_terms, not_fitted, terms_columns
from .panel_common import ENDPOINT_NAMES, Panel, model_mismatch, safe_name

KM_MESSAGE = {"coverage_gate": "PSI too sparse in this cohort\n(coverage gate)",
              "too_few_patients_or_events": "too few patients or events\nfor a log-rank test",
              "low_psi_variance": "PSI nearly constant\nin this cohort"}
NAN = float("nan")


def _as_list(x) -> list[str]:
    return [str(x)] if isinstance(x, (str, int)) else [str(v) for v in x]


def _get(row, key, default=NAN):
    v = row.get(key, default) if hasattr(row, "get") else default
    return default if v is None else v


_safe = safe_name
W_DETAIL = 6.96                             # inches of text for the event-detail lines under the title


def _highlight(h, events, endpoint) -> tuple[dict, dict]:
    """(event, cohort) -> label, and label -> colour, from a table (event_id, cohort, label[, endpoint, colour])
    or a dict {(event, cohort): label}."""
    if h is None:
        return {}, {}
    if isinstance(h, dict):
        return {(str(e), str(c)): str(v) for (e, c), v in h.items() if str(e) in events}, {}
    df = h if isinstance(h, pd.DataFrame) else read_table(h)
    missing = [c for c in ("event_id", "cohort", "label") if c not in df.columns]
    if missing:
        raise InputError(f"highlight: missing column(s) {', '.join(missing)}")
    if "endpoint" in df.columns:
        ep = df.endpoint.astype(str).str.strip()
        df = df[ep.eq(endpoint) | df.endpoint.isna() | ep.eq("")]
    marks = {(str(r.event_id), str(r.cohort)): ("" if pd.isna(r.label) else str(r.label)) for r in df.itertuples()
             if str(r.event_id) in events}
    colors = {}
    if "colour" in df.columns or "color" in df.columns:
        col = "colour" if "colour" in df.columns else "color"
        colors = {str(r.label): str(getattr(r, col)) for r in df.itertuples() if pd.notna(getattr(r, col))}
    return marks, colors


def _model(gtf, ev: pd.DataFrame, geoms, s: Settings) -> GeneModel | None:
    if gtf is None:
        return None
    chrom, gene, gid = ev.chrom.iloc[0], ev.gene.iloc[0], ev.gene_id.iloc[0]
    lo = min(g.span[0] for g in geoms)
    hi = max(g.span[1] for g in geoms)
    if isinstance(gtf, pd.DataFrame):
        return gene_model(gtf, gene, chrom, (lo, hi), gid, s)
    win = (chrom, lo - s.gtf_flank, hi + s.gtf_flank)
    tab = read_gtf(gtf, [win], genes=[g for g in (gene, gid) if g])
    m = gene_model(tab, gene, chrom, (lo, hi), gid, s)
    if m is not None:
        rec = tab[tab.feature.eq("gene") & (tab.gene_name.eq(m.gene) | (tab.gene_id.eq(gid.split(".")[0]) if gid
                                                                         else False))]
        g0, g1 = int(rec.start.min()), int(rec.end.max())
        if g0 < win[1] or g1 > win[2]:              # a long host gene: read its whole span for nested genes
            tab = read_gtf(gtf, [(chrom, g0 - s.gtf_flank, g1 + s.gtf_flank)], genes=[g for g in (gene, gid) if g])
            m = gene_model(tab, gene, chrom, (lo, hi), gid, s)
    return m


def _and(items) -> str:
    items = list(items)
    return items[0] if len(items) == 1 else ", ".join(items[:-1]) + " and " + items[-1]


def _stacked_layout(W, has_ref, marks, label_w):
    """One row per cell: comparison view, KM, and the cell's model at the right; the forest below at full width."""
    km_x0, km_w = (0.47 + 1.45 + 0.60, 1.30) if has_ref else (0.95, 2.30)
    mdl_x0 = km_x0 + km_w + 0.50
    lab_x = 0.12 + label_w + 0.04
    mark_x = lab_x + 0.11
    x_a = (mark_x + 0.13) if marks else lab_x + 0.10
    avail = W - 0.15 - x_a
    if has_ref:
        aw = min(2.6, (avail - 0.45) / 2)
        axL, axR = (x_a, aw), (x_a + aw + 0.45, aw)
    else:
        axL, axR = None, (x_a, min(3.2, avail))
    return dict(tn=(0.47, 1.45) if has_ref else None, km=(km_x0, km_w), model=(mdl_x0, W - 0.12 - mdl_x0),
                lab_x=lab_x, mark_x=mark_x, axL=axL, axR=axR)


def _legend_lines(W, first, second):
    """Split legend items into lines that fit the figure width (inches)."""
    lines = []
    for items in (first, second):
        cur, x = [], 0.12
        for kind, text in items:
            w = 0.21 + S.text_width(text, 6) + 0.2
            if cur and x + w > W - 0.05:
                lines.append(cur)
                cur, x = [], 0.12
            cur.append((kind, text))
            x += w
        lines.append(cur)
    return lines


def _layout(W, has_ref, marks, label_w):
    """Horizontal positions (inches) of the left panels and the forest for a figure of width W."""
    axR_w = 0.95 if has_ref else 1.6
    x_axR = W - 0.13 - axR_w
    x_axL = x_axR - 0.21 - 0.85 if has_ref else None
    first = x_axL if has_ref else x_axR
    mark_x = first - 0.13
    lab_x = mark_x - 0.11 if marks else first - 0.08
    tn_x0, tn_w = (0.47, 1.55) if has_ref else (None, 0.0)
    km_x0 = 0.47 + 1.55 + 0.62 if has_ref else 0.95
    km_w = lab_x - label_w - 0.50 - km_x0
    return dict(axR=(x_axR, axR_w), axL=(x_axL, 0.85) if has_ref else None, mark_x=mark_x, lab_x=lab_x,
                tn=(tn_x0, tn_w) if has_ref else None, km=(km_x0, km_w))


def event_panel(ds: Dataset, events, cohorts, endpoint: str, *, settings: Settings | None = None,
                model: CoxModel | None = None, gtf=None, highlight=None, highlight_title: str | None = None,
                highlight_legend: str | None = None, focus=None, cox_name: str = "Cox",
                endpoint_label: str | None = None, out_dir=None, stem: str | None = None,
                results: Results | None = None, detail=None, layout: str = "auto", proteins=None,
                detail_results: Results | None = None, part: tuple | None = None, shade_marks: bool = True,
                row_notes: dict | None = None, title_note: str | None = None) -> Panel:
    """Draw (and, with out_dir, save) the panel of one or two events of one gene for one endpoint.

    events           one or two event IDs of the same gene
    cohorts          the cohorts that get a comparison view and a KM panel (one row per event x cohort)
    model            the Cox model (default: PSI + host expression when an expression table is given)
    gtf              a GTF path (.gtf or .gtf.gz) or a table from annotation.read_gtf; without it the gene track and
                     nested small RNAs are omitted
    highlight        optional table or dict marking cells in the forest with a short label (e.g. a tier letter):
                     columns event_id, cohort, label, optional endpoint and colour
    focus            explicit [(event, cohort)] rows instead of every event x cohort
    results          precomputed analyze() results covering these events, every cohort and this endpoint; its
                     q values (BH within each gene) are printed beside the p values, so for gene-wide q it should
                     cover all of the gene's events (computed so when not given)
    detail_results   the same for the model rows' model (computed when not given)
    detail           add the full Cox model of each cell shown, with every term. True uses `model`; a CoxModel uses
                     that one (e.g. model.with_clinical(["age", "stage"]))
    layout           "side": the forest right of the cells and the models in a band below; "stacked": one row per
                     cell with its model at the right and the forest below at full width; "auto": stacked when models
                     are shown for three or more cells
    proteins         a protein cache (protein.ProteinCache or its folder) or {event: ProteinChange}: adds a band
                     under the schematic with the suggested protein of each form, when there is a good match
    part             (i, n) for page i of n when the cohorts are split over several pages (`page_parts`): the title
                     says so and the default stem ends in _p<i>of<n>; each page keeps the full forest
    shade_marks      False: the forest shades only the cells shown at left, not every highlighted one
    row_notes        {(event, cohort): text} written after each shown cell's heading (e.g. its lines of evidence)
    title_note       a few words added to the title (e.g. how its cohorts were chosen)
    """
    import matplotlib
    from matplotlib.figure import Figure
    from matplotlib.lines import Line2D
    from matplotlib.patches import Rectangle

    s = settings or Settings()
    events, cohorts = _as_list(events), _as_list(cohorts)
    if not 1 <= len(events) <= 2:
        raise InputError("a panel shows one or two events")
    for e in events:
        if e not in ds.events.index or e not in ds.psi.index:
            raise InputError(f"event {e} is not in the data")
    ev = ds.events.loc[events]
    for col in ("gene", "chrom", "strand"):
        if ev[col].nunique() > 1:
            raise InputError(f"the events of one panel must share {col}: {', '.join(ev[col])}")
    chrom, strand, gene = ev.chrom.iloc[0], ev.strand.iloc[0], ev.gene.iloc[0]
    if not chrom or strand not in ("+", "-"):
        raise InputError(f"event {events[0]}: chrom and strand are needed to draw it")
    geoms = [geometry(e, ev.loc[e]) for e in events]
    qty = EV.quantity(geoms[0].event_type)                  # PSI, or HIT index for HIT events
    if any(EV.quantity(g.event_type) != qty for g in geoms):
        raise InputError("the events of one page must measure the same quantity (PSI, or the HIT index)")
    q_lim = (-1.05, 1.05) if qty != "PSI" else None
    s_ev = event_settings(geoms[0].event_type, settings or Settings())     # its hit threshold
    if endpoint not in ds.endpoints:
        raise InputError(f"endpoint {endpoint} is not in the survival table")
    for c in cohorts:
        if c not in ds.cohorts:
            raise InputError(f"cohort {c} is not in the samples table")
    focus = [(e, c) for e in events for c in cohorts] if focus is None else [(str(e), str(c)) for e, c in focus]
    same = ds.events.gene.isin(set(ev.gene))                # the q family: the gene's events of this quantity
    gene_ids = [e for e in ds.events.index[same] if EV.quantity(ds.events.at[e, "event_type"]) == qty]
    hit_apart = qty == "PSI" and bool(ds.events.event_type[same].map(EV.opt_in).any())   # HIT has its own family
    if results is None:
        results = analyze(ds, events=gene_ids, endpoints=[endpoint], settings=s, model=model)
    elif model is not None and model != results.model:            # the page would name a model it did not draw
        raise InputError(model_mismatch(model, results.model))
    missing = [e for e in events if e not in set(results.groups.event_id)]
    if missing:
        raise InputError(f"results= does not cover {', '.join(missing)}"
                         + (" (HIT-index events are analyzed only on request: include_hit=True)"
                            if any(EV.opt_in(ds.events.at[e, "event_type"]) for e in missing) else "")
                         + "; leave results out to analyze them here")
    model = results.model if model is None else model
    tn = results.groups[results.groups.event_id.isin(events)].set_index(["event_id", "cohort"])
    sv = results.survival[results.survival.event_id.isin(events) & results.survival.endpoint.eq(endpoint)]
    sv = sv.set_index(["event_id", "cohort"])
    marks, mark_colors = _highlight(highlight, events, endpoint)
    absent = sorted({c for _, c in marks if c not in set(ds.cohorts)})
    known = ((ds.notes or {}).get("subset") or {}).get("all_cohorts") or ds.cohorts    # before any subset
    typos = [c for c in absent if c not in set(known)]
    if typos:                                                   # an error, as a --cohort typo is
        raise InputError(f"highlight: cohort(s) {', '.join(typos)} not in the data (cohorts: {', '.join(known)})")
    if absent:                                                  # cohorts the --keep / --where subset left out
        warnings.warn(f"highlight: cohort(s) {', '.join(absent)} not in the subset; their rows are not drawn",
                      stacklevel=2)
        marks = {k: v for k, v in marks.items() if k[1] not in absent}
    lacking = sorted(c for c in {c for _, c in focus} | {c for _, c in marks}
                     if any((e, c) not in tn.index or (e, c) not in sv.index for e in events))
    if lacking:
        raise InputError(f"results= does not cover cohort(s) {', '.join(lacking)} for {endpoint}; leave results out "
                         "to analyze them here")
    shade = (set(marks) if shade_marks else set()) | set(focus)
    labels = [ev.at[e, "label"] for e in events]
    glab = dict(case=s.case_label or ds.labels["case"], reference=s.reference_label or ds.labels["reference"])
    ylabel = endpoint_label or ENDPOINT_NAMES.get(endpoint, endpoint)
    has_ref = ds.has_reference
    use_expr = ds.expression is not None and model.expression and \
        any(ds.events.at[e, "expression_gene"] in ds.expression.index for e in events)   # else PSI alone
    model_note = f"Cox: {model.describe(use_expr).replace('PSI', qty, 1)}" + \
        ridge_text(results.settings, bool(model.covariates), [qty] + (["host expression"] if use_expr else []))

    # ------------------------------------------------------------------ forest records
    fc = {c for _, c in focus} | {c for (_, c) in marks}
    for e in events:
        t_, v_ = tn.loc[e], sv.loc[e]
        fc |= set(t_.index[t_.paired_status.eq("tested") | t_.unpaired_status.eq("tested")])
        fc |= set(v_.index[v_.cox_status.isin(["tested", "failed"])])
    fcoh = sorted(fc)
    recs_by = {}
    for c in fcoh:
        recs = []
        for e in events:
            t_, v_ = tn.loc[(e, c)], sv.loc[(e, c)]
            d = dict(event_id=e, mark=marks.get((e, c), ""), shade=(e, c) in shade, cox_status=v_.cox_status)
            for design in ("paired", "unpaired"):
                st = t_[f"{design}_status"]
                d[f"{design}_status"] = st
                d[f"{design}_delta"] = float(_get(t_, f"{design}_delta_median")) if st == "tested" else NAN
                d[f"{design}_p"] = float(_get(t_, f"{design}_p")) if st == "tested" else NAN
            d["paired_n"] = _get(t_, "paired_n_pairs")
            ok = v_.cox_status == "tested"
            u = s.psi_hr_unit                               # the HR per SD, or per IQR
            for k_out, k_in in (("hr", f"hr_per_{u}"), ("lo", f"ci_low_{u}"), ("hi", f"ci_high_{u}"),
                                ("cox_p", "cox_p"), ("cox_q", "cox_q"), ("cox_beta", "cox_beta"), ("cox_se", "cox_se"),
                                ("ph_p", "ph_p")):
                d[k_out] = float(_get(v_, k_in)) if ok else NAN
            for k in ("cox_n", "cox_events"):
                d[k] = _get(v_, k)
            d["spread"] = _get(v_, f"psi_{u}")
            d["cox_model"] = _get(v_, "cox_model", "")
            recs.append(d)
        recs_by[c] = recs
    two = len(events) == 2
    alone = [(f"{ev.at[r['event_id'], 'label']} " if two else "") + c for c in fcoh for r in recs_by[c]
             if r["cox_status"] == "tested" and "host expression" not in str(r["cox_model"])]
    if use_expr and alone:                                  # fits without usable host expression: PSI alone there
        model_note += (f"; {qty} alone in " + (", ".join(alone) if len(alone) <= 2 else
                                                f"{len(alone)} {'fits' if two else 'cohorts'}"))

    # ------------------------------------------------------------------ layout (inches)
    gm = _model(gtf, ev, geoms, s)
    n_ev, n_rows = len(events), len(focus)
    label_w = max(S.text_width(c, 6.0, weight="bold") for c in fcoh)
    detail_model = (model if detail is True else detail) if detail else None
    if layout not in ("auto", "side", "stacked"):
        raise InputError('layout must be "auto", "side" or "stacked"')
    stacked = detail_model is not None and (layout == "stacked" or (layout == "auto" and n_rows >= 3))
    details = [(lab, EV.describe(g, str(ev.at[e, "chrom"]), str(ev.at[e, "strand"])))
               for e, lab, g in zip(events, labels, geoms)]
    detail_lines = [(lab, line) for lab, d in details
                    for line in PR._wrap(f"{lab}: {d}", W_DETAIL, 6.2)]  # under the title, one per event (wrapped)
    top, sch_h, gap1 = 0.30 + 0.14 * len(detail_lines), SCH.height(n_ev, gm is not None), 0.14
    row_title, tab_h = (0.52, 0.66) if has_ref else (0.46, 0.62)
    qcols = [c for c in ("paired_q", "unpaired_q") if c in tn.columns]
    any_q = bool(any(np.isfinite(float(_get(tn.loc[f], c))) for f in focus for c in qcols if f in tn.index) or
                 any(np.isfinite(float(_get(sv.loc[f], "km_q"))) for f in focus if f in sv.index))
    if any_q:
        row_title += 0.10                                   # a line for q above the views
    if _any_km_ph([_get(sv.loc[f], "km_ph_p") for f in focus if f in sv.index], s):
        row_title += 0.10                                   # a line for a KM proportional-hazards note
    forest_h = (0.106 if n_ev == 1 else 0.18) * len(fcoh)
    W = 7.2
    if detail_model is not None and detail_results is None:
        same = detail_model == results.model and results.settings == s
        detail_results = results if same else analyze(ds, events=gene_ids, endpoints=[endpoint], settings=s,
                                                      model=detail_model, cox_only=True)
    elif detail_model is not None and detail_results.model != detail_model:
        raise InputError(model_mismatch(detail_model, detail_results.model, "detail_results="))
    band = _model_band(ds, focus, endpoint, s, detail_model, W, detail_results) if detail_model is not None else None
    if stacked:
        L = _stacked_layout(W, has_ref, marks, label_w)
        n_terms = max([len(c["disp"]) for c in band["cells"]] + [3])
        band["stack_extra"] = _head_extra(band["cells"], L["model"][1], 5.6, lambda c: _cell_model_text(c["row"], band),
                                          2)
        row_h = max(row_title + 1.02 + tab_h, 0.50 + band["stack_extra"] + BAND_PITCH * n_terms + BAND_AXIS + 0.10)
        ax_h = row_h - row_title - tab_h
        main_h = n_rows * row_h
        after_h = 0.20 + 0.25 + forest_h + 0.85               # the forest below the rows
    else:
        L = _layout(W, has_ref, marks, label_w)
        if L["km"][1] < 1.35:                               # long cohort names: widen rather than squeeze the KM
            W += 1.35 - L["km"][1]
            L = _layout(W, has_ref, marks, label_w)
            if band is not None:
                band = _model_band(ds, focus, endpoint, s, detail_model, W, detail_results)
        ax_h = min(1.9, max(1.02, (forest_h + 0.35 - n_rows * (row_title + tab_h)) / n_rows))
        row_h = row_title + ax_h + tab_h
        main_h = max(n_rows * row_h, forest_h + 0.62 + 0.22)
        after_h = band["height"] if band else 0.0
    nested_label = "snoRNAs" if set(s.nested_biotypes) <= {"snoRNA", "scaRNA"} else "Nested genes"
    col_ev = S.TYPE_COLOR.get(geoms[0].event_type, S.TYPE_DEFAULT)
    paired_drawn = has_ref and any(_get(tn.loc[f], "paired_status", "") == "tested" for f in focus)
    all_drawn = has_ref and any(_get(tn.loc[f], "unpaired_status", "") == "tested" for f in focus)
    first = [("exon", "Exon"), ("var", f"Region measured by {qty}")]
    if any(g.other for g in geoms):
        first.append(("other", "Other exon"))
    span = (min(g.span[0] for g in geoms), max(g.span[1] for g in geoms))
    win = SCH.window_for(span, gm)
    nested_shown = gm is not None and len(gm.nested if not win else gm.nested[(gm.nested.end > win[0]) &
                                                                                (gm.nested.start < win[1])])
    if nested_shown:
        first.append(("sno", nested_label[:-1] if nested_label.endswith("s") else nested_label))
    if any(g.psi_arcs or g.other_arcs for g in geoms):      # junction arcs are drawn (none for AFE, ALE or HIT)
        first += [("arc_up", "Junction of the PSI form"), ("arc_dn", "of the other form")]
    cl = S.in_sentence(glab["case"])
    second = [("up", f"Higher in {cl}"), ("down", f"Lower in {cl}")] if paired_drawn else []
    if all_drawn:
        second.append(("split", "KM split"))
    if has_ref:
        second += [("pd", "Paired Δ"), ("ud", "Unpaired Δ")]
    mark_text = highlight_legend or ((f"{highlight_title} cell (letter: {highlight_title.lower()})" if highlight_title
                                      else "Highlighted cell") if marks else "Cohort shown at left")
    band_text = mark_text if shade_marks else "Cohort shown at left"
    second += [("hs", f"{cox_name} p < {s.alpha:g}"), ("hn", f"{cox_name} p ≥ {s.alpha:g}"), ("band", band_text)]
    if any(marks.values()) and not shade_marks:             # the letters mark other cells than the shading
        second.append(("mark", mark_text))
    if band is not None and any(c["disp"] for c in band["cells"]):     # the model rows' other terms
        second += [("ts", f"Model term p < {s.alpha:g}"), ("tn", f"p ≥ {s.alpha:g}")]
    if _q_marked(focus, fcoh, tn, sv, band, s):
        second.append(("q", f"q < {s.q_mark_below:g} (Benjamini–Hochberg)"))
    if _ph_marked(focus, fcoh, sv, band, s):
        second.append(("ph", f"non-proportional hazards (p < {s.ph_note_below:g})"))
    if any(np.isfinite(r["hr"]) and F.low_power(r["cox_events"], s.cox_low_power_events)
           for c in fcoh for r in recs_by[c]):
        second.append(("low", f"fewer than {s.cox_low_power_events} events: low power"))
    lines = _legend_lines(W, first, second)
    legend_h = 0.42 + 0.15 * (len(lines) - 2)
    fam = "HIT-index events" if qty != "PSI" else ("PSI events (HIT index apart)" if hit_apart else "events")
    q_lines = PR._wrap(_q_note(gene, endpoint, tn, sv, detail_results, s, fam), W - 0.26, 5.4)
    legend_h += 0.05 + 0.115 * len(q_lines)
    if band is not None:                                    # one column layout for every model forest (aligned)
        band["cols"] = terms_columns(L["model"][1] if stacked else band["width"],
                                     [c["disp"] for c in band["cells"] if c["disp"]] or [[]], 5.8 if stacked else 5.9,
                                     s.ph_note_below, s.q_mark_below)
    changes = _protein_changes(proteins, ds, events, geoms)
    pbands = [b for b in (PR.prepare(changes.get(e), W, col_ev, lab if n_ev == 2 else "")
                          for e, lab in zip(events, labels)) if b is not None]
    prot_h = sum(b["height"] for b in pbands) + (0.06 if pbands else 0.0)
    H = top + sch_h + prot_h + gap1 + main_h + after_h + legend_h
    fx = lambda x: x / W                                                     # noqa: E731
    fy = lambda y: 1 - y / H                                                 # y in inches from the top  # noqa: E731
    rows: list[dict] = []

    with matplotlib.rc_context(S.rc()):
        fig = Figure(figsize=(W, H))
        # title
        fig.text(fx(0.12), fy(0.10), gene, fontsize=8.5, fontweight="bold", style="italic", ha="left", va="top")
        tw = 0.12 + S.text_width(gene, 8.5, style="italic", weight="bold") + 0.05
        shown = list(dict.fromkeys(c for _, c in focus))
        sub = (ds.notes or {}).get("subset")
        tail = [endpoint] + ([f"subset: {_subset_name(sub)}"] if sub else []) + (
            [f"page {part[0]} of {part[1]}"] if part else []) + ([title_note] if title_note else [])
        title_rest = " · ".join([_and(labels), _and(shown)] + tail)
        if tw + S.text_width(title_rest, 8.5, weight="bold") > W - 0.15:     # too many to list: count them
            title_rest = " · ".join([_and(labels), f"{len(shown)} cohorts"] + tail)
        fig.text(fx(tw), fy(0.10), title_rest, fontsize=8.5, fontweight="bold", ha="left", va="top")
        for i, (_, line) in enumerate(detail_lines):           # what the event is, in 1-based coordinates
            fig.text(fx(0.12), fy(0.30 + 0.14 * i), line, fontsize=6.2, color=S.INK2, ha="left", va="top")
            rows.append(dict(panel="event_detail", text=line))
        # schematic
        ax_s = fig.add_axes([fx(1.05), fy(top + sch_h), fx(W - 1.05 - 0.15), sch_h / H])
        span = (min(g.span[0] for g in geoms), max(g.span[1] for g in geoms))
        SCH.draw(ax_s, fig, [dict(event_id=e, label=lab, geometry=g) for e, lab, g in zip(events, labels, geoms)],
                 chrom, strand, gm, rows, nested_label, window=SCH.window_for(span, gm))
        for e, g in zip(events, geoms):
            rows.append(dict(panel="printed", tag=e, what="variable region length (nt)", value=g.variable_length,
                             source="events.variable"))
        # ------------------------------------------------------------------ the protein band(s)
        yp = top + sch_h + 0.06
        for b in pbands:
            PR.draw(fig, W, H, yp, b, rows)
            yp += b["height"]
        # ------------------------------------------------------------------ left: one row per cell shown
        y_main = top + sch_h + prot_h + gap1
        printed = []
        for ri, (e, c) in enumerate(focus):
            y_row = y_main + ri * row_h
            t_, v_ = tn.loc[(e, c)], sv.loc[(e, c)]
            mark = marks.get((e, c), "")
            head = ([ev.at[e, "label"]] if n_ev == 2 else []) + [c]
            if mark:
                head.append(f"{highlight_title} {mark}" if highlight_title else mark)
            hx = 0.12
            if n_ev == 2:                                   # the event's shade in the forest
                fig.add_artist(Line2D([fx(0.17)], [fy(y_row + 0.085)], marker="o", ms=4.2, ls="",
                                      color=(S.INK, S.SECOND)[events.index(e)]))
                hx = 0.27
            fig.text(fx(hx), fy(y_row + 0.02), " · ".join(head), fontsize=7.2, fontweight="bold", ha="left", va="top")
            if (row_notes or {}).get((e, c)):
                fig.text(fx(hx + S.text_width(" · ".join(head), 7.2, weight="bold") + 0.12), fy(y_row + 0.03),
                         row_notes[(e, c)], fontsize=6.2, color=S.INK2, ha="left", va="top")
                rows.append(dict(panel="printed", tag=f"{e}|{c}", what="row note", value=row_notes[(e, c)],
                                 source="row_notes"))
            tag = f"{e}|{c}"
            v = ds.psi.loc[e]
            if has_ref:
                ax_tn = fig.add_axes([fx(L["tn"][0]), fy(y_row + row_title + ax_h), fx(L["tn"][1]), ax_h / H])
                _group_view(fig, ax_tn, ds, e, c, v, t_, v_, s, glab, rows, printed, tag, quantity=qty, lim=q_lim)
            ax_km = fig.add_axes([fx(L["km"][0]), fy(y_row + row_title + ax_h), fx(L["km"][1]), ax_h / H])
            _km_view(fig, ax_km, ds, e, c, v, v_, endpoint, ylabel, s, rows, printed, tag, what=qty)
            if stacked:                                     # this cell's model, right of its KM
                rows += _draw_cell_model(fig, W, H, y_row + 0.06, L["model"], band, band["cells"][ri], s)
            rows.append(dict(panel="focus_cell", tag=tag, endpoint=endpoint, event_id=e, cohort=c, label=mark,
                             paired_hit=_get(t_, "paired_hit"), unpaired_hit=_get(t_, "unpaired_hit"),
                             within_patient_support=_get(t_, "within_patient_support"),
                             km_p=_get(v_, "km_p"), cox_p=_get(v_, "cox_p"),
                             **{f"hr_per_{s.psi_hr_unit}": _get(v_, f"hr_per_{s.psi_hr_unit}")}))
        # ------------------------------------------------------------------ forest
        f_top = y_main + (main_h + 0.20 if stacked else 0.0)       # below the rows when stacked
        fy0, fh = fy(f_top + 0.25 + forest_h), forest_h / H
        axL = fig.add_axes([fx(L["axL"][0]), fy0, fx(L["axL"][1]), fh]) if has_ref else None
        axR = fig.add_axes([fx(L["axR"][0]), fy0, fx(L["axR"][1]), fh])
        lim, xlo, xhi, n_clip = F.axis_limits(recs_by)
        F.set_axes(axL, axR, lim, xlo, xhi, s_ev.min_abs_delta, glab, model_note, quantity=qty,
                   hr_label=f"HR per {s.psi_hr_unit.upper()}\nof {qty} (95% CI)")
        F.draw(fig, axL, axR, fx(L["lab_x"]), fx(L["mark_x"]), fcoh, recs_by, n_ev, s.alpha, cox_name, mark_colors,
               ph_below=s.ph_note_below, q_below=s.q_mark_below, low_below=s.cox_low_power_events)
        rows.append(dict(panel="forest_axis", hr_axis_low=xlo, hr_axis_high=xhi, n_ci_clipped=n_clip,
                         cox_model=model_note))
        u = s.psi_hr_unit
        for c in fcoh:
            for r in recs_by[c]:
                drawn = bool(np.isfinite(r["hr"]))          # the marks follow a drawn CI (forest.draw)
                rows.append(dict(panel="forest", endpoint=endpoint, event_id=r["event_id"], cohort=c, label=r["mark"],
                                 paired_status=r["paired_status"], paired_n_pairs=r["paired_n"],
                                 paired_delta_median=r["paired_delta"], paired_p=r["paired_p"],
                                 unpaired_status=r["unpaired_status"], unpaired_delta_median=r["unpaired_delta"],
                                 unpaired_p=r["unpaired_p"], cox_status=r["cox_status"], cox_n=r["cox_n"],
                                 cox_events=r["cox_events"], cox_beta=r["cox_beta"], cox_se=r["cox_se"],
                                 **{f"psi_{u}": r["spread"], f"hr_per_{u}": r["hr"], f"ci_low_{u}": r["lo"],
                                    f"ci_high_{u}": r["hi"]}, cox_p=r["cox_p"], cox_q=r["cox_q"],
                                 q_marked=drawn and bool(s.q_mark_below > 0 and r["cox_q"] < s.q_mark_below),
                                 ph_p=r["ph_p"],
                                 ph_marked=drawn and bool(s.ph_note_below > 0 and r["ph_p"] < s.ph_note_below),
                                 low_power_marked=drawn and F.low_power(r["cox_events"], s.cox_low_power_events),
                                 cox_model=r["cox_model"]))
        head_y = fy(f_top + 0.25 - 0.04)
        if marks and highlight_title:
            fig.text(fx(L["mark_x"]), head_y, highlight_title, fontsize=5.6, color=S.INK2, ha="center", va="bottom")
        if n_ev == 2:
            x0_ = (L["axL"] or L["axR"])[0]
            w_ = (L["axL"] or L["axR"])[1]
            for i, (lab, col) in enumerate(zip(labels, (S.INK, S.SECOND))):
                fig.text(fx(x0_ + w_ * (0.02 + 0.5 * i)), head_y, f"● {lab}", color=col, fontsize=5.8, ha="left",
                         va="bottom")
        # ------------------------------------------------------------------ the model band
        if band and not stacked:
            rows += _draw_band(fig, W, H, y_main + main_h, band, s, n_ev)
        elif band:
            fig.text(fx(0.12), fy(y_main + main_h + 0.02), "Cox models: " + SCALE_NOTE,
                     fontsize=5.6, color=S.MUTED, ha="left", va="top")
        # ------------------------------------------------------------------ legend
        for i, line in enumerate(q_lines):                 # the q footnote, at the very bottom
            fig.text(fx(0.12), fy(H - 0.08 - 0.115 * (len(q_lines) - 1 - i)), line, fontsize=5.4, color=S.MUTED,
                     ha="left", va="bottom")
        for li, items in enumerate(lines):
            x = 0.12
            yy = fy(H - legend_h + 0.10 + li * 0.15)
            for kind, text in items:
                if kind in ("exon", "var"):
                    fig.add_artist(Rectangle((fx(x), yy - 0.035 / H), fx(0.16), 0.07 / H, lw=0,
                                             facecolor=S.EXON if kind == "exon" else col_ev))
                elif kind == "other":
                    fig.add_artist(Rectangle((fx(x), yy - 0.03 / H), fx(0.16), 0.06 / H, facecolor="white",
                                             edgecolor=col_ev, lw=0.8))
                elif kind == "sno":
                    fig.add_artist(Rectangle((fx(x), yy - 0.022 / H), fx(0.16), 0.044 / H, facecolor=S.NESTED, lw=0))
                elif kind in ("arc_up", "arc_dn"):
                    t = np.linspace(0, 1, 30)
                    sgn = 1 if kind == "arc_up" else -1
                    fig.add_artist(Line2D(fx(x + 0.16 * t), yy - sgn * 0.035 / H + sgn * 0.09 / H * 4 * t * (1 - t),
                                          color=col_ev if kind == "arc_up" else S.EXON_LINE, lw=0.8))
                elif kind in ("up", "down"):
                    fig.add_artist(Line2D([fx(x), fx(x + 0.16)], [yy, yy], color=S.UP if kind == "up" else S.DOWN,
                                          lw=1.1))
                elif kind == "split":
                    fig.add_artist(Line2D([fx(x), fx(x + 0.16)], [yy, yy], color=S.INK, lw=0.6, ls=(0, (2, 1.5))))
                elif kind == "ph":
                    fig.text(fx(x + 0.08), yy, F.PH_MARK, fontsize=7, color=S.INK, ha="center", va="center")
                elif kind == "low":
                    fig.text(fx(x + 0.08), yy, F.LOW_MARK, fontsize=7, color=S.INK, ha="center", va="center")
                elif kind == "q":
                    fig.text(fx(x + 0.08), yy, S.Q_MARK, fontsize=7.5, color=S.INK, ha="center", va="center")
                elif kind == "band":
                    fig.add_artist(Rectangle((fx(x), yy - 0.05 / H), fx(0.16), 0.10 / H, facecolor=S.PAPER, lw=0))
                    if shade_marks and any(marks.values()):     # the first label present, in sorted order
                        lab0 = sorted(v for v in marks.values() if v)[0]
                        fig.add_artist(Rectangle((fx(x + 0.04), yy - 0.04 / H), fx(0.08), 0.08 / H, lw=0,
                                                 facecolor=mark_colors.get(lab0, S.HIGHLIGHT.get(
                                                     lab0, S.HIGHLIGHT_DEFAULT))))
                elif kind == "mark":                            # a letter as the forest draws it
                    lab0 = sorted(v for v in marks.values() if v)[0]
                    fig.add_artist(Rectangle((fx(x + 0.03), yy - 0.05 / H), fx(0.10), 0.10 / H, lw=0,
                                             facecolor=mark_colors.get(lab0, S.HIGHLIGHT.get(lab0,
                                                                                             S.HIGHLIGHT_DEFAULT))))
                    fig.text(fx(x + 0.08), yy, lab0, fontsize=5.0, ha="center", va="center", fontweight="bold",
                             color=S.INK if lab0 in S.DARK_LABELS else "white")
                else:
                    mk = "o" if kind in ("pd", "ud") else "s" if kind in ("ts", "tn") else "D"
                    fig.add_artist(Line2D([fx(x + 0.08)], [yy], marker=mk, ms=3.6 if mk == "o" else 3.2, ls="",
                                          markerfacecolor=S.INK if kind in ("pd", "hs", "ts") else "white",
                                          markeredgecolor=S.INK, mew=0.7))
                fig.text(fx(x + 0.21), yy, text, fontsize=6, color=S.INK2, ha="left", va="center")
                x += 0.21 + S.text_width(text, 6) + 0.2
    rows += printed
    table = pd.DataFrame(rows)
    shown = list(dict.fromkeys(c for _, c in focus))
    named = ([_safe(c) for c in shown] if len(shown) <= 4 else           # many cohorts: count them, plus a short hash
             [f"{len(shown)}cohorts", hashlib.sha1(",".join(shown).encode()).hexdigest()[:6]])
    stem = stem or "_".join([_safe(gene)] + [_safe(lab) for lab in labels] + named + [_safe(endpoint)]
                            + ([f"p{part[0]}of{part[1]}"] if part else []))
    prov = _provenance(ds, events, gene_ids, endpoint, s, model, detail_model, highlight,
                       dict(events=events, cohorts=cohorts, endpoint=endpoint, focus=focus, cox_name=cox_name,
                            gtf=None if gtf is None or isinstance(gtf, pd.DataFrame) else str(gtf),
                            detail_model=None if detail_model is None else detail_model.to_dict(),
                            proteins=_protein_source(proteins),
                            protein_status={e: changes[e].status for e in events if e in changes}))
    paths = {}
    if out_dir is not None:
        paths = S.save(fig, out_dir, stem, s.formats, s.dpi, table, prov)
    return Panel(figure=fig, table=table, paths=paths, provenance=prov, stem=stem)


# ============================================================================================ the model band
BAND_PITCH, BAND_HEAD, BAND_AXIS = 0.165, 0.46, 0.44


def _model_band(ds, focus, endpoint, s, model, W, results=None) -> dict:
    """The full Cox model of each cell shown: terms and layout (two cells per line)."""
    use_expr = ds.expression is not None and model.expression
    cells = []
    for e, c in focus:
        try:
            row, terms = model_terms(ds, e, c, endpoint, model, s, results)
            q = float(_get(row, "cox_q"))
            cells.append(dict(event=e, cohort=c, row=row, terms=terms, disp=display_rows(terms, q), error="", q=q,
                              qty=EV.quantity(ds.events.at[e, "event_type"])))
        except InputError as err:
            cells.append(dict(event=e, cohort=c, row=None, terms=None, disp=[], error=str(err)))
    ncol = 1 if len(cells) == 1 else 2
    nline = (len(cells) + ncol - 1) // ncol
    n_terms = max([len(x["disp"]) for x in cells] + [3])
    width = min(4.6, W - 0.24) if ncol == 1 else (W - 0.24 - 0.35) / 2
    describe = model.describe(use_expr)
    head_extra = _head_extra(cells, width, 5.8, lambda c: f"{_cell_model_text(c['row'], {'describe': describe})} · "
                                                          f"{int(c['row'].cox_n)} patients, {int(c['row'].cox_events)} events",
                             1)
    line_h = BAND_HEAD + head_extra + BAND_PITCH * n_terms + BAND_AXIS
    shown = [d["r"] for x in cells for d in x["disp"] if d["r"] is not None]
    xlim = _range(shown) if shown else None                   # one HR axis for every forest of the band
    return dict(cells=cells, ncol=ncol, nline=nline, line_h=line_h, width=width, xlim=xlim, head_extra=head_extra,
                height=0.12 + nline * line_h + 0.16, describe=describe, model=model)


def _subset_name(sub: dict) -> str:
    """How a subset is named on a page: the keep file's name and the where conditions, e.g. Her2_samples."""
    return sub.get("name") or f"{sub.get('patients')} patients"


def _family(frame, col):
    """The size of a q family (the same on every row of a gene), or None."""
    if frame is None or col not in frame.columns:
        return None
    v = frame[col].dropna()
    return int(v.max()) if len(v) else None


def _q_note(gene, endpoint, tn, sv, detail_results, s, fam: str = "events") -> str:
    """The footnote naming the q families, with their sizes as computed. `fam`: the events in them ('events', or
    'HIT-index events', which form families of their own)."""
    adj = None
    if detail_results is not None and len(detail_results.survival):
        d = detail_results.survival
        adj = _family(d[d.gene.eq(gene) & d.endpoint.eq(endpoint)], "cox_q_tests")
    fams = [("within patients", _family(tn, "paired_q_tests")), ("all samples", _family(tn, "unpaired_q_tests")),
            ("KM", _family(sv, "km_q_tests")), ("Cox, forest", _family(sv, "cox_q_tests")),
            ("model rows", adj)]
    sizes = ", ".join(f"{name} {n}" for name, n in fams if n)
    note = (f"q: Benjamini–Hochberg within {gene} for {endpoint}, one family per kind of test over all its {fam} × "
            f"cohorts ({sizes or 'no tests'} tests); none for families under {s.fdr_min_family}.")
    return note + (f" Settings changed from the defaults: {s.changed_text()}." if s.changed() else "")


def _q(cell, s) -> str:
    """' · PSI q 0.031 *' (or 'HIT index q'; * below q_mark_below) for a model whose tested term has a q, else ''."""
    q = cell.get("q", NAN)
    return f" · {cell.get('qty', 'PSI')} {S.q_text(q, s.q_mark_below)}" if np.isfinite(q) else ""


def _q_marked(focus, fcoh, tn, sv, band, s) -> bool:
    """Whether the page marks a q with *: a group view's q, a KM header's q, the forest (base model) or a model
    header."""
    if not s.q_mark_below > 0:
        return False
    low = (lambda v: v is not None and np.isfinite(float(v)) and float(v) < s.q_mark_below)
    if any(low(_get(tn.loc[f], c)) for f in focus if f in tn.index for c in ("paired_q", "unpaired_q")):
        return True
    if any(low(_get(sv.loc[f], "km_q")) for f in focus if f in sv.index):
        return True
    if any(low(_get(r, "cox_q")) for (e, c), r in sv.iterrows() if c in set(fcoh) and _get(r, "cox_status") == "tested"):
        return True
    return bool(band) and any(low(cell.get("q")) for cell in band["cells"])


def _ph_marked(focus, fcoh, sv, band, s) -> bool:
    """Whether the page draws a proportional-hazards dagger: a KM panel, the forest (base model) or a model row."""
    if not s.ph_note_below > 0:
        return False
    low = (lambda v: v is not None and np.isfinite(float(v)) and float(v) < s.ph_note_below)
    if any(low(_get(sv.loc[f], "km_ph_p")) for f in focus if f in sv.index):
        return True
    if any(low(_get(r, "ph_p")) for (e, c), r in sv.iterrows() if c in set(fcoh) and _get(r, "cox_status") == "tested"):
        return True
    return bool(band) and any(d["r"] is not None and low(getattr(d["r"], "ph_p", None))
                              for cell in band["cells"] for d in cell["disp"])


def _any_km_ph(values, s) -> bool:
    """Whether any KM panel prints the proportional-hazards note (its p below ph_note_below)."""
    return s.ph_note_below > 0 and any(np.isfinite(float(v)) and float(v) < s.ph_note_below for v in values
                                       if v is not None)


def _cell_model_text(r, band) -> str:
    """The model fitted in this cell, with its notes (variables left out, merged levels, flags)."""
    txt = r.cox_model if isinstance(r.get("cox_model"), str) else band["describe"]
    notes = r.get("cox_notes")
    return txt + (f" ({notes})" if isinstance(notes, str) and notes else "")


MODEL_LINE, MODEL_LINES = 0.085, 3                  # line pitch (inches) and most lines of a model's header text


def _model_lines(text: str, width: float, size: float) -> list[str]:
    """The model text wrapped to its box, at most MODEL_LINES lines; a cut ends in '…' (the CSV keeps it whole)."""
    lines = PR._wrap(text, width, size)
    if len(lines) > MODEL_LINES:
        lines = lines[:MODEL_LINES]
        lines[-1] = lines[-1].rstrip(" ;,(") + " …"
    return lines


def _head_extra(cells, width: float, size: float, text, fits: int) -> float:
    """Extra header height (inches) when a model's text wraps past the lines its header holds (`fits`)."""
    n = max([len(_model_lines(text(c), width, size)) for c in cells if c["row"] is not None] + [fits])
    return (n - fits) * MODEL_LINE


def _draw_cell_model(fig, W, H, top, box, band, cell, s) -> list[dict]:
    """One cell's model in the stacked layout: header, then the term forest on the band's shared axis."""
    fx = lambda x: x / W                                                     # noqa: E731
    fy = lambda y: 1 - y / H                                                 # noqa: E731
    x0, width = box
    tag = f"{cell['event']}|{cell['cohort']}"
    if cell["row"] is None:
        fig.text(fx(x0), fy(top), "Cox model", fontsize=6.6, fontweight="bold", ha="left", va="top")
        fig.text(fx(x0), fy(top + 0.18), cell["error"].split(": ", 1)[-1], fontsize=5.8, color=S.MUTED, ha="left",
                 va="top")
        return [dict(panel="cox_detail_not_fitted", tag=tag, message=cell["error"])]
    r = cell["row"]
    fig.text(fx(x0), fy(top), f"Cox model · {int(r.cox_n)} patients, {int(r.cox_events)} events" + _q(cell, s),
             fontsize=6.6, fontweight="bold", ha="left", va="top")
    for k, line in enumerate(_model_lines(_cell_model_text(r, band), width, 5.6)):
        fig.text(fx(x0), fy(top + 0.15 + k * MODEL_LINE), line, fontsize=5.6, color=S.INK2, ha="left", va="top")
    drawn = draw_terms(fig, W, H, top + 0.50 + band.get("stack_extra", 0.0), x0, width, cell["disp"], s.alpha,
                       BAND_PITCH, fs=5.8, tag=tag, ph_below=s.ph_note_below, q_below=s.q_mark_below,
                       xlim=band["xlim"], cols=band["cols"])
    return [dict(panel="cox_detail", cox_model=_cell_model_text(r, band), cox_n=int(r.cox_n),
                 cox_events=int(r.cox_events), **d) for d in drawn]


def _draw_band(fig, W, H, y0, band, s, n_ev) -> list[dict]:
    from matplotlib.lines import Line2D

    fx = lambda x: x / W                                                     # noqa: E731
    fy = lambda y: 1 - y / H                                                 # noqa: E731
    y0 += 0.12
    fig.add_artist(Line2D([fx(0.12), fx(W - 0.12)], [fy(y0 - 0.06)] * 2, color=S.GRID, lw=0.6))
    rows = []
    for i, cell in enumerate(band["cells"]):
        col, line = i % band["ncol"], i // band["ncol"]
        x0 = 0.12 + col * (band["width"] + 0.35)
        top = y0 + line * band["line_h"]
        head = cell["cohort"] if n_ev == 1 else f"{cell['event']} · {cell['cohort']}"
        fig.text(fx(x0), fy(top + 0.02), f"{head} · Cox model" + _q(cell, s), fontsize=7.0, fontweight="bold",
                 ha="left", va="top")
        tag = f"{cell['event']}|{cell['cohort']}"
        if cell["row"] is None:
            fig.text(fx(x0), fy(top + 0.24), cell["error"].split(": ", 1)[-1], fontsize=6.0, color=S.MUTED,
                     ha="left", va="top")
            rows.append(dict(panel="cox_detail_not_fitted", tag=tag, message=cell["error"]))
            continue
        r = cell["row"]
        text = f"{_cell_model_text(r, band)} · {int(r.cox_n)} patients, {int(r.cox_events)} events"
        for k, line in enumerate(_model_lines(text, band["width"], 5.8)):
            fig.text(fx(x0), fy(top + 0.19 + k * MODEL_LINE), line, fontsize=5.8, color=S.INK2, ha="left", va="top")
        drawn = draw_terms(fig, W, H, top + BAND_HEAD + band["head_extra"], x0, band["width"], cell["disp"], s.alpha,
                           BAND_PITCH, fs=5.9, ph_below=s.ph_note_below, q_below=s.q_mark_below,
                           tag=tag, xlim=band["xlim"], cols=band["cols"])
        rows += [dict(panel="cox_detail", cox_model=_cell_model_text(r, band), cox_n=int(r.cox_n),
                      cox_events=int(r.cox_events), **d) for d in drawn]
    fig.text(fx(0.12), fy(y0 + band["nline"] * band["line_h"] + 0.02), "Cox models: " + SCALE_NOTE,
             fontsize=5.6, color=S.MUTED, ha="left", va="top")
    return rows


# ============================================================================================ host-gene expression
GEX_HEAD = 0.48


def _expression_model(model: CoxModel) -> CoxModel:
    """The Cox model of the expression rows: expression itself plus the given model's clinical terms."""
    return CoxModel(expression=False, covariates=model.covariates, categorical=model.categorical, strata=model.strata,
                    scale=model.scale, baseline=dict(model.baseline))


def _gex_section(ds, gene, coh, endpoint, s, model, W, has_ref, marks, label_w, row_title, tab_h,
                 res=None) -> dict | None:
    """The expression rows: statistics of the host gene in each cohort given, and their layout (one row per cohort,
    as the stacked layout); None without expression or survival data for the gene. `res` may hold the gene's
    expression results already (any superset of the cohorts)."""
    if ds.expression is None or gene not in ds.expression.index or ds.survival is None:
        return None
    gm = _expression_model(model)
    coh = list(dict.fromkeys(coh))
    if res is None:
        res = analyze_expression(ds, genes=[gene], endpoints=[endpoint], cohorts=coh, settings=s, model=gm)
    cells = []
    for c in coh:
        v = res.survival[res.survival.cohort.eq(c)].iloc[0]
        t = res.cox_terms[res.cox_terms.cohort.eq(c)] if len(res.cox_terms) else pd.DataFrame()
        ok = v.cox_status == "tested"
        cells.append(dict(event=gene, cohort=c, g=res.groups[res.groups.cohort.eq(c)].iloc[0], v=v,
                          row=v if ok else None, disp=display_rows(t) if ok else [],
                          error="" if ok else f"not fitted: {not_fitted(v, res.settings)}"))
    L = _stacked_layout(W, has_ref, marks, label_w)
    n_terms = max([len(x["disp"]) for x in cells] + [3])
    describe = gm.describe(False).replace("PSI", "expression", 1)
    extra = _head_extra(cells, L["model"][1], 5.6, lambda c: _cell_model_text(c["row"], {"describe": describe}), 2)
    if _any_km_ph([c["v"].get("km_ph_p") for c in cells], s):
        row_title += 0.10                                   # a line for a KM proportional-hazards note
    row_h = max(row_title + 1.02 + tab_h, 0.50 + extra + BAND_PITCH * n_terms + BAND_AXIS + 0.10)
    shown = [d["r"] for x in cells for d in x["disp"] if d["r"] is not None]
    smp = ds.samples[ds.samples.cohort.isin(coh) & ds.samples.role.isin(["case", "reference"])]
    vals = ds.expression.loc[gene].reindex(smp.sample_id).to_numpy(float)
    vals = vals[np.isfinite(vals)]
    lo, hi = (float(vals.min()), float(vals.max())) if len(vals) else (0.0, 1.0)
    pad = 0.06 * (hi - lo) if hi > lo else 0.5
    return dict(gene=gene, model=gm, res=res, cells=cells, L=L, row_h=row_h, lim=(lo - pad, hi + pad),
                xlim=_range(shown) if shown else None, describe=describe, stack_extra=extra, row_title=row_title,
                cols=terms_columns(L["model"][1], [x["disp"] for x in cells if x["disp"]] or [[]], 5.8),
                height=GEX_HEAD + len(coh) * row_h + 0.12)


def _split_words(rule) -> str:
    """'its median', 'its mean' or '4.5 (set)' for a KM split rule."""
    return f"its {rule}" if isinstance(rule, str) else f"{rule:g} (set)"


def expression_panel(ds: Dataset, gene: str, cohorts, endpoint: str, *, settings: Settings | None = None,
                     model: CoxModel | None = None, out_dir=None, stem: str | None = None,
                     part: tuple | None = None, endpoint_label: str | None = None) -> Panel | None:
    """The host gene's own expression on one page, drawn once per gene (after its splicing pages in a probe).

    A forest of every cohort with a test (case-vs-reference delta of expression, and the Cox HR per SD of expression),
    then one row per cohort given: group view, KM split, and Cox on expression plus `model`'s clinical terms (the
    tested term is expression itself). None when the gene has no expression values or there is no survival table.
    Expression is not adjusted for multiple testing.
    """
    import matplotlib
    from matplotlib.figure import Figure
    from matplotlib.lines import Line2D
    from matplotlib.patches import Rectangle

    s = settings or Settings()
    model = model or CoxModel()
    if ds.expression is None or gene not in ds.expression.index or ds.survival is None:
        return None
    cohorts = [c for c in dict.fromkeys(_as_list(cohorts)) if c in set(ds.cohorts)]
    if not cohorts:
        raise InputError(f"{gene}: none of the cohorts is in the data")
    gm = _expression_model(model)
    res = analyze_expression(ds, genes=[gene], endpoints=[endpoint], settings=s, model=gm)
    gs = res.groups.set_index("cohort")
    vs = res.survival.set_index("cohort")
    fcoh = [c for c in ds.cohorts if c in set(cohorts) or (c in vs.index and (
        _get(vs.loc[c], "cox_status") in ("tested", "failed") or _get(gs.loc[c], "unpaired_status") == "tested"))]
    glab = dict(case=s.case_label or ds.labels["case"], reference=s.reference_label or ds.labels["reference"])
    ylabel = endpoint_label or ENDPOINT_NAMES.get(endpoint, endpoint)
    has_ref = ds.has_reference
    W = 7.2
    label_w = max(S.text_width(c, 6.0, weight="bold") for c in fcoh)
    row_title, tab_h = (0.52, 0.66) if has_ref else (0.46, 0.62)
    gx = _gex_section(ds, gene, cohorts, endpoint, s, model, W, has_ref, {}, label_w, row_title, tab_h, res=res)
    L = _stacked_layout(W, has_ref, {}, label_w)
    recs = {}
    for c in fcoh:
        g_, v_ = gs.loc[c], vs.loc[c]
        ok = _get(v_, "cox_status") == "tested"
        d = dict(event_id=gene, mark="", shade=c in set(cohorts), cox_status=_get(v_, "cox_status"))
        for design in ("paired", "unpaired"):
            st = _get(g_, f"{design}_status")
            d[f"{design}_status"] = st
            d[f"{design}_delta"] = float(_get(g_, f"{design}_delta_median")) if st == "tested" else NAN
            d[f"{design}_p"] = float(_get(g_, f"{design}_p")) if st == "tested" else NAN
        for k_out, k_in in (("hr", "hr_per_sd"), ("lo", "ci_low_sd"), ("hi", "ci_high_sd"), ("cox_p", "cox_p"),
                            ("ph_p", "ph_p")):
            d[k_out] = float(_get(v_, k_in)) if ok else NAN
        d["cox_events"] = _get(v_, "cox_events")
        recs[c] = [d]
    top, sub_h = 0.30, 0.20
    forest_h = 0.106 * len(fcoh)
    forest_block = 0.25 + forest_h + 0.85
    first = []
    second = []
    paired_drawn = has_ref and any(_get(gs.loc[c], "paired_status") == "tested" for c in cohorts if c in gs.index)
    cl = S.in_sentence(glab["case"])
    if paired_drawn:
        second += [("up", f"Higher in {cl}"), ("down", f"Lower in {cl}")]
    if has_ref:
        second += [("split", "KM split"), ("pd", "Paired Δ"), ("ud", "Unpaired Δ")]
    second += [("hs", f"Cox p < {s.alpha:g}"), ("hn", f"Cox p ≥ {s.alpha:g}"), ("band", "Cohort drawn below")]
    if gx and gx["cells"]:                                  # the model rows' other terms
        second += [("ts", f"Model term p < {s.alpha:g}"), ("tn", f"p ≥ {s.alpha:g}")]
    low = (lambda v: s.ph_note_below > 0 and v is not None and np.isfinite(float(v)) and float(v) < s.ph_note_below)
    if any(low(r[0]["ph_p"]) for r in recs.values()) or (gx and any(
            low(_get(cell["v"], "km_ph_p")) or any(d["r"] is not None and low(getattr(d["r"], "ph_p", None))
                                                   for d in cell["disp"]) for cell in gx["cells"])):
        second.append(("ph", f"non-proportional hazards (p < {s.ph_note_below:g})"))
    if any(np.isfinite(r[0]["hr"]) and F.low_power(r[0]["cox_events"], s.cox_low_power_events) for r in recs.values()):
        second.append(("low", f"fewer than {s.cox_low_power_events} events: low power"))
    lines = _legend_lines(W, first, second)
    lines = [ln for ln in lines if ln]
    note = (f"Host-gene expression of {gene} on the expression table's scale; not adjusted for multiple testing. A "
            f"{S.in_sentence(glab['case'])}–{S.in_sentence(glab['reference'])} hit needs p < {s.alpha:g} and "
            f"|Δ median| > {s.gex_min_abs_delta:g}. "
            f"The HR is per SD of expression in each cohort, with {', '.join(gm.covariates) or 'no other terms'}."
            + (f" Settings changed from the defaults: {s.changed_text()}." if s.changed() else ""))
    note_lines = PR._wrap(note, W - 0.26, 5.4)
    legend_h = 0.42 + 0.15 * (len(lines) - 1) + 0.05 + 0.115 * len(note_lines)
    rows_h = gx["height"] if gx else 0.0
    H = top + sub_h + forest_block + rows_h + legend_h
    fx = lambda x: x / W                                                     # noqa: E731
    fy = lambda y: 1 - y / H                                                 # noqa: E731
    rows: list[dict] = []
    with matplotlib.rc_context(S.rc()):
        fig = Figure(figsize=(W, H))
        title = f"{gene} expression"
        fig.text(fx(0.12), fy(0.10), gene, fontsize=8.5, fontweight="bold", style="italic", ha="left", va="top")
        tw = 0.12 + S.text_width(gene, 8.5, style="italic", weight="bold") + 0.05
        sub_ = (ds.notes or {}).get("subset")
        tail = [endpoint] + ([f"subset: {_subset_name(sub_)}"] if sub_ else []) + (
            [f"page {part[0]} of {part[1]}"] if part else [])
        rest = " · ".join(["expression", _and(cohorts)] + tail)
        if tw + S.text_width(rest, 8.5, weight="bold") > W - 0.15:
            rest = " · ".join(["expression", f"{len(cohorts)} cohorts"] + tail)
        fig.text(fx(tw), fy(0.10), rest, fontsize=8.5, fontweight="bold", ha="left", va="top")
        fig.text(fx(0.12), fy(0.30), f"The forest shows {gene} expression in every cohort with a test; the rows below "
                 "show the cohorts drawn on its splicing pages", fontsize=6.0, color=S.INK2, ha="left", va="top")
        # forest of every cohort
        f_top = top + sub_h
        fy0, fh = fy(f_top + 0.25 + forest_h), forest_h / H
        axL = fig.add_axes([fx(L["axL"][0]), fy0, fx(L["axL"][1]), fh]) if has_ref else None
        axR = fig.add_axes([fx(L["axR"][0]), fy0, fx(L["axR"][1]), fh])
        lim, xlo, xhi, n_clip = F.axis_limits(recs)
        F.set_axes(axL, axR, max(lim, 1.2 * s.gex_min_abs_delta), xlo, xhi, s.gex_min_abs_delta, glab,
                   f"Cox: {gm.describe(False).replace('PSI', 'expression', 1)}"
                   + ridge_text(s, bool(gm.covariates), ["expression"]), quantity="expression",
                   hr_label="HR per SD of\nexpression (95% CI)")
        F.draw(fig, axL, axR, fx(L["lab_x"]), fx(L["mark_x"]), fcoh, recs, 1, s.alpha, "Cox",
               ph_below=s.ph_note_below, low_below=s.cox_low_power_events)
        rows.append(dict(panel="gex_forest_axis", hr_axis_low=xlo, hr_axis_high=xhi, n_ci_clipped=n_clip))
        for c in fcoh:
            r = recs[c][0]
            rows.append(dict(panel="gex_forest", endpoint=endpoint, gene=gene, cohort=c, shade=r["shade"],
                             paired_status=r["paired_status"], paired_delta_median=r["paired_delta"],
                             paired_p=r["paired_p"], unpaired_status=r["unpaired_status"],
                             unpaired_delta_median=r["unpaired_delta"], unpaired_p=r["unpaired_p"],
                             cox_status=r["cox_status"], hr_per_sd=r["hr"], ci_low_sd=r["lo"], ci_high_sd=r["hi"],
                             cox_p=r["cox_p"], ph_p=r["ph_p"], cox_events=r["cox_events"],
                             low_power_marked=bool(np.isfinite(r["hr"]))
                             and F.low_power(r["cox_events"], s.cox_low_power_events)))
        # one row per cohort drawn in full
        if gx:
            rows += _draw_gex(fig, W, H, f_top + forest_block, gx, ds, endpoint, ylabel, s, glab, has_ref,
                              row_title, tab_h)
        # legend and footnote
        for i, line in enumerate(note_lines):
            fig.text(fx(0.12), fy(H - 0.08 - 0.115 * (len(note_lines) - 1 - i)), line, fontsize=5.4, color=S.MUTED,
                     ha="left", va="bottom")
        for li, items in enumerate(lines):
            x = 0.12
            yy = fy(H - legend_h + 0.10 + li * 0.15)
            for kind, text in items:
                if kind in ("up", "down"):
                    fig.add_artist(Line2D([fx(x), fx(x + 0.16)], [yy, yy], color=S.UP if kind == "up" else S.DOWN,
                                          lw=1.1))
                elif kind == "split":
                    fig.add_artist(Line2D([fx(x), fx(x + 0.16)], [yy, yy], color=S.INK, lw=0.6, ls=(0, (2, 1.5))))
                elif kind == "band":
                    fig.add_artist(Rectangle((fx(x), yy - 0.05 / H), fx(0.16), 0.10 / H, facecolor=S.PAPER, lw=0))
                elif kind == "ph":
                    fig.text(fx(x + 0.08), yy, F.PH_MARK, fontsize=7, color=S.INK, ha="center", va="center")
                elif kind == "low":
                    fig.text(fx(x + 0.08), yy, F.LOW_MARK, fontsize=7, color=S.INK, ha="center", va="center")
                else:
                    mk = "o" if kind in ("pd", "ud") else "s" if kind in ("ts", "tn") else "D"
                    fig.add_artist(Line2D([fx(x + 0.08)], [yy], marker=mk, ms=3.6 if mk == "o" else 3.2, ls="",
                                          markerfacecolor=S.INK if kind in ("pd", "hs", "ts") else "white",
                                          markeredgecolor=S.INK, mew=0.7))
                fig.text(fx(x + 0.21), yy, text, fontsize=6, color=S.INK2, ha="left", va="center")
                x += 0.21 + S.text_width(text, 6) + 0.2
    table = pd.DataFrame(rows)
    named = ([_safe(c) for c in cohorts] if len(cohorts) <= 4 else
             [f"{len(cohorts)}cohorts", hashlib.sha1(",".join(cohorts).encode()).hexdigest()[:6]])
    stem = stem or "_".join([_safe(gene), "expression"] + named + [_safe(endpoint)]
                            + ([f"p{part[0]}of{part[1]}"] if part else []))
    sv = ds.survival[ds.survival.endpoint.eq(endpoint)]
    cl = None
    if ds.clinical is not None and gm.clinical_columns:
        cl = ds.clinical[[c for c in gm.clinical_columns if c in ds.clinical.columns]]
    prov = P.record("expression_panel", s, dict(expression=ds.expression.loc[[gene]], samples=ds.samples,
                                                  pairs=ds.pairs, survival=sv, clinical=cl),
                    dict(gene=gene, cohorts=cohorts, endpoint=endpoint, model=gm.to_dict(), title=title))
    paths = {}
    if out_dir is not None:
        paths = S.save(fig, out_dir, stem, s.formats, s.dpi, table, prov)
    return Panel(figure=fig, table=table, paths=paths, provenance=prov, stem=stem)


def _gex_row(r: dict) -> dict:
    """A plotted-values row of the expression section: panels prefixed gex_, PSI columns named value."""
    r = dict(r, panel=f"gex_{r['panel']}")
    for k in ("psi_reference", "psi_case", "psi"):
        if k in r:
            r[k.replace("psi", "value")] = r.pop(k)
    return r


def _draw_gex(fig, W, H, y0, gx, ds, endpoint, ylabel, s, glab, has_ref, row_title, tab_h) -> list[dict]:
    from matplotlib.lines import Line2D

    fx = lambda x: x / W                                                     # noqa: E731
    fy = lambda y: 1 - y / H                                                 # noqa: E731
    gene, L = gx["gene"], gx["L"]
    row_title = gx.get("row_title", row_title)              # with the expression rows' own PH line
    fig.add_artist(Line2D([fx(0.12), fx(W - 0.12)], [fy(y0 + 0.06)] * 2, color=S.GRID, lw=0.6))
    fig.text(fx(0.12), fy(y0 + 0.12), "Host-gene expression", fontsize=7.6, fontweight="bold", ha="left", va="top")
    fig.text(fx(0.12 + S.text_width("Host-gene expression", 7.6, weight="bold") + 0.06), fy(y0 + 0.12),
             f"· {gene} itself: {S.in_sentence(glab['case'])} vs {S.in_sentence(glab['reference'])}, a KM split at "
             f"{_split_words(s.km_split_expression)}, and Cox on {gx['describe']}", fontsize=6.0, color=S.INK2,
             ha="left", va="top")
    v_all = ds.expression.loc[gene]
    rows = []
    for ri, cell in enumerate(gx["cells"]):  # noqa: B007
        c = cell["cohort"]
        y_row = y0 + GEX_HEAD + ri * gx["row_h"]
        ax_h = gx["row_h"] - row_title - tab_h
        tag = f"{gene}|{c}|expression"
        fig.text(fx(0.12), fy(y_row + 0.02), f"{c} · {gene} expression", fontsize=7.2, fontweight="bold", ha="left",
                 va="top")
        tmp, printed = [], []
        if has_ref:
            ax = fig.add_axes([fx(L["tn"][0]), fy(y_row + row_title + ax_h), fx(L["tn"][1]), ax_h / H])
            _group_view(fig, ax, ds, gene, c, v_all, cell["g"], cell["v"], s, glab, tmp, printed, tag,
                        quantity="Expression", lim=gx["lim"])
        ax_km = fig.add_axes([fx(L["km"][0]), fy(y_row + row_title + ax_h), fx(L["km"][1]), ax_h / H])
        _km_view(fig, ax_km, ds, gene, c, v_all, cell["v"], endpoint, ylabel, s, tmp, printed, tag, what="expression")
        tmp += _draw_cell_model(fig, W, H, y_row + 0.06, L["model"], gx, cell, s)
        rows += [_gex_row(r) for r in tmp + printed]
        rows.append(dict(panel="gex_cell", tag=tag, gene=gene, cohort=c, endpoint=endpoint,
                         paired_p=_get(cell["g"], "paired_p"), unpaired_p=_get(cell["g"], "unpaired_p"),
                         km_p=_get(cell["v"], "km_p"), cox_p=_get(cell["v"], "cox_p"),
                         hr_per_sd=_get(cell["v"], "hr_per_sd")))
    return rows


# ============================================================================================ the two left panels
def _group_view(fig, ax, ds, e, c, v, t_, v_, s, glab, rows, printed, tag, quantity="PSI", lim=None):
    st_p, st_u = t_["paired_status"], t_["unpaired_status"]
    pr = ds.pairs[ds.pairs.cohort.eq(c)]
    pr_r, pr_c = v.reindex(pr.reference_sample).to_numpy(float), v.reindex(pr.case_sample).to_numpy(float)
    ok = np.isfinite(pr_r) & np.isfinite(pr_c)
    n_pairs = int(ok.sum())
    smp = ds.samples[ds.samples.cohort.eq(c)]
    case, ref = smp.role.eq("case"), smp.role.eq("reference")
    xt = patient_values(v, smp.sample_id[case], smp.patient_id[case])      # one value per patient, as tested
    yt = patient_values(v, smp.sample_id[ref], smp.patient_id[ref])
    x, y = xt.value.to_numpy(float), yt.value.to_numpy(float)
    cut = float(_get(v_, "cutoff"))
    if st_u == "tested":
        _check(len(x) == int(t_["unpaired_n_case"]) and len(y) == int(t_["unpaired_n_reference"]), "group sizes", e, c)
    qp, qu = float(_get(t_, "paired_q")), float(_get(t_, "unpaired_q"))
    all_stats = (f"Δ {S.fd(t_['unpaired_delta_median'])}\np {S.fp(t_['unpaired_p'])}" if st_u == "tested"
                 else "not tested\n") + ("\n" + S.q_text(qu, s.q_mark_below) if st_u == "tested" and np.isfinite(qu)
                                          else "")
    if st_p == "tested":
        _check(n_pairs == int(t_["paired_n_pairs"]), "pair count", e, c)
        _check(abs(np.median(pr_c[ok] - pr_r[ok]) - t_["paired_delta_median"]) < 1e-12, "paired delta", e, c)
        (mn, mt), meds = T.paired_and_all(
            fig, ax, pr_r[ok], pr_c[ok], y, x, cut, glab,
            f"{n_pairs} pairs\nΔ {S.fd(t_['paired_delta_median'])}\np {S.fp(t_['paired_p'])}"
            + ("\n" + S.q_text(qp, s.q_mark_below) if np.isfinite(qp) else ""),
            f"all samples\n{all_stats}", quantity=quantity, lim=lim)
        rows += [dict(panel="pairs", tag=tag, patient_id=pid, psi_reference=a, psi_case=b)
                 for pid, a, b in zip(pr.patient_id.to_numpy()[ok], pr_r[ok], pr_c[ok])]
        rows.append(dict(panel="pairs_median", tag=tag, psi_reference=mn, psi_case=mt))
        printed += [dict(panel="printed", tag=tag, what="pairs", value=n_pairs, source="groups.paired_n_pairs"),
                    dict(panel="printed", tag=tag, what="paired delta median", value=t_["paired_delta_median"],
                         source="groups.paired_delta_median"),
                    dict(panel="printed", tag=tag, what="paired p (exact signed-rank)", value=t_["paired_p"],
                         source="groups.paired_p"),
                    dict(panel="printed", tag=tag, what="paired q (BH within gene)", value=qp, source="groups.paired_q")]
    elif st_u == "tested":
        if st_p == "all_differences_zero":
            note = "all within-patient\ndifferences are zero"
        elif n_pairs == 0 and len(pr) == 0:
            note = "no matched pairs:\nno within-patient test"
        else:
            note = f"{n_pairs} pairs (< {s.min_pairs}):\nno within-patient test"
        meds = T.all_samples(fig, ax, y, x, cut, glab, f"all samples\n{all_stats}", note, quantity=quantity, lim=lim)
        printed.append(dict(panel="printed", tag=tag, what="pairs available", value=n_pairs,
                            source="groups.paired_n_pairs"))
    else:
        reason = {"no_reference_samples": f"no {S.in_sentence(glab['reference'])} samples",
                  "constant": f"{quantity} constant in both groups"}.get(
            st_u, f"fewer than {max(s.min_group, 1)} patient{'s' if s.min_group > 1 else ''} per group")
        if st_u == "too_few_samples" and len(y) == 0:
            reason = f"no {S.in_sentence(glab['reference'])} values"
        T.empty(ax, f"no {S.in_sentence(glab['case'])} vs {S.in_sentence(glab['reference'])}\ntest: {reason}",
                quantity=quantity,
                lim=lim)
        rows.append(dict(panel="groups_not_tested", tag=tag, paired_status=st_p, unpaired_status=st_u))
        return
    rows += [dict(panel="all_samples", tag=tag, group=g, patient_id=pid, n_samples=int(k), psi=float(val))
             for g, tab in (("reference", yt), ("case", xt)) for pid, val, k in tab.itertuples(index=False)]
    rows.append(dict(panel="all_samples_median", tag=tag, psi_reference=meds["reference"], psi_case=meds["case"]))
    if st_u == "tested":
        printed += [dict(panel="printed", tag=tag, what="reference patients", value=len(y),
                         source="groups.unpaired_n_reference"),
                    dict(panel="printed", tag=tag, what="case patients", value=len(x), source="groups.unpaired_n_case"),
                    dict(panel="printed", tag=tag, what="unpaired delta median", value=t_["unpaired_delta_median"],
                         source="groups.unpaired_delta_median"),
                    dict(panel="printed", tag=tag, what="unpaired p (Mann-Whitney)", value=t_["unpaired_p"],
                         source="groups.unpaired_p"),
                    dict(panel="printed", tag=tag, what="unpaired q (BH within gene)", value=qu,
                         source="groups.unpaired_q")]


def _km_view(fig, ax, ds, e, c, v, v_, endpoint, ylabel, s, rows, printed, tag, what="PSI"):
    if v_["km_status"] != "tested":
        K.empty(ax, KM_MESSAGE.get(v_["km_status"], str(v_["km_status"])), ylabel)
        rows.append(dict(panel="km_not_tested", tag=tag, km_status=v_["km_status"]))
        return
    base = ds.samples[ds.samples.cohort.eq(c) & ds.samples.survival_cohort]
    sv = ds.survival[ds.survival.endpoint.eq(endpoint)].set_index("patient_id")
    base = base[base.patient_id.isin(sv.index)]
    x = v.reindex(base.sample_id).to_numpy(float)
    ok = np.isfinite(x)
    pid = base.patient_id.to_numpy()[ok]
    x, t, ev = x[ok], sv.time.reindex(pid).to_numpy(float), sv.event.reindex(pid).to_numpy(int)
    high = km_high(x, v_["cutoff"], ties_high(v_))                       # above the cutoff, or at it (see km_cut)
    _check(high.sum() == v_["n_high"] and (~high).sum() == v_["n_low"], "KM arms", e, c)
    _check(ev[high].sum() == v_["events_high"] and ev[~high].sum() == v_["events_low"], "KM events", e, c)
    K.draw(fig, ax, s.years(t), ev, high, v_, ylabel, rows, tag, s, what=what)
    printed += [dict(panel="printed", tag=tag, what=f"KM split ({what})", value=v_["cutoff"], source="survival.cutoff"),
                dict(panel="printed", tag=tag, what="log-rank HR high vs low", value=v_["logrank_hr"],
                     source="survival.logrank_hr"),
                dict(panel="printed", tag=tag, what="log-rank p", value=v_["km_p"], source="survival.km_p"),
                dict(panel="printed", tag=tag, what="log-rank q (BH within gene)", value=_get(v_, "km_q"),
                     source="survival.km_q"),
                dict(panel="printed", tag=tag, what="KM proportional-hazards p", value=_get(v_, "km_ph_p"),
                     source="survival.km_ph_p"),
                dict(panel="printed", tag=tag, what="events high arm", value=int(v_["events_high"]),
                     source="survival.events_high"),
                dict(panel="printed", tag=tag, what="events low arm", value=int(v_["events_low"]),
                     source="survival.events_low"),
                dict(panel="km_arms", tag=tag, n_low=int(v_["n_low"]), n_high=int(v_["n_high"]),
                     events_low=int(v_["events_low"]), events_high=int(v_["events_high"]))]


def _protein_changes(proteins, ds, events, geoms) -> dict:
    """{event: ProteinChange} from a cache, a cache folder, or a ready dict; {} without proteins."""
    if proteins is None:
        return {}
    if isinstance(proteins, dict):
        return proteins
    from ..protein import ProteinCache, protein_change
    cache = proteins if isinstance(proteins, ProteinCache) else ProteinCache.load(proteins)
    return {e: protein_change(cache, e, ds.events.loc[e].to_dict(), g) for e, g in zip(events, geoms)}


def _protein_source(proteins):
    if proteins is None or isinstance(proteins, dict):
        return None if proteins is None else "precomputed"
    return getattr(proteins, "source", None) or str(proteins)


def _check(ok: bool, what: str, e: str, c: str):
    if not ok:
        raise RuntimeError(f"internal inconsistency in {what} for {e} in {c}; please report this")


def _provenance(ds, events, family, endpoint, s, model, detail_model, highlight, call) -> dict:
    """The page's record: `model` is the forest's model (call.model), `detail_model` the model rows' (call.detail_model,
    when shown). Hashed: the rows of every input the drawn numbers depend on, including the PSI of the other events of
    the q families (`family`: the gene's events of the same quantity), whose tests set the printed q values."""
    ev = ds.events.loc[events]
    others = [e for e in family if e not in set(events) and e in ds.psi.index]
    ex = None
    if ds.expression is not None:
        hosts = ds.events.loc[list(events) + others, "expression_gene"].unique()
        ex = ds.expression.loc[[g for g in hosts if g in ds.expression.index]]
    sv = None if ds.survival is None else ds.survival[ds.survival.endpoint.eq(endpoint)]
    cols = list(dict.fromkeys(model.clinical_columns + ([] if detail_model is None else detail_model.clinical_columns)))
    cl = None
    if ds.clinical is not None and cols:
        cl = ds.clinical[[c for c in cols if c in ds.clinical.columns]]
    hl = highlight if isinstance(highlight, pd.DataFrame) else (None if highlight is None else pd.DataFrame(
        [dict(event_id=k[0], cohort=k[1], label=v) for k, v in highlight.items()]) if isinstance(highlight, dict)
        else read_table(highlight))
    return P.record("event_panel", s, dict(psi=ds.psi.loc[events], samples=ds.samples, pairs=ds.pairs, survival=sv,
                                           events=ev, expression=ex, clinical=cl, highlight=hl,
                                           psi_q_family=ds.psi.loc[others] if others else None),
                    dict(call, model=model.to_dict()))
