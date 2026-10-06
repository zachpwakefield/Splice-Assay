"""The event schematic: collapsed host-gene model, nested small RNAs, and one row per event on one genomic scale.

Minus-strand genes are drawn 5' to 3' (coordinates decrease left to right). A gene much longer than the events is
drawn only around them (a window), with marks where it continues. On each event row the constant exons are
grey, the region PSI measures is coloured by event type (for MXE the other exon is outlined), the arcs above the row
are the junctions of the form PSI counts and the arc below is the other form. Regions shorter than 25 nt get a caret.
The region's length is written above it, lifted above a junction arc that would cross the label.
"""
from __future__ import annotations

import numpy as np
from matplotlib import transforms
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path as MPath

from ..events import PSI_MEANING
from . import style as S

TINY = 25                        # variable regions shorter than this get a caret
FIRST_ROW = {True: 0.84, False: 0.42}   # y (inches from the schematic top) of the first event row, with/without gene
PITCH = 0.42


def height(n_events: int, has_model: bool) -> float:
    return FIRST_ROW[has_model] + PITCH * (n_events - 1) + 0.19


def arc(ax, x0, x1, y, h, **kw):
    """Junction arc from x0 to x1 at baseline y with apex h (data units; y grows downwards)."""
    verts = [(x0, y), ((x0 + x1) / 2, y - 2 * h), (x1, y)]
    ax.add_patch(PathPatch(MPath(verts, [MPath.MOVETO, MPath.CURVE3, MPath.CURVE3]), fill=False, **kw))


def arc_rise(arcs, lo: float, hi: float, h: float) -> float:
    """How far the junction arcs (apex h) rise above their baseline anywhere over [lo, hi] (data units): an arc from
    a to b rises 4 h t (1 - t) at t = (x - a) / (b - a)."""
    top = 0.0
    for a, b in arcs:
        a, b = min(a, b), max(a, b)
        s, e = max(lo, a), min(hi, b)
        if s >= e:
            continue
        m = (a + b) / 2
        t = ((m if s <= m <= e else (e if e < m else s)) - a) / (b - a)
        top = max(top, 4 * h * t * (1 - t))
    return top


def box(ax, x0, x1, yc, h, **kw):
    ax.add_patch(Rectangle((x0, yc - h / 2), x1 - x0, h, **{"linewidth": 0, **kw}))


def _kb(n: float) -> str:
    return f"{n / 1000:.0f} kb" if n >= 10000 else (f"{n / 1000:.1f} kb" if n >= 1000 else f"{int(n)} nt")


def scale_bar(span: float) -> tuple[int, str]:
    if span < 1500:
        bar = 100
    elif span < 3000:
        bar = 250
    elif span < 6000:
        bar = 500
    else:
        bar = max([c for c in (1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000)
                   if c <= span / 6] or [1000])
    return bar, (f"{bar / 1000:g} kb" if bar >= 1000 else f"{bar} nt")


def window_for(events_span, model, min_gene=20000, ratio=3.0, margin=0.12):
    """A drawing window around the events when the gene is much longer than them (else None)."""
    if model is None:
        return None
    lo, hi = events_span
    g_lo, g_hi = model.blocks[0][0], model.blocks[-1][1]
    if g_hi - g_lo <= max(ratio * (hi - lo), min_gene):
        return None
    pad = max(margin * (hi - lo), 300)
    return int(lo - pad), int(hi + pad)


def draw(ax, fig, events, chrom: str, strand: str, model, rows: list, nested_label: str = "snoRNAs", window=None):
    """events: [dict(event_id, label, geometry)]; model: annotation.GeneModel or None; window: (start, end) to draw
    only part of a long gene. Appends plotted values."""
    H = ax.get_position().height * fig.get_figheight()
    pieces = [p for ev in events for p in ev["geometry"].constant + ev["geometry"].variable + ev["geometry"].other]
    blocks = model.blocks if model is not None else []
    gene_span = (blocks[0][0], blocks[-1][1]) if blocks else None
    if window is not None and blocks:
        blocks = [(max(a, window[0]), min(b, window[1]), n) for a, b, n in blocks if b > window[0] and a < window[1]]
    lo = min([b[0] for b in blocks] + [p[0] for p in pieces] + ([window[0]] if window else []))
    hi = max([b[1] for b in blocks] + [p[1] for p in pieces] + ([window[1]] if window else []))
    pad = 0.015 * (hi - lo)
    x0, x1 = lo - pad, hi + pad
    ax.set_xlim((x1, x0) if strand == "-" else (x0, x1))                   # 5' on the left
    ax.set_ylim(H, 0)
    ax.axis("off")
    lab = transforms.blended_transform_factory(ax.transAxes, ax.transData)
    # header: drawn span (1-based) and strand; scale bar at the right end
    a0, a1 = (lo, hi) if window else ((blocks[0][0], blocks[-1][1]) if blocks else (lo, hi))
    note = "minus strand, drawn 5′→3′" if strand == "-" else "plus strand"
    if window and gene_span:
        note += f" · view {_kb(a1 - a0)} of a {_kb(gene_span[1] - gene_span[0])} gene"
    ax.text(0, 0.07, f"{chrom}:{a0 + 1:,}–{a1:,} · {note}", transform=lab, ha="left", va="center", fontsize=6.2,
            color=S.INK2)
    bar, bar_text = scale_bar(hi - lo)
    xs = x1 - pad if strand != "-" else x0 + pad
    sgn = -1 if strand != "-" else 1
    ax.plot([xs, xs + sgn * bar], [0.07, 0.07], color=S.INK2, lw=0.8, solid_capstyle="butt")
    ax.text(xs + sgn * bar / 2, 0.02, bar_text, ha="center", va="bottom", fontsize=6, color=S.INK2)
    rows.append(dict(panel="schematic_header", chrom=chrom, start_1based=a0 + 1, end=a1, strand=strand,
                     scale_bar_nt=bar, windowed=bool(window)))
    # gene track and nested small RNAs
    if model is not None:
        yg = 0.27
        g0 = lo if window and gene_span[0] < lo else (blocks[0][0] if blocks else lo)
        g1 = hi if window and gene_span[1] > hi else (blocks[-1][1] if blocks else hi)
        ax.plot([g0, g1], [yg, yg], color=S.EXON_LINE, lw=0.6, zorder=1)
        if window:                                   # the gene continues beyond the window
            for side, x_ in (("lo", lo), ("hi", hi)):
                beyond = gene_span[0] < lo if side == "lo" else gene_span[1] > hi
                if beyond:
                    left = (side == "lo") == (strand != "-")      # drawn at the left edge?
                    ax.text(x_, yg, "‹" if left else "›", ha="center", va="center", fontsize=9, color=S.EXON_LINE)
        for a, b, n in blocks:
            box(ax, a, b, yg, 0.13, facecolor=S.EXON, zorder=2)
            rows.append(dict(panel="schematic_gene_block", gene=model.gene, start_0based=a, end=b, n_transcript_keys=n))
        ax.text(-0.005, yg, model.gene, transform=lab, ha="right", va="center", fontsize=7, style="italic",
                color=S.INK)
        ys = 0.445
        nested = model.nested if not window else model.nested[(model.nested.end > lo) & (model.nested.start < hi)]
        for r in nested.itertuples():
            box(ax, r.start, r.end, ys, 0.075, facecolor=S.NESTED, zorder=2)
            ax.text((r.start + r.end) / 2, ys + 0.06, r.gene_name, ha="center", va="top", fontsize=5.8,
                    style="italic", color=S.INK)
            rows.append(dict(panel="schematic_nested_gene", gene=r.gene_name, start_0based=r.start, end=r.end))
        if len(nested):
            ax.text(-0.005, ys, nested_label, transform=lab, ha="right", va="center", fontsize=6.2, color=S.INK2)
        rows.append(dict(panel="schematic_model_rule", gene=model.gene, n_transcripts_used=model.n_transcripts,
                         min_transcripts=model.min_transcripts))
    # events
    for i, ev in enumerate(events):
        g = ev["geometry"]
        y = FIRST_ROW[model is not None] + i * PITCH
        col = S.TYPE_COLOR.get(g.event_type, S.TYPE_DEFAULT)
        pp = g.constant + g.variable + g.other
        ax.plot([min(p[0] for p in pp), max(p[1] for p in pp)], [y, y], color=S.EXON_LINE, lw=0.5, zorder=1)
        for a, b in g.constant:
            box(ax, a, b, y, 0.13, facecolor=S.EXON, zorder=2)
        for a, b in g.variable:
            box(ax, a, b, y, 0.07 if g.event_type == "RI" else 0.13, facecolor=col, zorder=3)
            if b - a < TINY:
                ax.plot([(a + b) / 2], [y - 0.115], marker="v", ms=3.0, color=col, mew=0, zorder=4)
        for a, b in g.other:
            box(ax, a, b, y, 0.13 - 0.02, facecolor="white", edgecolor=col, linewidth=0.8, zorder=3)
        for a, b in g.psi_arcs:
            arc(ax, a, b, y - 0.07, 0.075, edgecolor=col, lw=0.8, zorder=2)
        for a, b in g.other_arcs:
            arc(ax, a, b, y + 0.07, -0.075, edgecolor=S.EXON_LINE, lw=0.8, zorder=2)
        vlen = g.variable_length
        vc = np.mean([(a + b) / 2 for a, b in g.variable])
        text, ly = f"{vlen} nt", y - 0.12 - (0.03 if vlen < TINY else 0)          # ly: the label's bottom
        half = S.text_width(text, 5.8) / 2 * (x1 - x0) / (ax.get_position().width * fig.get_figwidth())
        top = y - 0.07 - arc_rise(g.psi_arcs, vc - half, vc + half, 0.075)        # the arcs' top under the label
        if top < ly + 0.015:                             # a junction arc rises into (or up to) the label: lift it
            ly = top - 0.015
        ax.text(vc, ly, text, ha="center", va="bottom", fontsize=5.8, color=col)
        ax.text(-0.005, y - 0.045, ev["label"], transform=lab, ha="right", va="center", fontsize=6.6, color=col,
                fontweight="bold")
        ax.text(-0.005, y + 0.085, PSI_MEANING.get(g.event_type, "PSI"), transform=lab, ha="right", va="center",
                fontsize=5.8, color=S.INK2)
        for part, lst in (("constant_exon", g.constant), ("variable_region", g.variable), ("other_exon", g.other),
                          ("junction_psi_form", g.psi_arcs), ("junction_other_form", g.other_arcs)):
            for a, b in lst:
                rows.append(dict(panel="schematic_event", event_id=ev["event_id"], part=part, start_0based=a, end=b))
