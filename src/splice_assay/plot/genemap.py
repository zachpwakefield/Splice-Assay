"""The probe's gene map: one figure per gene with its collapsed model and every probed event, 5' to 3', each row
beside its Cox result in every cohort (the overview's cells; the gene's own row holds its expression's) and, with
correlations, the median over the cohorts of the Spearman rho between the rows (a lower triangle).

The event rows are drawn as on the assay pages: constant exons grey, the region PSI measures in the event type's
colour (for MXE the other exon outlined), the junctions of the PSI form arched above and of the other form below.
Pale columns carry the gene's exons down through the event rows, so an event's exons line up with the gene's. Long
introns, and exons far longer than the gene's typical exon, are drawn shortened (see `squeeze`).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from ..config import Settings
from ..dataset import Dataset, InputError
from ..events import PSI_MEANING, geometry, quantity
from . import grid as G
from . import schematic as SCH
from . import style as S

MAX_EVENTS = 40
PITCH = 0.32                   # inches per event row
GENE_ROW, NESTED_ROW, NESTED_ROW2 = 0.38, 0.24, 0.10     # the gene's row; its nested genes; a second label line
EXON_H, ARC_H, CELL_H = 0.11, 0.055, 0.21
BAR_W = 0.9                    # the narrowest colour bar
EXON_BAND = "#f5f4f0"          # the gene's exons, carried down through the event rows
INTRON_SHARE = 1.0             # introns are drawn shortened when longer than this many times the exons together
EXON_CAP, EXON_POWER = 3.0, 0.3   # an exonic segment longer than EXON_CAP x the median one grows as l^EXON_POWER
EXON_MIN_CAP = 300               # nt: shorter segments keep their length, however short the median one


class MapNotDrawn(UserWarning):
    """A gene map that could not be drawn; the probe goes on without it."""


def map_events(ds: Dataset, gene: str, ranked: pd.DataFrame, max_events: int = MAX_EVENTS,
               keep=None) -> tuple[list, dict]:
    """The events a gene's map draws, in genomic order 5' to 3': its best-ranked events that can be drawn (at most
    max_events, on the chromosome and strand of the best-ranked one; only those in `keep` when given), and how many
    of its events were left out: {"beyond": past max_events, "rarely observed": not in `keep` (probe.observed),
    "no coordinates": none given, "elsewhere": another chromosome or strand, "unreadable": coordinates that do not
    form the event}."""
    left = dict.fromkeys(("beyond", "rarely observed", "no coordinates", "elsewhere", "unreadable"), 0)
    span = {}
    for e in ranked.event_id[ranked.gene.eq(gene)]:
        if keep is not None and e not in keep:
            left["rarely observed"] += 1
            continue
        r = ds.events.loc[e]
        if not (r.variable and r.strand in ("+", "-") and r.chrom):
            left["no coordinates"] += 1
            continue
        try:
            span[e] = geometry(e, r).span
        except InputError:                               # the pages would refuse it too; the map leaves it out
            left["unreadable"] += 1
    ok = list(span)
    if not ok:
        return [], left
    chrom, strand = ds.events.at[ok[0], "chrom"], ds.events.at[ok[0], "strand"]
    same = [e for e in ok if ds.events.at[e, "chrom"] == chrom and ds.events.at[e, "strand"] == strand]
    left["elsewhere"] = len(ok) - len(same)
    keep = same[:max_events]
    left["beyond"] = len(same) - len(keep)
    key = (lambda e: (span[e][0], span[e][1])) if strand == "+" else (lambda e: (-span[e][1], -span[e][0]))
    return sorted(keep, key=key), left


def squeeze(exonic, lo: int, hi: int, share: float = INTRON_SHARE) -> tuple[np.ndarray, np.ndarray, str]:
    """A piecewise-linear map from genomic coordinates (lo to hi) to drawn ones. Exonic segments keep their length,
    except those longer than EXON_CAP times the median segment (and than EXON_MIN_CAP), drawn c (l / c)^EXON_POWER
    long (c the cap). The gaps between segments (introns) are shortened: a gap of l nt is drawn min(l, sqrt(l k))
    long, with k such that the gaps together take `share` times the exonic drawn length; nothing is shortened when
    they take less already. Returns the genomic knots, the drawn knots (np.interp maps between them) and what was
    shortened, in words ("" for nothing)."""
    segs: list[list[int]] = []
    for a, b in sorted((max(a, lo), min(b, hi)) for a, b in exonic if b > lo and a < hi):
        if segs and a <= segs[-1][1]:
            segs[-1][1] = max(segs[-1][1], b)
        else:
            segs.append([a, b])
    edges = sorted({lo, hi, *[x for sg in segs for x in sg]})
    parts = [(b - a, any(s0 <= a and b <= s1 for s0, s1 in segs)) for a, b in zip(edges[:-1], edges[1:])]
    lens = [n for n, ex in parts if ex]
    cap = max(EXON_CAP * float(np.median(lens)), EXON_MIN_CAP) if lens else 0.0
    exon_len = (lambda n: n if n <= cap else cap * (n / cap) ** EXON_POWER)      # noqa: E731
    long_exons = any(n > cap for n in lens) and len(lens) > 2
    if not long_exons:
        exon_len = (lambda n: n)                                                    # noqa: E731
    ex_len = sum(exon_len(n) for n in lens)
    gaps = [n for n, ex in parts if not ex]
    short_introns = bool(gaps) and ex_len > 0 and sum(gaps) > share * ex_len
    k = np.inf
    if short_introns:
        over = (lambda k_: sum(min(n, np.sqrt(n * k_)) for n in gaps) - share * ex_len)   # rises with k
        k0, k1 = 0.0, float(max(gaps))
        for _ in range(60):
            k0, k1 = ((k0 + k1) / 2, k1) if over((k0 + k1) / 2) < 0 else (k0, (k0 + k1) / 2)
        k = (k0 + k1) / 2
    drawn = [exon_len(n) if ex else min(n, np.sqrt(n * k)) for n, ex in parts]
    note = ("long introns and the longest exons shortened" if long_exons and short_introns else
            "the longest exons shortened" if long_exons else
            "long introns shortened, exons to scale" if short_introns else "")
    return np.array(edges, float), np.concatenate([[0.0], np.cumsum(drawn)]), note


def _medians(corr: pd.DataFrame | None, ids: list) -> dict:
    """(event, event) -> the median over the cohorts tested of the two events' rho."""
    if corr is None or not len(corr):
        return {}
    t = corr[corr.corr_status.eq("tested") & corr.kind.eq("event") & corr.event_id.isin(ids) & corr.partner.isin(ids)]
    pairs = {}
    for (a, b), v in t.groupby(["event_id", "partner"], sort=False).rho.median().items():
        pairs[(a, b)] = pairs[(b, a)] = v
    return pairs


def _not_drawn(left: dict, s: Settings) -> str:
    """' · 2 events not drawn (no coordinates, another chromosome)' or ''."""
    why = {"rarely observed": f"observed in under {s.observed_frac:.0%} of every cohort's survival samples",
           "no coordinates": "no coordinates", "elsewhere": "another chromosome or strand",
           "unreadable": "coordinates that do not form the event"}
    n = sum(left[k] for k in why)
    return (f" · {n} event{'s' if n != 1 else ''} not drawn ({', '.join(w for k, w in why.items() if left[k])})"
            if n else "")


def gene_map(ds: Dataset, gene: str, ranked: pd.DataFrame, cells: pd.DataFrame, endpoint: str,
             settings: Settings | None = None, *, gtf=None, gex_cells: pd.DataFrame | None = None, gex_model: str = "",
             corr: pd.DataFrame | None = None, fit_note: str = "", max_events: int = MAX_EVENTS, keep=None):
    """The gene map of `gene` (a matplotlib Figure), or None when none of its events can be drawn.

    `ranked`: the probe's ranked events (probe.rank_events); `cells`: its event x cohort table; `gtf`: a GTF path or
    the records read from one (annotation.read_gtf), for the gene model; `gex_cells`: the host gene's statistics
    (expression_cells.csv) for the gene's row, with `gex_model` the name of its Cox model; `corr`: the correlations
    (correlation.correlations) for the triangle of the events' pairwise correlation; `fit_note`: which fit the cells
    show, as the overview's caption says it (default: the adjusted model where fitted); `keep`: the events present
    enough to draw (probe.observed; default all)."""
    import matplotlib
    from matplotlib.figure import Figure
    from matplotlib.patches import Rectangle

    from .panel import _model

    s = settings or Settings()
    ids, left = map_events(ds, gene, ranked, max_events, keep)
    if not ids:
        return None
    E = ds.events.loc[ids]
    chrom, strand = E.chrom.iloc[0], E.strand.iloc[0]
    geoms = [geometry(e, E.loc[e]) for e in ids]
    model = _model(gtf, E, geoms, s) if gtf is not None else None
    host = E.expression_gene.mode().iloc[0]                 # the gene's row: the expression its events adjust for
    rank = ranked.set_index("event_id")["rank"].to_dict()
    cohorts = sorted(cells.cohort.unique())
    nc, n = len(cohorts), len(ids)
    idx = cells.set_index(["event_id", "cohort"])
    gx = None
    if gex_cells is not None and len(gex_cells):
        g = gex_cells[gex_cells.gene.eq(host)]
        g = g[g.endpoint.eq(endpoint)] if "endpoint" in g else g
        gx = g.set_index("cohort") if len(g) else None
    rho_e = _medians(corr, ids)
    items = ids if rho_e and len(ids) > 1 else []             # the triangle's rows and columns
    U = s.psi_hr_unit.upper()
    v = " or ".join(dict.fromkeys(sorted(E.event_type.map(quantity), key=lambda q: q != "PSI")))
    # ------------------------------------------------------------------ the drawn span
    pieces = [p for g in geoms for p in g.constant + g.variable + g.other]
    e_lo, e_hi = min(p[0] for p in pieces), max(p[1] for p in pieces)
    blocks = list(model.blocks) if model is not None else []
    gene_span = (blocks[0][0], blocks[-1][1]) if blocks else None
    win = SCH.window_for((e_lo, e_hi), model)
    if win is not None and blocks:
        blocks = [(max(a, win[0]), min(b, win[1]), k) for a, b, k in blocks if b > win[0] and a < win[1]]
    lo = min([b[0] for b in blocks] + [e_lo] + ([win[0]] if win else []))
    hi = max([b[1] for b in blocks] + [e_hi] + ([win[1]] if win else []))
    nested = pd.DataFrame(columns=["gene_name", "start", "end"])
    if model is not None:
        nested = model.nested[(model.nested.end > lo) & (model.nested.start < hi)]
    exonic = [(a, b) for a, b, _ in blocks] + [p for g in geoms for p in g.constant + g.other] + [
        p for g in geoms if g.event_type != "RI" for p in g.variable]     # a retained intron is an intron
    knots, drawn, scale_note = squeeze(exonic, lo, hi)
    # ------------------------------------------------------------------ horizontal layout (inches)
    labels = {e: str(E.at[e, "label"]) for e in ids}
    rank_w = S.text_width(str(max(rank.get(e, 0) for e in ids)), 6.0) + 0.04
    lab_w = max([S.text_width(labels[e], 6.6, weight="bold") for e in ids]
                + [S.text_width(PSI_MEANING.get(g.event_type, "PSI"), 5.6) for g in geoms]
                + [S.text_width(gene, 7.2, weight="bold", style="italic")]
                + ([S.text_width("snoRNAs", 6.0)] if len(nested) else []))
    x_rank = 0.12 + rank_w                                   # rank numbers end here; labels start after
    x_lab = x_rank + 0.08
    x_s0 = x_lab + lab_w + 0.22                              # the structure
    cw = min(0.22, max(0.09, 1.76 / max(nc, 1)))             # one survival cell
    strip_head = f"Cox HR per {U} · {endpoint}"
    strip_w = max(nc * cw, S.text_width(strip_head, 6.2, weight="bold"), BAR_W)
    mw = min(0.26, max(0.12, 2.0 / max(len(items), 1))) if items else 0.0
    mat_w = max(len(items) * mw, S.text_width("Spearman ρ, median of cohorts", 5.6), BAR_W) if items else 0.0
    fixed = x_s0 + 0.34 + strip_w + (0.34 + mat_w if items else 0.0) + 0.14
    W = max(7.2, fixed + 3.2)
    struct_w = W - fixed
    x_strip = x_s0 + struct_w + 0.34
    x_mat = x_strip + strip_w + 0.34
    span_d = drawn[-1]
    pad = 0.015 * span_d
    x0, x1 = -pad, span_d + pad
    per_in = (x1 - x0) / struct_w                            # drawn units per inch
    X = (lambda c: float(np.interp(c, knots, drawn)))                   # noqa: E731
    page_x = (lambda d: (d - x0) / per_in if strand != "-" else (x1 - d) / per_in)   # noqa: E731  (inches)
    # nested genes' names: on one line where they fit, else on a second, else left out
    name_row, ends = {}, [-np.inf, -np.inf]
    for r in nested.sort_values("start").itertuples() if strand != "-" else \
            nested.sort_values("end", ascending=False).itertuples():
        c_in = page_x(X((r.start + r.end) / 2))
        half = S.text_width(r.gene_name, 5.4, style="italic") / 2
        row = next((k for k in (0, 1) if c_in - half > ends[k] + 0.04), None)
        if row is not None:
            name_row[r.gene_name, r.start] = row
            ends[row] = c_in + half
    two_lines = any(row == 1 for row in name_row.values())
    # ------------------------------------------------------------------ vertical layout (inches from the top)
    head_y = 0.56                                            # the blocks' headers
    lab_h = max(S.text_width(c, 5.6) for c in cohorts) + 0.04 if cohorts else 0.2
    y0 = head_y + 0.14 + lab_h                               # the gene's row starts
    gene_row = model is not None or gx is not None              # else the title names the gene
    gene_h = (GENE_ROW if gene_row else 0.0) + (NESTED_ROW + (NESTED_ROW2 if two_lines else 0.0) if len(nested)
                                                else 0.0)
    yg = y0 + GENE_ROW / 2
    y_ev = y0 + gene_h
    yc = {e: y_ev + PITCH * (i + 0.5) for i, e in enumerate(ids)}
    y_end = y_ev + PITCH * n
    caption = (f"Rows: the gene and its probed events, 5′ to 3′, numbered by probe rank"
               + (f" (the best-ranked {n} of {n + left['beyond']})" if left["beyond"] else "")
               + ("; pale columns: the gene's exons (the GTF's transcripts collapsed)" if model is not None else
                  "; no gene model (no GTF, or the gene is not in it)")
               + (f"; {scale_note}" if scale_note else "")
               + f"; cells: Cox HR per {U} of {v} per cohort ({fit_note or 'adjusted model where fitted'})"
               + (f"; the gene's row: HR per SD of {'its' if host == gene else host} expression (Cox on "
                  f"{gex_model or 'expression'}), no q" if gx is not None else "")
               + f"; dot: p < {s.alpha:g} (large: < 0.01)"
               + (f"; {S.Q_MARK}: q < {s.q_mark_below:g}" if s.q_mark_below > 0 else "")
               + "; frame: group hit; grey: not tested"
               + (f"; ρ between events: Spearman correlation of two events' {v} in each cohort's survival samples, "
                  "median over the cohorts (each cohort in correlation.png)"
                  + (f"; frame: |ρ| ≥ {s.corr_note_above:g}" if s.corr_note_above > 0 else "") if items else ""))
    cap = G.clauses(caption, W - 0.24, 5.6)
    H = y_end + 0.62 + 0.115 * len(cap) + 0.08
    hr_scale, rho_scale = G.hr_scale(), G.rho_scale()
    with matplotlib.rc_context(S.rc()):
        fig = Figure(figsize=(W, H))
        fx = (lambda x: x / W)                                            # noqa: E731
        fy = (lambda y: 1 - y / H)                                        # noqa: E731

        def axes(x, w, top):
            ax = fig.add_axes([fx(x), fy(y_end), fx(w), (y_end - top) / H])
            ax.set_ylim(y_end, top)
            ax.axis("off")
            return ax
        # title and subtitle
        fig.text(fx(0.12), fy(0.10), gene, fontsize=8.5, fontweight="bold", style="italic", va="top")
        tw = 0.12 + S.text_width(gene, 8.5, weight="bold", style="italic") + 0.05
        fig.text(fx(tw), fy(0.10), f"gene map · {n} event{'s' if n != 1 else ''} · {endpoint}", fontsize=8.5,
                 fontweight="bold", va="top")
        a0, a1 = (lo, hi) if win else (gene_span or (lo, hi))
        sub = (f"{chrom}:{a0 + 1:,}–{a1:,} · " + ("minus strand, drawn 5′→3′" if strand == "-" else "plus strand")
               + (f" · view {SCH._kb(a1 - a0)} of a {SCH._kb(gene_span[1] - gene_span[0])} gene"
                  if win and gene_span else "") + _not_drawn(left, s))
        fig.text(fx(0.12), fy(0.30), sub, fontsize=6.0, color=S.INK2, va="top")
        # ---------------------------------------------------------------- structure
        ax = axes(x_s0, struct_w, head_y - 0.12)
        ax.set_xlim((x1, x0) if strand == "-" else (x0, x1))
        x_end = x1 - pad if strand != "-" else x0 + pad      # the right end, 3' (drawn)
        sgn = -1 if strand != "-" else 1                     # leftwards on the page
        if scale_note:
            ax.text(x_end, head_y, scale_note, ha="right", va="center", fontsize=5.8, color=S.INK2)
        else:
            bar, bar_text = SCH.scale_bar(hi - lo)
            ax.plot([x_end, x_end + sgn * bar], [head_y, head_y], color=S.INK2, lw=0.8, solid_capstyle="butt")
            ax.text(x_end + sgn * bar / 2, head_y - 0.04, bar_text, ha="center", va="bottom", fontsize=5.8,
                    color=S.INK2)

        def bx(a, b, y, h, min_in=0.016, **kw):              # a box at least min_in wide
            xa, xb = X(a), X(b)
            if xb - xa < min_in * per_in:
                c = (xa + xb) / 2
                xa, xb = c - min_in * per_in / 2, c + min_in * per_in / 2
            SCH.box(ax, xa, xb, y, h, **kw)
            return xb - xa

        def arc(a, b, y, up, **kw):                          # flatter over short spans
            h = min(ARC_H, 0.35 * abs(X(b) - X(a)) / per_in)
            SCH.arc(ax, X(a), X(b), y, h if up else -h, **kw)
        for a, b, _ in blocks:                               # the gene's exons, down through the event rows
            ax.add_patch(Rectangle((X(a), yg), X(b) - X(a), y_end - yg, facecolor=EXON_BAND, lw=0, zorder=0))
        if model is not None:
            g0 = lo if (win and gene_span[0] < lo) or not blocks else blocks[0][0]
            g1 = hi if (win and gene_span[1] > hi) or not blocks else blocks[-1][1]
            ax.plot([X(g0), X(g1)], [yg, yg], color=S.EXON_LINE, lw=0.6, zorder=1)
            for a, b, _ in blocks:
                bx(a, b, yg, 0.14, facecolor=S.EXON, zorder=2)
            if win:
                for x_, beyond_ in ((lo, gene_span[0] < lo), (hi, gene_span[1] > hi)):
                    if beyond_:
                        left_edge = (x_ == lo) == (strand != "-")
                        ax.text(X(x_), yg, "‹" if left_edge else "›", ha="center", va="center", fontsize=9,
                                color=S.EXON_LINE)
            ys = yg + 0.20
            for r in nested.itertuples():
                bx(r.start, r.end, ys, 0.075, min_in=0.03, facecolor=S.NESTED, zorder=2)
                row = name_row.get((r.gene_name, r.start))
                if row is not None:
                    ax.text(X((r.start + r.end) / 2), ys + 0.055 + 0.09 * row, r.gene_name, ha="center", va="top",
                            fontsize=5.4, style="italic", color=S.INK)
        elif gene_row:
            ax.plot([X(lo), X(hi)], [yg, yg], color=S.GRID, lw=0.8, ls=(0, (2, 2)), zorder=1)
        for e, g in zip(ids, geoms):
            y = yc[e]
            col = S.TYPE_COLOR.get(g.event_type, S.TYPE_DEFAULT)
            pp = g.constant + g.variable + g.other
            ax.plot([X(min(p[0] for p in pp)), X(max(p[1] for p in pp))], [y, y], color=S.EXON_LINE, lw=0.5,
                    zorder=1)
            for a, b in g.constant:
                bx(a, b, y, EXON_H, facecolor=S.EXON, zorder=2)
            for a, b in g.variable:
                w = bx(a, b, y, 0.06 if g.event_type == "RI" else EXON_H, facecolor=col, zorder=3)
                if w < 0.03 * per_in:                        # too narrow to see: a caret above it
                    ax.plot([X((a + b) / 2)], [y - 0.095], marker="v", ms=2.6, color=col, mew=0, zorder=4)
            for a, b in g.other:
                bx(a, b, y, EXON_H - 0.02, facecolor="white", edgecolor=col, linewidth=0.7, zorder=3)
            for a, b in g.psi_arcs:
                arc(a, b, y - EXON_H / 2, True, edgecolor=col, lw=0.7, zorder=2)
            for a, b in g.other_arcs:
                arc(a, b, y + EXON_H / 2, False, edgecolor=S.EXON_LINE, lw=0.7, zorder=2)
        # ---------------------------------------------------------------- row labels
        if gene_row:
            fig.text(fx(x_lab), fy(yg), gene, fontsize=7.2, fontweight="bold", style="italic", va="center")
        if len(nested):
            fig.text(fx(x_lab), fy(yg + 0.20), "snoRNAs" if set(s.nested_biotypes) <= {"snoRNA", "scaRNA"}
                     else "Nested genes", fontsize=6.0, color=S.INK2, va="center")
        for e, g in zip(ids, geoms):
            col = S.TYPE_COLOR.get(g.event_type, S.TYPE_DEFAULT)
            fig.text(fx(x_rank), fy(yc[e] - 0.045), str(rank.get(e, "")), fontsize=6.0, color=S.MUTED, ha="right",
                     va="center")
            fig.text(fx(x_lab), fy(yc[e] - 0.045), labels[e], fontsize=6.6, fontweight="bold", color=col,
                     va="center")
            fig.text(fx(x_lab), fy(yc[e] + 0.075), PSI_MEANING.get(g.event_type, "PSI"), fontsize=5.6,
                     color=S.INK2, va="center")
        # ---------------------------------------------------------------- survival cells
        cmap, norm = hr_scale
        ax = axes(x_strip, nc * cw, y0)
        ax.set_xlim(0, nc)
        fig.text(fx(x_strip), fy(head_y), strip_head, fontsize=6.2, fontweight="bold", va="center")
        for j, c in enumerate(cohorts):
            ax.text(j + 0.5, y0 - 0.04, c, rotation=90, ha="center", va="bottom", fontsize=5.6, color=S.INK2,
                    clip_on=False)

        def cell(j, y, tested, hr, p, q, hit):
            face = cmap(norm(np.log2(hr))) if tested and np.isfinite(hr) else G.UNTESTED
            ax.add_patch(Rectangle((j + 0.05, y - CELL_H / 2), 0.9, CELL_H, facecolor=face, lw=0))
            if hit:
                ax.add_patch(Rectangle((j + 0.05, y - CELL_H / 2), 0.9, CELL_H, fill=False, edgecolor=S.INK,
                                       lw=0.7))
            if tested and s.q_mark_below > 0 and np.isfinite(q) and q < s.q_mark_below:
                ax.text(j + 0.5, y + 0.035, S.Q_MARK, fontsize=8.5, color=G.ink_on(face), ha="center", va="center",
                        fontweight="bold")
            elif tested and np.isfinite(p) and p < s.alpha:
                ax.plot([j + 0.5], [y], "o", ms=2.6 if p >= 0.01 else 3.6, color=G.ink_on(face), mew=0)
        for j, c in enumerate(cohorts):
            if gx is not None and c in gx.index:
                x = gx.loc[c]
                cell(j, yg, x.get("cox_status") == "tested", x.get("hr_per_sd", np.nan), x.get("cox_p", np.nan),
                     np.nan, G.flag(x.get("group_hit")))
            for e in ids:
                if (e, c) in idx.index:
                    x = idx.loc[(e, c)]
                    cell(j, yc[e], *G.shown_fit(x, s), G.flag(x.get("group_hit")))
        G.colorbar(fig, x_strip, y_end + 0.20, min(1.4, max(nc * cw, 0.9)), W, H, hr_scale, [-1, 0, 1],
                   ["0.5", "1", "2"], f"HR per {U}" + (" (expression: per SD)" if gx is not None and U != "SD"
                                                       else ""))
        # ---------------------------------------------------------------- correlation triangle
        if items:
            cmap_r, norm_r = rho_scale
            ax = axes(x_mat, len(items) * mw, y0)
            ax.set_xlim(0, len(items))
            fig.text(fx(x_mat), fy(head_y), "ρ between events", fontsize=6.2, fontweight="bold", va="center")
            fig.text(fx(x_mat), fy(head_y + 0.13), "Spearman ρ, median of cohorts", fontsize=5.6, color=S.INK2,
                     va="center")
            show = mw >= 0.19                                 # room to print the value
            for i, a in enumerate(items):
                y = yc[a]
                ax.text(i + 0.5, y, str(rank.get(a, "")), ha="center", va="center", fontsize=5.6, color=S.MUTED)
                for j, b in enumerate(items[:i]):
                    r = rho_e.get((a, b))
                    ok = r is not None and np.isfinite(r)
                    face = cmap_r(norm_r(r)) if ok else G.UNTESTED
                    ax.add_patch(Rectangle((j + 0.05, y - CELL_H / 2), 0.9, CELL_H, facecolor=face, lw=0))
                    if ok and s.corr_note_above > 0 and abs(r) >= s.corr_note_above:
                        ax.add_patch(Rectangle((j + 0.05, y - CELL_H / 2), 0.9, CELL_H, fill=False,
                                               edgecolor=S.INK, lw=0.7))
                    if ok and show:
                        ax.text(j + 0.5, y, G.rho_text(r), ha="center", va="center", fontsize=4.9,
                                color=G.ink_on(face))
            G.colorbar(fig, x_mat, y_end + 0.20, min(1.4, max(len(items) * mw, 0.9)), W, H, rho_scale,
                       [-1, 0, 1], ["−1", "0", "1"], "Spearman ρ")
        for k, line in enumerate(cap):
            fig.text(fx(0.12), fy(y_end + 0.62 + 0.115 * k), line, fontsize=5.6, color=S.INK2, va="top")
    return fig
