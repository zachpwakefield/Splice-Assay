"""The protein band of an event panel: the suggested protein of each form, its features, and the event's residues.

One row per protein-coding isoform (the PSI form first), on a common residue axis aligned with the schematic. The
backbone is a thin rounded bar. The features drawn are:
  - named InterPro entries: rounded boxes, named inside when the name fits (else in the legend); of overlapping
    entries, one is shown;
  - disordered regions (MobiDB-lite): a line under the backbone;
  - transmembrane helices and signal peptides: small dark boxes;
  - ELM motifs: lollipops above the backbone.
The event's own residues are drawn as in the schematic: the backbone in the event colour, a tint over the domains it
touches and a bracket with the residue numbers. On an isoform that lacks the region, a triangle marks where it would
sit; an event that encodes no residue gets a triangle before residue 1 (in the 5' UTR) or after the last residue (in
the 3' UTR, or holding only the stop codon). Above the rows, one sentence says what the event does. The band is drawn
only when there is a suggestion (protein.ProteinChange.ok).
"""
from __future__ import annotations

from ..protein import domain_clusters, form_names
from . import style as S

DOMAIN_COLORS = ("#9dbde0", "#ebc07f", "#a9d3a2", "#e0a9b2", "#c1b3df", "#d3bf97", "#92cec6", "#e9b08b",
                 "#bfcd86", "#d8aad3")
BACKBONE, DISORDER, TM, SIGNAL, MOTIF = "#cfcdc4", "#a9a79e", "#4a4a4a", "#c98a2b", "#6d6b65"
ROW_H, AXIS_H, LINE_H, HEAD_PAD = 0.42, 0.26, 0.13, 0.06
DOM_H, BONE_W, LABEL_PT = 0.15, 3.0, 5.4       # domain box height (in), backbone width (pt), label size (pt)
X0 = 1.05                                      # left edge of the rows (inches), as the schematic
SHORT = (("Immunoglobulin", "Ig"),)


def _wrap(text: str, width: float, size: float) -> list[str]:
    lines, cur = [], ""
    for w in text.split(" "):
        nxt = f"{cur} {w}".strip()
        if cur and S.text_width(nxt, size) > width:
            lines.append(cur)
            cur = w
        else:
            cur = nxt
    return lines + ([cur] if cur else [])


def _ticks(n: int) -> list[int]:
    step = next(s for s in (50, 100, 250, 500, 1000, 2500, 5000) if n / s <= 8)
    return list(range(0, n + 1, step))


def _in_box(label: str, width_in: float) -> str | None:
    """The domain's name to write inside a box this wide (inches): the name, a shorter form, or None."""
    short = label
    for a, b in SHORT:
        short = short.replace(a, b)
    short = short.removesuffix(" domain")
    for text in dict.fromkeys((label, short)):
        if S.text_width(text, LABEL_PT) + 0.08 <= width_in:
            return text
    return None


def prepare(pc, W: float, color: str, label: str = "") -> dict | None:
    """Everything the band needs, and its height (inches); None when there is no suggestion to show."""
    if pc is None or not pc.ok:
        return None
    isos = list(pc.shown)
    if not isos:
        return None
    lines = _wrap(pc.effect["text"][:1].upper() + pc.effect["text"][1:] + ".", W - 0.27, 6.2)
    n_max = max(len(i.protein) for i in isos)
    inch_per_res = (W - X0 - 0.15) / (1.02 * n_max)
    doms = {i.form: [(lab, a, b, _in_box(lab, (b - a + 1) * inch_per_res)) for lab, a, b in
                     domain_clusters(i.features)] for i in isos}
    names = list(dict.fromkeys(d[0] for i in isos for d in doms[i.form]))
    colors = {n: DOMAIN_COLORS[k % len(DOMAIN_COLORS)] for k, n in enumerate(names)}
    unnamed = [n for n in names if not any(d[3] for i in isos for d in doms[i.form] if d[0] == n)]
    kinds = set().union(*[set(i.features.database[i.features.named]) for i in isos])
    legend = [("dom", n) for n in unnamed]
    legend += [(k, t) for k, t in (("mobidblite", "Disordered"), ("tmhmm", "Transmembrane helix"),
                                   ("signalp", "Signal peptide"), ("elm", "ELM motif")) if k in kinds]
    if any(i.event_aa for i in isos):
        legend.append(("event", "Event residues"))
    where = next((_outside(i) for i in isos if _outside(i)), "")
    if where:
        legend.append(("utr", OUTSIDE[where]))
    if any(i.insert_after for i in isos):
        legend.append(("insert", "Where the other form's region sits"))
    leg_lines, cur, x = [], [], 0.12
    for k, t in legend:
        w = 0.21 + S.text_width(t, 5.8) + 0.18
        if cur and x + w > W - 0.05:
            leg_lines.append(cur)
            cur, x = [], 0.12
        cur.append((k, t))
        x += w
    leg_lines.append(cur)
    h = HEAD_PAD + LINE_H * (1 + len(lines)) + 0.04 + ROW_H * len(isos) + AXIS_H + LINE_H * len(leg_lines) + 0.04
    return dict(pc=pc, isos=isos, label=label, lines=lines, doms=doms, colors=colors, legend=leg_lines,
                n_max=n_max, color=color, height=h)


OUTSIDE = {"5′ UTR": "Event: in the 5′ UTR, before residue 1", "3′ UTR": "Event: in the 3′ UTR, after the last residue",
           "stop codon": "Event: holds the stop codon, after the last residue"}


def _outside(iso) -> str:
    """Where the event sits when it encodes none of this isoform's residues (and has no insertion point): '5′ UTR',
    '3′ UTR' or 'stop codon'; '' otherwise."""
    if iso.event_aa or iso.insert_after:
        return ""
    return iso.utr or ("stop codon" if iso.codon == "stop codon" else "")


def _rbox(ax, x0, x1, yc, h, xpi, ypi, r_in=0.025, **kw):
    """A box with round corners (radius r_in inches) on axes scaled xpi, ypi data units per inch."""
    from matplotlib.patches import FancyBboxPatch
    w = x1 - x0
    r = max(min(r_in * xpi, 0.5 * w, 0.5 * h * xpi / ypi), 1e-9)
    ax.add_patch(FancyBboxPatch((x0, yc - h / 2), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
                                mutation_aspect=ypi / xpi, **{"linewidth": 0, **kw}))


def draw(fig, W: float, H: float, y0: float, band: dict, rows: list) -> None:
    """Draw a prepared band with its top at y0 (inches from the top of the figure); append its values to rows."""
    from matplotlib.lines import Line2D
    from matplotlib.patches import Rectangle

    fx = lambda x: x / W                                     # noqa: E731
    fy = lambda y: 1 - y / H                                 # noqa: E731
    pc, isos, color = band["pc"], band["isos"], band["color"]
    y = y0 + HEAD_PAD
    head = (f"{band['label']} · " if band["label"] else "") + "Protein"
    fig.text(fx(0.12), fy(y), head, fontsize=6.6, fontweight="bold", color=S.INK, ha="left", va="top")
    fig.text(fx(0.12 + S.text_width(head, 6.6, weight="bold") + 0.04), fy(y),
             "· suggested from annotation, not measured", fontsize=6.2, color=S.INK2, ha="left", va="top")
    for i, line in enumerate(band["lines"]):
        fig.text(fx(0.12), fy(y + (i + 1) * LINE_H), line, fontsize=6.2, color=S.INK, ha="left", va="top")
    y += LINE_H * (1 + len(band["lines"])) + 0.04
    rows_h = ROW_H * len(isos)
    ax_w = W - X0 - 0.15
    ax = fig.add_axes([fx(X0), fy(y + rows_h), fx(ax_w), rows_h / H])
    n = band["n_max"]
    ax.set_xlim(-0.01 * n, n * 1.01)
    ax.set_ylim(len(isos) - 0.5, -0.5)
    ax.set_yticks([])
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)
    ax.set_xticks(_ticks(n))
    ax.set_xlabel("Residue", fontsize=6, labelpad=1.5, color=S.INK2)
    ax.tick_params(axis="x", labelsize=5.8, length=2, color=S.MUTED)
    ax.spines["bottom"].set_color(S.MUTED)
    ypi = len(isos) / rows_h                                 # data units per inch: rows, residues
    xpi = 1.02 * n / ax_w
    rows.append(dict(panel="protein", what="effect", event_id=pc.event_id, value=pc.effect["text"],
                     source="annotation"))
    for k, iso in enumerate(isos):
        L = len(iso.protein)
        yc = k
        form = form_names(pc.event_type)[0 if iso.form == "INC" else 1]
        fig.text(fx(0.12), fy(y + k * ROW_H + ROW_H / 2 - 0.05), form[:1].upper() + form[1:], fontsize=6.2,
                 fontweight="bold", ha="left", va="center")
        fig.text(fx(0.12), fy(y + k * ROW_H + ROW_H / 2 + 0.07), f"{iso.name} · {L} aa", fontsize=5.6,
                 color=S.INK2, ha="left", va="center")
        ax.add_line(Line2D([1, L], [yc, yc], lw=BONE_W, color=BACKBONE, solid_capstyle="round", zorder=1))
        f = iso.features[iso.features.named]
        for r in f[f.database.eq("mobidblite")].itertuples():
            ax.add_line(Line2D([r.start, r.stop], [yc + 0.10 * ypi] * 2, lw=1.4, color=DISORDER,
                               solid_capstyle="round", zorder=1))
        for a, b in iso.event_aa:                            # the event: halo, coloured backbone, bracket
            _rbox(ax, a - 0.5, b + 0.5, yc, (DOM_H + 0.08) * ypi, xpi, ypi, r_in=0.03, facecolor=color,
                  alpha=0.2, zorder=2.8)                    # a tint over the domains it touches
            ax.add_line(Line2D([a, b], [yc, yc], lw=BONE_W, color=color, solid_capstyle="butt", zorder=1.5))
            yb = yc - (DOM_H / 2 + 0.07) * ypi
            ax.add_line(Line2D([a - 0.5, a - 0.5, b + 0.5, b + 0.5], [yb + 0.025 * ypi, yb, yb, yb + 0.025 * ypi],
                               lw=0.7, color=color, zorder=6))
            ax.text((a + b) / 2, yb - 0.015 * ypi, f"{a}–{b}" if b > a else f"{a}", fontsize=5.4, color=color,
                    ha="center", va="bottom", zorder=6)
            rows.append(dict(panel="protein", what="event_residues", event_id=pc.event_id, isoform=iso.transcript,
                             start=a, stop=b))
        for lab, a, b, text in band["doms"][iso.form]:
            _rbox(ax, a - 0.5, b + 0.5, yc, DOM_H * ypi, xpi, ypi, facecolor=band["colors"][lab], zorder=2)
            if text:
                ax.text((a + b) / 2, yc, text, fontsize=LABEL_PT, color=S.INK, ha="center", va="center", zorder=3,
                        clip_on=True)
            rows.append(dict(panel="protein", what="domain", event_id=pc.event_id, isoform=iso.transcript,
                             label=lab, start=a, stop=b))
        for db, col in (("tmhmm", TM), ("signalp", SIGNAL)):
            for r in f[f.database.eq(db)].itertuples():
                _rbox(ax, r.start - 0.5, r.stop + 0.5, yc, 0.11 * ypi, xpi, ypi, r_in=0.012, facecolor=col,
                      zorder=2.5)
        for r in f[f.database.eq("elm")].itertuples():
            xm = (r.start + r.stop) / 2
            top = yc - (DOM_H / 2 + 0.045) * ypi
            ax.add_line(Line2D([xm, xm], [yc - 0.02 * ypi, top], color=MOTIF, lw=0.5, zorder=1.2))
            ax.add_line(Line2D([xm], [top], marker="o", ms=2.2, ls="", color=MOTIF, mew=0, zorder=1.2))
        if iso.insert_after:
            xi = iso.insert_after + 0.5
            ax.add_line(Line2D([xi, xi], [yc - 0.05 * ypi, yc + 0.05 * ypi], color=color, lw=0.9, zorder=4))
            ax.add_line(Line2D([xi], [yc - (DOM_H / 2 + 0.05) * ypi], marker="v", ms=3.4, ls="", color=color,
                               zorder=6))
            rows.append(dict(panel="protein", what="insert_after", event_id=pc.event_id, isoform=iso.transcript,
                             value=iso.insert_after))
        where = _outside(iso)
        if where:                                            # before residue 1 (5' UTR) or after the last residue
            xu = 0.5 if where == "5′ UTR" else L + 0.5
            ax.add_line(Line2D([xu], [yc - (DOM_H / 2 + 0.05) * ypi], marker="v", ms=3.4, ls="", color=color,
                               zorder=6))
            rows.append(dict(panel="protein", what="event_outside", event_id=pc.event_id, isoform=iso.transcript,
                             value=where, position=xu))
        rows.append(dict(panel="protein", what="isoform", event_id=pc.event_id, isoform=iso.transcript,
                         label=iso.name, form=iso.form, value=L, biotype=iso.biotype, tsl=iso.tsl,
                         match_context=iso.context, source="annotation"))
    # legend
    y_leg = y + rows_h + AXIS_H
    for li, items in enumerate(band["legend"]):
        x = 0.12
        yy = fy(y_leg + 0.06 + li * LINE_H)
        for kind, text in items:
            if kind == "dom":
                fig.add_artist(Rectangle((fx(x), yy - 0.045 / H), fx(0.16), 0.09 / H, lw=0,
                                         facecolor=band["colors"][text]))
            elif kind == "mobidblite":
                fig.add_artist(Line2D([fx(x + 0.01), fx(x + 0.15)], [yy, yy], lw=1.4, color=DISORDER,
                                      solid_capstyle="round"))
            elif kind in ("tmhmm", "signalp"):
                fig.add_artist(Rectangle((fx(x + 0.04), yy - 0.04 / H), fx(0.08), 0.08 / H, lw=0,
                                         facecolor=TM if kind == "tmhmm" else SIGNAL))
            elif kind == "elm":
                fig.add_artist(Line2D([fx(x + 0.08)] * 2, [yy - 0.035 / H, yy + 0.025 / H], color=MOTIF, lw=0.5))
                fig.add_artist(Line2D([fx(x + 0.08)], [yy + 0.025 / H], marker="o", ms=2.2, ls="", color=MOTIF,
                                      mew=0))
            elif kind == "event":
                fig.add_artist(Rectangle((fx(x), yy - 0.05 / H), fx(0.16), 0.10 / H, facecolor=color, alpha=0.2,
                                         lw=0))
                fig.add_artist(Line2D([fx(x), fx(x + 0.16)], [yy, yy], lw=BONE_W, color=color,
                                      solid_capstyle="butt"))
            elif kind in ("insert", "utr"):
                fig.add_artist(Line2D([fx(x + 0.08)], [yy], marker="v", ms=3.4, ls="", color=color))
            fig.text(fx(x + 0.21), yy, text, fontsize=5.8, color=S.INK2, ha="left", va="center")
            x += 0.21 + S.text_width(text, 5.8) + 0.18
