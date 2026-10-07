"""Pieces shared by the probe's summary figures (the overview, the gene maps and the correlation grid): their colour
scales, the fit a survival cell shows and how it is drawn, colour bars, keys to the marks, and captions wrapped between
their clauses."""
from __future__ import annotations

from typing import NamedTuple

import numpy as np

from . import style as S

HR_COLORS = ["#2f5f98", "#9ebbd9", "#f4f3ee", "#e6a88a", "#b0412c"]    # log2 HR, -HR_LIM to HR_LIM
RHO_COLORS = ["#5e3c99", "#b2abd2", "#f4f3ee", "#a6dba0", "#1b7837"]   # Spearman rho, -1 (purple) to 1 (green)
HR_LIM = 1.5
UNTESTED = "#e8e7e1"
NEUTRAL = HR_COLORS[2]
PALE = 0.82                         # an imprecise cell's colour: this share of the way to the neutral middle
HATCH = "#b3b1a9"                   # and its hatching
KEY_W, KEY_H, KEY_ROW = 0.17, 0.13, 0.17   # inches: a key's sample cell, and one row of the key


def hr_scale():
    """(colormap, norm) of the HR cells: log2 HR from -1.5 to 1.5."""
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    return LinearSegmentedColormap.from_list("hr", HR_COLORS), Normalize(-HR_LIM, HR_LIM)


def rho_scale():
    """(colormap, norm) of the correlation cells: rho from -1 to 1."""
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    return LinearSegmentedColormap.from_list("rho", RHO_COLORS), Normalize(-1, 1)


def flag(v) -> bool:
    """A True/False cell that may be missing (NaN, None), a numpy bool, or text read back from a table."""
    if isinstance(v, str):
        return v.strip().lower() in ("true", "1", "1.0")
    return False if v is None or (isinstance(v, float) and np.isnan(v)) else bool(v)


def _luminance(rgb) -> float:
    """WCAG relative luminance of an sRGB colour (components 0-1)."""
    c = np.asarray(rgb[:3], float)
    c = np.where(c <= 0.03928, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)
    return float(0.2126 * c[0] + 0.7152 * c[1] + 0.0722 * c[2])


def ink_on(face) -> str:
    """The mark colour (white or ink) with the higher contrast on a cell of colour `face`."""
    from matplotlib.colors import to_rgb
    lum = _luminance(to_rgb(face))
    ink = _luminance(to_rgb(S.INK))
    return "white" if (1.05 / (lum + 0.05)) > ((lum + 0.05) / (ink + 0.05)) else S.INK


class Shown(NamedTuple):
    """The fit a survival cell shows."""
    tested: bool
    hr: float
    p: float
    q: float
    imprecise: bool


def imprecise(lo, hi, low_power, s) -> bool:
    """Whether a fit's HR is too uncertain to be read from its colour: its 95% CI (lo, hi) spans more than
    Settings.imprecise_ci_ratio, or it is low power. Never when that setting is 0."""
    if not s.imprecise_ci_ratio:
        return False
    if flag(low_power):
        return True
    try:
        lo, hi = float(lo), float(hi)
    except (TypeError, ValueError):
        return False
    return not (np.isnan(lo) or np.isnan(hi)) and hi > lo * s.imprecise_ci_ratio


def shown_fit(x, s) -> Shown:
    """The fit a survival cell shows: the adjusted model where it was fitted, else the base model; the HR per
    Settings.psi_hr_unit. `x`: a row of the probe's cells."""
    adj = "adj_cox_status" in x and x.adj_cox_status == "tested"
    u, pre = s.psi_hr_unit, "adj_" if adj else ""
    hr, p = (x[f"adj_hr_per_{u}"], x.adj_cox_p) if adj else (x.get(f"hr_per_{u}"), x.get("cox_p"))
    q = x.get("adj_cox_q", np.nan) if adj else x.get("cox_q", np.nan)
    tested = bool(adj or x.cox_status == "tested")
    return Shown(tested, hr, p, q, tested and imprecise(x.get(f"{pre}ci_low_{u}"), x.get(f"{pre}ci_high_{u}"),
                                                        x.get(f"{pre}cox_low_power"), s))


def hr_face(hr, tested: bool, pale: bool, scale):
    """The colour of an HR cell: the scale's colour of log2 HR, most of the way to the neutral middle when the fit is
    imprecise; grey when untested."""
    from matplotlib.colors import to_rgba
    if not tested or hr is None or not np.isfinite(hr):
        return UNTESTED
    cmap, norm = scale
    c = np.asarray(cmap(norm(np.log2(hr))))
    return tuple(c + PALE * (np.asarray(to_rgba(NEUTRAL)) - c)) if pale else tuple(c)


def hatch(ax, rect, step: float = 0.34) -> None:
    """Thin diagonal lines inside `rect` (a Rectangle already on `ax`): the mark of an imprecise cell."""
    from matplotlib.collections import LineCollection
    x, y, w, h = rect.get_x(), rect.get_y(), rect.get_width(), rect.get_height()
    lines = LineCollection([[(x + k * w, y + h), (x + (k + 1) * w, y)] for k in np.arange(-1.0, 1.0 + 1e-9, step)],
                           colors=HATCH, linewidths=0.4, capstyle="butt")
    ax.add_collection(lines, autolim=False)
    lines.set_clip_path(rect)


def key_rows(items, width: float, size: float = 5.6) -> list[list[int]]:
    """The key's items (sample, label) laid out in rows at most `width` inches wide: lists of item indices. An item
    None starts a new row (between groups of items)."""
    rows, used = [[]], 0.0
    for i, item in enumerate(items):
        if item is None:
            if rows[-1]:
                rows.append([])
                used = 0.0
            continue
        w = KEY_W + 0.05 + S.text_width(item[1], size) + 0.16
        if rows[-1] and used + w > width:
            rows.append([])
            used = 0.0
        rows[-1].append(i)
        used += w
    return [r for r in rows if r]


def key(fig, x: float, y: float, width: float, W: float, H: float, items, size: float = 5.6) -> float:
    """A key to the cells' marks, drawn as the grids draw them: for each item a sample cell, then its label, in rows
    at most `width` inches wide from (x, y) inches from the figure's top left. `items`: (sample, label), the sample a
    dict of face (colour; None: no cell), dot ("small" or "large"), q (the q mark), frame, hatch (imprecise), value (a
    printed value and its colour) and draw (a function drawing on the sample's axes, 0-1 across and down). Returns
    the height used, in inches."""
    from matplotlib.patches import Rectangle
    rows = key_rows(items, width, size)
    for r, row in enumerate(rows):
        xi, yi = x, y + r * KEY_ROW
        for i in row:
            spec, label = items[i]
            ax = fig.add_axes([xi / W, 1 - (yi + KEY_H) / H, KEY_W / W, KEY_H / H])
            ax.set_xlim(0, 1)
            ax.set_ylim(1, 0)
            ax.axis("off")
            face = spec.get("face", NEUTRAL)
            if face is not None:
                rect = Rectangle((0, 0), 1, 1, facecolor=face, lw=0)
                ax.add_patch(rect)
                if spec.get("hatch"):
                    hatch(ax, rect)
            if spec.get("draw"):
                spec["draw"](ax)
            if spec.get("frame"):
                ax.add_patch(Rectangle((0, 0), 1, 1, fill=False, edgecolor=S.INK, lw=0.7, clip_on=False))
            ink = ink_on(face if face is not None else "white")
            if spec.get("dot"):
                ax.plot([0.5], [0.5], "o", ms=3.6 if spec["dot"] == "large" else 2.6, color=ink, mew=0)
            if spec.get("q"):
                ax.text(0.5, 0.62, S.Q_MARK, fontsize=8.5, color=ink, ha="center", va="center", fontweight="bold")
            if spec.get("value"):
                text, colour = spec["value"]
                ax.text(0.5, 0.5, text, fontsize=4.9, color=colour, ha="center", va="center")
            fig.text((xi + KEY_W + 0.05) / W, 1 - (yi + KEY_H / 2) / H, label, fontsize=size, color=S.INK2,
                     va="center")
            xi += KEY_W + 0.05 + S.text_width(label, size) + 0.16
    return len(rows) * KEY_ROW


def p_dots(face, s) -> list:
    """The key's items for the p dots: small for p < alpha, large for p < 0.01 (only the large one when alpha is
    0.01 or less, as then every dot is large)."""
    large = (dict(face=face, dot="large"), f"p < {min(s.alpha, 0.01):g}")
    return [large] if s.alpha <= 0.01 else [(dict(face=face, dot="small"), f"p < {s.alpha:g}"), large]


def hr_key(s, *, group: bool = True, q: bool = True, extra=()) -> list:
    """The key to the HR cells: p and q marks, group-hit frame, untested, and imprecise (with its rule)."""
    cmap, norm = hr_scale()
    mid = cmap(norm(-0.6))                              # a sample cell: a modest HR below 1
    items = p_dots(mid, s)
    if q and s.q_mark_below > 0:
        items.append((dict(face=mid, q=True), f"q < {s.q_mark_below:g}"))
    if group:
        items.append((dict(face=mid, frame=True), "group hit"))
    items.append((dict(face=UNTESTED), "not tested"))
    if s.imprecise_ci_ratio:
        power = f", or under {s.cox_low_power_events} events" if s.cox_low_power_events > 0 else ""
        items.append((dict(face=hr_face(2 ** -1.5, True, True, (cmap, norm)), hatch=True),
                      f"imprecise: 95% CI over {s.imprecise_ci_ratio:g}-fold{power}"))
    return items + list(extra)


def rho_text(v: float) -> str:
    """A correlation to print in a cell: two decimals without the leading zero, typographic minus ('−.42')."""
    t = f"{v:.2f}"
    t = {"1.00": "1", "-1.00": "−1", "-0.00": ".00"}.get(t, t.replace("0.", ".", 1))
    return t.replace("-", "−")


def colorbar(fig, x: float, y: float, w: float, W: float, H: float, scale, ticks, labels, title: str, ends=None):
    """A horizontal colour bar w inches wide, its top left at (x, y) inches from the figure's top left; `ends`: what
    its low and high ends mean, written above them."""
    cmap, norm = scale
    cax = fig.add_axes([x / W, 1 - (y + 0.08) / H, w / W, 0.08 / H])
    cax.imshow(np.linspace(norm.vmin, norm.vmax, 256)[None, :], aspect="auto", cmap=cmap, norm=norm,
               extent=(norm.vmin, norm.vmax, 0, 1))
    cax.set_yticks([])
    if w < BAR_FULL and len(ticks) > 3:                 # a narrow bar: the inner ticks only, so labels do not touch
        ticks, labels = ticks[1:-1], labels[1:-1]
    cax.set_xticks(ticks, labels, fontsize=5.6)
    cax.tick_params(length=2, pad=1.5)
    for sp in cax.spines.values():
        sp.set_visible(False)
    cax.set_xlabel(title, fontsize=5.8, labelpad=1, color=S.INK2)
    if ends and w >= BAR_FULL:
        for at, text, ha in ((0.0, f"← {ends[0]}", "left"), (1.0, f"{ends[1]} →", "right")):
            cax.text(at, 1.35, text, transform=cax.transAxes, ha=ha, va="bottom", fontsize=5.2, color=S.INK2)
    return cax


HR_TICKS = ([-HR_LIM, -1, 0, 1, HR_LIM], ["≤0.35", "0.5", "1", "2", "≥2.8"])   # the HR bar's ticks (log2) and labels
BAR_FULL = 1.4                      # inches: a colour bar at least this wide shows its end ticks and what its ends mean
HR_ENDS = ("lower hazard", "higher hazard")


class Legend(NamedTuple):
    """Where a grid figure's colour bar and key go, below the grid: the key beside the bar when it fits, else below."""
    bottom: float                       # inches below the grid, down to the figure's bottom edge
    key_x: float
    beside: bool


def legend_layout(W: float, x_bar: float, bar_w: float, items) -> Legend:
    """The space a colour bar at x_bar (bar_w wide) and the key to `items` need under a grid (figure W wide)."""
    key_x = x_bar + bar_w + 0.45
    beside = W - key_x - 0.12 >= 2.4
    key_x = key_x if beside else 0.12
    n = len(key_rows(items, W - key_x - 0.12))
    return Legend(max(0.62, 0.29 + n * KEY_ROW) if beside else 0.66 + n * KEY_ROW + 0.06, key_x, beside)


def bar_and_key(fig, W: float, H: float, lay: Legend, x_bar: float, bar_w: float, items, scale, ticks, labels,
                title: str, ends=None) -> None:
    """Draw the colour bar and the key where legend_layout put them (the grid ends lay.bottom above the bottom)."""
    g = H - lay.bottom                                  # the grid's bottom edge, in inches from the top
    colorbar(fig, x_bar, g + 0.29, bar_w, W, H, scale, ticks, labels, title, ends=ends)
    key(fig, lay.key_x, g + (0.23 if lay.beside else 0.66), W - lay.key_x - 0.12, W, H, items)


def clauses(text: str, width: float, size: float) -> list[str]:
    """Caption lines at most `width` inches wide, broken between clauses ('; ' outside parentheses), and between
    words only where a clause is wider than a line. Lines are measured on the glyph outlines; a raster at low dpi
    draws small text a little wider, so 5% is kept in reserve."""
    width *= 0.95
    parts, depth, cur = [], 0, ""
    for i, ch in enumerate(text):
        depth += (ch == "(") - (ch == ")")
        if ch == ";" and depth == 0 and text[i + 1:i + 2] == " ":
            parts.append(cur + ";")
            cur = ""
            continue
        cur += ch
    parts.append(cur)
    parts = [p.strip() for p in parts if p.strip()]
    lines, line = [], ""
    for p in parts:
        nxt = f"{line} {p}".strip()
        if line and S.text_width(nxt, size) > width:
            lines.append(line)
            line = p
        else:
            line = nxt
        while S.text_width(line, size) > width and " " in line:     # one clause wider than a line: by words
            words, head = line.split(" "), ""
            for k, w_ in enumerate(words):
                if head and S.text_width(f"{head} {w_}", size) > width:
                    lines.append(head)
                    line = " ".join(words[k:])
                    break
                head = f"{head} {w_}".strip()
            else:
                break
    return lines + ([line] if line else [])
