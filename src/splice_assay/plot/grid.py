"""Pieces shared by the probe's summary figures (the overview, the gene maps and the correlation grid): their colour
scales, the fit a survival cell shows, colour bars, and captions wrapped between their clauses."""
from __future__ import annotations

import numpy as np

from . import style as S

HR_COLORS = ["#2f5f98", "#9ebbd9", "#f4f3ee", "#e6a88a", "#b0412c"]    # log2 HR, -HR_LIM to HR_LIM
RHO_COLORS = ["#5e3c99", "#b2abd2", "#f4f3ee", "#a6dba0", "#1b7837"]   # Spearman rho, -1 (purple) to 1 (green)
HR_LIM = 1.5
UNTESTED = "#e8e7e1"


def hr_scale():
    """(colormap, norm) of the HR cells: log2 HR from -1.5 to 1.5."""
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    return LinearSegmentedColormap.from_list("hr", HR_COLORS), Normalize(-HR_LIM, HR_LIM)


def rho_scale():
    """(colormap, norm) of the correlation cells: rho from -1 to 1."""
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    return LinearSegmentedColormap.from_list("rho", RHO_COLORS), Normalize(-1, 1)


def flag(v) -> bool:
    """A True/False cell that may be missing (NaN, None) or a numpy bool."""
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


def shown_fit(x, s) -> tuple[bool, float, float, float]:
    """(tested, hr, p, q) of the fit a survival cell shows: the adjusted model where it was fitted, else the base
    model; the HR per Settings.psi_hr_unit. `x`: a row of the probe's cells."""
    adj = "adj_cox_status" in x and x.adj_cox_status == "tested"
    u = s.psi_hr_unit
    hr, p = (x[f"adj_hr_per_{u}"], x.adj_cox_p) if adj else (x.get(f"hr_per_{u}"), x.get("cox_p"))
    q = x.get("adj_cox_q", np.nan) if adj else x.get("cox_q", np.nan)
    return bool(adj or x.cox_status == "tested"), hr, p, q


def rho_text(v: float) -> str:
    """A correlation to print in a cell: two decimals without the leading zero, typographic minus ('−.42')."""
    t = f"{v:.2f}"
    t = {"1.00": "1", "-1.00": "−1", "-0.00": ".00"}.get(t, t.replace("0.", ".", 1))
    return t.replace("-", "−")


def colorbar(fig, x: float, y: float, w: float, W: float, H: float, scale, ticks, labels, title: str):
    """A horizontal colour bar w inches wide, its top left at (x, y) inches from the figure's top left."""
    cmap, norm = scale
    cax = fig.add_axes([x / W, 1 - (y + 0.08) / H, w / W, 0.08 / H])
    cax.imshow(np.linspace(norm.vmin, norm.vmax, 256)[None, :], aspect="auto", cmap=cmap, norm=norm,
               extent=(norm.vmin, norm.vmax, 0, 1))
    cax.set_yticks([])
    cax.set_xticks(ticks, labels, fontsize=5.6)
    cax.tick_params(length=2, pad=1.5)
    for sp in cax.spines.values():
        sp.set_visible(False)
    cax.set_xlabel(title, fontsize=5.8, labelpad=1, color=S.INK2)
    return cax


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
