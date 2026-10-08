"""Palette, font, number formats and deterministic saving (SVG with live text, PDF with embedded fonts, PNG)."""
from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

INK, INK2, MUTED, GRID, PAPER = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#f4f3ee"
TUMOUR, NORMAL = "#2a78d6", "#eb6834"
UP, DOWN = "#2a78d6", "#a8a69f"                  # paired lines: higher / lower in tumour
KM_LOW, KM_HIGH = "#86b6ef", "#184f95"
SECOND = MUTED                                  # the second event of a two-event figure
EXON, EXON_LINE, NESTED = "#a3a19a", "#8f8d86", "#1f1f1f"
HIGHLIGHT = {"A": "#184f95", "B": "#6c9bd2", "C": "#c3d6ee",          # tier letters, and the probe's evidence grades
             "A+": "#184f95", "B+": "#6c9bd2", "C+": "#c3d6ee", "D": "#dcd8cc", "E": "#e8b9b1"}
HIGHLIGHT_DEFAULT = "#6c9bd2"
DARK_LABELS = {"C", "C+", "D", "E"}              # highlight letters drawn in ink rather than white
TYPE_COLOR = {"SE": "#3b6fb6", "MXE": "#7b5ea7", "RI": "#b8474d", "A3SS": "#2e8b57", "A5SS": "#c98a2b",
              "AFE": "#178f8f", "ALE": "#b5527f", "HIT": "#5b6b7f"}
TYPE_DEFAULT = "#555555"
FONT_CANDIDATES = ("Arial", "Liberation Sans", "Helvetica", "Nimbus Sans", "DejaVu Sans")


@lru_cache(maxsize=1)
def font_family() -> str:
    """Arial when installed, else the closest metric-compatible sans serif."""
    from matplotlib import font_manager

    for name in FONT_CANDIDATES:
        try:
            font_manager.findfont(font_manager.FontProperties(family=name), fallback_to_default=False)
            return name
        except ValueError:
            continue
    return "DejaVu Sans"


def rc() -> dict:
    return {
        "font.family": font_family(), "font.size": 7, "axes.labelsize": 7, "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5, "axes.linewidth": 0.6, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5, "axes.spines.top": False, "axes.spines.right": False,
        "axes.edgecolor": INK2, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
        "svg.fonttype": "none", "pdf.fonttype": 42, "legend.frameon": False, "svg.hashsalt": "splice-assay",
        "path.simplify": True, "figure.dpi": 100,
    }


def in_sentence(name: str) -> str:
    """A group name inside a sentence: 'Tumour' -> 'tumour', but 'IDH-mutant', 'PD-L1 high' or 'Non-Responder'
    keep their capitals (only a lone leading capital is lowered)."""
    name = str(name)
    return name if any(ch.isupper() for ch in name[1:]) else name[:1].lower() + name[1:]


def text_width(text: str, size: float, **kw) -> float:
    """Width of a text in inches: the extent of its glyph outlines, as Path.get_extents gives it, with each curve's
    extremes in closed form and no per-segment objects (get_extents is ~30x slower; the layout measures thousands of
    texts). Blank text: 0."""
    text = str(text)
    return _text_width(text, float(size), font_family(), tuple(sorted(kw.items()))) if text.strip() else 0.0


@lru_cache(maxsize=16384)
def _text_width(text: str, size: float, family: str, kw: tuple) -> float:
    from matplotlib.font_manager import FontProperties
    from matplotlib.textpath import TextPath

    path = TextPath((0, 0), text, size=size, prop=FontProperties(family=family, **dict(kw)))
    xs = _x_points(path.vertices[:, 0].tolist(), path.codes.tolist() if path.codes is not None else None)
    return (max(xs) - min(xs)) / 72 if xs else 0.0


def _x_points(x: list, codes) -> list:
    """The x of a path's segment ends and of its curves' extremes (quadratic and cubic Beziers)."""
    from matplotlib.path import Path as MPath

    if codes is None:
        return list(x)
    out, i, prev = [], 0, 0.0
    while i < len(x):
        c = codes[i]
        if c in (MPath.MOVETO, MPath.LINETO):
            prev = x[i]
            out.append(prev)
            i += 1
        elif c == MPath.CURVE3:                         # p0 (prev), control p1, end p2
            p0, p1, p2 = prev, x[i], x[i + 1]
            den = p0 - 2 * p1 + p2
            if den:
                s = (p0 - p1) / den
                if 0 < s < 1:
                    out.append((1 - s) ** 2 * p0 + 2 * (1 - s) * s * p1 + s * s * p2)
            prev = p2
            out.append(prev)
            i += 2
        elif c == MPath.CURVE4:                         # p0 (prev), controls p1 p2, end p3
            p0, p1, p2, p3 = prev, x[i], x[i + 1], x[i + 2]
            a, b, d = -p0 + 3 * p1 - 3 * p2 + p3, 2 * (p0 - 2 * p1 + p2), p1 - p0   # x'(s)/3 = a s^2 + b s + d
            roots = ([-d / b] if b else []) if not a else (
                [(-b + sg * (b * b - 4 * a * d) ** 0.5) / (2 * a) for sg in (1, -1)] if b * b >= 4 * a * d else [])
            for s in roots:
                if 0 < s < 1:
                    out.append((1 - s) ** 3 * p0 + 3 * (1 - s) ** 2 * s * p1 + 3 * (1 - s) * s * s * p2 + s ** 3 * p3)
            prev = p3
            out.append(prev)
            i += 3
        else:                                           # CLOSEPOLY, STOP: no point of the outline
            i += 1
    return out


def fp(p) -> str:
    """p to print: 2 decimals from 0.1, 3 below, 1.2e-5 style below 0.001."""
    if p is None or not np.isfinite(p):
        return "–"
    if p < 1e-3:
        m, e = f"{p:.1e}".split("e")
        return f"{m}e{int(e)}"
    return f"{p:.3f}" if p < 0.1 else f"{p:.2f}"


Q_MARK = "*"                                    # follows a q value below Settings.q_mark_below


def q_text(q, below: float) -> str:
    """'q 0.003 *' (the mark when q is below `below`), 'q 0.21', or '' when there is no q."""
    if q is None or not np.isfinite(float(q)):
        return ""
    return f"q {fp(float(q))}" + (f" {Q_MARK}" if below > 0 and float(q) < below else "")


def fd(v) -> str:
    """A delta PSI to print: up to 4 decimals, typographic minus."""
    t = f"{v:+.4f}".rstrip("0").rstrip(".")
    return t.replace("-", "−")


def fcut(v) -> str:
    return f"{v:.4f}".rstrip("0").rstrip(".")


META = {"svg": {"Date": None}, "pdf": {"CreationDate": None, "ModDate": None}, "png": {}}


def save(fig, out_dir, stem: str, formats=("svg", "pdf", "png"), dpi: int = 400, table: pd.DataFrame | None = None,
         provenance: dict | None = None) -> dict[str, Path]:
    """Write the figure in each format, the plotted values as CSV and the provenance as JSON; byte-for-byte
    reproducible for the same inputs and library versions."""
    import matplotlib

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    paths = {}
    with matplotlib.rc_context(rc()):
        for ext in formats:
            p = out / f"{stem}.{ext}"
            fig.savefig(p, format=ext, metadata=META[ext], **({"dpi": dpi} if ext == "png" else {}))
            paths[ext] = p
    if table is not None:
        paths["csv"] = out / f"{stem}.csv"
        table.to_csv(paths["csv"], index=False, lineterminator="\n")
    if provenance is not None:
        paths["json"] = out / f"{stem}.provenance.json"
        paths["json"].write_text(json.dumps(provenance, indent=2, sort_keys=True, default=str) + "\n")
    return paths
