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
HIGHLIGHT = {"A": "#184f95", "B": "#6c9bd2", "C": "#c3d6ee"}
HIGHLIGHT_DEFAULT = "#6c9bd2"
DARK_LABELS = {"C"}                              # highlight letters drawn in ink rather than white
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
    """Width of a text in inches."""
    from matplotlib.font_manager import FontProperties
    from matplotlib.textpath import TextPath

    return TextPath((0, 0), text, size=size, prop=FontProperties(family=font_family(), **kw)).get_extents().width / 72


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
