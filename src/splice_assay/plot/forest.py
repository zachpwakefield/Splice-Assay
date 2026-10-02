"""Twin forest: per cohort, the case-vs-reference delta (paired filled, unpaired open) beside the Cox HR per IQR (or
per SD) of PSI. Without reference samples in the data only the HR side is drawn.

Rows are cohorts in which a plotted event had a tested comparison or a Cox fit (tested or failed), plus the cohorts
shown at left. With two events each cohort has two sub-rows (first event ink, second grey).
"""
from __future__ import annotations

import numpy as np
from matplotlib import ticker as mticker
from matplotlib.patches import Rectangle

from . import style as S

NICE_LO, NICE_HI = (0.125, 0.25, 0.5), (2.0, 4.0, 8.0)
PH_MARK = "†"                               # a proportional-hazards test below Settings.ph_note_below


def axis_limits(recs_by_cohort: dict) -> tuple[float, float, float, int]:
    """Delta axis half-width; HR axis limits: every shaded cell's CI in full, other CIs to their 10th/90th
    percentile, snapped to 1/8 ... 8; beyond that CIs are clipped with an arrowhead (counted)."""
    recs = [r for rs in recs_by_cohort.values() for r in rs]
    vals = np.array([r[k] for r in recs for k in ("paired_delta", "unpaired_delta")], float)
    lim = max(0.3, float(np.nanmax(np.abs(vals))) * 1.12) if np.isfinite(vals).any() else 0.3
    his = [r["hi"] for r in recs if np.isfinite(r["hi"])]
    los = [r["lo"] for r in recs if np.isfinite(r["lo"])]
    his_t = [r["hi"] for r in recs if np.isfinite(r["hi"]) and r["shade"]]
    los_t = [r["lo"] for r in recs if np.isfinite(r["lo"]) and r["shade"]]
    lo_target = min([np.percentile(los, 10)] + los_t) if los else 0.5
    hi_target = max([np.percentile(his, 90)] + his_t) if his else 2.0
    xlo = max([v for v in NICE_LO if v <= lo_target] or [NICE_LO[0]])
    xhi = min([v for v in NICE_HI if v >= hi_target] or [NICE_HI[-1]])
    n_clip = int(sum(v < xlo for v in los) + sum(v > xhi for v in his))
    return lim, xlo, xhi, n_clip


def set_axes(axL, axR, lim, xlo, xhi, min_abs_delta, labels, model_note="", quantity="PSI",
             hr_label="HR per IQR\nof PSI (95% CI)"):
    """`quantity` names the delta axis (PSI, HIT index, expression); `hr_label` the HR axis."""
    if axL is not None:
        axL.set_xlim(-lim, lim)
        axL.xaxis.set_major_locator(mticker.MultipleLocator(0.2) if lim <= 0.6 else
                                    mticker.MaxNLocator(nbins=5, steps=[1, 2, 2.5, 5, 10]))
        for v in (-min_abs_delta, min_abs_delta):
            axL.axvline(v, color=S.MUTED, lw=0.5, ls=(0, (1.5, 1.5)), zorder=1)
        axL.axvline(0, color=S.INK2, lw=0.6, zorder=1)
        axL.set_xlabel(f"Δ median {quantity}\n({S.in_sentence(labels['case'])} − "
                       f"{S.in_sentence(labels['reference'])})", labelpad=1)
    axR.set_xscale("log")
    axR.set_xlim(xlo, xhi)
    tk = [t for t in (0.25, 0.5, 1, 2, 4, 8) if xlo <= t <= xhi]          # no 0.125 label (crowds 0.25)
    axR.set_xticks(tk, [f"{t:g}" for t in tk])
    axR.xaxis.set_minor_locator(mticker.NullLocator())
    axR.axvline(1, color=S.INK2, lw=0.6, zorder=1)
    axR.set_xlabel(hr_label, labelpad=1)
    if model_note:
        axR.annotate(model_note, xy=(0.5, 0), xycoords="axes fraction", xytext=(0, -30), textcoords="offset points",
                     ha="center", va="top", fontsize=5.4, color=S.MUTED, linespacing=1.15)


def draw(fig, axL, axR, lab_x, mark_x, cohorts, recs_by_cohort, n_ev, alpha, cox_name, colors=None,
         ph_below: float = 0.0, q_below: float = 0.0) -> dict:
    """recs_by_cohort[cohort][i]: values of event i (deltas, HR, CI, p, statuses, mark label; cox_q and ph_p
    optional). axL may be None. After a CI, * when its q is below q_below (the fill shows p < alpha) and a dagger when
    its proportional-hazards p is below ph_below; returns which of the two were drawn ({"q": bool, "ph": bool})."""
    colors = colors or {}
    marked = dict(q=False, ph=False)
    n = len(cohorts)
    offs = [0.0] if n_ev == 1 else [-0.24, 0.24]
    shades = [S.INK] if n_ev == 1 else [S.INK, S.SECOND]
    pos = axR.get_position()
    W = fig.get_figwidth()
    axes = [a for a in (axL, axR) if a is not None]
    for j, co in enumerate(cohorts):
        y = n - 1 - j
        recs = recs_by_cohort[co]
        shaded = any(r["shade"] for r in recs)
        if shaded:
            for ax in axes:
                ax.axhspan(y - 0.5, y + 0.5, color=S.PAPER, lw=0, zorder=0)
        yf = pos.y0 + (y + 0.5) / n * pos.height
        fig.text(lab_x, yf, co, fontsize=6.0, ha="right", va="center", color=S.INK if shaded else S.INK2,
                 fontweight="bold" if shaded else "normal")
        for i, r in enumerate(recs):
            yy = y - offs[i]
            col = shades[i]
            if r["mark"]:
                h = 0.8 / n if n_ev == 1 else 0.44 / n
                yb = pos.y0 + (yy + 0.5) / n * pos.height
                bw, bh = 0.09 / W, h * pos.height
                face = colors.get(r["mark"], S.HIGHLIGHT.get(r["mark"], S.HIGHLIGHT_DEFAULT))
                fig.add_artist(Rectangle((mark_x - bw / 2, yb - bh / 2), bw, bh, transform=fig.transFigure,
                                         facecolor=face, lw=0))
                fig.text(mark_x, yb, r["mark"], fontsize=4.9 if n_ev == 2 else 5.6, ha="center", va="center",
                         color=S.INK if r["mark"] in S.DARK_LABELS else "white", fontweight="bold")
            if axL is not None:
                for design, face in (("paired", col), ("unpaired", "white")):
                    v = r[f"{design}_delta"]
                    if np.isfinite(v):
                        axL.scatter([v], [yy], s=13 if n_ev == 1 else 10, facecolor=face, edgecolor=col, lw=0.7,
                                    zorder=3)
            if np.isfinite(r["hr"]):
                lo, hi = r["lo"], r["hi"]
                xlo, xhi = axR.get_xlim()
                axR.plot([max(lo, xlo), min(hi, xhi)], [yy, yy], color=col, lw=0.8, zorder=2, solid_capstyle="butt")
                for edge, cond, mk in ((xlo, lo < xlo, "<"), (xhi, hi > xhi, ">")):
                    if cond:
                        axR.plot([edge], [yy], marker=mk, ms=2.6, color=col, mew=0, zorder=3, clip_on=False)
                axR.scatter([min(max(r["hr"], xlo), xhi)], [yy], s=15 if n_ev == 1 else 11, marker="D",
                            facecolor=col if r["cox_p"] < alpha else "white", edgecolor=col, lw=0.7, zorder=3)
                q, ph = r.get("cox_q", np.nan), r.get("ph_p", np.nan)
                after = [m for m, on in ((S.Q_MARK, q_below > 0 and q < q_below),
                                         (PH_MARK, ph_below > 0 and ph < ph_below)) if on]
                if after:
                    axR.annotate(" ".join(after), (min(hi, xhi), yy), xytext=(3, 0), textcoords="offset points",
                                 fontsize=6.0, color=col, ha="left", va="center", annotation_clip=False)
                    marked["q"] |= S.Q_MARK in after
                    marked["ph"] |= PH_MARK in after
            elif r["cox_status"] == "failed":
                axR.text(1.0, yy, f"{cox_name} fit failed", fontsize=5.2, color=S.MUTED, ha="center", va="center")
    for ax in axes:
        ax.set_ylim(-0.5, n - 0.5)
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.xaxis.grid(True, color=S.GRID, lw=0.4)
        ax.set_axisbelow(True)
    return marked
