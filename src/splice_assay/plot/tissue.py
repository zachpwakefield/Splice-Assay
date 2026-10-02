"""Case-vs-reference views: matched pairs beside all samples, all samples alone, or the reason there is no test."""
from __future__ import annotations

import numpy as np
from matplotlib import transforms
from matplotlib.ticker import MaxNLocator

from . import style as S

YTICKS = ([0, 0.25, 0.5, 0.75, 1.0], ["0", "", "0.5", "", "1"])
P0, P1, A0, A1, XLIM = 0.0, 1.0, 2.3, 3.3, (-0.38, 3.78)     # paired half at 0/1, all-samples half at 2.3/3.3


def _frame(ax, ylabel=True, quantity="PSI", lim=None):
    """PSI on 0-1; another quantity (e.g. expression) on `lim` with about four ticks."""
    if lim is None:
        ax.set_ylim(-0.02, 1.02)
        ax.set_yticks(*YTICKS)
    elif quantity == "HIT index":                         # the signed -1..1 scale, in halves
        ax.set_ylim(*lim)
        ax.set_yticks([-1, -0.5, 0, 0.5, 1], ["−1", "−0.5", "0", "0.5", "1"])
    else:
        ax.set_ylim(*lim)
        ax.yaxis.set_major_locator(MaxNLocator(4))
    if ylabel:
        ax.set_ylabel(quantity, labelpad=1)
    ax.yaxis.grid(True, color=S.GRID, lw=0.4)
    ax.set_axisbelow(True)


def _pairs(ax, x0, x1, ref, case):
    for r_, c_ in zip(ref, case):
        d = round(c_ - r_, 9)
        ax.plot([x0, x1], [r_, c_], color=S.UP if d > 0 else S.DOWN, lw=0.6, alpha=0.75, zorder=2 if d > 0 else 1,
                solid_capstyle="round")
    mr, mc = float(np.median(ref)), float(np.median(case))
    ax.plot([x0, x1], [mr, mc], color=S.INK, lw=1.7, zorder=4, solid_capstyle="round")
    ax.scatter([x0, x1], [mr, mc], s=14, color=S.INK, edgecolor="white", lw=0.5, zorder=5)
    return mr, mc


def _groups(ax, xs, ref, case, width=0.5, size=3):
    rng = np.random.default_rng(1)
    meds = {}
    for x0, (g, v, col) in zip(xs, (("reference", ref, S.NORMAL), ("case", case, S.TUMOUR))):
        ax.scatter(x0 + rng.uniform(-0.32 * width, 0.32 * width, len(v)), v, s=size, color=col, alpha=0.45, lw=0,
                   zorder=2)
        ax.boxplot([v], positions=[x0], widths=width, showfliers=False, patch_artist=True, zorder=3,
                   boxprops=dict(facecolor="none", edgecolor=S.INK2, lw=0.7), whiskerprops=dict(color=S.INK2, lw=0.6),
                   capprops=dict(color=S.INK2, lw=0.6), medianprops=dict(lw=0))
        meds[g] = float(np.median(v))
        ax.scatter([x0], [meds[g]], marker="D", s=15, color=col, edgecolor="white", lw=0.5, zorder=5)
    return meds


def _ticks(ax, fig, positions, labels, size=6.0):
    """Tick labels, rotated when neighbours would touch."""
    ax.set_xticks(positions, labels, fontsize=size)
    spacing = min(np.diff(positions)) / (ax.get_xlim()[1] - ax.get_xlim()[0]) * ax.get_position().width * \
        fig.get_figwidth()
    widest = max(S.text_width(line, size) for lab in labels for line in lab.split("\n"))
    if widest > 0.92 * spacing:
        for t in ax.get_xticklabels():
            t.set_rotation(35)
            t.set_ha("right")
            t.set_rotation_mode("anchor")
    ax.tick_params(axis="x", length=0, pad=2)


def paired_and_all(fig, ax, pairs_ref, pairs_case, ref, case, cutoff, labels, paired_title, all_title,
                   quantity="PSI", lim=None):
    """Left half: one line per pair; right half: all samples (points, boxes, median diamonds) and the KM split."""
    mr, mc = _pairs(ax, P0, P1, pairs_ref, pairs_case)
    meds = _groups(ax, (A0, A1), ref, case, width=0.5, size=2.2)
    ax.axvline((P1 + A0) / 2, color=S.GRID, lw=0.7, zorder=0)
    if np.isfinite(cutoff):
        ax.plot([A0 - 0.45, XLIM[1]], [cutoff, cutoff], color=S.INK, lw=0.6, ls=(0, (2, 1.5)), zorder=4)
    ax.set_xlim(*XLIM)
    _frame(ax, quantity=quantity, lim=lim)
    rl, cl = labels["reference"], labels["case"]
    _ticks(ax, fig, [P0, P1, A0, A1], [rl, cl, f"{rl}\n{len(ref)}", f"{cl}\n{len(case)}"])
    tr = transforms.blended_transform_factory(ax.transData, ax.transAxes)
    for x0, text in (((P0 + P1) / 2, paired_title), ((A0 + A1) / 2, all_title)):
        ax.text(x0, 1.02, text, transform=tr, ha="center", va="bottom", fontsize=6.0, color=S.INK2, linespacing=1.15)
    return (mr, mc), meds


def all_samples(fig, ax, ref, case, cutoff, labels, title, note, quantity="PSI", lim=None):
    """All samples only (no paired test): points, boxes, median diamonds, the KM split and the reason."""
    meds = _groups(ax, (0, 1), ref, case)
    if np.isfinite(cutoff):
        ax.plot([-0.5, 1.5], [cutoff, cutoff], color=S.INK, lw=0.6, ls=(0, (2, 1.5)), zorder=4)
    ax.set_xlim(-0.6, 1.6)
    _frame(ax, quantity=quantity, lim=lim)
    _ticks(ax, fig, [0, 1], [f"{labels['reference']}\n{len(ref)}", f"{labels['case']}\n{len(case)}"], size=6.5)
    ax.text(0.5, 1.02, title, transform=ax.transAxes, ha="center", va="bottom", fontsize=6.0, color=S.INK2,
            linespacing=1.15)
    ax.annotate(note, xy=(0.5, 0), xycoords="axes fraction", xytext=(0, -24), textcoords="offset points",
                ha="center", va="top", fontsize=5.6, color=S.MUTED, linespacing=1.15)
    return meds


def empty(ax, message: str, quantity="PSI", lim=None):
    """No case-vs-reference test in this cohort: the axis carries the reason."""
    ax.set_xlim(0, 1)
    ax.set_xticks([])
    _frame(ax, quantity=quantity, lim=lim)
    ax.text(0.5, 0.5, message, transform=ax.transAxes, ha="center", va="center", fontsize=5.8, color=S.MUTED,
            linespacing=1.2)
