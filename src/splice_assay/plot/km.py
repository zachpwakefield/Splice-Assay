"""Kaplan-Meier panel: the high/low split (median, mean or a set value), 95% log-log bands, censoring ticks and an
at-risk table."""
from __future__ import annotations

import numpy as np
import pandas as pd
from matplotlib import transforms
from matplotlib.lines import Line2D

from ..stats.core import km_with_ci
from . import style as S


def _curve(t, e, z):
    c = km_with_ci(t, e, z)
    c = pd.DataFrame({"t": np.r_[0.0, c.time], "s": np.r_[1.0, c.surv], "lo": np.r_[1.0, c.lo], "hi": np.r_[1.0, c.hi]})
    return pd.concat([c, pd.DataFrame(dict(t=[t.max()], s=[c.s.iloc[-1]], lo=[c.lo.iloc[-1]], hi=[c.hi.iloc[-1]]))],
                     ignore_index=True)


SHORT = {"expression": "expr."}
PH_MARK = "†"                               # marks a proportional-hazards test below Settings.ph_note_below


def split_text(row, what: str) -> str:
    """'split at median PSI 0.7705' (or mean), or 'split at PSI 0.5 (set)' for a split given in the settings."""
    how = row.get("km_split") or "median"
    how = how if isinstance(how, str) else "median"
    if how == "set":
        return f"split at {what} {S.fcut(row['cutoff'])} (set)"
    return f"split at {how} {what} {S.fcut(row['cutoff'])}"


def draw(fig, ax, t_years, event, high, row: dict, ylabel: str, rows: list, tag: str, settings, what: str = "PSI"):
    """row: the survival row of the cell (cutoff, logrank_hr, km_p, events_high, events_low); `what` names the split
    variable (PSI, expression; the curve labels use a short form)."""
    s = settings
    t, e = np.asarray(t_years, float), np.asarray(event, int)
    tmax = float(min(s.km_max_years, s.km_tick_years * np.ceil(t.max() / s.km_tick_years)))
    ticks = [float(x) for x in np.arange(0, s.km_max_years + 1e-9, s.km_tick_years) if x <= tmax + 1e-9]
    risk = {}
    low_lab, high_lab = f"Low {SHORT.get(what, what)}", f"High {SHORT.get(what, what)}"
    for arm, col, lab in ((False, S.KM_LOW, low_lab), (True, S.KM_HIGH, high_lab)):
        ta, ea = t[high == arm], e[high == arm]
        cu = _curve(ta, ea, s.ci_z)
        # arrays, not Series: with matplotlib 3.7.0 and numpy 1.24, fill_between fails on a pandas >= 2.1 Series
        ax.fill_between(cu.t.to_numpy(), cu.lo.to_numpy(), cu.hi.to_numpy(), step="post", color=col, alpha=0.13,
                        lw=0, zorder=1)
        ax.step(cu.t.to_numpy(), cu.s.to_numpy(), where="post", color=col, lw=1.3, zorder=3)
        cens = ta[(ea == 0) & (ta <= tmax)]
        ax.scatter(cens, [cu.s[cu.t <= c].iloc[-1] for c in cens], marker="|", s=7, color=col, lw=0.5, zorder=4)
        send = float(cu.s[cu.t <= tmax].iloc[-1])
        risk[lab] = (col, [int((ta >= x).sum()) for x in ticks], send)
        rows += [dict(panel="km_curve", tag=tag, arm=lab, time_years=a, surv=b, ci_low=c, ci_high=d)
                 for a, b, c, d in cu[["t", "s", "lo", "hi"]].itertuples(index=False)]
        rows += [dict(panel="km_at_risk", tag=tag, arm=lab, time_years=x, n_at_risk=n) for x, n in zip(ticks, risk[lab][1])]
    # curve-end labels, nudged apart when close
    ends = {lab: v[2] for lab, v in risk.items()}
    dn, up = sorted(ends, key=ends.get)                 # two labels even when the curves end level (e.g. both at 0)
    if ends[up] - ends[dn] < 0.09:
        mid = (ends[up] + ends[dn]) / 2
        ends[up], ends[dn] = mid + 0.045, mid - 0.045
    if ends[dn] < 0.04:
        ends[up], ends[dn] = ends[up] + 0.04 - ends[dn], 0.04
    for lab in risk:
        ax.text(tmax * 1.02, ends[lab], lab, fontsize=6, color=S.INK2, va="center", ha="left", clip_on=False)
    ax.set_xlim(0, tmax)
    ax.set_ylim(0, 1.02)
    ax.set_xticks(ticks, [f"{x:g}" for x in ticks])
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0], ["0", "", "0.5", "", "1"])
    ax.yaxis.grid(True, color=S.GRID, lw=0.4)
    ax.set_axisbelow(True)
    ax.set_xlabel("Years", labelpad=1)
    ax.set_ylabel(ylabel, labelpad=1)
    q = row.get("km_q", np.nan)
    q = np.nan if q is None else float(q)
    ph = row.get("km_ph_p", np.nan)
    ph = np.nan if ph is None else float(ph)
    ph_line = s.ph_note_below > 0 and np.isfinite(ph) and ph < s.ph_note_below
    ax.text(0.0, 1.03, split_text(row, what) + "\n"
            f"HR {row['logrank_hr']:.2f} (high vs low) · log-rank p {S.fp(row['km_p'])}" + (f" {PH_MARK}" if ph_line
                                                                                         else "") + "\n"
            + (S.q_text(q, s.q_mark_below) + "\n" if np.isfinite(q) else "") +
            f"events {int(row['events_high'])} high, {int(row['events_low'])} low"
            + (f"\n{PH_MARK} non-proportional hazards (p {S.fp(ph)})" if ph_line else ""), transform=ax.transAxes,
            fontsize=6.0, va="bottom", ha="left", color=S.INK2, linespacing=1.2)
    # at-risk table under the axis
    pos = ax.get_position()
    W, H = fig.get_figwidth(), fig.get_figheight()
    tr = transforms.blended_transform_factory(ax.transData, fig.transFigure)
    lh = 0.105 / H
    ytab = pos.y0 - 0.36 / H
    xl = pos.x0 - 0.16 / W
    fig.text(xl, ytab + lh * 0.95, "At risk", fontsize=5.6, color=S.INK2, ha="right", va="center")
    for j, lab in enumerate((high_lab, low_lab)):
        col, nr, _ = risk[lab]
        y = ytab - j * lh
        fig.text(xl, y, lab.split()[0], fontsize=5.6, color=S.INK2, ha="right", va="center")
        fig.add_artist(Line2D([xl - 0.36 / W, xl - 0.24 / W], [y, y], color=col, lw=1.4, transform=fig.transFigure))
        for x, n in zip(ticks, nr):
            ax.text(x, y, str(n), transform=tr, fontsize=5.6, color=S.INK2, ha="center", va="center", clip_on=False)
    return risk


def empty(ax, message: str, ylabel: str):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.02)
    ax.set_xticks([])
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0], ["0", "", "0.5", "", "1"])
    ax.set_ylabel(ylabel, labelpad=1)
    ax.text(0.5, 0.5, message, transform=ax.transAxes, ha="center", va="center", fontsize=5.8, color=S.MUTED,
            linespacing=1.2)
