"""The Cox model of one event in one cohort: every term's hazard ratio with its 95% CI (a multivariable forest) beside
the numbers, and the coefficient table as CSV. `draw_terms` draws the same forest inside another figure (the event
panel's model band)."""
from __future__ import annotations

import numpy as np
import pandas as pd

from .. import events as EV
from .. import provenance as P
from ..analysis import Results, analyze
from ..config import Settings
from ..dataset import Dataset, InputError
from ..stats.survival import CoxModel
from . import style as S
from .panel_common import ENDPOINT_NAMES, Panel, model_mismatch, safe_name

MODEL_NOTE = "Filled: p < {alpha:g}. Numeric covariates per SD of the fit cohort; categories against the level named."
MAIN = ("psi_iqr", "psi_sd", "gex")                      # the tested term of a model: drawn bold, with a diamond


def model_terms(ds: Dataset, event: str, cohort: str, endpoint: str, model: CoxModel | None = None,
                settings: Settings | None = None, results: Results | None = None) -> tuple[pd.Series, pd.DataFrame]:
    """The survival row and the Cox terms of one cell (raises when the model was not fitted)."""
    res = results if results is not None else analyze(ds, events=[event], cohorts=[cohort], endpoints=[endpoint],
                                                     settings=settings, model=model)
    sv = res.survival
    row = sv[sv.event_id.eq(event) & sv.cohort.eq(cohort) & sv.endpoint.eq(endpoint)]
    if row.empty:
        raise InputError(f"no survival row for {event} in {cohort} ({endpoint})")
    row = row.iloc[0]
    if row.cox_status != "tested":
        raise InputError(f"the Cox model of {event} in {cohort} ({endpoint}) was not fitted: "
                         f"{not_fitted(row, res.settings)}")
    t = res.cox_terms
    t = t[t.event_id.eq(event) & t.cohort.eq(cohort) & t.endpoint.eq(endpoint)].reset_index(drop=True)
    return row, t


def not_fitted(row, s: Settings | None = None) -> str:
    """Why a Cox model was not fitted, in words, with the counts that decided it."""
    status, ev = str(row.get("cox_status", "")), row.get("cox_events", np.nan)
    if status == "too_few_events_per_term" and pd.notna(row.get("cox_n_terms", np.nan)):
        need = f"; needs {s.cox_min_events_per_term:g} per term" if s is not None else ""
        k = int(row["cox_n_terms"])
        return f"too few events per term ({int(ev)} events for {k} term{'s' if k != 1 else ''}{need})"
    n = row.get("cox_n", np.nan)
    return status.replace("_", " ") + ("" if pd.isna(n) else f" ({int(n)} patients, {int(ev)} events)")


def display_rows(terms: pd.DataFrame, q: float = np.nan) -> list[dict]:
    """The rows of a model forest: PSI per SD (or per IQR, Settings.psi_hr_unit), host expression and numeric
    covariates per SD, and each categorical covariate as a header with one indented row per level. `q`: the tested
    PSI term's q, kept on its row."""
    out = []
    for r in terms.itertuples():
        if r.kind == "psi":
            continue                                     # PSI is shown per IQR or SD (its per +0.10 row is in the CSV)
        if r.kind in ("psi_iqr", "psi_sd"):
            u = r.kind[4:].upper()
            out.append(dict(label=f"{r.term} (per {u}, {r.unit[len(u) + 2:-1]})", short=f"{r.term} (per {u})",
                            indent=False, r=r, q=q))
        elif r.kind == "expression":
            out.append(dict(label="Host expression (per SD)", short="Host expr. (per SD)", indent=False, r=r))
        elif r.kind == "gex":                            # the main term of an expression model
            out.append(dict(label="Expression (per SD)", short="Expr. (per SD)", indent=False, r=r))
        elif r.kind == "numeric":
            unit = "SD" if r.unit == "SD" else "unit"
            out.append(dict(label=f"{r.term} (per {unit})", short=f"{r.term} (per {unit})", indent=False, r=r))
        else:
            if not out or out[-1].get("head") != r.term:
                out.append(dict(label=f"{r.term} (vs {r.reference})", short=f"{r.term} (vs {r.reference})",
                                indent=False, r=None, head=r.term))
            out.append(dict(label=str(r.level), short=str(r.level), indent=True, r=r, head=r.term))
    return out


def axis_range(vals) -> tuple[float, float]:
    """An HR axis covering every estimate and most CIs (10th/90th percentile of their ends), snapped to powers of 2
    within 1/16..16; a CI beyond it ends in an arrowhead."""
    hr = [r.hr for r in vals if np.isfinite(r.hr)]
    los = [r.ci_low for r in vals if np.isfinite(r.ci_low)]
    his = [r.ci_high for r in vals if np.isfinite(r.ci_high)]
    lo = min([0.9] + hr + ([float(np.percentile(los, 10))] if los else []))
    hi = max([1.1] + hr + ([float(np.percentile(his, 90))] if his else []))
    return 2.0 ** max(-4.0, np.floor(np.log2(lo))), 2.0 ** min(4.0, np.ceil(np.log2(hi)))


_range = axis_range


def _ticks(xlo, xhi, plot_w, fs) -> list[float]:
    """Powers of 2 on the axis, every other one (keeping 1) while their labels would crowd."""
    ticks = [2.0 ** k for k in range(int(np.log2(xlo)), int(np.log2(xhi)) + 1)]
    step = 1
    while True:
        keep = [t for t in ticks if round(np.log2(t)) % step == 0]
        need = sum(S.text_width(f"{t:g}", fs) for t in keep) + 0.1 * (len(keep) - 1)
        if need <= 0.9 * plot_w or len(keep) <= 3:
            return keep
        step *= 2


def p_text(d: dict, ph_below: float = 0.0, q_below: float = 0.0) -> str:
    """A model row's p, then * when the tested term's q (d["q"], set for that row only) is below q_below, and a dagger
    when the term's proportional-hazards p is below ph_below."""
    r = d["r"]
    q, ph = d.get("q", np.nan), float(getattr(r, "ph_p", np.nan))
    after = [m for m, on in ((S.Q_MARK, q_below > 0 and q is not None and q < q_below),
                             ("†", ph_below > 0 and ph < ph_below)) if on]
    return S.fp(r.p) + (" " + " ".join(after) if after else "")


def terms_columns(width: float, disp_lists: list[list[dict]], fs: float, ph_below: float = 0.0,
                  q_below: float = 0.0) -> dict:
    """One column layout for one or several model forests (so that they align): which labels (full or short), the
    font size, and the widths of the label, plot, HR and p columns (the p column holds the marks after p)."""
    disp = [d for dl in disp_lists for d in dl]
    vals = [d["r"] for d in disp if d["r"] is not None]
    if not vals:
        return dict(key="label", fs=fs, lw=0.0, hw=0.0, pw=0.0, plot_w=0.45, g=0.08)
    hr_text = [f"{r.hr:.2f} ({r.ci_low:.2f}–{r.ci_high:.2f})" for r in vals] + ["HR (95% CI)"]
    g = 0.08
    for attempt in range(3):                         # labels, then shorter labels, then a smaller font
        key = "label" if attempt == 0 else "short"
        lw = max(S.text_width(d[key], fs, weight="bold" if d["r"] is not None and d["r"].kind in MAIN
                              else "normal") + (0.14 if d["indent"] else 0) for d in disp)
        hw = max(S.text_width(t, fs - 0.2) for t in hr_text)
        pw = max([S.text_width(p_text(d, ph_below, q_below), fs - 0.2) for d in disp if d["r"] is not None]
                 + [0.1])
        plot_w = width - lw - hw - pw - 3 * g
        if plot_w >= 0.85 or attempt == 2:
            break
        if attempt == 1:
            fs -= 0.4
    return dict(key=key, fs=fs, lw=lw, hw=hw, pw=pw, plot_w=min(max(plot_w, 0.45), 2.2), g=g)


def draw_terms(fig, W: float, H: float, top: float, x0: float, width: float, disp: list[dict], alpha: float,
               pitch: float = 0.19, fs: float = 6.4, tag: str = "", xlim=None, cols: dict | None = None,
               ph_below: float = 0.0, q_below: float = 0.0) -> list[dict]:
    """Draw a model forest into `fig`: term labels from x0, the CI plot, then 'HR (95% CI)' and p, within `width`
    inches; the column heads sit just above `top` (inches from the figure top). `xlim` shares an HR axis between
    several forests. Returns the plotted rows."""
    from matplotlib import ticker as mticker

    fx = lambda x: x / W                                                     # noqa: E731
    fy = lambda y: 1 - y / H                                                 # noqa: E731
    n = len(disp)
    vals = [d["r"] for d in disp if d["r"] is not None]
    c = cols or terms_columns(width, [disp], fs, ph_below, q_below)
    fs, plot_w, g = c["fs"], c["plot_w"], c["g"]
    px = x0 + width - c["pw"] - c["hw"] - 2 * g - plot_w
    hx, pr = px + plot_w + g, x0 + width
    labels = [d[c["key"]] for d in disp]
    xlo, xhi = xlim or axis_range(vals)
    ax = fig.add_axes([fx(px), fy(top + pitch * n), plot_w / W, pitch * n / H])
    ax.set_xscale("log")
    ax.set_xlim(xlo, xhi)
    ax.set_ylim(n - 0.5, -0.5)
    ax.axvline(1, color=S.INK2, lw=0.6, zorder=1)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.xaxis.grid(True, color=S.GRID, lw=0.4)
    ax.set_axisbelow(True)
    ticks = _ticks(xlo, xhi, plot_w, fs - 0.4)
    ax.set_xticks(ticks, [f"{t:g}" for t in ticks])
    ax.xaxis.set_minor_locator(mticker.NullLocator())
    ax.set_xlabel("Hazard ratio (95% CI)", labelpad=1)
    fig.text(fx(hx), fy(top - 0.05), "HR (95% CI)", fontsize=fs - 0.2, color=S.INK2, ha="left", va="bottom")
    fig.text(fx(pr), fy(top - 0.05), "p", fontsize=fs - 0.2, color=S.INK2, ha="right", va="bottom")
    rows = []
    for i, d in enumerate(disp):
        yy = fy(top + pitch * (i + 0.5))
        r = d["r"]
        main = r is not None and r.kind in MAIN
        fig.text(fx(x0 + (0.14 if d["indent"] else 0)), yy, labels[i], fontsize=fs, ha="left", va="center",
                 color=S.INK if r is None or main else S.INK2, fontweight="bold" if main else "normal")
        if r is None:
            continue
        ax.plot([max(r.ci_low, xlo), min(r.ci_high, xhi)], [i, i], color=S.INK, lw=0.8, solid_capstyle="butt")
        for edge, cond, mk in ((xlo, r.ci_low < xlo, "<"), (xhi, r.ci_high > xhi, ">")):
            if cond:
                ax.plot([edge], [i], marker=mk, ms=2.6, color=S.INK, mew=0, clip_on=False)
        ax.scatter([min(max(r.hr, xlo), xhi)], [i], marker="D" if main else "s", s=16 if pitch >= 0.18 else 12,
                   facecolor=S.INK if r.p < alpha else "white", edgecolor=S.INK, lw=0.7, zorder=3)
        fig.text(fx(hx), yy, f"{r.hr:.2f} ({r.ci_low:.2f}–{r.ci_high:.2f})", fontsize=fs - 0.2, ha="left",
                 va="center")
        ph = float(getattr(r, "ph_p", np.nan))
        bad = ph_below > 0 and ph < ph_below                 # the term's hazards are not proportional
        q = d.get("q", np.nan)
        fig.text(fx(pr), yy, p_text(d, ph_below, q_below), fontsize=fs - 0.2, ha="right", va="center")
        rows.append(dict(tag=tag, label=d["label"], term=r.term, kind=r.kind, level=r.level, reference=r.reference,
                         unit=r.unit, coef=r.coef, se=r.se, hr=r.hr, ci_low=r.ci_low, ci_high=r.ci_high, p=r.p,
                         ph_p=ph, ph_marked=bad, q=q, q_marked=bool(q_below > 0 and q is not None and q < q_below)))
    return rows


def cox_model_figure(ds: Dataset, event: str, cohort: str, endpoint: str, *, model: CoxModel | None = None,
                     settings: Settings | None = None, results: Results | None = None, out_dir=None,
                     stem: str | None = None) -> Panel:
    """Draw (and, with out_dir, save) the multivariable Cox forest of one event in one cohort for one endpoint."""
    import matplotlib
    from matplotlib.figure import Figure

    s = settings or Settings()
    if results is not None and model is not None and model != results.model:
        raise InputError(model_mismatch(model, results.model))
    model = model or (results.model if results is not None else CoxModel())
    row, terms = model_terms(ds, event, cohort, endpoint, model, s, results)
    disp = display_rows(terms, float(getattr(row, "cox_q", np.nan)))
    gene, label = ds.events.at[event, "gene"], ds.events.at[event, "label"]
    # the other events of the q family (results computed over the gene) set the printed q: hashed too
    qty = EV.quantity(ds.events.at[event, "event_type"])
    family = [] if results is None else [
        e for e in dict.fromkeys(results.survival.event_id) if e != event and e in ds.psi.index
        and ds.events.at[e, "gene"] == gene and EV.quantity(ds.events.at[e, "event_type"]) == qty]
    hosts = [h for h in dict.fromkeys(ds.events.loc[[event] + family, "expression_gene"])
             if ds.expression is not None and model.expression and h in ds.expression.index]
    W, top, pitch = 5.6, 0.72, 0.19
    H = top + pitch * len(disp) + 0.55
    fx = lambda x: x / W                                                     # noqa: E731
    fy = lambda y: 1 - y / H                                                 # noqa: E731
    with matplotlib.rc_context(S.rc()):
        fig = Figure(figsize=(W, H))
        fig.text(fx(0.12), fy(0.10), gene, fontsize=8.5, fontweight="bold", style="italic", ha="left", va="top")
        tw = 0.12 + S.text_width(gene, 8.5, style="italic", weight="bold") + 0.05
        fig.text(fx(tw), fy(0.10), f"{label} · {cohort} · {endpoint}", fontsize=8.5, fontweight="bold", ha="left",
                 va="top")
        lp = row.get("cox_low_power", False)
        low = " (low power)" if pd.notna(lp) and bool(lp) else ""
        fig.text(fx(0.12), fy(0.33), f"Cox: {row.cox_model} · {int(row.cox_n)} patients, "     # as fitted here
                 f"{int(row.cox_events)} events{low} · {ENDPOINT_NAMES.get(endpoint, endpoint).lower()}",
                 fontsize=6.2, color=S.INK2, ha="left", va="top")
        rows = draw_terms(fig, W, H, top, 0.12, 5.33, disp, s.alpha, pitch, ph_below=s.ph_note_below,
                          q_below=s.q_mark_below)
        fig.text(fx(0.12), fy(H - 0.08), MODEL_NOTE.format(alpha=s.alpha), fontsize=5.6, color=S.MUTED, ha="left",
                 va="bottom")
    head = pd.DataFrame([dict(panel="header", cox_model=row.cox_model, cox_n=int(row.cox_n),       # as printed
                              cox_events=int(row.cox_events), cox_low_power=bool(low), endpoint=endpoint)])
    table = pd.concat([head, pd.DataFrame(rows).assign(panel="term"), terms.assign(panel="cox_terms_all")],
                      ignore_index=True)
    stem = stem or "_".join([safe_name(gene), safe_name(label), safe_name(cohort), safe_name(endpoint), "cox"])
    prov = P.record("cox_model_figure", s, dict(psi=ds.psi.loc[[event]], samples=ds.samples, survival=ds.survival,
                                                 clinical=None if ds.clinical is None else ds.clinical[
                                                     [c for c in model.clinical_columns if c in ds.clinical.columns]],
                                                 expression=ds.expression.loc[hosts] if hosts else None,
                                                 psi_q_family=ds.psi.loc[family] if family else None),
                    dict(event=event, cohort=cohort, endpoint=endpoint, model=model.to_dict()))
    paths = S.save(fig, out_dir, stem, s.formats, s.dpi, table, prov) if out_dir is not None else {}
    return Panel(figure=fig, table=table, paths=paths, provenance=prov, stem=stem)
