"""Probe: every event of a gene (or every event given) in every cohort, ranked, with assay pages for the best-ranked
events (max_pages).

For one endpoint (default OS) and every event x cohort, the probe runs the group tests, the KM split and two Cox
models: the base model (PSI + host expression, the forest's model) and the adjusted model (by default + age + sex +
stage, found in the clinical table; see clinical.py). Then:

    cells.csv       one row per event x cohort: group tests, KM, base and adjusted Cox, q values (BH within each
                    gene, per kind of test)
    events.csv      one row per event, ranked (see RANKING)
    overview.png    events x cohorts: adjusted HR per IQR or SD (colour), adjusted p < 0.05 (dot), q < 0.05 (*),
                    group hit (frame); with an expression table, the host genes' own rows below the events
    gene_map_<GENE>.png   per gene: its model (with a GTF) and its probed events observed in enough samples, 5' to 3',
                    each row beside its cells of the overview (the gene's own row: its expression) and, with
                    correlation=True, the median rho between its events (see plot/genemap.py)
    pages/          one assay page per ranked, measurable event, at most max_pages (the event, its most promising
                    cohorts, their models, the forest); with an expression table, also the host gene's expression
                    page (<GENE>_expression)
    probe.pdf       the overview, the gene maps, the correlation figure (with correlation=True), every event page in
                    rank order, then the expression page
    report.md       what was run, the ranking, how to read it, and how to reproduce it
    proteins.csv    with a protein cache: the suggested protein change of every event (see protein.py), also
                    summarised in events.csv and drawn on the pages
    expression_cells.csv   with an expression table: the host gene's own statistics per cohort (case vs reference,
                    KM split, Cox on expression + the adjusted model's clinical terms), as on its expression page, and
                    one row per host gene in the overview
    agent_prompt.md the instructions and aggregate results an agent needs to write a narrative of the probe (see
                    agent.py; `splice-assay summarize` or `probe --agent-summary` hands it to Claude Code)
    correlations.csv, correlation.png   with correlation=True (--correlation): Spearman rho per cohort of each event
                    with its host gene's expression and of each pair of events of one gene (see correlation.py); rho
                    with expression also in cells.csv (expr_rho) and events.csv (best_expr_rho)

Everything is nominal: a probe ranks candidates, it does not test a hypothesis. The report states how many
p < 0.05 results chance alone would give.
"""
from __future__ import annotations

import argparse
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from . import __version__
from .analysis import Results, analyze
from .clinical import adjusted as clinical_adjusted
from .clinical import auto_clinical, found_summary
from .config import Settings
from .dataset import Dataset, InputError
from .events import opt_in, quantity
from .stats.survival import CoxModel, ridge_text

ALL = 10_000                     # `top` meaning every cohort with a test
MAX_MAPS = 20                    # gene maps: the genes of the overview's best-ranked events, at most this many
RANKING = ("measurable events first; then cohorts where the adjusted Cox p < 0.05, then where the base Cox p < 0.05, "
           "then cohorts with both a group hit and a survival hit, then cohorts where the KM p < 0.05, counting Cox "
           "fits only when they are not low power; then the low-power Cox hits (‡), adjusted, then base; then the p of "
           "the best cohort (adjusted Cox, else base Cox, else KM)")
ADJ_COLS = ["cox_status", "cox_model", "cox_n", "cox_events", "cox_n_dropped", "cox_notes", "hr_per_iqr",
            "ci_low_iqr", "ci_high_iqr", "hr_per_sd", "ci_low_sd", "ci_high_sd", "cox_p", "cox_q", "cox_q_tests",
            "cox_low_power", "cox_events_per_term", "psi_narrow", "ph_p"]


@dataclass
class ProbeResult:
    cells: pd.DataFrame
    events: pd.DataFrame
    paths: dict = field(default_factory=dict)
    pages: list = field(default_factory=list)
    correlations: pd.DataFrame | None = None     # with correlation=True: one row per pair x cohort (correlation.py)


def select_events(ds: Dataset, genes=None, events=None, include_hit: bool = False) -> list[str]:
    """The events to probe: those named; else every event of `genes`, else every event. HIT-index events are left
    out of the last two unless include_hit (a named HIT event is probed)."""
    if events:
        miss = [e for e in events if e not in ds.psi.index]
        if miss:
            raise InputError(f"event(s) not in the data: {', '.join(miss[:5])}")
        return list(dict.fromkeys(events))                # an event named twice is probed once
    if genes:
        want = {str(g).lower() for g in genes}
        ev = ds.events
        found = ev.index[ev.gene.str.lower().isin(want) | ev.gene_id.str.lower().isin(want)]
        if not len(found):
            raise InputError(f"no event of gene(s) {', '.join(genes)} in the data")
        ids = [e for e in ds.psi.index if e in set(found)]
    else:
        ids = list(ds.psi.index)
    out = [e for e in ids if include_hit or not opt_in(ds.events.at[e, "event_type"])]
    if ids and not out:
        raise InputError("only HIT-index events " + (f"for {', '.join(map(str, genes))}" if genes else "in the data")
                         + ": add --include-hit (include_hit=True) to probe them")
    return out


def _describe(model: CoxModel, ds: Dataset, events, s: Settings | None = None) -> str:
    """The model's terms as fitted for these events, with a ridge penalty that reaches them (describe_model)."""
    from .analysis import describe_model
    return describe_model(model, ds, events, s)


def default_endpoint(ds: Dataset, endpoint=None) -> str:
    if endpoint:
        return endpoint
    if ds.survival is None or not ds.endpoints:
        raise InputError("no survival table: a probe needs survival data")
    return "OS" if "OS" in ds.endpoints else ds.endpoints[0]


def combine(base: Results, adj: Results | None, endpoint: str) -> pd.DataFrame:
    """One row per event x cohort for one endpoint: group tests, KM, base Cox and adjusted Cox (adj_*)."""
    sv = base.survival[base.survival.endpoint.eq(endpoint)]
    cells = sv.merge(base.groups, on=["event_id", "gene", "cohort"], how="left")
    if adj is not None:
        a = adj.survival[adj.survival.endpoint.eq(endpoint)]
        a = a[["event_id", "cohort"] + [c for c in ADJ_COLS if c in a.columns]]
        cells = cells.merge(a.rename(columns={c: f"adj_{c}" for c in ADJ_COLS}), on=["event_id", "cohort"], how="left")
    return cells                                         # q values come from the results: BH within each gene


def _flag(v) -> bool:
    """A True/False cell that may be missing (NaN, None) or a numpy bool."""
    return bool(v) if v is not None and not (isinstance(v, float) and np.isnan(v)) else False


def _true(col: pd.Series) -> pd.Series:
    """Where a True/False column is True, a missing value counting as False (without pandas' fillna downcasting of an
    object column, deprecated)."""
    return col.map(lambda v: False if pd.isna(v) else bool(v)).astype(bool)


def _low(d: pd.DataFrame, col: str) -> pd.Series:
    """A low-power flag column (fewer events than Settings.cox_low_power_events), False where absent or missing."""
    return _true(d[col]) if col in d else pd.Series(False, index=d.index)


def _tier(p: pd.Series, low: pd.Series, alpha: float) -> np.ndarray:
    """The order of an event's cohorts: p < alpha in a fit that is not low power (0), p < alpha in a low-power fit
    (1), then the other fits (2) and the other low-power fits (3); a missing p counts as not below alpha."""
    with np.errstate(invalid="ignore"):
        hit = (p < alpha).to_numpy()
    lo = low.to_numpy(bool)
    return np.where(hit, np.where(lo, 1, 0), np.where(lo, 3, 2))


def observed(cells: pd.DataFrame, s: Settings) -> set:
    """The events observed in at least s.observed_frac of a cohort's survival samples, in some cohort: those the gene
    maps and correlations take."""
    f = cells.groupby("event_id").frac_obs.max()
    return set(f.index[f >= s.observed_frac])


def rank_events(cells: pd.DataFrame, ds: Dataset, alpha: float, unit: str = "sd") -> pd.DataFrame:
    """One row per event, ranked (RANKING): Cox hits in fits that are not low power first, low-power hits (fewer
    events than Settings.cox_low_power_events) as tie-breakers. `unit`: the HR reported, per "sd" or per "iqr"
    (Settings.psi_hr_unit)."""
    hr = f"hr_per_{unit}"
    rows = []
    for eid, d in cells.groupby("event_id", sort=False):
        t = d[d.cox_status.eq("tested")]
        sig = t[t.cox_p < alpha]
        has_adj = "adj_cox_p" in d
        ta = d[d.get("adj_cox_status", pd.Series(index=d.index, dtype=object)).eq("tested")] if has_adj else d.iloc[:0]
        sig_a = ta[ta.adj_cox_p < alpha] if has_adj else ta
        with np.errstate(invalid="ignore"):            # a survival hit: KM, or base Cox in a fit not low power
            surv_hit = (d.km_status.eq("tested") & (d.km_p < alpha)) | (
                d.cox_status.eq("tested") & (d.cox_p < alpha) & ~_low(d, "cox_low_power"))
        tk = d[d.km_status.eq("tested")]
        if len(ta):                                     # the best cohort: a hit not low power, then a low-power hit
            best = ta.assign(_t=_tier(ta.adj_cox_p, _low(ta, "adj_cox_low_power"), alpha)) \
                .sort_values(["_t", "adj_cox_p"], kind="stable")
            pcol, hcol, model = "adj_cox_p", f"adj_{hr}", "adjusted"
        elif len(t):
            best = t.assign(_t=_tier(t.cox_p, _low(t, "cox_low_power"), alpha)).sort_values(["_t", "cox_p"],
                                                                                         kind="stable")
            pcol, hcol, model = "cox_p", hr, "base"
        else:                                           # no Cox model anywhere (e.g. too few deaths): the KM test
            best, pcol, hcol, model = tk.sort_values("km_p"), "km_p", None, "KM"
        b = best.iloc[0] if len(best) else None
        rows.append(dict(
            event_id=eid, gene=ds.events.at[eid, "gene"], label=ds.events.at[eid, "label"],
            event_type=ds.events.at[eid, "event_type"], cohorts_cox_tested=len(t),
            cox_p05=len(sig), cox_p05_hr_up=int((sig[hr] > 1).sum()), cox_p05_hr_down=int((sig[hr] < 1).sum()),
            expected_by_chance=round(alpha * len(t), 2), adj_cohorts_tested=len(ta), adj_cox_p05=len(sig_a),
            adj_cox_p05_low_power=int(_low(sig_a, "adj_cox_low_power").sum()),
            cox_p05_low_power=int(_low(sig, "cox_low_power").sum()),
            group_hits=int(_true(d.group_hit).sum()),
            within_patient=int(_true(d.within_patient_support).sum()),
            group_and_survival=int((_true(d.group_hit) & surv_hit).sum()),
            share_hr_up=round(float((t[hr] > 1).mean()), 2) if len(t) else np.nan,
            km_tested=len(tk), km_p05=int((tk.km_p < alpha).sum()),
            best_cohort=None if b is None else b.cohort,
            **{f"best_{hr}": np.nan if b is None or hcol is None else b[hcol]},
            best_logrank_hr=np.nan if b is None or model != "KM" else b.logrank_hr,
            best_p=np.nan if b is None else b[pcol], best_model=model if b is not None else "",
            best_low_power=_flag(None if b is None else
                                 b.get({"adjusted": "adj_cox_low_power", "base": "cox_low_power"}.get(model, ""))),
            **{f"best_psi_{unit}": np.nan if b is None or model == "KM" else b[f"psi_{unit}"]},
            min_cox_q=float(t.cox_q.min()) if len(t) else np.nan))
    ev = pd.DataFrame(rows)
    if ev.empty:
        return ev
    ev["measurable"] = (ev.cohorts_cox_tested > 0) | (ev.km_tested > 0)
    ev = ev.assign(_adj=ev.adj_cox_p05 - ev.adj_cox_p05_low_power, _base=ev.cox_p05 - ev.cox_p05_low_power)
    ev = ev.sort_values(["measurable", "_adj", "_base", "group_and_survival", "km_p05", "adj_cox_p05_low_power",
                         "cox_p05_low_power", "best_p"], ascending=[False] * 7 + [True], na_position="last",
                        kind="stable").drop(columns=["_adj", "_base"]).reset_index(drop=True)
    ev.insert(0, "rank", np.arange(1, len(ev) + 1))
    return ev


def pick_cohorts(cells: pd.DataFrame, event: str, n: int = 3, alpha: float = 0.05) -> list[str]:
    """The event's most promising cohorts: Cox p < alpha in a fit that is not low power, then in a low-power fit, then
    the rest (the adjusted model where fitted, else the base model); within each, the smallest adjusted Cox p, then
    base Cox p, then KM p."""
    d = cells[cells.event_id.eq(event)].copy()
    d["_a"] = d.get("adj_cox_p", pd.Series(np.nan, index=d.index)).where(
        d.get("adj_cox_status", pd.Series("", index=d.index)).eq("tested"))
    d["_b"] = d.cox_p.where(d.cox_status.eq("tested")) if "cox_p" in d else np.nan
    d["_k"] = d.km_p.where(d.km_status.eq("tested")) if "km_p" in d else np.nan
    d = d[d[["_a", "_b", "_k"]].notna().any(axis=1)]
    has_a = d._a.notna()
    d["_t"] = _tier(d._a.where(has_a, d._b), _low(d, "adj_cox_low_power").where(has_a, _low(d, "cox_low_power")),
                    alpha)
    d = d.sort_values(["_t", "_a", "_b", "_k"], na_position="last", kind="stable")
    return d.cohort.head(n).tolist()


def parse_top(v) -> int:
    """--top: a number of cohorts, or 'all' (every cohort with a test)."""
    if str(v).strip().lower() == "all":
        return ALL
    try:
        n = int(v)
    except ValueError:
        raise argparse.ArgumentTypeError(f"--top expects a number or 'all', not {v!r}") from None
    if n < 1:
        raise argparse.ArgumentTypeError("--top must be at least 1")
    return n


def _gex_rows(gex_cells: pd.DataFrame | None, ev: pd.DataFrame, endpoint: str,
              hosts=None) -> list[tuple[str, pd.DataFrame]]:
    """The overview's expression rows: (gene, its cells by cohort) for the host genes of the events shown (`hosts`:
    event ID -> host gene, default the event's gene) that have expression statistics, in rank order; every gene of
    gex_cells when no event is shown."""
    if gex_cells is None or not len(gex_cells):
        return []
    g = gex_cells[gex_cells.endpoint.eq(endpoint)] if "endpoint" in gex_cells else gex_cells
    have = list(dict.fromkeys(g.gene))
    host = (lambda r: hosts.get(r.event_id, r.gene)) if hosts is not None else (lambda r: r.gene)
    shown = [h for h in dict.fromkeys(host(r) for r in ev.itertuples()) if h in set(have)]
    return [(h, g[g.gene.eq(h)].set_index("cohort")) for h in (shown if len(ev) else have)]


def _fit_note(cells: pd.DataFrame, ev: pd.DataFrame, s: Settings) -> str:
    """Which fit the survival cells of these events show, for a caption: "adjusted model where fitted", with a ridge
    penalty named when the fits shown have one."""
    shown = cells[cells.event_id.isin(ev.event_id)]                       # the fit each cell shows: adjusted, else base
    adj = shown["adj_cox_status"].eq("tested") if "adj_cox_model" in shown else pd.Series(False, index=shown.index)
    names = pd.concat([shown.loc[adj, "adj_cox_model"] if adj.any() else pd.Series(dtype=str),
                       shown.loc[~adj & shown.cox_status.eq("tested"), "cox_model"]]).astype(str)
    pen = names[names.str.contains("; ridge λ", regex=False)]
    ridge = "" if pen.empty else f"; ridge λ {s.cox_ridge_penalty:g} on " + {
        "clinical": "clinical terms", "all": "all terms",
        "molecular": _value_name(ev.event_type) + (" and host expression" if pen.str.contains(
            "host expression", regex=False).any() else "")
    }[s.cox_ridge]
    return f"adjusted model where fitted{ridge}"


def overview(cells: pd.DataFrame, events: pd.DataFrame, endpoint: str, s: Settings, max_rows: int = 60,
             gex_cells: pd.DataFrame | None = None, gex_model: str = "", hosts=None):
    """Events (rows, rank order) x cohorts: HR per SD or IQR (adjusted when fitted, else base) in colour, p < alpha as a
    dot, a group hit as a frame; grey = not tested.

    With `gex_cells` (the host genes' own statistics, as in expression_cells.csv), one row per host gene below the
    events (`hosts`: event ID -> host gene, e.g. Dataset.events.expression_gene): the HR per SD of expression from Cox
    on `gex_model` (e.g. "expression + age + sex + stage"), p < alpha as a dot (no q mark: expression is not adjusted
    for multiple testing), an expression group hit as a frame."""
    import matplotlib
    from matplotlib.figure import Figure
    from matplotlib.patches import Rectangle

    from .plot import grid as G
    from .plot import style as S

    ev = events[events.measurable].head(max_rows)
    gx = _gex_rows(gex_cells, ev, endpoint, hosts)
    cohorts = sorted(cells.cohort.unique())
    nr, nc = len(ev), len(cohorts)
    gap = 0.5 if nr and gx else 0.0                                       # between the events and the gene rows
    n_rows = max(nr + gap + len(gx), 1)
    cell_w, cell_h = min(0.2, 5.6 / max(nc, 1)), 0.17
    label_w = max([S.text_width(f"{r.label}  {r.gene}", 6.0) for r in ev.itertuples()]
                  + [S.text_width(f"{g} expression", 6.0) for g, _ in gx] + [0.8]) + 0.25
    v = _value_name(ev.event_type)
    caption = (f"colour: HR per {s.psi_hr_unit.upper()} of {v} ({_fit_note(cells, ev, s)})"
               + f"; dot: p < {s.alpha:g} (large: < 0.01)"
               + (f"; {S.Q_MARK}: q < {s.q_mark_below:g}" if s.q_mark_below > 0 else "")
               + "; frame: group hit; grey: not tested")
    gex_caption = (f"{'last row' if len(gx) == 1 else 'last rows'}: host-gene expression, HR per SD of expression "
                   f"(Cox on {gex_model or 'expression'}); frame: group hit (|Δ median| > {s.gex_min_abs_delta:g}); "
                   "not adjusted for multiple testing") if gx else ""
    W = max(4.5, label_w + nc * cell_w + 0.4)
    lines = [G.clauses(c, W - 0.24, 5.8) for c in (caption, gex_caption) if c]       # each caption on its own lines
    n_lines = sum(len(x) for x in lines)
    lab_h = max(S.text_width(c, 5.6) for c in cohorts) + 0.06 if cohorts else 0.2   # the cohort names, upright
    top = 0.30 + 0.12 * n_lines + 0.10 + lab_h
    H = top + n_rows * cell_h + 0.62
    scale = G.hr_scale()
    cmap, norm = scale
    idx = cells.set_index(["event_id", "cohort"])
    with matplotlib.rc_context(S.rc()):
        fig = Figure(figsize=(W, H))
        ax = fig.add_axes([label_w / W, 0.62 / H, nc * cell_w / W, n_rows * cell_h / H])
        ax.set_xlim(0, nc)
        ax.set_ylim(n_rows, 0)
        for k, (gene, g) in enumerate(gx):                                  # the host genes' own rows
            i = nr + gap + k
            for j, c in enumerate(cohorts):
                if c not in g.index:
                    continue
                x = g.loc[c]
                tested = x.get("cox_status") == "tested"
                hr, p = x.get("hr_per_sd", np.nan), x.get("cox_p", np.nan)
                face = cmap(norm(np.log2(hr))) if tested and np.isfinite(hr) else G.UNTESTED
                ax.add_patch(Rectangle((j + 0.04, i + 0.06), 0.92, 0.88, facecolor=face, lw=0))
                if _flag(x.get("group_hit", False)):
                    ax.add_patch(Rectangle((j + 0.04, i + 0.06), 0.92, 0.88, fill=False, edgecolor=S.INK, lw=0.7))
                if tested and np.isfinite(p) and p < s.alpha:
                    ax.plot([j + 0.5], [i + 0.5], "o", ms=2.6 if p >= 0.01 else 3.6, color=G.ink_on(face), mew=0)
        for i, r in enumerate(ev.itertuples()):
            for j, c in enumerate(cohorts):
                if (r.event_id, c) not in idx.index:
                    continue
                x = idx.loc[(r.event_id, c)]
                tested, hr, p, q = G.shown_fit(x, s)
                face = cmap(norm(np.log2(hr))) if tested and np.isfinite(hr) else G.UNTESTED
                ax.add_patch(Rectangle((j + 0.04, i + 0.06), 0.92, 0.88, facecolor=face, lw=0))
                if _flag(x.get("group_hit", False)):
                    ax.add_patch(Rectangle((j + 0.04, i + 0.06), 0.92, 0.88, fill=False, edgecolor=S.INK, lw=0.7))
                if tested and s.q_mark_below > 0 and np.isfinite(q) and q < s.q_mark_below:   # q, not just p
                    ax.text(j + 0.5, i + 0.62, S.Q_MARK, fontsize=9, color=G.ink_on(face), ha="center",
                            va="center", fontweight="bold")
                elif tested and np.isfinite(p) and p < s.alpha:
                    ax.plot([j + 0.5], [i + 0.5], "o", ms=2.6 if p >= 0.01 else 3.6, color=G.ink_on(face), mew=0)
        ax.set_xticks(np.arange(nc) + 0.5, cohorts, rotation=90, fontsize=5.6)
        ax.set_yticks(list(np.arange(nr) + 0.5) + [nr + gap + k + 0.5 for k in range(len(gx))],
                      [f"{r.rank}. {r.label}  {r.gene}" for r in ev.itertuples()] + [f"{g} expression" for g, _ in gx],
                      fontsize=6.0)
        for lab, r in zip(ax.get_yticklabels(), ev.itertuples()):         # each event in its type's colour
            lab.set_color(S.TYPE_COLOR.get(str(r.event_type).upper(), S.TYPE_DEFAULT))
        ax.tick_params(length=0)
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.xaxis.tick_top()
        fig.text(0.12 / W, 1 - 0.1 / H, f"Probe overview · {endpoint}", fontsize=8.5, fontweight="bold", va="top")
        for k, line in enumerate(x for block in lines for x in block):
            fig.text(0.12 / W, 1 - (0.30 + 0.12 * k) / H, line, fontsize=5.8, color=S.INK2, va="top")
        G.colorbar(fig, label_w, H - 0.33, min(1.6, max(nc * cell_w, 1.0)), W, H, scale, [-1, 0, 1],
                   ["0.5", "1", "2"], f"HR per {s.psi_hr_unit.upper()}"
                   + (" (expression: per SD)" if gx and s.psi_hr_unit != "sd" else ""))
    return fig


def gene_events(events_table: pd.DataFrame, genes) -> list[str]:
    """The event IDs of some genes (by name or gene ID, case-insensitive) in an events table."""
    want = {str(g).lower() for g in genes}
    e = events_table
    gid = e["gene_id"].astype(str).str.lower() if "gene_id" in e.columns else pd.Series("", index=e.index)
    return e.loc[e["gene"].astype(str).str.lower().isin(want) | gid.isin(want), "event_id"].astype(str).tolist()


def probe(ds: Dataset, genes=None, events=None, cohorts=None, endpoint=None, *, settings: Settings | None = None,
          model: CoxModel | None = None, adjusted="auto", baseline: dict | None = None, gtf=None, out_dir=None,
          top: int = 3, max_pages: int = 30, log=print, call: str = "", proteins=None, gex: bool = True,
          include_hit: bool = False, correlation: bool = False) -> ProbeResult:
    """Probe events (all of `genes`, the `events` given, or every event) in `cohorts` (default all) for one endpoint
    (default OS). `adjusted`: "auto" (age, sex and stage found in the clinical table), a CoxModel, or None.
    `proteins`: a protein cache (protein.ProteinCache or its folder) for the suggested protein changes. `gex`: the
    host gene's own expression statistics (expression_cells.csv, its own page and its row in the overview) when
    expression is given. `include_hit`: also probe HIT-index events (left out by default: one per exon, a far larger
    set). `correlation`: Spearman rho of each event with its host gene's expression and with the other events of its
    gene, per cohort (correlations.csv, correlation.png; see correlation.py)."""
    from matplotlib.backends.backend_pdf import PdfPages

    from .annotation import read_gtf
    from .events import geometry
    from .plot import event_panel, expression_panel
    from .plot import style as S
    from .plot.panel_common import page_parts, safe_name

    s = settings or Settings()
    ep = default_endpoint(ds, endpoint)
    ev = select_events(ds, genes, events, include_hit)
    hit_left = 0 if events or include_hit else len(select_events(ds, genes, None, True)) - len(ev)
    cohorts = None if not cohorts or [str(c).lower() for c in cohorts] == ["all"] else list(cohorts)
    base_model = model or CoxModel()
    found = {}
    if adjusted == "auto":                              # the base model plus the age, sex and stage found
        ds, auto, found = auto_clinical(ds, keep=base_model.clinical_columns)
        adj_model = None if auto is None else clinical_adjusted(base_model, auto, found)
        if adj_model is not None and baseline:
            adj_model = adj_model.with_clinical(baseline=baseline)
    else:
        adj_model = adjusted or None
    if s.changed():
        log(f"settings changed from the defaults: {s.changed_text()}")
    if hit_left:
        log(f"{hit_left} HIT-index event(s) left out (--include-hit adds them)")
    log(f"probe: {len(ev)} event(s) x {len(cohorts or ds.cohorts)} cohort(s), {ep}; base Cox "
        f"{_describe(base_model, ds, ev, s)}"
        + (f"; adjusted {_describe(adj_model, ds, ev, s)}" if adj_model else "; no adjusted model"))
    base = analyze(ds, events=ev, endpoints=[ep], cohorts=cohorts, settings=s, model=base_model)
    adj = analyze(ds, events=ev, endpoints=[ep], cohorts=cohorts, settings=s, model=adj_model) if adj_model else None
    cells = combine(base, adj, ep)
    keep = observed(cells, s)                          # present enough for the gene maps and correlations
    rare = [e for e in ev if e not in keep]
    if rare:
        log(f"{len(rare)} event(s) observed in under {s.observed_frac:.0%} of every cohort's survival samples: left "
            "out of the gene maps" + (" and correlations" if correlation else "") + " (--min-observed)")
    corr = None
    if correlation:
        from .correlation import CORR_FIRST, correlations
        kept = [e for e in ev if e in keep]
        corr = correlations(ds, events=kept, cohorts=cohorts, settings=s) if kept else pd.DataFrame(
            columns=CORR_FIRST)
        x = corr[corr.kind.eq("expression")]
        if len(x):                                     # rho with host expression beside each cell's tests
            cells = cells.merge(x[["event_id", "cohort", "corr_n", "rho", "corr_p", "corr_q"]].rename(columns=dict(
                corr_n="expr_rho_n", rho="expr_rho", corr_p="expr_rho_p", corr_q="expr_rho_q")),
                on=["event_id", "cohort"], how="left")
        elif ds.expression is not None:                # every event left out: the columns, empty
            cells = cells.assign(expr_rho_n=np.nan, expr_rho=np.nan, expr_rho_p=np.nan, expr_rho_q=np.nan)
        n_t = corr.corr_status.eq("tested")
        log(f"correlation: {int((n_t & corr.kind.eq('expression')).sum())} event x cohort cells with host expression, "
            f"{int((n_t & corr.kind.eq('event')).sum())} event pair x cohort cells tested (Spearman, at least "
            f"{s.corr_min_n} patients)")
    ranked = rank_events(cells, ds, s.alpha, s.psi_hr_unit)
    if "expr_rho" in cells and len(ranked):            # rho with expression in each event's best cohort
        rho = cells.set_index(["event_id", "cohort"]).expr_rho
        ranked["best_expr_rho"] = [rho.get((e, c), np.nan) for e, c in zip(ranked.event_id, ranked.best_cohort)]
    ptab, pchanges = None, {}
    if proteins is not None and len(ranked):
        from .protein import SHOWN, ProteinCache, protein_changes
        cache = proteins if isinstance(proteins, ProteinCache) else ProteinCache.load(proteins)
        ptab, pchanges = protein_changes(cache, ds, ranked.event_id.tolist())
        ranked = ranked.merge(ptab.reindex(columns=["event_id", "protein_status", "effect_short", "summary"])
                              .rename(columns={"effect_short": "protein_change", "summary": "protein_summary"}),
                              on="event_id", how="left")
        log(f"proteins: a suggestion for {int(ptab.protein_status.isin(SHOWN).sum())} of {len(ptab)} event(s) "
            f"({cache.source})")
    gex_cells, gex_model = None, ""
    genes_x = [g for g in dict.fromkeys(ds.events.loc[ev, "expression_gene"])
               if ds.expression is not None and g in ds.expression.index]
    if gex and genes_x:
        from .analysis import analyze_expression
        gm = CoxModel(expression=False, covariates=(adj_model or base_model).covariates,
                      categorical=(adj_model or base_model).categorical, strata=(adj_model or base_model).strata,
                      baseline=dict((adj_model or base_model).baseline))
        gex_cells = analyze_expression(ds, genes=genes_x, endpoints=[ep], cohorts=cohorts, settings=s,
                                       model=gm).cells()
        gex_model = gm.describe(False).replace("PSI", "expression", 1) + ridge_text(s, bool(gm.covariates),
                                                                                    ["expression"])
    out = ProbeResult(cells, ranked, correlations=corr)
    if out_dir is None:
        return out
    out_dir = Path(out_dir)
    (out_dir / "pages").mkdir(parents=True, exist_ok=True)
    cells.to_csv(out_dir / "cells.csv", index=False, lineterminator="\n")
    ranked_path = out_dir / "events.csv"
    gtf_table = None
    todo = ranked[ranked.measurable].head(max_pages) if len(ranked) else ranked
    drawable = [e for e in todo.event_id if ds.events.at[e, "variable"] and ds.events.at[e, "strand"] in ("+", "-")]
    from .plot.genemap import map_events
    best = ranked[ranked.measurable].head(60) if len(ranked) else ranked          # the overview's rows
    mapped = {g: map_events(ds, g, ranked, keep=keep)[0] for g in dict.fromkeys(best.gene)}
    mapped = dict([(g, ids) for g, ids in mapped.items() if ids][:MAX_MAPS])   # each gene's map: its events, 5' to 3'
    if gtf is not None and (drawable or mapped):
        wins, genes_ = [], set()
        for e in drawable:
            row = ds.events.loc[e]
            lo, hi = geometry(e, row).span
            wins.append((row.chrom, lo - s.gtf_flank, hi + s.gtf_flank))
            genes_ |= {g for g in (row.gene, row.gene_id) if g}
        for ids in mapped.values():                     # one window over each map's events
            spans = [geometry(e, ds.events.loc[e]).span for e in ids]
            row = ds.events.loc[ids[0]]
            wins.append((row.chrom, min(a for a, _ in spans) - s.gtf_flank, max(b for _, b in spans) + s.gtf_flank))
            genes_ |= {g for g in (row.gene, row.gene_id) if g}
        gtf_table = read_gtf(gtf, wins, genes=sorted(genes_))
    pages, page_col, expr_pages, map_pages = [], {}, [], []
    focus_of = {}                                       # the cohorts each event's page shows
    for r in todo.itertuples():
        if r.event_id not in drawable:
            continue
        focus = [c for c in (cohorts or []) if c in set(cells.cohort)][:top] if cohorts and len(cohorts) <= top \
            else pick_cohorts(cells, r.event_id, top, s.alpha)
        if focus:
            focus_of[r.event_id] = focus
    import matplotlib
    with matplotlib.rc_context(S.rc()), PdfPages(out_dir / "probe.pdf", metadata={"CreationDate": None,
                                                                                  "ModDate": None}) as pdf:
        fig = overview(cells, ranked, ep, s, gex_cells=gex_cells, gex_model=gex_model,
                       hosts=ds.events.expression_gene.to_dict())
        S.save(fig, out_dir, "overview", s.formats, s.dpi)
        fig.savefig(pdf, format="pdf")
        from .plot.genemap import MapNotDrawn, gene_map
        for gene, ids in mapped.items():                 # each gene with its events, 5' to 3'
            try:                                         # a summary figure: its failure does not stop the probe
                gm_fig = gene_map(ds, gene, ranked, cells, ep, s, gtf=gtf_table, gex_cells=gex_cells,
                                  gex_model=gex_model, corr=corr, keep=keep,
                                  fit_note=_fit_note(cells, ranked[ranked.event_id.isin(ids)], s))
            except Exception as e:                       # noqa: BLE001
                warnings.warn(f"gene map of {gene} not drawn: {type(e).__name__}: {e}", MapNotDrawn, stacklevel=2)
                log(f"  gene map: {gene} not drawn ({type(e).__name__}: {e})")
                continue
            if gm_fig is None:
                continue
            stem = f"gene_map_{safe_name(gene)}"
            S.save(gm_fig, out_dir, stem, s.formats, s.dpi)
            gm_fig.savefig(pdf, format="pdf")
            map_pages.append(out_dir / f"{stem}.png")
            log(f"  gene map: {gene} ({len(ids)} event{'s' if len(ids) != 1 else ''})")
        cf = None
        if corr is not None:
            from .correlation import correlation_figure
            cf = correlation_figure(corr, ranked, s)
            if cf is not None:
                S.save(cf, out_dir, "correlation", s.formats, s.dpi)
                cf.savefig(pdf, format="pdf")
        corr_drawn = cf is not None
        for r in todo.itertuples():
            focus = focus_of.get(r.event_id)
            if not focus:
                continue
            parts = page_parts(focus, s.cohorts_per_page)  # many cohorts: several pages, each with the full forest
            for i, chunk in enumerate(parts, 1):
                stem = f"{r.rank:03d}_{safe_name(r.event_id)}" + (f"_p{i}" if len(parts) > 1 else "")
                p = event_panel(ds, r.event_id, chunk, ep, settings=s, model=base_model, gtf=gtf_table,
                                results=_slice(base, r.event_id), detail=adj_model, out_dir=out_dir / "pages",
                                stem=stem, proteins={r.event_id: pchanges[r.event_id]} if r.event_id in pchanges
                                else None, detail_results=_slice(adj, r.event_id) if adj is not None else None,
                                part=(i, len(parts)) if len(parts) > 1 else None)
                p.figure.savefig(pdf, format="pdf")
                pages.append((r.event_id, stem, chunk))
                page_col.setdefault(r.event_id, f"pages/{stem}.png")
                log(f"  page {len(pages)}: {r.event_id} ({', '.join(chunk)})")
        if gex and ds.expression is not None:           # each gene's expression last, once, on its own page
            for gene in dict.fromkeys(ds.events.at[e, "expression_gene"] for e in focus_of):
                coh = list(dict.fromkeys(c for e, f in focus_of.items()
                                         if ds.events.at[e, "expression_gene"] == gene for c in f))
                parts = page_parts(coh, s.cohorts_per_page)
                for i, chunk in enumerate(parts, 1):     # unnumbered, so it also sorts after the event pages
                    stem = f"{safe_name(gene)}_expression" + (f"_p{i}" if len(parts) > 1 else "")
                    p = expression_panel(ds, gene, chunk, ep, settings=s, model=adj_model or base_model,
                                         out_dir=out_dir / "pages", stem=stem,
                                         part=(i, len(parts)) if len(parts) > 1 else None)
                    if p is None:
                        break
                    p.figure.savefig(pdf, format="pdf")
                    expr_pages.append(out_dir / "pages" / f"{stem}.png")
                    log(f"  expression page: {gene} ({', '.join(chunk)})")
    ranked["page"] = ranked.event_id.map(page_col).fillna("")
    ranked.to_csv(ranked_path, index=False, lineterminator="\n")
    report = _report(ds, cells, ranked, ep, base_model, adj_model, found, s, top, max_pages, call, ptab, gex_cells,
                     hit_left, corr, [m.stem for m in map_pages], len(rare))
    (out_dir / "report.md").write_text(report)
    from .agent import write_prompt                     # what an agent needs to summarise the probe (agent.py)
    agent_prompt = write_prompt(out_dir, report, ranked, cells, s, ds.events.expression_gene.to_dict(), gex_cells,
                                expression=bool(genes_x and base_model.expression
                                                and (adj_model is None or adj_model.expression)))
    out.events, out.pages = ranked, pages
    out.paths = dict(report=out_dir / "report.md", events=ranked_path, cells=out_dir / "cells.csv",
                     overview=out_dir / "overview.png", pdf=out_dir / "probe.pdf", pages=out_dir / "pages",
                     agent_prompt=agent_prompt)
    if expr_pages:
        out.paths["expression_pages"] = expr_pages[0] if len(expr_pages) == 1 else expr_pages
    if map_pages:
        out.paths["gene_maps"] = map_pages[0] if len(map_pages) == 1 else map_pages
    if ptab is not None:
        ptab.to_csv(out_dir / "proteins.csv", index=False, lineterminator="\n")
        out.paths["proteins"] = out_dir / "proteins.csv"
    if gex_cells is not None:
        gex_cells.to_csv(out_dir / "expression_cells.csv", index=False, lineterminator="\n")
        out.paths["expression"] = out_dir / "expression_cells.csv"
    if corr is not None:
        corr.to_csv(out_dir / "correlations.csv", index=False, lineterminator="\n")
        out.paths["correlations"] = out_dir / "correlations.csv"
        if corr_drawn:
            out.paths["correlation"] = out_dir / "correlation.png"
    return out


def _slice(res: Results, event: str) -> Results:
    """One event's rows of gene-wide results (its q values keep their gene-wide families)."""
    return Results(res.groups[res.groups.event_id.eq(event)].reset_index(drop=True),
                   res.survival[res.survival.event_id.eq(event)].reset_index(drop=True),
                   res.cox_terms[res.cox_terms.event_id.eq(event)].reset_index(drop=True) if len(res.cox_terms)
                   else res.cox_terms, settings=res.settings, model=res.model)


KM_RULES = {"km_split", "km_split_expression"}      # settings that move the KM split rather than relax a gate
GATES = {"min_pairs", "min_group", "coverage_frac", "min_off_modal", "km_min_group", "km_min_events", "cox_min_n",
         "cox_min_events", "cox_min_events_per_term", "fdr_min_family", "low_psi_variance_sd"}   # minimums of a test


def _ridge_line(cells: pd.DataFrame, s: Settings, v: str, base_clinical: bool) -> str | None:
    """The report's line on a ridge penalty, for the models it reached (from their fitted names), in the scope's
    terms; None when no fitted model was penalized. `base_clinical`: whether the base model has clinical terms."""
    models = (("base", "cox_status", "cox_model"), ("adjusted", "adj_cox_status", "adj_cox_model"))
    fits = {w: cells.loc[cells[st].eq("tested") & cells[mo].astype(str).str.contains("; ridge λ", regex=False), mo]
            .astype(str) if st in cells and mo in cells else pd.Series(dtype=str) for w, st, mo in models}
    which = [w for w, names in fits.items() if len(names)]
    if not which:
        return None
    expr = any(names.str.contains("host expression", regex=False).any() for names in fits.values())
    clin = s.cox_ridge == "all" and ("adjusted" in which or base_clinical)   # clinical terms among the penalized
    one = not (expr or clin or s.cox_ridge == "clinical")                  # one penalized term: the PSI term
    terms = "the clinical terms" if s.cox_ridge == "clinical" else f"the {v} term" if one else \
        f"the {v} and host-expression terms" if s.cox_ridge == "molecular" else "all terms"
    partial = (f"adjust {v} only partly: its HR stays closer to the HR without them, and an association that full "
               "adjustment would weaken can survive it")
    if s.cox_ridge == "clinical":
        what = f"they are shrunk jointly toward HR 1, so they {partial}"
    elif one:
        what = "its HR is shrunk toward 1"
    else:                     # jointly: one HR can move away from 1 when terms are correlated
        what = (f"they are shrunk jointly toward HR 1, though a single HR ({v}'s too) can move away from 1 when "
                "terms are correlated" + (f"; the clinical terms {partial}" if clin else ""))
    return (f"- **Penalty:** ridge λ {s.cox_ridge_penalty:g} on {terms} of the {' and '.join(which)} "
            f"model{'s' if len(which) > 1 else ''}: {what}. The CIs and p values of penalized terms are approximate.")


def _relaxed(s: Settings, correlation: bool = False) -> bool:
    """Whether a gate (a minimum a test needs) is below its default; corr_min_n counts when correlation ran."""
    gates = GATES | ({"corr_min_n"} if correlation else set())
    return any(k in gates and v < d for k, (v, d) in s.changed().items())


def _value_name(types) -> str:
    """What the events' values are called in text that covers them all: PSI, HIT index, or PSI or HIT index."""
    q = {quantity(t) for t in types}
    return " or ".join(n for n in ("PSI", "HIT index") if n in q) or "PSI"


def _have(n: int) -> str:
    return "has" if n == 1 else "have"


def _flag_lines(cells: pd.DataFrame, ds: Dataset, s: Settings, adj_model) -> list[str]:
    """The report's notes on the survival tests: narrow PSI ranges, few events per term and non-proportional hazards
    (notes only; every test ran)."""
    out = []
    t = cells[cells.cox_status.eq("tested")]
    ta = cells[cells.get("adj_cox_status", pd.Series(dtype=object, index=cells.index)).eq("tested")] if adj_model \
        else cells.iloc[:0]
    tk = cells[cells.km_status.eq("tested")]
    name = (lambda r: f"{ds.events.at[r.event_id, 'label']} {r.cohort}")
    v = _value_name(ds.events.loc[cells.event_id.unique(), "event_type"])
    if len(t) and s.narrow_psi_below > 0:
        nar = t[_true(t.psi_narrow)]
        sig = nar[nar.cox_p < s.alpha]
        out.append(f"- **Narrow {v} range** ({s.narrow_measure.upper()} of {v} below {s.narrow_psi_below:g} in "
                   f"the fit cohort, so the HR covers a few {v} points): {len(nar)} of {len(t)} base-model fits"
                   + (f"; with p < {s.alpha:g}: {', '.join(name(r) for r in sig.itertuples())}" if len(sig) else "")
                   + ".")
    if len(t) and s.cox_events_per_term > 0:
        few = (lambda d, col: d[d[col] < s.cox_events_per_term] if col in d else d.iloc[:0])
        fb, fa = few(t, "cox_events_per_term"), few(ta, "adj_cox_events_per_term")
        out.append(f"- **Overfit risk** (fewer than {s.cox_events_per_term:g} events per model term): {len(fb)} of "
                   f"{len(t)} base-model fits" + (f", {len(fa)} of {len(ta)} adjusted fits" if adj_model else "")
                   + (f" (cohorts: {', '.join(sorted(set(fb.cohort) | set(fa.cohort)))})" if len(fb) or len(fa) else "")
                   + ".")
    if len(t) and s.cox_low_power_events > 0:
        weak = (lambda d, col: d[d[col] < s.cox_low_power_events] if col in d else d.iloc[:0])
        wb, wa = weak(t, "cox_events"), weak(ta, "adj_cox_events")
        hit = wb.index[wb.cox_p < s.alpha].union(wa.index[wa.adj_cox_p < s.alpha] if len(wa) else wa.index)
        sig = cells[cells.index.isin(hit)]                  # p < alpha in the base or the adjusted model
        out.append(f"- **Low power** (fewer than {s.cox_low_power_events} events; ‡ in the forests and after p in the "
                   f"ranked table): {len(wb)} of {len(t)} base-model fits"
                   + (f", {len(wa)} of {len(ta)} adjusted fits" if adj_model else "")
                   + (f"; with p < {s.alpha:g}{' (base or adjusted model)' if adj_model else ''}: "
                      f"{', '.join(name(r) for r in sig.itertuples())}" if len(sig) else "")
                   + ". There, p ≥ 0.05 says little against an association.")
    gone = (lambda d, col: d[d[col].fillna("").astype(str).str.contains("host expression left out", regex=False)]
            if col in d else d.iloc[:0])
    lb, la = gone(t, "cox_notes"), gone(ta, "adj_cox_notes")
    if len(lb) or len(la):
        out.append(f"- **Host expression left out** (no values for the host gene in the cohort, recorded for fewer "
                   f"than {s.covariate_min_complete:.0%} of its patients, or constant; those fits have no expression "
                   "term): "
                   f"{len(lb)} of {len(t)} base-model fits"
                   + (f", {len(la)} of {len(ta)} adjusted fits" if adj_model else "")
                   + f" (cohorts: {', '.join(sorted(set(lb.cohort) | set(la.cohort)))}).")
    if s.cox_min_events_per_term > 0:
        cut = (lambda col: cells[cells[col].eq("too_few_events_per_term")] if col in cells else cells.iloc[:0])
        cb, ca = cut("cox_status"), cut("adj_cox_status")
        if len(cb) or len(ca):
            out.append(f"- **Not fitted: too few events per term** (fewer than {s.cox_min_events_per_term:g} events "
                       f"per estimated term, `cox_min_events_per_term`): {len(cb)} base-model and {len(ca)} adjusted "
                       f"fits (cohorts: {', '.join(sorted(set(cb.cohort) | set(ca.cohort)))}).")
    if (len(t) or len(tk)) and s.ph_note_below > 0:
        low = (lambda d, col: d[d[col] < s.ph_note_below] if col in d else d.iloc[:0])
        pb, pk = low(t, "ph_p"), low(tk, "km_ph_p")
        sig = pb[pb.cox_p < s.alpha]
        out.append(f"- **Non-proportional hazards** (Schoenfeld test on the Kaplan–Meier time scale, p below "
                   f"{s.ph_note_below:g}; the HR then averages an effect that changes over follow-up): the {v} term in "
                   f"{len(pb)} of {len(t)} base-model fits"
                   + (f" (with p < {s.alpha:g}: {', '.join(name(r) for r in sig.itertuples())})" if len(sig) else "")
                   + f"; the KM split in {len(pk)} of {len(tk)} log-rank tests. Look at the KM curves of those cells.")
    if not out:
        return out
    stopped = s.cox_min_events_per_term > 0 and any(
        cells.get(c, pd.Series(dtype=object, index=cells.index)).eq("too_few_events_per_term").any()
        for c in ("cox_status", "adj_cox_status"))
    return ["## Notes on the survival tests", "", *out,
            "- All are notes on the pages and in `cells.csv`; every test ran"
            + (" (except the fits the events-per-term minimum stopped)" if stopped else "")
            + " and no result is removed. Settings: "
            "`narrow_psi_below`, `narrow_psi_measure` (default: the HR's unit, `psi_hr_unit`), `cox_events_per_term`, "
            "`cox_low_power_events`, `ph_note_below`.", ""]


def _fmt(v, nd=2):
    return "–" if v is None or (isinstance(v, float) and not np.isfinite(v)) else f"{v:.{nd}f}"


def _gex_cox_line(gex_cells: pd.DataFrame, s: Settings) -> str:
    """The report's line on the host genes' own Cox fits: every cohort with p < alpha (by p), ‡ marking low power."""
    t = gex_cells[gex_cells.cox_status.eq("tested")] if "cox_status" in gex_cells else gex_cells.iloc[:0]
    sig = t[t.cox_p < s.alpha].sort_values("cox_p") if len(t) else t
    low = sig.cox_low_power.map(_flag) if "cox_low_power" in sig else pd.Series(False, index=sig.index)
    gene = (lambda r: f"{r.gene} " if gex_cells.gene.nunique() > 1 else "")
    names = [f"{gene(r)}{r.cohort} HR {r.hr_per_sd:.2f}" + (" ‡" if lo else "") for r, lo in zip(sig.itertuples(), low)]
    return (f"- **Cox:** {len(t)} cohort{'' if len(t) == 1 else 's'} tested; expression has p < {s.alpha:g} in "
            f"{len(sig)}" + (f" ({', '.join(names)})" if names else "")
            + (f"; ‡ fewer than {s.cox_low_power_events} events" if low.any() else "") + ".")


def _afe_ale_pairs(corr: pd.DataFrame, ds: Dataset) -> bool:
    """Whether any event pair is two alternative first exons, or two alternative last exons, of one gene."""
    d = corr[corr.kind.eq("event")]
    if not len(d):
        return False
    t = ds.events.event_type.astype(str).str.upper()
    a, b = d.event_id.map(t), d.partner.map(t)
    return bool((a.eq(b) & a.isin(["AFE", "ALE"])).any())


def _rho_text(v) -> str:
    return f"{v:.2f}".replace("-", "−")


def _corr_lines(corr: pd.DataFrame, ds: Dataset, s: Settings, gex: bool = False, rare: int = 0,
                n_events: int = 0) -> list[str]:
    """The report's section on correlation (correlation=True): what was computed, the tests with host expression
    against chance, the strong cells and pairs (|rho| >= corr_note_above), and how to read them. `gex`: whether the
    overview and the expression page show the host gene; `rare`: events left out as rarely observed, of `n_events`."""
    from .correlation import strong
    lab = (lambda e: ds.events.at[e, "label"])
    out = ["", "## Correlation (Spearman, within each cohort)", ""]
    if not len(corr):
        if rare and rare >= n_events:
            return out + [f"- **Nothing to correlate:** every event is observed in under {s.observed_frac:.0%} of "
                          "every cohort's survival samples (`--min-observed`)."]
        return out + ["- **Nothing to correlate:** no expression table, and no gene with two of the events "
                      + (f"kept (the {rare} rarely observed are left out, `--min-observed`)." if rare else "probed.")]
    cut = f"|ρ| ≥ {s.corr_note_above:g}"
    out.append("- **What:** Spearman ρ over each cohort's survival samples (one case sample per patient, the samples "
               f"the survival tests draw on), where at least {s.corr_min_n} patients have both values (`corr_min_n`).")
    x = corr[corr.kind.eq("expression")]
    if len(x):
        t = x[x.corr_status.eq("tested")]
        sig = int((t.corr_p < s.alpha).sum())
        st = strong(t, s)
        listed = ", ".join(f"{lab(r.event_id)} {r.cohort} (ρ {_rho_text(r.rho)})" for r in st.head(12).itertuples()) \
            + (f", and {len(st) - 12} more" if len(st) > 12 else "")
        out.append(f"- **With host-gene expression:** {len(t)} of {len(x)} event × cohort cells tested; {sig} "
                   f"{_have(sig)} p < {s.alpha:g}, against about {s.alpha * len(t):.0f} by chance"
                   + ((f"; {cut} in {len(st)}: {listed}" if len(st) else f"; none has {cut}")
                      if s.corr_note_above > 0 else "") + ".")
    e = corr[corr.kind.eq("event")]
    if len(e):
        t = e[e.corr_status.eq("tested")]
        n_coh = t.groupby(["event_id", "partner"], sort=False).size()
        items = []
        for (a, b), d in strong(t, s).groupby(["event_id", "partner"], sort=False):      # strongest pair first
            lo, hi = d.rho.min(), d.rho.max()
            rng = _rho_text(lo) if lo == hi else f"{_rho_text(lo)} to {_rho_text(hi)}"
            items.append(f"{lab(a)}–{lab(b)} in {len(d)} of {n_coh[(a, b)]} cohorts (ρ {rng})")
        listed = ", ".join(items[:12]) + (f", and {len(items) - 12} more" if len(items) > 12 else "")
        out.append(f"- **Between events of a gene:** {len(t)} of {len(e)} pair × cohort cells tested "
                   f"({e[['event_id', 'partner']].drop_duplicates().shape[0]} pairs)"
                   + ((f"; {cut} for {len(items)} pair{'' if len(items) == 1 else 's'}: {listed}" if items else
                       f"; no pair has {cut}") if s.corr_note_above > 0 else "")
                   + ". Events of one gene share patients and often reads, so these are not compared with chance.")
    return out + [
        "- **Reading:**",
        "  - An event strongly correlated with its host gene's expression carries much of the same information. The "
        "Cox models that adjust for host expression (the default) then estimate what the event adds beyond "
        "expression, less precisely (a wider CI). The KM split is not adjusted: a KM hit of such an event could be "
        "expression's" + (" (compare the gene's row in the overview and its expression page)." if gex else "."),
        "  - Strongly correlated events of one gene likely measure one isoform change (rMATS often lists an exon "
        "several times with different flanking exons): count their hits as one piece of evidence.",
        *(["  - The PSI values of a gene's alternative first (last) exons share its first (last) exon use. Two such "
           "exons sum to 1, so their ρ is near −1 by construction; with more of them, pairs tend to be negative but "
           "need not be."] if _afe_ale_pairs(corr, ds) else []),
        "  - ρ describes how two values move together across patients, not why. p is nominal; q (Benjamini–Hochberg "
        "within each gene and kind) is in `correlations.csv`.",
        "- **Where:** `correlations.csv` (every pair × cohort), `correlation.png` (each cohort, after the gene maps in "
        "`probe.pdf`), the gene maps' triangles (the median over the cohorts of each pair tested), `expr_rho` in "
        "`cells.csv`, and "
        "`best_expr_rho` (ρ with expression in the best cohort) in `events.csv`."]


def _report(ds, cells, ranked, ep, base_model, adj_model, found, s, top, max_pages, call, ptab=None,
            gex_cells=None, hit_left: int = 0, corr=None, maps=(), rare: int = 0) -> str:
    from .plot.style import fp

    from .plot.style import in_sentence
    cmp_ = f"{in_sentence(ds.labels['case'])}–{in_sentence(ds.labels['reference'])}"   # as named in the data
    n_cox = int(cells.cox_status.eq("tested").sum())
    n_km = int(cells.km_status.eq("tested").sum())
    n_km_sig = int((cells.km_status.eq("tested") & (cells.km_p < s.alpha)).sum())
    n_sig = int((cells.cox_status.eq("tested") & (cells.cox_p < s.alpha)).sum())
    n_adj = int(cells.get("adj_cox_status", pd.Series(dtype=object)).eq("tested").sum()) if adj_model else 0
    n_adj_sig = int((cells.get("adj_cox_status", pd.Series(dtype=object)).eq("tested") &
                     (cells.get("adj_cox_p", pd.Series(dtype=float)) < s.alpha)).sum()) if adj_model else 0
    meas = ranked[ranked.measurable] if len(ranked) else ranked
    v = _value_name(ranked.event_type) if len(ranked) else "PSI"
    has_hit = bool(len(ranked)) and bool(ranked.event_type.map(opt_in).any())
    U = s.psi_hr_unit.upper()
    L = [f"# splice-assay probe · {', '.join(sorted(set(ranked.gene)))[:120]} · {ep}", "",
         f"splice-assay {__version__}. **Everything here is nominal: a probe ranks candidates for a closer look; it "
         "does not test a hypothesis.**", "",
         "## What was run", "",
         f"- **Events:** {len(ranked)} ({len(meas)} measurable in at least one cohort)."
         + (f" {hit_left} HIT-index event{'s were' if hit_left != 1 else ' was'} left out: the HIT index covers every "
            "exon, a far larger set; `--include-hit` adds them." if hit_left else "")
         + (f" {rare} {'are' if rare != 1 else 'is'} observed in under {s.observed_frac:.0%} of every cohort's "
            "survival samples and left out of the gene maps" + (" and correlations" if corr is not None else "")
            + " (`--min-observed`)." if rare else ""),
         f"- **Cohorts:** {cells.cohort.nunique()}.",
         (f"- **Subset:** {ds.notes['subset']['patients']} of {ds.notes['subset']['of_patients']} patients "
          f"({ds.notes['subset']['samples']} samples, their normals included), selected by "
          f"{ds.notes['subset']['by']}."
          if (ds.notes or {}).get("subset") else None),
         f"- **Endpoint:** {ep}.",
         (f"- **Settings changed from the defaults:** {s.changed_text()}."
          + (" Results that pass only these relaxed gates rest on fewer samples or events, or less variation, than "
             "the defaults require." if _relaxed(s, corr is not None) else "")
          + (" The KM arms are split by the rule given, not at the median." if KM_RULES & set(s.changed()) else "")
          if s.changed() else None),
         f"- **Base model** (the forest on every page): Cox {_describe(base_model, ds, cells.event_id.unique(), s)}.",
         ("- **Adjusted model** (the model rows on every page and the ranking): Cox "
          f"{_describe(adj_model, ds, cells.event_id.unique(), s)}"
          + (f"; baselines {', '.join(f'{k} {v}' for k, v in adj_model.baseline)}" if adj_model.baseline else "")
          + (f"; read from {found_summary(ds.clinical, found)}" if found else "") + "."
          if adj_model else "- **Adjusted model:** none (no clinical table, or no age/sex/stage column found)."),
         _ridge_line(cells, s, v, bool(base_model.covariates)),
         f"- **Pages:** one per ranked, measurable event (at most {max_pages}). Each shows "
         + ("every cohort with a test" if top >= ALL else f"the event's {top} most promising cohorts")
         + f" (Cox p < {s.alpha:g} first, low-power fits after the others within each; then by adjusted Cox p) "
         "with their models, and every cohort in the forest. Change with `--top N` or `--top all`.",
         "",
         "## How much is chance",
         "",
         f"- **Base model:** {n_cox} Cox tests; {n_sig} {_have(n_sig)} p < {s.alpha:g}, where chance alone gives about "
         f"{s.alpha * n_cox:.0f}.",
         (f"- **Adjusted model:** {n_adj} tests; {n_adj_sig} {_have(n_adj_sig)} p < {s.alpha:g}, against about "
          f"{s.alpha * n_adj:.0f} "
          "by chance." if adj_model else ""),
         f"- **KM:** {n_km} log-rank tests; {n_km_sig} {_have(n_km_sig)} p < {s.alpha:g}, against about "
         f"{s.alpha * n_km:.0f} by "
         "chance.",
         (f"- **No Cox model could be fitted:** no cohort reaches {s.cox_min_n} patients and {s.cox_min_events} "
          f"events with this {v} measured (Settings cox_min_n, cox_min_events). The ranking uses the KM tests."
          if n_cox == 0 else None),
         "- **q values** are Benjamini–Hochberg within each gene, with each kind of test its own family over all "
         f"of the gene's events × cohorts tested for {ep}"
         + (" (HIT-index events form families of their own)" if has_hit else "") + ":",
         f"  - {cmp_} within patients (`paired_q`) and over all samples (`unpaired_q`);",
         "  - KM (`km_q`);",
         f"  - Cox on {v} in the base model (`cox_q`, the forest) and in the adjusted model (`adj_cox_q`, the model "
         "rows).",
         f"  - A family of fewer than {s.fdr_min_family} tests gets no q; `*_q_tests` gives each family's size. The "
         "pages print q beside p. Host-gene expression is not adjusted.",
         "- **What converging evidence looks like:**",
         "  - the same HR direction in several cohorts;",
         f"  - a {cmp_} (group) hit in the same cohort;",
         "  - an association that survives the clinical adjustment.",
         "",
         *_flag_lines(cells, ds, s, adj_model),
         f"## Ranked events (ranking: {RANKING})",
         "",
         "| Rank | Event | Type | Adj. p<.05 | Base p<.05 (↑/↓) | KM p<.05 | Expected | Group hits | Group + survival | "
         f"Best cohort | HR per {U} | p |" + (" Protein (suggested) |" if ptab is not None else "") + " Page |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|" + ("---|" if ptab is not None else "") + "---|"]
    flat = marked = False
    for r in meas.head(max(40, max_pages)).itertuples():
        page = f"[page]({r.page})" if r.page else ""
        spread = getattr(r, f"best_psi_{s.psi_hr_unit}")
        if np.isfinite(spread) and spread == 0:
            flat = True
            page += " §"
        prot = ""
        if ptab is not None:
            pc = getattr(r, "protein_change", "")
            prot = f" {pc if isinstance(pc, str) and pc else '–'} |"
        p_txt = (fp(r.best_p) + (" (KM)" if r.best_model == "KM" else "") + (" ‡" if r.best_low_power else "")) \
            if np.isfinite(r.best_p) else "–"
        low_a, low_b = r.adj_cox_p05_low_power, r.cox_p05_low_power
        marked = marked or bool(low_a or low_b or r.best_low_power)
        L.append(f"| {r.rank} | {r.label} ({r.gene}) | {r.event_type} | {r.adj_cox_p05}"
                 + (f" ({low_a}‡)" if low_a else "") + f" | {r.cox_p05} ({r.cox_p05_hr_up}/{r.cox_p05_hr_down}"
                 + (f"; {low_b}‡" if low_b else "") + f") | {r.km_p05} | {r.expected_by_chance:.1f} | {r.group_hits} | "
                 f"{r.group_and_survival} | {r.best_cohort or '–'} | {_fmt(getattr(r, f'best_hr_per_{s.psi_hr_unit}'))} "
                 f"| {p_txt} |{prot} {page} |")
    if marked:
        L += ["", f"‡ low power (fewer than {s.cox_low_power_events} Cox events): after p, the best cohort's fit; in "
                  "the counts, how many of the hits are in such fits. The ranking counts those hits after the others. "
                  "For the best cohort and the pages it prefers a hit in a fit that is not low power, then a "
                  "low-power hit, then the other fits (those not low power first)."]
    if flat:
        L += ["", f"§ {v} barely varies in the best cohort ({U} 0): the association rests on a few samples, and "
                  f"there is no HR per {U}."]
    unmeas = ranked[~ranked.measurable] if len(ranked) else ranked
    if len(unmeas):
        L += ["", f"Not measurable in any cohort (coverage or gates): {len(unmeas)} event(s), listed in events.csv."]
    if ptab is not None:
        from .protein import SHOWN
        n_ok = int(ptab.protein_status.isin(SHOWN).sum())
        L += ["", "## Protein changes (suggested from annotation)", "",
              f"- **What:** for {n_ok} of {len(ptab)} events, the annotated transcripts that best represent the "
              "inclusion and exclusion forms (SpliceImpactR's matching), and what the event does to the protein.",
              "- **How sure:** these are suggestions read from annotation, not measurements. A form without an "
              "annotated transcript is described from the coding coordinates of the other one. An event with no "
              "matching transcript, or only noncoding ones, has no suggestion.",
              "- **Where:** the Protein column above; `proteins.csv` (transcripts, how they matched, residues, "
              "features); the band under the schematic on each page."]
    if gex_cells is not None and len(gex_cells):
        L += ["", "## Host-gene expression", "",
              f"- **What:** the same tests for the host gene's own expression ({', '.join(sorted(set(gex_cells.gene)))}): "
              "case vs reference, a KM split, and Cox on expression per SD plus the adjusted model's clinical terms. "
              "They have their own page per gene (`pages/<gene>_expression.png`, last in `probe.pdf`): a forest of "
              "every cohort, then the cohorts of the gene's splicing pages in full.",
              _gex_cox_line(gex_cells, s),
              "- **Reading:** a splicing association that holds while expression itself is not prognostic (or the "
              "reverse) is easier to interpret; the splicing models adjust for expression wherever it is recorded (see "
              "the notes above).",
              "- **Where:** `expression_cells.csv`."]
    if corr is not None:
        L += _corr_lines(corr, ds, s, gex_cells is not None and len(gex_cells) > 0, rare, len(ranked))
    L += ["", "## Files", "",
          "| File | Content |", "|---|---|",
          "| `events.csv` | One row per event: counts across cohorts, best cohort, page |",
          "| `cells.csv` | One row per event × cohort: group tests, KM, base and adjusted Cox (`adj_*`), q values |",
          "| `overview.png` | Events × cohorts: HR colour, p < 0.05 dot, group-hit frame"
          + ("; the host gene's expression in its own row" if gex_cells is not None and len(gex_cells) else "") + " |",
          *([f"| `gene_map_<GENE>.png` ({', '.join(m.removeprefix('gene_map_') for m in maps)}) | Each gene's model "
             "and its probed events, 5′ to 3′, beside their cells of the overview"
             + (", and the median ρ between them" if corr is not None and len(corr) else "") + " |"] if maps else []),
          "| `pages/` | One assay page per ranked, measurable event, at most " + str(max_pages)
          + " (`--max-pages`; SVG, PDF, PNG, CSV of plotted values, provenance) |",
          "| `probe.pdf` | The overview, " + ("the gene maps, " if maps else "")
          + ("the correlation figure, " if corr is not None and len(corr) else "") + "then every page |",
          *(["| `proteins.csv` | One row per event: suggested transcripts, protein change, features |"]
            if ptab is not None else []),
          *(["| `expression_cells.csv` | One row per host gene × cohort: its own group tests, KM and Cox |"]
            if gex_cells is not None else []),
          "| `agent_prompt.md` | What an agent needs to write a narrative of this probe: instructions, rules, this "
          "report and the best-ranked events' notable cohorts (`splice-assay summarize` hands it to Claude Code) |",
          *(["| `correlations.csv` | One row per pair × cohort: Spearman ρ of an event with its host gene's expression "
             "or with another event of its gene |"] if corr is not None else []),
          *(["| `correlation.png` | Each cohort's ρ as a grid, gene by gene: ρ colour and value (black where p < 0.05), "
             "strong-ρ frame |"]
            if corr is not None and len(corr) else []),
          "", "## Reproduce", "", "```bash", call or "splice-assay probe <data> --gene <gene>", "```", ""]
    return "\n".join(x for x in L if x is not None) + "\n"
