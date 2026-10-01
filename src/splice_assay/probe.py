"""Probe: every event of a gene (or every event given) in every cohort, ranked, with one assay page per event.

For one endpoint (default OS) and every event x cohort, the probe runs the group tests, the KM split and two Cox
models: the base model (PSI + host expression, the forest's model) and the adjusted model (by default + age + sex +
stage, found in the clinical table; see clinical.py). Then:

    cells.csv       one row per event x cohort: group tests, KM, base and adjusted Cox, q values (BH within each
                    gene, per kind of test)
    events.csv      one row per event, ranked (see RANKING)
    overview.png    events x cohorts: adjusted HR per IQR (colour), adjusted p < 0.05 (dot), group hit (frame)
    pages/          one assay page per ranked event: the event, its most promising cohorts, their models, the forest
    probe.pdf       the overview, then every page, in rank order
    report.md       what was run, the ranking, how to read it, and how to reproduce it
    proteins.csv    with a protein cache: the suggested protein change of every event (see protein.py), also
                    summarised in events.csv and drawn on the pages
    expression_cells.csv   with an expression table: the host gene's own statistics per cohort (case vs reference,
                    KM on its median split, Cox on expression + the adjusted model's clinical terms); also drawn as
                    expression rows on the pages

Everything is nominal: a probe ranks candidates, it does not test a hypothesis. The report states how many
p < 0.05 results chance alone would give.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from . import __version__
from .analysis import Results, analyse
from .clinical import auto_clinical
from .config import Settings
from .dataset import Dataset, InputError
from .stats.survival import CoxModel

ALL = 10_000                     # `top` meaning every cohort with a test
RANKING = ("measurable events first; then cohorts where the adjusted Cox p < 0.05, then where the base Cox p < 0.05, "
           "then cohorts with both a group hit and a survival hit, then cohorts where the KM p < 0.05, then the "
           "smallest p (adjusted Cox, else base Cox, else KM)")
ADJ_COLS = ["cox_status", "cox_model", "cox_n", "cox_events", "cox_n_dropped", "cox_notes", "hr_per_iqr",
            "ci_low_iqr", "ci_high_iqr", "cox_p", "cox_q", "cox_q_tests", "cox_events_per_term", "psi_narrow"]


@dataclass
class ProbeResult:
    cells: pd.DataFrame
    events: pd.DataFrame
    paths: dict = field(default_factory=dict)
    pages: list = field(default_factory=list)


def select_events(ds: Dataset, genes=None, events=None) -> list[str]:
    if events:
        miss = [e for e in events if e not in ds.psi.index]
        if miss:
            raise InputError(f"event(s) not in the data: {', '.join(miss[:5])}")
        return list(events)
    if genes:
        want = {str(g).lower() for g in genes}
        ev = ds.events
        hit = ev.index[ev.gene.str.lower().isin(want) | ev.gene_id.str.lower().isin(want)]
        if not len(hit):
            raise InputError(f"no event of gene(s) {', '.join(genes)} in the data")
        return [e for e in ds.psi.index if e in set(hit)]
    return list(ds.psi.index)


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


def rank_events(cells: pd.DataFrame, ds: Dataset, alpha: float) -> pd.DataFrame:
    rows = []
    for eid, d in cells.groupby("event_id", sort=False):
        t = d[d.cox_status.eq("tested")]
        sig = t[t.cox_p < alpha]
        has_adj = "adj_cox_p" in d
        ta = d[d.get("adj_cox_status", pd.Series(index=d.index, dtype=object)).eq("tested")] if has_adj else d.iloc[:0]
        sig_a = ta[ta.adj_cox_p < alpha] if has_adj else ta
        surv_hit = d.survival_hit.fillna(False).astype(bool)
        tk = d[d.km_status.eq("tested")]
        if len(ta):
            best, pcol, hcol, model = ta.sort_values("adj_cox_p"), "adj_cox_p", "adj_hr_per_iqr", "adjusted"
        elif len(t):
            best, pcol, hcol, model = t.sort_values("cox_p"), "cox_p", "hr_per_iqr", "base"
        else:                                           # no Cox model anywhere (e.g. too few deaths): the KM test
            best, pcol, hcol, model = tk.sort_values("km_p"), "km_p", None, "KM"
        b = best.iloc[0] if len(best) else None
        rows.append(dict(
            event_id=eid, gene=ds.events.at[eid, "gene"], label=ds.events.at[eid, "label"],
            event_type=ds.events.at[eid, "event_type"], cohorts_cox_tested=len(t),
            cox_p05=len(sig), cox_p05_hr_up=int((sig.hr_per_iqr > 1).sum()), cox_p05_hr_down=int((sig.hr_per_iqr < 1).sum()),
            expected_by_chance=round(alpha * len(t), 2), adj_cohorts_tested=len(ta), adj_cox_p05=len(sig_a),
            group_hits=int(d.group_hit.fillna(False).astype(bool).sum()),
            within_patient=int(d.within_patient_support.fillna(False).astype(bool).sum()),
            group_and_survival=int((d.group_hit.fillna(False).astype(bool) & surv_hit).sum()),
            share_hr_up=round(float((t.hr_per_iqr > 1).mean()), 2) if len(t) else np.nan,
            km_tested=len(tk), km_p05=int((tk.km_p < alpha).sum()),
            best_cohort=None if b is None else b.cohort,
            best_hr_per_iqr=np.nan if b is None or hcol is None else b[hcol],
            best_logrank_hr=np.nan if b is None or model != "KM" else b.logrank_hr,
            best_p=np.nan if b is None else b[pcol], best_model=model if b is not None else "",
            best_psi_iqr=np.nan if b is None or model == "KM" else b.psi_iqr,
            min_cox_q=float(t.cox_q.min()) if len(t) else np.nan))
    ev = pd.DataFrame(rows)
    if ev.empty:
        return ev
    ev["measurable"] = (ev.cohorts_cox_tested > 0) | (ev.km_tested > 0)
    ev = ev.sort_values(["measurable", "adj_cox_p05", "cox_p05", "group_and_survival", "km_p05", "best_p"],
                        ascending=[False, False, False, False, False, True], na_position="last",
                        kind="stable").reset_index(drop=True)
    ev.insert(0, "rank", np.arange(1, len(ev) + 1))
    return ev


def pick_cohorts(cells: pd.DataFrame, event: str, n: int = 3) -> list[str]:
    """The event's most promising cohorts: smallest adjusted Cox p, then base Cox p, then KM p."""
    d = cells[cells.event_id.eq(event)].copy()
    d["_a"] = d.get("adj_cox_p", pd.Series(np.nan, index=d.index)).where(
        d.get("adj_cox_status", pd.Series("", index=d.index)).eq("tested"))
    d["_b"] = d.cox_p.where(d.cox_status.eq("tested")) if "cox_p" in d else np.nan
    d["_k"] = d.km_p.where(d.km_status.eq("tested")) if "km_p" in d else np.nan
    d = d[d[["_a", "_b", "_k"]].notna().any(axis=1)]
    d = d.sort_values(["_a", "_b", "_k"], na_position="last")
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


def overview(cells: pd.DataFrame, events: pd.DataFrame, endpoint: str, s: Settings, max_rows: int = 60):
    """Events (rows, rank order) x cohorts: HR per IQR (adjusted when fitted, else base) in colour, p < alpha as a
    dot, a group hit as a frame; grey = not tested."""
    import matplotlib
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    from matplotlib.figure import Figure
    from matplotlib.patches import Rectangle

    from .plot import style as S

    ev = events[events.measurable].head(max_rows)
    cohorts = sorted(cells.cohort.unique())
    nr, nc = len(ev), len(cohorts)
    cell_w, cell_h = min(0.2, 5.6 / max(nc, 1)), 0.17
    label_w = max([S.text_width(f"{r.label}  {r.gene}", 6.0) for r in ev.itertuples()] + [0.8]) + 0.25
    W = max(4.5, label_w + nc * cell_w + 0.4)
    H = 0.95 + max(nr, 1) * cell_h + 0.75
    cmap = LinearSegmentedColormap.from_list("hr", ["#2f5f98", "#9ebbd9", "#f4f3ee", "#e6a88a", "#b0412c"])
    norm = Normalize(-1.5, 1.5)
    idx = cells.set_index(["event_id", "cohort"])
    with matplotlib.rc_context(S.rc()):
        fig = Figure(figsize=(W, H))
        ax = fig.add_axes([label_w / W, 0.75 / H, nc * cell_w / W, max(nr, 1) * cell_h / H])
        ax.set_xlim(0, nc)
        ax.set_ylim(max(nr, 1), 0)
        for i, r in enumerate(ev.itertuples()):
            for j, c in enumerate(cohorts):
                if (r.event_id, c) not in idx.index:
                    continue
                x = idx.loc[(r.event_id, c)]
                adj = "adj_cox_status" in x and x.adj_cox_status == "tested"
                hr, p = (x.adj_hr_per_iqr, x.adj_cox_p) if adj else (x.get("hr_per_iqr"), x.get("cox_p"))
                tested = adj or x.cox_status == "tested"
                face = cmap(norm(np.log2(hr))) if tested and np.isfinite(hr) else "#e8e7e1"
                ax.add_patch(Rectangle((j + 0.04, i + 0.06), 0.92, 0.88, facecolor=face, lw=0))
                gh = x.get("group_hit", False)
                if pd.notna(gh) and bool(gh):
                    ax.add_patch(Rectangle((j + 0.04, i + 0.06), 0.92, 0.88, fill=False, edgecolor=S.INK, lw=0.7))
                if tested and np.isfinite(p) and p < s.alpha:
                    ax.plot([j + 0.5], [i + 0.5], "o", ms=2.6 if p >= 0.01 else 3.6, color=S.INK, mew=0)
        ax.set_xticks(np.arange(nc) + 0.5, cohorts, rotation=90, fontsize=5.6)
        ax.set_yticks(np.arange(nr) + 0.5, [f"{r.rank}. {r.label}  {r.gene}" for r in ev.itertuples()], fontsize=6.0)
        ax.tick_params(length=0)
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.xaxis.tick_top()
        fig.text(0.12 / W, 1 - 0.1 / H, f"Probe overview · {endpoint}", fontsize=8.5, fontweight="bold", va="top")
        fig.text(0.12 / W, 1 - 0.3 / H, "colour: HR per IQR of PSI (adjusted model where fitted); dot: p < "
                 f"{s.alpha:g} (large: < 0.01); frame: group hit; grey: not tested", fontsize=5.8, color=S.INK2,
                 va="top")
        cax = fig.add_axes([label_w / W, 0.25 / H, min(1.6, nc * cell_w) / W, 0.08 / H])
        grad = np.linspace(-1.5, 1.5, 256)[None, :]
        cax.imshow(grad, aspect="auto", cmap=cmap, norm=norm, extent=(-1.5, 1.5, 0, 1))
        cax.set_yticks([])
        cax.set_xticks([-1, 0, 1], ["0.5", "1", "2"], fontsize=5.6)
        cax.set_xlabel("HR per IQR", fontsize=5.8, labelpad=1)
    return fig


def gene_events(events_table: pd.DataFrame, genes) -> list[str]:
    """The event IDs of some genes (by name or gene ID, case-insensitive) in an events table."""
    want = {str(g).lower() for g in genes}
    e = events_table
    gid = e["gene_id"].astype(str).str.lower() if "gene_id" in e.columns else pd.Series("", index=e.index)
    return e.loc[e["gene"].astype(str).str.lower().isin(want) | gid.isin(want), "event_id"].astype(str).tolist()


def probe(ds: Dataset, genes=None, events=None, cohorts=None, endpoint=None, *, settings: Settings | None = None,
          model: CoxModel | None = None, adjusted="auto", baseline: dict | None = None, gtf=None, out_dir=None,
          top: int = 3, max_pages: int = 30, log=print, call: str = "", proteins=None, gex: bool = True) -> ProbeResult:
    """Probe events (all of `genes`, the `events` given, or every event) in `cohorts` (default all) for one endpoint
    (default OS). `adjusted`: "auto" (age, sex and stage found in the clinical table), a CoxModel, or None.
    `proteins`: a protein cache (protein.ProteinCache or its folder) for the suggested protein changes. `gex`: the
    host gene's own expression statistics (expression_cells.csv, and rows on the pages) when expression is given."""
    from matplotlib.backends.backend_pdf import PdfPages

    from .annotation import read_gtf
    from .events import geometry
    from .plot import event_panel
    from .plot import style as S
    from .plot.panel_common import page_parts, safe_name

    s = settings or Settings()
    ep = default_endpoint(ds, endpoint)
    ev = select_events(ds, genes, events)
    cohorts = None if not cohorts or [str(c).lower() for c in cohorts] == ["all"] else list(cohorts)
    found = {}
    if adjusted == "auto":
        ds, adj_model, found = auto_clinical(ds)
        if adj_model is not None and baseline:
            adj_model = adj_model.with_clinical(baseline=baseline)
    else:
        adj_model = adjusted or None
    base_model = model or CoxModel()
    if s.changed():
        log(f"settings changed from the defaults: {s.changed_text()}")
    log(f"probe: {len(ev)} event(s) x {len(cohorts or ds.cohorts)} cohort(s), {ep}; base Cox "
        f"{base_model.describe(ds.expression is not None)}"
        + (f"; adjusted {adj_model.describe(ds.expression is not None)}" if adj_model else "; no adjusted model"))
    base = analyse(ds, events=ev, endpoints=[ep], cohorts=cohorts, settings=s, model=base_model)
    adj = analyse(ds, events=ev, endpoints=[ep], cohorts=cohorts, settings=s, model=adj_model) if adj_model else None
    cells = combine(base, adj, ep)
    ranked = rank_events(cells, ds, s.alpha)
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
    gex_cells = None
    genes_x = [g for g in dict.fromkeys(ds.events.loc[ev, "expression_gene"])
               if ds.expression is not None and g in ds.expression.index]
    if gex and genes_x:
        from .analysis import analyse_expression
        gm = CoxModel(expression=False, covariates=(adj_model or base_model).covariates,
                      categorical=(adj_model or base_model).categorical, strata=(adj_model or base_model).strata,
                      baseline=dict((adj_model or base_model).baseline))
        gex_cells = analyse_expression(ds, genes=genes_x, endpoints=[ep], cohorts=cohorts, settings=s,
                                       model=gm).cells()
    out = ProbeResult(cells, ranked)
    if out_dir is None:
        return out
    out_dir = Path(out_dir)
    (out_dir / "pages").mkdir(parents=True, exist_ok=True)
    cells.to_csv(out_dir / "cells.csv", index=False, lineterminator="\n")
    ranked_path = out_dir / "events.csv"
    gtf_table = None
    todo = ranked[ranked.measurable].head(max_pages) if len(ranked) else ranked
    drawable = [e for e in todo.event_id if ds.events.at[e, "variable"] and ds.events.at[e, "strand"] in ("+", "-")]
    if gtf is not None and drawable:
        wins, genes_ = [], set()
        for e in drawable:
            row = ds.events.loc[e]
            lo, hi = geometry(e, row).span
            wins.append((row.chrom, lo - s.gtf_flank, hi + s.gtf_flank))
            genes_ |= {g for g in (row.gene, row.gene_id) if g}
        gtf_table = read_gtf(gtf, wins, genes=sorted(genes_))
    pages, page_col = [], {}
    import matplotlib
    with matplotlib.rc_context(S.rc()), PdfPages(out_dir / "probe.pdf", metadata={"CreationDate": None,
                                                                                  "ModDate": None}) as pdf:
        fig = overview(cells, ranked, ep, s)
        S.save(fig, out_dir, "overview", s.formats, s.dpi)
        fig.savefig(pdf, format="pdf")
        for r in todo.itertuples():
            if r.event_id not in drawable:
                continue
            focus = [c for c in (cohorts or []) if c in set(cells.cohort)][:top] if cohorts and len(cohorts) <= top \
                else pick_cohorts(cells, r.event_id, top)
            if not focus:
                continue
            parts = page_parts(focus, s.cohorts_per_page)  # many cohorts: several pages, each with the full forest
            for i, chunk in enumerate(parts, 1):
                stem = f"{r.rank:03d}_{safe_name(r.event_id)}" + (f"_p{i}" if len(parts) > 1 else "")
                p = event_panel(ds, r.event_id, chunk, ep, settings=s, model=base_model, gtf=gtf_table,
                                results=_slice(base, r.event_id), detail=adj_model, out_dir=out_dir / "pages",
                                stem=stem, proteins={r.event_id: pchanges[r.event_id]} if r.event_id in pchanges
                                else None, gex=gex, detail_results=_slice(adj, r.event_id) if adj is not None else None,
                                part=(i, len(parts)) if len(parts) > 1 else None)
                p.figure.savefig(pdf, format="pdf")
                pages.append((r.event_id, stem, chunk))
                page_col.setdefault(r.event_id, f"pages/{stem}.png")
                log(f"  page {len(pages)}: {r.event_id} ({', '.join(chunk)})")
    ranked["page"] = ranked.event_id.map(page_col).fillna("")
    ranked.to_csv(ranked_path, index=False, lineterminator="\n")
    report = _report(ds, cells, ranked, ep, base_model, adj_model, found, s, top, max_pages, call, ptab, gex_cells)
    (out_dir / "report.md").write_text(report)
    out.events, out.pages = ranked, pages
    out.paths = dict(report=out_dir / "report.md", events=ranked_path, cells=out_dir / "cells.csv",
                     overview=out_dir / "overview.png", pdf=out_dir / "probe.pdf", pages=out_dir / "pages")
    if ptab is not None:
        ptab.to_csv(out_dir / "proteins.csv", index=False, lineterminator="\n")
        out.paths["proteins"] = out_dir / "proteins.csv"
    if gex_cells is not None:
        gex_cells.to_csv(out_dir / "expression_cells.csv", index=False, lineterminator="\n")
        out.paths["expression"] = out_dir / "expression_cells.csv"
    return out


def _slice(res: Results, event: str) -> Results:
    """One event's rows of gene-wide results (its q values keep their gene-wide families)."""
    return Results(res.groups[res.groups.event_id.eq(event)].reset_index(drop=True),
                   res.survival[res.survival.event_id.eq(event)].reset_index(drop=True),
                   res.cox_terms[res.cox_terms.event_id.eq(event)].reset_index(drop=True) if len(res.cox_terms)
                   else res.cox_terms, settings=res.settings, model=res.model)


def _have(n: int) -> str:
    return "has" if n == 1 else "have"


def _flag_lines(cells: pd.DataFrame, ds: Dataset, s: Settings, adj_model) -> list[str]:
    """The report's flags on the Cox fits: narrow PSI ranges and few events per term (notes only; the fits ran)."""
    out = []
    t = cells[cells.cox_status.eq("tested")]
    ta = cells[cells.get("adj_cox_status", pd.Series(dtype=object, index=cells.index)).eq("tested")] if adj_model \
        else cells.iloc[:0]
    if not len(t) or not (s.narrow_psi_below > 0 or s.cox_events_per_term > 0):
        return out
    name = (lambda r: f"{ds.events.at[r.event_id, 'label']} {r.cohort}")
    if s.narrow_psi_below > 0:
        nar = t[t.psi_narrow.fillna(False).astype(bool)]
        sig = nar[nar.cox_p < s.alpha]
        out.append(f"- **Narrow PSI range** ({s.narrow_psi_measure.upper()} of PSI below {s.narrow_psi_below:g} in the "
                   f"fit cohort, so the HR covers a few PSI points): {len(nar)} of {len(t)} base-model fits"
                   + (f"; with p < {s.alpha:g}: {', '.join(name(r) for r in sig.itertuples())}" if len(sig) else "")
                   + ".")
    if s.cox_events_per_term > 0:
        few = (lambda d, col: d[d[col] < s.cox_events_per_term] if col in d else d.iloc[:0])
        fb, fa = few(t, "cox_events_per_term"), few(ta, "adj_cox_events_per_term")
        out.append(f"- **Overfit risk** (fewer than {s.cox_events_per_term:g} events per model term): {len(fb)} of "
                   f"{len(t)} base-model fits" + (f", {len(fa)} of {len(ta)} adjusted fits" if adj_model else "")
                   + (f" (cohorts: {', '.join(sorted(set(fb.cohort) | set(fa.cohort)))})" if len(fb) or len(fa) else "")
                   + ".")
    return ["## Flags on the Cox fits", "", *out,
            "- Both are notes on the page's model header; the fits ran. Settings: `narrow_psi_below`, "
            "`narrow_psi_measure`, `cox_events_per_term`.", ""]


def _fmt(v, nd=2):
    return "–" if v is None or (isinstance(v, float) and not np.isfinite(v)) else f"{v:.{nd}f}"


def _report(ds, cells, ranked, ep, base_model, adj_model, found, s, top, max_pages, call, ptab=None,
            gex_cells=None) -> str:
    from .plot.style import fp

    has_expr = ds.expression is not None
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
    L = [f"# splice-assay probe · {', '.join(sorted(set(ranked.gene)))[:120]} · {ep}", "",
         f"splice-assay {__version__}. **Everything here is nominal: a probe ranks candidates for a closer look; it "
         "does not test a hypothesis.**", "",
         "## What was run", "",
         f"- **Events:** {len(ranked)} ({len(meas)} measurable in at least one cohort).",
         f"- **Cohorts:** {cells.cohort.nunique()}.",
         (f"- **Subset:** {ds.notes['subset']['patients']} of {ds.notes['subset']['of_patients']} patients "
          f"({ds.notes['subset']['samples']} samples, their normals included), selected by "
          f"{ds.notes['subset']['by']}."
          if (ds.notes or {}).get("subset") else None),
         f"- **Endpoint:** {ep}.",
         (f"- **Settings changed from the defaults:** {s.changed_text()}. Results that pass only these relaxed gates "
          "rest on fewer samples or events than the defaults require." if s.changed() else None),
         f"- **Base model** (the forest on every page): Cox {base_model.describe(has_expr)}.",
         ("- **Adjusted model** (the model rows on every page and the ranking): Cox "
          f"{adj_model.describe(has_expr)}"
          + (f"; baselines {', '.join(f'{k} {v}' for k, v in adj_model.baseline)}" if adj_model.baseline else "")
          + (f"; columns found: {', '.join(f'{k} = {v}' for k, v in found.items())}" if found else "") + "."
          if adj_model else "- **Adjusted model:** none (no clinical table, or no age/sex/stage column found)."),
         f"- **Pages:** one per ranked, measurable event (at most {max_pages}). Each shows "
         + ("every cohort with a test" if top >= ALL else f"the event's {top} most promising cohorts")
         + " (ordered by adjusted Cox p) with their models, and every cohort in the forest. Change with `--top N` "
         "or `--top all`.",
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
          "events with this PSI measured (Settings cox_min_n, cox_min_events). The ranking uses the KM tests."
          if n_cox == 0 else None),
         "- **q values** are Benjamini–Hochberg within each gene, with each kind of test its own family over all "
         f"of the gene's events × cohorts tested for {ep}:",
         f"  - {cmp_} within patients (`paired_q`) and over all samples (`unpaired_q`);",
         "  - KM (`km_q`);",
         "  - Cox on PSI in the base model (`cox_q`, the forest) and in the adjusted model (`adj_cox_q`, the model "
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
         "Best cohort | HR per IQR | p |" + (" Protein (suggested) |" if ptab is not None else "") + " Page |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|" + ("---|" if ptab is not None else "") + "---|"]
    flat = False
    for r in meas.head(max(40, max_pages)).itertuples():
        page = f"[page]({r.page})" if r.page else ""
        if np.isfinite(r.best_psi_iqr) and r.best_psi_iqr == 0:
            flat = True
            page += " †"
        prot = ""
        if ptab is not None:
            pc = getattr(r, "protein_change", "")
            prot = f" {pc if isinstance(pc, str) and pc else '–'} |"
        p_txt = (fp(r.best_p) + (" (KM)" if r.best_model == "KM" else "")) if np.isfinite(r.best_p) else "–"
        L.append(f"| {r.rank} | {r.label} ({r.gene}) | {r.event_type} | {r.adj_cox_p05} | {r.cox_p05} "
                 f"({r.cox_p05_hr_up}/{r.cox_p05_hr_down}) | {r.km_p05} | {r.expected_by_chance:.1f} | {r.group_hits} | "
                 f"{r.group_and_survival} | {r.best_cohort or '–'} | {_fmt(r.best_hr_per_iqr)} | {p_txt} |{prot} {page} |")
    if flat:
        L += ["", "† PSI barely varies in the best cohort (IQR 0): the association rests on a few samples, and there "
                  "is no HR per IQR."]
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
        t = gex_cells[gex_cells.cox_status.eq("tested")] if "cox_status" in gex_cells else gex_cells.iloc[:0]
        sig = t[t.cox_p < s.alpha] if len(t) else t
        L += ["", "## Host-gene expression", "",
              f"- **What:** the same tests for the host gene's own expression ({', '.join(sorted(set(gex_cells.gene)))}): "
              "case vs reference, a KM split at its median, and Cox on expression per SD plus the adjusted model's "
              "clinical terms. Each page shows them under the splicing rows.",
              f"- **Cox:** {len(t)} cohort{'' if len(t) == 1 else 's'} tested; expression has p < {s.alpha:g} in "
              f"{len(sig)}"
              + (f" ({', '.join(f'{r.cohort} HR {r.hr_per_sd:.2f}' for r in sig.sort_values('cox_p').head(6).itertuples())})"
                 if len(sig) else "") + ".",
              "- **Reading:** a splicing association that holds while expression itself is not prognostic (or the "
              "reverse) is easier to interpret; the splicing models already adjust for expression.",
              "- **Where:** `expression_cells.csv`."]
    L += ["", "## Files", "",
          "| File | Content |", "|---|---|",
          "| `events.csv` | One row per event: counts across cohorts, best cohort, page |",
          "| `cells.csv` | One row per event × cohort: group tests, KM, base and adjusted Cox (`adj_*`), q values |",
          "| `overview.png` | Events × cohorts: HR colour, p < 0.05 dot, group-hit frame |",
          "| `pages/` | One assay page per ranked event (SVG, PDF, PNG, CSV of plotted values, provenance) |",
          "| `probe.pdf` | The overview, then every page |",
          *(["| `proteins.csv` | One row per event: suggested transcripts, protein change, features |"]
            if ptab is not None else []),
          *(["| `expression_cells.csv` | One row per host gene × cohort: its own group tests, KM and Cox |"]
            if gex_cells is not None else []),
          "", "## Reproduce", "", "```bash", call or "splice-assay probe <data> --gene <gene>", "```", ""]
    return "\n".join(x for x in L if x is not None) + "\n"
