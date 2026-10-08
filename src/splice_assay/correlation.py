"""Correlation of splicing events with their host gene's expression and with the other events of their gene (opt-in:
`probe --correlation`, or `correlations(ds, ...)` in Python).

In each cohort, over its survival samples (one case sample per patient, the samples the survival tests draw on),
Spearman's rho between

    expression   an event's values (PSI, or the HIT index) and its host gene's expression;
    event        two events of the same gene.

A pair is tested when at least Settings.corr_min_n patients have both values and each value differs from its modal
value in at least Settings.min_off_modal of them (the survival tests' gate), so one or two patients cannot make rho.
Every survival sample with both values counts, whether or not its patient enters an endpoint's Cox fit (a survival record,
complete covariates), so n can exceed the fit's.
q values: Benjamini-Hochberg within each gene, per kind, over its pairs x cohorts; pairs with a HIT index form families
of their own (as in analysis.py). A correlation describes whether two values move together across patients; it does
not say why.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .config import Settings
from .dataset import Dataset
from .events import quantity

CORR_FIRST = ["kind", "gene", "event_id", "partner", "cohort", "corr_status", "corr_n", "rho", "corr_p", "corr_q",
              "corr_q_tests"]


def _rho(x: np.ndarray, y: np.ndarray, s: Settings) -> tuple[str, int, float, float]:
    """(status, n, rho, p): Spearman's rho over the patients with both values."""
    from scipy.stats import spearmanr
    ok = np.isfinite(x) & np.isfinite(y)
    n = int(ok.sum())
    if n < s.corr_min_n:
        return "too_few_patients", n, np.nan, np.nan
    a, b = np.round(x[ok], s.round_decimals), np.round(y[ok], s.round_decimals)
    if not (np.ptp(a) > 0 and np.ptp(b) > 0):
        return "constant", n, np.nan, np.nan
    if min(n - np.unique(v, return_counts=True)[1].max() for v in (a, b)) < s.min_off_modal:
        return "few_off_modal", n, np.nan, np.nan          # rho would rest on a handful of patients
    rho, p = spearmanr(a, b)
    return "tested", n, float(rho), float(p)


def correlations(ds: Dataset, events=None, cohorts=None, settings: Settings | None = None, include_hit: bool = False,
                 expression: bool = True, within_gene: bool = True) -> pd.DataFrame:
    """One row per pair x cohort (CORR_FIRST first): each event with its host gene's expression (`expression`; needs
    an expression table) and each pair of events of one gene (`within_gene`), for the events given (default: every
    event, HIT-index events only with include_hit) in `cohorts` (default all). Status: tested, too_few_patients,
    constant, few_off_modal (fewer than min_off_modal patients away from a value's modal value), or no_expression
    (the expression table has no row for the host gene)."""
    from .analysis import _names, _order, default_events, gene_fdr
    s = settings or Settings()
    events = default_events(ds, include_hit) if events is None else \
        list(dict.fromkeys(_names(events, "event(s)", set(ds.psi.index))))
    cohorts = ds.cohorts if cohorts is None else _names(cohorts, "cohort(s)", set(ds.cohorts))
    ev = ds.events.loc[events]
    fam = ev.event_type.map(quantity)
    smp = {c: ds.samples.sample_id[ds.samples.cohort.eq(c) & ds.samples.survival_cohort].to_numpy() for c in cohorts}
    psi = {c: ds.psi.loc[events].reindex(columns=smp[c]).to_numpy(float) for c in cohorts}   # events x samples
    pos = {e: i for i, e in enumerate(events)}
    rows = []
    if expression and ds.expression is not None:
        hosts = [h for h in dict.fromkeys(ev.expression_gene) if h in ds.expression.index]
        expr = {c: ds.expression.loc[hosts].reindex(columns=smp[c]) for c in cohorts}   # host genes x samples
        for e in events:
            host = ev.at[e, "expression_gene"]
            rec = dict(kind="expression", gene=ev.at[e, "gene"], event_id=e, partner=host, _family=fam[e])
            for c in cohorts:
                if host not in ds.expression.index:
                    rows.append(dict(rec, cohort=c, corr_status="no_expression", corr_n=0))
                    continue
                st, n, r, p = _rho(psi[c][pos[e]], expr[c].loc[host].to_numpy(float), s)
                rows.append(dict(rec, cohort=c, corr_status=st, corr_n=n, rho=r, corr_p=p))
    if within_gene:
        for gene, ids in ev.groupby("gene", sort=False).groups.items():
            members = set(ids)
            ids = [e for e in events if e in members]                     # in the order given
            for i, a in enumerate(ids):
                for b in ids[i + 1:]:
                    rec = dict(kind="event", gene=gene, event_id=a, partner=b,
                               _family="HIT index" if "HIT index" in (fam[a], fam[b]) else "PSI")
                    for c in cohorts:
                        st, n, r, p = _rho(psi[c][pos[a]], psi[c][pos[b]], s)
                        rows.append(dict(rec, cohort=c, corr_status=st, corr_n=n, rho=r, corr_p=p))
    if not rows:
        return pd.DataFrame(columns=CORR_FIRST)
    return _order(gene_fdr(pd.DataFrame(rows), "corr", ["gene", "kind", "_family"], s.fdr_min_family)
                  .drop(columns="_family"), CORR_FIRST)


def strong(corr: pd.DataFrame, s: Settings) -> pd.DataFrame:
    """The tested pairs with |rho| >= Settings.corr_note_above (none when it is 0), strongest first."""
    if s.corr_note_above <= 0 or not len(corr):
        return corr.iloc[:0]
    t = corr[corr.corr_status.eq("tested")]
    return t[t.rho.abs() >= s.corr_note_above].sort_values("rho", key=np.abs, ascending=False, kind="stable")


def correlation_figure(corr: pd.DataFrame, events: pd.DataFrame, s: Settings, max_rows: int = 60):
    """Rows x cohorts, Spearman rho in colour (and printed where the cells have room), gene by gene in the order of
    their best rank: each event with its host gene's expression (in rank order), then each pair of the gene's events
    (by the better rank, then the other's). A printed value is black where p < alpha and grey otherwise (a dot marks
    p < alpha where no value is printed); |rho| >= corr_note_above is framed; grey cells were not tested. `events`:
    the ranked events (probe.rank_events). At most max_rows rows: the expression rows of the best-ranked events, at
    most half of them unless the pairs leave more (max(max_rows // 2, max_rows - pairs)), then the pairs among each
    gene's best-ranked events (adding its events one by one in rank order), the genes taking turns. None when there
    is nothing to draw."""
    import matplotlib
    from matplotlib.figure import Figure
    from matplotlib.lines import Line2D
    from matplotlib.patches import Rectangle

    from .plot import grid as G
    from .plot import style as S

    if not len(corr) or not len(events):
        return None
    ev = events[events.event_id.isin(set(corr.event_id) | set(corr.partner))]
    rank = dict(zip(ev.event_id, ev["rank"]))
    label = dict(zip(ev.event_id, ev.label.astype(str)))
    color = {e: S.TYPE_COLOR.get(str(t).upper(), S.TYPE_DEFAULT) for e, t in zip(ev.event_id, ev.event_type)}
    per_gene = {}                            # gene -> (its host, expression rows, pair rows), all in rank order
    for gene in dict.fromkeys(ev.gene):
        ids = ev.event_id[ev.gene.eq(gene)].tolist()
        x = corr[corr.kind.eq("expression") & corr.event_id.isin(ids)]
        pr = corr[corr.kind.eq("event") & corr.event_id.isin(ids) & corr.partner.isin(ids)]
        pairs = sorted({tuple(sorted((a, b), key=rank.get)) for a, b in zip(pr.event_id, pr.partner)},
                       key=lambda ab: (rank[ab[0]], rank[ab[1]]))
        rows_x = [(("rank", e), x[x.event_id.eq(e)].set_index("cohort")) for e in ids if e in set(x.event_id)]
        rows_p = [(("pair", a, b), pr[(pr.event_id.eq(a) & pr.partner.eq(b)) | (pr.event_id.eq(b) & pr.partner.eq(a))]
                   .set_index("cohort")) for a, b in pairs]
        if rows_x or rows_p:
            per_gene[gene] = (x.partner.iloc[0] if len(x) else "", rows_x, rows_p)
    n_all = sum(len(rx) + len(rp) for _, rx, rp in per_gene.values())
    n_pairs = sum(len(rp) for _, _, rp in per_gene.values())
    budget = max(max_rows // 2, max_rows - n_pairs)          # the expression rows: at most half while pairs wait
    by_rank = sorted((rank[k[1]], g, i) for g, (_, rx, _) in per_gene.items() for i, (k, _) in enumerate(rx))
    chosen = {(g, i) for _, g, i in by_rank[:max(budget, 0)]}             # the best-ranked events' expression rows
    keep_x = {g: [r for i, r in enumerate(rx) if (g, i) in chosen] for g, (_, rx, _) in per_gene.items()}
    keep_p = {g: [] for g in per_gene}
    budget = max_rows - len(chosen)
    queue = {g: sorted(rp, key=lambda r: (rank[r[0][2]], rank[r[0][1]]))   # the pairs among its best-ranked events:
             for g, (_, _, rp) in per_gene.items()}                     # by the worse rank, then the better
    while budget > 0 and any(queue.values()):                # then the pairs, gene by gene in turn
        for g in per_gene:
            if budget > 0 and queue[g]:
                keep_p[g].append(queue[g].pop(0))
                budget -= 1
    layout, n_data = [], 0                   # (kind, payload): gene and block heads, then rows of cells
    for g, (host, rx, rp) in per_gene.items():
        kept = [r for r in rp if any(r is k for k in keep_p[g])]       # back in rank order
        if not (keep_x[g] or kept):
            continue
        layout.append(("gene", g))
        for head, rows in ((f"with {host} expression", keep_x[g]), ("between its events", kept)):
            if rows:
                layout.append(("sub", head))
                layout += [("row", r) for r in rows]
                n_data += len(rows)
    if not n_data:
        return None
    cohorts = sorted(corr.cohort.unique())
    nc = len(cohorts)
    cw = min(0.26, max(0.12, 5.6 / max(nc, 1)))
    values = cw >= 0.15                                   # room to print rho in a cell
    size = 6.2

    def wide(text):
        return S.text_width(text, size, weight="bold")
    # a row's label in two columns around a mark: "rank | label" or "label × | label", the marks in one column
    keys = [p[0] for kind, p in layout if kind == "row"]
    right_w = max(wide(label[k[-1]]) for k in keys)
    left_w = max([wide(label[k[1]]) for k in keys if k[0] == "pair"] + [S.text_width("99", size)])
    label_w = max(left_w + 0.20 + right_w, max([S.text_width(p, 6.0) + 0.1 for kind, p in layout if kind == "sub"]
                                               + [0.8]))
    x_grid = 0.12 + label_w + 0.14
    x_right = x_grid - 0.10 - right_w                     # the right-hand labels start here
    bar_w = 1.6
    W = max(4.5, x_grid + max(nc * cw, bar_w) + 0.3)     # the colour bar may reach past a narrow grid
    x_all = sum(len(rx) for _, rx, _ in per_gene.values())
    x_kept, p_kept = len(chosen), sum(map(len, keep_p.values()))
    shown = ", then ".join(filter(None, (
        ("every expression row" if x_kept == x_all else "the best-ranked events' expression rows") if x_kept else "",
        ("every pair" if p_kept == n_pairs else "the pairs among each gene's best-ranked events") if p_kept else "")))
    more = f"; {n_data} of {n_all} rows: {shown} (all in correlations.csv)" if n_all > n_data else ""
    caption = (f"colour{' and value' if values else ''}: Spearman ρ in each cohort's survival samples (one case sample "
               "per patient, as in the survival tests)" + more)
    cap = G.clauses(caption, W - 0.24, 5.8)
    scale = G.rho_scale()
    cmap, norm = scale
    fair, strong_ = cmap(norm(0.42)), cmap(norm(0.8))                 # sample cells for the key
    items = ([(dict(face=fair, value=(".42", S.INK)), f"p < {s.alpha:g}"),
              (dict(face=fair, value=(".42", S.MUTED)), f"p ≥ {s.alpha:g}")] if values else
             G.p_dots(fair, s))
    if s.corr_note_above > 0:
        items.append((dict(face=strong_, frame=True, value=(".80", G.ink_on(strong_)) if values else None),
                      f"|ρ| ≥ {s.corr_note_above:g}"))
    items.append((dict(face=G.UNTESTED), f"not tested (under {s.corr_min_n} patients, or under {s.min_off_modal} "
                                         "off a value's mode)"))     # short enough for the key beside the bar
    lay = G.legend_layout(W, x_grid, bar_w, items)
    lab_h = max(S.text_width(c, 5.6) for c in cohorts) + 0.06
    top = 0.30 + 0.12 * len(cap) + 0.10 + lab_h
    pitch = {"gene": 0.30, "sub": 0.20, "row": 0.20}
    ys, y = [], 0.0
    for k, (kind, _) in enumerate(layout):
        y += 0.10 if kind == "gene" and k else 0.0        # a gap between genes
        ys.append(y)
        y += pitch[kind]
    n_h = y
    H = top + n_h + lay.bottom
    with matplotlib.rc_context(S.rc()):
        fig = Figure(figsize=(W, H))
        fx = (lambda x: x / W)                                            # noqa: E731
        fy = (lambda y_: 1 - y_ / H)                                      # noqa: E731
        ax = fig.add_axes([fx(x_grid), fy(top + n_h), fx(nc * cw), n_h / H])
        ax.set_xlim(0, nc)
        ax.set_ylim(n_h, 0)
        ax.axis("off")
        for j, c in enumerate(cohorts):
            ax.text(j + 0.5, -0.06, c, rotation=90, ha="center", va="bottom", fontsize=5.6, color=S.INK2,
                    clip_on=False)
        for (kind, payload), y0 in zip(layout, ys):
            if kind == "gene":
                fig.add_artist(Line2D([fx(0.12), fx(x_grid + nc * cw)], [fy(top + y0 + 0.03)] * 2, color=S.GRID,
                                      lw=0.6))
                fig.text(fx(0.12), fy(top + y0 + 0.17), payload, fontsize=7.0, fontweight="bold", style="italic",
                         va="center")
                continue
            if kind == "sub":
                fig.text(fx(0.12), fy(top + y0 + 0.11), payload, fontsize=6.0, color=S.INK2, va="center")
                continue
            key, d = payload
            yc = y0 + 0.10
            b = key[-1]
            fig.text(fx(x_right), fy(top + yc), label[b], fontsize=size, color=color[b], va="center",
                     fontweight="bold")
            if key[0] == "rank":
                fig.text(fx(x_right - 0.07), fy(top + yc), str(rank[b]), fontsize=size, color=S.MUTED, ha="right",
                         va="center")
            else:
                fig.text(fx(x_right - 0.10), fy(top + yc), "×", fontsize=size, color=S.INK2, ha="center",
                         va="center")
                fig.text(fx(x_right - 0.17), fy(top + yc), label[key[1]], fontsize=size, color=color[key[1]],
                         ha="right", va="center", fontweight="bold")
            for j, c in enumerate(cohorts):
                if c not in d.index:
                    continue
                r = d.loc[c]
                tested = r.corr_status == "tested" and np.isfinite(r.rho)
                ax.add_patch(Rectangle((j + 0.05, yc - 0.08), 0.9, 0.16, facecolor=cmap(norm(r.rho)) if tested
                                       else G.UNTESTED, lw=0))
                if not tested:
                    continue
                if s.corr_note_above > 0 and abs(r.rho) >= s.corr_note_above:
                    ax.add_patch(Rectangle((j + 0.05, yc - 0.08), 0.9, 0.16, fill=False, edgecolor=S.INK, lw=0.7))
                sig = np.isfinite(r.corr_p) and r.corr_p < s.alpha
                ink = G.ink_on(cmap(norm(r.rho)))
                if values:                                   # grey where p >= alpha, unless the cell needs white
                    ax.text(j + 0.5, yc, G.rho_text(r.rho), ha="center", va="center", fontsize=4.9,
                            color=ink if sig or ink == "white" else S.MUTED)
                elif sig:
                    ax.plot([j + 0.5], [yc], "o", ms=2.6 if r.corr_p >= 0.01 else 3.6, mew=0, color=ink)
        fig.text(fx(0.12), fy(0.10), "Probe correlation", fontsize=8.5, fontweight="bold", va="top")
        for k, line in enumerate(cap):
            fig.text(fx(0.12), fy(0.30 + 0.12 * k), line, fontsize=5.8, color=S.INK2, va="top")
        G.bar_and_key(fig, W, H, lay, x_grid, bar_w, items, scale, [-1, -0.5, 0, 0.5, 1], ["−1", "−.5", "0", ".5", "1"],
                      "Spearman ρ", ends=("opposite", "together"))
    return fig
