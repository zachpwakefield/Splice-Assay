"""Command line: splice-assay <command> ...  (run `splice-assay <command> -h` for the options)."""
from __future__ import annotations

import argparse
import os
import shlex
import sys
import warnings
from pathlib import Path

import pandas as pd

from . import __version__
from .config import Settings
from .dataset import EVENT_COLUMNS, Dataset, InputError, read_table
from .stats.survival import CoxModel


def _settings(args) -> Settings:
    s = Settings.from_json(args.settings) if getattr(args, "settings", None) else Settings()
    try:
        if getattr(args, "km_split", None):
            s = s.replace(km_split=args.km_split)
        if getattr(args, "km_split_expression", None):
            s = s.replace(km_split_expression=args.km_split_expression)
    except ValueError as e:
        raise InputError(str(e)) from None
    if getattr(args, "case_label", None) or getattr(args, "reference_label", None):
        s = s.replace(case_label=args.case_label or s.case_label, reference_label=args.reference_label or
                      s.reference_label)
    return s


def _pairs(items, what) -> dict:
    out = {}
    for it in items or []:
        if "=" not in it:
            raise InputError(f"{what} expects NAME=VALUE, got {it!r}")
        k, v = it.split("=", 1)
        out[k.strip()] = v.strip()
    return out


def _dataset(args, event_ids=None) -> Dataset:
    return Dataset.from_dir(args.data, event_ids=event_ids, case=args.case or None,
                            reference=args.reference or None, columns=_pairs(args.column, "--column") or None,
                            na_values=args.na_value, keep=args.keep, keep_column=args.keep_column,
                            where=args.where, **_pairs(args.table, "--table"))


def _model(args) -> CoxModel:
    return CoxModel(expression=not args.no_expression, covariates=tuple(args.covariate or ()),
                    categorical=tuple(args.categorical or ()), strata=tuple(args.strata or ()),
                    baseline=_pairs(getattr(args, "baseline", None), "--baseline"))


def _detail(args, model: CoxModel):
    """The model band's model: the panel's model plus --detail-* clinical terms; True for --detail alone."""
    extra = [getattr(args, f"detail_{k}", None) for k in ("covariate", "categorical", "strata")]
    if any(extra):
        return model.with_clinical(*(tuple(x or ()) for x in extra))
    return True if getattr(args, "detail", False) else None


def _with_gene_events(args, events) -> list[str]:
    """The events plus the other events of their genes in the same q families: q values are computed within each
    gene, and HIT-index events form families of their own (so a PSI page never analyses them, nor a HIT page PSI)."""
    from .dataset import events_table
    from .events import quantity
    tables = _pairs(args.table, "--table")
    ev = events_table(tables.get("events") or tables.get("psi") or args.data, _pairs(args.column, "--column") or None)
    miss = [e for e in events if e not in ev.index]
    if miss:
        raise InputError(f"event(s) not in the events table: {', '.join(miss[:5])}")
    kind = ev.get("event_type", pd.Series("", index=ev.index)).map(quantity)
    fam = {(ev.at[e, "gene"], kind[e]) for e in events}
    return list(dict.fromkeys(list(events) + [e for e, g, k in zip(ev.index, ev.gene, kind) if (g, k) in fam]))


def _top(v) -> int:
    from .probe import parse_top
    return parse_top(v)


def _split(v) -> list[str]:
    return [x.strip() for x in str(v).split(";") if x.strip()]


# ============================================================================================ commands
def cmd_validate(args) -> int:
    ds = _dataset(args, args.event or None)
    print(f"samples {len(ds.samples)} · events {len(ds.psi)} · cohorts {len(ds.cohorts)} · endpoints "
          f"{', '.join(ds.endpoints) or 'none'} · pairs {len(ds.pairs)} ({ds.notes['pairs_source']})")
    print(f"groups: case = {ds.labels['case']}, reference = {ds.labels['reference']}"
          + (f"; not tested: {ds.notes['other_groups']}" if ds.notes["other_groups"] else ""))
    for what, where in ds.notes.get("sources", {}).items():
        print(f"{what}: from the {where}")
    sub = ds.notes.get("subset")
    if sub:
        print(f"subset: {sub['patients']} of {sub['of_patients']} patients ({sub['samples']} samples), selected by "
              f"{sub['by'].replace('`', '')}")
    if ds.notes["survival_rows_dropped"]:
        print("survival rows dropped (missing or non-positive time, or event not 0/1):",
              ", ".join(f"{k} {v}" for k, v in ds.notes["survival_rows_dropped"].items()))
    if ds.clinical is not None:
        print(f"clinical: {len(ds.clinical)} patients; columns {', '.join(map(str, ds.clinical.columns))}")
    with pd.option_context("display.width", 200, "display.max_columns", 50):
        print(ds.summary().to_string(index=False))
    missing_geo = ds.events[(ds.events.event_type == "") | (ds.events.variable == "")]
    if len(missing_geo):
        print(f"note: {len(missing_geo)} event(s) have no geometry; they can be analysed but not drawn")
    from .events import opt_in
    n_hit = int(ds.events.event_type.map(opt_in).sum())
    if n_hit:
        print(f"note: {n_hit} HIT-index event(s); analyse and probe leave them out unless --include-hit (a HIT event "
              "named with --event is always analysed)")
    print("ok")
    return 0


def cmd_analyse(args) -> int:
    from .analysis import analyse
    ds = _dataset(args, args.event or None)
    res = analyse(ds, events=args.event or None, cohorts=args.cohort or None, endpoints=args.endpoint or None,
                  settings=_settings(args), model=_model(args), include_hit=args.include_hit)
    for k, p in res.write(args.out).items():
        print(f"{k}: {p}")
    return 0


def _gtf(args):
    return args.gtf or os.environ.get("SPLICE_ASSAY_GTF") or None


def _proteins(args):
    """The protein cache folder: --proteins, else $SPLICE_ASSAY_PROTEINS; None for neither (or --no-proteins)."""
    if getattr(args, "no_proteins", False):
        return None
    return getattr(args, "proteins", None) or os.environ.get("SPLICE_ASSAY_PROTEINS") or None


def cmd_panel(args) -> int:
    """One assay page. Without --cohort the most promising cohorts are picked; without --endpoint OS is used;
    without --detail-* the model rows use age, sex and stage found in the clinical table (--no-detail: none)."""
    from .analysis import analyse
    from .clinical import auto_clinical
    from .plot import event_panel
    from .plot import expression_panel
    from .plot.panel_common import page_parts
    from .probe import combine, default_endpoint, pick_cohorts
    ids = _with_gene_events(args, args.event)                # q values are within the gene: analyse all its events
    ds = _dataset(args, ids)
    s, model = _settings(args), _model(args)
    ep = default_endpoint(ds, args.endpoint)
    detail = None if args.no_detail else _detail(args, model)
    if detail is None and not args.no_detail:
        ds, detail, found = auto_clinical(ds)
        if detail is not None:
            if args.baseline:
                detail = detail.with_clinical(baseline=_pairs(args.baseline, "--baseline"))
            print(f"model rows: Cox {detail.describe(ds.expression is not None)} (found {', '.join(found.values())})")
    cohorts = [c for c in (args.cohort or []) if c.lower() != "all"]
    results = analyse(ds, events=ids, endpoints=[ep], settings=s, model=model)
    adj = analyse(ds, events=ids, endpoints=[ep], settings=s, model=detail) if isinstance(detail, CoxModel) else None
    if not cohorts:
        cohorts = pick_cohorts(combine(results, adj, ep), args.event[0], args.top)
        if not cohorts:
            raise InputError(f"{args.event[0]}: no cohort has a survival test for {ep}")
        print(f"cohorts shown (most promising of {len(ds.cohorts)}): {', '.join(cohorts)}")
    parts = page_parts(cohorts, s.cohorts_per_page)         # many cohorts: several pages, each with the full forest
    gtf = _gtf(args)
    for i, chunk in enumerate(parts, 1):
        part = (i, len(parts)) if len(parts) > 1 else None
        p = event_panel(ds, args.event, chunk, ep, settings=s, model=model, gtf=gtf, highlight=args.highlight,
                        highlight_title=args.highlight_title, out_dir=args.out,
                        stem=f"{args.stem}_p{i}" if args.stem and part else args.stem, detail=detail,
                        layout=args.layout, results=results, proteins=_proteins(args), detail_results=adj, part=part)
        if i == 1:
            status = p.provenance["call"].get("protein_status") or {}
            for e, st in status.items():
                print(f"protein {e}: " + (st if st not in ("pair", "inclusion_only", "exclusion_only") else "drawn"))
        for k, v in p.paths.items():
            print(f"{k}: {v}")
    if not args.no_gex:                                     # the host gene's expression, on its own page
        gene = ds.events.at[args.event[0], "expression_gene"]
        xm = detail if isinstance(detail, CoxModel) else model
        for i, chunk in enumerate(parts, 1):
            part = (i, len(parts)) if len(parts) > 1 else None
            stem = (f"{args.stem}_expression" + (f"_p{i}" if part else "")) if args.stem else None
            x = expression_panel(ds, gene, chunk, ep, settings=s, model=xm, out_dir=args.out, stem=stem, part=part)
            if x is None:
                break
            print(f"expression: {x.paths.get('png', x.stem)}")
    return 0


def cmd_probe(args) -> int:
    """Every event of a gene (or every event) in every cohort: ranked, one assay page per event."""
    from .dataset import events_table
    from .probe import gene_events, probe
    event_ids = args.event or None
    if args.gene and not event_ids:                       # read only the gene's events
        tables = _pairs(args.table, "--table")
        src = tables.get("events") or tables.get("psi") or args.data
        tab = events_table(src, _pairs(args.column, "--column") or None).reset_index()
        event_ids = gene_events(tab, args.gene)
        if not event_ids:
            raise InputError(f"no event of gene(s) {', '.join(args.gene)} in the events table")
    ds = _dataset(args, event_ids)
    s = _settings(args)
    base = CoxModel(expression=not args.no_expression)
    if args.no_adjust:
        adjusted = None
    elif args.covariate or args.categorical or args.strata:
        adjusted = _model(args)
    else:
        adjusted = "auto"
    eps = args.endpoint or [None]                          # one probe per endpoint, each in its own folder
    if [str(e).lower() for e in eps] == ["all"]:
        eps = ds.endpoints
    for ep_arg in eps:
        ep = ep_arg or ("OS" if ds.survival is not None and "OS" in set(ds.survival.endpoint) else None)
        if args.out:
            out = Path(args.out) / ep if len(eps) > 1 else args.out
        else:
            out = f"probe_{'_'.join(args.gene or ['events'])}_{ep or 'endpoint'}"
        if len(eps) > 1:
            print(f"== {ep}")
        res = probe(ds, genes=args.gene, events=args.event, cohorts=args.cohort, endpoint=ep_arg, settings=s,
                    model=base, adjusted=adjusted, baseline=_pairs(args.baseline, "--baseline") or None,
                    gtf=_gtf(args), out_dir=out, top=args.top, max_pages=args.max_pages, proteins=_proteins(args),
                    gex=not args.no_gex, include_hit=args.include_hit,
                    call="splice-assay " + " ".join(shlex.quote(a) for a in args._argv))
        for k, v in res.paths.items():
            print(f"{k}: {v}")
        top = res.events[res.events.measurable].head(5) if len(res.events) else res.events
        for r in top.itertuples():
            if r.best_model == "KM":                      # no Cox model could be fitted
                print(f"  {r.rank}. {r.label} ({r.gene}): no Cox model (too few events); KM p<.05 in {r.km_p05} of "
                      f"{r.km_tested} cohorts; best {r.best_cohort} log-rank HR {r.best_logrank_hr:.2f}, "
                      f"p {r.best_p:.3g}")
            else:
                print(f"  {r.rank}. {r.label} ({r.gene}): adjusted p<.05 in {r.adj_cox_p05}, base p<.05 in "
                      f"{r.cox_p05} of {r.cohorts_cox_tested} cohorts; best {r.best_cohort} HR {r.best_hr_per_iqr:.2f}")
    return 0


def cmd_panels(args) -> int:
    """Many figures from a spec table: events and cohorts separated by ';', one endpoint per row, optional stem."""
    from .analysis import analyse
    from .plot import event_panel
    spec = read_table(args.spec)
    miss = [c for c in ("events", "cohorts", "endpoint") if c not in spec.columns]
    if miss:
        raise InputError(f"{args.spec}: missing column(s) {', '.join(miss)}")
    events = sorted({e for v in spec.events for e in _split(v)})
    s, model = _settings(args), _model(args)
    ids = _with_gene_events(args, events)                   # q values are within each gene
    ds = _dataset(args, ids)
    res = analyse(ds, events=ids, settings=s, model=model)
    gtf = _gtf(args)
    if gtf:                                                 # read the GTF once for every figure
        from .annotation import read_gtf
        from .events import geometry
        wins, genes = [], set()
        for e in events:
            row = ds.events.loc[e]
            lo, hi = geometry(e, row).span
            wins.append((row.chrom, lo - s.gtf_flank, hi + s.gtf_flank))
            genes |= {g for g in (row.gene, row.gene_id) if g}
        gtf = read_gtf(_gtf(args), wins, genes=sorted(genes))
    cache = None
    if _proteins(args):                                     # read the protein cache once for every figure
        from .protein import ProteinCache
        cache = ProteinCache.load(_proteins(args))
    for r in spec.itertuples():
        stem = getattr(r, "stem", None)
        p = event_panel(ds, _split(r.events), _split(r.cohorts), str(r.endpoint).strip(), settings=s, gtf=gtf,
                        highlight=args.highlight, highlight_title=args.highlight_title, out_dir=args.out,
                        stem=None if stem is None or pd.isna(stem) else str(stem), results=res,
                        detail=_detail(args, model), layout=args.layout, proteins=cache)
        print(p.paths.get("png", p.stem))
    if not args.no_gex:                                     # one expression page per gene and endpoint
        from .plot import expression_panel
        detail = _detail(args, model)
        xm = detail if isinstance(detail, CoxModel) else model
        todo = {}
        for r in spec.itertuples():
            for e in _split(r.events):
                key = (ds.events.at[e, "expression_gene"], str(r.endpoint).strip())
                todo.setdefault(key, [])
                todo[key] += [c for c in _split(r.cohorts) if c not in todo[key]]
        for (gene, ep), coh in todo.items():
            x = expression_panel(ds, gene, coh, ep, settings=s, model=xm, out_dir=args.out)
            if x is not None:
                print(x.paths.get("png", x.stem))
    return 0


def cmd_protein_cache(args) -> int:
    """Export a SpliceImpactR annotation cache to a protein cache folder (once; needs Rscript, any R)."""
    from .protein import build_cache
    for k, v in build_cache(args.annotation_cache, args.out, rscript=args.rscript).items():
        print(f"{k}: {v}")
    print(f"use it with --proteins {args.out} (or export SPLICE_ASSAY_PROTEINS={args.out})")
    return 0


def cmd_proteins(args) -> int:
    """The suggested protein change of events: printed, and with --out written as CSV."""
    from .dataset import events_table
    from .probe import gene_events
    from .protein import SHOWN, ProteinCache, protein_changes
    cols = _pairs(args.column, "--column") or None
    ev = events_table(args.events, cols)
    ids = list(args.event or []) + (gene_events(ev.reset_index(), args.gene) if args.gene else [])
    if not args.event and not args.gene:
        ids = list(ev.index)
    miss = [e for e in ids if e not in ev.index]
    if miss:
        raise InputError(f"event(s) not in the events table: {', '.join(miss[:5])}")
    if not ids:
        raise InputError("no events selected")
    folder = _proteins(args)
    if not folder:
        raise InputError("give the protein cache with --proteins (or set SPLICE_ASSAY_PROTEINS)")
    cache = ProteinCache.load(folder)

    class _Events:                                          # protein_changes needs only .events
        events = ev
    tab, _ = protein_changes(cache, _Events, list(dict.fromkeys(ids)))
    for r in tab.itertuples():
        text = getattr(r, "summary", "") if r.protein_status in SHOWN else f"no suggestion ({r.protein_status})"
        print(f"{r.event_id}: {text}")
    if args.out:
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        tab.to_csv(out, index=False, lineterminator="\n")
        print(f"table: {out}")
    return 0


def cmd_cox(args) -> int:
    """The full Cox model of one event in one cohort: printed, and with --out also as CSV and a figure."""
    from .plot import cox_model_figure, model_terms
    ds = _dataset(args, [args.event])
    s, model = _settings(args), _model(args)
    row, terms = model_terms(ds, args.event, args.cohort, args.endpoint, model, s)
    use_expr = ds.expression is not None and model.expression
    print(f"{args.event} · {args.cohort} · {args.endpoint}: Cox {model.describe(use_expr)}")
    print(f"{int(row.cox_n)} patients, {int(row.cox_events)} events"
          + (f"; {int(row.cox_n_dropped)} left out for missing values" if row.get("cox_n_dropped", 0) else "")
          + (f"; {row.cox_notes}" if isinstance(row.get("cox_notes"), str) and row.cox_notes else ""))
    show = terms[["term", "level", "unit", "coef", "se", "hr", "ci_low", "ci_high", "p"]].copy()
    with pd.option_context("display.width", 200, "display.max_columns", 20, "display.float_format", "{:.4g}".format):
        print(show.to_string(index=False))
    if args.out:
        p = cox_model_figure(ds, args.event, args.cohort, args.endpoint, model=model, settings=s, out_dir=args.out)
        for k, v in p.paths.items():
            print(f"{k}: {v}")
    return 0


def cmd_example(args) -> int:
    from . import example
    from .analysis import analyse
    from .plot import cox_model_figure, event_panel, expression_panel
    out = Path(args.out)
    paths = example.write(out, separate=args.separate)
    ds = Dataset.from_dir(out / "data")
    res = analyse(ds)
    res.write(out / "results")
    clinical = CoxModel().with_clinical(("age", "sex", "stage"), baseline={"stage": "I", "sex": "female"})
    figs = [expression_panel(ds, "SYN1", ["COH1", "COH2"], "OS", model=clinical, out_dir=out / "figures"),
            event_panel(ds, "SYN1:SE:1", ["COH1", "COH2"], "OS", gtf=paths["gtf"], out_dir=out / "figures",
                        results=res, detail=clinical, proteins=paths["proteins"]),
            event_panel(ds, ["SYN1:SE:1", "SYN1:RI:1"], ["COH1"], "DSS", gtf=paths["gtf"], out_dir=out / "figures",
                        results=res),
            event_panel(ds, "SYN2:MXE:1", ["COH3"], "OS", gtf=paths["gtf"], out_dir=out / "figures", results=res),
            # res leaves out the HIT index (as analyse does by default), so this page analyses SYN3's HIT events
            event_panel(ds, "SYN3:HIT:0002", ["COH1"], "OS", gtf=paths["gtf"], out_dir=out / "figures"),
            cox_model_figure(ds, "SYN1:SE:1", "COH1", "OS", model=clinical, out_dir=out / "figures")]
    print(f"data:     {out / 'data'}")
    print(f"hitindex: {paths['hitindex']} (HITindex matrices of SYN3, for import-hitindex)")
    print(f"proteins: {paths['proteins']} (a synthetic protein cache)")
    print(f"results:  {out / 'results'}")
    for f in figs:
        print(f"figure:   {f.paths['png']}")
    return 0


def cmd_import_rmats(args) -> int:
    from .rmats import import_rmats
    names1 = args.names1.split(",") if args.names1 else None
    names2 = args.names2.split(",") if args.names2 else None
    types = [t.strip().upper() for t in args.types.split(",")] if args.types else ["SE", "RI", "A3SS", "A5SS", "MXE"]
    events, psi = import_rmats(args.rmats_dir, args.b1, args.b2, counting=args.counting, event_types=types,
                               prefix=args.prefix, names1=names1, names2=names2, mxe_psi_exon=args.mxe_psi_exon)
    _write_event_tables(Path(args.out), events, psi, not args.separate, args.append)
    return 0


def cmd_import_hitindex(args) -> int:
    from .hitindex import import_hitindex
    gtf = _gtf(args)
    if not gtf:
        raise InputError("import-hitindex needs --gtf (or $SPLICE_ASSAY_GTF): HITindex IDs carry no strand")
    events, psi = import_hitindex(args.matrices, gtf, genes=args.gene, prefix=args.prefix)
    n = events.event_type.value_counts().to_dict()
    print("imported: " + ", ".join(f"{n[k]} {k}" for k in ("AFE", "ALE", "HIT") if k in n))
    _write_event_tables(Path(args.out), events, psi, not args.separate, args.append)
    return 0


def _write_event_tables(out: Path, events: pd.DataFrame, psi: pd.DataFrame, one_table: bool, append: bool) -> None:
    """Write events and PSI as events.csv + psi.csv (long), or one psi.csv with the event columns. With `append`,
    add them to the tables already in `out` (e.g. HITindex events beside rMATS ones), in the layout found there:
    events.csv beside a long or wide psi.csv, or the event columns inside psi.csv (wide or long). The rows already
    there are rewritten exactly as read."""
    from .dataset import find_tables
    out.mkdir(parents=True, exist_ok=True)
    ev_f, psi_f = out / "events.csv", out / "psi.csv"
    found = find_tables(out)
    if append and found.get("psi") is not None:
        other = [found[t].name for t in ("psi", "events") if t in found and found[t].name != f"{t}.csv"]
        if other:
            raise InputError(f"--append adds to CSV tables; {out} has {', '.join(other)}")
        old = pd.read_csv(psi_f, dtype=str, keep_default_na=False, low_memory=False)      # kept exactly as written
        in_psi = "events" not in found                              # the event columns are inside psi.csv
        old_events = (old[[c for c in old.columns if c in EVENT_COLUMNS]] if in_psi
                      else pd.read_csv(ev_f, dtype=str, keep_default_na=False, low_memory=False))
        clash = sorted(set(old_events.get("event_id", [])) & set(events.event_id))
        if clash:
            raise InputError(f"--append: event(s) already in {out}: {', '.join(clash[:5])}")
        if "sample_id" in old.columns:                              # long: one row per event and sample
            add = psi.merge(events, on="event_id", how="left") if in_psi else psi
        else:                                                       # wide: one column per sample
            add = psi.pivot(index="event_id", columns="sample_id", values="psi")
            add = add.reindex(index=events.event_id, columns=list(dict.fromkeys(psi.sample_id))).reset_index()
            if in_psi:
                add = events.merge(add, on="event_id")
        table = pd.concat([old, add], ignore_index=True)
        table.to_csv(psi_f, index=False, lineterminator="\n")
        if not in_psi:
            ev_all = pd.concat([old_events, events], ignore_index=True)
            ev_all.to_csv(ev_f, index=False, lineterminator="\n")
            print(f"events: {ev_f} ({len(ev_all)}, {len(events)} added)")
        print(f"psi:    {psi_f} ({len(events)} events added to the "
              + ("long" if "sample_id" in old.columns else "wide") + " table)")
        return
    if one_table:                                       # the event columns, then one PSI column per sample
        wide = psi.pivot(index="event_id", columns="sample_id", values="psi")
        wide = wide.reindex(index=events.event_id, columns=list(dict.fromkeys(psi.sample_id))).reset_index()
        table = events.merge(wide, on="event_id")
        table.to_csv(psi_f, index=False, lineterminator="\n")
        print(f"psi:    {psi_f} ({len(table)} events with their event columns)")
    else:
        events.to_csv(ev_f, index=False, lineterminator="\n")
        psi.to_csv(psi_f, index=False, lineterminator="\n")
        print(f"events: {ev_f} ({len(events)})")
        print(f"psi:    {psi_f} ({psi.sample_id.nunique()} samples)")
    if "samples" not in found:
        print("next: add samples.csv (sample_id, patient_id, cohort, group; survival and clinical columns may go in it "
              "too) to the same folder")


def cmd_gtf_subset(args) -> int:
    from .annotation import subset_gtf
    from .dataset import events_table
    from .events import geometry
    ev = events_table(args.events)
    wins, genes = [], set()
    for e, row in ev.iterrows():
        if not row.variable:
            continue
        lo, hi = geometry(e, row).span
        wins.append((row.chrom, lo - args.flank, hi + args.flank))
        genes |= {g for g in (row.gene, row.gene_id) if g}
    n = subset_gtf(args.gtf, args.out, wins, sorted(genes))
    print(f"{n} GTF lines -> {args.out}")
    return 0


# ============================================================================================ parser
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="splice-assay", description=__doc__.split("\n")[0])
    p.add_argument("--version", action="version", version=f"splice-assay {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    def data_args(sp):
        sp.add_argument("data", help="folder with samples, psi, events[, survival, pairs, expression, clinical] tables")
        sp.add_argument("--table", action="append", metavar="NAME=PATH", help="use this file for one table")
        sp.add_argument("--column", action="append", metavar="LOGICAL=ACTUAL",
                        help="your column name for a logical one, e.g. cohort=cancer, sample_id=File.ID")
        sp.add_argument("--case", action="append",
                        help="group value of the case samples (default: the samples table's role column, else tumour)")
        sp.add_argument("--reference", action="append",
                        help="group value of the reference samples (default: the role column, else normal)")
        sp.add_argument("--case-label", help="display name of the case group")
        sp.add_argument("--reference-label", help="display name of the reference group")
        sp.add_argument("--na-value", action="append", metavar="CODE",
                        help='a code to read as missing in every table, e.g. missing or "[Not Available]"')
        sp.add_argument("--settings", help="JSON file of setting overrides")
        sp.add_argument("--km-split", metavar="median|mean|VALUE",
                        help="KM split of the event values, PSI or HIT index (default median; Settings km_split)")
        sp.add_argument("--km-split-expression", metavar="median|mean|VALUE",
                        help="KM split of host-gene expression (default median; Settings km_split_expression)")
        sp.add_argument("--keep", metavar="FILE",
                        help="only the patients listed in this table (all their samples, e.g. one subtype); IDs may be "
                             "patient or sample IDs, or barcodes that start with them")
        sp.add_argument("--keep-column", help="the column of --keep holding the IDs (default: the first)")
        sp.add_argument("--where", action="append", metavar="COLUMN=VALUE[,VALUE]",
                        help="only patients whose clinical row or samples have one of these values in COLUMN "
                             "(case-insensitive; repeat for more conditions, all must hold), e.g. Subtype=Her2")

    def model_args(sp):
        sp.add_argument("--covariate", action="append", help="clinical column to adjust for (repeat for more)")
        sp.add_argument("--categorical", action="append", help="treat this covariate as categorical")
        sp.add_argument("--strata", action="append", help="clinical column used as Cox strata")
        sp.add_argument("--baseline", action="append", metavar="COLUMN=LEVEL",
                        help="reference level of a categorical covariate, e.g. stage=I")
        sp.add_argument("--no-expression", action="store_true", help="do not adjust for host-gene expression")

    def detail_args(sp):
        sp.add_argument("--detail", action="store_true",
                        help="add the full Cox model of each cohort shown (the model band)")
        sp.add_argument("--detail-covariate", action="append",
                        help="model band only: clinical column added to the model (repeat for more)")
        sp.add_argument("--detail-categorical", action="append", help="model band only: categorical covariate")
        sp.add_argument("--detail-strata", action="append", help="model band only: strata column")
        sp.add_argument("--layout", default="auto", choices=["auto", "side", "stacked"],
                        help="with models: side = band below; stacked = one row per cohort with its model, forest "
                             "below (auto: stacked from three cohorts)")

    def layout_args(sp):
        sp.add_argument("--separate", action="store_true",
                        help="write events.csv and a long psi.csv instead of the default: one psi.csv holding the "
                             "event columns and one value column per sample")
        sp.add_argument("--one-table", action="store_true", help=argparse.SUPPRESS)   # the default now; kept for old
        # commands

    def hit_arg(sp):
        sp.add_argument("--include-hit", action="store_true",
                        help="also analyse HIT-index events (left out by default: the HIT index covers every exon, a "
                             "far larger set; a HIT event named with --event is always analysed)")

    sp = sub.add_parser("validate", help="check the input tables and print per-cohort counts")
    data_args(sp)
    sp.add_argument("--event", action="append", help="check only these events")
    sp.set_defaults(func=cmd_validate)

    sp = sub.add_parser("analyse", aliases=["analyze"], help="group and survival statistics as CSV")
    data_args(sp)
    model_args(sp)
    sp.add_argument("--out", required=True)
    sp.add_argument("--event", action="append")
    sp.add_argument("--cohort", action="append")
    sp.add_argument("--endpoint", action="append")
    hit_arg(sp)
    sp.set_defaults(func=cmd_analyse)

    def gex_args(sp):
        sp.add_argument("--no-gex", action="store_true",
                        help="no host-gene expression page (default: one per gene when an expression table is given)")

    def protein_args(sp):
        sp.add_argument("--proteins", metavar="DIR",
                        help="protein cache (from `splice-assay protein-cache`) for the suggested protein change "
                             "(default: $SPLICE_ASSAY_PROTEINS)")
        sp.add_argument("--no-proteins", action="store_true", help="ignore $SPLICE_ASSAY_PROTEINS")

    def fig_args(sp):
        sp.add_argument("--gtf", help="GTF (.gtf or .gtf.gz) for the gene track (default: $SPLICE_ASSAY_GTF)")
        sp.add_argument("--highlight", help="table of cells to mark: event_id, cohort, label[, endpoint, colour]")
        sp.add_argument("--highlight-title", help="header of the mark column, e.g. Tier")
        sp.add_argument("--out", required=True)

    sp = sub.add_parser("panel", help="one assay page: one or two events; cohorts, endpoint and model rows default "
                                      "to the most promising cohorts, OS and age + sex + stage")
    data_args(sp)
    model_args(sp)
    sp.add_argument("--event", action="append", required=True)
    sp.add_argument("--cohort", action="append", help="cohorts to show (default: the most promising; 'all' = same)")
    sp.add_argument("--endpoint", help="default OS")
    sp.add_argument("--top", type=_top, default=3,
                    help="how many cohorts to pick when none are given: a number or 'all' (default 3)")
    sp.add_argument("--no-detail", action="store_true", help="no model rows")
    sp.add_argument("--stem", help="output file name without extension")
    fig_args(sp)
    detail_args(sp)
    protein_args(sp)
    gex_args(sp)
    sp.set_defaults(func=cmd_panel)

    sp = sub.add_parser("probe", help="explore: every event of a gene (or every event) in every cohort, ranked, with "
                                      "one assay page per event")
    data_args(sp)
    model_args(sp)
    sp.add_argument("--gene", action="append", help="probe the events of this gene (repeat for more; default: all)")
    sp.add_argument("--event", action="append", help="probe only these events")
    sp.add_argument("--cohort", action="append", help="cohorts to probe (default: all)")
    sp.add_argument("--endpoint", action="append",
                    help="endpoint to probe, e.g. OS (repeat for more, or 'all'; each gets its own folder; default OS)")
    sp.add_argument("--gtf", help="GTF for the gene track (default: $SPLICE_ASSAY_GTF)")
    sp.add_argument("--top", type=_top, default=3,
                    help="cohorts shown per page: a number or 'all' (default 3; the forest always shows every cohort)")
    sp.add_argument("--max-pages", type=int, default=30, help="pages for the best-ranked events (default 30)")
    hit_arg(sp)
    sp.add_argument("--no-adjust", action="store_true", help="no adjusted model (default: age + sex + stage found)")
    sp.add_argument("--out", help="output folder (default probe_<gene>_<endpoint>)")
    protein_args(sp)
    gex_args(sp)
    sp.set_defaults(func=cmd_probe)

    sp = sub.add_parser("panels", help="many figures from a spec table (events, cohorts, endpoint[, stem])")
    data_args(sp)
    model_args(sp)
    sp.add_argument("--spec", required=True)
    fig_args(sp)
    detail_args(sp)
    protein_args(sp)
    gex_args(sp)
    sp.set_defaults(func=cmd_panels)

    sp = sub.add_parser("proteins", help="the suggested protein change of events (from a protein cache)")
    sp.add_argument("events", help="an events table (or a psi table with the event columns), or a data folder")
    sp.add_argument("--gene", action="append", help="events of this gene (repeat for more)")
    sp.add_argument("--event", action="append", help="this event (repeat for more; default: every event)")
    sp.add_argument("--column", action="append", metavar="LOGICAL=ACTUAL", help="your column name for a logical one")
    sp.add_argument("--out", help="CSV file for the table")
    protein_args(sp)
    sp.set_defaults(func=cmd_proteins)

    sp = sub.add_parser("protein-cache", help="a protein cache from a SpliceImpactR annotation cache (once; needs "
                                              "Rscript, any R)")
    sp.add_argument("annotation_cache", help="folder of SpliceImpactR objects: <name>.gtf.rds, <name>_sequences.rds "
                                             "and protein-feature .rds files (docs/annotation-cache.md)")
    sp.add_argument("--out", required=True, help="folder for the protein cache")
    sp.add_argument("--rscript", default="Rscript", help="Rscript to use (default: Rscript on the PATH)")
    sp.set_defaults(func=cmd_protein_cache)

    sp = sub.add_parser("cox", help="the full Cox model of one event in one cohort (table, and a figure with --out)")
    data_args(sp)
    model_args(sp)
    sp.add_argument("--event", required=True)
    sp.add_argument("--cohort", required=True)
    sp.add_argument("--endpoint", required=True)
    sp.add_argument("--out", help="folder for the CSV and figure")
    sp.set_defaults(func=cmd_cox)

    sp = sub.add_parser("example", help="write a synthetic dataset and example outputs")
    sp.add_argument("out")
    sp.add_argument("--separate", action="store_true",
                    help="write the data as separate tables (samples, psi, events, survival, clinical, expression) "
                         "instead of the default two (samples, psi) plus expression")
    sp.set_defaults(func=cmd_example)

    sp = sub.add_parser("import-rmats", help="psi.csv (with the event columns) from an rMATS output folder")
    sp.add_argument("rmats_dir")
    sp.add_argument("--b1", required=True, help="rMATS b1.txt (BAM paths of group 1)")
    sp.add_argument("--b2", help="rMATS b2.txt")
    sp.add_argument("--names1", help="comma-separated sample names instead of BAM names (group 1)")
    sp.add_argument("--names2", help="comma-separated sample names (group 2)")
    sp.add_argument("--counting", default="JCEC", choices=["JC", "JCEC"])
    sp.add_argument("--types", help="comma-separated event types (default SE,RI,A3SS,A5SS,MXE)")
    sp.add_argument("--prefix", default="", help="prefix for event IDs")
    sp.add_argument("--mxe-psi-exon", default="transcript_upstream", choices=["transcript_upstream", "first_listed"])
    sp.add_argument("--append", action="store_true", help="add to the events and psi tables already in --out")
    layout_args(sp)
    sp.add_argument("--out", required=True)
    sp.set_defaults(func=cmd_import_rmats)

    sp = sub.add_parser("import-hitindex", help="psi.csv (with the event columns) from HITindex matrices (AFE, ALE, "
                                                "HIT index)")
    sp.add_argument("matrices", nargs="+", help="afe/ale/hit matrices: rows '<gene_id>;<chrom>:<start>-<end>;<kind>' "
                                               "(1-based), one column per sample (CSV or TSV, .gz allowed)")
    sp.add_argument("--gtf", help="GTF giving each gene's strand and symbol (default: $SPLICE_ASSAY_GTF)")
    sp.add_argument("--gene", action="append", help="import only this gene (symbol or Ensembl ID; repeat for more)")
    sp.add_argument("--prefix", default="", help="prefix for event IDs")
    layout_args(sp)
    sp.add_argument("--append", action="store_true",
                    help="add to the events and psi tables already in --out (e.g. beside imported rMATS events)")
    sp.add_argument("--out", required=True)
    sp.set_defaults(func=cmd_import_hitindex)

    sp = sub.add_parser("gtf-subset", help="a small GTF with the records around the events")
    sp.add_argument("gtf")
    sp.add_argument("--events", required=True)
    sp.add_argument("--out", required=True)
    sp.add_argument("--flank", type=int, default=20000)
    sp.set_defaults(func=cmd_gtf_subset)
    return p


def main(argv=None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    args = build_parser().parse_args(argv)
    args._argv = argv
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("default")
            return args.func(args)
    except InputError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
