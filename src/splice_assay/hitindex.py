"""Import HITindex output: alternative first and last exons (AFE, ALE) and the HIT index of exons.

Each input is one matrix per kind: one row per exon, one column per sample (CSV or TSV, optionally gzipped). The
first column holds the exon's ID, '<gene_id>;<chrom>:<start>-<end>;<kind>' in HITindex's 1-based inclusive
coordinates, with kind afe, ale or hit (a matrix assembled from HITindex's per-sample outputs looks like this). The
values are AFE or ALE PSI (0..1: the exon's use among the gene's first or last exons) or the HIT index (-1..1: -1 for
an exon used as a first exon, 1 as a last exon, about 0 for an internal exon). The HIT index is not PSI and has its
own rules in the analysis (analysis.event_settings).

HITindex IDs carry no strand and no gene symbol, so a GTF gives both.

    events   event_id '<prefix><SYMBOL>:<TYPE>:<nnnn>', numbered per gene and type in transcript order (5' to 3');
             variable = the exon (0-based half-open); for AFE and ALE, constant = the gene's other first or last exons
             in the matrix (drawn grey); label '<TYPE>:<nnnn>'; source_id = the HITindex ID. When one symbol names
             several genes (gene IDs, or the X and Y copies of a pseudoautosomal gene), each is ordered and drawn on
             its own and the numbering continues across them
    psi      event_id, sample_id, psi (the AFE or ALE PSI, or the HIT index)
"""
from __future__ import annotations

import gzip
import io
import re
from pathlib import Path

import pandas as pd

from .annotation import ATTR
from .dataset import InputError
from .events import format_intervals

KINDS = {"afe": "AFE", "ale": "ALE", "hit": "HIT"}
_ID = re.compile(r"^(?P<gene>[^;]+);(?P<chrom>[^:;]+):(?P<start>\d+)-(?P<end>\d+);(?P<kind>[A-Za-z]+)$")


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def parse_id(text) -> tuple[str, str, int, int, str]:
    """'ENSG00000204371.13;chr6:31879759-31880123;afe' -> (unversioned gene_id, chrom, start0, end, 'AFE'):
    the exon in 0-based half-open coordinates."""
    m = _ID.match(str(text).strip())
    if not m or m["kind"].lower() not in KINDS:
        raise InputError(f"HITindex ID {str(text)[:80]!r}: expected <gene_id>;<chrom>:<start>-<end>;<afe|ale|hit>")
    start, end = int(m["start"]), int(m["end"])
    if end < start:
        raise InputError(f"HITindex ID {text!r}: end before start")
    return m["gene"].split(".")[0], m["chrom"], start - 1, end, KINDS[m["kind"].lower()]


def gtf_genes(gtf) -> pd.DataFrame:
    """Each gene of a GTF: symbol (gene_name), chrom and strand, indexed by the unversioned gene_id."""
    rows = []
    with _open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if len(f) < 9 or f[2] != "gene":
                continue
            a = dict(ATTR.findall(f[8]))
            rows.append((a.get("gene_id", "").split(".")[0], a.get("gene_name", ""), f[0], f[6]))
    g = pd.DataFrame(rows, columns=["gene_id", "gene_name", "chrom", "strand"])
    return g.drop_duplicates("gene_id").set_index("gene_id")


def read_matrix(path, gene_ids=None) -> pd.DataFrame:
    """The rows of one HITindex matrix (all, or those of the unversioned `gene_ids`), indexed by HITindex ID with one
    float column per sample. Large matrices are streamed, so filtering by gene reads little into memory."""
    want = None if gene_ids is None else set(gene_ids)
    with _open(path) as fh:
        header = fh.readline()
        sep = "\t" if "\t" in header else ","
        kept = [header]
        for line in fh:                  # the ID may be quoted (R's write.csv), as read_csv would read it
            if want is None or line[:line.find(sep)].strip().strip("\"'").split(";", 1)[0].split(".")[0] in want:
                kept.append(line)
    m = pd.read_csv(io.StringIO("".join(kept)), sep=sep, index_col=0, dtype={0: str}, low_memory=False)
    m.index = m.index.astype(str).str.strip()
    m.columns = [str(c).strip() for c in m.columns]
    return m.apply(pd.to_numeric, errors="coerce")


def import_hitindex(matrices, gtf, genes=None, prefix: str = "") -> tuple[pd.DataFrame, pd.DataFrame]:
    """Events and long PSI from HITindex matrices (any of afe, ale and hit; each file may mix kinds). `genes`: gene
    symbols or Ensembl IDs to keep (default all)."""
    info = gtf_genes(gtf)
    want = None
    if genes:
        names = {str(g).split(".")[0] for g in genes}
        want = set(info.index[info.index.isin(names) | info.gene_name.isin(names)])
        if not want:
            raise InputError(f"none of {', '.join(sorted(names))} is a gene of the GTF")
    parts = [read_matrix(p, want) for p in ([matrices] if isinstance(matrices, (str, Path)) else matrices)]
    m = pd.concat(parts) if parts else pd.DataFrame()
    if not len(m):
        raise InputError("no HITindex row" + (" for these genes" if genes else "") + " in the matrices")
    dup = m.index[m.index.duplicated()]
    if len(dup):
        raise InputError(f"HITindex ID listed twice: {dup[0]}")
    rows = []
    for hid in m.index:
        gid, chrom, a, b, kind = parse_id(hid)
        g = info.loc[gid] if gid in info.index else None
        rows.append(dict(source_id=hid, gene_id=gid, gene=(g.gene_name if g is not None and g.gene_name else gid),
                         chrom=chrom, strand=g.strand if g is not None else "", event_type=kind, start=a, end=b))
    ex = pd.DataFrame(rows)
    missing = ex.loc[ex.strand.eq(""), "gene_id"].unique()
    if len(missing):
        import warnings
        warnings.warn(f"{len(missing)} gene(s) not in the GTF (e.g. {missing[0]}): no strand, so their events are not "
                      "drawn", stacklevel=2)
    events = []
    for (gene, kind), d in ex.groupby(["gene", "event_type"], sort=True):
        n = 0                            # numbered per symbol (unique IDs); ordered, and other exons, per gene
        for _, g in d.groupby(["gene_id", "chrom"], sort=True):    # one symbol can name several genes, or PAR copies
            minus = g.strand.iloc[0] == "-"
            g = g.sort_values(["end", "start"], ascending=False) if minus else g.sort_values(["start", "end"])
            ivs = list(zip(g.start, g.end))
            for r in g.itertuples():
                n += 1
                others = sorted(iv for iv in ivs if iv != (r.start, r.end)) if kind in ("AFE", "ALE") else []
                events.append(dict(event_id=f"{prefix}{gene}:{kind}:{n:04d}", gene=gene, gene_id=r.gene_id,
                                   chrom=r.chrom, strand=r.strand, event_type=kind,
                                   constant=format_intervals(others), variable=f"{r.start}-{r.end}",
                                   label=f"{kind}:{n:04d}", source_id=r.source_id))
    events = pd.DataFrame(events)
    ids = dict(zip(events.source_id, events.event_id))
    psi = m.rename(index=ids).rename_axis("event_id").reset_index().melt(
        id_vars="event_id", var_name="sample_id", value_name="psi").dropna(subset=["psi"])
    order = {e: i for i, e in enumerate(events.event_id)}
    psi = psi.sort_values(["event_id", "sample_id"], key=lambda c: c.map(order) if c.name == "event_id" else c)
    return events, psi.reset_index(drop=True)
