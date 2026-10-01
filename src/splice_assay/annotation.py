"""Gene models from a GTF: read the records near some events, collapse the host gene's transcripts into exon blocks,
and find the small RNA genes nested in it.

Collapsed gene model (a drawing rule): the host gene's transcripts with >= 2 exons, minus excluded transcript types
(default retained_intron). An exon is keyed by its splice sites; a transcript-terminal end is free and drawn at the
median of the transcripts that share the key. Keys used by >= min_frac of these transcripts (and at least 2) are kept,
except exons that span an intron used by >= min_frac of them (intron-retaining variants). Overlapping kept exons are
merged into one block.
"""
from __future__ import annotations

import gzip
import re
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .config import Settings

ATTR = re.compile(r'(\S+) "([^"]*)"')
GTF_COLUMNS = ["chrom", "feature", "start", "end", "strand", "gene_id", "gene_name", "gene_type", "transcript_id",
               "transcript_type"]


def _aliases(chrom: str) -> set[str]:
    c = str(chrom)
    base = c[3:] if c.startswith("chr") else c
    out = {c, base, "chr" + base}
    if base in ("M", "MT"):
        out |= {"M", "MT", "chrM", "chrMT"}
    return out


def read_gtf(path, windows, genes=(), features=("gene", "transcript", "exon")) -> pd.DataFrame:
    """GTF records overlapping any (chrom, start0, end) window, plus every record of the named genes on those
    chromosomes. `genes` holds gene names or gene IDs (with or without version). Coordinates are returned 0-based,
    half-open; chromosome names are returned as spelled in `windows` ('chr1' and '1' match each other)."""
    wins, spell = {}, {}
    for chrom, s, e in windows:
        for a in _aliases(chrom):
            wins.setdefault(a, []).append((int(s), int(e)))
            spell[a] = str(chrom)
    tags = []
    for g in genes:
        g = str(g)
        tags += [f'gene_id "{g}"', f'gene_id "{g}.', f'gene_name "{g}"']
    opener = gzip.open if str(path).endswith(".gz") else open
    rows = []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            chrom = line[: line.find("\t")]
            w = wins.get(chrom)
            if w is None:
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] not in features:
                continue
            s1, e1 = int(f[3]), int(f[4])
            if not (any(s1 <= b and e1 > a for a, b in w) or any(t in f[8] for t in tags)):
                continue
            a = dict(ATTR.findall(f[8]))
            rows.append((spell[chrom], f[2], s1 - 1, e1, f[6], a.get("gene_id", "").split(".")[0],
                         a.get("gene_name", ""), a.get("gene_type", a.get("gene_biotype", "")),
                         a.get("transcript_id", ""), a.get("transcript_type", a.get("transcript_biotype", ""))))
    return pd.DataFrame(rows, columns=GTF_COLUMNS)


def collapsed_model(exons: pd.DataFrame, min_frac: float = 0.10, exclude=("retained_intron",)):
    """Exon blocks [(start, end, n_keys)] of the collapsed model, the number of transcripts used and the threshold."""
    ex = exons[~exons.transcript_type.isin(exclude)]
    keys, introns, ntx = {}, {}, 0
    for _, d in ex.groupby("transcript_id"):
        d = d.sort_values("start")
        s, e = d.start.to_numpy(int), d.end.to_numpy(int)
        m = len(s)
        if m < 2:
            continue
        ntx += 1
        for i in range(m):
            k = (None if i == 0 else int(s[i]), None if i == m - 1 else int(e[i]))
            keys.setdefault(k, []).append((int(s[i]), int(e[i])))
        for i in range(m - 1):
            introns[(int(e[i]), int(s[i + 1]))] = introns.get((int(e[i]), int(s[i + 1])), 0) + 1
    thr = max(2.0, min_frac * ntx)
    frequent = [k for k, c in introns.items() if c >= thr]
    kept = []
    for k, v in keys.items():
        if len(v) < thr:
            continue
        a = k[0] if k[0] is not None else int(round(np.median([x[0] for x in v])))
        b = k[1] if k[1] is not None else int(round(np.median([x[1] for x in v])))
        if any(a <= i0 and b >= i1 for i0, i1 in frequent):
            continue
        kept.append((a, b, len(v)))
    kept.sort()
    blocks = []
    for a, b, c in kept:
        if blocks and a <= blocks[-1][1]:
            blocks[-1] = [blocks[-1][0], max(blocks[-1][1], b), blocks[-1][2] + c]
        else:
            blocks.append([a, b, c])
    return [tuple(b) for b in blocks], ntx, thr


@dataclass
class GeneModel:
    gene: str
    chrom: str
    strand: str
    blocks: list            # [(start, end, n_keys)]
    n_transcripts: int
    min_transcripts: float
    nested: pd.DataFrame    # gene_name, start, end (small RNA genes inside the drawn span)


def gene_model(gtf: pd.DataFrame, gene: str, chrom: str, span, gene_id: str = "",
               settings: Settings | None = None) -> GeneModel | None:
    """The collapsed model of the host gene (matched by gene_id when given, else by name, on `chrom`, overlapping
    `span`) and its nested small RNA genes; None (with a warning) when the GTF has no such gene."""
    s = settings or Settings()
    gid = str(gene_id).split(".")[0] if gene_id else ""
    on = gtf[gtf.chrom.isin(_aliases(chrom))]
    host = on[on.gene_id.eq(gid)] if gid else on[on.gene_name.eq(gene)]
    rec = host[host.feature.eq("gene")]
    rec = rec[(rec.start < span[1]) & (rec.end > span[0])]
    if rec.empty:
        warnings.warn(f"GTF: no gene {gid or gene} on {chrom} overlapping the events; the gene track is omitted",
                      stacklevel=2)
        return None
    rec = rec.iloc[0]
    host = host[host.gene_id.eq(rec.gene_id)]
    blocks, ntx, thr = collapsed_model(host[host.feature.eq("exon")], s.gene_model_min_frac,
                                       s.exclude_transcript_types)
    if not blocks:
        warnings.warn(f"GTF: gene {rec.gene_name or rec.gene_id} has no multi-exon transcript; the gene track is "
                      "omitted", stacklevel=2)
        return None
    sn = on[on.feature.eq("gene") & on.gene_type.isin(s.nested_biotypes) & (on.start < rec.end) & (on.end > rec.start)]
    sn = sn.groupby("gene_name", as_index=False).agg(start=("start", "min"), end=("end", "max"))
    lo, hi = blocks[0][0], blocks[-1][1]
    sn = sn[(sn.end > lo) & (sn.start < hi)].sort_values("start").reset_index(drop=True)
    return GeneModel(gene=rec.gene_name or gene, chrom=chrom, strand=rec.strand, blocks=blocks, n_transcripts=ntx,
                     min_transcripts=thr, nested=sn)


def subset_gtf(path, out_path, windows, genes=()) -> int:
    """Write the GTF lines overlapping `windows` or naming `genes` to a smaller (optionally gzipped) GTF."""
    wins = {}
    for chrom, s, e in windows:
        for a in _aliases(chrom):
            wins.setdefault(a, []).append((int(s), int(e)))
    tags = [t for g in genes for t in (f'gene_id "{g}"', f'gene_id "{g}.', f'gene_name "{g}"')]
    opener = gzip.open if str(path).endswith(".gz") else open
    writer = gzip.open if str(out_path).endswith(".gz") else open
    n = 0
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with opener(path, "rt") as fh, writer(out_path, "wt") as out:
        for line in fh:
            if line.startswith("#"):
                continue
            chrom = line[: line.find("\t")]
            w = wins.get(chrom)
            if w is None:
                continue
            f = line.split("\t", 8)
            if len(f) < 9:
                continue
            s1, e1 = int(f[3]), int(f[4])
            if any(s1 <= b and e1 > a for a, b in w) or any(t in f[8] for t in tags):
                out.write(line)
                n += 1
    return n
