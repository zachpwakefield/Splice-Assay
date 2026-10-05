"""Protein consequences: the likely transcript of each splice form, and how the protein changes.

This is a suggestion read from annotation; nothing here is measured. It ports two parts of SpliceImpactR.

Matching (get_ranked_pairs, protocol 1). Each event has an inclusion form (INC, the form PSI counts) and an exclusion
form (EXC). Each form is a list of required intervals and a list of forbidden ones (1-based, closed):
    SE    INC  upstream, cassette, downstream               EXC  upstream, downstream; forbids the cassette
    RI    INC  upstream, intron, downstream (one exon)      EXC  upstream, downstream; forbids the intron
    A3SS  INC  long exon, flanking exon                     EXC  short exon, flanking exon; forbids the extension
    A5SS  as A3SS
    MXE   INC  upstream, PSI exon, downstream;              EXC  upstream, other exon, downstream;
          forbids the other exon                                 forbids the PSI exon
A transcript represents a form when all three of these hold:
  - every required interval lies inside one of its exons;
  - no exon touches a forbidden interval;
  - the splice structure agrees exactly. That means consecutive exons whose inner boundaries are the event's
    junctions. A retained intron must lie inside one exon. An A3SS/A5SS site must be an exon boundary with a
    neighbouring exon.
SpliceImpactR's rMATS converter leaves out the flanking exon of A3SS/A5SS events. It is supplied here, so that junction
is checked too.

rMATS often names flanking exons that are not the annotated neighbours (a third exon can lie between them). When no
transcript matches a form exactly, the form falls back to SpliceImpactR's definition from the event's own exons:
  - SE, MXE: the exon, with exactly its boundaries, as an internal exon;
  - RI: the intron inside one exon;
  - A3SS/A5SS: the alternative splice site as an exon boundary.
The SE and RI exclusion forms have no such definition. The `inc_context` / `exc_context` columns say which applied
(junction_chain, partial_event, retained_intron, splice_site_only).

Pairs (INC transcript, EXC transcript) of distinct transcripts are ranked by, in order:
  1. how many of the two are protein-coding;
  2. the worse transcript support tier (TSL 1-3, 4-5, unknown);
  3. the evidence used;
  4. the similarity of their coding sequence outside the event;
  5. the similarity of their whole coding sequence;
  6. the worse TSL, then the summed TSL;
  7. the transcript IDs.
Similarity is the Jaccard overlap of coding genomic positions, which is SpliceImpactR's fallback evidence.
SpliceImpactR aligns the sequences when it has them, so a near-tie can rank differently. When only one form has a
transcript, that transcript is chosen by protein-coding status, then TSL. `n_tied` counts the pairs that tie with
the chosen one.

The forms are named by event type: inclusion / skipping (SE), retention / splicing (RI), long / short form
(A3SS, A5SS), PSI exon / other exon (MXE).

The change (`effect`: kind, frame, a short label and a sentence):
  - When both chosen transcripts are coding and differ only at the event, the two proteins are compared directly.
    "Only at the event" means one of two things:
      - they share every coding position outside it;
      - they share at least 95% of those positions, the event keeps the frame, and every differing residue lies
        within the event's residues (two residues of slack for split codons).
  - Otherwise the event's residues and frame are read from the coding coordinates of the transcript that carries the
    event (or, when only the other form is coding, of the place where the event would sit in it). The frame is the
    length difference of the two forms modulo 3. A stop codon inside an unannotated insertion cannot be seen, and
    the sentence says so.
  - A pair sharing less than half of its coding positions outside the event is too different to compare. Only the
    transcript carrying the event is drawn.
  - Feature databases: InterPro, Pfam, CDD, MobiDB-lite, TMHMM, SignalP and ELM. A feature on one transcript's
    event exons is "only with" that form when the other protein has no feature of the same name whose genomic span
    overlaps it (get_domains).
  - Named features: InterPro families, superfamilies and entries whose pieces span 90% or more of the protein are
    not named. Overlapping InterPro entries are shown as one (preferring a domain or repeat).
There is no suggestion when no chosen transcript is protein-coding.

The protein cache is three tables: exons with coding coordinates, proteins, and features. `build_cache` makes them
from a SpliceImpactR annotation cache; this needs Rscript (any R, no packages), once. Other annotation can be supplied in
the same format (docs/proteins.md).
"""
from __future__ import annotations

import math
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, field
from functools import cached_property
from pathlib import Path

import numpy as np
import pandas as pd

from .dataset import InputError

EXON_COLS = ["transcript_id", "transcript_name", "transcript_type", "transcript_support_level", "gene_id",
             "gene_name", "chr", "strand", "start", "end", "exon_id", "cds_gen_start", "cds_gen_stop",
             "cds_rel_start", "cds_rel_stop"]
PROTEIN_COLS = ["transcript_id", "protein_seq"]
FEATURE_COLS = ["transcript_id", "database", "name", "start", "stop", "g_start", "g_end"]
READABLE = {"mobidblite": "disordered region", "tmhmm": "transmembrane helix", "signalp": "signal peptide"}
SUMMARY_DBS = ("interpro", "mobidblite", "tmhmm", "signalp", "elm")
SHOWN = ("pair", "inclusion_only", "exclusion_only")
_TSL = re.compile(r"^\s*([1-5])(?:\s.*)?$")

R_EXPORT = r"""
# base R only, so any R that can read the cache works (no data.table needed)
a <- commandArgs(TRUE); d <- a[1]; out <- a[2]; pre <- a[3]
rd <- function(f) { x <- readRDS(f); class(x) <- "data.frame"; x }
wr <- function(x, f) {
  con <- gzfile(file.path(out, f), "w")
  write.table(x, con, sep = "\t", quote = TRUE, qmethod = "double", row.names = FALSE, na = "")
  close(con)
}
g <- rd(file.path(d, paste0(pre, ".gtf.rds")))
cols <- c("transcript_id","transcript_name","transcript_type","transcript_support_level","gene_id","gene_name","chr",
          "strand","start","end","exon_id","exon_number","protein_id","cds_gen_start","cds_gen_stop","cds_rel_start",
          "cds_rel_stop")
wr(g[g$type == "exon", cols], "exons.tsv.gz")
s <- rd(file.path(d, paste0(pre, "_sequences.rds")))
s <- s[!is.na(s$protein_seq) & nzchar(s$protein_seq), c("transcript_id", "protein_id", "protein_seq")]
wr(s[!duplicated(s), ], "proteins.tsv.gz")
dbs <- c("interpro","pfam","cdd","elm","mobidblite","signalp","tmhmm")
f <- do.call(rbind, lapply(dbs, function(x) {
  p <- file.path(d, paste0(x, ".rds"))
  if (!file.exists(p)) return(NULL)
  y <- rd(p)
  m <- regmatches(y$name, regexec(";[^;:]+:([0-9]+)-([0-9]+)$", y$name))
  gs <- vapply(m, function(v) if (length(v) == 3) as.numeric(v[2]) else NA_real_, 0)
  ge <- vapply(m, function(v) if (length(v) == 3) as.numeric(v[3]) else NA_real_, 0)
  data.frame(transcript_id = y$ensembl_transcript_id, database = as.character(y$database),
             feature_id = as.character(y$feature_id), name = as.character(y$clean_name),
             start = as.integer(y$start), stop = as.integer(y$stop),
             g_start = as.integer(pmin(gs, ge)), g_end = as.integer(pmax(gs, ge)), stringsAsFactors = FALSE)
}))
if (is.null(f)) f <- data.frame(transcript_id = character(), database = character(), feature_id = character(),
                                name = character(), start = integer(), stop = integer(), g_start = integer(),
                                g_end = integer())
wr(f[!duplicated(f), ], "features.tsv.gz")
"""


# ============================================================================================ the cache
def _strip_version(s: pd.Series) -> pd.Series:
    return s.astype(str).str.replace(r"\.\d+(_PAR_Y)?$", "", regex=True)


@dataclass
class ProteinCache:
    exons: pd.DataFrame        # one row per exon (1-based, closed), with coding coordinates where coding
    proteins: pd.DataFrame     # transcript_id, protein_seq (protein_id optional)
    features: pd.DataFrame     # transcript_id, database, name, start, stop (residues), g_start, g_end (genomic)
    source: str = ""

    @classmethod
    def load(cls, folder) -> "ProteinCache":
        folder = Path(folder).expanduser()
        if not folder.is_dir():
            raise InputError(f"{folder}: no such protein cache folder (make one with `splice-assay protein-cache`)")

        def read(name, need):
            for ext in (".parquet", ".tsv.gz", ".tsv", ".csv"):
                p = folder / f"{name}{ext}"
                if p.exists():
                    df = (pd.read_parquet(p) if ext == ".parquet" else
                          pd.read_csv(p, sep="," if ext == ".csv" else "\t", low_memory=False))
                    miss = [c for c in need if c not in df.columns]
                    if miss:
                        raise InputError(f"{p}: missing column(s) {', '.join(miss)}")
                    return df
            raise InputError(f"{folder}: no {name} table (.parquet, .tsv.gz, .tsv or .csv)")

        ex = read("exons", EXON_COLS)
        for c in ("gene_id", "transcript_id"):
            ex[c] = _strip_version(ex[c])
        ex["chr"] = ex["chr"].astype(str)
        for c in ("start", "end", "cds_gen_start", "cds_gen_stop", "cds_rel_start", "cds_rel_stop"):
            ex[c] = pd.to_numeric(ex[c], errors="coerce")
        pr = read("proteins", PROTEIN_COLS)
        pr["transcript_id"] = _strip_version(pr["transcript_id"])
        ft = read("features", FEATURE_COLS)
        ft["transcript_id"] = _strip_version(ft["transcript_id"])
        ft["name"] = ft["name"].astype(str)
        return cls(ex, pr.drop_duplicates("transcript_id"), ft, str(folder))

    @cached_property
    def _by_gene(self) -> dict:
        return self.exons.groupby("gene_id", sort=False).indices

    @cached_property
    def _by_name(self) -> dict:
        return self.exons.groupby("gene_name", sort=False).indices

    @cached_property
    def _features_by_tx(self) -> dict:
        return self.features.groupby("transcript_id", sort=False).indices

    @cached_property
    def _seq(self) -> dict:
        return dict(zip(self.proteins.transcript_id, self.proteins.protein_seq))

    def sequence(self, transcript: str) -> str:
        s = self._seq.get(transcript, "")
        return s if isinstance(s, str) else ""

    def features_of(self, transcript: str) -> pd.DataFrame:
        """A protein's features, with a readable label and `named`: whether the feature is worth naming (not an
        InterPro family or superfamily, and not an entry whose pieces together span 90% or more of the protein)."""
        f = self.features.iloc[self._features_by_tx.get(transcript, [])].copy()
        f["label"] = [READABLE.get(d, n) for d, n in zip(f.database, f.name)]
        n = len(self.sequence(transcript)) or (int(f.stop.max()) if len(f) else 0)
        low = f.name.str.lower()
        family = f.database.eq("interpro") & (low.str.contains("superfamily") | low.str.endswith(" family")
                                              | low.str.contains(" family,"))
        span = f.groupby(["database", "name"]).stop.transform("max") - f.groupby(["database", "name"]).start.transform(
            "min") + 1 if len(f) else f.stop                # all pieces of an entry together
        whole = f.database.isin(["interpro", "pfam", "cdd"]) & (span >= 0.9 * max(n, 1))
        f["named"] = f.database.isin(SUMMARY_DBS) & ~family & ~whole
        return f.reset_index(drop=True)

    def transcripts_for(self, row, span: tuple[int, int]) -> pd.DataFrame:
        """Exons of the transcripts that could represent an event. These are the gene's transcripts (found by ID,
        else by name) on the event's chromosome and strand. When the gene is not found, every transcript on that
        strand that overlaps the event is used."""
        ex = self.exons
        chrom = str(row.get("chrom", "") or "")
        strand = str(row.get("strand", "") or "")
        chroms = {chrom, "chr" + chrom.removeprefix("chr"), chrom.removeprefix("chr")} - {"", "chr"}

        def here(df):
            keep = df.chr.isin(chroms) if chroms else pd.Series(True, index=df.index)
            if strand in ("+", "-"):
                keep &= df.strand.eq(strand)
            return df[keep]

        gid = str(row.get("gene_id", "") or "").split(".")[0]
        sub = here(ex.iloc[self._by_gene.get(gid, [])]) if gid and gid != "nan" else ex.iloc[:0]
        if sub.empty and row.get("gene"):
            sub = here(ex.iloc[self._by_name.get(str(row["gene"]), [])])
        if sub.empty:
            near = here(ex[(ex.start <= span[1]) & (ex.end >= span[0])])
            sub = ex[ex.transcript_id.isin(near.transcript_id.unique())]
        return sub


def cache_prefix(folder) -> str:
    """The annotation of a SpliceImpactR export folder: '<name>' for the one <name>.gtf.rds that has its
    <name>_sequences.rds beside it (e.g. human_gencode_v45; docs/annotation-cache.md shows how to write them)."""
    src = Path(folder).expanduser()
    names = sorted(p.name[:-len(".gtf.rds")] for p in src.glob("*.gtf.rds")
                   if (src / f"{p.name[:-len('.gtf.rds')]}_sequences.rds").exists())
    if not names:
        raise InputError(f"{src}: no <name>.gtf.rds with its <name>_sequences.rds, e.g. human_gencode_v45.gtf.rds "
                         "(docs/annotation-cache.md shows how to save them from SpliceImpactR)")
    if len(names) > 1:
        raise InputError(f"{src}: {len(names)} annotations ({', '.join(names)}); keep one per folder")
    return names[0]


def build_cache(annotation_cache, out_dir, rscript: str = "Rscript") -> dict[str, Path]:
    """Export a SpliceImpactR annotation folder to a protein cache folder. The folder holds <name>.gtf.rds and
    <name>_sequences.rds (e.g. human_gencode_v45) plus the protein-feature .rds files. The output is Parquet when
    pyarrow is installed, else gzipped TSV. Any R (3.5 or newer) works; no R package is needed."""
    found = shutil.which(rscript)
    if found is None:
        raise InputError(f"{rscript} was not found: Rscript is needed once, to read the SpliceImpactR .rds files "
                         "(give its path with --rscript)")
    src = Path(annotation_cache).expanduser()
    pre = cache_prefix(src)
    out = Path(out_dir).expanduser()
    out.mkdir(parents=True, exist_ok=True)
    paths = {}
    with tempfile.TemporaryDirectory() as tmp:
        script = Path(tmp) / "export.R"
        script.write_text(R_EXPORT)
        r = subprocess.run([found, str(script), str(src), tmp, pre], capture_output=True, text=True)
        if r.returncode:
            raise InputError(f"the R export with {found} failed (another R can be given with --rscript): "
                             f"{r.stderr.strip()[-800:]}")
        for name in ("exons", "proteins", "features"):
            df = pd.read_csv(Path(tmp) / f"{name}.tsv.gz", sep="\t", low_memory=False)
            try:
                p = out / f"{name}.parquet"
                df.to_parquet(p, index=False)
            except ImportError:
                p = out / f"{name}.tsv.gz"
                df.to_csv(p, sep="\t", index=False)
            paths[name] = p
    (out / "SOURCE.txt").write_text(f"Exported from {src} ({pre}) by splice-assay protein-cache.\n")
    return paths


# ============================================================================================ intervals
Iv = tuple[int, int]


def _merge(ivs) -> list[Iv]:
    out: list[list[int]] = []
    for a, b in sorted((int(a), int(b)) for a, b in ivs):
        if out and a <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return [(a, b) for a, b in out]


def _subtract(ivs, mask) -> list[Iv]:
    out = []
    for a, b in _merge(ivs):
        pieces = [(a, b)]
        for m0, m1 in mask:
            pieces = [q for p0, p1 in pieces for q in ((p0, min(p1, m0 - 1)), (max(p0, m1 + 1), p1)) if q[0] <= q[1]]
        out += pieces
    return out


def _length(ivs) -> int:
    return sum(b - a + 1 for a, b in ivs)


def _intersect(a, b) -> list[Iv]:
    return _merge([(max(x0, y0), min(x1, y1)) for x0, x1 in a for y0, y1 in b if max(x0, y0) <= min(x1, y1)])


def _jaccard(a, b) -> float:
    a, b = _merge(a), _merge(b)
    inter = _length(_intersect(a, b))
    total = _length(a) + _length(b) - inter
    return inter / total if total else float("nan")


def _sym_diff(a, b) -> list[Iv]:
    a, b = _merge(a), _merge(b)
    return _subtract(a, b) + _subtract(b, a)


# ============================================================================================ forms and matching
def _closed(iv) -> Iv:
    """0-based half-open (the events table) -> 1-based closed (annotation)."""
    return int(iv[0]) + 1, int(iv[1])


@dataclass(frozen=True)
class Form:
    inc: tuple[Iv, ...]        # required intervals, sorted
    exc: tuple[Iv, ...]        # forbidden intervals


def forms(geom, flanks: bool = True) -> dict[str, Form]:
    """The two forms, per SpliceImpactR's rules (see the module docstring). With flanks=False, the definitions from
    the event's own exons; the SE and RI exclusion forms have none. {} for other event types."""
    t = geom.event_type
    c = sorted(_closed(x) for x in geom.constant)
    v = [_closed(x) for x in geom.variable]
    o = [_closed(x) for x in geom.other]
    f = lambda inc, exc=(): Form(tuple(sorted(inc)), tuple(sorted(exc)))  # noqa: E731
    if t in ("SE", "RI") and len(c) == 2 and len(v) == 1:
        if not flanks:
            return {"INC": f([v[0]])}
        return {"INC": f([c[0], v[0], c[1]]), "EXC": f([c[0], c[1]], [v[0]])}
    if t in ("A3SS", "A5SS") and len(c) == 2 and len(v) == 1:
        short = next((x for x in c if x[0] == v[0][1] + 1 or x[1] == v[0][0] - 1), None)
        if short is None:
            return {}
        flank = c[1] if short == c[0] else c[0]
        long_ = (min(short[0], v[0][0]), max(short[1], v[0][1]))
        if not flanks:
            return {"INC": f([long_]), "EXC": f([short], [v[0]])}
        return {"INC": f([long_, flank]), "EXC": f([short, flank], [v[0]])}
    if t == "MXE" and len(c) == 2 and len(v) == 1 and len(o) == 1:
        if not flanks:
            return {"INC": f([v[0]], [o[0]]), "EXC": f([o[0]], [v[0]])}
        return {"INC": f([c[0], v[0], c[1]], [o[0]]), "EXC": f([c[0], o[0], c[1]], [v[0]])}
    return {}


def tsl_rank(x) -> int:
    m = _TSL.match(str(x))
    return int(m.group(1)) if m else 6


def _tier(r: int) -> int:
    return 0 if r <= 3 else (1 if r <= 5 else 2)


@dataclass
class Check:
    ok: bool
    reason: str
    context: str = ""
    exons: tuple[int, ...] = ()        # rows (of the transcript's start-sorted exons) holding each required interval


def compatible(tx: pd.DataFrame, form: Form, event_type: str, which: str, strand: str) -> Check:
    """Whether one transcript (its exons sorted by start) represents a form; a port of .si_rank_compatible."""
    inc, exc = form.inc, form.exc
    if not inc:
        return Check(False, "no_required_interval")
    if tx.empty:
        return Check(False, "no_exons")
    st, en = tx.start.to_numpy(), tx.end.to_numpy()
    for a, b in exc:
        if ((st <= b) & (en >= a)).any():
            return Check(False, "forbidden_interval_overlap")
    ix = []
    for lo_, hi_ in inc:
        w = np.flatnonzero((st <= lo_) & (en >= hi_))
        if not len(w):
            return Check(False, "required_interval_not_covered")
        ix.append(int(w[0]))
    ni, ne = len(inc), len(tx)
    lo = np.array([a for a, _ in inc])
    hi = np.array([b for _, b in inc])

    def chain_ok():
        if ni < 2 or np.any(np.diff(ix) != 1):
            return False
        return bool(np.all(en[ix[:-1]] == hi[:-1]) and np.all(st[ix[1:]] == lo[1:]) and np.all(lo[1:] > hi[:-1] + 1))

    context = "interval_only"
    if event_type in ("SE", "MXE"):
        expected = 2 if (event_type == "SE" and which == "EXC") else 3
        if ni == expected:
            if not chain_ok():
                return Check(False, "splice_junction_chain_mismatch")
            context = "junction_chain"
        elif ni == 1 and not (event_type == "SE" and which == "EXC"):
            j = ix[0]
            if j == 0 or j == ne - 1 or st[j] != lo[0] or en[j] != hi[0]:
                return Check(False, "internal_exon_boundary_mismatch")
            context = "partial_event"
        else:
            return Check(False, "insufficient_event_context")
    elif event_type == "RI":
        if which == "INC":
            contiguous = ni == 1 or bool(np.all(lo[1:] == hi[:-1] + 1))
            if not contiguous or len(set(ix)) != 1:
                return Check(False, "retained_intron_not_in_one_exon")
            context = "retained_intron" if ni == 3 else "partial_event"
        else:
            if ni != 2 or not chain_ok():
                return Check(False, "retained_intron_spliced_chain_mismatch")
            context = "junction_chain"
    elif event_type in ("A5SS", "A3SS"):
        if ni not in (1, 2):
            return Check(False, "insufficient_event_context")
        variable_end = (event_type == "A5SS" and strand == "+") or (event_type == "A3SS" and strand == "-")
        anchor = 0 if (ni == 1 or variable_end) else ni - 1
        j = ix[anchor]
        ok = (en[j] == hi[anchor] and j < ne - 1) if variable_end else (st[j] == lo[anchor] and j > 0)
        if not ok:
            return Check(False, "alternative_splice_boundary_mismatch")
        if ni == 2 and not chain_ok():
            return Check(False, "splice_junction_chain_mismatch")
        context = "junction_chain" if ni == 2 else "splice_site_only"
    return Check(True, "eligible", context, tuple(ix))


def _cds(tx: pd.DataFrame) -> list[Iv]:
    ok = (tx.cds_gen_start.notna() & tx.cds_gen_stop.notna() & (tx.cds_gen_start >= tx.start)
          & (tx.cds_gen_stop <= tx.end) & (tx.cds_gen_stop >= tx.cds_gen_start))
    return _merge(zip(tx.cds_gen_start[ok].astype(int), tx.cds_gen_stop[ok].astype(int)))


PAIR_KEYS = ["coding_count", "tsl_tier", "evidence_rank", "context_similarity", "orf_similarity", "tsl_worst",
             "tsl_sum"]


@dataclass
class Match:
    status: str                       # 'pair', 'inclusion_only', 'exclusion_only', or why there is no suggestion
    inc: str = ""                     # the chosen transcripts ('' when that form has none)
    exc: str = ""
    candidates: pd.DataFrame = field(default_factory=pd.DataFrame)     # every transcript x form: eligible, reason
    pairs: pd.DataFrame = field(default_factory=pd.DataFrame)          # ranked pairs (or single transcripts)
    checks: dict = field(default_factory=dict)                         # (form, transcript) -> the Check used
    mask: list = field(default_factory=list)                           # the event's own intervals (1-based)


def match_forms(exons: pd.DataFrame, geom, strand: str) -> Match:
    """Choose the transcripts of an event's two forms (see the module docstring)."""
    full, part = forms(geom), forms(geom, flanks=False)
    if not full:
        return Match("unsupported_event_type")
    ex = exons.sort_values(["transcript_id", "start", "end", "exon_id"], kind="stable")
    by_tx = {t: d.reset_index(drop=True) for t, d in ex.groupby("transcript_id", sort=True)}
    meta = ex.drop_duplicates("transcript_id").set_index("transcript_id")
    rows, checks, used = [], {}, {}
    for which in ("INC", "EXC"):
        for level, fm in (("exact", full), ("event_exons", part)):
            if which not in fm:
                continue
            cs = {t: compatible(d, fm[which], geom.event_type, which, strand) for t, d in by_tx.items()}
            for t, c in cs.items():
                rows.append(dict(form=which, level=level, transcript_id=t, eligible=c.ok, reason=c.reason,
                                 context=c.context, coding=meta.at[t, "transcript_type"] == "protein_coding",
                                 tsl=tsl_rank(meta.at[t, "transcript_support_level"])))
            if any(c.ok for c in cs.values()):
                used[which] = fm[which]
                checks.update({(which, t): c for t, c in cs.items() if c.ok})
                break
    cand = pd.DataFrame(rows, columns=["form", "level", "transcript_id", "eligible", "reason", "context", "coding",
                                       "tsl"])
    mask = _merge([*_sym_diff(full["INC"].inc, full["EXC"].inc), *full["INC"].exc, *full["EXC"].exc])
    out = Match("no_transcript", candidates=cand, checks=checks, mask=mask)
    ok = cand[cand.eligible].drop_duplicates(["form", "transcript_id"])
    ca, cb = ok[ok.form.eq("INC")], ok[ok.form.eq("EXC")]
    if ca.empty and cb.empty:
        return out
    if ca.empty or cb.empty:                                # one form only: the best single transcript
        one = (ca if len(ca) else cb).assign(tier=lambda d: d.tsl.map(_tier))
        one = one.sort_values(["coding", "tier", "tsl", "transcript_id"], ascending=[False, True, True, True],
                              kind="stable").reset_index(drop=True)
        best = one.iloc[0]
        one["tied_best"] = (one.coding.eq(best.coding) & one.tier.eq(best.tier) & one.tsl.eq(best.tsl))
        out.pairs = one
        if not best.coding:
            out.status = "noncoding"
            return out
        if len(ca):
            out.status, out.inc = "inclusion_only", str(best.transcript_id)
        else:
            out.status, out.exc = "exclusion_only", str(best.transcript_id)
        return out
    cds = {t: _cds(by_tx[t]) for t in set(ca.transcript_id) | set(cb.transcript_id)}
    exn = {t: _merge(zip(by_tx[t].start.astype(int), by_tx[t].end.astype(int))) for t in cds}
    pairs = []
    for a in ca.itertuples():
        for b in cb.itertuples():
            if a.transcript_id == b.transcript_id:
                continue
            if (compatible(by_tx[a.transcript_id], used["EXC"], geom.event_type, "EXC", strand).ok
                    and compatible(by_tx[b.transcript_id], used["INC"], geom.event_type, "INC", strand).ok):
                continue                                     # neither transcript tells the forms apart
            ra, rb = cds[a.transcript_id], cds[b.transcript_id]
            evidence, rank = "cds_coordinates", 1
            if not ra or not rb:
                ra, rb = exn[a.transcript_id], exn[b.transcript_id]
                evidence, rank = "exon_coordinates", 2
            context = _jaccard(_subtract(ra, mask), _subtract(rb, mask))
            if math.isnan(context):
                evidence, rank = "no_outside_event_context", 3
            pairs.append(dict(inc=a.transcript_id, exc=b.transcript_id, coding_count=int(a.coding) + int(b.coding),
                              tsl_tier=max(_tier(a.tsl), _tier(b.tsl)), evidence_rank=rank, evidence=evidence,
                              context_similarity=context, orf_similarity=_jaccard(ra, rb),
                              tsl_worst=max(a.tsl, b.tsl), tsl_sum=a.tsl + b.tsl))
    if not pairs:
        out.status = "no_distinct_pair"
        return out
    p = pd.DataFrame(pairs).sort_values(PAIR_KEYS + ["inc", "exc"],
                                        ascending=[False, True, True, False, False, True, True, True, True],
                                        na_position="last", kind="stable").reset_index(drop=True)
    best = p.iloc[0]
    p["tied_best"] = np.all([np.isclose(p[k].astype(float), float(best[k]), rtol=0, atol=1e-12, equal_nan=True)
                             for k in PAIR_KEYS], axis=0)
    out.pairs = p
    if best.coding_count == 0:
        out.status = "noncoding"
        return out
    out.status, out.inc, out.exc = "pair", str(best.inc), str(best.exc)
    return out


# ============================================================================================ the change
def _aa_spans(tx: pd.DataFrame, region: list[Iv]) -> list[Iv]:
    """Residue spans that genomic intervals encode in a transcript (residue = ceil(coding nt / 3))."""
    rel = []
    for e in tx.itertuples():
        if pd.isna(e.cds_gen_start) or pd.isna(e.cds_rel_start):
            continue
        for a, b in region:
            lo, hi = max(a, int(e.cds_gen_start)), min(b, int(e.cds_gen_stop))
            if lo > hi:
                continue
            if e.strand == "+":
                r0, r1 = e.cds_rel_start + (lo - e.cds_gen_start), e.cds_rel_start + (hi - e.cds_gen_start)
            else:
                r0, r1 = e.cds_rel_start + (e.cds_gen_stop - hi), e.cds_rel_start + (e.cds_gen_stop - lo)
            rel.append((int(r0), int(r1)))
    return [(math.ceil(a / 3), math.ceil(b / 3)) for a, b in _merge(rel)]


def _exon_aa(e) -> Iv | None:
    """An exon's residue span, as get_exon_features (inclusive)."""
    if pd.isna(e.cds_rel_start) or pd.isna(e.cds_rel_stop):
        return None
    return math.ceil(e.cds_rel_start / 3), math.ceil(e.cds_rel_stop / 3)


def _overlapping(f: pd.DataFrame, spans) -> pd.DataFrame:
    if f.empty or not spans:
        return f.iloc[:0]
    hit = np.zeros(len(f), dtype=bool)
    for a, b in spans:
        hit |= (f.start.to_numpy() <= b) & (f.stop.to_numpy() >= a)
    return f[hit]


def _only(on_exons: pd.DataFrame, other: pd.DataFrame) -> pd.DataFrame:
    """Features on one transcript's event exons such that the other protein has no feature of the same database and
    name whose genomic span overlaps them (SpliceImpactR .diff_domains)."""
    if on_exons.empty or other.empty:
        return on_exons
    keep = []
    for r in on_exons.itertuples():
        o = other[other.database.eq(r.database) & other.name.eq(r.name)]
        keep.append(not ((o.g_start <= r.g_end) & (o.g_end >= r.g_start)).any())
    return on_exons[np.array(keep, dtype=bool)]


def _compare(p: str, q: str, mxe: bool = False, subject: str = "inclusion") -> dict:
    """How protein p (the PSI form) differs from q (the other form), when the two differ only by the event."""
    n = 0
    while n < min(len(p), len(q)) and p[n] == q[n]:
        n += 1
    m = 0
    while m < min(len(p), len(q)) - n and p[-1 - m] == q[-1 - m]:
        m += 1
    xp, xq = p[n:len(p) - m], q[n:len(q) - m]
    span = f"residues {n + 1}–{n + len(xp)}"
    def junction(k):
        return f"; {k} residue{'s' if k > 1 else ''} at the junction change{'' if k > 1 else 's'}"
    if p == q:
        kind, text = "identical", "the same protein"
    elif mxe and m >= 5:
        kind, text = "exchange", (f"the PSI form has {len(xp)} aa ({span}) where the other form has {len(xq)} aa; "
                                  f"the rest of the protein is the same")
    elif not xq and m:
        kind, text = "insertion", f"{subject} adds {len(xp)} aa in frame ({span})"
    elif not xp and m:
        kind, text = "deletion", (f"{subject} removes {len(xq)} aa in frame (residues {n + 1}–{n + len(xq)} "
                                  "without it)")
    elif len(xp) == len(xq) <= 3 and not mxe:
        k = len(xp)
        kind, text = "substitution", (f"{subject} changes residue{'s' if k > 1 else ''} {n + 1}"
                                      + (f"–{n + k}" if k > 1 else "") + f" ({xq}→{xp})")
    elif m >= 5 and len(xq) <= 2 < len(xp):
        kind, text = "insertion", f"{subject} adds {len(xp) - len(xq)} aa in frame ({span}{junction(len(xq))})"
    elif m >= 5 and len(xp) <= 2 < len(xq):
        kind, text = "deletion", (f"{subject} removes {len(xq) - len(xp)} aa in frame (residues {n + 1}–"
                                  f"{n + len(xq)} without it{junction(len(xp))})")
    elif m >= 5:
        kind, text = "replacement", f"{subject} replaces {len(xq)} aa with {len(xp)} aa ({span})"
    elif n == 0:
        kind, text = "n_terminus", "the proteins start differently (another start codon)"
    else:
        kind, text = "c_terminus", f"the proteins differ from residue {n + 1} to the end (frameshift or new stop codon)"
    return dict(kind=kind, text=text, shared_prefix=n, shared_suffix=m, inc_block=len(xp), exc_block=len(xq))


def labels(f: pd.DataFrame | None) -> list[str]:
    """Readable names of the features worth naming (InterPro domains and sites, disorder, TM helices, signal
    peptides, ELM motifs; see ProteinCache.features_of)."""
    if f is None or f.empty:
        return []
    return sorted(set(f[f.named].label.astype(str)))


def domain_clusters(f: pd.DataFrame) -> list[tuple[str, int, int]]:
    """Named InterPro entries, overlapping ones (by half the shorter) grouped. Each group is shown as one member:
    the longest whose name says domain or repeat, else the longest; its own span is used."""
    if f is None or f.empty:
        return []
    d = f[f.named & f.database.eq("interpro")].sort_values(["start", "stop", "label"])
    groups: list[list] = []
    for r in d.itertuples():
        a, b = int(r.start), int(r.stop)
        for g in groups:
            if min(b, g[1]) - max(a, g[0]) + 1 >= 0.5 * min(b - a + 1, g[1] - g[0] + 1):
                g[0], g[1] = min(a, g[0]), max(b, g[1])
                g[2].append((r.label, a, b))
                break
        else:
            groups.append([a, b, [(r.label, a, b)]])

    def pick(members):
        low = lambda x: x.lower()  # noqa: E731
        return sorted(members, key=lambda m: (not ("domain" in low(m[0]) or "repeat" in low(m[0])),
                                              -(m[2] - m[1]), m[0]))[0]
    return [pick(g[2]) for g in groups]


def _named_at(f: pd.DataFrame, spans) -> list[str]:
    """Feature names at some residues: merged InterPro spans, and the other named features. A feature must cover
    at least 3 of the residues (all of them when there are fewer), so one that only touches the edge is not named."""
    if f is None or f.empty or not spans:
        return []

    def at(a, b):
        return any(min(b, y) - max(a, x) + 1 >= min(3, y - x + 1) for x, y in spans)
    out = [lab for lab, a, b in domain_clusters(f) if at(a, b)]
    other = f[f.named & ~f.database.eq("interpro")]
    out += sorted({r.label for r in other.itertuples() if at(int(r.start), int(r.stop))})
    return list(dict.fromkeys(out))


def _list(xs, k):
    """Feature names joined with ' / ' (InterPro names contain commas)."""
    xs = list(xs)
    return " / ".join(xs[:k]) + (f" and {len(xs) - k} more" if len(xs) > k else "")


def _pct(x: float) -> str:
    return f"{x:.1%}" if 0.99 <= x < 1 else f"{x:.0%}"


BIOTYPE_TEXT = {"protein_coding_CDS_not_defined": "has no annotated coding sequence",
                "nonsense_mediated_decay": "is annotated as a nonsense-mediated decay target",
                "non_stop_decay": "is annotated as a non-stop decay target",
                "retained_intron": "is annotated as a retained-intron (noncoding) transcript"}
NOUN = {"SE": "exon", "RI": "retained intron", "A3SS": "extension", "A5SS": "extension", "MXE": "exon"}
FORM_NAMES = {"SE": ("inclusion", "skipping"), "RI": ("retention", "splicing"), "A3SS": ("long form", "short form"),
              "A5SS": ("long form", "short form"), "MXE": ("PSI exon", "other exon")}


def form_names(event_type: str) -> tuple[str, str]:
    """What the PSI form and the other form are called for an event type (inclusion / skipping, ...)."""
    return FORM_NAMES.get(event_type, ("PSI form", "other form"))
WEAK_PAIR = 0.5                    # a pair sharing less of its coding sequence outside the event is not drawn


@dataclass
class Isoform:
    form: str                          # INC or EXC
    transcript: str
    name: str
    biotype: str
    tsl: int
    context: str                       # how it matched (junction_chain, partial_event, ...)
    coding: bool
    protein: str                       # sequence ('' when not coding)
    event_aa: list                     # residue spans the event's own intervals encode
    insert_after: int | None           # residue after which the absent region would sit (when it has none)
    utr: str                           # "5′ UTR" / "3′ UTR" when the event lies outside the coding sequence
    codon: str                         # "start codon" / "stop codon" when the event's region holds one
    features: pd.DataFrame


@dataclass
class ProteinChange:
    event_id: str
    status: str                        # 'pair', 'inclusion_only', 'exclusion_only', or why there is no suggestion
    event_type: str = ""
    match: Match | None = None
    inc: Isoform | None = None
    exc: Isoform | None = None
    effect: dict = field(default_factory=dict)          # kind, frame, text
    shown: list = field(default_factory=list)          # the isoforms to draw
    only_inc: pd.DataFrame = field(default_factory=pd.DataFrame)
    only_exc: pd.DataFrame = field(default_factory=pd.DataFrame)
    in_event: list = field(default_factory=list)       # feature labels at the event's residues

    @property
    def ok(self) -> bool:
        return self.status in SHOWN

    @property
    def n_tied(self) -> int:
        m = self.match
        return int(m.pairs.tied_best.sum()) if m is not None and len(m.pairs) else 0

    @property
    def context_similarity(self) -> float:
        m = self.match
        return float(m.pairs.context_similarity.iloc[0]) if self.status == "pair" else float("nan")

    @property
    def isoforms(self) -> list[Isoform]:
        return [i for i in (self.inc, self.exc) if i is not None]

    def summary(self) -> str:
        if not self.ok:
            return ""
        parts = [self.effect["text"]]
        if self.in_event:
            parts.append("at the event: " + _list(self.in_event, 4))
        oi, oe = labels(self.only_inc), labels(self.only_exc)
        if oi or oe:
            a, b = form_names(self.event_type)
            parts.append(f"features only with {a}: " + (_list(oi, 3) or "none")
                         + f"; only with {b}: " + (_list(oe, 3) or "none"))
        return "; ".join(parts)

    def row(self) -> dict:
        r = dict(event_id=self.event_id, protein_status=self.status)
        if self.ok:
            for k, iso in (("inc", self.inc), ("exc", self.exc)):
                if iso is None:
                    continue
                r.update({f"{k}_transcript": iso.transcript, f"{k}_name": iso.name, f"{k}_biotype": iso.biotype,
                          f"{k}_tsl": iso.tsl, f"{k}_context": iso.context, f"{k}_length": len(iso.protein),
                          f"{k}_event_residues": ";".join(f"{a}-{b}" for a, b in iso.event_aa),
                          f"{k}_insert_after": iso.insert_after})
            p = self.match.pairs
            r.update(n_candidates=len(p), n_tied=self.n_tied)
            if self.status == "pair":
                r.update(evidence=p.evidence.iloc[0], context_similarity=round(self.context_similarity, 6))
            r.update(effect=self.effect.get("kind", ""), frame=self.effect.get("frame", ""),
                     effect_short=self.effect.get("short", ""), effect_text=self.effect.get("text", ""),
                     features_at_event="; ".join(self.in_event),
                     features_only_inc="; ".join(labels(self.only_inc)),
                     features_only_exc="; ".join(labels(self.only_exc)), summary=self.summary())
        return r


def _cds_ends(tx: pd.DataFrame, strand: str) -> tuple[int, int] | None:
    cds = _cds(tx)
    if not cds:
        return None
    lo, hi = cds[0][0], cds[-1][1]
    return (lo, hi) if strand == "+" else (hi, lo)            # (first coding nt, last coding nt), transcript order


def _stop_codon(exons: list[Iv], last: int, strand: str) -> list[int]:
    """The genomic positions of the stop codon: the three nucleotides after the last coding one, in transcript order
    (the cache's coding coordinates leave the stop codon out, the GENCODE convention), across an exon junction if
    need be; [] when the transcript ends first (no annotated stop)."""
    out = []
    for a, b in (exons if strand == "+" else exons[::-1]):
        pos = range(max(a, last + 1), b + 1) if strand == "+" else range(min(b, last - 1), a - 1, -1)
        for x in pos:
            out.append(x)
            if len(out) == 3:
                return out
    return []


def _isoform(cache, tx, which, check, mask, strand) -> Isoform:
    t = str(tx.transcript_id.iloc[0])
    seq = cache.sequence(t)
    coding = tx.transcript_type.iloc[0] == "protein_coding" and len(seq) > 0 and bool(_cds(tx))
    exn = _merge(zip(tx.start.astype(int), tx.end.astype(int)))
    own = _intersect(mask, exn)
    aa, after, utr, codon = [], None, "", ""
    if coding:
        first, last = _cds_ends(tx, strand)
        lo_c, hi_c = min(first, last), max(first, last)
        if own:
            aa = _aa_spans(tx, own)
            ref = own
            inside = lambda x: any(a <= x <= b for a, b in own)  # noqa: E731
            stop = _stop_codon(exn, last, strand)
            codon = "start codon" if inside(first) else ("stop codon" if any(inside(x) for x in stop) else "")
        else:                                                   # the region is absent: where would it sit?
            up = ([b for _, b in exn if b < mask[0][0]] if strand == "+" else
                  [a for a, _ in exn if a > mask[-1][1]])
            ref = [(max(up), max(up))] if strand == "+" and up else ([(min(up), min(up))] if up else [])
            sp = _aa_spans(tx, ref) if ref else []
            after = sp[0][1] if sp else None
        if not aa and after is None and ref and not codon:     # no coding residue and no stop codon: a UTR
            a, b = ref[0][0], ref[-1][1]
            before = (b < lo_c) if strand == "+" else (a > hi_c)
            utr = "5′ UTR" if before else "3′ UTR"
    return Isoform(which, t, str(tx.transcript_name.iloc[0]), str(tx.transcript_type.iloc[0]),
                   tsl_rank(tx.transcript_support_level.iloc[0]), check.context, coding, seq if coding else "",
                   aa, after, utr, codon, cache.features_of(t) if coding else cache.features_of(""))


def _effect(pc: ProteinChange, geom) -> dict:
    """What the event does to the protein, in words (see the module docstring)."""
    fm = forms(geom)
    mxe = geom.event_type == "MXE"
    d_inc = _length(fm["INC"].inc[1:2]) if mxe else _length(_subtract(fm["INC"].inc, fm["EXC"].inc))
    d_exc = _length(fm["EXC"].inc[1:2]) if mxe else 0
    shift = (d_inc - d_exc) % 3
    frame = "in frame" if shift == 0 else "frameshift"
    inc, exc, sim = pc.inc, pc.exc, pc.context_similarity
    noun = NOUN.get(geom.event_type, "region")
    psi_form, other_form = form_names(geom.event_type)
    both = inc is not None and exc is not None and inc.coding and exc.coding
    if both and inc.protein == exc.protein:
        if inc.event_aa or exc.event_aa:        # the event is coding, yet both forms translate the same (e.g. an MXE)
            return dict(kind="identical", frame="", short="no protein change",
                        text=f"{inc.name} and {exc.name} encode the same protein, although the event lies in the "
                             "coding sequence")
        return dict(kind="utr", frame="", short="no protein change",
                    text=f"{inc.name} and {exc.name} encode the same protein: the event lies outside the coding "
                         f"sequence")
    if both:                                # the pair differs only at the event: compare the two proteins
        c = _compare(inc.protein, exc.protein, mxe, psi_form)
        a0 = inc.event_aa[0][0] if inc.event_aa else (exc.insert_after or 0)
        b0 = inc.event_aa[-1][1] if inc.event_aa else (exc.insert_after or 0)
        first, last = c["shared_prefix"] + 1, c["shared_prefix"] + c["inc_block"]      # differing residues (INC)
        confined = a0 - 2 <= first and last <= b0 + 2
        if sim == 1.0 or (sim >= 0.95 and shift == 0 and confined):
            pair = f"{inc.name} ({len(inc.protein)} aa) vs {exc.name} ({len(exc.protein)} aa)"
            k, ki, ke = c["kind"], c["inc_block"], c["exc_block"]
            short = {"insertion": f"+{ki - ke} aa, in frame", "deletion": f"−{ke - ki} aa, in frame",
                     "substitution": f"{ki} aa changed", "replacement": f"{ke} → {ki} aa",
                     "exchange": "exon swap, in frame",
                     "n_terminus": "new N-terminus", "c_terminus": "new C-terminus"}.get(k, k)
            if mxe and inc.event_aa and exc.event_aa and k == "exchange":
                (a, b), (a2, b2) = (inc.event_aa[0][0], inc.event_aa[-1][1]), (exc.event_aa[0][0], exc.event_aa[-1][1])
                return dict(kind="exchange", frame=frame, short=short, text=(
                    f"the PSI exon (residues {a}–{b} of {inc.name}) takes the place of the other exon (residues "
                    f"{a2}–{b2} of {exc.name}) in frame; the rest of the protein is the same: {pair}"))
            text = c["text"]
            if k in ("insertion", "substitution", "replacement") and inc.event_aa:   # locate it as drawn
                text = re.sub(r"\(residues \d+–\d+", f"(the {noun} spans residues {a0}–{b0} of {inc.name}", text,
                              count=1)
            return dict(kind=k, frame=frame, short=short, text=text + f"; {pair}")
    ref = inc if (inc is not None and inc.coding) else exc
    other = exc if ref is inc else inc
    name = ref.name
    if ref.utr:
        kind, frame, short = "utr", "", ref.utr
        text = f"the event lies in the {ref.utr} of {name}: the protein is unchanged"
    elif ref.codon:
        kind, short = ("start" if ref.codon == "start codon" else "end"), ref.codon
        where = (f"residues {ref.event_aa[0][0]}–{ref.event_aa[-1][1]}" if ref.event_aa
                 else f"after residue {len(ref.protein)}")      # the stop codon alone, after the last residue
        text = (f"the {noun} holds the {ref.codon} of {name} ({where}), so the other form changes the "
                f"protein's {'start' if kind == 'start' else 'end'}")
    elif mxe:
        a, b = ref.event_aa[0][0], ref.event_aa[-1][1]
        mine, theirs = (d_inc, d_exc) if ref.form == "INC" else (d_exc, d_inc)
        kind, short = "exchange", f"exon swap, {frame}"
        text = (f"the {'PSI' if ref.form == 'INC' else 'other'} exon ({mine} nt) encodes residues {a}–{b} of {name}; "
                f"the {'other' if ref.form == 'INC' else 'PSI'} exon is {theirs} nt, so the exchange "
                f"{'keeps' if shift == 0 else 'shifts'} the frame")
    elif ref.form == "INC":
        a, b = ref.event_aa[0][0], ref.event_aa[-1][1]
        kind = "insertion" if shift == 0 else "frameshift"
        short = f"+{d_inc // 3} aa, in frame" if shift == 0 else "frameshift"
        text = (f"{psi_form} adds {d_inc // 3} aa in frame (the {noun} spans residues {a}–{b} of {name})"
                if shift == 0
                else f"the {noun} is {d_inc} nt, not a multiple of 3: without it the frame shifts after residue "
                     f"{a - 1} of {name}")
    else:
        kind = "insertion" if shift == 0 else "frameshift"
        short = f"+{d_inc // 3} aa if no stop" if shift == 0 else "frameshift"
        text = (f"{psi_form} would add {d_inc} nt after residue {ref.insert_after} of {name}: "
                + (f"{d_inc // 3} aa in frame, unless they hold a stop codon" if shift == 0 else "a frameshift"))
    if both and sim >= WEAK_PAIR:
        text += (f"; the annotated pair, {inc.name} ({len(inc.protein)} aa) and {exc.name} ({len(exc.protein)} aa), "
                 f"also differs outside the event ({_pct(sim)} of coding positions shared)")
    elif both:
        text += (f"; the closest {other_form if ref is inc else psi_form} transcript ({other.name}) shares "
                 + (f"only {_pct(sim)} of its" if sim > 0 else "none of its") + " coding positions with it, so it is "
                 "not drawn")
    elif other is not None:
        text += (f"; the {other_form if ref is inc else psi_form} transcript ({other.name}) "
                 f"{BIOTYPE_TEXT.get(other.biotype, 'is ' + other.biotype.replace('_', ' ') + ' (noncoding)')}")
    return dict(kind=kind, frame=frame, short=short, text=text)


def protein_change(cache: ProteinCache, event_id: str, row, geom) -> ProteinChange:
    """The suggested transcripts of an event's two forms and how the protein changes."""
    full = forms(geom)
    if not full:
        return ProteinChange(event_id, "unsupported_event_type", geom.event_type)
    span = (min(a for f in full.values() for a, _ in f.inc), max(b for f in full.values() for _, b in f.inc))
    ex = cache.transcripts_for(row, span)
    if ex.empty:
        return ProteinChange(event_id, "gene_not_in_cache", geom.event_type)
    strand = str(row.get("strand", "") or "")
    strand = strand if strand in ("+", "-") else str(ex.strand.iloc[0])
    m = match_forms(ex, geom, strand)
    if m.status not in SHOWN:
        return ProteinChange(event_id, m.status, geom.event_type, match=m)
    pc = ProteinChange(event_id, m.status, geom.event_type, match=m)
    on_exons = {}
    for which, t in (("INC", m.inc), ("EXC", m.exc)):
        if not t:
            continue
        tx = ex[ex.transcript_id.eq(t)].sort_values(["start", "end", "exon_id"], kind="stable").reset_index(drop=True)
        iso = _isoform(cache, tx, which, m.checks[(which, t)], m.mask, strand)
        setattr(pc, which.lower(), iso)
        spans = [s for s in (_exon_aa(tx.iloc[i]) for i in sorted(set(m.checks[(which, t)].exons))) if s]
        on_exons[which] = _overlapping(iso.features, spans)
    coding = [i for i in pc.isoforms if i.coding]
    if not coding:
        pc.status = "noncoding"
        return pc
    pc.effect = _effect(pc, geom)
    close = len(coding) == 2 and pc.context_similarity >= WEAK_PAIR
    pc.shown = coding if close else coding[:1]
    if close:
        pc.only_inc = _only(on_exons["INC"], pc.exc.features)
        pc.only_exc = _only(on_exons["EXC"], pc.inc.features)
    at: list[str] = []
    for iso in pc.shown:
        spans = iso.event_aa or ([(iso.insert_after, iso.insert_after + 1)] if iso.insert_after else [])
        at += _named_at(iso.features, spans)
    pc.in_event = list(dict.fromkeys(at))
    return pc


def protein_changes(cache: ProteinCache, ds, events) -> tuple[pd.DataFrame, dict]:
    """One row per event (protein_status in SHOWN when there is a suggestion), and the ProteinChange objects."""
    from .events import geometry
    out, objs = [], {}
    ev = ds.events
    if "event_id" in ev.columns:
        ev = ev.drop_duplicates("event_id").set_index("event_id")
    for e in events:
        row = ev.loc[e].to_dict()
        try:
            g = geometry(e, row)
        except (InputError, KeyError) as err:
            objs[e] = ProteinChange(e, "no_geometry")
            out.append(dict(event_id=e, protein_status="no_geometry", note=str(err)))
            continue
        pc = protein_change(cache, e, row, g)
        objs[e] = pc
        out.append(pc.row())
    return pd.DataFrame(out), objs
