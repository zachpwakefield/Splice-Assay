"""Import rMATS output: events (with geometry) and per-sample PSI.

Events come from fromGTF.<TYPE>.txt (every detected event) or <TYPE>.MATS.<JC|JCEC>.txt (tested events); PSI comes
from the IncLevel1/IncLevel2 columns of <TYPE>.MATS.<JC|JCEC>.txt, whose values follow the BAM order of b1.txt and
b2.txt. rMATS coordinates are 0-based, half-open, so they are used as they are.

Event IDs are '<prefix><TYPE>:<rMATS ID>' (rMATS numbers events per type).

MXE: rMATS lists the two exons in genomic order. In the data this package was developed on, MXE PSI measured the
transcript-upstream exon: the first-listed exon on the plus strand and the second-listed on the minus strand
(checked on FGFR1-3 IIIb/IIIc). That is the default (`mxe_psi_exon="transcript_upstream"`). If your rMATS version
may differ, check a known event and use `mxe_psi_exon="first_listed"` if needed.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .dataset import InputError
from .events import EVENT_TYPES, format_intervals

COORDS = {
    "SE": ["exonStart_0base", "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
    "RI": ["riExonStart_0base", "riExonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE"],
    "MXE": ["1stExonStart_0base", "1stExonEnd", "2ndExonStart_0base", "2ndExonEnd", "upstreamES", "upstreamEE",
            "downstreamES", "downstreamEE"],
    "A3SS": ["longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"],
    "A5SS": ["longExonStart_0base", "longExonEnd", "shortES", "shortEE", "flankingES", "flankingEE"],
}


def _unquote(s: pd.Series) -> pd.Series:
    return s.astype(str).str.strip().str.strip('"')


def _read(path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype={"GeneID": str, "geneSymbol": str, "chr": str, "strand": str,
                                              "IncLevel1": str, "IncLevel2": str}, low_memory=False)


def _find(folder: Path, event_type: str, source: str) -> Path:
    names = {"fromGTF": [f"fromGTF.{event_type}.txt"], "JC": [f"{event_type}.MATS.JC.txt"],
             "JCEC": [f"{event_type}.MATS.JCEC.txt"]}
    order = ["fromGTF", "JCEC", "JC"] if source == "auto" else [source]
    for s in order:
        for n in names[s]:
            if (folder / n).exists():
                return folder / n
    raise InputError(f"{folder}: no rMATS {event_type} table ({', '.join(n for s in order for n in names[s])})")


def geometry_columns(df: pd.DataFrame, event_type: str, mxe_psi_exon: str = "transcript_upstream") -> pd.DataFrame:
    """constant and variable columns (0-based, half-open) from the rMATS coordinate columns of one event type."""
    t = event_type.upper()
    missing = [c for c in COORDS[t] if c not in df.columns]
    if missing:
        raise InputError(f"rMATS {t}: missing column(s) {', '.join(missing)}")
    c = {k: df[k].to_numpy(np.int64) for k in COORDS[t]}
    const, var = [], []
    for i in range(len(df)):
        if t == "SE":
            k = sorted([(c["upstreamES"][i], c["upstreamEE"][i]), (c["downstreamES"][i], c["downstreamEE"][i])])
            v = [(c["exonStart_0base"][i], c["exonEnd"][i])]
        elif t == "RI":
            k = sorted([(c["upstreamES"][i], c["upstreamEE"][i]), (c["downstreamES"][i], c["downstreamEE"][i])])
            v = [(k[0][1], k[1][0])]
        elif t == "MXE":
            k = sorted([(c["upstreamES"][i], c["upstreamEE"][i]), (c["downstreamES"][i], c["downstreamEE"][i])])
            first = (c["1stExonStart_0base"][i], c["1stExonEnd"][i])
            second = (c["2ndExonStart_0base"][i], c["2ndExonEnd"][i])
            if mxe_psi_exon == "first_listed" or df["strand"].iloc[i] == "+":
                v = [first, second]
            elif mxe_psi_exon == "transcript_upstream":
                v = [second, first]
            else:
                raise ValueError("mxe_psi_exon must be 'transcript_upstream' or 'first_listed'")
        else:                                                   # A3SS / A5SS
            short = (c["shortES"][i], c["shortEE"][i])
            flank = (c["flankingES"][i], c["flankingEE"][i])
            lo, hi = c["longExonStart_0base"][i], c["longExonEnd"][i]
            v = [(lo, short[0])] if lo < short[0] else [(short[1], hi)]
            k = sorted([short, flank])
        const.append(format_intervals((int(a), int(b)) for a, b in k))
        var.append(format_intervals((int(a), int(b)) for a, b in v))
    return pd.DataFrame({"constant": const, "variable": var}, index=df.index)


def read_events(folder, event_types=EVENT_TYPES, source: str = "auto", prefix: str = "",
                mxe_psi_exon: str = "transcript_upstream") -> pd.DataFrame:
    """The events table (event_id, gene, gene_id, chrom, strand, event_type, constant, variable, label)."""
    folder = Path(folder)
    out = []
    for t in event_types:
        t = t.upper()
        if t not in COORDS:
            raise InputError(f"rMATS event type {t} is not one of {', '.join(COORDS)}")
        df = _read(_find(folder, t, source))
        if df.empty:
            continue
        geo = geometry_columns(df, t, mxe_psi_exon)
        gene = _unquote(df["geneSymbol"])
        gid = _unquote(df["GeneID"])
        gene = gene.where(~gene.isin(["NA", "nan", ""]), gid)
        out.append(pd.DataFrame({"event_id": prefix + t + ":" + df["ID"].astype(str), "gene": gene.to_numpy(),
                                 "gene_id": gid.to_numpy(), "chrom": df["chr"].astype(str).to_numpy(),
                                 "strand": df["strand"].astype(str).to_numpy(), "event_type": t,
                                 "constant": geo.constant.to_numpy(), "variable": geo.variable.to_numpy(),
                                 "label": t + ":" + df["ID"].astype(str)}))
    if not out:
        raise InputError(f"{folder}: no rMATS events found")
    return pd.concat(out, ignore_index=True)


def sample_names(b_file) -> list[str]:
    """Sample names from an rMATS b1.txt/b2.txt: BAM file names without directory and '.bam'."""
    text = Path(b_file).read_text().strip()
    return [Path(p.strip()).name.removesuffix(".bam") for p in text.split(",") if p.strip()]


def read_psi(mats_file, names1, names2=None, prefix: str = "") -> pd.DataFrame:
    """Long PSI (event_id, sample_id, psi) from the IncLevel1/IncLevel2 columns of a <TYPE>.MATS.*.txt file."""
    path = Path(mats_file)
    t = path.name.split(".")[0].upper()
    df = _read(path)
    rows = []
    for col, names in (("IncLevel1", names1), ("IncLevel2", names2)):
        if not names:
            continue
        vals = df[col].fillna("").astype(str).str.split(",")
        lens = vals.map(len)
        if (lens != len(names)).any():
            bad = int(lens[lens != len(names)].iloc[0])
            raise InputError(f"{path.name}: {col} has {bad} values per event but {len(names)} sample names were given")
        mat = pd.DataFrame(vals.tolist(), columns=list(names))
        mat = mat.apply(lambda s: pd.to_numeric(s.where(~s.str.strip().isin(["NA", ""])), errors="coerce"))
        mat.insert(0, "event_id", prefix + t + ":" + df["ID"].astype(str))
        rows.append(mat.melt(id_vars="event_id", var_name="sample_id", value_name="psi"))
    if not rows:
        raise InputError("read_psi: no sample names given")
    return pd.concat(rows, ignore_index=True)


def import_rmats(folder, b1, b2=None, counting: str = "JCEC", event_types=EVENT_TYPES, prefix: str = "",
                 names1=None, names2=None, mxe_psi_exon: str = "transcript_upstream"):
    """(events, psi) for the event types present in an rMATS output folder, from its <TYPE>.MATS.<counting>.txt."""
    folder = Path(folder)
    n1 = list(names1) if names1 is not None else sample_names(b1)
    n2 = list(names2) if names2 is not None else (sample_names(b2) if b2 else None)
    dup = set(n1) & set(n2 or [])
    if dup or len(set(n1)) < len(n1) or (n2 and len(set(n2)) < len(n2)):
        raise InputError("rMATS sample names must be unique across b1 and b2")
    ev, ps, seen = [], [], []
    for t in event_types:
        f = folder / f"{t.upper()}.MATS.{counting}.txt"
        if not f.exists():
            continue
        seen.append(t.upper())
        try:
            no_rows = pd.read_csv(f, sep="\t", nrows=1).empty
        except pd.errors.EmptyDataError:                 # not even a header line: a broken output
            raise InputError(f"{f}: the file is empty (rMATS writes a header line even when it finds no event)") \
                from None
        if no_rows:                                      # rMATS writes every type's file, only a header when it
            continue                                     # found no event of that type
        ev.append(read_events(folder, [t], source=counting, prefix=prefix, mxe_psi_exon=mxe_psi_exon))
        ps.append(read_psi(f, n1, n2, prefix))
    if not seen:
        raise InputError(f"{folder}: no *.MATS.{counting}.txt files")
    if not ev:
        raise InputError(f"{folder}: no rMATS events found in the {counting} files of {', '.join(seen)} (each holds "
                         "only a header)")
    return pd.concat(ev, ignore_index=True), pd.concat(ps, ignore_index=True)
