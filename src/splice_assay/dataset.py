"""The input contract: the tables splice-assay reads, their columns, and the checks applied to them.

Every table is a CSV, TSV (.tsv or .txt) or Parquet file, or a pandas DataFrame. Column names below are the logical
names; tables may use their own names through `columns={"sample_id": "File.ID", "cohort": "cancer", ...}`.

samples     sample_id, patient_id, cohort, group (e.g. tumour / normal; `tissue` is accepted for `group`); optional
            survival_cohort (true/false: which case samples enter survival analysis; default every case sample,
            which must then be one per patient and cohort)
psi         long: event_id, sample_id, psi; or wide: event_id plus one column per sample. PSI on a 0-1 scale;
            empty = not measured
events      event_id, gene; for figures also chrom, strand, event_type, constant, variable (see events.py); optional
            label, gene_id, expression_gene, psi_junctions, other_junctions
survival    long: patient_id, endpoint, time, event (1 = event, 0 = censored); or wide, one pair of columns per
            endpoint: <EP>.time (or <EP>_time) with <EP>, <EP>.event or <EP>.status (e.g. OS.time and OS). Rows
            with a missing or non-positive time or an event other than 0/1 are dropped (and counted). sample_id may
            stand in for patient_id
pairs       cohort, patient_id, case_sample, reference_sample (tumour_sample / normal_sample are accepted). Optional:
            without it, a patient with exactly one case and one reference sample in a cohort forms a pair
expression  long: gene, sample_id, value; or wide: gene plus one column per sample. Optional host-gene expression
            (e.g. log2(TPM + 1)); when given, the Cox model adjusts PSI for it
clinical    patient_id plus any columns (age, sex, stage, ...). Optional covariates for the Cox model

Groups: any two conditions can be compared, and their names (as spelled in `group`) are used in every table and
figure. Which group is the case and which the reference (effects are case minus reference) comes from, in order:
  - `case` and `reference` (--case / --reference): the group values compared (matching ignores case, "tumor" =
    "tumour");
  - an optional `role` column in samples: case or reference (control is accepted for reference) for each sample;
  - otherwise "tumour" and "normal".
Samples of any other group are kept aside and not tested.

Fewer tables: the samples table may carry the patient-level columns, and the psi table the event columns.
  samples   + survival columns as in a wide survival table (OS.time with OS, ...), used when no survival table is
            given; + any other columns as clinical covariates, used when no clinical table is given; + an optional
            pair_id (samples sharing a pair_id within a cohort form a pair), used when no pairs table is given.
            Survival and clinical values are read from each patient's survival sample (see survival_cohort).
  psi       + event columns (gene, chrom, strand, event_type, constant, variable, ...) beside the PSI columns (wide:
            one column per sample; long: repeated on each row), used when no events table is given.
So two tables (samples, psi) are enough; separate tables are still read, and take precedence.

Subsets: `keep` (--keep FILE) restricts the run to the patients listed (e.g. one subtype). A listed ID selects a
patient when it is that patient's patient_id or one of its sample_ids; any other ID selects the patient of its longest
leading part, up to a "-", that is one (so TCGA aliquot barcodes select their patient). IDs are read as text; a list
without a header line works too (its first line is read as an ID when it names a patient or sample). `where`
(--where COLUMN=VALUE[,VALUE...], repeatable) selects the patients whose clinical row, or any of whose samples, has
one of the values in that column (case-insensitive text, or numbers by value: 1 matches 1.0; a missing value matches
nothing); every condition must hold, also two on one column, and with keep both apply. All samples of a selected
patient are kept, so its normals stay. Tables keyed by sample_id (survival, clinical) may list samples that the subset
leaves out.
"""
from __future__ import annotations

import re
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

TABLES = ("samples", "psi", "events", "survival", "pairs", "expression", "clinical")
SUFFIXES = (".csv", ".tsv", ".txt", ".csv.gz", ".tsv.gz", ".txt.gz", ".parquet")
TEXT_COLUMNS = {"sample_id", "patient_id", "cohort", "group", "tissue", "survival_cohort", "event_id", "gene",
                "gene_id", "expression_gene", "label", "chrom", "strand", "event_type", "constant", "variable",
                "psi_junctions", "other_junctions", "endpoint", "case_sample", "reference_sample", "tumour_sample",
                "tumor_sample", "normal_sample"}
ALIASES = {"group": ("tissue",), "case_sample": ("tumour_sample", "tumor_sample"),
           "reference_sample": ("normal_sample",)}
TRUE, FALSE = {"true", "t", "1", "yes", "y"}, {"false", "f", "0", "no", "n"}
EVENT_OPTIONAL = ("label", "gene_id", "expression_gene", "chrom", "strand", "event_type", "constant", "variable",
                  "psi_junctions", "other_junctions")
TIME_COLUMN = re.compile(r"^(?P<ep>.+?)[._ ]?time$", re.IGNORECASE)
SAMPLE_OWN = {"sample_id", "patient_id", "cohort", "group", "tissue", "survival_cohort", "pair_id", "role"}
EVENT_COLUMNS = ("event_id", "gene") + EVENT_OPTIONAL


class InputError(ValueError):
    """An input table does not meet the contract; the message says what to change."""


# ============================================================================================ reading
def read_table(path, columns: dict | None = None, na_values=None, usecols=None) -> pd.DataFrame:
    """Read a CSV, TSV (.tsv, .txt) or Parquet table. ID columns stay text; numbers parse exactly (round trip).
    `na_values` are extra codes read as missing (e.g. "missing", "[Not Available]"), besides pandas' own; `usecols`
    reads only the named columns that exist."""
    p = Path(path)
    name = p.name.lower()
    header = table_columns(p)
    use = None if usecols is None else [c for c in header if c in set(usecols)]
    if name.endswith(".parquet"):
        return with_missing(pd.read_parquet(p, columns=use), na_values)
    sep = "\t" if name.endswith((".tsv", ".txt", ".tsv.gz", ".txt.gz")) else ","
    text = TEXT_COLUMNS | set((columns or {}).values())
    dtype = {c: str for c in header if c in text}
    return pd.read_csv(p, sep=sep, dtype=dtype, float_precision="round_trip", low_memory=False,
                       na_values=list(na_values or []), usecols=use)


def table_columns(path) -> list[str]:
    """The column names of a table file, without reading its rows."""
    p = Path(path)
    if p.name.lower().endswith(".parquet"):
        import pyarrow.parquet as pq
        return list(pq.read_schema(p).names)
    sep = "\t" if p.name.lower().endswith((".tsv", ".txt", ".tsv.gz", ".txt.gz")) else ","
    return list(pd.read_csv(p, sep=sep, nrows=0).columns)


def events_table(source, columns: dict | None = None, na_values=None) -> pd.DataFrame:
    """The events table (normalised) of a data folder or file: its events table, or else the event columns of its psi
    table (read without the PSI values)."""
    p = Path(source)
    if p.is_dir():
        found = find_tables(p)
        p = found.get("events") or found.get("psi")
        if p is None:
            raise InputError(f"{source}: no events or psi table")
    want = set(EVENT_COLUMNS) | set((columns or {}).values())
    return normalise_events(read_table(p, columns, na_values, usecols=want), columns=columns)


def with_missing(df: pd.DataFrame, na_values) -> pd.DataFrame:
    """The table with the given codes replaced by missing values (text cells only)."""
    if not na_values:
        return df
    codes = set(map(str, na_values))
    obj = [c for c in df.columns if df[c].dtype == object]
    out = df.copy()
    for c in obj:
        out[c] = out[c].where(~out[c].astype(str).str.strip().isin(codes))
    return out


def check_table_names(names) -> None:
    """Every name must be one of TABLES: a misspelt one would leave the folder's table of that name in use."""
    unknown = sorted(set(names) - set(TABLES))
    if unknown:
        raise InputError(f"unknown table name{'s' if len(unknown) > 1 else ''} "
                         f"{', '.join(map(repr, unknown))}: expected one of {', '.join(TABLES)}")


def find_tables(folder) -> dict[str, Path]:
    """Files named <table>.<csv|tsv|txt|parquet>[.gz] in a folder, e.g. samples.csv, psi.parquet."""
    folder = Path(folder)
    if not folder.is_dir():
        raise InputError(f"not a folder: {folder}")
    out = {}
    for t in TABLES:
        hits = [folder / f"{t}{s}" for s in SUFFIXES if (folder / f"{t}{s}").exists()]
        if len(hits) > 1:
            raise InputError(f"{folder}: more than one {t} table ({', '.join(h.name for h in hits)})")
        if hits:
            out[t] = hits[0]
    return out


# ============================================================================================ small helpers
def _map_columns(df: pd.DataFrame, columns: dict | None) -> pd.DataFrame:
    """Rename the user's columns to logical names (mapping, then aliases such as tissue -> group)."""
    ren = {}
    for logical, actual in (columns or {}).items():
        if actual in df.columns and logical not in df.columns and actual != logical:
            ren[actual] = logical
    df = df.rename(columns=ren)
    for logical, alts in ALIASES.items():
        if logical not in df.columns:
            for a in alts:
                if a in df.columns:
                    df = df.rename(columns={a: logical})
                    break
    return df


def _require(df: pd.DataFrame, table: str, cols) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise InputError(f"{table}: missing column(s) {', '.join(missing)} (found: {', '.join(map(str, df.columns))}). "
                         "Name them so, or map them, e.g. columns={'cohort': 'cancer'} or --column cohort=cancer")


def _examples(values, k: int = 3) -> str:
    v = list(dict.fromkeys(map(str, values)))
    return ", ".join(v[:k]) + (f" and {len(v) - k} more" if len(v) > k else "")


def _ids(s: pd.Series, table: str, col: str) -> pd.Series:
    if s.isna().any():
        raise InputError(f"{table}: column {col} has {int(s.isna().sum())} empty value(s)")
    out = s.astype(str).str.strip()
    if (out == "").any():
        raise InputError(f"{table}: column {col} has {int((out == '').sum())} empty value(s)")
    return out


def _num(s: pd.Series, table: str, col: str) -> pd.Series:
    out = pd.to_numeric(s, errors="coerce")
    bad = s[out.isna() & s.notna() & (s.astype(str).str.strip() != "")]
    if len(bad):
        raise InputError(f"{table}: column {col} must be numeric; found {_examples(bad)}")
    return out.astype(float)


def _bool(s: pd.Series, table: str, col: str) -> np.ndarray:
    if s.dtype == bool:
        return s.to_numpy()
    t = s.astype(str).str.strip().str.lower()
    ok = t.isin(TRUE | FALSE)
    if not ok.all():
        raise InputError(f"{table}: column {col} must be true or false; found {_examples(s[~ok])}")
    return t.isin(TRUE).to_numpy()


def _text(df: pd.DataFrame, col: str) -> pd.Series:
    if col not in df.columns:
        return pd.Series([""] * len(df), index=df.index, dtype=object)
    return df[col].map(lambda v: "" if pd.isna(v) else str(v).strip())


def _canon(v) -> str:
    return str(v).strip().lower().replace("tumor", "tumour")


def _as_tuple(v) -> tuple:
    return (v,) if isinstance(v, str) else tuple(v)


def _label(values: pd.Series, fallback: str) -> str:
    """Display label of a group: its most common spelling in the data, capitalised when all lower case."""
    v = values.value_counts().index[0] if len(values) else fallback
    return v[:1].upper() + v[1:] if v == v.lower() else v


# ============================================================================================ tables
ROLE_NAMES = {"case": "case", "reference": "reference", "control": "reference"}


def normalise_samples(df: pd.DataFrame, case=None, reference=None, columns=None) -> tuple[pd.DataFrame, dict]:
    """Samples with a `role` (case / reference / other) and `survival_cohort`; also the display labels (the group
    names as spelled in the data). The roles come from `case` / `reference`, else a `role` column, else tumour /
    normal (see the module docstring)."""
    df = _map_columns(df, columns)
    _require(df, "samples", ["sample_id", "patient_id", "cohort", "group"])
    s = pd.DataFrame({c: _ids(df[c], "samples", c).to_numpy() for c in ("sample_id", "patient_id", "cohort")})
    dup = s.sample_id[s.sample_id.duplicated()]
    if len(dup):
        raise InputError(f"samples: sample_id must be unique; repeated: {_examples(dup)}")
    s["group"] = _ids(df["group"], "samples", "group").to_numpy()
    if case is None and reference is None and "role" in df.columns:          # the samples table says which is which
        r = [ROLE_NAMES.get(str(v).strip().lower(), "other") if pd.notna(v) else "other" for v in df["role"]]
        s["role"] = r
        both = s[s.role.ne("other")].groupby("group").role.nunique()
        if (both > 1).any():
            raise InputError(f"samples: group {_examples(both.index[both > 1])} has both roles; each group needs one "
                             "role (case or reference)")
        case = sorted(set(s.group[s.role.eq("case")])) or ["case"]
        reference = sorted(set(s.group[s.role.eq("reference")])) or ["reference"]
        if not (s.role == "case").any():
            raise InputError(f"samples: the role column names no case sample (values found: "
                             f"{_examples(df['role'].dropna().unique())}; use case and reference)")
    else:
        case = "tumour" if case is None else case
        reference = "normal" if reference is None else reference
        cases, refs = {_canon(v) for v in _as_tuple(case)}, {_canon(v) for v in _as_tuple(reference)}
        if cases & refs:
            raise InputError("case and reference groups must differ")
        g = s.group.map(_canon)
        s["role"] = np.where(g.isin(cases), "case", np.where(g.isin(refs), "reference", "other"))
        if not (s.role == "case").any():
            raise InputError(f"samples: no sample of the case group {_examples(_as_tuple(case))}; groups found: "
                             f"{_examples(s.group.unique())}. Name the groups with --case and --reference (case= and "
                             "reference= in Python), or add a role column (case / reference) to the samples table")
    labels = dict(case=_label(s.group[s.role == "case"], _as_tuple(case)[0]),
                  reference=_label(s.group[s.role == "reference"], _as_tuple(reference)[0]))
    others = {str(k): int(v) for k, v in s.group[s.role == "other"].value_counts().items()}
    if others:
        warnings.warn(f"samples: {sum(others.values())} sample(s) of other groups are not tested "
                      f"({_examples(others)})", stacklevel=3)
    if "survival_cohort" in df.columns:
        sc = _bool(df["survival_cohort"], "samples", "survival_cohort")
        bad = sc & (s.role != "case").to_numpy()
        if bad.any():
            raise InputError(f"samples: survival_cohort is true for non-case sample(s) {_examples(s.sample_id[bad])}")
        s["survival_cohort"] = sc
    else:
        s["survival_cohort"] = (s.role == "case").to_numpy()
    base = s[s.survival_cohort]
    dup = base[base.duplicated(["cohort", "patient_id"], keep=False)]
    if len(dup):
        raise InputError(f"samples: {dup.patient_id.nunique()} patient(s) have more than one {labels['case']} sample "
                         f"in the survival cohort of one cohort (e.g. {_examples(dup.patient_id)}). Add a column "
                         "survival_cohort that is true for exactly one of them per patient.")
    return s, dict(labels=labels, other_groups=others)


def normalise_psi(df: pd.DataFrame, samples: pd.DataFrame, event_ids=None, columns=None,
                  subset: bool = False, signed=()) -> pd.DataFrame:
    """Wide float matrix: one row per event_id, one column per sample of the samples table (its order). Values lie in
    [0, 1], except for the events in `signed` (HIT index events, -1 to 1)."""
    df = _map_columns(df, columns)
    if "event_id" not in df.columns and df.index.name == "event_id":
        df = df.reset_index()
    keep = None if event_ids is None else set(map(str, event_ids))
    if {"event_id", "sample_id", "psi"} <= set(df.columns):
        d = pd.DataFrame({"event_id": _ids(df["event_id"], "psi", "event_id").to_numpy(),
                          "sample_id": _ids(df["sample_id"], "psi", "sample_id").to_numpy(),
                          "psi": _num(df["psi"], "psi", "psi").to_numpy()})
        if keep is not None:
            d = d[d.event_id.isin(keep)]
        dup = d[d.duplicated(["event_id", "sample_id"])]
        if len(dup):
            raise InputError(f"psi: {len(dup)} repeated (event_id, sample_id) row(s), e.g. "
                             f"{dup.event_id.iloc[0]} / {dup.sample_id.iloc[0]}")
        wide = d.pivot(index="event_id", columns="sample_id", values="psi")
    elif "event_id" in df.columns:
        w = df.copy()
        w["event_id"] = _ids(w["event_id"], "psi", "event_id")
        if keep is not None:
            w = w[w.event_id.isin(keep)]
        dup = w.event_id[w.event_id.duplicated()]
        if len(dup):
            raise InputError(f"psi: event_id must be unique in a wide table; repeated: {_examples(dup)}")
        wide = w.set_index("event_id")
        wide.columns = [str(c).strip() for c in wide.columns]
        wide = wide.apply(lambda c: _num(c, "psi", c.name))
    else:
        raise InputError("psi: expected columns event_id, sample_id, psi (long) or event_id plus one column per "
                         "sample (wide)")
    known = set(samples.sample_id)
    extra = [c for c in wide.columns if c not in known]
    if len(wide.columns) and len(extra) == len(wide.columns):
        raise InputError(f"psi: no sample matches samples.sample_id (psi has e.g. {_examples(extra)})")
    if extra and not subset:
        warnings.warn(f"psi: {len(extra)} sample(s) not in the samples table are ignored (e.g. {_examples(extra)})",
                      stacklevel=3)
    wide = wide.reindex(columns=samples.sample_id.to_numpy()).astype(float)
    v = wide.to_numpy()
    if np.isinf(v).any():
        raise InputError("psi: infinite values; use an empty cell for a missing PSI")
    hit = wide.index.isin(set(signed))
    for rows, lo, what in ((~hit, 0.0, "[0, 1]"), (hit, -1.0, "[-1, 1] (HIT index)")):
        fin = v[rows][np.isfinite(v[rows])]
        if fin.size and (fin.min() < lo - 1e-9 or fin.max() > 1 + 1e-9):
            hint = " Values above 1 look like percentages: divide by 100." if fin.max() > 1 else ""
            raise InputError(f"psi: values must lie in {what}; found {fin.min():g} to {fin.max():g}.{hint}")
        if fin.size and (fin.min() < lo or fin.max() > 1):
            wide.loc[rows] = wide.loc[rows].clip(lo, 1.0)
    wide.index.name, wide.columns.name = "event_id", "sample_id"
    return wide


def _keep_list(keep, column: str | None = None) -> tuple[list[str], str | None]:
    """The IDs of a keep list and, for a table read by its first column, that column's name (None otherwise). Every
    cell is read as text, so an ID such as 00123 keeps its leading zeros."""
    head = None
    if isinstance(keep, (str, Path)):
        p = Path(keep)
        try:
            if p.name.lower().endswith(".parquet"):
                t = pd.read_parquet(p)
            else:
                sep = "\t" if p.name.lower().endswith((".tsv", ".txt", ".tsv.gz", ".txt.gz")) else ","
                t = pd.read_csv(p, sep=sep, dtype=str, keep_default_na=False)
        except pd.errors.EmptyDataError:                     # no line at all (or only blank ones)
            return [], None
        except OSError as e:
            raise InputError(f"keep: cannot read {keep} ({e.strerror or e})") from None
        col = column or t.columns[0]
        if col not in t.columns:
            raise InputError(f"{keep}: no column {col} (found: {', '.join(map(str, t.columns))})")
        head = None if column else str(col).strip()
        keep = t[col]
    return [str(v).strip() for v in keep if pd.notna(v) and str(v).strip()], head


def keep_ids(keep, column: str | None = None) -> list[str]:
    """The IDs of a keep list: a list of IDs, or a table (CSV/TSV/Parquet) and its `column` (default the first). IDs
    are read as text."""
    ids, _ = _keep_list(keep, column)
    if not ids and isinstance(keep, (str, Path)):
        raise InputError(f"{keep}: the keep list is empty")
    return ids


def where_conditions(where) -> list[tuple[str, list[str]]]:
    """[(column, [values])] from 'COLUMN=VALUE[,VALUE...]' strings or a dict. Every condition must hold, each by one
    of its values, so a column named in two conditions must meet both (stage=I,II and stage=II,III select II)."""
    if isinstance(where, dict):
        return [(str(k), [str(v)] if isinstance(v, str) else [str(x) for x in v]) for k, v in where.items()]
    out: list[tuple[str, list[str]]] = []
    for w in [where] if isinstance(where, str) else where:
        col, sep, vals = str(w).partition("=")
        values = [v.strip() for v in vals.split(",") if v.strip()]
        if not sep or not col.strip() or not values:
            raise InputError(f"where: expected COLUMN=VALUE[,VALUE...], not {w!r}")
        out.append((col.strip(), values))
    return out


def _where_hits(values: pd.Series, wanted) -> pd.Series:
    """Cells equal to one of the wanted values: as text, ignoring case and surrounding spaces, or as numbers when both
    are numbers (so 1 matches a column read as 1.0). A missing cell matches nothing."""
    want = {str(v).strip().lower() for v in wanted}
    hit = values.astype(str).str.strip().str.lower().isin(want)
    exact = 2.0 ** 53                                       # beyond it floats are not exact: long codes match as text
    nums = {float(x) for x in pd.to_numeric(pd.Series(sorted(want)), errors="coerce").dropna() if abs(x) < exact}
    if nums:
        num = pd.to_numeric(values, errors="coerce")
        hit |= num.isin(nums) & (num.abs() < exact)
    return hit & values.notna()


def _where_mask(raw: pd.DataFrame, clinical: pd.DataFrame | None, conds: list[tuple[str, list[str]]]) -> np.ndarray:
    """Rows (samples) of patients meeting every condition, read from the samples table (any of the patient's samples)
    or else the clinical table (one row per patient); values compare as `_where_hits` says."""
    pid = raw["patient_id"].astype(str).str.strip()
    chosen = set(pid)
    for col, vals in conds:
        if col in raw.columns:
            hit = pid[_where_hits(raw[col], vals).to_numpy()]
        elif clinical is not None and col in clinical.columns and "patient_id" in clinical.columns:
            cp = clinical["patient_id"].astype(str).str.strip()
            hit = cp[_where_hits(clinical[col], vals).to_numpy()]
        else:
            cols = list(raw.columns) + ([] if clinical is None else list(clinical.columns))
            raise InputError(f"where: no column {col} in the samples or clinical table (found: "
                             f"{', '.join(map(str, dict.fromkeys(cols)))})")
        chosen &= set(hit)
    return pid.isin(chosen).to_numpy()


def _keep_mask(raw: pd.DataFrame, ids) -> np.ndarray:
    """Rows (samples) of patients selected by the keep list. A listed ID that is a patient or sample ID selects that
    patient; any other selects the patient of its longest leading '-'-separated part that is one (a TCGA aliquot
    barcode TCGA-XX-0001-01A-11R-... selects TCGA-XX-0001), so an ID such as 1-2 never also selects patient 1."""
    sid = raw["sample_id"].astype(str).str.strip()
    pid = raw["patient_id"].astype(str).str.strip()
    owner = dict(zip(sid, pid))
    patients = set(pid)
    chosen = set()
    for i in ids:
        parts = str(i).split("-")
        for k in range(len(parts), 0, -1):                 # the ID itself first, then shorter leading parts
            key = "-".join(parts[:k])
            if key in patients or key in owner:
                chosen |= {key} & patients | ({owner[key]} if key in owner else set())
                break
    return pid.isin(chosen).to_numpy()


def split_psi_events(df: pd.DataFrame, samples: pd.DataFrame, columns=None, table: str = "psi"):
    """(psi, events) from one table holding both: the PSI columns (wide: the samples' IDs; long: sample_id and psi)
    and the event columns (all others; one value per event)."""
    df = _map_columns(df, columns)
    if "event_id" not in df.columns and df.index.name == "event_id":
        df = df.reset_index()
    _require(df, table, ["event_id"])
    if "gene" not in df.columns:
        raise InputError(f"no events table was given, so the {table} table needs the event columns too: at least gene "
                         "(and chrom, strand, event_type, constant and variable to draw); or add an events table")
    if {"sample_id", "psi"} <= set(df.columns):                       # long: event columns repeated on each row
        ann = df.drop(columns=["sample_id", "psi"]).drop_duplicates()
        dup = ann.event_id[ann.event_id.duplicated()]
        if len(dup):
            raise InputError(f"{table}: the event columns differ between rows of one event, e.g. {_examples(dup)}")
        return df[["event_id", "sample_id", "psi"]], ann
    known = set(samples.sample_id)
    sample_cols = [c for c in df.columns if str(c).strip() in known]
    if not sample_cols:
        raise InputError(f"{table}: no column is a sample of the samples table (expected one PSI column per "
                         f"sample, e.g. {_examples(samples.sample_id)})")
    ann = [c for c in df.columns if c not in sample_cols]
    return df[["event_id"] + sample_cols], df[ann]


def _from_samples(raw: pd.DataFrame, s: pd.DataFrame, endpoints=None) -> dict:
    """The survival, clinical and pairs tables carried by the samples table (each None when it carries none)."""
    out = dict(survival=None, clinical=None, pairs=None, endpoint_columns=[])
    base = s.survival_cohort.to_numpy()
    pid = s.patient_id.to_numpy()[base]
    eps = endpoints or wide_endpoints([c for c in raw.columns if c not in SAMPLE_OWN])
    if eps:
        rows = raw[base]
        sv = pd.concat([pd.DataFrame({"patient_id": pid, "endpoint": ep, "time": rows[t].to_numpy(),
                                      "event": rows[e].to_numpy()}) for ep, (t, e) in eps.items()], ignore_index=True)
        out["survival"] = sv[~(sv.time.isna() & sv.event.isna())]           # both cells empty: no record
        out["endpoint_columns"] = [c for te in eps.values() for c in te]
    extra = [c for c in raw.columns if c not in SAMPLE_OWN and c not in out["endpoint_columns"]]
    if extra:
        out["clinical"] = raw.loc[base, extra].assign(patient_id=pid).drop_duplicates()
    if "pair_id" in raw.columns:
        d = s.assign(pair_id=[("" if pd.isna(v) else str(v).strip()) for v in raw["pair_id"]])
        d = d[d.pair_id.ne("") & d.role.isin(["case", "reference"])]
        rows, bad = [], []
        for (cohort, pair), g in d.groupby(["cohort", "pair_id"], sort=False):
            c, r = g[g.role.eq("case")], g[g.role.eq("reference")]
            if len(c) == 1 and len(r) == 1:
                rows.append(dict(cohort=cohort, patient_id=c.patient_id.iloc[0], case_sample=c.sample_id.iloc[0],
                                 reference_sample=r.sample_id.iloc[0]))
            else:
                bad.append(f"{pair} in {cohort} ({len(c)} case, {len(r)} reference)")
        if bad:
            raise InputError(f"samples: a pair_id must name one case and one reference sample of a cohort; "
                             f"{_examples(bad)}")
        out["pairs"] = pd.DataFrame(rows, columns=["cohort", "patient_id", "case_sample", "reference_sample"])
    return out


def normalise_events(df: pd.DataFrame, event_ids=None, columns=None) -> pd.DataFrame:
    """One row per event, indexed by event_id; optional text columns filled with ''."""
    df = _map_columns(df, columns)
    _require(df, "events", ["event_id", "gene"])
    e = pd.DataFrame({"event_id": _ids(df["event_id"], "events", "event_id").to_numpy(),
                      "gene": _ids(df["gene"], "events", "gene").to_numpy()})
    for c in EVENT_OPTIONAL:
        e[c] = _text(df, c).to_numpy()
    dup = e.event_id[e.event_id.duplicated()]
    if len(dup):
        raise InputError(f"events: event_id must be unique; repeated: {_examples(dup)}")
    e["label"] = np.where(e.label == "", e.event_id, e.label)
    e["expression_gene"] = np.where(e.expression_gene == "", e.gene, e.expression_gene)
    e["event_type"] = e.event_type.str.upper()
    bad = e[~e.strand.isin(["", "+", "-"])]
    if len(bad):
        raise InputError(f"events: strand must be + or -; found {_examples(bad.strand)} "
                         f"(e.g. event {bad.event_id.iloc[0]})")
    if event_ids is not None:
        e = e[e.event_id.isin(set(map(str, event_ids)))]
    return e.set_index("event_id")


def _patients(df: pd.DataFrame, samples: pd.DataFrame, table: str) -> pd.DataFrame:
    """patient_id from sample_id when a table is keyed by sample."""
    if "patient_id" in df.columns or "sample_id" not in df.columns:
        return df
    m = samples.set_index("sample_id").patient_id
    pid = df["sample_id"].astype(str).str.strip().map(m)
    if pid.isna().any():
        raise InputError(f"{table}: sample_id not in samples: {_examples(df.sample_id[pid.isna()])}")
    return df.assign(patient_id=pid.to_numpy()).drop(columns="sample_id")


def wide_endpoints(columns) -> dict[str, tuple[str, str]]:
    """Detect endpoint column pairs: <EP>.time with <EP>, <EP>.event, <EP>_event, <EP>.status or <EP>_status."""
    cols = list(map(str, columns))
    out = {}
    for c in cols:
        m = TIME_COLUMN.match(c)
        if not m or not m.group("ep"):
            continue
        ep = m.group("ep").rstrip("._ ")
        for cand in (ep, f"{ep}.event", f"{ep}_event", f"{ep}.status", f"{ep}_status", f"{ep}.Event"):
            if cand in cols and cand != c:
                out[ep] = (c, cand)
                break
    return out


def normalise_survival(df: pd.DataFrame, samples: pd.DataFrame, columns=None, endpoints=None) -> tuple[pd.DataFrame,
                                                                                                        dict]:
    """Long table (patient_id, endpoint, time, event) of valid rows, and the dropped counts per endpoint."""
    df = _patients(_map_columns(df, columns), samples, "survival")
    if not {"endpoint", "time", "event"} <= set(df.columns):
        pairs = endpoints or wide_endpoints(df.columns)
        if not pairs:
            raise InputError("survival: expected columns patient_id, endpoint, time, event (long), or endpoint "
                             "column pairs such as OS.time and OS (wide)")
        _require(df, "survival", ["patient_id"])
        df = pd.concat([pd.DataFrame({"patient_id": df["patient_id"], "endpoint": ep, "time": df[t], "event": df[e]})
                        for ep, (t, e) in pairs.items()], ignore_index=True)
    _require(df, "survival", ["patient_id", "endpoint", "time", "event"])
    s = pd.DataFrame({"patient_id": _ids(df["patient_id"], "survival", "patient_id").to_numpy(),
                      "endpoint": _ids(df["endpoint"], "survival", "endpoint").to_numpy(),
                      "time": _num(df["time"], "survival", "time").to_numpy(),
                      "event": _num(df["event"], "survival", "event").to_numpy()})
    ok = s.event.isin([0, 1]) & np.isfinite(s.time) & (s.time > 0)
    dropped = {str(k): int(v) for k, v in s[~ok].groupby("endpoint").size().items()}
    s = s[ok].drop_duplicates()
    s["event"] = s.event.astype(int)
    dup = s[s.duplicated(["patient_id", "endpoint"], keep=False)]
    if len(dup):
        raise InputError(f"survival: one row per patient and endpoint; conflicting rows for "
                         f"{_examples(dup.patient_id + ' / ' + dup.endpoint)}")
    return s.reset_index(drop=True), dropped


def normalise_clinical(df: pd.DataFrame, samples: pd.DataFrame, columns=None) -> pd.DataFrame:
    """Patient-level covariates indexed by patient_id (identical repeated rows are merged)."""
    df = _patients(_map_columns(df, columns), samples, "clinical")
    _require(df, "clinical", ["patient_id"])
    c = df.copy()
    c["patient_id"] = _ids(c["patient_id"], "clinical", "patient_id")
    c = c.drop_duplicates()
    dup = c.patient_id[c.patient_id.duplicated()]
    if len(dup):
        raise InputError(f"clinical: one row per patient; conflicting rows for {_examples(dup)}")
    return c.set_index("patient_id")


def normalise_pairs(df: pd.DataFrame, samples: pd.DataFrame, columns=None) -> pd.DataFrame:
    df = _map_columns(df, columns)
    cols = ["cohort", "patient_id", "case_sample", "reference_sample"]
    _require(df, "pairs", cols)
    p = pd.DataFrame({c: _ids(df[c], "pairs", c).to_numpy() for c in cols})
    idx = samples.set_index("sample_id")
    problems = []
    for col, role in (("case_sample", "case"), ("reference_sample", "reference")):
        missing = ~p[col].isin(idx.index)
        if missing.any():
            problems.append(f"{col} not in samples: {_examples(p[col][missing])}")
            continue
        ref = idx.loc[p[col]]
        for what, got, want in (("group", ref.role.to_numpy(), role), ("cohort", ref.cohort.to_numpy(), p.cohort),
                                ("patient_id", ref.patient_id.to_numpy(), p.patient_id)):
            bad = np.asarray(got != np.asarray(want))
            if bad.any():
                problems.append(f"{col} with a different {what} in samples: {_examples(p[col][bad])}")
    if problems:
        raise InputError("pairs: " + "; ".join(problems))
    dup = p[p.duplicated(["cohort", "patient_id"], keep=False)]
    if len(dup):
        raise InputError(f"pairs: one pair per patient and cohort; repeated: {_examples(dup.patient_id)}")
    for col in ("case_sample", "reference_sample"):
        d = p[col][p[col].duplicated()]
        if len(d):
            raise InputError(f"pairs: a sample can be in one pair only; repeated {col}: {_examples(d)}")
    return p


def derive_pairs(samples: pd.DataFrame) -> pd.DataFrame:
    """A pair for every patient with exactly one case and one reference sample in a cohort; more is an error."""
    rows, ambiguous = [], []
    s = samples[samples.role.isin(["case", "reference"])]
    for (cohort, pid), d in s.groupby(["cohort", "patient_id"], sort=False):
        c = d.sample_id[d.role.eq("case")].tolist()
        r = d.sample_id[d.role.eq("reference")].tolist()
        if c and r:
            if len(c) == 1 and len(r) == 1:
                rows.append((cohort, pid, c[0], r[0]))
            else:
                ambiguous.append(f"{pid} in {cohort} ({len(c)} case, {len(r)} reference)")
    if ambiguous:
        raise InputError(f"samples: {len(ambiguous)} patient(s) cannot be paired unambiguously, e.g. "
                         f"{_examples(ambiguous)}. Supply a pairs table (cohort, patient_id, case_sample, "
                         "reference_sample).")
    return pd.DataFrame(rows, columns=["cohort", "patient_id", "case_sample", "reference_sample"])


def normalise_expression(df: pd.DataFrame, samples: pd.DataFrame, genes=None, columns=None) -> pd.DataFrame:
    """Wide float matrix: one row per gene, one column per sample of the samples table (its order)."""
    df = _map_columns(df, columns)
    if "gene" not in df.columns and df.index.name == "gene":
        df = df.reset_index()
    if {"gene", "sample_id", "value"} <= set(df.columns):
        d = pd.DataFrame({"gene": _ids(df["gene"], "expression", "gene").to_numpy(),
                          "sample_id": _ids(df["sample_id"], "expression", "sample_id").to_numpy(),
                          "value": _num(df["value"], "expression", "value").to_numpy()})
        if genes is not None:
            d = d[d.gene.isin(set(genes))]
        dup = d[d.duplicated(["gene", "sample_id"])]
        if len(dup):
            raise InputError(f"expression: repeated (gene, sample_id) row(s), e.g. {dup.gene.iloc[0]} / "
                             f"{dup.sample_id.iloc[0]}")
        wide = d.pivot(index="gene", columns="sample_id", values="value")
    elif "gene" in df.columns:
        w = df.copy()
        w["gene"] = _ids(w["gene"], "expression", "gene")
        if genes is not None:
            w = w[w.gene.isin(set(genes))]
        dup = w.gene[w.gene.duplicated()]
        if len(dup):
            raise InputError(f"expression: gene must be unique in a wide table; repeated: {_examples(dup)}")
        wide = w.set_index("gene")
        wide.columns = [str(c).strip() for c in wide.columns]
        wide = wide.apply(lambda c: _num(c, "expression", c.name))
    else:
        raise InputError("expression: expected columns gene, sample_id, value (long) or gene plus one column per "
                         "sample (wide)")
    extra = [c for c in wide.columns if c not in set(samples.sample_id)]
    if len(wide.columns) and len(extra) == len(wide.columns):
        raise InputError(f"expression: no sample matches samples.sample_id (e.g. {_examples(extra)})")
    wide = wide.reindex(columns=samples.sample_id.to_numpy()).astype(float)
    if np.isinf(wide.to_numpy()).any():
        raise InputError("expression: infinite values")
    wide.index.name, wide.columns.name = "gene", "sample_id"
    return wide


# ============================================================================================ the dataset
@dataclass
class Dataset:
    """Validated inputs. Build it with `Dataset.from_tables(...)` or `Dataset.from_dir(folder)`."""

    samples: pd.DataFrame
    psi: pd.DataFrame
    events: pd.DataFrame
    survival: pd.DataFrame | None = None
    pairs: pd.DataFrame | None = None
    expression: pd.DataFrame | None = None
    clinical: pd.DataFrame | None = None
    labels: dict = field(default_factory=lambda: dict(case="Tumour", reference="Normal"))
    notes: dict = field(default_factory=dict)

    @classmethod
    def from_tables(cls, samples, psi=None, events=None, survival=None, pairs=None, expression=None, clinical=None,
                    event_ids=None, case=None, reference=None, columns: dict | None = None,
                    endpoints: dict | None = None, na_values=None, keep=None, keep_column: str | None = None,
                    where=None):
        """Check and normalise the tables (DataFrames or paths).

        event_ids   keep only these events (saves memory on large PSI tables)
        case, reference   the group values compared (a string or a list of strings each); default: the samples
                    table's role column, else tumour and normal
        columns     logical name -> the tables' own column name, e.g. {"cohort": "cancer", "sample_id": "File.ID"}
        endpoints   explicit wide survival columns, e.g. {"OS": ("OS.time", "OS")}
        na_values   extra codes to read as missing in every table, e.g. ["missing", "[Not Available]"]
        keep        only the patients listed (IDs, or a table and keep_column): see the module docstring
        where       only the patients meeting conditions, e.g. ["Subtype=Her2"] or {"Subtype": "Her2"}: see the module
                    docstring
        """
        def load(x):
            if x is None:
                return None
            return with_missing(x, na_values) if isinstance(x, pd.DataFrame) else read_table(x, columns, na_values)
        raw_samples = load(samples)
        subset, id_map = None, None
        if keep is not None or where:                       # only the selected patients (and all their samples)
            raw_m = _map_columns(raw_samples, columns)
            _require(raw_m, "samples", ["sample_id", "patient_id"])
            # every sample's patient, before the subset: tables keyed by sample_id may list samples left out
            id_map = pd.DataFrame({c: _ids(raw_m[c], "samples", c).to_numpy() for c in ("sample_id", "patient_id")})
            dup = id_map.sample_id[id_map.sample_id.duplicated()]
            if len(dup):
                raise InputError(f"samples: sample_id must be unique; repeated: {_examples(dup)}")
            mask, by, name, listed = np.ones(len(raw_m), bool), [], [], {}
            if keep is not None:
                ids, head = _keep_list(keep, keep_column)
                if head and _keep_mask(raw_m, [head]).any():   # a list without a header: its first line is an ID
                    warnings.warn(f"keep: {keep} has no header; its first line ({head}) is read as an ID",
                                  stacklevel=2)
                    ids = [head] + ids
                if not ids:
                    raise InputError(f"{keep}: the keep list is empty")
                m = _keep_mask(raw_m, ids)
                if not m.any():
                    raise InputError(f"keep: none of the {len(ids)} listed IDs (e.g. {_examples(ids)}) is a patient or "
                                     f"sample of the samples table (e.g. {_examples(raw_m.patient_id)})")
                mask &= m
                src = str(keep) if isinstance(keep, (str, Path)) else ""
                listed = dict(listed=len(ids), source=src or "a list")
                by.append(f"the {len(ids)} IDs in `{src}`" if src else f"a list of {len(ids)} IDs")
                name.append(Path(src).stem if src else f"{len(ids)} IDs")
            if where:
                conds = where_conditions(where)
                clin = None if clinical is None else _patients(_map_columns(load(clinical), columns), id_map,
                                                                "clinical")
                mask &= _where_mask(raw_m, clin, conds)
                by.append(" and ".join(f"{c} = {' or '.join(v)}" for c, v in conds))
                name.append(", ".join(f"{c}={'/'.join(v)}" for c, v in conds))
            if not mask.any():
                raise InputError(f"subset: no patient is selected by {' and '.join(by)}")
            subset = dict(patients=int(raw_m.patient_id[mask].nunique()), samples=int(mask.sum()),
                          of_patients=int(raw_m.patient_id.nunique()), by=" and ".join(by), name=", ".join(name),
                          all_cohorts=sorted(set(raw_m["cohort"].dropna().astype(str).str.strip()))
                          if "cohort" in raw_m.columns else [], **listed)
            raw_samples = raw_samples[mask].reset_index(drop=True)
        s, meta = normalise_samples(raw_samples, case, reference, columns)
        sources = {}
        P, E = load(psi), load(events)
        if P is None and E is None:
            raise InputError("no psi table was given")
        if P is None or E is None:                     # one table holds the PSI values and the event columns
            which = "psi" if P is not None else "events"
            P, E = split_psi_events(P if P is not None else E, s, columns, which)
            sources["events"] = f"{which} table"
        ev = normalise_events(E, event_ids, columns)
        p = normalise_psi(P, s, event_ids, columns, subset=subset is not None,
                          signed=ev.index[ev.event_type.astype(str).str.upper().eq("HIT")])
        if event_ids is not None:
            want = list(dict.fromkeys(map(str, event_ids)))
            miss = [e for e in want if e not in p.index or e not in ev.index]
            if miss:
                raise InputError(f"event(s) missing from the psi or events table: {_examples(miss)}")
            p = p.loc[want]
        no_row = [e for e in p.index if e not in ev.index]
        if no_row:
            raise InputError(f"events: {len(no_row)} event(s) in psi have no events row (e.g. {_examples(no_row)})")
        ev = ev.loc[p.index]
        own = _from_samples(_map_columns(raw_samples, columns).reset_index(drop=True), s,
                            endpoints if survival is None else None)
        by_sample = s if id_map is None else id_map          # maps the sample_id of sample-keyed tables to patients
        if survival is not None:
            sv, dropped = normalise_survival(load(survival), by_sample, columns, endpoints)
        elif own["survival"] is not None:
            sv, dropped = normalise_survival(own["survival"], s)
            sources["survival"] = f"samples table ({', '.join(sorted(sv.endpoint.unique()))})"
        else:
            sv, dropped = None, {}
        if pairs is not None:
            P2 = _map_columns(load(pairs), columns)
            if subset is not None and {"case_sample", "reference_sample"} <= set(P2.columns):
                kept = set(s.sample_id)
                P2 = P2[P2.case_sample.astype(str).str.strip().isin(kept)
                        & P2.reference_sample.astype(str).str.strip().isin(kept)]
            pr, pairs_source = normalise_pairs(P2, s), "table"
        elif own["pairs"] is not None:
            pr, pairs_source = normalise_pairs(own["pairs"], s), "samples.pair_id"
        else:
            pr, pairs_source = derive_pairs(s), "derived"
        genes = None if event_ids is None else set(ev.expression_gene)
        ex = None if expression is None else normalise_expression(load(expression), s, genes, columns)
        if clinical is not None:
            cl = normalise_clinical(load(clinical), by_sample, columns)
        elif own["clinical"] is not None:
            try:
                cl = normalise_clinical(own["clinical"], s)
            except InputError as err:
                raise InputError(f"samples table, clinical columns: {err}") from None
            sources["clinical"] = f"samples table ({', '.join(map(str, cl.columns[:6]))}" + (
                ", ..." if len(cl.columns) > 6 else "") + ")"
        else:
            cl = None
        notes = dict(survival_rows_dropped=dropped, pairs_source=pairs_source, other_groups=meta["other_groups"],
                     sources=sources, subset=subset)
        return cls(samples=s, psi=p, events=ev, survival=sv, pairs=pr, expression=ex, clinical=cl,
                   labels=meta["labels"], notes=notes)

    @classmethod
    def from_dir(cls, folder, event_ids=None, case=None, reference=None, columns=None, endpoints=None,
                 na_values=None, keep=None, keep_column=None, where=None, **override):
        """Read <table>.<csv|tsv|txt|parquet> files from a folder; keyword arguments replace single tables."""
        check_table_names(override)
        found = find_tables(folder)
        if "samples" not in found and override.get("samples") is None:
            raise InputError(f"{folder}: no samples table (expected e.g. samples.csv)")
        if not ({"psi", "events"} & set(found)) and override.get("psi") is None and override.get("events") is None:
            raise InputError(f"{folder}: no psi table (expected e.g. psi.csv, with or without the event columns)")
        tables = {k: override.get(k) if override.get(k) is not None else found.get(k) for k in TABLES}
        return cls.from_tables(event_ids=event_ids, case=case, reference=reference, columns=columns,
                               endpoints=endpoints, na_values=na_values, keep=keep, keep_column=keep_column,
                               where=where, **tables)

    # ------------------------------------------------------------------ summaries
    @property
    def cohorts(self) -> list[str]:
        return sorted(self.samples.cohort.unique())

    @property
    def endpoints(self) -> list[str]:
        return [] if self.survival is None else sorted(self.survival.endpoint.unique())

    @property
    def has_reference(self) -> bool:
        return bool((self.samples.role == "reference").any())

    def summary(self) -> pd.DataFrame:
        """Per cohort: case and reference samples, pairs, survival samples, and patients and events per endpoint."""
        s = self.samples
        out = s.groupby("cohort").agg(case=("role", lambda x: int((x == "case").sum())),
                                      reference=("role", lambda x: int((x == "reference").sum())),
                                      survival_samples=("survival_cohort", "sum"))
        out = out.rename(columns={"case": self.labels["case"], "reference": self.labels["reference"]})
        out["pairs"] = self.pairs.groupby("cohort").size().reindex(out.index).fillna(0).astype(int)
        base = s[s.survival_cohort]
        for ep in self.endpoints:
            sv = self.survival[self.survival.endpoint.eq(ep)].set_index("patient_id")
            b = base.assign(ev=base.patient_id.map(sv.event))
            out[f"{ep}_patients"] = b.groupby("cohort").ev.count().reindex(out.index).fillna(0).astype(int)
            out[f"{ep}_events"] = b.groupby("cohort").ev.sum().reindex(out.index).fillna(0).astype(int)
        return out.reset_index()
