"""Default clinical adjustment: find age, sex and stage in the clinical table, clean them, and build the model.

Detection (column names, case-insensitive; the first match wins):
    age     age, age_years, age_at_diagnosis, age_at_index, age_at_initial_pathologic_diagnosis, diagnosis_age,
            or any column named age_* / *_age
    sex     sex, gender
    stage   stage, ajcc_pathologic_stage, pathologic_stage, ajcc_pathologic_tumor_stage, tumor_stage, clinical_stage,
            figo_stage, ajcc_clinical_stage, overall_stage, or any other column containing "stage" that is not a T,
            N or M stage (a t / n / m part of the name, e.g. pathologic_t_stage) or a summary stage (SEER)
Cleaning (only for the detected columns):
    age     numeric; anything else is missing (units do not matter: age enters per SD)
    sex     female / male ('f', 'woman', 'm', 'man' accepted); anything else is missing
    stage   the Roman numeral of the overall stage: 'Stage IIA' -> II, 'IIIc' -> III, 'Stage IV' -> IV, 'Stage 0' -> 0;
            in the named columns above, stage numbers too: 2, 2.0, '02', '2B', 'Stage 4' -> II; anything else
            ('missing', '[Not Available]', 'unknown', 'X', 5, ...) is missing
A column the base model already uses (`keep`) is left as it is: the base model's own columns are never rewritten.
Model: PSI + host expression (when an expression table is given) + age + sex + stage, stage against I and sex
against female. A variable recorded for too few patients of a cohort is left out of that cohort's model (Settings.
covariate_min_complete), and the table says so. `adjusted` builds it on top of a base model, so that it keeps the base
model's terms and choices (e.g. no host expression, or strata).
"""
from __future__ import annotations

import re
from dataclasses import replace

import numpy as np
import pandas as pd

from .dataset import Dataset
from .stats.survival import CoxModel

NAMES = {
    "age": ("age", "age_years", "age_at_diagnosis", "age_at_index", "age_at_initial_pathologic_diagnosis",
            "diagnosis_age", "age_at_diagnosis_years"),
    "sex": ("sex", "gender"),
    "stage": ("stage", "ajcc_pathologic_stage", "pathologic_stage", "ajcc_pathologic_tumor_stage", "tumor_stage",
              "clinical_stage", "figo_stage", "ajcc_clinical_stage", "overall_stage"),
}
STAGE = re.compile(r"^(?:stage\s*)?(IV|III|II|I|0)(?:[A-C]\d?)?$", re.IGNORECASE)
NOT_OVERALL = re.compile(r"(^|_)[tnm](_|$)|summary|seer", re.IGNORECASE)      # a T, N or M stage, or a summary stage
STAGE_NUMBER = re.compile(r"^(?:stage\s*)?([0-4])(?:[A-C]\d?)?$", re.IGNORECASE)      # 2, 2B, Stage 4
ROMAN = {0: "0", 1: "I", 2: "II", 3: "III", 4: "IV"}
SEX = {"female": "female", "f": "female", "woman": "female", "male": "male", "m": "male", "man": "male"}


def detect(columns) -> dict[str, str]:
    """{'age': column, 'sex': column, 'stage': column} for the columns found."""
    low = {str(c).lower(): str(c) for c in columns}
    out = {}
    for role, names in NAMES.items():
        hit = next((low[n] for n in names if n in low), None)
        if hit is None and role == "age":
            hit = next((c for lc, c in low.items() if re.search(r"(^|_)age(_|$)", lc)), None)
        if hit is None and role == "stage":
            cands = [c for lc, c in low.items() if "stage" in lc and not NOT_OVERALL.search(lc)]
            hit = sorted(cands, key=lambda c: (0 if "patholog" in c.lower() else 1, len(c)))[0] if cands else None
        if hit is not None:
            out[role] = hit
    return out


def _stage_number(x: float):
    return ROMAN.get(int(x), np.nan) if np.isfinite(x) and float(x).is_integer() else np.nan


def clean_stage(v, numbers: bool = True) -> float | str:
    """The overall stage as a Roman numeral (0, I-IV), from Roman numerals ('Stage IIA' -> II, 'IIIc' -> III) or, with
    `numbers`, from stage numbers (2, 2.0, '02', '2B', 'Stage 4' -> II, II, II, II, IV); anything else is missing."""
    if v is None or isinstance(v, bool) or (isinstance(v, float) and np.isnan(v)):
        return np.nan
    if isinstance(v, (int, float, np.integer, np.floating)):              # a numeric column
        return _stage_number(float(v)) if numbers else np.nan
    text = str(v).strip()
    try:
        x = float(text)                                                     # '2', '2.0', '02'
    except ValueError:
        x = None
    if x is not None:
        return _stage_number(x) if numbers else np.nan
    m = STAGE.match(text)
    if m:
        return m.group(1).upper()
    m = STAGE_NUMBER.match(text)
    return ROMAN[int(m.group(1))] if m and numbers else np.nan


def clean_sex(v) -> float | str:
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return np.nan
    return SEX.get(str(v).strip().lower(), np.nan)


def auto_clinical(ds: Dataset, keep=()) -> tuple[Dataset, CoxModel | None, dict]:
    """The dataset with cleaned age / sex / stage columns, the default adjusted model (None without a clinical table
    or without any of the three), and what was found (role -> column). `keep`: the columns the base model already
    uses (covariates, strata); a variable named or found in one of them is left as it is and out of this model."""
    if ds.clinical is None:
        return ds, None, {}
    keep = set(keep)
    found = {r: c for r, c in detect(ds.clinical.columns).items() if r not in keep and c not in keep}
    if not found:
        return ds, None, {}
    cl = ds.clinical.copy()
    for role, col in found.items():
        v = cl[col]
        if role == "age":
            cl[role] = pd.to_numeric(v, errors="coerce")
        elif role == "sex":
            cl[role] = v.map(clean_sex)
        else:                                               # stage numbers only from a column named as a stage
            numbers = str(col).lower() in NAMES["stage"]
            cl[role] = v.map(lambda x: clean_stage(x, numbers))
    covs = [r for r in ("age", "sex", "stage") if r in found and cl[r].notna().any()]
    if not covs:
        return ds, None, found
    baseline = {}
    if "stage" in covs and (cl.stage == "I").any():
        baseline["stage"] = "I"
    if "sex" in covs and (cl.sex == "female").any():
        baseline["sex"] = "female"
    model = CoxModel(covariates=tuple(covs), categorical=tuple(c for c in covs if c != "age"), baseline=baseline)
    return replace(ds, clinical=cl), model, found


def found_summary(clinical: pd.DataFrame, found: dict) -> str:
    """What the automatic adjustment read, per variable: its column and the cleaned values, e.g. 'stage =
    ajcc_pathologic_stage (I 40, II 55, III 30, IV 12; 8 missing)', for checking the cleaning."""
    out = []
    for role, col in found.items():
        v = clinical[role]
        miss = int(v.isna().sum())
        if role == "age":
            what = f"{int(v.notna().sum())} recorded"
        else:
            n = v.value_counts()
            order = ["0", "I", "II", "III", "IV"] if role == "stage" else sorted(n.index)
            what = ", ".join(f"{k} {int(n[k])}" for k in order if k in n.index) or "none recorded"
        out.append(f"{role} = {col} ({what}" + (f"; {miss} missing" if miss else "") + ")")
    return "; ".join(out)


def adjusted(base: CoxModel, auto: CoxModel, found: dict) -> CoxModel:
    """The default adjusted model: the base model (its host-expression choice, covariates, strata and baselines) plus
    the clinical terms of `auto` (from `auto_clinical`, with `found`: role -> the column it was found in). A variable
    the base model already has, as a covariate or a stratum, under its role or its column (e.g. gender for sex), is not
    added again; a baseline the base model names wins over the default one."""
    have = set(base.covariates) | set(base.strata)
    covs = tuple(c for c in auto.covariates if c not in have and found.get(c, c) not in have)
    return base.with_clinical(covariates=covs, categorical=tuple(c for c in auto.categorical if c in covs),
                              baseline={**{k: v for k, v in auto.baseline if k in covs}, **dict(base.baseline)})
