"""Default clinical adjustment: find age, sex and stage in the clinical table, clean them, and build the model.

Detection (column names, case-insensitive; the first match wins):
    age     age, age_years, age_at_diagnosis, age_at_index, age_at_initial_pathologic_diagnosis, diagnosis_age,
            or any column named age_* / *_age
    sex     sex, gender
    stage   stage, ajcc_pathologic_stage, pathologic_stage, ajcc_pathologic_tumor_stage, tumor_stage, clinical_stage,
            figo_stage, or any other column containing "stage"
Cleaning (only for the detected columns):
    age     numeric; anything else is missing (units do not matter: age enters per SD)
    sex     female / male ('f', 'woman', 'm', 'man' accepted); anything else is missing
    stage   the Roman numeral of the overall stage: 'Stage IIA' -> II, 'IIIc' -> III, 'Stage IV' -> IV, 'Stage 0' -> 0;
            anything else ('missing', '[Not Available]', 'unknown', 'X', ...) is missing
Model: PSI + host expression (when an expression table is given) + age + sex + stage, stage against I and sex
against female. A variable recorded for too few patients of a cohort is left out of that cohort's model (Settings.
covariate_min_complete), and the table says so.
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
              "clinical_stage", "figo_stage", "ajcc_clinical_stage"),
}
STAGE = re.compile(r"^(?:stage\s*)?(IV|III|II|I|0)(?:[A-C]\d?|\d)?$", re.IGNORECASE)
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
            cands = [c for lc, c in low.items() if "stage" in lc]
            hit = sorted(cands, key=lambda c: (0 if "patholog" in c.lower() else 1, len(c)))[0] if cands else None
        if hit is not None:
            out[role] = hit
    return out


def clean_stage(v) -> float | str:
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return np.nan
    m = STAGE.match(str(v).strip())
    return m.group(1).upper() if m else np.nan


def clean_sex(v) -> float | str:
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return np.nan
    return SEX.get(str(v).strip().lower(), np.nan)


def auto_clinical(ds: Dataset) -> tuple[Dataset, CoxModel | None, dict]:
    """The dataset with cleaned age / sex / stage columns, the default adjusted model (None without a clinical table
    or without any of the three), and what was found."""
    if ds.clinical is None:
        return ds, None, {}
    found = detect(ds.clinical.columns)
    if not found:
        return ds, None, {}
    cl = ds.clinical.copy()
    for role, col in found.items():
        v = cl[col]
        if role == "age":
            cl[role] = pd.to_numeric(v, errors="coerce")
        elif role == "sex":
            cl[role] = v.map(clean_sex)
        else:
            cl[role] = v.map(clean_stage)
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
