"""Pieces shared by the event panel and the Cox model figure."""
from __future__ import annotations

import re
from dataclasses import dataclass

import pandas as pd

ENDPOINT_NAMES = {"OS": "Overall survival", "DSS": "Disease-specific survival", "PFI": "Progression-free interval",
                  "DFI": "Disease-free interval", "PFS": "Progression-free survival", "DFS": "Disease-free survival",
                  "RFS": "Relapse-free survival", "EFS": "Event-free survival"}


@dataclass
class Panel:
    figure: object          # matplotlib.figure.Figure
    table: pd.DataFrame     # every plotted value (the CSV)
    paths: dict             # written files by format
    provenance: dict
    stem: str


def page_parts(cohorts, n: int) -> list[list]:
    """Cohorts split into balanced pages of at most n each (n = 0: one page): 27 by 6 gives 6, 6, 5, 5, 5."""
    cohorts = list(cohorts)
    if n <= 0 or len(cohorts) <= n:
        return [cohorts]
    k = -(-len(cohorts) // n)
    size, extra = divmod(len(cohorts), k)
    out, i = [], 0
    for j in range(k):
        m = size + (j < extra)
        out.append(cohorts[i:i + m])
        i += m
    return out


def safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(text)).strip("_")
