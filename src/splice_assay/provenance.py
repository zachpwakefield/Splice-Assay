"""Provenance: what produced an output. Every figure and table can be traced to the exact input values (SHA-256 of a
canonical serialisation of the rows used), the settings, and the software versions."""
from __future__ import annotations

import hashlib
import io
import platform
import sys
from importlib import metadata

import pandas as pd

from ._version import __version__

LIBRARIES = ("numpy", "pandas", "scipy", "matplotlib", "lifelines")


def table_sha256(df: pd.DataFrame | None) -> str | None:
    """SHA-256 of a table written as CSV with every float in round-trip precision (index included)."""
    if df is None:
        return None
    buf = io.StringIO()
    df.to_csv(buf, float_format="%.17g", lineterminator="\n")
    return hashlib.sha256(buf.getvalue().encode("utf-8")).hexdigest()


def environment() -> dict:
    libs = {}
    for name in LIBRARIES:
        try:
            libs[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            libs[name] = None
    return dict(splice_assay=__version__, python=platform.python_version(), libraries=libs,
                platform=sys.platform)


def record(what: str, settings, inputs: dict, call: dict) -> dict:
    """inputs: name -> DataFrame of the rows used (hashed); call: the arguments that selected them."""
    return dict(output=what, call=call, settings=settings.to_dict(),
                inputs_sha256={k: table_sha256(v) for k, v in inputs.items()},
                input_rows={k: (None if v is None else int(len(v))) for k, v in inputs.items()},
                environment=environment())
