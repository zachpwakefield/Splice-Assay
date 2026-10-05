"""splice-assay: case-vs-reference and survival panels for alternative-splicing events, from plain tables.

    import splice_assay as sa
    ds = sa.Dataset.from_dir("my_data/")                 # samples, psi, events, survival, [pairs, expression, clinical]
    res = sa.analyze(ds)                                 # per event x cohort (x endpoint) statistics
    sa.event_panel(ds, "GENE1:SE:1", ["COH1"], "OS", gtf="annotation.gtf.gz", out_dir="figures/")
    sa.cox_model_figure(ds, "GENE1:SE:1", "COH1", "OS", model=sa.CoxModel(covariates=["age", "stage"]))
"""
from ._version import __version__
from .analysis import Results, analyze
from .config import Settings
from .dataset import Dataset, InputError
from .events import Geometry, geometry
from .plot import Panel, cox_model_figure, event_panel, expression_panel, model_terms
from .stats.survival import CoxModel

analyse = analyze                                       # the earlier spelling, still accepted

__all__ = ["CoxModel", "Dataset", "Geometry", "InputError", "Panel", "Results", "Settings", "__version__", "analyze",
           "cox_model_figure", "event_panel", "expression_panel", "geometry", "model_terms"]
