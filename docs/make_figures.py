"""Regenerate the figures in docs/ from the synthetic example (no real data is used).

    python docs/make_figures.py

It writes docs/example_panel.png (the README figure), the reading guide's crops in docs/reading/ and the probe overview
used by the guide. The crops are pixel boxes of two pages drawn at 200 dpi (an assay page and the gene's expression
page); the script stops when either page changes size, because the boxes then need adjusting (blank horizontal bands
separate the pages' sections).
"""
from __future__ import annotations

import shutil
import tempfile
from pathlib import Path

from PIL import Image

import splice_assay as sa
from splice_assay import example
from splice_assay.clinical import auto_clinical
from splice_assay.probe import probe

DOCS = Path(__file__).resolve().parent
OUT = DOCS / "reading"
PAGE_SIZE = (1440, 2567)                    # the guide page at 200 dpi, in pixels
CROPS = {                                   # name: (top, bottom) in pixels of the guide page
    "01_schematic": (0, 300),
    "02_protein": (300, 640),
    "03_cohort_row": (640, 1116),
    "04_cohort_without_normals": (1588, 2060),
    "05_forest": (2130, 2388),
    "06_legend_footnote": (2400, 2567),
}
EXPRESSION_SIZE = (1440, 2014)              # the gene's expression page for the same cohorts
EXPRESSION_CROPS = {"07_expression_page": (0, 925)}   # title, forest and the first cohort


def main() -> None:
    OUT.mkdir(exist_ok=True)
    s = sa.Settings(formats=("png",), dpi=200)
    with tempfile.TemporaryDirectory() as tmp:
        example.write(tmp)
        data, gtf, proteins = Path(tmp) / "data", Path(tmp) / "data" / "annotation.gtf", Path(tmp) / "proteins"
        ds = sa.Dataset.from_dir(data)
        # the README figure: two cohorts, the model band below
        res = sa.analyze(ds)
        model = sa.CoxModel().with_clinical(("age", "sex", "stage"), baseline={"stage": "I", "sex": "female"})
        sa.event_panel(ds, "SYN1:SE:1", ["COH1", "COH2"], "OS", gtf=gtf, results=res, detail=model,
                       proteins=proteins, out_dir=tmp, stem="example_panel", settings=s)
        shutil.copy(Path(tmp) / "example_panel.png", DOCS / "example_panel.png")
        # the guide page: three cohorts as a probe draws them (pairs; too few pairs; no normal samples)
        dsa, adjusted, _ = auto_clinical(ds)
        gene = list(dsa.events.index[dsa.events.gene.eq("SYN1")])
        res = sa.analyze(dsa, events=gene, endpoints=["OS"], settings=s)
        sa.event_panel(dsa, "SYN1:SE:1", ["COH1", "COH2", "COH4"], "OS", gtf=gtf, results=res, detail=adjusted,
                       proteins=proteins, out_dir=tmp, stem="guide_page", settings=s)
        sa.expression_panel(dsa, "SYN1", ["COH1", "COH2", "COH4"], "OS", model=adjusted, out_dir=tmp,
                            stem="guide_expression", settings=s)
        for stem, size, crops in (("guide_page", PAGE_SIZE, CROPS),
                                  ("guide_expression", EXPRESSION_SIZE, EXPRESSION_CROPS)):
            page = Image.open(Path(tmp) / f"{stem}.png")
            if page.size != size:
                raise SystemExit(f"{stem} is {page.size} px, not {size}: adjust the crops in {__file__}")
            for name, (top, bottom) in crops.items():
                page.crop((0, top, page.size[0], bottom)).save(OUT / f"{name}.png", optimize=True)
        # the probe's overview, as the guide shows it
        r = probe(ds, genes=["SYN1"], settings=s, gtf=gtf, out_dir=Path(tmp) / "probe", max_pages=1,
                  log=lambda *_: None)
        shutil.copy(r.paths["overview"], OUT / "08_overview.png")
    print(f"wrote {DOCS / 'example_panel.png'} and {len(CROPS) + len(EXPRESSION_CROPS) + 1} figures in {OUT}")


if __name__ == "__main__":
    main()
