# Changelog

## 0.1.0 (unreleased)

First version, released under GPL-3.0-only (as SpliceImpactR, whose transcript matcher it ports).
- **Input contract.** Plain tables (CSV, TSV or Parquet): samples, psi (long or wide), events, survival (long or wide,
  e.g. `OS.time` + `OS`), and optional pairs, expression and clinical.
  - Your own column names are mapped (`--column cohort=cancer`).
  - Checks give actionable error messages.
- **Probe.** `splice-assay probe DATA --gene G` runs every event × cohort, fits base and adjusted models, ranks the
  events, and writes one assay page per event, an overview, `probe.pdf` and `report.md`.
- **Defaults.**
  - `panel` needs only `--event`: it picks the most promising cohorts, uses OS, and adds models adjusted for age, sex
    and stage, found and cleaned automatically.
  - The GTF can come from `$SPLICE_ASSAY_GTF`.
- **Agent guide.** `AGENT_GUIDE.md` is an operating manual for agents.
- **Missing-value codes.** `--na-value missing` (repeatable) reads such codes as missing in every table.
- **Groups.** Any two groups can be compared (`--case`, `--reference`; default tumour and normal). Figures and tables
  use their names. Data without a reference group are supported.
- **Statistics.**
  - Case vs reference per cohort: paired exact signed-rank and unpaired Mann–Whitney, with a hit rule that needs
    robustness to PSI 0/1. Within-patient support and composition diagnostics are also reported.
  - Survival per cohort and endpoint: median-split KM and log-rank.
  - Cox on PSI adjusted for host expression and any clinical covariates (numeric, categorical, strata), with gates,
    failure rules and every model term reported.
- **Event panel.**
  - Schematic: a collapsed gene model from a GTF (windowed for long genes) and nested snoRNAs.
  - Group view: matched pairs beside all samples.
  - KM with an at-risk table.
  - Forest across cohorts.
- **Cox model figure.** Every term of one cell's model.
- **Model band.** `event_panel(detail=...)` or `panel --detail` puts the model of each cohort shown into the event
  panel: one assay figure per event.
- **Baseline levels.** `--baseline stage=I` sets the reference level of categorical covariates.
- **Rare categories.** A level with fewer than `level_min_patients` (10) patients in a cohort, or no events, joins
  its neighbour before the fit (stage I into "I–II"), so the adjusted model no longer fails where stage I is rare
  (all BLCA and SKCM fits in the EHMT2 probe).
- **Flags on Cox fits.** Notes only, with the estimates unchanged:
  - fewer than `cox_events_per_term` (10) events per model term (an overfit risk);
  - a PSI spread below `narrow_psi_below` (0.05), measured by IQR or SD (`narrow_psi_measure`).

  They appear in `cox_notes`, the page's model header (which now wraps) and a section of the probe report.
- **`--where COLUMN=VALUE[,VALUE]`.** Selects patients by a clinical or samples column (case-insensitive;
  conditions combine with each other and with `--keep`).
- **Several endpoints in one probe.** `--endpoint OS --endpoint DSS` (or `all`), one folder per endpoint.
- **Proportional-hazards notes.**
  - Every Cox term, and the KM high/low split, gets a Schoenfeld test (`ph_p`, `cox_terms.ph_p`, `km_ph_p`).
  - A p below `ph_note_below` (0.05) adds a note to the model header, the KM header, `cox_notes` and `km_notes`,
    and to the probe report's "Notes on the survival tests". No result is removed.
- **Guides.** [Reading the results, panel by panel](docs/reading-results.md) and
  [getting the annotation cache from SpliceImpactR](docs/annotation-cache.md).
  - `docs/make_figures.py` regenerates their figures and the README figure.
- **`protein-cache`** reads any `<name>.gtf.rds` with its `<name>_sequences.rds` (any species or release), and works
  without feature files.
- **Long pages are split.** More than `cohorts_per_page` (6) cohorts are drawn over balanced pages, each with the
  full forest (a `--top all` page of 27 cohorts was 140 inches tall).
- **Stacked layout.** For three or more cohorts, one row per cohort holds its comparison, KM and Cox model; the
  forest moves below.
- **Outputs.** SVG, PDF and PNG, a CSV of every plotted value, and a provenance JSON. Outputs are byte-identical across
  runs.
- **rMATS import.** Events with geometry, and PSI from IncLevel.
- **Subsets.** `--keep FILE` runs any command on the listed patients only. TCGA barcodes are matched to their
  patients, and normals are kept.
  - Results tables keep their full column set when nothing was tested.
  - The probe ranks by KM when no Cox model can be fitted.
- **Condition names from the samples table.** A `role` column (case / reference) names the comparison without
  flags. Group names keep inner capitals in sentences, and the probe report uses them too.
- **q values within each gene.** Benjamini–Hochberg per kind of test (within patients, all samples, KM, Cox PSI
  term; per endpoint and model) over all of the gene's events × cohorts.
  - Shown beside p on the pages, with a footnote naming the family.
  - No q for families of fewer than 10 tests (`fdr_min_family`).
  - `panel` and `panels` now analyse every event of the gene so that the family is complete.
- **Two tables are enough.**
  - `samples` may carry survival (wide `OS.time` + `OS`, …), clinical columns and `pair_id`.
  - `psi` may carry the event columns.
  - Separate tables still work and take precedence.
- **Host-gene expression rows.** With an expression table, each page shows the host gene's own tumour vs normal,
  KM on its median split, and Cox on expression + the clinical terms (HR per SD) for every cohort shown.
  - `analyse_expression` computes the same statistics; the probe writes `expression_cells.csv`.
  - `--no-gex` turns them off.
- **`--top all`.** Probe pages and `panel` can show every cohort with a test. Long titles and file names are
  shortened.
- **Protein consequences.** With a protein cache (`splice-assay protein-cache`, from a SpliceImpactR annotation
  cache), each event gets the annotated transcripts of its two forms and a suggested protein change.
  - Matching ports SpliceImpactR's ranked-pair matcher, with a fallback to the event's own exons.
  - The change: residues, frame, and features gained or lost.
  - Where it appears: a band on every page, `proteins.csv` and a Protein column in probes, and the `proteins`
    command.
  - See docs/proteins.md.
- **Command line.** `validate`, `analyse`, `panel`, `panels`, `probe`, `cox`, `proteins`, `protein-cache`, `example`,
  `import-rmats` and `gtf-subset`.
- **Synthetic example.** A simulated dataset with a small GTF and clinical covariates.
