# Changelog

## 0.2.0 (2026-10-04)

- **Naming.** `analyze` is the command and the Python function (`sa.analyze`, `analyze_expression`); the
  earlier spelling `analyse` still works as an alias.
- **Options.**
  - `--ridge clinical|molecular|all` (`cox_ridge`) adds an opt-in ridge penalty λ/2·β² to the clinical terms, to PSI
    and host expression, or to all terms (β per SD, per level for categories); `--ridge-penalty` (`cox_ridge_penalty`,
    default 1) sets λ. Penalized models are named so (tables, pages, the overview's caption), and the probe report
    says which terms the penalty reached. Penalized CIs and p values are approximate, and with `clinical` or `all`
    the shrunk covariates adjust PSI only partly (docs/methods.md).
  - `--min-events-per-term N` (`cox_min_events_per_term`, default off) fits a Cox model only with at least N events per
    estimated term (`too_few_events_per_term`; the probe report counts such fits, and a model row says why it is
    missing).
- **Changes in results.**
  - The PSI (and HIT-index) hazard ratio is shown per SD of the cohort by default (`psi_hr_unit = "sd"`;
    `--hr-unit iqr` for the HR per IQR). The narrow-range note follows the unit (SD below 0.05). The tables hold both
    units, as before.
  - Cox needs 10 events (`cox_min_events`, was 20). A fit with fewer than 20 events (`cox_low_power_events`) runs
    and is marked low power: a note, `cox_low_power` in the tables, ‡ after its forest CI, a line in the probe
    report, ‡ after p in its ranked table, and `best_low_power` in `events.csv`. The probe ranks events by Cox hits
    in fits that are not low power; low-power hits only break ties, before the smallest p (`adj_cox_p05_low_power`,
    `cox_p05_low_power`, "4 (1‡)" in the report), and the best cohort and the cohorts on each page prefer fits that
    are not low power. A cohort with 10–19 deaths cannot push an event ahead of one with more evidence from
    better-powered fits.
  - Host expression follows the rule for clinical variables: recorded for fewer than 80% of a cohort's fit patients
    (none when the expression table has no row for the host gene), or constant, it is left out and the Cox model is
    PSI alone there, with a note ("host expression left out (…)"); such cells were not fitted (`no_host_expression`,
    `constant_expression`). Reports, printouts and the forest name the models as fitted.
  - The automatic clinical adjustment reads stage numbers (1–4, "02", "2B", 2.0) as stages I–IV in a column named as
    an overall stage; they were missing. A T, N or M stage or a summary stage (SEER) is not taken for the stage.
    `validate`, the probe report and `panel` print what it read ("stage = stage (I 254, II 169, III 107)"). It never
    rewrites a column the page's or probe's own model uses, so the forest is the same with or without model rows.
  - The unpaired (all-samples) test compares one value per patient: a patient's several samples in one group (e.g.
    replicate aliquots) enter with their mean. `min_group` counts patients; new columns `unpaired_n_case_samples` and
    `unpaired_n_reference_samples` count the samples behind them. Data with one sample per patient and group give the
    same results as before.
  - Mann–Whitney ranks and the constant check use values rounded to 12 decimals, as every other comparison does.
  - The 0/1 check is `fail`, not `untestable`, when enough values remain but none differ.
  - A most common category without events merges into the next most common one, so its Cox fit no longer fails.
  - The default adjusted model is the base model plus the age, sex and stage found: `--no-expression` now applies to
    it too (`probe` and `panel`), and `panel`'s model rows keep the page's covariates and strata (a variable the page
    already has, e.g. `--covariate gender`, is not added again as `sex`). `panel --detail` keeps the adjustment, and
    `panels --detail` now draws the same model rows as `panel` (it showed the page's model alone).
- **Fixes.**
  - Keep lists are read as text (`00123` no longer selects patient `123`); a list without a header keeps its first
    ID; an empty or unreadable keep file is an input error.
  - Repeated `--where` conditions on one column must all hold, as documented (they were combined into one).
    `--where` matches numbers by value (`grade=1` matches a column read as 1.0), and a missing value never matches.
  - A keep-list ID that is itself a patient or sample ID selects only that patient (`1-2` no longer also selects
    patient `1`); another ID selects the patient of its longest leading part that is one (TCGA barcodes).
  - Survival and clinical tables keyed by `sample_id` work with `--keep` and `--where`.
  - A Cox cohort without a patient complete for every covariate is a gate, not a crash.
  - `import-rmats` skips an event type whose MATS file has no events instead of stopping (a file without even a
    header line is an input error).
  - `import-hitindex`: genes sharing a symbol are ordered and drawn each on its own; `--gene` reads quoted IDs.
  - `panel`, `panels` and `probe` with `--table psi=FILE` read the folder's events table.
  - A highlight row for a cohort not in the data is an input error, not a crash; a cohort that `--keep` or
    `--where` left out is skipped with a warning.
  - A bad `--settings` file is an input error, and `validate` checks the settings too.
  - `--categorical COL` (and `--detail-categorical`) adds COL as a covariate; on its own it was ignored. A column
    named both as a covariate and a stratum is an input error, not a traceback.
  - Kaplan–Meier panels draw with the oldest versions allowed (matplotlib 3.7.0, numpy 1.24, pandas 2.1), where
    they failed with "ufunc 'isfinite' not supported".
  - `event_panel` and `cox_model_figure` refuse a `model` other than the one their `results` were computed with
    (the page would name one model and draw another's numbers), and `event_panel` names the cohorts its `results`
    do not cover instead of failing with a KeyError.
  - `min_group` or `min_pairs` of 0 no longer crash on an event without values.
  - Kaplan–Meier curve labels are kept apart when both curves end at the same survival (e.g. both at 0); they were
    drawn on each other.
  - The probe report's host-expression line names every cohort with p < 0.05 (it stopped at six), ‡ marking low
    power; its low-power note names fits with p < 0.05 in the adjusted model too, as the ranking uses them.
  - The probe no longer raises pandas' FutureWarning on downcasting in `fillna` (real data, pandas ≥ 2.2).
- **Labels and records.**
  - The probe report and the `cox` figure and printout name each model as fitted (no host expression with
    `--no-expression`; a covariate left out of a cohort is not named; the HIT index is called so).
  - The probe report calls settings "relaxed gates" only when a minimum a test needs is lowered (now including
    `low_psi_variance_sd`); a ridge, another HR unit or a stricter minimum is listed without that caution.
  - `analyze` writes `analysis.json` (model and settings), and `Results.read` restores them.
  - Provenance: `call.model` is the forest's model, the clinical columns of both models are hashed, and so are the
    PSI of the gene's other events (`psi_q_family`); the `cox` figure hashes the expression it adjusts for, and its
    q family when it is given gene-wide results.
  - The case-vs-reference view's CSV names the patient behind each point (`patient_id`, `n_samples`), and `validate`
    notes patients with several samples of one group.
  - The `cox` figure's CSV has a `header` row with the model, patients, events and low power it prints.
  - The forest's CSV holds `ph_p` and `ph_marked`; `q_marked` and `ph_marked` follow a drawn CI.
  - Protein changes: a region holding the stop codon (the three nucleotides after the last coding one) is not
    called 3′ UTR; identical proteins from a coding event are not said to lie outside the coding sequence.
  - The gene-track warning gives the real reason when one transcript is not enough for the collapsed model.
- **Packaging.** `lifelines>=0.29.0` (0.27.8 and 0.28 fail to import with SciPy ≥ 1.14) and `pandas>=2.1`. CI runs
  the tests on macOS (Python 3.12) and checks that the minimum versions (`.github/minimum-versions.txt`, Python 3.10)
  install and draw the example; the publish workflow checks that the built wheel installs and draws the example
  before uploading.
- **Docs.** The panel and `cox` examples name `--out` and `--endpoint`; page suffixes and the stacked layout's rule
  are described as they are. The agent guide asks to compare `adj_cox_n` with `cox_n` before reading a change in p
  after adjustment.

## 0.1.1 (2026-10-02)

- **Install guide.** The README's install section covers the PyPI release: a virtual environment, checking the
  install, the Parquet extra (`pip install "splice-assay[parquet]"`), upgrading, and installing from a clone to run
  the tests. No code changes.

## 0.1.0 (2026-10-02)

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
  - Survival per cohort and endpoint: KM and log-rank, split at the median, the mean or a set value (`km_split`,
    `--km-split`; `km_split_expression` for expression). The KM header names the split ("split at median PSI 0.7705").
  - Cox on PSI adjusted for host expression and any clinical covariates (numeric, categorical, strata), with gates,
    failure rules and every model term reported. The PSI HR is shown per IQR, or per SD with `--hr-unit sd`
    (`psi_hr_unit`); the tables hold both.
- **Event panel.**
  - Under the title, one line per event: what it is, in 1-based coordinates, and what its value measures.
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
  - In the figures a dagger (†) marks each such test: after the KM log-rank p, the forest CI and the p of a model
    row, with a legend entry.
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
- **HITindex import.** `import-hitindex` reads HITindex matrices of alternative first exons (AFE), last exons (ALE)
  and the HIT index, with a GTF for strand and symbol.
  - AFE and ALE are PSI.
  - The HIT index (−1 to 1) has its own rules: no 0/1 check, a hit needs |Δ| > `hit_min_abs_delta` (0.20), and it
    is labelled "HIT index" on the pages and in the reports.
  - HIT-index events are left out of `analyse` and `probe` unless `--include-hit` (`include_hit=True`), because
    there is one per exon; a HIT event named with `--event` is always analysed. They form q families of their own.
  - `--append` (both importers) adds events to tables already in the folder.
  - The synthetic example has a third gene, SYN3, with AFE, ALE and HIT events, and writes its HITindex matrices.
- **Subsets.** `--keep FILE` runs any command on the listed patients only. TCGA barcodes are matched to their
  patients, and normals are kept.
  - Results tables keep their full column set when nothing was tested.
  - The probe ranks by KM when no Cox model can be fitted.
- **Condition names from the samples table.** A `role` column (case / reference) names the comparison without
  flags. Group names keep inner capitals in sentences, and the probe report uses them too.
- **q values within each gene.** Benjamini–Hochberg per kind of test (within patients, all samples, KM, Cox PSI
  term; per endpoint and model) over all of the gene's events × cohorts.
  - Shown beside p on the pages, with a footnote naming the family.
  - A `*` marks each q below `q_mark_below` (0.05): after the q, the forest's CI, the tested term's p and in the
    probe overview. The filled markers stay on p.
  - No q for families of fewer than 10 tests (`fdr_min_family`).
  - `panel` and `panels` now analyse every event of the gene so that the family is complete.
- **Two tables are the default input.**
  - `samples` carries survival (wide `OS.time` + `OS`, …), clinical columns and an optional `pair_id`.
  - `psi` carries the event columns beside one PSI column per sample.
  - `example`, `import-rmats` and `import-hitindex` write this form; `--separate` writes separate tables.
  - Separate tables still work and take precedence.
- **Host-gene expression page.** With an expression table, the host gene gets one page of its own, drawn once after
  its splicing pages so they do not repeat it (`expression_panel`).
  - A forest of every cohort with a test (Δ expression, Cox HR per SD), then one row per cohort of the gene's pages:
    tumour vs normal, KM, and Cox on expression + the clinical terms.
  - `analyse_expression` computes the same statistics; the probe writes `expression_cells.csv` and
    `pages/GENE_expression.png` (last in `probe.pdf`).
  - `--no-gex` turns it off.
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
  `import-rmats`, `import-hitindex` and `gtf-subset`.
- **Synthetic example.** A simulated dataset with a small GTF and clinical covariates.
