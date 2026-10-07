# splice-assay

Case-vs-reference and survival panels for alternative-splicing events, from plain tables.

Give it PSI values, a sample sheet, survival data and the coordinates of your events. For any event it computes:
- tests between two groups (tumour vs normal, or any two groups you name), within patients and across all samples;
- a Kaplan–Meier test per cohort, split at the median (or the mean, or a value you set);
- a Cox model per cohort, adjusted for host-gene expression and any clinical variables you choose.

It then draws one figure per event: the event on its gene with its coordinates, the group comparison, the survival
curves, a forest of every cohort, and optionally the full Cox model of each cohort shown. That is one "assay" figure
per event. With expression data, the host gene's own expression gets a page of its own, drawn once per gene after
its splicing pages, with the same views (tumour vs normal, KM, Cox on expression + clinical).

When you are not sure where to look, `probe` runs every event of a gene in every cohort and ranks them. It draws an
overview of events × cohorts and a map of each gene with its events, can add how events move together and with their
gene's expression (`--correlation`), and writes a prompt from which an agent can draft a short narrative of the
results (`summarize`).

Events can be rMATS types (SE, RI, A3SS, A5SS, MXE) or HITindex's alternative first and last exons and HIT index
(AFE, ALE, HIT).

![Example panel](https://raw.githubusercontent.com/zachpwakefield/Splice-Assay/main/docs/example_panel.png)

*The example uses simulated data: `splice-assay example demo/`.*

## Guides

- [Reading the results, panel by panel](https://github.com/zachpwakefield/Splice-Assay/blob/main/docs/reading-results.md): every part of a page, the probe's report and
  tables, and a reading order.
- [Getting the annotation cache from SpliceImpactR](https://github.com/zachpwakefield/Splice-Assay/blob/main/docs/annotation-cache.md): what the protein suggestions need,
  and how to build it.
- [Methods](https://github.com/zachpwakefield/Splice-Assay/blob/main/docs/methods.md): every statistic and threshold.
- [Protein suggestions](https://github.com/zachpwakefield/Splice-Assay/blob/main/docs/proteins.md): how events are matched to transcripts, and what the band shows.
- [Agent guide](https://github.com/zachpwakefield/Splice-Assay/blob/main/AGENT_GUIDE.md): the steps from data to report, for agents and new users.

## Install

splice-assay is on [PyPI](https://pypi.org/project/splice-assay/) and needs Python 3.10 or newer. A virtual
environment keeps it apart from your other projects:

```bash
python3 -m venv splice-env
source splice-env/bin/activate
pip install splice-assay
```

`splice-assay --version` checks the install, and `splice-assay example demo/` builds simulated data to try it on
(see the quick start below). The dependencies (numpy, pandas, scipy, matplotlib and lifelines) come with it.

To read Parquet tables, add pyarrow (CSV and TSV need nothing extra):

```bash
pip install "splice-assay[parquet]"
```

To update to a new release:

```bash
pip install --upgrade splice-assay
```

Changes on GitHub that are not yet released:

```bash
pip install "splice-assay @ git+https://github.com/zachpwakefield/Splice-Assay"
```

To develop or run the tests, install from a clone:

```bash
git clone https://github.com/zachpwakefield/Splice-Assay
cd Splice-Assay
pip install -e ".[test]"
pytest
```

The command is `splice-assay` and the Python package is `splice_assay`.

## Quick start

```bash
splice-assay example demo/                     # simulated data, results and six figures
splice-assay validate demo/data                # check your tables; per-cohort counts
export SPLICE_ASSAY_GTF=demo/data/annotation.gtf
export SPLICE_ASSAY_PROTEINS=demo/proteins     # optional: the suggested protein change on every page

# not sure where to look? probe every event of a gene in every cohort: ranked, a page per top event
splice-assay probe demo/data --gene SYN1       # -> probe_SYN1_OS/report.md, probe.pdf, pages/
splice-assay probe demo/data --gene SYN3       # HITindex events: alternative first/last exons (AFE, ALE)
splice-assay probe demo/data --gene SYN3 --include-hit   # ... and the HIT index of each exon (left out by default)
splice-assay probe demo/data --gene SYN1 --correlation   # ... and how its events move together and with expression
splice-assay summarize probe_SYN1_OS           # optional: a narrative drafted by Claude Code (see Agent summary)

# one event: the most promising cohorts, OS, and models adjusted for age + sex + stage are the defaults
splice-assay panel demo/data --event SYN1:SE:1 --out figures/

# or choose everything yourself
splice-assay panel demo/data --event SYN1:SE:1 --cohort COH1 --cohort COH2 --endpoint DSS \
    --detail-covariate age --detail-covariate stage --baseline stage=I --out figures/
splice-assay analyze demo/data --out results/  # statistics for every event x cohort x endpoint, as CSV
```

An agent (or a new user) should read [AGENT_GUIDE.md](https://github.com/zachpwakefield/Splice-Assay/blob/main/AGENT_GUIDE.md): the steps from data to report, how to read a
probe, and what not to claim.

In Python:

```python
import splice_assay as sa

ds = sa.Dataset.from_dir("demo/data")          # or sa.Dataset.from_tables(samples=df, psi=df, ...)
res = sa.analyze(ds)                           # res.groups, res.survival, res.cox_terms
clinical = sa.CoxModel().with_clinical(["age", "sex", "stage"], baseline={"stage": "I"})
panel = sa.event_panel(ds, "SYN1:SE:1", ["COH1", "COH2"], "OS", gtf="demo/data/annotation.gtf",
                       detail=clinical, out_dir="figures/")    # panel.figure is a matplotlib Figure
sa.cox_model_figure(ds, "SYN1:SE:1", "COH1", "OS", model=clinical, out_dir="figures/")
```

## Input tables

Two tables are enough, plus an optional expression table. Put them in one folder, named `<table>.csv` (or `.tsv`,
`.txt`, `.parquet`), or pass DataFrames. IDs are read as text.

| Table | Columns |
|---|---|
| `samples` | One row per sample: `sample_id`, `patient_id`, `cohort`, `group`; the survival columns (`OS.time` + `OS`, `DSS.time` + `DSS`, …); any clinical columns (age, sex, stage, …); optional `role` (case/reference), `survival_cohort` (true/false) and `pair_id` |
| `psi` | One row per event: the event columns (`event_id`, `gene`; to draw it also `chrom`, `strand`, `event_type`, `constant`, `variable`; optional `label`, `gene_id`, `expression_gene`, `psi_junctions`, `other_junctions`), then one PSI column per sample |
| `expression` (optional) | long: `gene`, `sample_id`, `value`; or wide: `gene` + one column per sample |

How they are read:
- **Survival and clinical values** come from each patient's survival sample (one tumour per patient). They may be
  repeated on the patient's other rows or left empty there.
- **Pairs.** Samples sharing a `pair_id` within a cohort form a pair. Without `pair_id`, a patient's one tumour and
  one normal sample are paired.
- **This is the default form.** `splice-assay example` writes its data this way, and `import-rmats` and
  `import-hitindex` write one `psi.csv` holding the event columns.

### Separate tables (also accepted)

Each part can instead have its own table. A separate table takes precedence over the same columns in `samples` or
`psi`, and `validate` says where each table came from. `splice-assay example DIR --separate` writes this form.

`--table NAME=PATH` reads one table from a file elsewhere instead of the folder, for example
`--table clinical=/shared/clinical.tsv`; it replaces the folder's table of that name. NAME is one of the tables above
(`samples`, `psi`, `events`, `survival`, `pairs`, `expression`, `clinical`; any other is an error), and a relative
PATH is read from the current directory, not the data folder. Repeat it for more tables. It works in either form,
on `validate`, `analyze`, `panel`, `panels`, `probe` and `cox`.

| Table | Required columns | Optional columns |
|---|---|---|
| `samples` | `sample_id`, `patient_id`, `cohort`, `group` | `role` (case/reference), `survival_cohort` (true/false) |
| `psi` | long: `event_id`, `sample_id`, `psi`; or wide: `event_id` + one column per sample | |
| `events` | `event_id`, `gene`; to draw also `chrom`, `strand`, `event_type`, `constant`, `variable` | `label`, `gene_id`, `expression_gene`, `psi_junctions`, `other_junctions` |
| `survival` | long: `patient_id`, `endpoint`, `time`, `event`; or wide: `OS.time` + `OS`, `DSS.time` + `DSS`, … | |
| `pairs` | `cohort`, `patient_id`, `case_sample`, `reference_sample` | |
| `expression` | long: `gene`, `sample_id`, `value`; or wide: `gene` + one column per sample | |
| `clinical` | `patient_id` | any covariates: age, sex, stage, … |

### Subsets (e.g. one subtype)

`--keep FILE` restricts any command to the patients listed in a table. The first column holds the IDs, or name the
column with `--keep-column`.
- **Which IDs match.** Patient IDs, sample IDs, or barcodes that start with either. TCGA aliquot barcodes such as
  `TCGA-XX-0001-01A-11R-0000-07` select patient `TCGA-XX-0001`. An ID that is itself a patient or sample ID selects
  only that one, so `1-2` does not also select patient `1`.
- **IDs are text.** `00123` stays `00123`. A plain list without a header line also works: its first line is read as
  an ID when it names a patient or sample (with a warning).
- **What is kept.** All samples of a selected patient, so its normal samples stay for the tumour–normal views.
- **Where it shows.** `validate` reports the subset, the probe report states it, and page titles end with
  "subset: <file>".

```bash
splice-assay probe data/ --gene EHMT2 --keep Her2_samples.tsv
```

`--where COLUMN=VALUE` selects by a column of the clinical or samples table instead, for example
`--where Subtype=Her2`.
- **Matching.** Values match case-insensitively, and numbers by value (`--where grade=1` matches a column read as
  1.0); `--where stage=II,III` accepts either value. A missing value never matches.
- **Several conditions.** Repeat `--where`; all conditions must hold (`--where stage=I,II --where stage=II,III`
  keeps stage II), and `--keep` can be added.
- **Tables keyed by sample.** Survival and clinical tables may be keyed by `sample_id`; they may list samples the
  subset leaves out.
- **What is kept.** All samples of a selected patient, as with `--keep`.

Small subsets often have too few deaths. Cox needs 30 patients and 10 events by default, and the log-rank test 10
patients per arm and 10 events; a subset that misses one gate gets only the other test (the probe then ranks by it),
and one with fewer than 10 deaths gets neither. A Cox fit with fewer than 20 events runs but is marked low power (‡).

### Your own column names

Nothing needs renaming. Map the logical names above to your columns:

```bash
splice-assay validate data/ --column sample_id=File.ID --column patient_id=Case.ID --column cohort=cancer \
    --column group=sample_type --case "Primary Tumor" --reference "Solid Tissue Normal"
```

In Python this is `Dataset.from_dir(..., columns={...}, case=..., reference=...)`.

### Missing-value codes

Empty cells, `NA` and `NaN` are read as missing. Clinical files often use other codes, such as "missing" or
"[Not Available]". Declare those codes with `--na-value missing` (repeatable) or `na_values=[...]`. Otherwise a code
becomes a category of its own in the Cox model.

### Groups

- **Any two conditions.** `group` holds each sample's condition, with whatever names you use (`tissue` is accepted
  too): Responder / Non-responder, Metastasis / Primary, IDH-mutant / IDH-wildtype, …
- **Which is the case.** Effects are case minus reference. The two groups are named, in order of precedence:
  - by `--case` and `--reference` (`case=`, `reference=` in Python);
  - else by an optional `role` column in the samples table: `case` or `reference` (`control` is accepted too);
  - else as tumour and normal.
- **Matching.** Matching ignores case, and "tumor" equals "tumour".
- **Labels.** Figures, tables and the probe report use the names as spelled in `group`, for example "Higher in
  responder" and "Δ median PSI (metastasis − primary)". Names with inner capitals (IDH-mutant) keep them.
  `--case-label` and `--reference-label` change the displayed names only.
- **Other groups.** Samples of any other group are set aside, with a warning.

### Survival

- **One tumour per patient.** Survival analysis uses one case sample per patient. If a patient has several, mark one
  with `survival_cohort`.
- **Wide or long.** Survival can be long, or wide with one time/event column pair per endpoint (e.g. `OS.time` and
  `OS`). Wide pairs are detected automatically.
- **Units.** Time is in days by default (`Settings(time_unit=...)`).
- **Dropped rows.** Rows with a missing or non-positive time are dropped and counted.

### Without reference samples

If your data have no reference group, for example a tumour-only cohort with outcomes, figures leave out the group
view and the Δ forest and show only survival.

### Event geometry

Coordinates are 0-based and half-open (BED and rMATS convention), written `start-end` and joined with `;`.

| `event_type` | `constant` | `variable` (the region PSI measures) |
|---|---|---|
| `SE` | the two flanking exons | the cassette exon |
| `RI` | the two exons around the intron | the intron |
| `A3SS`, `A5SS` | the short form of the alternative exon, and the flanking exon | the extension of the long form |
| `MXE` | the two flanking exons | the exon PSI measures, then the other exon |
| `AFE`, `ALE` | the gene's other first (last) exons, drawn grey (may be empty) | the first (last) exon whose use PSI measures |
| `HIT` | empty | the exon whose HIT index is given |

The junctions of both forms follow from these coordinates. For any other event type, give `psi_junctions` and
`other_junctions` (introns as `donor-acceptor`) to draw arcs. Genes much longer than their events (over 3× and over
20 kb) are drawn around the events only, with marks where the gene continues.

### From rMATS

```bash
splice-assay import-rmats rmats_out/ --b1 b1.txt --b2 b2.txt --counting JCEC --out my_data/
```

This writes one `psi.csv`: the event columns (with geometry), then the PSI values from `IncLevel1`/`IncLevel2`, one
column per sample. Add `samples.csv`, with the survival and clinical columns, to the same folder.
- **Counts.** `--counting JCEC` (junction and exon-body reads, the default) or `JC` (junction reads only) picks the
  rMATS files that are read.
- **Sample names.** The columns are named after the BAM files, or after `--names1` and `--names2` (comma-separated,
  in the order of `b1.txt` and `b2.txt`).
- **Event IDs.** `TYPE:ID`, with rMATS's ID within each type. `--prefix` puts text in front (`--prefix run2_` gives
  `run2_SE:17`), and `--types SE,RI` imports only those types.
- **Layout.** `--separate` writes `events.csv` and a long `psi.csv` instead. `--append` adds the events to the tables
  already in `--out`, as for HITindex below.
- **MXE.** PSI is taken to measure the transcript-upstream exon, as in the data this was checked on (FGFR1–3
  IIIb/IIIc). Check one known event in your data; `--mxe-psi-exon first_listed` takes the exon rMATS lists first.

### From HITindex

```bash
splice-assay import-hitindex afe.csv ale.csv hit.csv --gtf gencode.gtf.gz --out my_data/ --append
```

Each matrix has one row per exon and one column per sample. The first column is HITindex's exon ID,
`gene_id;chrom:start-end;afe` (or `ale`, `hit`), in its 1-based coordinates. A matrix assembled from HITindex's
per-sample outputs looks like this; `splice-assay example` writes three in `demo/hitindex/`.
- **The GTF** gives each gene's strand and symbol, which the IDs lack.
- **Events** are named `GENE:AFE:0001`, numbered per gene and type from 5′ to 3′; `--prefix` puts text in front.
  `source_id` keeps the HITindex ID.
- **`--gene`** keeps the listed genes (symbols or Ensembl IDs). The matrices are streamed, so a large file costs
  little memory.
- **`--append`** adds the events to the tables already in the folder (for example after `import-rmats`), in the
  layout it finds there; an event ID already present is an error. Without it, the events go into a new `psi.csv`
  with their event columns (`--separate`: `events.csv` and a long `psi.csv`).
- **Values.** AFE and ALE values are PSI: the exon's use among the gene's first (last) exons. The HIT index (−1 to 1)
  is not PSI and has its own rules: no 0/1 check, and a hit needs |Δ| > 0.20 (`hit_min_abs_delta`). Pages say "HIT
  index" where they would say PSI. See [docs/methods.md](https://github.com/zachpwakefield/Splice-Assay/blob/main/docs/methods.md#event-types-and-their-values).
- **The HIT index is left out by default.** HITindex reports it for every exon of every gene, a far larger set than
  the splicing events, so `analyze` and `probe` skip HIT events unless `--include-hit` (`include_hit=True`). A HIT
  event named with `--event` is always analyzed and drawn. HIT-index events form q families of their own, so
  including them changes no PSI q value.

## Probing: when you are not sure what to look at

`splice-assay probe DATA --gene GENE` runs every event of the gene (or `--event ...`, or every event in the data) in
every cohort for one endpoint (OS unless `--endpoint`). Repeat `--endpoint` (or use `--endpoint all`) to probe
several endpoints. Each gets its own folder (`probe_GENE_OS`, `probe_GENE_DSS`, or `OUT/OS` with `--out`), because
OS and DSS are never mixed in one figure. It fits two Cox models per cell:
- the base model: PSI + host expression;
- the adjusted model: + age + sex + stage (see below).

It writes:

| File | Content |
|---|---|
| `report.md` | What was run, how many p < 0.05 chance would give, and the ranked events with links to their pages |
| `events.csv` | One row per event, ranked: significant cohorts (adjusted, base, with directions), group hits, best cohort |
| `cells.csv` | One row per event × cohort: every statistic, base and adjusted, with BH q values within each gene |
| `overview.png` | Events × cohorts: HR colour, p < 0.05 dot, group-hit frame; with an expression table, the host gene's own expression in a row of its own |
| `gene_map_GENE.png` | Per gene of the best-ranked events (at most 20): its model (with `--gtf`) and its probed events observed in enough samples, 5′ to 3′, each beside its cells of the overview; with `--correlation`, the median ρ between its events |
| `pages/`, `probe.pdf` | One assay page per ranked, measurable event, at most 30 (`--max-pages N`), then the gene's expression page (`GENE_expression.png`, with an expression table); all of them in one PDF, after the overview and the gene maps (and the correlation figure) |
| `correlations.csv`, `correlation.png` | With `--correlation`: Spearman ρ per cohort of each event with its host gene's expression and with the gene's other events (see Correlation below) |
| `agent_prompt.md` | What an agent needs to write a narrative of the probe: instructions, the rules for reading it, the report and the best-ranked events' notable cohorts (see Agent summary below) |

Events measurable in at least one cohort come first. The ranking then orders events by these criteria in turn, with
Cox counted only in fits that are not low power (at least 20 events):
1. cohorts where the adjusted Cox p < 0.05;
2. cohorts where the base Cox p < 0.05;
3. cohorts with both a group hit and a survival hit;
4. cohorts where the KM p < 0.05;
5. the Cox hits in low-power fits (‡), adjusted, then base;
6. the p of the best cohort: adjusted Cox, else base Cox, else KM (when no Cox model can be fitted).

Each page shows that event's most promising cohorts (`--top`, default 3; `--top all` shows every cohort with a
test, one row each) and every cohort in the forest. The best cohort and the page's cohorts are those with Cox
p < 0.05 in a fit that is not low power, then in a low-power fit, then the other fits (those not low power first),
each by adjusted Cox p.

Pages go to the best-ranked measurable events, at most 30; an event among them without coordinates to draw gets none.
`--max-pages N` draws more or fewer, and a large N draws every one. `events.csv` and `cells.csv` hold every event
either way.

More than 6 cohorts are split over balanced pages (`_p1`, `_p2`, …, "page 1 of 3" in the title), each with the
full forest. Change the split with the `cohorts_per_page` setting (0 = one page).

**Gene maps** show where on the gene the signal sits. Each event is drawn as on its page, 5′ to 3′ under the gene's
collapsed model, with pale columns carrying the gene's exons down through the rows, so events that share an exon line
up. Long introns are drawn shortened. Beside each row are its overview cells, and with `--correlation` a triangle gives
the median ρ between each two events. Events observed in under half of every cohort's survival samples are not drawn
(`--min-observed FRAC`; in TCGA data most of a gene's annotated events can be nearly absent); the subtitle counts them.

![Gene map](https://raw.githubusercontent.com/zachpwakefield/Splice-Assay/main/docs/reading/09_gene_map.png)

A probe ranks candidates; it does not test a hypothesis.

## Protein consequences (optional)

Give splice-assay a protein cache and it suggests how each event could change the protein. It chooses the annotated
transcripts of the two forms with SpliceImpactR's matching, then compares their proteins. Each page gets a protein
band (both isoforms, domains, the event's residues), and probes get a Protein column.

```bash
splice-assay protein-cache ~/splice_assay_annotation --out protein_cache/   # once (docs/annotation-cache.md)
export SPLICE_ASSAY_PROTEINS=protein_cache/                          # or --proteins on panel, panels and probe
splice-assay proteins data/ --gene FNBP1      # e.g. "inclusion adds 61 aa in frame (residues 330–390 of FNBP1-202)"
```

`protein-cache` runs R once (any version, no packages; `--rscript PATH` when `Rscript` is not on the
PATH).

The protein changes are suggestions read from annotation; an event without a good match is not shown. The example
writes a small synthetic cache (`demo/proteins/`). [docs/annotation-cache.md](https://github.com/zachpwakefield/Splice-Assay/blob/main/docs/annotation-cache.md) shows how to build the real
one with SpliceImpactR. [docs/proteins.md](https://github.com/zachpwakefield/Splice-Assay/blob/main/docs/proteins.md) covers the matching, the cache format for other
annotation, and checks against SpliceImpactR.

## Defaults (nothing needs to be chosen)

| Setting | Default | Change with |
|---|---|---|
| Cohorts | all; pages show the most promising 3 | `--cohort`, `--top N` or `--top all` |
| Endpoint | OS (`cox` needs one named) | `--endpoint` (repeatable in `probe`, or `all`) |
| Groups | tumour vs normal | `--case`, `--reference` |
| Subset | every patient | `--keep FILE`, `--where COLUMN=VALUE` |
| Adjustment | host expression, plus age, sex and stage found in the clinical table | `--covariate`, `--no-adjust`, `--no-detail`; `--no-expression` drops expression from both models |
| Baselines | stage I, female (a level with < 10 patients or no deaths joins its neighbour) | `--baseline stage=II`; settings `level_min_patients` |
| Model notes | under 20 events (low power, ‡); under 10 events per term; PSI SD under 0.05 | settings `cox_low_power_events`, `cox_events_per_term`, `narrow_psi_below` |
| KM split | median | `--km-split mean` or `--km-split 0.5`; `--km-split-expression` for expression |
| HIT-index events | left out of `analyze` and `probe` (one per exon, a far larger set) | `--include-hit`, or `--event` for one |
| q mark | `*` for q < 0.05 (filled markers stay p < 0.05) | settings `q_mark_below` |
| PSI hazard ratio | per SD of the cohort's PSI | `--hr-unit iqr` (per IQR) |
| Cohorts per page | 6 (more are split over pages) | settings `cohorts_per_page` |
| Pages (probe) | the best-ranked measurable events, at most 30 | `--max-pages N` |
| GTF | `$SPLICE_ASSAY_GTF` | `--gtf` |
| Protein cache | `$SPLICE_ASSAY_PROTEINS` (none: no protein band) | `--proteins`, `--no-proteins` |
| Expression page | one per gene, after its splicing pages, when an expression table is given | `--no-gex` |
| Correlation (probe) | off | `--correlation` |
| Events in the gene maps and correlations (probe) | observed in at least `coverage_frac` (half) of a cohort's survival samples, in some cohort | `--min-observed FRAC` (0 = all) |
| Agent summary (probe) | off (`agent_prompt.md` is always written) | `--agent-summary`, or `splice-assay summarize PROBE_DIR` |

**How age, sex and stage are found.**
- **Names.** Columns named like `age`, `age_at_diagnosis`, `sex`, `gender`, `stage` or `ajcc_pathologic_stage` are
  used.
- **Cleaning.**
  - Stage is collapsed to its Roman numeral ("Stage IIA" becomes II); stage numbers work too (2 or "2B" becomes II).
  - Sex becomes female or male.
  - Other codes, such as "[Not Available]" or "missing", become missing.
- **Variables missing in a cohort.** A variable recorded for fewer than 80% of a cohort's patients is left out of that
  cohort's model rather than shrinking it; the page says so. This is typical of stage in brain tumours.
- **Rare categories.** A level with fewer than 10 patients in a cohort, or with no deaths, cannot be estimated, so it
  joins its neighbour: stage I with 2 patients becomes part of "I–II", which is then the reference. Other categories
  join the most common level; the most common level itself, when it has no deaths, joins the next most common. The
  page lists each merge. Change the limit with `level_min_patients` (0 = never).
- **Notes.** These warn without changing any result:
  - fewer than 20 events (low power; ‡ after the CI in the forest);
  - fewer than 10 events per model term (an overfit risk);
  - a PSI SD (or IQR, with `--hr-unit iqr`) below 0.05 in the fit cohort, so the HR covers a few PSI points;
  - non-proportional hazards: a term's Schoenfeld test, or the KM split's, has p < 0.05, so the effect changes over
    follow-up. Pages mark such tests with † after the p (KM, model rows) or the CI (forest).

  Settings `cox_low_power_events`, `cox_events_per_term`, `narrow_psi_below`, `narrow_psi_measure` (`sd` or `iqr`)
  and `ph_note_below` control them.

## The Cox model

- **Default model.** PSI (per 0.10) plus host-gene expression, z-scored within the cohort, whenever an expression
  table is given; PSI alone otherwise. Like a clinical variable, expression is left out of a cohort's model (PSI
  alone, with a note) when it is recorded for fewer than 80% of the cohort's patients (none when the table has no row
  for the host gene), or is constant.
- **Adding clinical variables.** Name clinical columns as covariates:
  - numeric columns enter per SD of the cohort;
  - other columns enter as categories against their most common level, or against the level you name with
    `--baseline stage=I`;
  - `--categorical` enters a column as categories even when it is numeric (no need to repeat it with `--covariate`);
  - `--strata` gives a separate baseline hazard per level;
  - `--no-expression` drops the expression term. The adjusted model (age, sex and stage found) is built on this
    model, so it has no expression term either, and it keeps the strata.

```bash
splice-assay analyze data/ --out results/ --covariate age --covariate stage --strata sex
```

- **Regularization (opt-in).** `--ridge clinical` adds a ridge (L2) penalty λ/2·β² to the clinical terms;
  `--ridge molecular` to PSI and host expression; `--ridge all` to every term. β is the log HR per SD (per level for
  categories), and `--ridge-penalty` sets λ (default 1, modest: like a normal prior with SD 1 on each log HR, though a
  sparse category level, which the data inform little, is shrunk more). It steadies models with few deaths or sparse
  categories. Penalized terms are shrunk jointly toward HR 1 (a single HR can move away from 1 when terms are
  correlated), and their CIs and p values are approximate; the model's name says so ("…; ridge λ 1 (clinical terms)").
  With `clinical` or `all`, the shrunk covariates adjust PSI only partly: its HR stays closer to the HR without them,
  so an association that full adjustment would weaken can survive it (docs/methods.md has the numbers).
- **Deaths per term (opt-in).** `--min-events-per-term 5` fits a Cox model only with at least 5 events per estimated
  term (below it: `too_few_events_per_term`). With few deaths per term the PSI p of a model with several clinical
  terms is anti-conservative; whether to require more is your choice (also `cox_min_events`).

- **What the forest shows.** The HR per SD of PSI from this model, with the model written under the axis. The HR per
  SD is the hazard ratio for PSI one standard deviation (of that cohort's PSI) higher, with the other variables held
  fixed.
- **Per IQR instead.** `--hr-unit iqr` (setting `psi_hr_unit`) shows the HR per IQR on the pages, in the probe and its
  overview: the hazard ratio between a patient at the 75th and one at the 25th percentile of PSI in that cohort. The
  tables hold both (`hr_per_sd`, `hr_per_iqr`, with their CIs).
- **Low power.** A fit with 10 to 19 events (`cox_min_events` 10, `cox_low_power_events` 20) runs, but it is marked:
  ‡ after its CI in the forest, "low power" in its model's notes and `cox_low_power` in the tables. A p ≥ 0.05
  there says little against an association.
- **Terms table.** `cox_terms.csv` holds every term of every fitted model.
- **The model band.** `panel` adds the full model of each cohort shown to the figure, below the KM panels: the
  page's model plus the age, sex and stage found in the clinical table. Each term gets an HR, a 95% CI and p, and the
  forests of the band share one axis.
  - `--no-detail` leaves it out. `--detail` shows the page's model alone when no age, sex or stage column is found.
  - `--detail-covariate` (with `--detail-categorical` and `--detail-strata`) adds clinical terms to the band only.
  - The cross-cohort forest then keeps the model every cohort supports. Clinical variables are often missing or
    coded differently between cohorts.
- **Any number of cohorts.**
  - With one or two rows at left (one per event and cohort), the models sit in a band below them.
  - From three rows, each row carries its own model at the right, aligned on one axis, and the cross-cohort forest
    moves below at full width.
  - `--layout side` or `--layout stacked` forces either.
- **One model on its own.** `splice-assay cox` prints and draws one model with all its terms.

## Multiple testing

Each splicing p value has a Benjamini–Hochberg q within its gene, with each kind of test as its own family:
- within-patient tests;
- all-samples tests;
- KM;
- the Cox PSI term (per endpoint and model).

The pages print q beside p and name the family in a footnote. Families of fewer than 10 tests get no q. HIT-index
events, when included, form families of their own. A `*` marks each q below 0.05: after the q itself, after the
forest's CI and after the tested term's p in the model rows, and in the probe overview. It is separate from the
filled markers, which show p < 0.05. The `q_mark_below` setting changes the threshold (0 = no marks). Hit rules, filled markers and the probe ranking stay on p. Details:
[docs/methods.md](https://github.com/zachpwakefield/Splice-Assay/blob/main/docs/methods.md#multiple-testing-q-values).

## Host-gene expression page

When an expression table is given, the host gene's own expression gets a page of its own. It is drawn once per gene,
after the gene's splicing pages, so those pages do not repeat it:
- a forest of every cohort with a test: Δ expression (case vs reference) and the Cox HR per SD, with the cohorts
  drawn below shaded;
- for each cohort drawn (those of the gene's splicing pages): case vs reference (pairs and all samples), a KM split
  (median by default; `--km-split-expression`), and the Cox model OS ~ expression + the same clinical terms as the
  model rows (HR per SD).

This answers whether the gene's level itself is prognostic or shifted, beside the splicing results. The splicing
models already adjust PSI for expression.
- **Where it goes.** `panel` writes it as `<GENE>_expression_<cohorts>_<endpoint>`. The probe writes
  `pages/<GENE>_expression.png`, last in `probe.pdf`, and the statistics of every cohort to
  `expression_cells.csv`.
- **`--no-gex`** leaves it out.
- **Matching.** Expression comes from the expression table, matched by the event's `gene` (or `expression_gene`).
  Normals need values too for the case-vs-reference view.
- **In the probe overview.** Below the events, one row per host gene: the HR per SD of expression (Cox on expression
  + the adjusted model's clinical terms), a dot for p < 0.05 and a frame for an expression group hit. No q mark,
  because expression is not adjusted for multiple testing.

## Correlation (optional)

`splice-assay probe ... --correlation` (`probe(..., correlation=True)`) adds Spearman correlations, computed in each
cohort over its survival samples (one tumour per patient, the samples the survival tests draw on):
- **each event with its host gene's expression.** A strong correlation means PSI and expression carry overlapping
  information. The Cox models adjust PSI for expression, so PSI's HR there is what it adds beyond expression,
  estimated less precisely; the KM split is not adjusted, so a KM hit of such an event could be expression's.
- **each pair of events of one gene.** Strongly correlated events (often the same exon listed with different
  flanking exons) are one piece of evidence, not several. Two alternative first (or last) exons of a gene share
  its first (last) exon use and sum to 1, so their ρ is near −1 by construction.

It writes `correlations.csv` (one row per pair × cohort: `rho`, `corr_p`, `corr_q`, `corr_n`, `corr_status`) and
`correlation.png` (each cohort's ρ, gene by gene, after the gene maps in `probe.pdf`), adds the median over the
cohorts of each pair of events as a triangle to each gene map, adds `expr_rho` to `cells.csv` and `best_expr_rho` (ρ in
the best cohort) to `events.csv`, and summarises the strong pairs (|ρ| ≥ 0.7, `corr_note_above`) in the report. A pair
needs 20 patients with both values (`corr_min_n`). q values are Benjamini–Hochberg within each gene and kind. Events
observed in under half of every cohort's survival samples are left out, as they are from the gene maps
(`--min-observed`). In Python, `splice_assay.correlation.correlations(ds, events=...)` returns the table on its own.

## Agent summary (optional)

Every probe writes `agent_prompt.md`: instructions for an agent, the rules for reading a probe and what not to claim
(as in the agent guide), the report, and the notable cohorts of the best-ranked events. It holds aggregate statistics
only, never the pages' per-patient values or your tables.

`splice-assay summarize PROBE_DIR`, or `probe ... --agent-summary`, gives it to Claude Code (`claude -p`, with no
tools, run from an empty folder) and writes:
- `agent_summary.md`: a short narrative (bottom line, candidates, caveats, next checks), marked as machine-written.
  Below it, a number check lists the decimal numbers the narrative quotes that the results do not print, and each HR
  quoted with a CI that the results do not print with that CI, p and q for the event and cohort (and model, where
  named) beside it. Whole numbers are not checked.
- `agent_summary.json`: the agent, the model, the date, a hash of the prompt, and the cost Claude Code reports.

Before you use it:
- **What is sent.** The prompt goes to Anthropic through your Claude Code login. Results from controlled-access data
  (TCGA clinical data, for example) may fall under a data use agreement: check it first. `--dry-run` shows what would
  be sent, and to which command.
- **It is a draft.** Read it against `report.md` and the pages. Its wording differs from run to run; the probe's other
  outputs are unchanged.
- **Setup.** Install Claude Code and sign in once (`claude`, then `/login`), or set `ANTHROPIC_API_KEY`.
  `--agent-model opus` (or `sonnet`) picks the model; by default Claude Code's own is used. Both routes wait up to
  10 minutes for the answer; `summarize --timeout SECONDS` changes that. Your own Claude Code
  settings (`~/.claude/CLAUDE.md`, hooks, output style) still apply and can change the summary's style. If the
  summary fails, the probe's outputs stand; `splice-assay summarize PROBE_DIR` tries again.
- **Cost.** Signed in with a Claude subscription, a summary counts toward the plan's usage limits and is not charged
  separately; the cost in `agent_summary.json` is Claude Code's estimate at API prices. With `ANTHROPIC_API_KEY`
  set, Claude Code may use the key instead, which is billed per token (`claude`, then `/status`, shows which). A
  summary is a small task: `--agent-model sonnet` uses less of a plan than Opus.
- **Other agents.** Any agent can be given `agent_prompt.md` by hand, without access to the probe's folder: the pages'
  CSV files hold per-patient values.

## Outputs

- **Statistics.** `analyze` writes three tables. Every column is defined in [docs/methods.md](https://github.com/zachpwakefield/Splice-Assay/blob/main/docs/methods.md).
  - `group_tests.csv`: one row per event × cohort.
  - `survival.csv`: one row per event × cohort × endpoint.
  - `cox_terms.csv`: one row per model term.
  - `analysis.json` beside them records the Cox model and the settings, so `Results.read` restores both.
- **Protein changes.** With a protein cache, `proteins --out` and the probe's `proteins.csv` hold one row per event
  (see [docs/proteins.md](https://github.com/zachpwakefield/Splice-Assay/blob/main/docs/proteins.md)).
- **Figures.** `panel` and `cox` write four kinds of file:
  - `<stem>.svg`, `.pdf` and `.png` (400 dpi);
  - `<stem>.csv`, with every plotted value;
  - `<stem>.provenance.json`, with input hashes, settings, the models (the forest's and the model rows') and
    versions. The hashes cover every input the drawn numbers depend on, including the PSI of the gene's other
    events, whose tests set the q values.
- **Reproducibility.** The same inputs, settings and library versions give byte-identical files.

### Figure options

- **Highlighting.** `--highlight table.csv` marks forest cells with a short label, such as a tier letter. Its columns
  are `event_id`, `cohort`, `label`, and optionally `endpoint` and `colour`. `--highlight-title Tier` names the
  column. A cohort that is not in the data is an error; with `--keep` or `--where` its rows are skipped, with a
  warning.
- **Batch figures.** `splice-assay panels data/ --spec panels.csv` makes many figures from a table.
  - Columns: `events` and `cohorts` (both separated by `;`), `endpoint`, and optionally `stem`.
  - The statistics are computed and the GTF is read only once.
  - No model band unless asked: `--detail` adds it as `panel` draws it (the page's model plus the age, sex and stage
    found), and `--detail-covariate` adds the terms you name.
- **Two events.** Two events of the same gene can share a figure, for example two retained introns.
- **Smaller GTF.** `splice-assay gtf-subset gencode.gtf.gz --events my_data/psi.csv --out small.gtf.gz` keeps only the
  records of your events' genes and those within 20 kb of the events (`--flank NT`).

## Settings

Every threshold has a default taken from the analysis these figures were designed for. You can change any of them
with a JSON file (`--settings my.json`) or in Python (`sa.Settings(min_pairs=5)`). The main defaults:

- **Groups.** 10 pairs for the paired test; 10 case and 10 reference patients for the unpaired test (a patient's
  several samples in one group, e.g. replicate aliquots, enter as their mean). A hit needs
  p < 0.05 and |Δ median PSI| > 0.10, and must be robust to PSI values of exactly 0 or 1. For the HIT index, a hit
  needs |Δ| > 0.20 and there is no 0/1 check.
- **Survival.**
  - PSI observed in at least 50% of the survival samples, with at least 10 values away from the mode.
  - KM: split at the median (`km_split`: `median`, `mean` or a value), 10 patients per arm and 10 events.
  - Cox: 30 patients and 10 events; a fit with fewer than 20 events is marked low power (‡).
- **Other analysis settings.**
  - `time_unit` (days) and `days_per_year` (365.25): the survival times' unit, converted to years for the KM axes;
  - `psi_step` (0.10): the PSI step of `cox_beta` and `hr_per_step` in `survival.csv` and of the `psi` row in
    `cox_terms.csv` (the HR per SD or IQR does not depend on it);
  - `ci_z` (1.96, for every 95% interval);
  - `ph_test` (the Schoenfeld tests) and `robust_01` (a hit's 0/1 check), both on;
  - `round_decimals` (12): PSI and differences are rounded before comparisons.
- **Drawing.**
  - `formats` (svg, pdf, png) and `dpi` (400, for the PNG files);
  - the KM axes: `km_max_years` (10) and a tick every `km_tick_years` (2);
  - `case_label` and `reference_label`, as `--case-label` and `--reference-label`;
  - the gene model: it leaves out transcripts of the `exclude_transcript_types` (retained_intron), draws genes of the
    `nested_biotypes` (snoRNA, scaRNA) in a track under it, and reads the GTF within `gtf_flank` (20000 nt) of the
    events.
- **Every setting.** `sa.Settings().to_dict()` lists them all with their values, and
  [config.py](https://github.com/zachpwakefield/Splice-Assay/blob/main/src/splice_assay/config.py) describes each in a
  comment. Reports and figures list every analysis setting you changed from its default.

## Scope and limits

- **Figure scope.** One endpoint per figure, and one or two events of one gene, measured on one scale (PSI or HIT
  index).
- **p values.** Nominal. q values (Benjamini–Hochberg within each gene, one family per kind of test) are printed
  beside them; hits, filled markers and the probe ranking use p.
- **Causality.** Cox models describe association; they do not establish causation.
- **Patient-level data.** None is shipped; the example is simulated.
- **Protein changes.** Suggestions from annotated isoforms, not measurements (see [docs/proteins.md](https://github.com/zachpwakefield/Splice-Assay/blob/main/docs/proteins.md)).

## Citation and licence

- **Licence.** GNU General Public License, version 3 only ([LICENSE](https://github.com/zachpwakefield/Splice-Assay/blob/main/LICENSE)), as for SpliceImpactR, whose
  transcript matcher the protein layer ports.
- **Citation.** A citation file will follow; until then, cite this repository.
