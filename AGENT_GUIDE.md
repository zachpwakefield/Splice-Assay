# splice-assay: operating guide for agents

This guide is for an agent (or a person) who runs splice-assay for a user: from data to a probe, to assay pages, to
an honest summary. Developer notes are at the end.

## What the tool answers

For an alternative-splicing event it answers three questions per cohort (e.g. per cancer type). The event's value is
PSI per sample for rMATS events (SE, RI, A3SS, A5SS, MXE) and HITindex's alternative first and last exons (AFE, ALE),
or the HIT index of an exon (HIT, −1 to 1).

1. **Does the event differ between two groups?**
   - Within patients: matched pairs, exact signed-rank test.
   - Across all samples: Mann–Whitney test.
2. **Does it go with survival?**
   - Kaplan–Meier / log-rank test, split at the median (or the mean, or a set value).
   - Cox HR per IQR of PSI (or per SD).
3. **Does the survival association hold after adjustment?**
   - Cox with host-gene expression.
   - By default also with age, sex and stage.

Every p value is nominal. The tool ranks and displays evidence; it does not declare discoveries.

## Step 1: get the data into shape

Put two tables in one folder (`.csv`, `.tsv` or `.parquet`), plus an optional expression table:

| Table | Columns |
|---|---|
| samples | `sample_id`, `patient_id`, `cohort`, `group`; survival columns `OS.time` + `OS`, …; clinical columns (age, sex and stage are found automatically); optional `pair_id`, `role` |
| psi | the event columns (`event_id`, `gene`; to draw: `chrom`, `strand`, `event_type`, `constant`, `variable`), then one PSI column (0–1) per sample |
| expression | `gene`, `sample_id`, `value` (e.g. log2 TPM + 1), or wide |

- **This is the default form.** `example`, `import-rmats` and `import-hitindex` write it.
  - Separate `events`, `survival`, `clinical` and `pairs` tables (and a long `psi`) also work, and take precedence.
  - `validate` reports where each table came from ("survival: from the samples table").
- **Subsets:** `--keep FILE` keeps only the patients listed (IDs, or barcodes starting with them), with all their
  samples. `--where COLUMN=VALUE[,VALUE]` keeps the patients with one of the values in a clinical or samples
  column; it ignores case, and repeated conditions must all hold.
  - Check the counts that `validate` prints. A subset of about 80 patients rarely has the 20 deaths Cox needs; then
    only KM runs, and the report says so.
  - Do not lower the gates silently: say so when reporting.
- **Different column names:** do not rename them; map them with `--column logical=actual` (e.g.
  `--column cohort=cancer --column sample_id=File.ID`).
- **Different group names:** any two conditions work, and their names flow into every figure and table.
  - Say which is the case either with `--case "Responder" --reference "Non-responder"` or with a `role` column in
    samples (`case` / `reference`).
  - Effects are case minus reference.
- **rMATS output:** `splice-assay import-rmats` writes one `psi.csv` holding the event columns and the PSI values
  (`--separate`: `events.csv` and a long `psi.csv`).
- **HITindex output** (alternative first and last exons, HIT index): `splice-assay import-hitindex afe.csv ale.csv
  hit.csv --gtf GTF --out DATA --append` adds them beside the rMATS events. The GTF supplies the strand, which the
  HITindex IDs lack.
  - AFE and ALE values are PSI.
  - The HIT index (−1 to 1) is not PSI: it has no 0/1 check, and a hit needs |Δ| > 0.20. Say "HIT index", not PSI,
    when reporting it.
  - `analyse` and `probe` leave HIT events out unless `--include-hit`: there is one per exon, a far larger set. A
    HIT event named with `--event` is always analysed. HIT-index q values come from families of their own.

Then validate, and read the whole output:

```bash
splice-assay validate DATA [--column ...] [--case ...] [--reference ...]
```

What to check:
- **Groups.** The case/reference line names the right groups. "Not tested" lists any other groups.
- **Counts.** The per-cohort counts are plausible: reference samples and pairs per cohort, and patients and events
  per endpoint.
- **Dropped survival rows.** They are reported; a large number means time or event columns were misread.
- **Missing-value codes.** Look at the clinical columns for codes such as `missing` or `[Not Available]`. Pass each
  one as `--na-value CODE`. Age, sex and stage found automatically are cleaned anyway.

Errors are `error: <table>: <what to change>`. Fix the input and rerun; do not work around the checks.

## Step 2: probe (when you do not yet know what to look at)

```bash
splice-assay probe DATA --gene GENE [--gtf annotation.gtf.gz] [--na-value missing]
```

- **What it covers.** Every event of the gene, in every cohort; HIT-index events only with `--include-hit` (the
  report says how many were left out). The endpoint is OS unless `--endpoint` is given. Repeat `--endpoint` (or use
  `--endpoint all`) for several; each gets its own folder.
- **Two Cox models:**
  - base: PSI + host expression;
  - adjusted: + age + sex + stage, found in the clinical table, against stage I and female.
- **Missing variables.** A variable recorded for fewer than 80% of a cohort's patients is left out of that cohort's
  model, and the page says so.
- **Rare levels.** A category with fewer than 10 patients in a cohort, or no deaths, joins its neighbour (stage I
  into "I–II"); the page lists the merge.
- **Output.** `probe_GENE_OS/`:

| File | Read it for |
|---|---|
| `report.md` | Start here: what was run, how many p < 0.05 chance would give, and the ranked events with links to pages ([how to read a page](docs/reading-results.md)) |
| `events.csv` | One row per event (see below) |
| `cells.csv` | One row per event × cohort: every statistic, base and adjusted (`adj_*`), and BH q within each gene |
| `overview.png` | Events × cohorts at a glance: HR colour, p < 0.05 dot (a `*` instead when q < 0.05), group-hit frame |
| `pages/NNN_EVENT.png` | One assay page per ranked event |
| `pages/GENE_expression.png` | The host gene's expression page (with an expression table) |
| `probe.pdf` | The overview, the event pages in rank order, then the expression page |
| `expression_cells.csv`, `proteins.csv` | The expression statistics per cohort; the suggested protein changes (with a protein cache) |

- **Ranking** (events.csv, in this order):
  - `measurable`: events with a Cox fit or a log-rank test in at least one cohort come first;
  - `adj_cox_p05`: cohorts where the adjusted Cox p < 0.05;
  - `cox_p05`: the same for the base model;
  - `group_and_survival`: cohorts with both a group hit and a survival hit;
  - `km_p05`: cohorts where the KM p < 0.05;
  - `best_p`: the smallest p of `best_model` (adjusted Cox, else base Cox, else KM).
- **Other columns:**
  - `best_cohort`, `best_hr_per_iqr` (`best_hr_per_sd` with `--hr-unit sd`) and `min_cox_q`;
  - `expected_by_chance` = 0.05 × cohorts tested;
  - `cox_p05_hr_up` / `cox_p05_hr_down`: directions of the significant cohorts;
  - `share_hr_up`: the share of tested cohorts with HR > 1;
  - `measurable`: false when the event was not testable anywhere.

### Notes on the survival tests

The model header on a page can carry notes in brackets, and a KM header a last line. The probe report sums them up
under "Notes on the survival tests". Each is a reason for caution, not a failed test:
- **`n events per term (< 10)`:** an overfit risk. Read the HR's CI, and prefer the base model when the
  adjusted one is flagged. Typical of small cohorts and subsets.
- **`narrow PSI range (IQR … < 0.05)`:** the HR per IQR covers a few PSI points (with `--hr-unit sd` the note
  checks the SD). It is common for low-inclusion events (retained introns, PSI near 0). Say so when reporting.
- **`stage I merged into II (…)`:** a level too rare to estimate joined its neighbour, so the reference may read
  "I–II".
- **`unstable: …`:** a covariate term with se > 3.
- **`non-proportional hazards: PSI (p …)`** (model header) **or `non-proportional hazards (p …)`** (KM header):
  the Schoenfeld test says the hazard ratio changes over follow-up, so the HR or log-rank result is an average.
  Look at the KM curves; crossing curves are the usual cause. `cells.csv` holds every test's p (`ph_p`,
  `adj_ph_p`, `km_ph_p`).
  - In the figures a † follows each such test: the KM log-rank p, the forest CI, and the p of a model row.

### Host-gene expression

With an expression table, the host gene's expression has a page of its own, drawn once per gene after the splicing
pages (`pages/GENE_expression.png`, last in `probe.pdf`), and the probe writes `expression_cells.csv`.
- At the top, a forest of every cohort: Δ expression and the Cox HR per SD.
- Below it, one row per cohort of the gene's splicing pages: tumour vs normal, KM (median split by default;
  `--km-split-expression`), and Cox on expression + age + sex + stage, with the HR per SD.

Use them when judging an event:
- **Probably an expression echo:** a splicing association in a cohort where expression itself is equally
  prognostic. The splicing model adjusts for expression, but check that the PSI term keeps its p there.
- **Easier to interpret:** a splicing association where expression is not prognostic.
- **Normals:** they need expression values for the tumour–normal expression view. Without them the view says "no
  normal values".

### Optional: protein changes

With a protein cache (`--proteins DIR` or `$SPLICE_ASSAY_PROTEINS`; build one once as in
[docs/annotation-cache.md](docs/annotation-cache.md)), the probe also suggests how each event could change the
protein:
- `protein_change` in events.csv, e.g. `+61 aa, in frame`, `frameshift`, `exon swap, in frame` or `5′ UTR`;
- `proteins.csv`, with the transcripts, how they matched, residues and features;
- a protein band on each page.

`splice-assay proteins DATA --gene GENE` prints the suggestions without a probe. An empty `protein_change` means no
annotated transcript matched well; nothing is shown then.

## Step 3: judge the probe (do this before you say anything to the user)

**An event is worth a closer look when several of these hold.**
- `cox_p05` is clearly above `expected_by_chance`.
- The significant cohorts agree in direction (all ↑ or all ↓).
- The association survives adjustment (`adj_cox_p05` close to `cox_p05`).
- The same cohorts show a group hit (`group_and_survival` > 0), ideally within patients (`within_patient`).
- The per-cohort HRs lean one way overall (`share_hr_up` far from 0.5).

**Treat an event with caution when any of these holds.**
- A single cohort at p ≈ 0.04 among 20 tested.
- Directions that flip between cohorts.
- An association that disappears after adjustment.
- The adjusted model dropped stage in that cohort (see the page note).
- Few events (a KM or Cox gate note on the page).

**Other checks:**
- **Where signal lives.** Open `overview.png`: a whole column of dots means a cohort-wide effect (many events
  associated in one cancer), which is less specific to the event.
- **Duplicate events.** rMATS often lists the same exon several times with different flanking exons, and their PSI
  values are correlated. Count them as one piece of evidence, not several.
- **Coverage.** Low coverage (`frac_obs` in cells.csv) makes an event noisy even when it passes the gates.

## Step 4: make the final pages

```bash
splice-assay panel DATA --event EVENT_ID [--cohort A --cohort B ...] [--gtf ...] [--na-value missing]
```

- **Defaults:**
  - the most promising cohorts (smallest adjusted p; `--top N`, default 3, or `--top all`);
  - endpoint OS;
  - model rows with age + sex + stage found in the clinical table (`--no-detail` for none).
- **Layout.** With three or more cohorts the page stacks: one row per cohort (comparison, KM, model), and the forest
  of every cohort below.
- **What the page says about the event.** Under the title, one line per event gives what it is (type, exon or
  intron, 1-based coordinates, length, strand) and what its value measures. Check it against the event you meant.
- **KM split.** The median by default; the KM header says where it fell ("split at median PSI 0.7705"). Use
  `--km-split mean` or `--km-split 0.3` (a set value) only for a stated reason, and say so when reporting.
- **Expression page.** `panel` also writes `<GENE>_expression_<cohorts>_<endpoint>` for the same cohorts
  (`--no-gex` leaves it out).
- **Many cohorts.** More than 6 are split over pages (`_p1`, `_p2`, …), each with the full forest; set
  `cohorts_per_page` to change this.
- **Explicit models.** `--detail-covariate COL` (repeatable), `--detail-strata COL`, `--baseline COL=LEVEL`.
- **One model on its own:** `splice-assay cox DATA --event E --cohort C --covariate age ...` prints the full
  coefficient table.
- **Outputs.** Each page writes `.png`, `.pdf`, `.svg`, a `.csv` of every plotted value, and a `.provenance.json`
  (input hashes, settings, model, versions). Quote numbers from the CSV, not from the image.

## Step 5: report to the user

- **Name the model.** Say which model a number comes from (base or adjusted) and the endpoint, e.g. "HR per IQR
  1.44 (95% CI 1.11–1.87), p 0.005, Cox with host expression, age, sex and stage, LUAD, OS".
- **What HR per IQR means.** It is the hazard ratio between a patient at the 75th and one at the 25th percentile of
  PSI in that cohort, other terms fixed. HR > 1 means higher PSI goes with worse survival. With `--hr-unit sd` the
  pages show the HR per SD of PSI instead; say which unit a number uses.
- **Markers.** Filled markers mean p < 0.05; a `*` means q < 0.05 within the gene's family; a † means the
  proportional-hazards test failed (p < 0.05). Say which one a claim rests on.
- **Chance.** State how many tests were run and how many p < 0.05 chance predicts. Quote q only with its family:
  e.g. "q 0.03, BH within EHMT2 over the 105 adjusted-model Cox tests for OS".
  - Each kind of test is its own family within the gene: `paired_q`, `unpaired_q`, `km_q`, `cox_q`, `adj_cox_q`.
    HIT-index events form families of their own.
  - `*_q_tests` gives the family size.
  - Families under 10 tests have no q.
- **Wording.** Call results "associations" or "candidates", never discoveries or causes. TCGA-style cohorts are
  observational, and PSI can reflect tumour purity or composition.
- **Composition.** Prefer within-patient evidence for tumour–normal claims. An all-samples shift without a paired
  shift can come from which tumours have matched normals (`composition_sensitive` in cells.csv).
- **Protein changes.** Say they are suggested from annotation: "the annotated isoforms suggest the exon adds 61 aa
  in frame, in a disordered linker". When the PSI form or the other form came from the fallback (`partial_event`,
  `splice_site_only` in proteins.csv) or only one form is annotated, say so. Never present a protein change as
  measured.

## Defaults at a glance

| Setting | Default | Change with |
|---|---|---|
| Cohorts | all (pages pick the most promising 3) | `--cohort`, `--top N` or `--top all` |
| Endpoint | OS | `--endpoint` (repeatable in `probe`, or `all`) |
| Subset | every patient | `--keep FILE`, `--where COLUMN=VALUE` |
| Cohorts per page | 6 (more are split over pages) | settings `cohorts_per_page` |
| Groups | tumour (case) vs normal (reference) | `--case`, `--reference` |
| Base model | PSI + host expression | `--no-expression` |
| Adjusted model | + age + sex + stage (found, cleaned) | `--covariate`/`--strata`, `--no-adjust` |
| Baselines | stage I, female | `--baseline COL=LEVEL` |
| Missing covariate | left out of a cohort's model below 80% recorded | settings `covariate_min_complete` |
| Rare category | merged into its neighbour below 10 patients or with no deaths | settings `level_min_patients` |
| Overfit-risk note | below 10 events per model term | settings `cox_events_per_term` |
| PSI hazard ratio | per IQR | `--hr-unit sd` (per SD; the tables hold both) |
| Narrow-PSI note | PSI spread below 0.05 in the fit cohort, in the HR's unit | settings `narrow_psi_below`, `narrow_psi_measure` |
| Tests | 10 pairs; 10 per group; KM 10 per arm and 10 events; Cox 30 patients and 20 events | `--settings file.json` |
| GTF | none (no gene track) | `--gtf` or `$SPLICE_ASSAY_GTF` |
| Protein cache | none (no protein band) | `--proteins` or `$SPLICE_ASSAY_PROTEINS` |
| KM split | median | `--km-split mean` or a value; `--km-split-expression` |
| HIT-index events | left out of `analyse` and `probe` | `--include-hit`, or `--event` for one |
| q mark | `*` for q < 0.05 | settings `q_mark_below` |
| Proportional-hazards mark | † for p < 0.05 | settings `ph_note_below` |
| Input form | two tables (samples, psi) plus expression | separate tables are also read |
| HIT-index hit | \|Δ\| > 0.20, no 0/1 check | settings `hit_min_abs_delta` |
| Expression page | on when expression is given (one per gene, last) | `--no-gex` |

## Pitfalls

- **MXE orientation.** MXE PSI from rMATS is taken to measure the transcript-upstream exon. Check one known event
  (e.g. FGFR2 IIIb/IIIc) in new data.
- **Endpoints.** Do not mix them: every figure and table is for one endpoint.
- **Patient data.** Never copy patient-level tables into a public place. The pages contain per-patient values in
  their CSVs.
- **Slow first run.** A large GTF makes the first run slow. `splice-assay gtf-subset` makes a small one.
- **The HIT index is large.** It covers every exon, so `--include-hit` multiplies the work: on one real gene the
  probe went from about 3 to 11 minutes, and 25 of its 30 pages were HIT events. Ask for it only when wanted.
- **Sparse AFE and ALE.** For a gene with one dominant first or last exon, most AFE and ALE events are observed in a
  handful of samples and are not measurable. That is expected, not a failure.

## Developer notes

- **Layout.** Python ≥ 3.10; code in `src/splice_assay/`; tests in `tests/`. Run the tests with `pip install -e
  ".[test]"` followed by `pytest`.
- **Statistics.** Definitions are in `docs/methods.md`. Change them only together with a test against an independent
  implementation.
- **Figures.** Figure code must stay deterministic (byte-identical outputs are tested). `docs/make_figures.py`
  regenerates the README figure and the reading guide's crops.
- **Commits.** Commits are authored by the maintainer only: do not add AI co-author trailers or attribution.
