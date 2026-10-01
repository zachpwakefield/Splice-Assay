# Methods

Every number splice-assay prints or draws is defined here. The defaults are the settings in `Settings`; each
threshold can be changed there. PSI is on a 0–1 scale, and every effect is case minus reference (for example
tumour minus normal).

## Cohorts, groups and samples

- **Cohort.** A group of patients analysed on its own, such as a cancer type or a study arm. Nothing is pooled across
  cohorts.
- **Groups.** `case` and `reference` name the two sample groups compared; the defaults are tumour and normal. Samples
  of other groups are set aside.
- **Pairs.** A case and a reference sample of the same patient in one cohort. They come from a `pairs` table, or, when
  none is given, from patients with exactly one sample of each group. A patient with more must be paired explicitly.
- **Survival samples.** One case sample per patient (`samples.survival_cohort`; by default every case sample, which
  must then be unique per patient and cohort).
- **Endpoint cohort.** Survival samples whose patient has a valid survival row for the endpoint: time finite and
  greater than 0, event 0 or 1. Every endpoint is analysed and drawn separately.

## Case vs reference (per event and cohort)

| | Paired (within patients) | Unpaired (all samples) |
|---|---|---|
| Data | pairs with both PSI values observed | every case and every reference sample with PSI observed |
| Minimum | `min_pairs` = 10 pairs | `min_group` = 10 of each |
| Test | exact two-sided Wilcoxon signed-rank, conditional on ties, zeros dropped | Mann–Whitney, normal approximation with continuity and tie correction |
| Δ | median of the differences | median(case) − median(reference) |
| HL | one-sample Hodges–Lehmann (median of Walsh averages) | two-sample Hodges–Lehmann (median of pairwise differences) |

- **Exact signed-rank p.** Doubled average ranks are integers. A dynamic program over the probability mass of the
  positive-rank sum enumerates every sign assignment exactly, so ties do not need a normal approximation.
  - Differences are rounded to 12 decimals first, so floating-point noise cannot create or break a tie or a zero.
  - p = min(1, 2·P(W ≤ min(W⁺, W_total − W⁺))).
- **Robustness to PSI 0/1.** The test is repeated without PSI values of exactly 0 or 1. For pairs, a pair is dropped
  when either value is 0 or 1. The result has four states:
  - `not_applicable`: nothing was dropped;
  - `untestable`: the reduced set is below the minimum;
  - `pass`: p < α, the same sign, and |Δ| > `min_abs_delta`;
  - `fail`: otherwise.
- **Hit.** p < `alpha` (0.05), round(|Δ|, 12) > `min_abs_delta` (0.10), sign(HL) = sign(Δ), and a 0/1 state of
  `pass` or `not_applicable`. Otherwise the hit status names the first rule that failed, in this order:
  - `not_tested`
  - `not_significant`
  - `small_effect`
  - `direction_ambiguous`
  - `robustness_fail` or `robustness_untestable`
- **`group_hit`.** The paired or the unpaired design is a hit.
- **Composition diagnostics** (unpaired side). These are descriptive and enter no rule.
  - `matched_*`: the case samples of paired patients against all reference samples.
  - `other_*`: the remaining case samples against the matched ones.
  - `composition_sensitive`: the paired and unpaired Δ differ by more than `min_abs_delta` or in sign.
- **Within-patient support.** The paired test ran, and either the paired design is a hit or the matched comparison
  is a hit in the unpaired direction.
- **No within-patient check.** A hit in a cohort where the paired test was not possible, because there were too few
  pairs or all differences were zero.

## Survival (per event, cohort and endpoint)

- **Coverage gate.** Applied on the survival samples. PSI must be observed in at least `coverage_frac` (50%) of them,
  and at least `min_off_modal` (10) observed values must differ from the modal value. Otherwise neither test runs
  (`coverage_gate`).
- **KM split.** The median PSI over the survival samples, fixed per cohort and shared by every endpoint.
  - PSI ≤ median is the low arm; PSI > median is the high arm.
  - The log-rank test needs `km_min_group` (10) patients per arm and `km_min_events` (10) events.
  - The log-rank HR is (O/E of the high arm) / (O/E of the low arm).
  - Bands are 95% log-log (exponential Greenwood) intervals. The at-risk table counts patients still followed at each
    tick.
- **Cox model.** lifelines `CoxPHFitter`, with Efron ties and no penalty.
  - Model: h(t) = h₀(t)·exp(β·PSI/0.10 + γ·z(host expression) + covariates) (`CoxModel`).
  - Host expression enters when an expression table is given. It is z-scored over the fit cohort.
  - Numeric clinical covariates enter per SD of the fit cohort, or per unit with `scale=False`.
  - Other covariates enter as indicator variables against a baseline level. That is the level named in
    `CoxModel(baseline=...)` when it occurs in the cohort, else the cohort's most common level.
  - Rare levels are merged before fitting. A level with fewer than `level_min_patients` (10) patients in the fit
    cohort, or with no events, has no stable estimate (with no events its likelihood has no maximum).
    - Stage numerals and numbers merge into the adjacent level with fewer patients, the rarest first. Stage I with 2
      patients joins II as "I–II", and IV with 3 joins III as "III–IV".
    - Other categories merge into the most common level ("white+asian").
    - A named baseline refers to its merged level ("I–II"). `cox_notes` lists each merge with its patients and
      events. `level_min_patients` = 0 turns merging off.
  - Strata give each level its own baseline hazard.
  - Patients missing any model variable are left out, and the count is reported (`cox_n_dropped`). A covariate that
    is constant in the cohort is dropped, with a note (`cox_notes`).
  - Requirements: `cox_min_n` (30) patients, `cox_min_events` (20) events, and `min_off_modal` PSI values away from
    the mode.
  - HR per IQR = exp(β·IQR/0.10), where IQR is the interquartile range of PSI in the fit cohort. Its 95% CI is
    exp((β ± 1.96·se)·IQR/0.10). It is left empty when the IQR is 0.
  - The PH test is lifelines' `proportional_hazard_test` with the KM time transform, for unstratified models. It is
    reported, not used.
  - `cox_terms` holds every term: coefficient, se, HR (per its unit), 95% CI and p.
- **Low variance.** When SD(PSI) is below `low_psi_variance_sd` (0.002), neither test runs (`low_psi_variance`).
- **Failed fits.** A fit is marked `failed` and reports no estimate when any of these occur: a lifelines convergence
  warning, a numeric warning inside lifelines, an exception, a non-finite estimate, or se ≤ 0. A Newton–Raphson
  failure is first refitted once with step size `nr_refit_step_size` (0.5), and `cox_fit_note` records the refit.
- **Survival hit.** KM p < α or Cox p < α, among the tests that ran.

p values are nominal. q values within each gene are defined in the next section.

## Multiple testing (q values)

Every p value for splicing gets a Benjamini–Hochberg q within its gene. Each kind of test is its own family:

| Family | Tests in it (one gene) | Column |
|---|---|---|
| Within patients | the paired signed-rank tests of every event × cohort | `paired_q` |
| All samples | the Mann–Whitney tests of every event × cohort | `unpaired_q` |
| KM | the log-rank tests of every event × cohort, per endpoint | `km_q` |
| Cox | the PSI terms of every event × cohort, per endpoint and model | `cox_q` (`adj_cox_q` for the probe's adjusted model) |

- **What counts.** Only tests that ran (status `tested`) are in a family. `*_q_tests` gives the family size. A family
  of fewer than `fdr_min_family` tests (default 10) gets no q.
- **Where q is shown.** The pages print q under or beside the p values: the tumour–normal titles, the KM header,
  and the Cox model headers. A footnote names the family.
- **How the family is built.** `panel`, `panels` and `probe` analyse all of a gene's events so the family is the
  whole gene; `analyse` uses the events it is given.
- **What stays on p.** The hit rules, the forest's filled markers and the probe's ranking use p, as in the screen.
  q is reported beside p.
- **Not adjusted.** Host-gene expression statistics.
- **Why BH within a gene.** The tests of one gene are correlated, mostly positively: rMATS lists the same exon with
  several flanks, and one cohort's events share patients. BH controls the false discovery rate under such positive
  dependence. Treat q as a guide to how much of the gene's signal could be chance, not as a per-test guarantee.

## Host-gene expression (per gene, cohort and endpoint)

With an expression table, the host gene's own expression gets the same three analyses as PSI (`analyse_expression`;
the expression rows of each page; the probe's `expression_cells.csv`).
- **Case vs reference.** The same paired signed-rank and Mann–Whitney tests and composition diagnostics, on the
  expression table's scale.
  - There is no 0/1 robustness check, because that is specific to PSI.
  - A hit needs |Δ median| > `gex_min_abs_delta`: 1.0 by default, which is two-fold on a log2 scale.
- **KM.** A median split of expression over the cohort's survival samples (expression ≤ cut is the low arm), with the
  same coverage and size gates as PSI.
- **Cox.** h(t) = h0(t) exp(b · z(expression) + clinical terms). The clinical terms and strata are those of the
  page's model rows (by default age, sex and stage), with the same completeness rule. The HR is per SD of
  expression in the fit cohort. Status `constant_expression` when expression does not vary.
- **Interpretation.** The splicing models already adjust PSI for expression. The expression rows show whether the
  gene's level itself goes with survival or differs between the groups, beside the splicing result.

## Default clinical adjustment and the probe

- **Found variables.** Age, sex and stage are found by column name and cleaned (`clinical.py`):
  - stage becomes its overall Roman numeral;
  - sex becomes female or male;
  - other codes become missing.
- **Adjusted model.** Base model + age (per SD) + sex (against female) + stage (against I). A stage too rare to
  estimate in a cohort joins its neighbour (see the Cox model above), so stage I with 2 patients makes the
  reference I–II.
- **Variables missing in a cohort.** Below `covariate_min_complete` (80%) recorded in a cohort, a variable is left out
  of that cohort's model (`cox_notes`).
- **Unstable terms.** A covariate term with se > 3 (a sparse level, or near separation) is reported and flagged as
  unstable in `cox_notes`.
- **Flags on a fit.** Two notes go in `cox_notes` and on the page's model header. Both are notes only: the fit runs
  and its estimates are unchanged.
  - **Overfit risk:** fewer than `cox_events_per_term` (10) events per estimated term. Terms are the PSI, expression
    and covariate coefficients, one per level beyond the reference; strata do not count. The value is
    `cox_events_per_term` in the tables.
  - **Narrow PSI range:** the PSI spread in the fit cohort is below `narrow_psi_below` (0.05), so the HR per IQR
    describes a change of a few PSI points. The spread is the IQR (the HR's unit) or the SD
    (`narrow_psi_measure`). The flag is `psi_narrow` in the tables.
  - The probe report counts both flags and names the flagged fits with p < α.
- **Probe ranking.** For one endpoint, events measurable in at least one cohort (a Cox fit or a log-rank test ran)
  come first. Events are then ranked by these criteria in turn:
  1. the number of cohorts with an adjusted Cox p < α;
  2. the same for the base Cox model;
  3. the number of cohorts with both a group hit and a survival hit (KM or base Cox p < α);
  4. the number of cohorts with a KM p < α;
  5. the smallest p: adjusted Cox, else base Cox, else KM (the fallback when no Cox model can be fitted).
- **Probe q values.** The probe analyses every event of each gene it probes, so its q values are the gene-wide
  families of the previous section; `adj_cox_q` is the adjusted model's Cox family.
- **Cohorts shown per page.** They are those with the smallest adjusted Cox p, then base Cox p, then KM p.

## The figures

- **Event panel: schematic.**
  - **Collapsed gene model.** Built from the host gene's GTF transcripts that have at least 2 exons, excluding
    `retained_intron`.
    - Each exon is keyed by its splice sites. A transcript-terminal end is drawn at the median of the transcripts
      that share the key.
    - A key is kept when at least `gene_model_min_frac` (10%, and at least 2) of the transcripts use it. Exons that
      span an intron used that often are dropped.
    - Overlapping exons are merged into one block.
  - **Nested genes.** snoRNA and scaRNA genes inside the host gene are drawn in black.
  - **Event rows.** Each event row shows constant exons in grey and the region PSI measures in the event-type colour.
    For MXE, the other exon is outlined. The arcs above the row are the junctions of the form PSI counts; the arcs
    below are those of the other form.
  - **Orientation.** Minus-strand genes are drawn 5′→3′.
  - **Long genes.** A gene more than 3× and more than 20 kb longer than its events is drawn in a window around the
    events, with marks where it continues.
- **Event panel: group view.**
  - With at least `min_pairs` pairs: one line per pair (left), beside all samples (right, with the KM split).
  - Otherwise: all samples, with the reason there is no paired test.
  - With no reference samples in the cohort: the reason only.
  - With no reference samples anywhere in the data: no group view.
- **Event panel: KM.** The split, bands, censoring ticks and the at-risk table. The log-rank HR, p and events are
  printed above the axes.
- **Event panel: forest.** One row per cohort where a plotted event had a tested comparison or a Cox fit, plus the
  cohorts shown at left.
  - Left: the paired Δ (filled) and unpaired Δ (open).
  - Right: the HR per IQR with its 95% CI, filled when Cox p < α.
  - The model is written under the axis.
  - The HR axis snaps to 1/8…8. It shows every shaded cell's CI in full; other CIs beyond the axis end in an
    arrowhead.
- **Event panel: model band** (`detail`). The full Cox model of each cell shown, on one shared HR axis and one
  column layout, so all models align.
  - **Side layout** (default for one or two cells): the models sit two per line below the cells.
  - **Stacked layout** (default from three): each cell's model sits at the right of its KM, and the forest moves
    below at full width.
  - Its model can add clinical terms to the panel's own model; the forest above keeps the panel's model.
  - A cell whose model could not be fitted shows the reason.
- **Cox model figure.** Every term of one cell's model, with HR, 95% CI and p. PSI is shown per IQR; host expression
  and numeric covariates per SD; categories against their baseline level.

## Reproducibility

- **Identical outputs.** The same inputs, settings and library versions give byte-identical SVG, PDF, PNG, CSV and
  JSON (tested).
- **Plotted values.** Every figure writes a CSV of each value it draws.
- **Provenance.** Each figure also writes a `.provenance.json` with:
  - the SHA-256 of the exact input rows used;
  - the settings;
  - the model;
  - the call;
  - the software versions.

## Validation

- **Tests.** The test suite checks each statistic against an independent implementation:
  - brute-force enumeration for the signed-rank p;
  - scipy for the Mann–Whitney p;
  - lifelines for log-rank, KM and Cox, including models with covariates and strata.
- **Acceptance run.** An acceptance run outside this repository (it uses controlled-access data) fed the package plain
  tables exported from the locked analysis these figures were designed for.
  - It reproduced that analysis's tables: group tests, coverage, KM, and Cox with and without clinical covariates.
  - It also reproduced the plotted values of its figures.
  - The largest difference was 4.4e-16.

## Protein consequences

The optional protein suggestions (transcript matching, effect, features) are described in
[proteins.md](proteins.md).
