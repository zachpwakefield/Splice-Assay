# Methods

Every number splice-assay prints or draws is defined here. The defaults are the settings in `Settings`; each
threshold can be changed there. PSI is on a 0–1 scale, the HIT index on a −1 to 1 scale (see Event types), and every
effect is case minus reference (for example tumour minus normal).

## Cohorts, groups and samples

- **Cohort.** A group of patients analyzed on its own, such as a cancer type or a study arm. Nothing is pooled across
  cohorts.
- **Groups.** `case` and `reference` name the two sample groups compared; the defaults are tumour and normal. Samples
  of other groups are set aside.
- **Pairs.** A case and a reference sample of the same patient in one cohort. They come from a `pairs` table, or, when
  none is given, from patients with exactly one sample of each group. A patient with more must be paired explicitly.
- **Survival samples.** One case sample per patient (`samples.survival_cohort`; by default every case sample, which
  must then be unique per patient and cohort).
- **Endpoint cohort.** Survival samples whose patient has a valid survival row for the endpoint: time finite and
  greater than 0, event 0 or 1. Every endpoint is analyzed and drawn separately.

## Event types and their values

| Type | Source | Value |
|---|---|---|
| SE, RI, A3SS, A5SS, MXE | rMATS (`import-rmats`) or your own table | PSI: inclusion of the event's variable region (SE exon, retained intron, long form, first MXE exon) |
| AFE, ALE | HITindex (`import-hitindex`) | PSI: the use of this first (last) exon among the gene's first (last) exons |
| HIT | HITindex (`import-hitindex`) | HIT index of the exon: −1 when it is used only as a first exon, 1 only as a last exon, about 0 when internal |

- **AFE and ALE** are PSI and follow every PSI rule. Their schematic draws the gene's other first (last) exons in
  grey.
- **The HIT index is not PSI** and has its own rules (`analysis.event_settings`):
  - values must lie in [−1, 1] (PSI must lie in [0, 1]);
  - no robustness check at 0 or 1, because those are not boundary values here (its state is `not_applicable`);
  - a hit needs |Δ| > `hit_min_abs_delta` (0.20) instead of `min_abs_delta`, because the index spans two units;
  - the tests themselves are the same: signed-rank, Mann–Whitney, KM and Cox, with the HR per SD of the HIT index;
  - the narrow-range note uses the same `narrow_psi_below` on the HIT index's spread.
- **Labels.** Pages, tables and reports say "HIT index" wherever they would say PSI (the Cox term is `HIT index`). A
  page draws events of one quantity: a HIT event cannot share a page with a PSI event.
- **Left out by default.** HITindex reports the HIT index for every exon of every gene, a far larger set than the
  splicing events. `analyze` and `probe` therefore skip HIT events unless `include_hit` (`--include-hit`); a HIT
  event named explicitly is always analyzed.
- **q families.** AFE and ALE events share their gene's PSI families. HIT-index events form families of their own
  (the same kinds of test over the gene's HIT-index events), so including them changes no PSI q value.
- **Event details.** Under each page's title, one line per event says what it is, in 1-based inclusive coordinates
  (e.g. "cassette exon chr7:103,001–103,150 (150 nt) between exons …; PSI = its inclusion · plus strand").

## Case vs reference (per event and cohort)

| | Paired (within patients) | Unpaired (all samples) |
|---|---|---|
| Data | pairs with both PSI values observed | every case and every reference patient with PSI observed, one value per patient and group: a patient with several samples in a group (e.g. replicate aliquots) enters with their mean, so that each patient counts once per group |
| Minimum | `min_pairs` = 10 pairs | `min_group` = 10 patients of each |
| Test | exact two-sided Wilcoxon signed-rank, conditional on ties, zeros dropped | Mann–Whitney, normal approximation with continuity and tie correction |
| Δ | median of the differences | median(case) − median(reference) |
| HL | one-sample Hodges–Lehmann (median of Walsh averages) | two-sample Hodges–Lehmann (median of pairwise differences) |

- **Counts.** `unpaired_n_case` and `unpaired_n_reference` count patients (the values tested);
  `unpaired_n_case_samples` and `unpaired_n_reference_samples` count the samples behind them. With one sample per
  patient and group the two agree. `validate` notes patients with several samples of one group.
- **Replicates.** The mean is meant for technical replicates (aliquots of one tissue). Samples of different kinds,
  such as a primary tumour and a metastasis, belong in separate groups. The paired test keeps the pairs' own samples;
  the unpaired test and the composition diagnostics below use the per-patient values. A patient with a sample in
  each group still contributes to both; the Mann–Whitney test ignores that pairing, which makes it conservative when
  a patient's values are positively correlated (the usual case).
- **Rounding.** Values are rounded to 12 decimals before they are ranked or compared, so floating-point noise
  cannot create or break a tie, and two groups equal up to that noise are `constant` (not tested).

- **Exact signed-rank p.** Doubled average ranks are integers. A dynamic program over the probability mass of the
  positive-rank sum enumerates every sign assignment exactly, so ties do not need a normal approximation.
  - Differences are rounded to 12 decimals first, so floating-point noise cannot create or break a tie or a zero.
  - p = min(1, 2·P(W ≤ min(W⁺, W_total − W⁺))).
- **Robustness to PSI 0/1.** The test is repeated without PSI values of exactly 0 or 1. For pairs, a pair is dropped
  when either value is 0 or 1. Unpaired, the sample values are dropped before each patient's mean is taken. The
  result has four states:
  - `not_applicable`: nothing was dropped;
  - `untestable`: the reduced set is below the minimum;
  - `pass`: p < α, the same sign, and |Δ| > `min_abs_delta`;
  - `fail`: otherwise, also when enough values remain but none differ (all differences zero, or constant values).
- **Hit.** p < `alpha` (0.05), round(|Δ|, 12) > `min_abs_delta` (0.10), sign(HL) = sign(Δ), and a 0/1 state of
  `pass` or `not_applicable`. Otherwise the hit status names the first rule that failed, in this order:
  - `not_tested`
  - `not_significant`
  - `small_effect`
  - `direction_ambiguous`
  - `robustness_fail` or `robustness_untestable`
- **`group_hit`.** The paired or the unpaired design is a hit.
- **Composition diagnostics** (unpaired side). These are descriptive and enter no rule.
  - `matched_*`: the case values of paired patients against all reference values (one value per patient, as above).
  - `other_*`: the remaining case patients against the matched ones.
  - `composition_sensitive`: the paired and unpaired Δ differ by more than `min_abs_delta` or in sign.
- **Within-patient support.** The paired test ran, and either the paired design is a hit or the matched comparison
  is a hit in the unpaired direction.
- **No within-patient check.** A hit in a cohort where the paired test was not possible, because there were too few
  pairs or all differences were zero.

## Survival (per event, cohort and endpoint)

- **Coverage gate.** Applied on the survival samples. PSI must be observed in at least `coverage_frac` (50%) of them,
  and at least `min_off_modal` (10) observed values must differ from the modal value. Otherwise neither test runs
  (`coverage_gate`).
- **KM split.** Set by `km_split` (`--km-split`): the median (default) or the mean of PSI over the survival
  samples, or a value you give (e.g. 0.5). The split is fixed per cohort and shared by every endpoint.
  - PSI ≤ the split is the low arm; PSI > the split is the high arm.
  - A given value applies to every event and cohort. Where it leaves an arm below `km_min_group`, no log-rank test
    runs.
  - `km_split` in the tables is `median`, `mean` or `set`, and `cutoff` is the value.
  - The log-rank test needs `km_min_group` (10) patients per arm and `km_min_events` (10) events.
  - The log-rank HR is (O/E of the high arm) / (O/E of the low arm).
  - Bands are 95% log-log (exponential Greenwood) intervals. The at-risk table counts patients still followed at each
    tick.
- **Cox model.** lifelines `CoxPHFitter`, with Efron ties and no penalty (unless the opt-in ridge below is set).
  - Model: h(t) = h₀(t)·exp(β·PSI/0.10 + γ·z(host expression) + covariates) (`CoxModel`).
  - Host expression enters when an expression table is given. It is z-scored over the fit cohort. Like a clinical
    variable, it is left out of a cohort's model when it is recorded for fewer than `covariate_min_complete` (80%) of
    the fit patients (none when the table has no row for the host gene) or is constant: the model is then PSI alone,
    with a note (`cox_notes`, and the fitted `cox_model`).
  - Numeric clinical covariates enter per SD of the fit cohort, or per unit with `scale=False`.
  - Other covariates enter as indicator variables against a baseline level. That is the level named in
    `CoxModel(baseline=...)` when it occurs in the cohort, else the cohort's most common level.
  - Rare levels are merged before fitting. A level with fewer than `level_min_patients` (10) patients in the fit
    cohort, or with no events, has no stable estimate (with no events its likelihood has no maximum).
    - Stage numerals and numbers merge into the adjacent level with fewer patients, the rarest first. Stage I with 2
      patients joins II as "I–II", and IV with 3 joins III as "III–IV".
    - Other categories merge into the most common level ("white+asian"). When the most common level itself has no
      events, it merges into the next most common ("black+white").
    - A named baseline refers to its merged level ("I–II"). `cox_notes` lists each merge with its patients and
      events. `level_min_patients` = 0 turns merging off.
  - Strata give each level its own baseline hazard.
  - Patients missing any model variable are left out, and the count is reported (`cox_n_dropped`). A covariate that
    is constant in the cohort is dropped, with a note (`cox_notes`).
  - Requirements: `cox_min_n` (30) patients, `cox_min_events` (10) events, and `min_off_modal` PSI values away from
    the mode.
  - **Low power.** A fit with fewer than `cox_low_power_events` (20) events runs, and is marked: a note in
    `cox_notes`, `cox_low_power` in the tables, and ‡ after its CI in the forest. A p ≥ 0.05 there says little
    against an association.
  - **Events per term (opt-in).** `cox_min_events_per_term` (default 0, off): a fit with fewer events per estimated
    term is not run (`too_few_events_per_term`). With few events per term the PSI term of a model with several
    clinical terms is anti-conservative: in two simulations with no true effect, the adjusted model (age, sex, stage)
    gave PSI p < 0.05 in about 7–10% of 40-patient cohorts with 10–14 deaths (the base model about 5%).
  - **Ridge (opt-in).** `cox_ridge` = `clinical`, `molecular` (PSI or the HIT index, and host expression) or `all`
    adds the penalty λ/2·β² (`cox_ridge_penalty`, default λ = 1) to those terms, where β is a log HR per SD of the fit
    cohort (per level, against the baseline, for categories): like a normal prior with SD 1/√λ on each, whatever the
    cohort's size. Strata are not penalized. The penalized terms are shrunk jointly toward HR 1: a single HR, PSI's
    included, can move away from 1 when terms are correlated, and an unpenalized term (PSI under `clinical`) can move
    either way. SEs, CIs and p values come from the penalized information matrix (lifelines), so they are
    approximate. Around the shrunk estimate the Wald p is conservative (the penalized variance is never smaller than
    the estimate's sandwich variance); but when confounders are penalized, the estimate itself keeps part of their
    association, and PSI's p can be anti-conservative (below). A penalized term's CI reads best as an approximate 95%
    posterior (credible) interval for its log HR per SD or per level under that normal prior. A penalized model's
    name says so (`cox_model`, the pages, the probe report). The test suite checks the fit against an independent
    penalized partial likelihood, with and without strata.
    - Trade-offs. `clinical` and `all` shrink the confounders, so they adjust PSI only partly: PSI's HR stays closer
      to the HR without them, and an association that full adjustment would weaken can survive it. In a simulation
      with no PSI effect, 100 patients (about 28 deaths) and PSI 0.24 higher per stage, about 60% (`clinical`) and
      55% (`all`) of the confounding stayed in the PSI estimate, and PSI p < 0.05 went from 6.5% of fits without the
      ridge to 11.5% and 9.8%. With 40 patients (about 13 deaths) and no such correlation, `clinical` brought it from
      6.9% to 5.1%, near the nominal 5%. `molecular` leaves the confounders unpenalized and kept PSI p < 0.05 near or
      below 5%: 4.9% and 5.3% in these simulations, and 3–4% with PSI and host expression correlated (r 0.8).
    - The proportional-hazards tests of penalized terms are conservative (lifelines scales their residuals by the
      penalized variance), so a † is rarer on them.
    - For categories the penalty is on each level against the baseline, so it depends on which level is the
      baseline. λ = 1 is mild for a term the data inform well, but a sparse level leans on the prior: with 40
      patients (about 13 deaths) it shrank the top stage level's log HR by about 38% (median). The prior also keeps
      the se of a penalized term per SD or per level below 1/√λ (1 at the default), so the `unstable` note (se > 3)
      cannot flag a sparse level or a covariate per SD there unless λ < 1/9; a covariate per unit
      (`CoxModel(scale=False)`), whose se is per unit, can still be flagged.
  - HR per SD = exp(β·SD/0.10), where SD is the standard deviation of PSI in the fit cohort (`hr_per_sd`,
    `ci_low_sd`, `ci_high_sd`). Its 95% CI is exp((β ± 1.96·se)·SD/0.10).
  - HR per IQR is the same with the interquartile range of PSI in the fit cohort; it is left empty when the IQR is 0.
    The tables hold both; `psi_hr_unit` ("sd", the default, or "iqr") picks the one the pages, the probe and its
    overview show, and the model term in `cox_terms` (`psi_sd` or `psi_iqr`).
  - **Proportional hazards.** Every fit gets lifelines' `proportional_hazard_test`: the Schoenfeld-residual test on
    the Kaplan–Meier time scale, one test per term, stratified models included.
    - `ph_p` is the PSI term's p; `cox_terms.ph_p` holds every term's.
    - The KM comparison gets the same test on a Cox model of the high-arm indicator (`km_ph_p`).
    - The tests are reported, never used to drop a result (see Notes on a fit).
    - On the pages a dagger (†) follows each test whose p is below `ph_note_below` (0.05): the KM log-rank p, the
      forest's CI and the p of each model row. The legend explains it.
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
- **HIT index apart.** The tests of HIT-index events form families of their own, kind by kind, so including them
  changes no PSI q value. A page's footnote says which events its families cover.
- **Where q is shown.** The pages print q under or beside the p values: the tumour–normal titles, the KM header,
  and the Cox model headers. A footnote names the family.
- **The q mark.** A `*` follows each q below `q_mark_below` (0.05), and the forest's CI and the tested term's p in
  the model rows carry it too. In the probe overview it replaces the p dot. It is separate from the filled markers,
  which stay on p.
- **How the family is built.** `panel` and `panels` analyze all of a gene's events of the page's kind (PSI events,
  or HIT-index events), and `probe` all of the gene's events it probes, so the family is the whole gene; `analyze`
  uses the events it is given.
- **What stays on p.** The hit rules, the forest's filled markers and the probe's ranking use p, as in the screen.
  q is reported beside p.
- **Not adjusted.** Host-gene expression statistics.
- **Why BH within a gene.** The tests of one gene are correlated, mostly positively: rMATS lists the same exon with
  several flanks, and one cohort's events share patients. BH controls the false discovery rate under such positive
  dependence. Treat q as a guide to how much of the gene's signal could be chance, not as a per-test guarantee.

## Host-gene expression (per gene, cohort and endpoint)

With an expression table, the host gene's own expression gets the same three analyses as PSI (`analyze_expression`;
the gene's expression page; the probe's `expression_cells.csv`).
- **Case vs reference.** The same paired signed-rank and Mann–Whitney tests and composition diagnostics, on the
  expression table's scale.
  - There is no 0/1 robustness check, because that is specific to PSI.
  - A hit needs |Δ median| > `gex_min_abs_delta`: 1.0 by default, which is two-fold on a log2 scale.
- **KM.** A split of expression over the cohort's survival samples (expression ≤ the split is the low arm), with the
  same coverage and size gates as PSI. The split is set by `km_split_expression` (`--km-split-expression`): the
  median (default), the mean, or a value on the expression table's scale.
- **Cox.** h(t) = h0(t) exp(b · z(expression) + clinical terms). The clinical terms and strata are those of the
  page's model rows (by default age, sex and stage), with the same completeness rule. The HR is per SD of
  expression in the fit cohort. Status `constant_expression` when expression does not vary.
- **Interpretation.** The splicing models already adjust PSI for expression. The expression page shows whether the
  gene's level itself goes with survival or differs between the groups, beside the splicing results.

## Default clinical adjustment and the probe

- **Found variables.** Age, sex and stage are found by column name and cleaned (`clinical.py`):
  - stage becomes its overall Roman numeral, from Roman numerals ("Stage IIA") or stage numbers (2, "2B");
  - sex becomes female or male;
  - other codes become missing.
- **Adjusted model.** Base model + age (per SD) + sex (against female) + stage (against I). A stage too rare to
  estimate in a cohort joins its neighbour (see the Cox model above), so stage I with 2 patients makes the
  reference I–II.
- **Variables missing in a cohort.** Below `covariate_min_complete` (80%) recorded in a cohort, a variable is left out
  of that cohort's model (`cox_notes`).
- **Unstable terms.** A covariate term with se > 3 (a sparse level, or near separation) is reported and flagged as
  unstable in `cox_notes`.
- **Notes on a fit.** These notes go in `cox_notes` and on the page's model header (the KM one in `km_notes` and on
  the KM header). All are notes only: the test runs and its estimates are unchanged.
  - **Overfit risk:** fewer than `cox_events_per_term` (10) events per estimated term. Terms are the PSI, expression
    and covariate coefficients, one per level beyond the reference; strata do not count. The value is
    `cox_events_per_term` in the tables.
  - **Narrow PSI range:** the PSI spread in the fit cohort is below `narrow_psi_below` (0.05), so the HR
    describes a change of a few PSI points. The spread is the HR's unit (SD, or IQR with `psi_hr_unit = "iqr"`)
    unless `narrow_psi_measure` names one. The flag is `psi_narrow` in the tables.
  - **Non-proportional hazards:** a term's proportional-hazards p is below `ph_note_below` (0.05), so its effect
    changes over follow-up and the HR is an average over it. The KM header says the same when the high/low split
    fails the test. Crossing curves are the usual cause.
  - The probe report counts each note and names the noted fits with p < α.
- **Probe ranking.** For one endpoint, events measurable in at least one cohort (a Cox fit or a log-rank test ran)
  come first. Events are then ranked by these criteria in turn, counting a Cox fit only when it is not low power
  (`cox_low_power_events`, 20 events):
  1. the number of cohorts with an adjusted Cox p < α;
  2. the same for the base Cox model;
  3. the number of cohorts with both a group hit and a survival hit (KM or base Cox p < α);
  4. the number of cohorts with a KM p < α;
  5. the Cox hits in low-power fits (‡): adjusted, then base. They break ties, so a cohort with 10–19 events
     cannot lift an event above one with the same evidence from better-powered fits;
  6. the p of the best cohort: adjusted Cox, else base Cox, else KM (the fallback when no Cox model can be fitted).
     The best cohort is the one with the smallest p < α in a fit that is not low power, else in a low-power fit,
     else the smallest p in a fit that is not low power, then in a low-power fit.

  KM hits count in full whatever a cohort's events: the log-rank test stays valid with few events (`km_min_events`,
  10), while it is the Cox model with several clinical terms that becomes anti-conservative.
- **Probe q values.** The probe analyzes every event of each gene it probes, so its q values are the gene-wide
  families of the previous section; `adj_cox_q` is the adjusted model's Cox family.
- **Cohorts shown per page.** As for the best cohort: Cox p < α in a fit that is not low power (the adjusted model
  where fitted, else the base model), then in a low-power fit, then the other fits that are not low power, then the
  other low-power fits; within each, the smallest adjusted Cox p, then base Cox p, then KM p.

## The figures

- **Title and event details.** The title names the gene, the events, the cohorts and the endpoint. Under it, one
  line per event gives its details (see Event types).
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
- **Event panel: KM.** The split, bands, censoring ticks and the at-risk table. Above the axes:
  - the split, e.g. "split at median PSI 0.7705" or "split at PSI 0.5 (set)";
  - the log-rank HR and p (with † for non-proportional hazards);
  - q, when the family is large enough;
  - the events per arm;
  - with †, a last line giving the proportional-hazards p.
- **Event panel: forest.** One row per cohort where a plotted event had a tested comparison or a Cox fit, plus the
  cohorts shown at left.
  - Left: the paired Δ (filled) and unpaired Δ (open).
  - Right: the HR per SD (or IQR) with its 95% CI, filled when Cox p < α; after the CI a `*` when its q < 0.05, a †
    when its proportional-hazards p < 0.05, and a ‡ when the fit has fewer than 20 events (low power).
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
- **Cox model figure.** Every term of one cell's model, with HR, 95% CI and p. PSI is shown per SD (or IQR); host expression
  and numeric covariates per SD; categories against their baseline level.
- **Expression page.** The host gene's expression is drawn once per gene, after its splicing pages, so the splicing
  pages do not repeat it.
  - At the top, a forest of every cohort with an expression test: Δ of expression (left) and the Cox HR per SD
    (right). The cohorts drawn below are shaded.
  - Then one row per cohort drawn: the group view, the KM split and the Cox model of expression plus the page's
    clinical terms.
  - `panel` and `probe` draw the cohorts of the gene's splicing pages (their union in a probe). `--no-gex` leaves the
    page out.

## Reproducibility

- **Identical outputs.** The same inputs, settings and library versions give byte-identical SVG, PDF, PNG, CSV and
  JSON (tested).
- **Plotted values.** Every figure writes a CSV of each value it draws.
- **Provenance.** Each figure also writes a `.provenance.json` with:
  - the SHA-256 of the exact input rows used, including the PSI of the gene's other events (`psi_q_family`), whose
    tests set the printed q values;
  - the settings;
  - the models: the forest's (`call.model`) and the model rows' (`call.detail_model`);
  - the call;
  - the software versions.
- **Analysis tables.** `analyze` writes `analysis.json` beside its CSV tables: the Cox model and the settings.

## Validation

- **Tests.** The test suite checks each statistic against an independent implementation:
  - brute-force enumeration for the signed-rank p;
  - scipy for the Mann–Whitney p, also on per-patient means computed independently (pandas);
  - lifelines for log-rank, KM and Cox, including models with covariates and strata.
- **Acceptance run.** An acceptance run outside this repository (it uses controlled-access data) fed the package plain
  tables exported from the locked analysis these figures were designed for.
  - It reproduced that analysis's tables: group tests, coverage, KM, and Cox with and without clinical covariates.
  - It also reproduced the plotted values of its figures.
  - The largest difference was 4.4e-16.
  - It predates the rule that a patient's several samples in one group enter the unpaired test as their mean. That
    rule changes nothing when every patient has one sample per group.

## Protein consequences

The optional protein suggestions (transcript matching, effect, features) are described in
[proteins.md](proteins.md).
