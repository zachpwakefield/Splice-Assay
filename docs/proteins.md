# Protein consequences (suggested from annotation)

splice-assay can suggest how an event could change the protein. It finds the annotated transcripts that best
represent the event's two forms, compares their proteins, and lists the protein features at the event. It draws a
protein band under the schematic of each page.

This is a suggestion read from annotation, not a measurement. An event with no good match gets no band and no
suggestion.

## Quick start

```bash
# once: export a SpliceImpactR annotation cache (needs R, any version from 3.5, no packages; about 30 s)
splice-assay protein-cache ~/annotation_cache --out protein_cache/
export SPLICE_ASSAY_PROTEINS=protein_cache/        # or pass --proteins protein_cache/ to each command

splice-assay proteins data/ --gene FNBP1             # the suggestion for every event of a gene, printed
splice-assay probe data/ --gene FNBP1                # pages get the band; events.csv and report.md a Protein column
splice-assay panel data/ --event FNBP1:SE:0020 --out figures/
```

`splice-assay example demo/` writes a small synthetic cache to `demo/proteins/` for trying this out.

In Python:

```python
from splice_assay.protein import ProteinCache, protein_changes

cache = ProteinCache.load("protein_cache/")
table, changes = protein_changes(cache, ds, ["FNBP1:SE:0020"])
sa.event_panel(ds, "FNBP1:SE:0020", ["LUAD"], "OS", proteins=cache, out_dir="figures/")
```

## What you get

| Where | What |
|---|---|
| Page band | One bar per protein-coding form on one residue axis, with the event's residues outlined. Shown: named InterPro domains (one of each set of overlapping entries), disordered regions, TM helices, signal peptides and ELM motifs. A triangle marks where an absent region would sit. |
| `proteins.csv` (probe), `proteins --out` | One row per event: the status, both transcripts (name, biotype, TSL, how they matched), protein lengths, event residues, the effect and the features. |
| `events.csv` (probe) | `protein_status`, `protein_change` (a short label such as `+61 aa, in frame`), `protein_summary`. |
| Page CSV | Rows with `panel = protein`: every drawn isoform, domain, event span and the effect text. |

The **status** is one of these:
- `pair`: both forms have a transcript (at least one coding);
- `inclusion_only` or `exclusion_only`: only one form has an annotated transcript, and it is coding;
- `noncoding`: the transcripts found are not protein-coding;
- `no_transcript`: no transcript represents either form;
- `gene_not_in_cache` or `unsupported_event_type`.

Only the first three are shown.

## How the transcripts are chosen

The matching ports SpliceImpactR's ranked-pair matcher (`get_ranked_pairs`, protocol 1). Each form is a list of
intervals that a transcript must contain and a list it must not touch (1-based, closed):

| Type | PSI form (INC) | Other form (EXC) |
|---|---|---|
| SE | upstream, cassette, downstream | upstream, downstream; not the cassette |
| RI | upstream, intron, downstream, in one exon | upstream, downstream; not the intron |
| A3SS / A5SS | long exon, flanking exon | short exon, flanking exon; not the extension |
| MXE | upstream, PSI exon, downstream; not the other exon | upstream, other exon, downstream; not the PSI exon |

**When a transcript represents a form.** All of these must hold:
- every required interval lies inside one of its exons;
- no exon touches a forbidden interval;
- the splice structure agrees exactly: consecutive exons whose inner boundaries are the event's junctions.

SpliceImpactR's rMATS converter leaves out the flanking exon of A3SS/A5SS events. It is supplied here, so that
junction is checked too.

**Fallback to the event's own exons.** rMATS often names flanking exons that are not the annotated neighbours. When no
transcript matches a form exactly, SpliceImpactR's definition from the event's own exons is used:
- SE and MXE: the exon, with exactly its boundaries, as an internal exon;
- RI: the intron inside one exon;
- A3SS/A5SS: the alternative site as an exon boundary.

The SE and RI exclusion forms have no such definition. `inc_context` and `exc_context` say how each transcript
matched: `junction_chain`, `retained_intron`, `partial_event` or `splice_site_only`.

**Ranking.** Pairs of distinct transcripts are ranked by, in order:
1. how many are protein-coding;
2. the worse transcript-support tier (TSL 1–3, 4–5, unknown);
3. the evidence;
4. the similarity of the coding sequence outside the event;
5. the similarity of the whole coding sequence;
6. the worse TSL, then the summed TSL;
7. the transcript IDs.

Similarity is the Jaccard overlap of coding genomic positions, which is SpliceImpactR's fallback evidence.
SpliceImpactR aligns the sequences when it has them, so a near-tie can rank differently. When only one form has
transcripts, the coding one with the best support is used.

## How the change is described

- **Direct comparison.** When both transcripts are coding and differ only at the event, the two proteins are compared
  directly. "Only at the event" means one of two things:
  - they share every coding position outside it;
  - they share at least 95% of them, the event keeps the frame, and every differing residue lies at the event. Example: "inclusion adds 43 aa in frame (residues 223–266; 1 residue at the
  junction changes)".
- **From the coding coordinates.** Otherwise the residues and the frame come from the coding coordinates of the
  transcript carrying the event. Example: "inclusion adds 61 aa in frame (the exon spans residues 330–390 of
  FNBP1-202)".
  - When only the other form is coding, the coordinates give where the event would sit in it. Example (the synthetic data):
    "retention would add 850 nt after residue 117 of SYN1-201: a frameshift".
  - The frame is the length difference of the two forms modulo 3.
  - A stop codon inside an unannotated insertion cannot be seen, and the sentence says so.
- **Other outcomes.** The event can lie in a UTR (the protein is unchanged) or hold the start or stop codon.
- **Pairs that differ elsewhere.** When the pair also differs elsewhere, the sentence says how much coding sequence
  they share. Below 50%, only the transcript carrying the event is drawn.
- **Features.** A feature on one transcript's event exons is "only with" that form when the other protein has no
  feature of the same name whose genomic span overlaps it. This is SpliceImpactR's `get_domains` rule.
  - InterPro families, superfamilies and whole-protein entries (counting all pieces of an entry) are not named.
  - Pfam and CDD entries are kept in the cache but not named, because InterPro integrates them.

The forms are named by event type:

| Type | PSI form | Other form |
|---|---|---|
| SE | inclusion | skipping |
| RI | retention | splicing |
| A3SS / A5SS | long form | short form |
| MXE | PSI exon | other exon |

## Checks

These were run on 818 rMATS events of eight genes (FNBP1, FGFR2, PKM, BCL2L1, TPM1, ACTN1, CD44, MKNK2) against
GENCODE v45.
- **Against SpliceImpactR.**
  - The compatibility calls agree with SpliceImpactR's `get_ranked_pairs` for all 30,882 transcript × form
    combinations, including the reason given for each rejection.
  - With both forms matched exactly, the same pair is chosen for 43 of 44 events.
  - The one difference is a near-tie that SpliceImpactR breaks by sequence alignment: FGFR2:A5SS:0001, where the
    coordinate similarities are 0.860 and 0.856.
- **Known answers.**
  - FGFR2 IIIb/IIIc (FGFR2:MXE:0006): FGFR2-215 (822 aa) vs FGFR2-206 (821 aa). These are UniProt's KGFR (IIIb) and
    BEK (IIIc) lengths. The exons are swapped in frame inside Ig-like domain III.
  - PKM1/PKM2 (exons 9 and 10): both proteins are 531 aa, swapped in frame at residues 381–436.
  - FNBP1's 183-nt exon (FNBP1:SE:0020): 61 aa in frame at residues 330–390 of FNBP1-202, in a disordered linker.
    No annotated coding transcript skips it.
- **Coverage.** 666 of the 818 events got a suggestion:
  - 154 with both forms;
  - 439 from the PSI form only;
  - 73 from the other form only.

  Of the other 152, 86 are MXE pairs that annotation does not support (both exons in every transcript, or
  neither).

## The cache

`protein-cache` reads a SpliceImpactR annotation cache:
- `human_gencode_v45.gtf.rds` and `human_gencode_v45_sequences.rds`;
- the protein-feature files `interpro.rds`, `pfam.rds`, `cdd.rds`, `elm.rds`, `mobidblite.rds`, `signalp.rds` and
  `tmhmm.rds`.

It writes three tables, as Parquet (with pyarrow) or gzipped TSV. You can supply other annotation in the same format,
as `.parquet`, `.tsv(.gz)` or `.csv`:

| Table | Columns |
|---|---|
| `exons` | `transcript_id`, `transcript_name`, `transcript_type`, `transcript_support_level`, `gene_id`, `gene_name`, `chr`, `strand`, `start`, `end`, `exon_id`, `cds_gen_start`, `cds_gen_stop`, `cds_rel_start`, `cds_rel_stop` |
| `proteins` | `transcript_id`, `protein_seq` (optional `protein_id`) |
| `features` | `transcript_id`, `database`, `name`, `start`, `stop`, `g_start`, `g_end` (optional `feature_id`) |

**Coordinates.**
- Exon coordinates are 1-based and closed.
- `cds_gen_*` is the coding part of the exon, empty when the exon is not coding.
- `cds_rel_*` is the same part in coding-sequence coordinates: 1 is the first coding nucleotide, counted in
  transcript orientation. The stop codon is not part of the coding sequence (the GENCODE convention).
- Feature `start` and `stop` are residues. `g_start` and `g_end` are the feature's genomic span; they tell a
  feature shared by both proteins from one that is gained.

**Matching events to genes.** Transcripts are matched to an event by the events table's `gene_id` (version ignored),
then by `gene`. As a last resort, any transcript on the same strand that overlaps the event is used.

## Limits

- **Annotated isoforms only.** The two forms are matched to annotated isoforms. A form that annotation lacks is
  described from the other form's coordinates, never invented.
- **Exon choice.** The match does not know which exon the reads came from. When annotation offers several exons at
  the same place, the best-supported transcript is used.
- **NMD.** Nonsense-mediated decay is reported only as annotated (biotype `nonsense_mediated_decay`), not
  predicted.
- **Predicted features.** Features are those of the annotation cache: InterPro member databases, MobiDB-lite, TMHMM,
  SignalP and ELM predictions.
