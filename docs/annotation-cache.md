# Getting the annotation cache from SpliceImpactR

The protein suggestions need three things: annotated transcripts, their protein sequences, and protein features.
They appear in two places: the band under the schematic, and the Protein column of a probe.
[SpliceImpactR](https://github.com/fiszbein-lab/SpliceImpactR) downloads and prepares this annotation, and
splice-assay converts it once into its own protein cache.

Without a cache, everything else works and the protein band is left out. How the matching works and what the band
shows is in [proteins.md](proteins.md).

## 1. Install SpliceImpactR (R)

It installs from Bioconductor or GitHub (see its
[installation notes](https://github.com/fiszbein-lab/SpliceImpactR#installation)):

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("SpliceImpactR")
# or: devtools::install_github("fiszbein-lab/SpliceImpactR")
```

## 2. Prepare the annotation and save it (R)

This step needs internet access and takes a while. It downloads GENCODE once (about 110 MB of GTF and FASTA) and
queries Ensembl BioMart, and ELM for motifs. The script saves the objects splice-assay reads into one folder:

```r
library(SpliceImpactR)
out <- path.expand("~/splice_assay_annotation")    # any folder
dir.create(out, showWarnings = FALSE)

ann <- get_annotation(load = "link", species = "human", release = 45, base_dir = out)
saveRDS(ann$annotations, file.path(out, "human_gencode_v45.gtf.rds"))
saveRDS(ann$sequences, file.path(out, "human_gencode_v45_sequences.rds"))

for (db in c("interpro", "pfam", "cdd", "elm", "mobidblite", "signalp", "tmhmm")) {
  saveRDS(get_protein_features(db, ann$annotations, base_dir = out), file.path(out, paste0(db, ".rds")))
}
```

- **Why save explicitly.** SpliceImpactR keeps its own processed objects inside a BiocFileCache database. The
  `saveRDS` lines write them under plain names that splice-assay can read with any R.
- **File names.**
  - The two annotation files must share a prefix: `<name>.gtf.rds` and `<name>_sequences.rds`.
  - Any name works, e.g. `mouse_gencode_vM34` with `species = "mouse"`.
  - Keep one annotation per folder.
- **Feature databases.**
  - splice-assay reads whichever of these seven files are present, and works without any (no feature boxes).
  - `get_protein_features()` uses Ensembl release 109 unless told otherwise. Pass `release = 111` to match
    GENCODE 45.
- **Transcripts kept.** By default SpliceImpactR keeps transcripts with support level 1–3. `filter_tsl = NULL` keeps
  all of them (see SpliceImpactR's
  [annotation notes](https://github.com/fiszbein-lab/SpliceImpactR#load-reference-resources)).
- **Genome build.** Events are matched by coordinates, so the annotation must be on your events' genome build
  (GRCh38 for current human GENCODE releases). The release itself may differ from the one used to quantify PSI.
- **Later sessions.** `get_annotation(load = "cached", base_dir = out)` reloads without downloading.
- **Older caches.** A folder from an earlier SpliceImpactR version that already holds
  `human_gencode_v45.gtf.rds`, `human_gencode_v45_sequences.rds` and the feature files works as it is.

## 3. Convert it for splice-assay (shell)

```bash
splice-assay protein-cache ~/splice_assay_annotation --out ~/splice_assay_proteins
```

This needs `Rscript` once (any R, no packages; `--rscript PATH` when it is not on the PATH). It writes three tables
(exons, proteins, features) as Parquet with pyarrow, else gzipped TSV, plus `SOURCE.txt`. The human GENCODE 45 cache
is about 70 MB.

## 4. Use it

```bash
export SPLICE_ASSAY_PROTEINS=~/splice_assay_proteins
```

Or pass `--proteins ~/splice_assay_proteins` to `panel`, `panels`, `probe` or `proteins`.
