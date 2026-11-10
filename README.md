# Mechanical and chemical macrophage activation differ by arm, not by amplitude

Analysis code and result tables for the transcriptomic meta-analysis asking whether
mechanically induced macrophage M1-like states are the same transcriptional program as
chemically induced (LPS + IFN-γ) M1.

Zhong C, Chen W, Zeng T, Zhang S, Zhu J, Cong Z.
Department of Urology and Stony Brook Cancer Center, Renaissance School of Medicine,
Stony Brook University, Stony Brook, NY, USA.
Correspondence: Zhaoqing.Cong@stonybrookmedicine.edu

## What this repository contains

| folder | contents |
|---|---|
| `analysis_code/` | scoring library, bulk and single-cell pipelines, manuscript build script |
| `analysis_code/diagnostics/` | the 33 scripts used for the reproduction audit and the recomputation under the single scoring convention |
| `provenance/` | contrast specifications, anchor catalogue, gene-set definitions, reproduction checks, expression-matrix cache manifest |
| `result_tables/` | intermediate and final result tables behind every reported number |
| `supplementary_tables/` | Supplementary Tables S1-S15 as deposited with the manuscript |

No primary data were generated. Every dataset analysed is public; accessions and the
per-contrast design are listed in `provenance/B1_contrast_specification.csv` and
Supplementary Table S3.

## What is not here

Processed expression matrices (about 2.1 GB across 130 files) are too large for a code
repository. `provenance/cache_MANIFEST.csv` lists every matrix with its accession, source
URL and file size, so the cache can be rebuilt from the public repositories;
`provenance/cache_README.md` documents how. E-MTAB-7900 has no processed matrix and was
quantified from the 48 single-end ENA runs with salmon against Ensembl GRCr8.

## Pipeline

Module scores use singscore: `simpleScore(rankGenes(X), upSet = module, centerScore = TRUE)`,
with a module required to have at least 5 detected genes. The study effect is the mean of
per-stratum (or per-donor) differences. Module scores are computed on the unfiltered gene
universe; differential-expression counts are computed on the detection-filtered matrix,
because 18-36% of genes on an unfiltered matrix have exactly zero variance and limma's
moderated variances are then unreliable. `score_lib_nofilter.R` is the current library;
`score_lib.R` is the earlier filtered-universe version kept for the reproduction audit.

Meta-analysis is DerSimonian-Laird random effects with confidence intervals on t with K-1
degrees of freedom, two-one-sided-tests equivalence, Kish effective study count, and an
exhaustive 2^K study-level sign-flip permutation.

## Running the code

The scripts locate their inputs through two environment variables, each falling back to a
directory of the same name next to the script tree:

| variable | default | contents |
|---|---|---|
| `M1_CACHE` | `cache_matrices` | processed expression matrices and the gene-set definitions (`genesets.csv`, `genesets_effector.csv`) |
| `M1_PROVENANCE` | `provenance` | contrast specifications and anchor catalogue (included in this repository) |

`genesets.csv` is also included under `provenance/`, so pointing `M1_CACHE` at a directory
holding the rebuilt matrices and copying that file across is enough to run the scoring
pipeline:

```bash
export M1_CACHE=/path/to/cache_matrices
Rscript analysis_code/run_bulk_nofilter.R
```

Rebuild the matrix cache from the public repositories using
`provenance/cache_MANIFEST.csv` and `provenance/cache_README.md`.


## Software

R 4.4.3 (singscore 1.26.0, limma 3.62.1, GSEABase 1.68.0, msigdbr, babelgene,
org.Mm.eg.db, org.Hs.eg.db, affy); Python 3.11 (numpy, scipy, pandas, statsmodels,
matplotlib, scanpy); salmon 2.9.1.

## License

Code is released under the MIT License. Result tables are released under CC BY 4.0.
