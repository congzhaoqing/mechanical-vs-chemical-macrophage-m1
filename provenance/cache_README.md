# Expression-matrix cache — Project 2 (mechanical vs chemical M1)

Written 2026-11-02. Purpose: never re-download or re-quantify these inputs again.
The sandbox workspace is swept when idle; this directory is on the local disk and is not.

## Why this exists
The processed expression matrices from the original analysis were lost when the execution
kernels were reset. The derived results survived as artifacts (Tables S4-S11), but the
matrices they were computed from did not, which means any NEW module (e.g. the NF-kB/AP-1
effector arm) would have required re-fetching every dataset. This cache removes that cost.

## Coverage
All 20 accessions used anywhere in the manuscript. See MANIFEST.csv for the per-file index
(file, kind, size, accession, role).

- main pool (K = 10): 11 accessions
- pre-specified sensitivity set: GSE283034, GSE299513
- chemical and mechanical anchors: GSE69607, GSE253423, GSE274020, GSE46903, GSE99445, GSE99446
- modality projection: GSE275018 (three-arm magnetomechanics, Figure 8)

## File kinds
- `processed matrix` - the deposited count/TPM matrix, used as-is
- `raw archive` (*_RAW.tar) - per-sample files; GSE36878 / GSE69605 / GSE99445 / GSE99446 are CEL
  and need normalisation before use
- `series matrix` - GSE99445 / GSE99446 were processed from the series matrix with the custom CDF
  hugene10st_Hs_ENTREZG (probe ids are Entrez gene ids, mapped via NCBI Homo_sapiens.gene_info)
- `GSE299513_kallisto/` - 43 author-deposited kallisto abundance files; `GSE299513_gene_tpm.csv`
  is the aggregated gene-level matrix derived from them
- `E-MTAB-7900_run_urls.csv` + `E-MTAB-7900.sdrf.txt` - the 48 single-end ENA runs (6.1 GB);
  no processed matrix is deposited for this study, it must be quantified from reads

## IMPORTANT - E-MTAB-7900 assembly caveat
The published ISG effect for this study (delta = +0.0498, the largest positive value in the main
pool) was computed from a salmon quantification whose transcriptome release is not recorded.
The re-quantification stored here uses Ensembl GRCr8 (current). Numbers from the two
quantifications must NOT be mixed. If the effector arm is added, re-score BOTH the ISG modules
and the effector modules on the same matrix so the paired within-study comparison is internally
consistent - which the paired-contrast design requires anyway.

## Processing conventions used throughout (reproduce these)
counts -> log2 CPM; TPM/FPKM and normalised array intensities -> log2(x+1); already-logged
matrices used as supplied. Drop genes detected in fewer than max(2, floor(0.25 * n_samples))
samples. Collapse duplicate symbols to the row of highest mean expression. No batch correction
anywhere - the singscore statistic is rank-based and single-sample.

`genesets.csv` and `score_lib.R` are the gene-set definitions and the R scoring helpers
(singscore 1.26.0 simpleScore, centerScore = TRUE).

## Validation of the E-MTAB-7900 re-quantification (performed 2026-11-02)

The matrix in this cache was re-quantified from the 48 ENA runs against Ensembl GRCr8, whereas
the published numbers came from a salmon run whose transcriptome release was not recorded. The
regenerated matrix reproduces the published dose trend on all four modules:

| module | published slope/decade | regenerated (GRCr8) | published endpoint | regenerated endpoint |
|---|---|---|---|---|
| ISGF3        | +0.0142 (SE 0.0062) | +0.0135 (SE 0.0052) | +0.0498 | +0.0485 |
| Hallmark IFN-alpha | +0.0125        | +0.0117             | -       | +0.0433 |
| Hallmark IFN-gamma | +0.0106        | +0.0077             | -       | +0.0297 |
| YAP/TAZ      | -0.0006             | -0.0013             | -       | -0.0213 |

Sign and magnitude agree throughout; the ISGF3 endpoint agrees to 2.6%. Spearman rho for ISGF3
is 0.618 here versus 0.558 published.

CAVEAT ON THE p VALUE. The validation fit above is a plain lm of score on log10(stiffness) over
the 16 samples (df = 14) and returns p = 0.022 for ISGF3. The published p_trend = 0.108 comes
from the donor-aware test over 4 donors. The smaller p here is a degrees-of-freedom artefact of
the quick check, NOT a strengthening of the result: the donor structure must be respected in any
analysis that uses this matrix. Slope and endpoint, which do not depend on df, are what match.

Design recovered exactly: 4 stiffness levels (0.1 / 1 / 10 / 100 kPa) x 4 samples each = 16
samples, from 48 runs (3 sequencing runs per sample, summed at the read-count level before
normalisation). Median mapping rate 40.6% (range 34.2-51.1%). 22,632 genes; 84/118 ISGF3
orthologues recovered in rat.
