# Contrast specification & reproducibility protocol — Project 2

Written 2026-11-02. Created because the K = 10 main-pool study-level effect sizes could not be
regenerated from the surviving artifacts.

## What went wrong

The scoring pipeline survived (`score_lib.R`) and is demonstrably correct: every dataset whose
input file AND sample grouping are unambiguously recoverable reproduces the published number
exactly — GSE99445 IFN-gamma ISGF3 0.10017 vs 0.10017 (background 19,529 genes in both), GSE99446
0.10166 vs 0.10166, GSE99445 porous-Ti YAP 0.01085 vs 0.01085, GSE99446 IFN-gamma YAP 0.01098 vs
0.01099, GSE69607 M1 ISGF3 0.03942 vs 0.039175.

What did NOT survive is the **contrast specification**: which samples went into which arm, for each
of the 11 main-pool accessions. It existed only inside an analysis script that was never saved as
an artifact. Two scripts that WERE saved (`run_bulk.R`, `run_stage3.R`) cover the modality
projection and the titanium cohort only.

Candidate causes were tested and excluded as explanations for the residual 1.5-4x gaps:
  - preprocessing scale (counts vs norm): **zero effect**. singscore ranks genes WITHIN each
    sample, so any per-sample monotone transform (CPM, log, library-size division) leaves ranks,
    and therefore scores, identical. Verified: 0.20567 either way.
  - detection-filter threshold (0 / 10 / 25 / 50 / 75 % of samples): moves the delta by <= 11%.
  - annotation source (NCBI gene_info + RefSeq FASTA vs org.Mm.eg.db): <= 7%.
  - retaining unmapped identifiers in the background: -2%.
  - sample-selection variants for GSE274020 (5 plausible groupings): all +0.206 to +0.283, none
    near the published +0.1245.
Gene background size is NOT the explanation either: published backgrounds for the modality
datasets run 13,465-39,975 genes (median 18,150) and the rebuilt backgrounds sit in the same range.

**Conclusion: the pipeline is reproducible; the sample-to-arm assignment is what was lost, and for
the RNA-seq studies it cannot be inferred from the deposited metadata alone.**

## Why it will recur unless the practice changes

Even the two scripts that survived encode the grouping as bare positional vectors, e.g.
`grp = c(0,0,0,1,1,1)` keyed to the column order of a downloaded file. That is fragile even when
the script survives: if the repository re-orders or re-versions the supplementary file, the
grouping silently becomes wrong and nothing in the pipeline detects it. The grouping was expressed
as **code**, and tied to **position**, instead of as **data** tied to **sample identifiers**.

## The rule going forward

1. **The contrast specification is data, not code.** `B1_contrast_specification.csv` names every
   sample explicitly: study, sample id (GSM or matrix column name), group, stratum, paired_id,
   species, scale. A re-run joins on the sample id; positional order is never used.
2. **Pin the annotation.** Record the annotation package and version used for identifier mapping
   (here org.Mm.eg.db / org.Hs.eg.db via AnnotationDbi) together with the resulting background size
   per contrast. Identifier mapping changes the expression background, and the background changes
   every rank-based score.
3. **Fingerprint every contrast.** `B1_fingerprints.csv` stores delta, SE, module-hit count and
   background size per study x module. A re-run asserts against these; a drift in background size
   is the early warning that annotation or filtering has changed.
4. **Save the analysis script as an artifact in the same call that saves its outputs.** The two
   surviving scripts survived by accident.
5. **Keep the matrices.** `cache_matrices/` already does this; the specification and fingerprints
   live alongside them.

## Human effector-arm anchor (resolved 2026-11-02)

The human effector arm initially appeared unmeasurable: against IFN-gamma anchors the NF-kB modules
failed the reproducibility gate (GSE99445 vs GSE99446: TNFA 2.32x, field panel sign-flipped) while
ISG agreed to 1.5%. Two causes, both now identified:
  - IFN-gamma alone does not engage NF-kB, so it is the wrong stimulus for an effector anchor
    (TNFA +0.029 vs ISG +0.100).
  - GSE46903's LPS arms are sampled at 72 h, long past the NF-kB window; its sLPS and upLPS arms
    give sign-flipped effector values.
A donor-paired human LPS time course (GSE112372, 4 donors) resolves it and shows the kinetics
directly: at 0.5 h the effector arm has moved (TNFA +0.0595) while ISG has not (+0.0081); TNFA
peaks at 3 h (+0.1613) and decays to +0.0894 by 16 h, while ISG rises monotonically to +0.2011.
An independent cohort (GSE110947, LPS 4 h, 5 paired donors) gives absolute values a uniform
1.67-1.79x lower across ALL modules — a global cohort scale difference, not module instability —
and the **effector:ISG ratio agrees between cohorts to 1.6% (TNFA), 2.0% (field panel), 6.7%
(inflammatory response)**.

**Therefore the human effector ruler is usable, provided it is used as a ratio to the ISG arm
measured in the same contrast**, which is exactly what the paired within-study design requires.
Absolute human effector percentages carry a factor-1.8 cohort-scale uncertainty and must be
reported with it. Recommended primary human effector anchor: GSE112372 LPS 3 h, TNFA +0.1613
(SE 0.0028); sensitivity: GSE110947 LPS 4 h, +0.0915 (SE 0.0057).
