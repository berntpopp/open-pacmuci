# Benchmark Results: Soft QUAL Scoring (v0.3.0)

Date: 2026-04-07

## Overview

Full pipeline benchmark on 44 simulated MucOneUp samples with the new soft QUAL
scoring system (`min_qual=5.0`, continuous QUAL-to-confidence mapping, boundary
repeat penalty).

## Test Data Generation

All samples generated with [MucOneUp](https://github.com/berntpopp/muconeup)
v0.44.2 using `scripts/generate_testdata.py`. Reads are simulated PacBio HiFi
CCS amplicon reads via pbsim3.

**Default settings:** 200x coverage, PacBio HiFi platform, 2 haplotypes per sample.

### Sample Definitions (26 samples from generate_testdata.py)

| Sample | Alleles | Mutation | Targets | Seed | Coverage | Purpose |
|--------|---------|----------|---------|------|----------|---------|
| sample_dupc_60_80 | 60/80 | dupC | 1,25 | 1001 | 200x | Baseline dupC |
| sample_dupa_60_80 | 60/80 | dupA | 1,25 | 1002 | 200x | Baseline dupA |
| sample_insg_60_80 | 60/80 | insG | 1,25 | 1003 | 200x | Baseline insG |
| sample_dupcccc_60_80 | 60/80 | insCCCC | 1,25 | 1004 | 200x | 4bp insertion |
| sample_del_60_80 | 60/80 | del18_31 | 1,25 | 1005 | 200x | 14bp deletion |
| sample_normal_60_80 | 60/80 | normal | -- | 1006 | 200x | Wild-type control |
| sample_homozygous_60_60 | 60/60 | dupC | 1,25 | 1007 | 200x | Same-length alleles |
| sample_asymmetric_25_140 | 25/140 | dupC | 1,10 | 1008 | 200x | Extreme PCR bias |
| sample_short_25_30 | 25/30 | dupC | 1,10 | 1009 | 200x | Short alleles |
| sample_long_120_140 | 120/140 | dupC | 1,50 | 1010 | 200x | Long alleles |
| sample_dupc_60_80_s2 | 60/80 | dupC | 1,25 | 2001 | 200x | QUAL variance seed 2 |
| sample_dupc_60_80_s3 | 60/80 | dupC | 1,25 | 3001 | 200x | QUAL variance seed 3 |
| sample_dupc_60_80_s4 | 60/80 | dupC | 1,25 | 4001 | 200x | QUAL variance seed 4 |
| sample_dupc_60_80_s5 | 60/80 | dupC | 1,25 | 5001 | 200x | QUAL variance seed 5 |
| sample_dupc_40_50 | 40/50 | dupC | 1,15 | 1011 | 200x | Short allele lengths |
| sample_dupc_80_100 | 80/100 | dupC | 1,35 | 1012 | 200x | Medium allele lengths |
| sample_dupc_100_120 | 100/120 | dupC | 1,45 | 1013 | 200x | Long allele lengths |
| sample_dupc_50_55 | 50/55 | dupC | 1,20 | 1014 | 200x | Close pair (gap=5) |
| sample_dupc_50_57 | 50/57 | dupC | 1,20 | 1015 | 200x | Close pair (gap=7) |
| sample_dupc_50_60 | 50/60 | dupC | 1,20 | 1016 | 200x | Close pair (gap=10) |
| sample_normal_50_55 | 50/55 | normal | -- | 1017 | 200x | Close pair control |
| sample_normal_50_60 | 50/60 | normal | -- | 1018 | 200x | Close pair control |
| sample_dupa_100_120 | 100/120 | dupA | 1,45 | 1019 | 200x | Long allele dupA |
| sample_insg_100_120 | 100/120 | insG | 1,45 | 1020 | 200x | Long allele insG |
| sample_normal_100_120 | 100/120 | normal | -- | 1021 | 200x | Long allele control |
| sample_dupc_60_80_cov50 | 60/80 | dupC | 1,25 | 1022 | 50x | Low coverage |

### Additional samples (from previous benchmarking)

| Sample | Alleles | Mutation | Purpose |
|--------|---------|----------|---------|
| sample_bench_5000-5004 | 60/80 | dupC | Manuscript benchmark (5 seeds) |
| sample_gap3-12_* | various | normal | Indel-valley splitting stress |
| sample_close_51_58 | 51/58 | dupC | Close allele pair |

## Results Summary

| Metric | Count | Rate |
|--------|-------|------|
| TP (template match) | 24 | -- |
| TP (partial) | 1 | -- |
| **Total detected** | **25/28** | **89.3%** |
| FN (missed) | 3 | -- |
| FP (false alarm) | 1 | -- |
| **TN (correct normal)** | **15/16** | **93.8%** |
| Errors | 0 | -- |

**Sensitivity:** 89.3% (25/28 mutant samples)
**Specificity:** 93.8% (15/16 normal samples)

## Per-Sample Results

### True Positives (24 template matches + 1 partial)

| Sample | Mutation | Status | Notes |
|--------|----------|--------|-------|
| sample_dupc_60_80 | dupC | TP | Baseline |
| sample_dupa_60_80 | dupA | TP | |
| sample_insg_60_80 | insG | TP | |
| sample_del_60_80 | del18_31 | TP | |
| sample_asymmetric_25_140 | dupC | TP | Despite extreme PCR bias |
| sample_short_25_30 | dupC | TP | Short alleles |
| sample_dupc_60_80_s2-s5 | dupC | TP | All 4 seeds |
| sample_dupc_40_50 | dupC | TP | Short allele combo |
| sample_dupc_80_100 | dupC | TP | Medium allele combo |
| sample_dupc_50_55 | dupC | TP | Close pair (gap=5) |
| sample_dupc_50_57 | dupC | TP | Close pair (gap=7) |
| sample_dupc_50_60 | dupC | TP | Close pair (gap=10) |
| sample_dupc_60_80_cov50 | dupC | TP | At only 50x coverage |
| sample_dupa_100_120 | dupA | TP | Long allele |
| sample_insg_100_120 | insG | TP | Long allele |
| sample_bench_5000-5004 | dupC | TP | **5004 previously FN (QUAL=11.23), now rescued** |
| sample_close_51_58 | dupC | TP | Close pair |
| sample_dupcccc_60_80 | insCCCC | TP_partial | Detected as frameshift, not via insCCCC template |

### False Negatives (3)

| Sample | Mutation | Root Cause |
|--------|----------|------------|
| sample_homozygous_60_60 | dupC | Indel-valley splitting incorrectly splits the single 60/60 peak into two sub-peaks (contig_51 + contig_54 with 723/276 reads). Each sub-allele Clair3 run finds 0 variants. **Bug in allele detection for true homozygous samples.** |
| sample_dupc_100_120 | dupC | Clair3 produces 0 variants on 91-repeat and 111-repeat contigs (~6-7kb). A single 1bp insertion is invisible to Clair3 in very long tandem repeats. **Clair3 limitation.** |
| sample_long_120_140 | dupC | Same as above: 0 variants on 111-repeat and 131-repeat contigs. **Clair3 limitation on long VNTRs.** |

### False Positives (1)

| Sample | Details |
|--------|---------|
| sample_gap4_40_44 | 5 spurious mutations on allele_2, none template-matched. **4/5 have confidence < 0.3** (scoring system correctly penalizes). One mutation at repeat 41 (closest=6p, vcf_qual=19.25, conf=0.96) is a classification error at the after-repeat boundary. |

## Confidence Scoring Effectiveness

The new continuous QUAL scoring and boundary penalty work as designed:

- **sample_bench_5004 rescued:** QUAL=11.23 dupC now passes min_qual=5.0, detected with
  moderate confidence weight (~0.71) rather than being hard-filtered
- **4/5 FP mutations penalized:** confidence < 0.3 due to no VCF support (vcf_score=0.3)
- **Boundary penalty working:** FP repeat 42 (last repeat) correctly flagged as boundary

## Known Limitations

1. **Homozygous same-length alleles:** Indel-valley splitting can incorrectly split
   a single peak, preventing Clair3 from calling hom-alt variants
2. **Long VNTRs (>100 repeats):** Clair3 cannot reliably detect 1bp insertions in
   contigs >6kb of tandem repeats
3. **insCCCC template:** 4bp insertion detected as generic frameshift rather than
   exact template match (template catalog issue)
4. **After-repeat boundary FPs:** Classification at the after-repeat boundary (6p/7/8/9)
   can produce false mutations not caught by the current boundary penalty
