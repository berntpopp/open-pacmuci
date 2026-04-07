# Known Limitations

## Allele Length Detection

### Homozygous Same-Length Alleles

When both alleles have the same repeat count (e.g., 60/60), the indel-valley splitting algorithm may incorrectly separate reads into two clusters. The pipeline reports `same_length: true` and uses a disambiguation strategy, but accuracy is reduced compared to heterozygous samples with distinct allele lengths.

### Extreme PCR Bias (Asymmetric Alleles)

For highly asymmetric allele pairs (e.g., 25/140), the shorter allele amplifies much more efficiently during PCR. The longer allele may have very low read coverage (<10x), making allele detection unreliable. The pipeline requires a minimum coverage threshold (default: 10 reads) per allele.

## Variant Calling

### Long VNTR Alleles (>100 repeats)

Clair3 has reduced sensitivity for detecting 1bp insertions within very long tandem repeat regions (>6kb). For alleles with >100 repeat units, single-nucleotide frameshifts (e.g., dupC) may not be detected by Clair3, leading to false negatives.

**Affected test samples:**
- `sample_dupc_100_120`: dupC mutation in 100-repeat allele not detected
- `sample_long_120_140`: mutations in 120+ repeat alleles not detected

### Boundary Repeat Mutations

Mutations detected in the last 3 repeat units of an allele receive a boundary penalty (0.5x confidence multiplier) because these positions are prone to alignment artifacts at the VNTR boundary.

## Classification

### Novel Repeat Types

Repeat sequences with >2 base substitutions from any known type are classified as `novel_repeat`. These may represent true biological variation or alignment/sequencing errors. The edit distance and identity percentage are reported for manual review.

### Bidirectional Classification Fallback

When the forward classification pass encounters a large indel (>30bp offset from the expected 60bp unit length), the algorithm falls back to bidirectional classification. This backward anchoring from the 3' end recovers downstream repeats but may produce different results than purely forward classification for the ambiguous region.
