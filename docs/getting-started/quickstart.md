# Quick Start

Get started with open-pacmuci in under 5 minutes. This tutorial walks through the full analysis pipeline.

---

## Prerequisites

- open-pacmuci installed ([Installation Guide](installation.md))
- External tools on PATH: minimap2, samtools, bcftools, Clair3
- PacBio HiFi CCS reads from a MUC1 VNTR PCR amplicon

---

## Full Pipeline (Recommended)

Run all five stages in a single command:

```bash
open-pacmuci run \
  --input reads.fastq \
  --output results/ \
  --clair3-model /path/to/clair3/models/hifi \
  --threads 8
```

**What happens:**

1. Generates a synthetic reference ladder (150 contigs, 1-150 repeat units)
2. Maps reads to the ladder with minimap2
3. Detects two allele lengths from the read count distribution
4. Calls variants per allele with Clair3
5. Builds consensus and classifies each 60bp repeat unit

---

## Step-by-Step Execution

For more control, run each stage individually:

### 1. Generate Reference Ladder

```bash
open-pacmuci ladder --output reference_ladder.fa
```

This creates a FASTA with 150 contigs, each containing 1-150 canonical X repeats plus flanking sequences.

### 2. Map Reads

```bash
open-pacmuci map \
  --input reads.fastq \
  --reference reference_ladder.fa \
  --output-dir results/ \
  --threads 8
```

### 3. Detect Alleles

```bash
open-pacmuci alleles \
  --input results/mapping.bam \
  --output-dir results/
```

### 4. Call Variants

```bash
open-pacmuci call \
  --input results/mapping.bam \
  --reference reference_ladder.fa \
  --alleles-json results/alleles.json \
  --output-dir results/ \
  --clair3-model /path/to/clair3/models/hifi
```

### 5. Build Consensus and Classify

```bash
open-pacmuci consensus \
  --input results/mapping.bam \
  --reference reference_ladder.fa \
  --alleles-json results/alleles.json \
  --output-dir results/

open-pacmuci classify \
  --input results/consensus_allele_1.fa \
  --output-dir results/
```

---

## Input Format

- **PacBio HiFi CCS reads** as FASTQ or BAM
- Reads should be from a PCR amplicon spanning the MUC1 VNTR (primers per Wenzel et al. 2018)
- Minimum recommended coverage: **10x per allele** (higher coverage improves confidence)

!!! tip "Coverage recommendation"
    For clinical-grade results, aim for 50-100x per allele. The `--min-coverage` flag (default: 10) controls the minimum read count threshold for allele detection.

---

## Output Files

| File | Description |
|------|-------------|
| `alleles.json` | Detected allele lengths, read counts, contig assignments |
| `allele_1/variants.vcf.gz` | Filtered Clair3 variants for allele 1 |
| `allele_2/variants.vcf.gz` | Filtered Clair3 variants for allele 2 |
| `consensus_allele_1.fa` | Consensus FASTA (VNTR region only) |
| `consensus_allele_2.fa` | Consensus FASTA (VNTR region only) |
| `repeats.json` | Per-repeat classification with confidence scores |
| `repeats.txt` | Human-readable VNTR structure string |
| `summary.json` | Combined allele + classification + mutation report |

### Example Output

```
results/
├── alleles.json
├── allele_1/
│   └── variants.vcf.gz
├── allele_2/
│   └── variants.vcf.gz
├── consensus_allele_1.fa
├── consensus_allele_2.fa
├── repeats.json
├── repeats.txt
└── summary.json
```

### Reading Results

```bash
# View allele lengths
cat results/alleles.json | python -m json.tool

# View VNTR structure
cat results/repeats.txt
# allele_1: 1 2 3 4 5 C F X X X X:dupC B A A B X X X V 6 7 8 9
# allele_2: 1 2 3 4 5 C F X X X X X B A A B X X X V 6 7 8 9

# View confidence scores
python -c "import json; d=json.load(open('results/repeats.json')); print(d['allele_1']['allele_confidence'])"
```

---

## Next Steps

- **[Core Concepts](concepts.md)** -- Understand the pipeline architecture
- **[Deviations from PacMUCI](deviations.md)** -- What changed and why
- **[CLI Reference](../reference/cli.md)** -- All command options
- **[Benchmarking](../guides/benchmarking.md)** -- Validate with simulated data
