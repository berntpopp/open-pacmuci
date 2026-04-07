# open-pacmuci

An open-source reconstruction of the PacMUCI bioinformatics pipeline for analyzing PacBio HiFi long-read sequencing data of the MUC1 VNTR region. The original pipeline was described in [Vrbacka et al. 2025](https://doi.org/10.1101/2025.09.06.673538) but its source code was never released. This project rebuilds it from the published methods, extending it with several algorithmic improvements.

## Background

Autosomal Dominant Tubulointerstitial Kidney Disease caused by MUC1 mutations (ADTKD-MUC1) is a rare genetic kidney disease caused by frameshift mutations in the variable number tandem repeat (VNTR) region of the MUC1 gene. The VNTR consists of 20-125 copies of a degenerate 60-bp repeat unit with extremely high GC content (>80%), making it inaccessible to standard short-read sequencing. Long-read sequencing (PacBio SMRT) combined with specialized bioinformatics is required to resolve the full VNTR structure and locate pathogenic frameshift mutations.

The most common mutation, 59dupC, duplicates a cytosine in the heptanucleotide C-tract of a canonical X repeat unit, producing a +1 frameshift that leads to expression of the toxic MUC1fs protein. At least 9 distinct frameshift mutation types have been identified across affected families.

## How It Works

The pipeline takes PacBio HiFi CCS reads (FASTQ or BAM) from a PCR amplicon spanning the MUC1 VNTR and produces a per-allele VNTR structure report with mutation calls and confidence scores.

### Pipeline Stages

```
Input reads (FASTQ/BAM)
  |
  v
1. LADDER MAPPING ---- Map reads to synthetic reference (150 contigs,
  |                     each with 1-150 canonical X repeats + flanking)
  v
2. ALLELE DETECTION --- Identify two allele lengths from read count
  |                     distribution across contigs (with indel-valley
  |                     splitting for close alleles)
  v
3. VARIANT CALLING ---- Clair3 deep neural network variant caller
  |                     per allele, with VCF quality filtering
  v
4. CONSENSUS ---------- bcftools consensus to build per-allele FASTA,
  |                     anchor-based flanking trim
  v
5. CLASSIFICATION ----- Classify each 60bp repeat unit using Vrbacka
                        nomenclature, detect frameshift mutations via
                        mutation template matching + confidence scoring
```

### Output

For each allele, the pipeline reports:

- **VNTR structure** as a sequence of repeat type identifiers (e.g., `1 2 3 4 5 C F X X X X:dupC B A A B X X X V 6 7 8 9`)
- **Mutation calls** with exact repeat position and mutation type (e.g., `X:59dupC at repeat 25`)
- **Confidence scores** per repeat (1.0 = exact template match) and per allele (mean confidence)

## Deviations from the Original PacMUCI Pipeline

open-pacmuci follows the architecture described in Vrbacka et al. 2025 (Figure 1) but makes several deliberate changes:

| Aspect | PacMUCI (Vrbacka et al.) | open-pacmuci | Rationale |
|--------|--------------------------|--------------|-----------|
| **Read mapper** | bwa-mem 0.7.16a | minimap2 2.28+ | minimap2 is the standard for long-read alignment; faster and more accurate for HiFi data |
| **Reference ladder** | 120 contigs (1-120 X repeats) | 150 contigs (1-150 X repeats) | Covers the full observed VNTR length range (Vrbacka reports alleles up to 125 repeats) |
| **Variant caller** | Clair 2.0.7 + Clair3 1.0.10 | Clair3 only (latest HiFi model) | Clair is deprecated; Clair3 subsumes its functionality with better accuracy |
| **Read filtering** | Manual IGV inspection for chimeric/incomplete reads; MAPQ > 4 filter | Automated read extraction per cluster + remapping to peak contig | Reproducible; no manual step required |
| **Allele detection** | `peaks.py` (idxstats peaks, details unpublished) | Gap-based clustering + **indel-valley splitting** | The indel-valley algorithm resolves close allele pairs (3-9 repeats apart) that simple peak-finding merges into one |
| **Mutation detection** | Custom Python scripts (unpublished) | Pre-computed mutation template catalog (13 known mutations) with exact-match-first probing | O(1) lookup for known mutations; no edit distance needed for common cases |
| **Flanking trim** | Fixed-position trim (assumes flanking is invariant) | Anchor-based trim using flanking/pre-repeat boundary sequences | Resilient to indels in flanking region caused by Clair3 false positive calls |
| **Same-length alleles** | Treated as homozygous (Clair3 skipped) | Disambiguated using Clair3 het genotype detection | Detects compound heterozygotes where both alleles have the same repeat count |
| **Confidence scoring** | Not reported | Per-repeat and per-allele confidence with VCF cross-validation | Enables automated QC; flags low-confidence calls |
| **Source code** | Proprietary (unpublished) | MIT license, fully open source | Reproducibility and community contribution |

### Allele Detection: Indel-Valley Splitting

The key algorithmic innovation in open-pacmuci is **indel-valley allele splitting**, which addresses a fundamental limitation of the ladder-mapping approach.

When two alleles differ by fewer than ~10 repeats, reads from both alleles map to overlapping sets of contigs, creating a single broad peak in the idxstats read count distribution instead of two distinct peaks. Gap-based clustering (used by the original PacMUCI) cannot separate them.

open-pacmuci solves this by analyzing the CIGAR indel lengths of individual reads:

- Reads aligned to the **correct-length contig** have near-zero indels (the read length matches the reference contig length).
- Reads aligned to a **wrong-length contig** have large indels proportional to the length mismatch (~60bp per repeat of difference).

By computing the mean indel length per contig, the algorithm finds two **local minima (valleys)** in the indel series. These valleys correspond to the two true allele lengths, even when the read count distribution shows a single merged peak.

**Benchmark:** 20/20 simulated allele pairs resolved correctly, including pairs as close as 3 repeats apart (where simple peak-finding fails).

## Installation

### From Source (Recommended for Development)

```bash
git clone https://github.com/berntpopp/open-pacmuci.git
cd open-pacmuci
make dev  # installs with uv in editable mode
```

### External Tool Dependencies

open-pacmuci requires the following bioinformatics tools on PATH:

| Tool | Version | Purpose | Install |
|------|---------|---------|---------|
| minimap2 | >= 2.28 | Long-read alignment | `conda install -c bioconda minimap2` |
| samtools | >= 1.21 | BAM processing | `conda install -c bioconda samtools` |
| bcftools | >= 1.10 | VCF filtering + consensus | `conda install -c bioconda bcftools` |
| Clair3 | >= 1.0.10 | Variant calling (HiFi model) | See [Clair3 docs](https://github.com/HKU-BAL/Clair3) |

### Conda Environment

```bash
conda env create -f conda/environment.yml
conda activate open-pacmuci-tools
pip install -e .
```

### Docker

```bash
docker pull ghcr.io/berntpopp/open-pacmuci:latest
docker run ghcr.io/berntpopp/open-pacmuci run --input reads.bam --output results/
```

### Important: Clair3 Python Environment

Clair3 requires its own Python environment (typically Python 3.9 with TensorFlow). When running open-pacmuci via `uv run`, the virtualenv Python may shadow the Clair3 conda Python. open-pacmuci handles this automatically by stripping `.venv/bin` from PATH when calling external tools. However, ensure the Clair3 conda environment's `bin/` directory is on your system PATH:

```bash
export PATH="/path/to/conda/envs/env_clair3/bin:$PATH"
```

## Usage

### Full Pipeline (Recommended)

```bash
open-pacmuci run \
  --input reads.fastq \
  --output results/ \
  --clair3-model /path/to/clair3/models/hifi \
  --threads 8
```

### Step-by-Step Execution

```bash
# 1. Generate reference ladder (once)
open-pacmuci ladder --output reference_ladder.fa

# 2. Map reads
open-pacmuci map --input reads.fastq --reference reference_ladder.fa --output-dir results/

# 3. Detect alleles
open-pacmuci alleles --input results/mapping.bam --output-dir results/

# 4. Call variants
open-pacmuci call --input results/mapping.bam --reference reference_ladder.fa \
  --alleles-json results/alleles.json --output-dir results/ \
  --clair3-model /path/to/models/hifi

# 5. Build consensus
open-pacmuci consensus --input results/mapping.bam --reference reference_ladder.fa \
  --alleles-json results/alleles.json --output-dir results/

# 6. Classify repeats
open-pacmuci classify --input results/consensus_allele_1.fa --output-dir results/
```

### Input Format

- **PacBio HiFi CCS reads** as FASTQ or BAM
- Reads should be from a PCR amplicon spanning the MUC1 VNTR (primers per Wenzel et al. 2018)
- Minimum recommended coverage: 10x per allele (higher coverage improves confidence)

### Output Files

| File | Description |
|------|-------------|
| `alleles.json` | Detected allele lengths, read counts, contig assignments |
| `allele_1/variants.vcf.gz` | Filtered Clair3 variants for allele 1 |
| `consensus_allele_1.fa` | Consensus FASTA (VNTR region only) |
| `repeats.json` | Per-repeat classification with confidence scores |
| `repeats.txt` | Human-readable VNTR structure string |
| `summary.json` | Combined allele + classification + mutation report |

## Generating Test Data with MucOneUp

[MucOneUp](https://github.com/berntpopp/muconeup) is a companion tool that generates simulated MUC1 VNTR haplotypes and PacBio HiFi amplicon reads with known ground truth, for benchmarking and validation.

### Prerequisites

- MucOneUp >= 0.44.0 (`pip install muc-one-up` or from [GitHub](https://github.com/berntpopp/muconeup))
- pbsim3 (via conda: `conda install -c bioconda pbsim3`)
- minimap2, samtools
- MucOneUp `config.json` (from the MucOneUp repository)

### Generate All Test Samples

```bash
# Ensure pbsim3 is available
export PATH="/path/to/conda/envs/env_pacbio/bin:$PATH"

# Generate 10 samples covering dupC, dupA, insG, insCCCC, del18_31,
# normal, homozygous, asymmetric, short, and long allele pairs
python scripts/generate_testdata.py
```

This creates simulated BAM files in `tests/data/generated/` with:

- 200x template coverage (PacBio HiFi amplicon)
- Realistic PCR bias (shorter alleles amplified more)
- Known VNTR structure and mutation positions as ground truth

### Run Integration Tests

```bash
export PATH="/path/to/conda/envs/env_pacbio/bin:$PATH"
uv run pytest tests/integration/ -v --no-cov
```

## Repeat Nomenclature

open-pacmuci uses the Vrbacka nomenclature (extended from Kirby et al. 2013 and Wenzel et al. 2018):

- **Pre-repeats:** `1`, `2`, `3`, `4`, `4'`, `5`, `5C` (precede the VNTR array)
- **Canonical repeats:** `X` (most common), `A`-`Z`, `aA`-`aZ` (variant types differing by 1-3 substitutions from X)
- **After-repeats:** `6`, `6'`, `7`, `8`, `9` (follow the VNTR array)
- **Mutations:** `X:59dupC`, `X:60dupA`, `X:58_59insG`, `C:42_57dupGGGCTCCACCGCCCCC`, etc.

A complete allele structure reads left-to-right as: `1 2 3 4 5 C F X X X ... X V 6 7 8 9`

## Known Mutations Detected

The pipeline includes a catalog of 13 known MUC1 frameshift mutations with pre-computed sequence templates for exact matching:

| Mutation | Repeat(s) | Effect | Citation |
|----------|-----------|--------|----------|
| 59dupC | X | +1bp (C duplication in 7C tract) | Kirby et al. 2013 |
| 60dupA | X | +1bp (A duplication) | Olinger et al. 2020 |
| 58_59insG | X | +1bp (G insertion) | Olinger et al. 2020 |
| 56_59dupCCCC | X | +4bp (CCCC duplication) | Vrbacka et al. 2025 |
| 18_31del | X | -14bp (deletion) | Vrbacka et al. 2025 |
| 42_57dupGGGCTCCACCGCCCCC | C | +16bp (16bp duplication) | Vrbacka et al. 2025 |
| 31ins25bp | A,B,J,K,N,S,X | +25bp (25bp insertion) | Saei et al. 2023 |
| delinsAT | X | net -1bp (delete 2, insert 2) | Olinger et al. 2020 |
| delGCCCA | multiple | -5bp (5bp deletion at start) | Saei et al. 2023 |
| insC_pos23 | A, E | +1bp | Vrbacka et al. 2025 |
| insG_pos58 | B, X | +1bp | Vrbacka et al. 2025 |
| insG_pos54 | B, J | +1bp | Vrbacka et al. 2025 |
| insA_pos54 | A, H | +1bp | Vrbacka et al. 2025 |

Novel mutations not in the catalog are still detected via edit-distance fallback, but with lower confidence.

## Validation Results

Tested against MucOneUp-simulated PacBio HiFi amplicon data:

| Metric | Result |
|--------|--------|
| VNTR classification (given correct sequence) | 20/20 haplotypes perfect (100% exact match) |
| Allele length detection (from BAMs, gap >= 5) | Exact match |
| Allele length detection (from BAMs, gap 3-4) | Within +/- 2 repeats |
| Full pipeline with Clair3 (dupC sample) | dupC detected at correct position, 100% confidence |
| Close allele pairs (gap 3-12 repeats) | 12/12 resolved by indel-valley splitting |

## License

MIT

## Citation

If you use open-pacmuci, please cite:

> Vrbacka A, Pristoupilova A, Kidd KO, et al. Long-Read Sequencing of the MUC1 VNTR: Genomic Variation, Mutational Landscape, and Its Impact on ADTKD Diagnosis and Progression. bioRxiv 2025. doi: [10.1101/2025.09.06.673538](https://doi.org/10.1101/2025.09.06.673538)

For the simulation toolkit:

> MucOneUp: MUC1 VNTR diploid reference simulator. [https://github.com/berntpopp/muconeup](https://github.com/berntpopp/muconeup)
