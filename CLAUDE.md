# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Open-source reconstruction of the PacMUCI pipeline (Vrbacka et al. 2025, bioRxiv: 10.1101/2025.09.06.673538) for analyzing PacBio HiFi sequencing data of PCR amplicons spanning the MUC1 VNTR region. The original source code was never published; this project rebuilds it from the methods description.

## Pipeline Architecture

The pipeline has 7 stages executed sequentially:

1. **Reference Ladder Generation** -- Build synthetic FASTA with 131 contigs (20-150 repeat units each), structured as `[flanking] [pre-repeats 1-5] [N * canonical 60bp repeat] [after-repeats 6-9] [flanking]`
2. **Read Mapping** -- bwa-mem alignment of HiFi reads to the ladder reference
3. **Allele Length Determination** -- samtools idxstats counts per contig, peak detection to find the two allele lengths (min 10x coverage)
4. **Read Filtering** -- Remove chimeric/incomplete reads
5. **Variant Calling** -- Clair3 with PacBio HiFi model to detect frameshift mutations
6. **Consensus Construction** -- bcftools consensus to build per-allele FASTA
7. **Repeat Unit Classification** -- Split into 60bp units, classify using Vrbacka nomenclature (leveraging [muc1repeats](https://github.com/pristanna/muc1repeats))

## Tech Stack

- **Language:** Python 3.10+
- **CLI:** Click 8.0+
- **External tools:** bwa-mem, samtools, Clair3, bcftools
- **Packaging:** pyproject.toml, Docker/Singularity, Conda
- **Entry point:** `open-pacmuci run --input reads.bam --output results/`

## Repository Structure

```
src/open_pacmuci/     # Main package
  cli.py              # Click CLI entry point
  ladder.py           # Reference ladder generation
  mapping.py          # bwa-mem mapping + allele detection
  alleles.py          # Allele length determination (peak detection)
  calling.py          # Clair3 variant calling
  consensus.py        # bcftools consensus
  classify.py         # Repeat unit classification (Vrbacka nomenclature)
  tools.py            # External tool wrappers
  config.py           # Configuration management
  version.py          # Package version
data/repeats/         # Repeat type definitions (repeats.json)
data/reference/       # Generated ladder reference (gitignored)
tests/unit/           # Unit tests (no external tools required)
tests/integration/    # Integration tests (require bwa, samtools, Clair3, bcftools)
docker/               # Dockerfile
conda/                # environment.yml
scripts/              # Helper scripts (e.g., generate_testdata.sh)
```

## Development Commands

```bash
make dev               # Install with dev dependencies
make test              # Run all tests with coverage
make test-fast         # Unit tests only, no coverage
make lint              # Run ruff linter
make format            # Format code
make type-check        # Run mypy
make check             # All quality checks
make generate-testdata # Generate MucOneUp test data
```

Run a single test:

```bash
uv run pytest tests/unit/test_classify.py::TestClassifyRepeat::test_exact_match -v --no-cov
```

## Testing Strategy

Test data is generated using [MucOneUp](https://github.com/berntpopp/muconeup) (v0.43.2+) to simulate PacBio HiFi amplicon reads with known ground truth. See `.planning/TESTING_WITH_MUCONEUP.md` for detailed test data generation commands and validation approach. Key test sets cover: basic validation, multiple mutation types (13 known mutations including dupC, dupA, insG, etc.), coverage sensitivity, extreme VNTR lengths, and cross-platform (ONT) validation.

## Key Domain Concepts

- **VNTR:** Variable Number Tandem Repeat -- the MUC1 gene contains a region of 60bp repeat units (typically 20-150 copies per allele)
- **Alleles:** Each person has two alleles (maternal/paternal) that may differ in repeat count
- **Frameshift mutations:** Insertions/deletions not divisible by 3 within repeat units (e.g., dupC) cause disease (ADTKD-MUC1)
- **Vrbacka nomenclature:** Classification system for repeat unit types: pre-repeats (1-5), canonical (A-Z, aA-aZ...), after-repeats (6-9)
- **PCR bias:** Shorter alleles amplify more efficiently, producing unequal read counts between alleles
