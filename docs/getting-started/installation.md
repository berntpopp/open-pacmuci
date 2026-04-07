# Installation

open-pacmuci supports installation from source with external bioinformatics tool dependencies.

---

## Prerequisites

**System Requirements:**

- **Python:** 3.10, 3.11, or 3.12
- **RAM:** 8GB recommended (Clair3 variant calling is memory-intensive)
- **Operating System:** Linux (primary), macOS (partial -- Clair3 Linux only)

---

## Quick Install (From Source)

=== "uv (Recommended)"

    ```bash
    # Clone repository
    git clone https://github.com/berntpopp/open-pacmuci.git
    cd open-pacmuci

    # Install with uv in editable mode
    make dev

    # Verify installation
    open-pacmuci --version
    ```

=== "pip"

    ```bash
    # Clone repository
    git clone https://github.com/berntpopp/open-pacmuci.git
    cd open-pacmuci

    # Install in editable mode
    pip install -e .

    # Verify installation
    open-pacmuci --version
    ```

**What this installs:**

- open-pacmuci CLI tool
- Core Python dependencies (Click, PyYAML)
- Entry point for `open-pacmuci` command

**Does NOT include:** External bioinformatics tools (install separately -- see below).

---

## External Tool Dependencies

open-pacmuci requires the following bioinformatics tools on PATH:

| Tool | Version | Purpose | Install |
|------|---------|---------|---------|
| minimap2 | >= 2.28 | Long-read alignment to ladder | `conda install -c bioconda minimap2` |
| samtools | >= 1.21 | BAM processing and indexing | `conda install -c bioconda samtools` |
| bcftools | >= 1.10 | VCF filtering and consensus | `conda install -c bioconda bcftools` |
| Clair3 | >= 1.0.10 | Variant calling (HiFi model) | See [Clair3 docs](https://github.com/HKU-BAL/Clair3) |

!!! note "All tools required for full pipeline"
    The `open-pacmuci run` command requires all four tools. Individual subcommands (`ladder`, `classify`) can run without external tools.

---

## Conda Environment

Install all external tools in a single conda environment:

```bash
conda env create -f conda/environment.yml
conda activate open-pacmuci-tools
pip install -e .
```

---

## Docker

```bash
# Pull pre-built image
docker pull ghcr.io/berntpopp/open-pacmuci:latest

# Run full pipeline
docker run --rm \
  -v $(pwd)/data:/data \
  ghcr.io/berntpopp/open-pacmuci:latest \
  run --input /data/reads.bam --output /data/results/
```

---

## Clair3 Python Environment

!!! warning "Clair3 requires its own Python environment"
    Clair3 typically runs under Python 3.9 with TensorFlow. When running open-pacmuci via `uv run`, the virtualenv Python may shadow the Clair3 conda Python.

open-pacmuci handles this automatically by stripping `.venv/bin` from PATH when calling external tools. However, ensure the Clair3 conda environment's `bin/` directory is on your system PATH:

```bash
export PATH="/path/to/conda/envs/env_clair3/bin:$PATH"
```

---

## Verify Installation

### Test Core Functionality

```bash
# Check version
open-pacmuci --version

# View available commands
open-pacmuci --help

# Generate reference ladder (no external tools needed)
open-pacmuci ladder --output test_ladder.fa
```

### Test Full Pipeline

```bash
# Ensure external tools are available
which minimap2 samtools bcftools run_clair3.sh

# Run pipeline (requires all tools + Clair3 model)
open-pacmuci run \
  --input reads.fastq \
  --output results/ \
  --clair3-model /path/to/clair3/models/hifi
```

---

## Development Setup

For contributors and developers:

```bash
# Install development environment
make dev

# Run all quality checks
make check

# Individual commands
make test          # Run pytest with coverage
make test-fast     # Unit tests only, no coverage
make lint          # Check code quality (ruff)
make format        # Auto-format code
make type-check    # Run mypy
```

**Development Commands:**

| Command | Action |
|---------|--------|
| `make dev` | Install with dev dependencies (editable mode) |
| `make test` | Run all tests with coverage |
| `make test-fast` | Unit tests only, no coverage |
| `make lint` | Check code quality (ruff) |
| `make format` | Auto-format code |
| `make type-check` | Run mypy type checker |
| `make check` | All quality checks (test + lint + type) |

---

## Troubleshooting

**Issue:** `open-pacmuci: command not found`

**Solution:** Ensure installation completed and the package is on PATH:

```bash
pip install -e .
# or if using uv:
uv pip install -e .
```

---

**Issue:** `run_clair3.sh: command not found`

**Solution:** Activate the Clair3 conda environment and ensure it is on PATH:

```bash
export PATH="/path/to/conda/envs/env_clair3/bin:$PATH"
which run_clair3.sh
```

---

**Issue:** `minimap2` or `samtools` not found

**Solution:** Install via conda:

```bash
conda install -c bioconda minimap2 samtools bcftools
```

---

## Next Steps

- **[Quick Start](quickstart.md)** -- Run your first analysis
- **[Core Concepts](concepts.md)** -- Understand the pipeline architecture
- **[Deviations from PacMUCI](deviations.md)** -- What changed and why
