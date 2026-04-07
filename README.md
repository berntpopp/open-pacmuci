# open-pacmuci

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Tests](https://github.com/berntpopp/open-pacmuci/workflows/Test%20%26%20Quality/badge.svg)](https://github.com/berntpopp/open-pacmuci/actions)
[![Documentation](https://img.shields.io/badge/docs-MkDocs%20Material-blue)](https://berntpopp.github.io/open-pacmuci/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Open-source MUC1 VNTR analysis pipeline for PacBio HiFi amplicon data.

---

## Overview

open-pacmuci reconstructs the **PacMUCI bioinformatics pipeline** (Vrbacka et al. 2025) for analyzing PacBio HiFi long-read sequencing data of the MUC1 VNTR region. The original source code was never published; this project rebuilds it from the methods and extends it with algorithmic improvements.

**Key Capabilities:**

- Detect frameshift mutations in the MUC1 VNTR that cause ADTKD-MUC1 kidney disease
- Resolve close allele pairs (3-9 repeats apart) via indel-valley splitting
- Pre-computed mutation template catalog (13 known mutations) with O(1) lookup
- Per-repeat and per-allele confidence scoring with VCF cross-validation
- Fully automated -- no manual IGV inspection steps

---

## Quick Start

```bash
open-pacmuci run \
  --input reads.fastq \
  --output-dir results/ \
  --clair3-model /path/to/clair3/models/hifi \
  --threads 8
```

**Documentation:** https://berntpopp.github.io/open-pacmuci/

---

## Installation

```bash
git clone https://github.com/berntpopp/open-pacmuci.git
cd open-pacmuci
make dev
open-pacmuci --version
```

**Requires:** minimap2, samtools, bcftools, Clair3 on PATH. See [Installation Guide](https://berntpopp.github.io/open-pacmuci/getting-started/installation/).

---

## Citation

If you use open-pacmuci, please cite:

> Vrbacka A, Pristoupilova A, Kidd KO, et al. Long-Read Sequencing of the MUC1 VNTR. *bioRxiv.* 2025. doi: [10.1101/2025.09.06.673538](https://doi.org/10.1101/2025.09.06.673538)

---

## License

MIT License -- see [LICENSE](LICENSE) for details.

---

**Maintained by:** [Bernt Popp](https://github.com/berntpopp)
