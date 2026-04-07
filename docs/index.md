# open-pacmuci

**Open-source MUC1 VNTR analysis pipeline for PacBio HiFi amplicon data**

---

## What is open-pacmuci?

open-pacmuci is an open-source reconstruction of the PacMUCI bioinformatics pipeline (Vrbacka et al. 2025) for analyzing PacBio HiFi long-read sequencing data of the MUC1 VNTR region. The original pipeline source code was never released; this project rebuilds it from the published methods and extends it with several algorithmic improvements including **indel-valley allele splitting**, **mutation template matching**, and **per-repeat confidence scoring**.

### Why open-pacmuci?

**Detect frameshift mutations** in the MUC1 VNTR that cause ADTKD-MUC1 kidney disease
**Resolve allele pairs** even when they differ by only 3-9 repeat units (indel-valley splitting)
**Score confidence** per repeat unit and per allele with VCF cross-validation
**Classify repeats** using the Vrbacka nomenclature with a pre-computed mutation catalog
**Fully open source** -- MIT licensed, reproducible, no manual inspection steps

---

## Scientific Context

**Autosomal Dominant Tubulointerstitial Kidney Disease caused by MUC1 mutations (ADTKD-MUC1)** is a rare genetic kidney disease caused by frameshift mutations in the Variable Number Tandem Repeat (VNTR) region of the MUC1 gene. The VNTR consists of 20-125 copies of a degenerate 60-bp repeat unit with extremely high GC content (>80%), making it inaccessible to standard short-read sequencing.

The most common mutation, **59dupC**, duplicates a cytosine in the heptanucleotide C-tract of a canonical X repeat unit, producing a +1 frameshift that leads to expression of the toxic MUC1fs protein. **Long-read sequencing** (PacBio SMRT) combined with specialized bioinformatics is required to resolve the full VNTR structure and locate pathogenic frameshift mutations.

---

## Key Features

### Mutation Catalog

Pre-computed sequence templates for **13 known MUC1 frameshift mutations** enable O(1) lookup for common mutations. Novel mutations are detected via edit-distance fallback.

### Indel-Valley Allele Splitting

Resolves close allele pairs (3-9 repeats apart) that simple peak-finding merges into a single cluster. Analyzes CIGAR indel lengths to find two local minima corresponding to the true allele lengths.

### Confidence Scoring

Per-repeat and per-allele confidence scores with **VCF cross-validation** -- flags low-confidence calls for automated QC.

### VCF Validation

Clair3 variant calls are cross-referenced against repeat classifications to adjust confidence scores and detect compound heterozygotes.

---

## Quick Example

```bash
open-pacmuci run \
  --input reads.fastq \
  --output-dir results/ \
  --clair3-model /path/to/clair3/models/hifi \
  --threads 8
```

**Output:** Per-allele VNTR structure, mutation calls with exact repeat position, and confidence scores.

---

## Documentation

<div class="grid cards" markdown>

-  **[Getting Started](getting-started/installation.md)**
  Installation, quick start tutorial, and core concepts

-  **[Guides](guides/benchmarking.md)**
  Benchmarking with MucOneUp simulated data

-  **[Reference](reference/cli.md)**
  CLI commands, known mutations, and repeat nomenclature

-  **[About](about/citation.md)**
  Citation guide, license, and changelog

</div>

---

## Citation

If you use open-pacmuci, please cite:

```bibtex
@article{vrbacka2025pacmuci,
  author = {Vrbacka, A. and Pristoupilova, A. and Kidd, K.O. and others},
  title = {Long-Read Sequencing of the MUC1 VNTR: Genomic Variation,
           Mutational Landscape, and Its Impact on ADTKD Diagnosis
           and Progression},
  journal = {bioRxiv},
  year = {2025},
  doi = {10.1101/2025.09.06.673538}
}
```

See [Citation Guide](about/citation.md)

---

**Development Status:** Active | **License:** MIT | **Maintained by:** [Bernt Popp](https://github.com/berntpopp)
