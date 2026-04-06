# open-pacmuci

Open-source reconstruction of the PacMUCI pipeline for MUC1 VNTR analysis from long-read sequencing data.

## Background

PacMUCI was described in [Vrbacka et al. 2025](https://doi.org/10.1101/2025.09.06.673538) as a pipeline for analyzing PacBio SMRT sequencing data of PCR amplicons spanning the MUC1 VNTR. The original source code was not published. This project reconstructs the pipeline from the published methods as an open-source, containerized tool.

## What it does

1. **Allele length determination** -- Maps long reads to a synthetic reference ladder (contigs with 20-150 repeat units) and identifies the two allele lengths from mapping peaks
2. **Variant calling** -- Uses Clair3 to detect frameshift mutations within the VNTR
3. **Repeat classification** -- Classifies each 60bp repeat unit using the Vrbacka nomenclature system

## Status

Under development. See `.planning/PLAN.md` for implementation details.

## License

MIT
