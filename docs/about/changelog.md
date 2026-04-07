# Changelog

All notable changes to open-pacmuci will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

---

## [0.2.0] - 2026-04-06

### Added
- Indel-valley allele splitting for close allele pairs (3-9 repeats apart)
- Resolves allele pairs that gap-based clustering merges into a single peak
- Analyzes CIGAR indel lengths to find two local minima corresponding to true allele lengths

### Changed
- Allele detection now uses indel-valley splitting as fallback when gap-based clustering finds only one peak

---

## [0.1.2] - 2026-04-06

### Fixed
- Strip virtualenv PATH from subprocess calls for external tools (Clair3 compatibility)

---

## [0.1.1] - 2026-04-06

### Fixed
- Handle empty Clair3 VCFs and remove INFO/DP filter (fixes #8)

---

## [0.1.0] - 2026-04-05

### Added
- Initial pipeline implementation with 5 stages
- Reference ladder generation (150 contigs, 1-150 repeat units)
- Read mapping with minimap2
- Allele detection with gap-based clustering
- Variant calling with Clair3 (per-allele)
- Consensus construction with bcftools and anchor-based flanking trim
- Repeat classification using Vrbacka nomenclature
- Pre-computed mutation template catalog (13 known mutations)
- Per-repeat and per-allele confidence scoring
- VCF-backed mutation validation with confidence adjustment
- Same-length allele disambiguation using Clair3 het genotypes
- Bidirectional classification fallback for novel large mutations
- Full pipeline `run` command and individual stage commands
- Click CLI with subcommands: `ladder`, `map`, `alleles`, `call`, `consensus`, `classify`, `run`

---

[Unreleased]: https://github.com/berntpopp/open-pacmuci/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/berntpopp/open-pacmuci/releases/tag/v0.2.0
[0.1.2]: https://github.com/berntpopp/open-pacmuci/releases/tag/v0.1.2
[0.1.1]: https://github.com/berntpopp/open-pacmuci/releases/tag/v0.1.1
[0.1.0]: https://github.com/berntpopp/open-pacmuci/releases/tag/v0.1.0
