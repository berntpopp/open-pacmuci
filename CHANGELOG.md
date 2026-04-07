# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.3.0] - 2026-04-07

### Added
- Soft QUAL scoring with continuous confidence and boundary penalty

### Fixed
- Address review feedback on soft QUAL scoring
- Address Copilot review comments on soft QUAL scoring PR

## [0.2.0] - 2026-04-07

### Added
- MkDocs Material documentation site with GitHub Pages deploy
- Comprehensive README with deviations, usage, and validation

## [0.1.2] - 2026-04-07

### Added
- Indel-valley allele splitting for close allele pairs

## [0.1.1] - 2026-04-06

### Fixed
- Strip virtualenv PATH from subprocess calls for external tools
- Handle empty Clair3 VCFs and remove INFO/DP filter

## [0.1.0] - 2026-04-06

### Added
- Initial release of open-pacmuci pipeline
- Reference ladder generation (20-150 repeat units)
- Read mapping with minimap2
- Allele length detection with peak finding
- Variant calling with Clair3
- Consensus building with bcftools
- Repeat unit classification with Vrbacka nomenclature
- Full pipeline CLI with individual subcommands
- Docker image with conda-based tool stack
- Unit and integration test suite

## [0.0.0] - 2026-04-06

### Added
- Project scaffolding with uv, ruff, mypy, pytest, CI
- Initial pipeline implementation

[Unreleased]: https://github.com/berntpopp/open-pacmuci/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/berntpopp/open-pacmuci/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/berntpopp/open-pacmuci/compare/v0.1.2...v0.2.0
[0.1.2]: https://github.com/berntpopp/open-pacmuci/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/berntpopp/open-pacmuci/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/berntpopp/open-pacmuci/compare/v0.0.0...v0.1.0
[0.0.0]: https://github.com/berntpopp/open-pacmuci/releases/tag/v0.0.0
