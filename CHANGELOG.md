# Changelog

All notable changes to **riboseq-flow** will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

---

## [v1.2.0] - 2025-11-23
### Added
- Added CITATION.cff for GitHub “Cite this repository” support.
- Added Zenodo concept DOI badge and citation badge to README.
- Added Nextflow version badge (≥ 21.10.3).
- Added CODEOWNERS file for maintainership.
- Added NOTICE file for attribution.
- Added CHANGELOG.md.
- Added additional Nextflow versions tested in documentation.
- Added `base.config` for resource management.

### Changed
- Updated workflow metadata, licensing to GPL-3.0.
- Added MultiQC configuration file as an input to the `MULTIQC` process, closes [#98](https://github.com/iraiosub/riboseq-flow/issues/98).
- Cap resources for profile `test` testing and override base profile labels.
- Added one more sample in the samplesheet for for profile `test` testing.
- Change references path in CI.

### Fixed
- Removed `tag` directive from processes not requiring it, enabling compatibility with Nextflow `25.10.0`, closes [#99](https://github.com/iraiosub/riboseq-flow/issues/99).
- Bug in `PREPARE_RIBOSEQ_REFERENCE`, `params.gtf` now replaced with `genome_gtf`.

---

## [v1.1.1] - 2024-01-23

### Fixed
- Fixed PCA `if` statement and colors.

---

## [v1.1.0] - 2024-01-18

### Added
- Add option to count reads over CDS.

### Changed
- Updated PCA display in MultiQC.
- Bump STAR, samtools, MultiQC versions.
- Made saving of trimmed files and saving indexes for alignment optional.
- Tidy-up `CUTADAPT` module.
- Removed cutadapt TSO arguments.
- Use all adaptor-trimmed reads for mapping length analysis and useful reads calculations (instead of length-filtered).
- Renamed contaminants input param name.
- Documentation updates.


---

## [v1.0.2] - 2023-12-19

### Added
- Add read length filter qc to MultiQC.

### Changed
- Changed order of sections in MultiQC.
- Updated location of coverage tracks in docs.
- Clarifications in docs and MultiQC section descriptions.
- PCA performed only if multiple samples provided.

---

## [v1.0.1] - 2023-12-13

### Fixed
- Fixed P-site BED output.

### Changed
- Updated run instructions in documentation.

---
## [v1.0.0] - 2023-12-13

- The first main release of **riboseq-flow**.


[v1.2.0]: https://github.com/iraiosub/riboseq-flow/releases/tag/v1.2.0
[v1.1.1]: https://github.com/iraiosub/riboseq-flow/releases/tag/v1.1.1
[v1.1.0]: https://github.com/iraiosub/riboseq-flow/releases/tag/v1.1.0
[v1.0.2]: https://github.com/iraiosub/riboseq-flow/releases/tag/v1.0.2
[v1.0.1]: https://github.com/iraiosub/riboseq-flow/releases/tag/v1.0.1
[v1.0.0]: https://github.com/iraiosub/riboseq-flow/releases/tag/v1.0.0
