# Changelog

All notable changes to **riboseq-flow** will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/).

---
## [v1.4.0] - 2026-11-25
### Added
- Added `--psite_length_range` to select P-site read lengths independently of the expected length range used for read QC.
- Added `--psite_per_sample` to run P-site inference separately for each sample and merge offset and coverage tables. Individual tables and compact diagnostics are published in sample-named folders, with merged tables and per-read P-site tables at the root of `psites`.
- Added `--psite_transcripts all_cds` to analyse all transcripts observed in the transcriptome BAM that have an annotated CDS. Representative transcript selection remains the default; RUST, P-site PCA and genomic P-site tracks use the matching reference.
- Added strand-specific genomic P-site BigWig or bedGraph tracks with `--genomic_psite_coverage`, retaining `--genomic_codon_coverage` as an alias. Repeated isoform assignments of a read count once; reads with conflicting genomic P-sites are excluded and reported in mapping statistics.
- Added a riboWaltz P-site QC section to MultiQC with region and frame distributions, start and stop codon profiles, retained read lengths and transcript selection details.
- Added compact RUST-inspired QC to MultiQC, with a heatmap of supported sample/read-length combinations and profiles for the most frequent and most periodic read lengths. Displayed lengths require at least 1,000 CDS P-sites after subsampling and dominant-frame selection; full PDF and TSV outputs retain all analysed lengths.
- Published the underlying RUST plot data and read-length support tables in `riboseq_qc/rust_analysis`.
- Added QC for the most periodic supported expected read length, a selected-length summary table, and warnings when the most abundant useful read length falls outside the expected range.
- Added a read 5' end profile around the stop codon to MultiQC, alongside the start codon profile.
- Added software versions grouped by process to MultiQC and a single consolidated `pipeline_info/software_versions.yml`, including R/Python dependencies and pipeline/Nextflow metadata from the current run.

### Changed
- Limited published genomic P-site outputs to coverage tracks and read mapping-statistics tables; BED files remain in the Nextflow work directory.
- Expanded the combined QC PDF to include dominant-frame percentages and dominant / least abundant frame ratios for the expected range, most abundant length and most periodic expected length.
- Adjusted CPU and memory requests for preprocessing, QC, reporting and track generation.

### Fixed
- Released genomic P-site mapping intermediates sooner and increased all-CDS track conversion memory to 72 GB with retry scaling, retaining representative-mode resource settings and read-deduplication rules.
- Reduced RUST memory use by reading only required P-site columns and joining transcript sequences after subsampling; increased its default resource allocation to 36 GB with retry scaling.
- Quoted the P-site inference script path so it can be launched from repository directories containing spaces.
- Fixed merging of gzipped per-sample P-site tables while preserving sample/transcript identifiers and coverage counts.
- Quoted UMI barcode patterns, including embedded apostrophes, so regex extraction no longer fails during shell parsing.
- Quoted Cutadapt adapter arguments so per-adapter settings containing semicolons reach Cutadapt intact.

---
## [v1.3.0] - 2026-11-23
### Added
- Added frame distributions and dominant-frame percentages for the most abundant useful read length.
- Added dominant / least abundant frame ratios for the expected read length range and the most abundant useful read length.
- Exported frame counts by read length and distances from start and stop codons in `riboseq_qc/qc_tables`, and combined QC metrics in `qc_summary.tsv.gz`.

### Changed
- Updated project licensing to MIT.
- Clarified the MultiQC frame plot title to identify the expected read length range.

### Fixed
- Filled missing UTR counts with zero so the MultiQC Region section remains a stacked bar chart, defaulting to percentages that sum to 100% per sample.

---
## [v1.2.3] - 2026-05-25
### Changed
- Pinned custom Docker container images to immutable SHA digests for reproducible execution.
- Updated Singularity pre-download instructions to use digest-pinned container references.
- Refined README wording for the "Why use riboseq-flow" section.

### Fixed
- Fixed `--skip_premap` runs so `--contaminants_fasta` is not required or checked when premapping is disabled.
- Replaced non-idiomatic bitwise `&` with logical `&&` in the `RUST_QC` skip condition.
- Pinned the GitHub Actions CI Nextflow version to `25.10.0`.

---
## [v1.2.2] - 2026-03-28
### Added
- Added .devcontainer for GitHub Codespaces support.

---

## [v1.2.1] - 2025-12-15
### Fixed
- Added `first()` to ensure the GTF channel can be reused for multiple samples for gene-level quantification with featureCounts (closes #89).

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
- Updated workflow metadata and licensing.
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

[v1.4.0]: https://github.com/iraiosub/riboseq-flow/releases/tag/v1.4.0
[v1.3.0]: https://github.com/iraiosub/riboseq-flow/releases/tag/v1.3.0
[v1.2.3]: https://github.com/iraiosub/riboseq-flow/releases/tag/v1.2.3
[v1.2.2]: https://github.com/iraiosub/riboseq-flow/releases/tag/v1.2.2
[v1.2.1]: https://github.com/iraiosub/riboseq-flow/releases/tag/v1.2.1
[v1.2.0]: https://github.com/iraiosub/riboseq-flow/releases/tag/v1.2.0
[v1.1.1]: https://github.com/iraiosub/riboseq-flow/releases/tag/v1.1.1
[v1.1.0]: https://github.com/iraiosub/riboseq-flow/releases/tag/v1.1.0
[v1.0.2]: https://github.com/iraiosub/riboseq-flow/releases/tag/v1.0.2
[v1.0.1]: https://github.com/iraiosub/riboseq-flow/releases/tag/v1.0.1
[v1.0.0]: https://github.com/iraiosub/riboseq-flow/releases/tag/v1.0.0
