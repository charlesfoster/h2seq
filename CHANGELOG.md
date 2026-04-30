# charlesfoster/h2seq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.3.0 - [2026-04-30]

### `Added`

- Added minimap2 competitive reference selection as the default reference-selection method
- Added sample-level QC failure reasons for missing primary mapped reads, low selected-reference coverage, and consensus genomes that are insufficient for HCV-GLUE reporting
- Added run-summary reporting to mark genotype/subtype assignments as unreliable for QC-failed samples

### `Changed`

- Raised the minimum supported Nextflow version to `26.04.0`

## v0.2.0 - [2026-03-24]

### `Added`

- Support for `fastp` QC-only operation via `discard_trimmed_pass`, allowing reports to be generated without writing passing trimmed reads
- An explicit consensus-generation path based on variant calling, variant filtering/masking, and `bcftools consensus`, rather than relying only on `samtools consensus`
- Stub-run compatibility improvements for FAIDX- and LoFreq-dependent paths so the Illumina test profile completes successfully under `-stub-run`
- Module patch tracking entries for additional nf-core components in `modules.json`

### `Fixed`

- Reworked BWA index channel handling in short-read mapping and primer-mapping subworkflows to match the current `BWA_INDEX` output contract
- Updated the short-read workflow call to pass the full `FASTP` input signature explicitly
- Improved the consensus workflow wiring around reference indexing, low-depth masking, and downstream variant-to-consensus steps
- Corrected MultiQC section ordering so the custom H2Seq summary renders before the software versions section

### `Dependencies`

- Refreshed local nf-core module patch metadata for `bwa/index`, `fastp`, `fastqc`, `kallisto/quant`, and related patched components tracked in `modules.json`

### `Deprecated`

## v0.1.1 - [2025-07-29]

Initial release of charlesfoster/h2seq, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Test dataset (and associated CI actions) for short read (Illumina) data

### `Fixed`

- Ensured the 'master' branch is the default

### `Dependencies`

### `Deprecated`

## v0.1.0 - [2025-07-28]

Initial release of charlesfoster/h2seq, created with the [nf-core](https://nf-co.re/) template.

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`
