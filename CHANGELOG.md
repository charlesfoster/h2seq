# charlesfoster/h2seq: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.6.0 - [Unreleased]

### `Added`

- Optional duplicate marking (`--mark_duplicates`) for Illumina alignments using `samtools fixmate`/`markdup`, with duplicate statistics published to `qc_metrics`
- Separate thresholds for writing IUPAC ambiguity codes into the IUPAC consensus: `--illumina_iupac_min_af`, `--ont_iupac_min_af`, `--illumina_iupac_min_alt_reads` and `--ont_iupac_min_alt_reads`
- `*.final_consensus.vcf.gz` output, the subset of `*.final.vcf.gz` applied to the IUPAC consensus
- `twist` profile for Twist hybrid-capture Illumina data (skips primer trimming and marks duplicates)

### `Changed`

- Updated Clair3 from 1.1.1 to 2.0.3 (Docker, Singularity and conda) and changed the default `--clair3_model_path` to `/opt/models/r1041_e82_400bps_sup_v520`. Custom Clair3 models must now be v2 (PyTorch) model directories containing `pileup.pt` and `full_alignment.pt`
- Low-frequency mixed Illumina SNVs (below `--illumina_iupac_min_af` or supported by fewer than `--illumina_iupac_min_alt_reads` reads) are still reported in `*.final.vcf.gz` but no longer produce IUPAC codes in the consensus

### `Fixed`

- Disabled the catch-all `publishDir` default so processes without their own publish rule no longer create stray top-level directories in `--outdir` (e.g. `compress/` containing the prepared LoFreq VCF)
- Run summary and MultiQC section builders now receive their inputs as staged files instead of reading the published output directory, which was invisible inside containers when `--outdir` was outside the project directory, producing empty summaries

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
