# charlesfoster/h2seq: Output

## Overview

`h2seq` writes per-sample outputs under the top-level results directory and also creates several run-level summary folders. The exact set of files depends on:

- long-read versus short-read input
- whether reference selection is enabled
- whether primer trimming is enabled
- whether `--run_hcv_glue` is enabled
- whether a sample has enough data to proceed through each stage

## Per-sample Structure

Each sample is organised by read type:

- `<sample>/long_reads/`
- `<sample>/short_reads/`

Within those directories you may see some or all of the following subdirectories.

### `qc_metrics`

Typical contents:

- raw and post-QC `seqkit stats` tables
- `FastQC` reports for short reads
- read-stat summaries used during reference selection
- `NANOQ` statistics for long reads

### `clean_reads`

Typical contents:

- trimmed and filtered reads, if saving cleaned reads is enabled
- `fastp` JSON and HTML reports for short reads
- filtered long-read FASTQ files, if saving cleaned reads is enabled

### `reference_selection`

Present when automatic reference selection is used, and also populated with synthetic reference metadata when `--skip_reference_selection` is set.

Typical contents:

- competitive BAM/BAI and `*.reference_selection_ranking.tsv` files from default `minimap2` reference selection
- abundance estimates when `--reference_selection_tool salmon` or `--reference_selection_tool kallisto` is used
- `*.best_reference.tsv`
- `*.best_reference.txt`
- selected reference FASTA files

### `read_mapping`

Contains read-mapping intermediates where publishing is enabled for the relevant aligner.

### `primer_clipping`

Present when primer trimming is used. Contains intermediate files related to primer localisation and clipping.

### `coverage`

Typical contents:

- `*.coverage_summary.tsv`
- `*.hcv_glue_coverage.tsv` parsed from the HCV-GLUE HTML report for HCV runs with `--run_hcv_glue`
- `*.mapped_read_count.txt`

### `consensus`

Typical contents:

- split consensus FASTA files such as `*.consensus_main.fa`
- alternate consensus FASTA files when multiple candidate references are retained
- mapped BAM/CSI files when `--save_mapped_bam` is enabled

### `variant_calling`

Contains intermediate and final files from the variant-calling and consensus workflow.

If `--majority_allele_consensus` is enabled, the final consensus contains no IUPAC ambiguity codes. In that mode, only majority variants with allele frequency `>= 0.5` are incorporated into the consensus, regardless of lower user-supplied SNV or indel thresholds.

### `reporting`

Present when optional reporting steps are enabled, most notably `HCV-GLUE`.

## Run-level Outputs

### `multiqc/`

Contains the aggregated MultiQC report. The report now includes native modules for supported tool outputs such as:

- `FastQC`
- `fastp`
- `seqkit stats`

It also includes an `h2seq Run Summary` custom section built from the final run summary table.

### `run_summary/`

Contains:

- `run_summary.csv`
- `run_summary_mqc.json`

`run_summary.csv` is the main per-sample summary table. Depending on the run, columns can include:

- designated genotype and subtype
- chosen mapping reference
- genome coverage and mean depth
- HCV polyprotein and region coverage
- raw read counts and post-QC read counts
- mapped read counts
- consensus count and main consensus path
- HCV-GLUE report path
- notes describing skipped or failed sample-level outcomes

Rows are emitted even for samples that stop early because they have no usable reads.

### `all_consensus_genomes/`

Collects consensus FASTA outputs from the entire run into one folder. Filenames are prefixed with `long.` or `short.` to reduce collisions.

### `all_hcv_glue_reports/`

When `--run_hcv_glue` is enabled, all HTML reports are collected into one folder. Filenames are prefixed with `long.` or `short.`.

### `pipeline_info/`

Contains standard Nextflow and nf-core run metadata, such as:

- execution reports and traces
- validated samplesheet
- parameter snapshot
- software version summary used by MultiQC

## Notes On Missing Outputs

Not every sample will produce every file. Negative controls and failed libraries can legitimately stop at QC, mapping, or consensus generation. In those cases:

- the pipeline continues
- downstream sample outputs may be absent
- the reason should be visible in `combined_results_summary.csv`
