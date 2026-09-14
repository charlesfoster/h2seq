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
- `*.mixed_assignment.tsv`, containing uniquely assigned, ambiguous, and unassigned fragment counts after competitive mapping
- `*.hcv_glue_coverage.tsv` parsed from the HCV-GLUE HTML report for HCV runs with `--run_hcv_glue`
- `*.mapped_read_count.txt`

For samples with more than one retained reference, coverage is reported separately
for every reference. Fragments are assigned using the primary alignments of the read
or read pair. Fragments below `--mixed_assignment_min_mapq`, or pairs whose primary
alignments point to different references, are counted as ambiguous and do not
contribute to genotype-specific coverage, variant calling, or consensus generation.

### `consensus`

Typical contents:

- split IUPAC consensus FASTA files such as `*.consensus_main.fa`
- split simple majority-allele consensus FASTA files such as `*.simple.consensus_main.fa`
- alternate consensus FASTA files when multiple candidate references are retained
- mapped BAM/CSI files when `--save_mapped_bam` is enabled

### `variant_calling`

Contains intermediate and final files from the variant-calling and consensus workflow.

The workflow always builds both an IUPAC ambiguity consensus and a simple majority-allele consensus. The simple consensus incorporates variants with allele frequency `>= 0.5` and is used for downstream reporting steps such as HCV-GLUE.

### `reporting`

Present when optional reporting steps are enabled, most notably `HCV-GLUE`.

## Run-level Outputs

### `multiqc/`

Contains the aggregated MultiQC report. The report now includes native modules for supported tool outputs such as:

- `FastQC`
- `fastp`
- `seqkit stats`

It also includes an `h2seq Run Summary` custom section built from the final run summary table.

### Run summary tables

The top-level results directory contains:

- `combined_results_summary.csv`
- `reference_component_summary.csv`

`combined_results_summary.csv` has one row per sample and read type. For a mixed
infection, its coverage and mean depth describe the selected main reference only;
`mixed_infection`, `secondary_subtypes`, and `component_fractions` describe the
final competitively assigned components. The preliminary selection-stage mixture
score remains available only in `*.best_reference.tsv` for diagnostic purposes.
Depending on the run, columns can include:

- designated genotype and subtype
- chosen mapping reference
- genome coverage and mean depth
- HCV polyprotein and region coverage
- raw read counts and post-QC read counts
- mapped read counts
- consensus count and main consensus path
- HCV-GLUE report path
- notes describing skipped or failed sample-level outcomes

Rows are emitted even for samples that stop early because they have no usable reads. Samples are marked `qc_fail` when they have no primary mapped reads, when less than `--qc_min_ref_coverage_pct` of the selected reference is covered at `--consensus_min_depth`, or when the main consensus is all `N` and cannot be sent to HCV-GLUE. Genotype/subtype calls for QC-failed samples are blanked in the combined summary and flagged as unreliable.

`reference_component_summary.csv` has one row per retained reference component.
It identifies each component as `main` or `secondary` and reports its own genome
coverage, mean depth, HCV polyprotein coverage where available, uniquely assigned
fragment count and fraction of confidently assigned fragments, together with the
sample-level ambiguous and unassigned fragment counts. This is the table to use when
comparing the coverage of genotype 1a and genotype 3a in the same sample.

Mixed infections produce a separate component PDF for the main and every retained
secondary reference. Filenames include the reference and component role. Each PDF
contains component-specific depth, coverage, assignment statistics, consensus
interpretation, and HCV-GLUE feature coverage where available.

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
