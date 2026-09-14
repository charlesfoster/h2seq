# charlesfoster/h2seq: Usage

> Documentation for pipeline parameters is generated automatically from `nextflow_schema.json`. This page focuses on the parts of usage that are specific to `h2seq`.

## Overview

`h2seq` accepts ONT long reads, Illumina short reads, or a mixture of both in the same run. The pipeline can either:

- choose a reference automatically from a multifasta or preset reference set
- skip reference selection and map directly to a user-supplied reference

The `hcv` preset adds HCV-oriented defaults and enables optional HCV-specific reporting such as `HCV-GLUE`.

## Samplesheet

The input samplesheet is a CSV with these columns:

| Column          | Required | Description                                                                               |
| --------------- | -------- | ----------------------------------------------------------------------------------------- |
| `sample`        | yes      | Sample identifier. Spaces are not allowed.                                                |
| `long_reads`    | no       | Path to an ONT FASTQ or FASTQ.GZ file.                                                    |
| `short_reads_1` | no       | Path to Illumina read 1 FASTQ.GZ file.                                                    |
| `short_reads_2` | no       | Path to Illumina read 2 FASTQ.GZ file. Required when paired-end short reads are provided. |

Each row must contain at least `long_reads` or `short_reads_1`.

Example mixed-modality samplesheet:

```csv
sample,long_reads,short_reads_1,short_reads_2
sample_long_only,/data/sample_long.fastq.gz,,
sample_short_only,,/data/sample_short_R1.fastq.gz,/data/sample_short_R2.fastq.gz
sample_both,/data/sample_both.fastq.gz,/data/sample_both_R1.fastq.gz,/data/sample_both_R2.fastq.gz
negative_control,/data/neg.fastq.gz,,
```

Repeated `sample` IDs are allowed when the same sample has been sequenced multiple times. Those rows are kept as distinct inputs in the current workflow, so sample IDs should still be chosen carefully and consistently.

## Typical Commands

Automatic reference selection with the HCV preset:

```bash
nextflow run charlesfoster/h2seq \
  -profile docker \
  --input samplesheet.csv \
  --outdir results \
  --virus_preset hcv
```

By default, automatic reference selection uses competitive `minimap2` mapping against the reference panel.

Automatic reference selection from a user-supplied multifasta:

```bash
nextflow run charlesfoster/h2seq \
  -profile docker \
  --input samplesheet.csv \
  --outdir results \
  --possible_references references.fasta
```

Competitive minimap2 reference selection from a user-supplied multifasta:

```bash
nextflow run charlesfoster/h2seq \
  -profile docker \
  --input samplesheet.csv \
  --outdir results \
  --possible_references references.fasta \
  --reference_selection_tool minimap2
```

Skip reference selection and map directly to one reference:

```bash
nextflow run charlesfoster/h2seq \
  -profile docker \
  --input samplesheet.csv \
  --outdir results \
  --skip_reference_selection \
  --reference_fasta reference.fasta
```

## HCV-GLUE

`--run_hcv_glue` is optional and should be treated as a site-specific integration.

It requires:

- Docker on the execution host
- a running `gluetools-mysql` container
- an already installed HCV GLUE database/project in that container

This step does not follow the normal nf-core-style container handling used by the rest of the pipeline.

## Negative Controls And Failed Samples

Samples with no usable data are handled explicitly so that they do not fail the entire run. In practice this means the pipeline can skip samples that have:

- no raw reads
- no reads after QC
- no mapped reads

Those samples still appear in the final `combined_results_summary.csv` with QC status and failure-reason fields describing what happened. Unexpected process failures still stop the pipeline.

## Mixed Infections

When reference selection retains multiple genotype references, the second mapping
pass is competitive. Each fragment (an individual long read or an Illumina read
pair) contributes to only one reference when its primary alignment is sufficiently
confident. Fragments with mapping quality below
`--mixed_assignment_min_mapq` (default `10`), and pairs whose mates map primarily to
different references, are excluded from genotype-specific coverage, variant calling,
and consensus generation.

The existing `combined_results_summary.csv` remains one row per sample and read type
and reports the main selected component. Its `component_fractions` field is derived
from the final confidently assigned fragments rather than the preliminary reference-
selection score. Use `reference_component_summary.csv` for one row per main or
secondary reference component, including component-specific coverage, mean depth,
assignment counts, and the number of ambiguous fragments. A separate PDF report is
generated for every retained component.

## Profiles

Containerised execution is the intended mode of operation:

- `docker`
- `singularity`
- `apptainer`
- `podman`

The bundled `test` profiles are for development and CI. Running without a container profile is possible only if the required tools are already installed on `PATH`.

## Reproducibility

Use a tagged release with `-r` for production analyses:

```bash
nextflow run charlesfoster/h2seq -r <tag> -profile docker --input samplesheet.csv --outdir results
```

Resume interrupted runs with:

```bash
nextflow run charlesfoster/h2seq -resume -profile docker --input samplesheet.csv --outdir results
```
