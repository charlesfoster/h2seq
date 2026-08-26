[![GitHub Actions CI Status](https://github.com/charlesfoster/h2seq/actions/workflows/ci.yml/badge.svg)](https://github.com/charlesfoster/h2seq/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/charlesfoster/h2seq/actions/workflows/linting.yml/badge.svg)](https://github.com/charlesfoster/h2seq/actions/workflows/linting.yml)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)
[![nf-core template](https://img.shields.io/badge/nf--core%20template-3.5.2-%2304B7B4?logo=nf-core&logoColor=white)](https://nf-co.re/)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A526.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/charlesfoster/h2seq)

## Introduction

**charlesfoster/h2seq** is a bioinformatics pipeline for viral sequencing data developed for the H2Seq study. It is currently focused on HCV-oriented analysis, but it is intended to work with other small viral genomes given an appropriate reference set. The workflow supports ONT long reads, Illumina short reads, or mixed-modality runs, and can be used for tiled amplicon data as well as untargeted sequencing data where primer trimming is skipped.

### What The Pipeline Does

At a high level, the workflow performs:

1. Input validation and modality-aware preprocessing for long and/or short reads
2. Read QC and filtering/trimming
3. Optional automatic reference selection from a multifasta reference set
4. Reference-guided mapping
5. Optional primer coordinate inference and primer soft-clipping for amplicon data
6. Coverage and mapped-read summarisation
7. Explicit variant calling, filtering, annotation, and consensus sequence generation
8. Per-sample summary outputs, plots, and MultiQC reporting

### Tools Used

The current workflow uses the following primary tools:

- QC and read preprocessing: [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [`fastp`](https://github.com/OpenGene/fastp), [`NanoQ`](https://github.com/esteinig/nanoq), [`seqkit`](https://bioinf.shenwei.me/seqkit/)
- Reference selection: [`salmon`](https://github.com/COMBINE-lab/salmon), [`kallisto`](https://github.com/pachterlab/kallisto), or competitive [`minimap2`](https://github.com/lh3/minimap2) mapping
- Read alignment and alignment processing: [`minimap2`](https://github.com/lh3/minimap2), [`bwa`](https://github.com/lh3/bwa), [`samtools`](http://www.htslib.org/), [`bedtools`](https://github.com/arq5x/bedtools2)
- Coverage analysis: [`mosdepth`](https://github.com/brentp/mosdepth)
- Variant calling:
  - ONT: [`Clair3`](https://github.com/HKU-BAL/Clair3)
  - Illumina: [`LoFreq`](https://csb5.github.io/lofreq/)
- Variant preparation and consensus generation: [`bcftools`](https://samtools.github.io/bcftools/)
- Reporting: [`MultiQC`](https://multiqc.info/)

### Workflow Summary

#### Read QC and preprocessing

- Long reads are summarised with `seqkit stats`, filtered with `NanoQ`, and then summarised again after QC.
- Short reads are summarised with `seqkit stats`, assessed with `FastQC`, filtered/trimmed with `fastp`, and then reassessed with `seqkit stats` and `FastQC`.
- Samples with no usable reads, no reads after QC, or no mapped reads are handled explicitly and are retained in run-level summaries with explanatory notes instead of crashing the whole run.

#### Reference selection

- The pipeline can automatically select a best-matching reference from a multifasta, either from a user-supplied `--possible_references` file or from the bundled HCV reference set when `--virus_preset hcv` is used.
- Reference selection uses competitive `minimap2` mapping by default, with `salmon` and `kallisto` available via `--reference_selection_tool`. The `minimap2` mode maps reads to the whole reference panel and ranks primary alignment evidence rather than abundance.
- Alternatively, automatic reference selection can be skipped entirely with `--skip_reference_selection` and a fixed `--reference_fasta`.

#### Mapping and primer trimming

- Long reads and short reads are mapped against the selected reference using dedicated long-read and short-read mapping subworkflows.
- Primer trimming is optional and is intended for tiled amplicon data. Primer coordinates can be inferred from a primer FASTA using `bwa` and `bedtools`, and primer regions are then soft-clipped with `samtools ampliconclip`.
- For untargeted data such as shotgun, capture-probe, or metagenomic-style viral sequencing, primer trimming should normally be skipped with `--skip_primer_trimming`.

#### Variant calling and consensus generation

- Variant calling is platform-specific:
  - ONT data are processed with `Clair3`
  - Illumina data are processed with `LoFreq`, including an indel-quality preparation step
- The workflow then prepares the resulting VCFs into a harmonised format, applies explicit filtering and annotation, generates consensus mask regions, and builds the final consensus sequence with `bcftools consensus`.
- This means the current consensus path is variant-driven rather than relying on `samtools consensus` alone.

#### Coverage, reporting, and HCV-specific outputs

- Coverage is calculated with `mosdepth`, and the pipeline produces depth summaries, mapped-read counts, and per-sample coverage metrics.
- For mixed infections, competitively mapped fragments are assigned to one retained reference; low-confidence or mate-conflicting fragments are excluded from component-specific analysis. The combined summary retains one row per sample/read type with final component fractions, while `reference_component_summary.csv` and reference-specific PDF reports describe every retained component.
- Samples with no primary mapped reads, low selected-reference coverage, or an all-`N` main consensus before HCV-GLUE are retained in the combined summary as `qc_fail`; genotype/subtype assignments for those rows are flagged as unreliable.
- Custom summary tables and plots are generated for integration into `MultiQC`.
- When enabled, the optional `HCV-GLUE` integration adds HCV-specific reporting and genomic-region coverage summaries for final consensus genomes.

For ONT data, the default `--ont_min_snv_af` is `0.15`. Values below this are not recommended because Clair3 models were trained in a setting where lower-frequency calls are less reliable.

> [!IMPORTANT]
> Additional options have been included over time, and this documentation will be updated accordingly at some stage. For now, just view all possible options by running the `--help` command (see: 'Usage' section below).

### Metagenomic Sequencing

There is no separate metagenomics-specific branch of the workflow. Instead, untargeted viral sequencing data should generally be run through the standard workflow with primer trimming disabled using `--skip_primer_trimming`.

### Specialised Modules

**HCV-GLUE**

Given the initial focus of this pipeline for the H2Seq study, a module has been included to run the excellent ([`HCV-GLUE`](https://github.com/giffordlabcvr/HCV-GLUE)) tool. While the tool can do _a lot_, in this case the use is for generating reports based on HCV consensus genomes. `HCV-GLUE` is an intentionally non-standard optional step in this pipeline because the upstream Docker installation model relies on a persistent MySQL-backed GLUE setup, rather than a single self-contained command line tool invocation. Accordingly, unlike the other modules in the pipeline, successful use of `HCV-GLUE` requires manual installation and preparation outside the normal nf-core container/conda flow.

If you use `--run_hcv_glue`, you should assume the following prerequisites:

- Docker must be available on the host running the task.
- A `gluetools-mysql` container must already be running.
- The HCV GLUE project/database must already be installed in that container.
- This optional step is therefore less portable than the rest of the pipeline and is best treated as a site-specific integration.

Please see the upstream Docker installation instructions for `HCV-GLUE`: <https://github.com/giffordlabcvr/HCV-GLUE/wiki/Docker-Installation>

### Future considerations:

- _placement of input samples into a phylogenetic tree_
- _host filtration_

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that (minimally) looks as follows:

`samplesheet.csv`:

```csv
sample,long_reads,short_reads_1,short_reads_2
sample2,sample1.long.fastq.gz,sample1_S1_L002_R1_001.fastq.gz,sample1_S1_L002_R2_001.fastq.gz
```

Each row points the workflow to fastq files associated with a sample. Each sample can have one long reads file, and one short read fastq file (single-end) or a pair of short read fastq files (paired end).

> [!IMPORTANT]
> Currently single-end short read functionality is not tested and might not work.

Extra columns can be added to the spreadsheet as required for local purposes (e.g., tracking barcodes/serial numbers etc.), but must occur _after_ the mandatory minimal columns and cannot use any of the mandatory minimal column names.

Now, you can minimally run the pipeline using:

```bash
nextflow run charlesfoster/h2seq \
   -profile <docker/apptainer/...> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!IMPORTANT]
> The pipeline is primarily developed and tested with containerised execution using `-profile docker`, `-profile apptainer`, or `-profile singularity`. When running with `-profile conda` it is unlikely but possible that you will run into build errors in a Linux environment, and certain if you are working on a Mac. _Please_ use `-profile docker` on a Mac, or run within a Linux VM.

> [!NOTE]
> The optional `--run_hcv_glue` step is an exception to the usual dependency model described above. It shells out to Docker directly and expects a pre-existing GLUE database/container setup as described in the `HCV-GLUE` section.

Available parameters to be configured can be viewed by running:

```bash
nextflow run charlesfoster/h2seq --help
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

For more details and further functionality, please refer to the local documentation in [docs/usage.md](docs/usage.md) and [docs/output.md](docs/output.md).

## Credits

charlesfoster/h2seq was originally written by Charles Foster.

We thank the nf-core community for the reusable workflow structure and module ecosystem that this pipeline builds on.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
