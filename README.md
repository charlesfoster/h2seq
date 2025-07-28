[![GitHub Actions CI Status](https://github.com/charlesfoster/h2seq/actions/workflows/ci.yml/badge.svg)](https://github.com/charlesfoster/h2seq/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/charlesfoster/h2seq/actions/workflows/linting.yml/badge.svg)](https://github.com/charlesfoster/h2seq/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.10.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/charlesfoster/h2seq)

## Introduction

**charlesfoster/h2seq** is a bioinformatics pipeline that has been designed to analyse molecular sequencing data of viruses for the H2Seq study. Accordingly, it has been designed with HCV and HIV in mind, but in theory should work with any (most?) viruses. The workflow will handle data from both long-read (ONT) and/or short-read (Illumina) (to be added) sequencing platforms, and can handle tiled amplicon sequencing and/or shotgun/metagenomic/capture probe sequencing.

### Quality Control

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Read filtering/trimming
   - Long reads: ([`NanoQ`](https://github.com/esteinig/nanoq))
   - Short reads: ([`fastp`](https://github.com/OpenGene/fastp))
3. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))

### Amplicon Sequencing

1. Selection of closest reference genome
   - Choice of ([`kallisto`](https://github.com/pachterlab/kallisto)) or ([`salmon`](https://github.com/COMBINE-lab/salmon))
2. Alignment of reads against closest reference genome using([`minimap2`](https://github.com/lh3/minimap2))
3. Masking of amplicon primer sequences
   - Determination of primer coordinates using ([`bwa`](https://github.com/lh3/bwa)) and ([`bedtools`](https://github.com/arq5x/bedtools2))
   - Soft clipping of primer regions with ([`samtools ampliconclip`](http://www.htslib.org/doc/samtools-ampliconclip.html))
4. Consensus genome generation with ([`samtools consensus`](http://www.htslib.org/doc/samtools-consensus.html))

> [!IMPORTANT]
> Additional options have been included over time, and this documentation will be updated accordingly at some stage. For now, just view all possible options by running the `--help` command (see: 'Usage' section below).

### Metagenomic Sequencing

Currently there are no 'specialised' modules for metagenomics data. Just run the pipeline as if your reads are derived from amplicon sequencing, but use the `--skip_primer_trimming` option (see: 'Usage' section below).

### Specialised Modules

**HCV-GLUE**

Given the initial focus of this pipeline for the H2Seq study, a module has been included to run the excellent ([`HCV-GLUE`](https://github.com/giffordlabcvr/HCV-GLUE)) tool. While the tool can do _a lot_, in this case the use is for generating reports based on HCV consensus genomes. The usage of `HCV-GLUE`, by design, needs SQL-like syntax and relies on a Docker daemon running in the background. While it makes sense for 'vanilla' usage of `HCV-GLUE`, unfortunately it does not make its incorporation into a pipeline strictly following the `nf-core` framework straightforward. Therefore, unlike all other modules in the pipeline, in which tools are automatically sourced from containers/conda, successful use of `HCV-GLUE` requires manual installation of it and its dependencies. Please see the ([installation instructions](http://hcv-glue.cvr.gla.ac.uk/#/aboutGlueProject)) for `HCV-GLUE` if you wish to run this module.

### Future considerations:

- _variant calling to produce a VCF file instead of/in addition to direct consensus generation_ (currently available with long reads only)
- _estimation of the amino acid consequences of SNPs/indels to aid with drug resistance analysis_ (currently only provided via the integrated HCV Glue workflow)
- _placement of input samples into a phylogenetic tree_
- _host filtration_
- _use `pycoQC` for Nanopore QC instead of fastQC_

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   -->

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
   -profile <docker/singularity/...> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!IMPORTANT]
> Development has focused on dependencies being handled by Docker or Singularity, i.e. by including `-profile docker` or `-profile singularity`. Currently `-profile conda` will _NOT_ work, but will work in the future.

Available parameters to be configured can be viewed by running:

```bash
nextflow run charlesfoster/h2seq --help
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/mag/usage).

## Credits

charlesfoster/h2seq was originally written by Charles Foster.

We thank the following people for their extensive assistance in the development of this pipeline:

- No one else yet!

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use charlesfoster/h2seq for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
