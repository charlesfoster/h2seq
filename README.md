[![GitHub Actions CI Status](https://github.com/charlesfoster/h2seq/actions/workflows/ci.yml/badge.svg)](https://github.com/charlesfoster/h2seq/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/charlesfoster/h2seq/actions/workflows/linting.yml/badge.svg)](https://github.com/charlesfoster/h2seq/actions/workflows/linting.yml)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)
[![nf-core template](https://img.shields.io/badge/nf--core%20template-3.5.2-%2304B7B4?logo=nf-core&logoColor=white)](https://nf-co.re/)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A523.10.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/charlesfoster/h2seq)

## Introduction

**charlesfoster/h2seq** is a bioinformatics pipeline that has been designed to analyse molecular sequencing data of viruses for the H2Seq study. Accordingly, it has been designed with HCV and HIV in mind, but in theory should work with any (most?) viruses. The workflow handles both long-read (ONT) and short-read (Illumina) sequencing data, and can handle tiled amplicon sequencing and/or shotgun/metagenomic/capture probe sequencing.

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
4. Variant-based consensus genome generation with `Clair3` (ONT) or `LoFreq` (Illumina), followed by [`bcftools consensus`](https://samtools.github.io/bcftools/bcftools.html#consensus)

For ONT data, the default `--ont_min_snv_af` is `0.15`. Values below this are not recommended because Clair3 was trained on human data with allele frequencies in the range of 15%-100%.

> [!IMPORTANT]
> Additional options have been included over time, and this documentation will be updated accordingly at some stage. For now, just view all possible options by running the `--help` command (see: 'Usage' section below).

### Metagenomic Sequencing

Currently there are no 'specialised' modules for metagenomics data. Just run the pipeline as if your reads are derived from amplicon sequencing, but use the `--skip_primer_trimming` option (see: 'Usage' section below).

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

- _estimation of the amino acid consequences of SNPs/indels to aid with drug resistance analysis_ (currently only provided via the integrated HCV Glue workflow)
- _placement of input samples into a phylogenetic tree_
- _host filtration_
- _use `pycoQC` for Nanopore QC instead of fastQC_

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
> Development has focused on dependencies being handled by Docker or Apptainer, i.e. by including `-profile docker` or `-profile apptainer`. Currently `-profile conda` will _NOT_ work, but will work in the future.

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
