/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// nf-core modules //

include { FASTQC as FASTQC_RAW_SHORT              } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIMMED_SHORT          } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_RAW_LONG               } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIMMED_LONG           } from '../modules/nf-core/fastqc/main'
include { FASTP                                   } from '../modules/nf-core/fastp/main'
include { MULTIQC                                 } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                        } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                  } from '../subworkflows/local/utils_nfcore_h2seq_pipeline'
include { NANOQ                                   } from '../modules/nf-core/nanoq/main'
include { KALLISTO_INDEX                          } from '../modules/nf-core/kallisto/index/main'
include { KALLISTO_QUANT                          } from '../modules/nf-core/kallisto/quant/main'
include { SALMON_INDEX                            } from '../modules/nf-core/salmon/index/main'
include { SALMON_QUANT as SALMON_QUANT_LONG       } from '../modules/nf-core/salmon/quant/main'
include { SALMON_QUANT as SALMON_QUANT_SHORT      } from '../modules/nf-core/salmon/quant/main'
include { MINIMAP2_ALIGN                          } from '../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_SALMON } from '../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_AMPLICONCLIP                   } from '../modules/nf-core/samtools/ampliconclip/main'
include { SAMTOOLS_SORT                           } from '../modules/nf-core/samtools/sort/main'
include { MOSDEPTH                                } from '../modules/nf-core/mosdepth/main'
include { BWA_INDEX                               } from '../modules/nf-core/bwa/index/main'
include { BWA_MEM as MAP_PRIMERS                  } from '../modules/nf-core/bwa/mem/main'
include { BWA_MEM                                 } from '../modules/nf-core/bwa/mem/main'
include { SEQKIT_GREP                             } from '../modules/nf-core/seqkit/grep/main'
include { BEDTOOLS_BAMTOBED                       } from '../modules/nf-core/bedtools/bamtobed/main'

// local modules //
include { SAMTOOLS_CONSENSUS        } from '../modules/local/samtools/consensus/main'
include { CALCULATE_READ_STATS      } from '../modules/local/custom/calculate_read_stats/main'
include { SELECT_BEST_REFERENCE     } from '../modules/local/custom/select_best_reference/main'
include { REMOVE_EMPTY_SEQUENCES    } from '../modules/local/custom/remove_empty_sequences/main'
include { HCV_GLUE                  } from '../modules/local/custom/hcv_glue/main'
include { SPLIT_CONSENSUS_GENOMES   } from '../modules/local/custom/split_consensus_genomes/main'

// local subworkflows //
include { LONG_READ_MAPPING        } from '../subworkflows/local/long_read_mapping'
include { SHORT_READ_MAPPING       } from '../subworkflows/local/short_read_mapping'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//TODO: make sure that the 'ghost' short reads folder isn't created when running long reads only
//TODO: update README to fully describe the pipeline, its parameters, and its installation
workflow H2SEQ {

    take:
    ch_raw_long_reads
    ch_raw_short_reads

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        SET UP FILE PATHS
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    if (!params.skip_reference_selection){
        if ( params.virus_preset == "hcv" ) {
            possible_references =  "${projectDir}/assets/reference_data/hcv_references.fasta"
        } else {
            possible_references = params.possible_references
        }

        if (possible_references) {
            ch_reference_fasta = Channel.fromPath(possible_references)
            ch_reference_fasta = ch_reference_fasta
                .map { fasta ->
                    def meta = [id: 'possible_references']
                    return [meta, fasta]
                }
        } else {
            ch_reference_fasta = Channel.empty()
        }
    } else {
        ch_reference_fasta = Channel.fromPath(params.reference_fasta)
        ch_reference_fasta = ch_reference_fasta
            .map { fasta ->
                def meta = [id: 'best_reference', ref_type: 'BEST']
                return [meta.id, meta, fasta]
            }
    }

    ch_fastp_adapter_path = params.fastp_adapter_path ? file(params.fastp_adapter_path) : []

    /*
    ================================================================================
                                    Preprocessing and QC for long reads
    ================================================================================
    */

    // TODO: add an appropriate long read QC tool before trimming

    // TODO: add nanoq out to multiqc
    NANOQ (
        ch_raw_long_reads,
        "fastq"
    )
    ch_clean_reads_long = NANOQ.out.reads
    ch_versions = ch_versions.mix(NANOQ.out.versions)

    // TODO: add an appropriate long read QC tool after trimming

    /*
    ================================================================================
                                    Preprocessing and QC for short reads
    ================================================================================
    */

    FASTQC_RAW_SHORT (
        ch_raw_short_reads
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW_SHORT.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC_RAW_SHORT.out.versions.first())

    FASTP (
        ch_raw_short_reads,
        ch_fastp_adapter_path,
        false,
        false
    )

    ch_clean_reads_short = FASTP.out.reads
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    FASTQC_TRIMMED_SHORT (
        ch_clean_reads_short
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMMED_SHORT.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC_TRIMMED_SHORT.out.versions.first())


    /*
    ================================================================================
                                    Reference selection section
    ================================================================================
    */

    ch_clean_reads_combined = ch_clean_reads_short.mix(ch_clean_reads_long)

    if (!params.skip_reference_selection){
        if (params.reference_selection_tool == "kallisto") {
            CALCULATE_READ_STATS(
                ch_clean_reads_combined
            )
            ch_reads_and_stats = CALCULATE_READ_STATS.out.reads_and_stats
            ch_versions = ch_versions.mix(CALCULATE_READ_STATS.out.versions)

            KALLISTO_INDEX (
                ch_reference_fasta
            )
            ch_reference_fasta_index = KALLISTO_INDEX.out.index
            ch_versions = ch_versions.mix(KALLISTO_INDEX.out.versions)

            ch_quant_input = ch_reads_and_stats
                .combine(ch_reference_fasta_index)

            KALLISTO_QUANT (
                ch_quant_input
            )
            ch_abundance_tsv = KALLISTO_QUANT.out.tsv
            ch_versions = ch_versions.mix(KALLISTO_QUANT.out.versions)
        } else if (params.reference_selection_tool == "salmon") {
            // Focus on (potential) long reads first
            ch_map_for_salmon_long = ch_clean_reads_long
                .combine(ch_reference_fasta)

            MINIMAP2_ALIGN_SALMON (
                ch_map_for_salmon_long,
                true, // output in bam format
                false, //sort output
                "bai",
                false,
                false
            )
            ch_versions = ch_versions.mix(MINIMAP2_ALIGN_SALMON.out.versions)
            ch_ref_for_salmon = ch_reference_fasta
                .map { _meta, fasta ->
                    [fasta]
                }

            ch_input_for_salmon_long = MINIMAP2_ALIGN_SALMON.out.bam
                .combine(ch_ref_for_salmon)

            SALMON_QUANT_LONG (
                ch_input_for_salmon_long,
                [],
                true,
                "A"
            )
            ch_versions = ch_versions.mix(SALMON_QUANT_LONG.out.versions)
            ch_abundance_tsv_long = SALMON_QUANT_LONG.out.tsv

            // Focus on (potential) short reads second

            SALMON_INDEX (
                ch_ref_for_salmon
            )

            ch_versions = ch_versions.mix(SALMON_INDEX.out.versions)
            ch_salmon_idx = SALMON_INDEX.out.index.first()

            ch_input_for_salmon_short = ch_clean_reads_short
                .combine(ch_ref_for_salmon)

            SALMON_QUANT_SHORT (
                ch_input_for_salmon_short,
                ch_salmon_idx,
                false,
                "A"
            )
            ch_versions = ch_versions.mix(SALMON_QUANT_SHORT.out.versions)
            ch_abundance_tsv_short = SALMON_QUANT_SHORT.out.tsv

            // combine the results
            ch_abundance_tsv = ch_abundance_tsv_long
                .mix(ch_abundance_tsv_short)
        }

        SELECT_BEST_REFERENCE (
            ch_abundance_tsv
        )

        ch_best_ref_tsv = SELECT_BEST_REFERENCE.out.best_ref_tsv
        // TODO: GET THE BEST REF TXT EVEN WHEN SKIPPING REFERENCE SELECTION
        ch_best_ref_txt = SELECT_BEST_REFERENCE.out.best_ref_txt
        ch_alt_ref_txt = SELECT_BEST_REFERENCE.out.alt_ref_txt
        ch_versions = ch_versions.mix(SELECT_BEST_REFERENCE.out.versions)

        ch_alt_seqkit_input = ch_reference_fasta
            .combine(ch_alt_ref_txt)

        SEQKIT_GREP (
            ch_alt_seqkit_input
        )

        ch_versions = ch_versions.mix(SEQKIT_GREP.out.versions)

        ch_best_ref_fasta = SEQKIT_GREP.out.filter
            .map{ meta, fasta ->
                [meta.id, meta, fasta]
            }

        ch_best_ref_long = ch_best_ref_fasta
            .filter { _id, meta, _fasta -> meta.long_reads == true }

        ch_best_ref_short = ch_best_ref_fasta
            .filter { _id, meta, _fasta -> meta.long_reads == false }
    } else {
        ch_best_ref_long = ch_reference_fasta
        ch_best_ref_short = ch_reference_fasta
        ch_best_ref_short.view()
        ch_best_ref_long.view()
    }

    /*
    ================================================================================
                                    Mapping and primer clipping
    ================================================================================
    */
    ch_consensus_bam_long = Channel.empty()
    ch_consensus_bam_short = Channel.empty()

    // Note: here we have separate subworkflows for short and long reads
    //       Might have been able to avoid this with careful multi-key combines (see Consensus section below),
    //       but for sustained development this seemed like a better choice.
    LONG_READ_MAPPING  ( ch_best_ref_long, ch_clean_reads_long )
    SHORT_READ_MAPPING ( ch_best_ref_short, ch_clean_reads_short )

    ch_versions = ch_versions.mix(LONG_READ_MAPPING.out.versions)
    ch_versions = ch_versions.mix(SHORT_READ_MAPPING.out.versions)

    ch_consensus_bam_long = LONG_READ_MAPPING.out.consensus_bam
    ch_consensus_bam_short = SHORT_READ_MAPPING.out.consensus_bam
    ch_consensus_bam = ch_consensus_bam_long.mix(ch_consensus_bam_short)

    /*
    ================================================================================
                                    Consensus generation
    ================================================================================
    */

    SAMTOOLS_CONSENSUS (
        ch_consensus_bam
    )

    ch_consensus_fa = SAMTOOLS_CONSENSUS.out.fasta
        .map { meta, fasta ->
            return [meta.id, meta.long_reads, meta, fasta]
        }

    ch_versions = ch_versions.mix(SAMTOOLS_CONSENSUS.out.versions)

    ch_best_ref_txt_keyed = ch_best_ref_txt
        .map { meta, txt ->
            return [meta.id, meta.long_reads, meta, txt]
        }

    ch_split_input = ch_consensus_fa
        .combine(ch_best_ref_txt_keyed, by:[0,1])
        .map {_id, _is_long, meta, fasta, _meta2, txt ->
            [meta, fasta, txt]
        }

    // split potential multiple consensuses into individual files
    // the best will be in "*.main.fa"
    // the others will be in "*.alt#.fa" (e.g. alt1, alt2, ...)

    SPLIT_CONSENSUS_GENOMES (
        ch_split_input
    )

    ch_split_consensuses = SPLIT_CONSENSUS_GENOMES.out.fastas

    /*
    ================================================================================
                                    HCV Analysis
    ================================================================================
    */

    if ( params.run_hcv_glue ){

        REMOVE_EMPTY_SEQUENCES (
            ch_split_consensuses
        )

        // Handy hint: transpose operator “transposes” each tuple from a source channel
        //      by flattening any nested list in each tuple, emitting each nested item separately.
        // Practical example: if the channel has [[meta], fasta1, fasta2], then transpose will lead
        //      to [[[meta], fasta1], [[meta],fasta2]]
        ch_glue_fa = REMOVE_EMPTY_SEQUENCES.out.fasta
            .transpose()

        HCV_GLUE (
            ch_glue_fa
        )
    }

    /*
    ================================================================================
                                    Version parsing and MultiQC
    ================================================================================
    */

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
