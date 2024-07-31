/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// nf-core modules //

include { FASTQC as FASTQC_RAW                    } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIMMED                } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                 } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                        } from 'plugin/nf-validation'
include { paramsSummaryMultiqc                    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                  } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                  } from '../subworkflows/local/utils_nfcore_h2seq_pipeline'
include { NANOQ                                   } from '../modules/nf-core/nanoq/main'                         
include { KALLISTO_INDEX                          } from '../modules/nf-core/kallisto/index/main'       
include { KALLISTO_QUANT                          } from '../modules/nf-core/kallisto/quant/main'       
include { SALMON_INDEX                            } from '../modules/nf-core/salmon/index/main'
include { SALMON_QUANT                            } from '../modules/nf-core/salmon/quant/main'
include { MINIMAP2_ALIGN                          } from '../modules/nf-core/minimap2/align/main'       
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_SALMON } from '../modules/nf-core/minimap2/align/main'       
include { SAMTOOLS_AMPLICONCLIP                   } from '../modules/nf-core/samtools/ampliconclip/main'
include { SAMTOOLS_SORT                           } from '../modules/nf-core/samtools/sort/main'
include { MOSDEPTH                                } from '../modules/nf-core/mosdepth/main'                   
include { BWA_INDEX                               } from '../modules/nf-core/bwa/index/main'                 
include { BWA_MEM                                 } from '../modules/nf-core/bwa/mem/main'                     
include { SEQKIT_GREP                             } from '../modules/nf-core/seqkit/grep/main'             
include { BEDTOOLS_BAMTOBED                       } from '../modules/nf-core/bedtools/bamtobed/main'

// local modules //
include { SAMTOOLS_CONSENSUS        } from '../modules/local/samtools/consensus/main'
include { CALCULATE_READ_STATS      } from '../modules/local/custom/calculate_read_stats/main'
include { SELECT_BEST_REFERENCE     } from '../modules/local/custom/select_best_reference/main'
include { REMOVE_EMPTY_SEQUENCES    } from '../modules/local/custom/remove_empty_sequences/main'
include { HCV_GLUE                  } from '../modules/local/custom/hcv_glue/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SET UP FILE PATHS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

if (!params.skip_reference_selection){
    if (params.possible_references) {
        ch_reference_fasta = Channel.fromPath(params.possible_references)
        ch_reference_fasta = ch_reference_fasta
            .map { fasta ->
                def meta = [id: 'possible_references']
                return [meta, fasta]
            }
    } else {
        ch_reference_fasta = Channel.empty()
    }
}

if (params.primer_fasta) {
    ch_primer_fasta = Channel.fromPath(params.primer_fasta)
    ch_primer_fasta = ch_primer_fasta
        .map { fasta ->
            def meta = [id: 'primers']
            return [meta, fasta]
        }
} else {
    ch_primer_fasta = Channel.empty()
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow H2SEQ {

    take:
    ch_raw_long_reads
    ch_raw_short_reads

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    /*
    ================================================================================
                                    Preprocessing and QC for long reads
    ================================================================================
    */

    FASTQC_RAW (
        ch_raw_long_reads
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first())

    // TODO nf-core: add nanoq out to multiqc
    NANOQ (
        ch_raw_long_reads,
        "fastq"
    )
    ch_clean_reads = NANOQ.out.reads
    ch_versions = ch_versions.mix(NANOQ.out.versions)

    FASTQC_TRIMMED (
        ch_clean_reads
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMMED.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC_TRIMMED.out.versions.first())

    /*
    ================================================================================
                                    Reference selection section
    ================================================================================
    */
    if (!params.skip_reference_selection){
        if (params.reference_selection_tool == "kallisto") {
            CALCULATE_READ_STATS( 
                ch_clean_reads
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
            ch_map_for_salmon = ch_clean_reads
                .combine(ch_reference_fasta)

            MINIMAP2_ALIGN_SALMON (
                ch_map_for_salmon,
                true, // output in bam format
                false, //sort output
                "bai",
                false,
                false
            )
            ch_versions = ch_versions.mix(MINIMAP2_ALIGN_SALMON.out.versions)
            ch_ref_for_salmon = ch_reference_fasta
                .map { meta, fasta ->
                    [fasta]
                }
            ch_alignment_for_salmon = MINIMAP2_ALIGN_SALMON.out.bam
                .combine(ch_ref_for_salmon)

            SALMON_QUANT (
                ch_alignment_for_salmon,
                [],
                true,
                "A"
            )
            ch_versions = ch_versions.mix(SALMON_QUANT.out.versions)
            ch_abundance_tsv = SALMON_QUANT.out.tsv
        }
        SELECT_BEST_REFERENCE (
            ch_abundance_tsv
        )
        ch_best_ref_tsv = SELECT_BEST_REFERENCE.out.best_ref_tsv
        ch_best_ref_txt = SELECT_BEST_REFERENCE.out.best_ref_txt
        ch_versions = ch_versions.mix(SELECT_BEST_REFERENCE.out.versions)

        ch_seqkit_input = ch_reference_fasta
            .combine(ch_best_ref_txt)
        
        SEQKIT_GREP (
            ch_seqkit_input
        )

        ch_versions = ch_versions.mix(SEQKIT_GREP.out.versions)
    
        ch_reads_for_minimap = ch_clean_reads
            .map{ meta, reads ->
                [meta.id, meta, reads]
            }

        ch_best_ref_fasta = SEQKIT_GREP.out.filter
            .map{ meta, fasta ->
                [meta.id, meta, fasta]
            }
        
        ch_fasta_for_samtools_sort = ch_best_ref_fasta
       
        ch_minimap2_input = ch_reads_for_minimap
            .combine(ch_best_ref_fasta, by:0)
            .map{ id, meta1, reads, meta2, fasta ->
                [meta1, reads, meta2, fasta]
            }
    } else {
        ch_best_ref_fasta = Channel.fromPath(params.reference_fasta)
            .map { fasta ->
                def meta = [id: 'best_reference']
                return [meta.id, meta, fasta]
            }
        ch_fasta_for_samtools_sort = ch_best_ref_fasta
        .map { meta, fasta ->
            [meta.id,meta,fasta]
        }
        ch_minimap2_input = ch_clean_reads
            .combine(ch_best_ref_fasta)
    }

    /*
    ================================================================================
                                    Mapping and primer clipping
    ================================================================================
    */

    MINIMAP2_ALIGN (
        ch_minimap2_input,
        true,
        true,
        "bai",
        false,
        false
    )
    ch_minimap2_bam = MINIMAP2_ALIGN.out.bam
    ch_minimap2_bam_idx = MINIMAP2_ALIGN.out.index
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

    if (!params.skip_primer_trimming){
        if (!params.primer_bed){
            
            ch_bwa_input_fasta = ch_best_ref_fasta
                .map{id,meta,fasta ->
                    [meta, fasta]
                }

            BWA_INDEX (
                ch_bwa_input_fasta
            )
            ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
            
            ch_bwa_mem_input = ch_primer_fasta
                .combine(BWA_INDEX.out.fasta_and_index)

            BWA_MEM (
                ch_bwa_mem_input,
                false
            )
            ch_primer_bam = BWA_MEM.out.bam
            ch_versions = ch_versions.mix(BWA_MEM.out.versions)

            BEDTOOLS_BAMTOBED (
                ch_primer_bam
            )
            ch_primer_bed = BEDTOOLS_BAMTOBED.out.bed
                .map { meta, bed ->
                    [ meta.id, meta, bed ]
                }

            ch_samtools_ampliconclip_input = ch_minimap2_bam
                .map { meta, bam ->
                    [meta.id, meta, bam]
                }
                .combine(ch_primer_bed, by:0)
                .map { id, meta, bam, meta2, bed ->
                    [meta, bam, bed]
                }
            ch_versions = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions)
        } else {
            ch_primer_bed = file(params.primer_bed, checkIfExists: true)
                .map { bed ->
                    def meta = [id: 'primer_bed']
                return [meta, bed]
            }

            ch_samtools_ampliconclip_input = ch_minimap2_bam
                .combine(ch_primer_bed)
        }
        
        SAMTOOLS_AMPLICONCLIP (
            ch_samtools_ampliconclip_input,
            true,
            false
        )

        ch_clipped_bam = SAMTOOLS_AMPLICONCLIP.out.bam
        ch_versions = ch_versions.mix(SAMTOOLS_AMPLICONCLIP.out.versions)
        
        ch_samtools_sort_input = ch_clipped_bam
            .map {meta, bam ->
                [meta.id, meta, bam]
            }
            .combine(ch_fasta_for_samtools_sort, by:0)
            .map {id, meta, bam, meta2, fasta ->
                [meta, bam, fasta]
            }

        SAMTOOLS_SORT (
            ch_samtools_sort_input
        )

        ch_consensus_bam = SAMTOOLS_SORT.out.bam
        ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)
       
    } else {
        ch_primer_bed = Channel.empty()
        ch_consensus_bam = ch_minimap2_bam
    }

    /*
    ================================================================================
                                    Consensus generation
    ================================================================================
    */

    SAMTOOLS_CONSENSUS (
        ch_consensus_bam
    ) 

    ch_consensus_fa = SAMTOOLS_CONSENSUS.out.fasta

    /*
    ================================================================================
                                    HCV Analysis
    ================================================================================
    */

    if ( params.run_hcv_glue ){
        REMOVE_EMPTY_SEQUENCES (
            ch_consensus_fa
        )

        ch_glue_fa = REMOVE_EMPTY_SEQUENCES.out.fasta

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
