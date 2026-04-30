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
include { SEQKIT_STATS as SEQKIT_STATS_RAW_LONG   } from '../modules/nf-core/seqkit/stats/main'
include { SEQKIT_STATS as SEQKIT_STATS_CLEAN_LONG } from '../modules/nf-core/seqkit/stats/main'
include { SEQKIT_STATS as SEQKIT_STATS_RAW_SHORT  } from '../modules/nf-core/seqkit/stats/main'
include { SEQKIT_STATS as SEQKIT_STATS_CLEAN_SHORT } from '../modules/nf-core/seqkit/stats/main'
include { KALLISTO_INDEX                          } from '../modules/nf-core/kallisto/index/main'
include { KALLISTO_QUANT                          } from '../modules/nf-core/kallisto/quant/main'
include { SALMON_INDEX                            } from '../modules/nf-core/salmon/index/main'
include { SALMON_QUANT as SALMON_QUANT_LONG       } from '../modules/nf-core/salmon/quant/main'
include { SALMON_QUANT as SALMON_QUANT_SHORT      } from '../modules/nf-core/salmon/quant/main'
include { MINIMAP2_ALIGN                          } from '../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_SALMON } from '../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_AMPLICONCLIP                   } from '../modules/nf-core/samtools/ampliconclip/main'
include { SAMTOOLS_SORT                           } from '../modules/nf-core/samtools/sort/main'
include { MOSDEPTH as MOSDEPTH_GENOME             } from '../modules/nf-core/mosdepth/main'
include { BWA_INDEX                               } from '../modules/nf-core/bwa/index/main'
include { BWA_MEM as MAP_PRIMERS                  } from '../modules/nf-core/bwa/mem/main'
include { BWA_MEM                                 } from '../modules/nf-core/bwa/mem/main'
include { SEQKIT_GREP                             } from '../modules/nf-core/seqkit/grep/main'
include { BEDTOOLS_BAMTOBED                       } from '../modules/nf-core/bedtools/bamtobed/main'
include { SAMTOOLS_FAIDX                          } from '../modules/nf-core/samtools/faidx/main'

// local modules //
include { CALCULATE_READ_STATS      } from '../modules/local/custom/calculate_read_stats/main'
include { SELECT_BEST_REFERENCE     } from '../modules/local/custom/select_best_reference/main'
include { MINIMAP2_REFERENCE_INDEX  } from '../modules/local/custom/minimap2_reference_index/main'
include { MINIMAP2_COMPETITIVE_MAP  } from '../modules/local/custom/minimap2_competitive_map/main'
include { SELECT_REFERENCE_FROM_BAM } from '../modules/local/custom/select_reference_from_bam/main'
include { REMOVE_EMPTY_SEQUENCES    } from '../modules/local/custom/remove_empty_sequences/main'
include { HCV_GLUE                  } from '../modules/local/custom/hcv_glue/main'
include { SPLIT_CONSENSUS_GENOMES   } from '../modules/local/custom/split_consensus_genomes/main'
include { CREATE_PATTERN_FILE       } from '../modules/local/custom/create_pattern_file/main'
include { COVERAGE_METRICS          } from '../modules/local/custom/coverage_metrics/main'
include { BUILD_RUN_SUMMARY         } from '../modules/local/custom/build_run_summary/main'
include { BUILD_MULTIQC_SECTIONS    } from '../modules/local/custom/build_multiqc_sections/main'
include { COUNT_MAPPED_READS        } from '../modules/local/custom/count_mapped_reads/main'
include { REFERENCE_METADATA_FROM_FASTA } from '../modules/local/custom/reference_metadata_from_fasta/main'
include { GENERATE_WHOLE_GENOME_BED } from '../modules/local/custom/generate_whole_genome_bed/main'
include { LOFREQ_INDELQUAL          } from '../modules/local/custom/lofreq_indelqual/main'
include { LOFREQ_CALL               } from '../modules/local/custom/lofreq_call/main'
include { PREPARE_CLAIR3_VCF        } from '../modules/local/custom/prepare_clair3_vcf/main'
include { PREPARE_LOFREQ_VCF        } from '../modules/local/custom/prepare_lofreq_vcf/main'
include { COMPRESS_INDEX_VCF as COMPRESS_PREPARED_LOFREQ_VCF } from '../modules/local/custom/compress_index_vcf/main'
include { FILTER_VARIANTS           } from '../modules/local/custom/filter_variants/main'
include { ANNOTATE_VARIANTS         } from '../modules/local/custom/annotate_variants/main'
include { CREATE_CONSENSUS_MASK     } from '../modules/local/custom/create_consensus_mask/main'
include { BCFTOOLS_CONSENSUS        } from '../modules/local/custom/bcftools_consensus/main'
include { CLAIR3                    } from '../modules/local/custom/clair3/main'
include { PARSE_HCV_GLUE_COVERAGE   } from '../modules/local/custom/parse_hcv_glue_coverage/main'
include { PLOT_HCV_SUMMARY          } from '../modules/local/custom/plot_hcv_summary/main'
include { PLOT_DEPTH_SUMMARY        } from '../modules/local/custom/plot_depth_summary/main'
include { RENDER_HCV_REPORT         } from '../modules/local/custom/render_hcv_report/main'
include { RENDER_SUMMARY_REPORT     } from '../modules/local/custom/render_summary_report/main'

// local subworkflows //
include { LONG_READ_MAPPING        } from '../subworkflows/local/long_read_mapping'
include { SHORT_READ_MAPPING       } from '../subworkflows/local/short_read_mapping'

def parseSeqkitCount(statsFile) {
    def lines = statsFile.text.readLines().findAll { it?.trim() }
    if (!lines) {
        throw new IllegalStateException("Empty seqkit stats file: ${statsFile}")
    }

    def rows = lines.collect { it.split('\t', -1) as List }
    def header = rows[0]
    def countIdx = header.indexOf('num_seqs')

    if (countIdx == -1) {
        throw new IllegalStateException("Missing 'num_seqs' column in seqkit stats file: ${statsFile}")
    }
    if (rows.size() < 2) {
        throw new IllegalStateException("Expected at least one data row in seqkit stats file: ${statsFile}")
    }

    def counts = rows.tail().collect { cols ->
        if (cols.size() <= countIdx || !cols[countIdx].trim()) {
            throw new IllegalStateException("Invalid seqkit stats row in ${statsFile}: ${cols.join('\t')}")
        }
        try {
            cols[countIdx].trim() as Long
        } catch (NumberFormatException e) {
            throw new IllegalStateException("Non-numeric num_seqs value in ${statsFile}: ${cols[countIdx]}", e)
        }
    }

    if (counts.toSet().size() > 1) {
        throw new IllegalStateException("Inconsistent num_seqs values across seqkit stats rows in ${statsFile}: ${counts}")
    }

    return counts[0]
}

def parseMappedReadCount(countFile) {
    def raw = countFile.text.trim()
    if (!raw) {
        throw new IllegalStateException("Empty mapped read count file: ${countFile}")
    }

    try {
        return raw as Long
    } catch (NumberFormatException e) {
        throw new IllegalStateException("Non-numeric mapped read count in ${countFile}: ${raw}", e)
    }
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

    SEQKIT_STATS_RAW_LONG (
        ch_raw_long_reads
    )
    ch_long_raw_stats = SEQKIT_STATS_RAW_LONG.out.stats
    ch_versions = ch_versions.mix(SEQKIT_STATS_RAW_LONG.out.versions)

    ch_raw_long_counts = ch_long_raw_stats
        .map { meta, stats ->
            [meta.id, meta.long_reads, meta, parseSeqkitCount(stats)]
        }

    ch_raw_long_reads_ready = ch_raw_long_reads
        .map { meta, reads ->
            [meta.id, meta.long_reads, meta, reads]
        }
        .combine(ch_raw_long_counts, by: [0, 1])
        .filter { _id, _long_reads, _meta1, _reads, _meta2, count -> count > 0 }
        .map { _id, _long_reads, meta, reads, _meta2, _count ->
            [meta, reads]
        }

    NANOQ (
        ch_raw_long_reads_ready,
        "fastq"
    )
    ch_clean_reads_long = NANOQ.out.reads
    ch_versions = ch_versions.mix(NANOQ.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(NANOQ.out.stats.map { _meta, files -> files }.flatten())

    SEQKIT_STATS_CLEAN_LONG (
        ch_clean_reads_long
    )
    ch_long_clean_stats = SEQKIT_STATS_CLEAN_LONG.out.stats
    ch_versions = ch_versions.mix(SEQKIT_STATS_CLEAN_LONG.out.versions)

    ch_clean_long_counts = ch_long_clean_stats
        .map { meta, stats ->
            [meta.id, meta.long_reads, meta, parseSeqkitCount(stats)]
        }

    ch_clean_reads_long_ready = ch_clean_reads_long
        .map { meta, reads ->
            [meta.id, meta.long_reads, meta, reads]
        }
        .combine(ch_clean_long_counts, by: [0, 1])
        .filter { _id, _long_reads, _meta1, _reads, _meta2, count -> count > 0 }
        .map { _id, _long_reads, meta, reads, _meta2, _count ->
            [meta, reads]
        }

    // TODO: add an appropriate long read QC tool after trimming

    /*
    ================================================================================
                                    Preprocessing and QC for short reads
    ================================================================================
    */

    SEQKIT_STATS_RAW_SHORT (
        ch_raw_short_reads
    )
    ch_short_raw_stats = SEQKIT_STATS_RAW_SHORT.out.stats
    ch_versions = ch_versions.mix(SEQKIT_STATS_RAW_SHORT.out.versions)

    ch_raw_short_counts = ch_short_raw_stats
        .map { meta, stats ->
            [meta.id, meta.long_reads, meta, parseSeqkitCount(stats)]
        }

    ch_raw_short_reads_ready = ch_raw_short_reads
        .map { meta, reads ->
            [meta.id, meta.long_reads, meta, reads]
        }
        .combine(ch_raw_short_counts, by: [0, 1])
        .filter { _id, _long_reads, _meta1, _reads, _meta2, count -> count > 0 }
        .map { _id, _long_reads, meta, reads, _meta2, _count ->
            [meta, reads]
        }

    FASTQC_RAW_SHORT (
        ch_raw_short_reads_ready
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW_SHORT.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC_RAW_SHORT.out.versions.first())

    FASTP (
        ch_raw_short_reads_ready,
        ch_fastp_adapter_path,
        false,
        false,
        false
    )

    ch_clean_reads_short = FASTP.out.reads
    ch_versions = ch_versions.mix(FASTP.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect { it[1] })

    SEQKIT_STATS_CLEAN_SHORT (
        ch_clean_reads_short
    )
    ch_short_clean_stats = SEQKIT_STATS_CLEAN_SHORT.out.stats
    ch_versions = ch_versions.mix(SEQKIT_STATS_CLEAN_SHORT.out.versions)

    ch_clean_short_counts = ch_short_clean_stats
        .map { meta, stats ->
            [meta.id, meta.long_reads, meta, parseSeqkitCount(stats)]
        }

    ch_clean_reads_short_ready = ch_clean_reads_short
        .map { meta, reads ->
            [meta.id, meta.long_reads, meta, reads]
        }
        .combine(ch_clean_short_counts, by: [0, 1])
        .filter { _id, _long_reads, _meta1, _reads, _meta2, count -> count > 0 }
        .map { _id, _long_reads, meta, reads, _meta2, _count ->
            [meta, reads]
        }

    FASTQC_TRIMMED_SHORT (
        ch_clean_reads_short_ready
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMMED_SHORT.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC_TRIMMED_SHORT.out.versions.first())


    /*
    ================================================================================
                                    Reference selection section
    ================================================================================
    */

    ch_clean_reads_combined = ch_clean_reads_short_ready.mix(ch_clean_reads_long_ready)
    ch_long_sample_ids = ch_clean_reads_long_ready
        .map { meta, _reads -> [meta.id] }
        .unique()
    ch_short_sample_ids = ch_clean_reads_short_ready
        .map { meta, _reads -> [meta.id] }
        .unique()

    if (!params.skip_reference_selection){
        if (params.reference_selection_tool == "kallisto" || params.reference_selection_tool == "salmon") {
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
                ch_map_for_salmon_long = ch_clean_reads_long_ready
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

                ch_input_for_salmon_short = ch_clean_reads_short_ready
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
            ch_best_ref_txt = SELECT_BEST_REFERENCE.out.best_ref_txt
            ch_alt_ref_txt = SELECT_BEST_REFERENCE.out.alt_ref_txt
            ch_versions = ch_versions.mix(SELECT_BEST_REFERENCE.out.versions)
        } else if (params.reference_selection_tool == "minimap2") {
            MINIMAP2_REFERENCE_INDEX (
                ch_reference_fasta
            )
            ch_versions = ch_versions.mix(MINIMAP2_REFERENCE_INDEX.out.versions)

            ch_reference_panel_for_minimap2 = ch_reference_fasta
                .map { meta, fasta ->
                    [meta.id, meta, fasta]
                }
                .combine(
                    MINIMAP2_REFERENCE_INDEX.out.index.map { meta, index ->
                        [meta.id, meta, index]
                    },
                    by: 0
                )
                .map { _id, reference_meta, reference_fasta, _index_meta, reference_index ->
                    [reference_meta, reference_fasta, reference_index]
                }

            ch_competitive_map_input = ch_clean_reads_combined
                .combine(ch_reference_panel_for_minimap2)
                .map { meta, reads, reference_meta, reference_fasta, reference_index ->
                    [meta, reads, reference_meta, reference_fasta, reference_index]
                }

            MINIMAP2_COMPETITIVE_MAP (
                ch_competitive_map_input
            )
            ch_versions = ch_versions.mix(MINIMAP2_COMPETITIVE_MAP.out.versions)

            SELECT_REFERENCE_FROM_BAM (
                MINIMAP2_COMPETITIVE_MAP.out.bam
            )
            ch_best_ref_tsv = SELECT_REFERENCE_FROM_BAM.out.best_ref_tsv
            ch_best_ref_txt = SELECT_REFERENCE_FROM_BAM.out.best_ref_txt
            ch_alt_ref_txt = SELECT_REFERENCE_FROM_BAM.out.alt_ref_txt
            ch_versions = ch_versions.mix(SELECT_REFERENCE_FROM_BAM.out.versions)
        }

        ch_alt_seqkit_input = ch_reference_fasta
            .combine(ch_alt_ref_txt)

        SEQKIT_GREP (
            ch_alt_seqkit_input
        )

        ch_versions = ch_versions.mix(SEQKIT_GREP.out.versions)

        ch_best_ref_fasta = SEQKIT_GREP.out.filter
            .map{ meta, fasta ->
                [meta.id, meta.long_reads, meta, fasta]
            }

        ch_best_ref_long = ch_best_ref_fasta
            .filter { _sample_id, long_reads, _meta, _fasta -> long_reads }
            .map { sample_id, _long_reads, meta, fasta ->
                [sample_id, meta, fasta]
            }

        ch_best_ref_short = ch_best_ref_fasta
            .filter { _sample_id, long_reads, _meta, _fasta -> !long_reads }
            .map { sample_id, _long_reads, meta, fasta ->
                [sample_id, meta, fasta]
            }
    } else {
        // When skipping reference selection, we need to create reference channels
        // for each sample with the proper structure [sample_id, meta, fasta]
        // Create reference channels by combining sample IDs with the reference fasta
        // Structure: [sample_id, ref_meta, fasta]
        ch_best_ref_long = ch_long_sample_ids
            .combine(ch_reference_fasta)
            .map { sample_id, ref_id, ref_meta, fasta ->
                def new_meta = ref_meta + [id: sample_id, long_reads: true]
                return [sample_id, new_meta, fasta]
            }

        ch_best_ref_short = ch_short_sample_ids
            .combine(ch_reference_fasta)
            .map { sample_id, ref_id, ref_meta, fasta ->
                def new_meta = ref_meta + [id: sample_id, long_reads: false]
                return [sample_id, new_meta, fasta]
            }

        REFERENCE_METADATA_FROM_FASTA (
            ch_best_ref_long.mix(ch_best_ref_short)
        )
        ch_best_ref_tsv = REFERENCE_METADATA_FROM_FASTA.out.best_ref_tsv
        ch_best_ref_txt = REFERENCE_METADATA_FROM_FASTA.out.best_ref_txt
        ch_versions = ch_versions.mix(REFERENCE_METADATA_FROM_FASTA.out.versions)
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
    LONG_READ_MAPPING  ( ch_best_ref_long, ch_clean_reads_long_ready )
    SHORT_READ_MAPPING ( ch_best_ref_short, ch_clean_reads_short_ready )

    ch_versions = ch_versions.mix(LONG_READ_MAPPING.out.versions)
    ch_versions = ch_versions.mix(SHORT_READ_MAPPING.out.versions)

    ch_consensus_bam_long = LONG_READ_MAPPING.out.consensus_bam
    ch_consensus_bam_short = SHORT_READ_MAPPING.out.consensus_bam
    ch_consensus_bam_idx_long = LONG_READ_MAPPING.out.consensus_bam_idx
    ch_consensus_bam_idx_short = SHORT_READ_MAPPING.out.consensus_bam_idx
    ch_consensus_bam = ch_consensus_bam_long.mix(ch_consensus_bam_short)
    ch_consensus_bam_idx = ch_consensus_bam_idx_long.mix(ch_consensus_bam_idx_short)

    COUNT_MAPPED_READS (
        ch_consensus_bam
    )
    ch_versions = ch_versions.mix(COUNT_MAPPED_READS.out.versions)

    ch_mapped_read_counts = COUNT_MAPPED_READS.out.txt
        .map { meta, countTxt ->
            def count = parseMappedReadCount(countTxt)
            [meta.id, meta.long_reads, meta, count]
        }

    ch_consensus_bam_ready = ch_consensus_bam
        .map { meta, bam ->
            [meta.id, meta.long_reads, meta, bam]
        }
        .combine(ch_mapped_read_counts, by: [0, 1])
        .filter { _id, _long_reads, _meta1, _bam, _meta2, count -> count > 0 }
        .map { _id, _long_reads, meta, bam, _meta2, _count ->
            [meta, bam]
        }

    ch_best_ref_all = ch_best_ref_long.mix(ch_best_ref_short)

    SAMTOOLS_FAIDX (
        ch_best_ref_all
            .map { _id, meta, fasta ->
                [meta, fasta]
            }
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    ch_reference_fai = SAMTOOLS_FAIDX.out.fa_and_idx
        .map { meta, fasta, fai ->
            [meta.id, meta.long_reads, meta + [reference_name: fai.text.readLines()[0].split('\t')[0]], fasta, fai]
        }

    GENERATE_WHOLE_GENOME_BED (
        ch_reference_fai.map { _id, _long_reads, meta, _fasta, fai ->
            [meta, fai]
        }
    )
    ch_versions = ch_versions.mix(GENERATE_WHOLE_GENOME_BED.out.versions)

    ch_mosdepth_input = ch_consensus_bam
        .map { meta, bam ->
            [meta.id, meta.long_reads, meta, bam]
        }
        .combine(
            ch_consensus_bam_idx.map { meta, idx ->
                [meta.id, meta.long_reads, meta, idx]
            },
            by: [0, 1]
        )
        .combine(
            GENERATE_WHOLE_GENOME_BED.out.bed.map { meta, bed ->
                [meta.id, meta.long_reads, meta, bed]
            },
            by: [0, 1]
        )
        .map { _id, _long_reads, meta, bam, _meta2, idx, refMeta, bed ->
            [meta + [reference_name: refMeta.reference_name], bam, idx, bed]
        }

    ch_mosdepth_reference = ch_reference_fai
        .map { _id, _long_reads, meta, fasta, _fai ->
            [meta, fasta]
        }

    MOSDEPTH_GENOME (
        ch_mosdepth_input,
        ch_mosdepth_reference
    )
    ch_versions = ch_versions.mix(MOSDEPTH_GENOME.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH_GENOME.out.summary_txt.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH_GENOME.out.global_txt.collect { it[1] })

    ch_coverage_input = MOSDEPTH_GENOME.out.regions_bed
        .map { meta, genome_bed_gz ->
            [meta.id, meta.long_reads, meta, genome_bed_gz]
        }
        .combine(
            MOSDEPTH_GENOME.out.per_base_bed.map { meta, per_base_bed_gz ->
                [meta.id, meta.long_reads, meta, per_base_bed_gz]
            },
            by: [0, 1]
        )
        .map { _id, _long_reads, meta, genome_bed_gz, _meta2, per_base_bed_gz ->
            [meta, genome_bed_gz, per_base_bed_gz]
        }

    COVERAGE_METRICS (
        ch_coverage_input
    )
    ch_versions = ch_versions.mix(COVERAGE_METRICS.out.versions)

    /*
    ================================================================================
                                    Variant calling and consensus generation
    ================================================================================
    */

    ch_variant_call_input = ch_consensus_bam_ready
        .map { meta, bam ->
            [meta.id, meta.long_reads, meta, bam]
        }
        .combine(
            ch_consensus_bam_idx.map { meta, idx ->
                [meta.id, meta.long_reads, meta, idx]
            },
            by: [0, 1]
        )
        .combine(
            ch_reference_fai.map { _id, _long_reads, meta, fasta, fai ->
                [meta.id, meta.long_reads, meta, fasta, fai]
            },
            by: [0, 1]
        )
        .map { _id, _long_reads, meta, bam, _meta2, bam_idx, refMeta, fasta, fai ->
            [meta + [reference_name: refMeta.reference_name], bam, bam_idx, fasta, fai]
        }

    ch_clair3_input = ch_variant_call_input
        .filter { meta, _bam, _bam_idx, _fasta, _fai -> meta.long_reads }

    CLAIR3 (
        ch_clair3_input
    )
    ch_versions = ch_versions.mix(CLAIR3.out.versions)

    PREPARE_CLAIR3_VCF (
        CLAIR3.out.vcf
    )
    ch_versions = ch_versions.mix(PREPARE_CLAIR3_VCF.out.versions)

    ch_lofreq_indelqual_input = ch_variant_call_input
        .filter { meta, _bam, _bam_idx, _fasta, _fai -> !meta.long_reads }
        .map { meta, bam, _bam_idx, fasta, fai ->
            [meta, bam, fasta, fai]
        }

    LOFREQ_INDELQUAL (
        ch_lofreq_indelqual_input
    )
    ch_versions = ch_versions.mix(LOFREQ_INDELQUAL.out.versions)

    ch_lofreq_call_input = LOFREQ_INDELQUAL.out.bam
        .map { meta, bam, bai ->
            [meta.id, meta.long_reads, meta, bam, bai]
        }
        .combine(
            ch_reference_fai.map { _id, _long_reads, meta, fasta, fai ->
                [meta.id, meta.long_reads, meta, fasta, fai]
            },
            by: [0, 1]
        )
        .map { _id, _long_reads, meta, bam, bai, _meta2, fasta, fai ->
            [meta, bam, bai, fasta, fai]
        }

    LOFREQ_CALL (
        ch_lofreq_call_input
    )
    ch_versions = ch_versions.mix(LOFREQ_CALL.out.versions)

    PREPARE_LOFREQ_VCF (
        LOFREQ_CALL.out.vcf
    )
    ch_versions = ch_versions.mix(PREPARE_LOFREQ_VCF.out.versions)

    COMPRESS_PREPARED_LOFREQ_VCF (
        PREPARE_LOFREQ_VCF.out.vcf
    )
    ch_versions = ch_versions.mix(COMPRESS_PREPARED_LOFREQ_VCF.out.versions)

    ch_prepared_variants = PREPARE_CLAIR3_VCF.out.vcf
        .mix(COMPRESS_PREPARED_LOFREQ_VCF.out.vcf)

    FILTER_VARIANTS (
        ch_prepared_variants
    )
    ch_versions = ch_versions.mix(FILTER_VARIANTS.out.versions)

    ANNOTATE_VARIANTS (
        FILTER_VARIANTS.out.vcf
    )
    ch_versions = ch_versions.mix(ANNOTATE_VARIANTS.out.versions)

    CREATE_CONSENSUS_MASK (
        MOSDEPTH_GENOME.out.per_base_bed
    )
    ch_versions = ch_versions.mix(CREATE_CONSENSUS_MASK.out.versions)

    ch_bcftools_consensus_input = FILTER_VARIANTS.out.vcf
        .map { meta, vcf, vcf_idx ->
            [meta.id, meta.long_reads, meta, vcf, vcf_idx]
        }
        .combine(
            ch_reference_fai.map { _id, _long_reads, meta, fasta, _fai ->
                [meta.id, meta.long_reads, meta, fasta]
            },
            by: [0, 1]
        )
        .combine(
            CREATE_CONSENSUS_MASK.out.bed.map { meta, mask_bed ->
                [meta.id, meta.long_reads, meta, mask_bed]
            },
            by: [0, 1]
        )
        .map { _id, _long_reads, meta, vcf, vcf_idx, refMeta, fasta, _meta3, mask_bed ->
            [meta + [reference_name: refMeta.reference_name], vcf, vcf_idx, fasta, mask_bed]
        }

    BCFTOOLS_CONSENSUS (
        ch_bcftools_consensus_input
    )
    ch_versions = ch_versions.mix(BCFTOOLS_CONSENSUS.out.versions)

    ch_consensus_fa = BCFTOOLS_CONSENSUS.out.fasta
        .map { meta, fasta ->
            return [meta.id, meta.long_reads, meta, fasta]
        }

    if (!params.skip_reference_selection){
        // need to (potentially) split the consensuses into multiple files
        // only relevant when reference selection is not skipped
        // can make robust later
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
    } else {
        // When reference selection is skipped, still run SPLIT_CONSENSUS_GENOMES
        // First create pattern files from the reference FASTA
        ch_reference_for_pattern = ch_reference_fasta
            .map { ref_id, ref_meta, fasta ->
                return [ref_meta, fasta]
            }

        CREATE_PATTERN_FILE (
            ch_reference_for_pattern
        )
        // ch_versions = ch_versions.mix(CREATE_PATTERN_FILE.out.versions)

        // Create pattern files for each sample
        ch_sample_patterns = ch_consensus_fa
            .map { meta_id, meta_long_reads, meta, fasta ->
                return [meta, fasta]
            }
            .combine(CREATE_PATTERN_FILE.out.pattern)
            .map { sample_meta, consensus_fasta, pattern_meta, pattern_file ->
                return [sample_meta, consensus_fasta, pattern_file]
            }

        SPLIT_CONSENSUS_GENOMES (
            ch_sample_patterns
        )

        ch_split_consensuses = SPLIT_CONSENSUS_GENOMES.out.fastas
    }

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

        ch_hcv_reports = HCV_GLUE.out.report
    } else {
        ch_hcv_reports = Channel.empty()
    }

    ch_hcv_main_reports = ch_hcv_reports
        .filter { _meta, html -> html.name.contains("consensus_main") }

    ch_depth_plot_input = MOSDEPTH_GENOME.out.per_base_bed
        .map { meta, per_base_bed_gz ->
            [meta.id, meta.long_reads, meta, per_base_bed_gz]
        }
        .combine(
            COVERAGE_METRICS.out.summary.map { meta, summary ->
                [meta.id, meta.long_reads, meta, summary]
            },
            by: [0, 1]
        )
        .map { _id, _long_reads, meta, per_base_bed_gz, _meta2, summary ->
            [meta, per_base_bed_gz, summary]
        }

    PLOT_DEPTH_SUMMARY (
        ch_depth_plot_input
    )
    ch_versions = ch_versions.mix(PLOT_DEPTH_SUMMARY.out.versions)

    if (params.virus_preset == "hcv" && params.run_hcv_glue) {
        PARSE_HCV_GLUE_COVERAGE (
            ch_hcv_main_reports
        )
        ch_versions = ch_versions.mix(PARSE_HCV_GLUE_COVERAGE.out.versions)

        ch_hcv_plot_input = COVERAGE_METRICS.out.summary
            .map { meta, summary ->
                [meta.id, meta.long_reads, meta, summary]
            }
            .combine(
                PARSE_HCV_GLUE_COVERAGE.out.tsv.map { meta, tsv ->
                    [meta.id, meta.long_reads, meta, tsv]
                },
                by: [0, 1]
            )
            .map { _id, _long_reads, meta, summary, _meta3, hcv_tsv ->
                [meta, summary, hcv_tsv]
            }

        PLOT_HCV_SUMMARY (
            ch_hcv_plot_input
        )
        ch_versions = ch_versions.mix(PLOT_HCV_SUMMARY.out.versions)
    }

    ch_report_base = ch_best_ref_tsv
        .map { meta, best_ref_tsv ->
            [meta.id, meta.long_reads, meta, best_ref_tsv]
        }
        .combine(
            COVERAGE_METRICS.out.summary.map { meta, summary ->
                [meta.id, meta.long_reads, meta, summary]
            },
            by: [0, 1]
        )
        .combine(
            PLOT_DEPTH_SUMMARY.out.depth.map { meta, depth_plot ->
                [meta.id, meta.long_reads, meta, depth_plot]
            },
            by: [0, 1]
        )
        .map { _id, _long_reads, meta, best_ref_tsv, _meta2, summary, _meta3, depth_plot ->
            [meta.id, meta.long_reads, meta, best_ref_tsv, summary, depth_plot]
        }

    if (params.virus_preset == "hcv" && params.run_hcv_glue) {
        ch_report_grouped = ch_report_base
            .map { id, long_reads, meta, best_ref_tsv, summary, depth_plot ->
                [[id, long_reads], [meta, best_ref_tsv, summary, depth_plot], null]
            }
            .mix(
                PLOT_HCV_SUMMARY.out.feature.map { meta, feature_plot ->
                    [[meta.id, meta.long_reads], null, feature_plot]
                }
            )
            .groupTuple()

        ch_report_input_with_feature = ch_report_grouped
            .map { _key, payloads, feature_plots ->
                def payload = payloads.find { it != null }
                def featurePlot = feature_plots.find { it != null }
                [payload, featurePlot]
            }
            .filter { payload, featurePlot -> payload != null && featurePlot != null }
            .map { payload, featurePlot ->
                def (meta, best_ref_tsv, summary, depth_plot) = payload
                [meta, best_ref_tsv, summary, depth_plot, featurePlot]
            }

        ch_report_input_without_feature = ch_report_grouped
            .map { _key, payloads, feature_plots ->
                def payload = payloads.find { it != null }
                def featurePlot = feature_plots.find { it != null }
                [payload, featurePlot]
            }
            .filter { payload, featurePlot -> payload != null && featurePlot == null }
            .map { payload, _featurePlot ->
                def (meta, best_ref_tsv, summary, depth_plot) = payload
                [meta, best_ref_tsv, summary, depth_plot]
            }

        RENDER_HCV_REPORT (
            ch_report_input_with_feature,
            file("${projectDir}/assets/h2seq_logo.png"),
            workflow.manifest.version ?: ""
        )

        RENDER_SUMMARY_REPORT (
            ch_report_input_without_feature,
            file("${projectDir}/assets/h2seq_logo.png"),
            workflow.manifest.version ?: ""
        )

        ch_versions = ch_versions
            .mix(RENDER_HCV_REPORT.out.versions)
            .mix(RENDER_SUMMARY_REPORT.out.versions)
    } else {
        ch_report_input_without_feature = ch_report_base
            .map { _id, _long_reads, meta, best_ref_tsv, summary, depth_plot ->
                [meta, best_ref_tsv, summary, depth_plot]
            }

        RENDER_SUMMARY_REPORT (
            ch_report_input_without_feature,
            file("${projectDir}/assets/h2seq_logo.png"),
            workflow.manifest.version ?: ""
        )

        ch_versions = ch_versions.mix(RENDER_SUMMARY_REPORT.out.versions)
    }

    ch_summary_triggers = COVERAGE_METRICS.out.summary
        .mix(ch_split_consensuses)
        .mix(ch_hcv_reports)
        .mix(NANOQ.out.stats)
        .mix(FASTP.out.json)
        .mix(ch_long_raw_stats)
        .mix(ch_long_clean_stats)
        .mix(ch_short_raw_stats)
        .mix(ch_short_clean_stats)
        .mix(COUNT_MAPPED_READS.out.txt)

    ch_summary_triggers = ch_summary_triggers.mix(ch_best_ref_tsv)
    ch_summary_triggers = ch_summary_triggers
        .mix(PLOT_DEPTH_SUMMARY.out.depth)
        .mix(RENDER_SUMMARY_REPORT.out.pdf)
    if (params.virus_preset == "hcv" && params.run_hcv_glue) {
        ch_summary_triggers = ch_summary_triggers
            .mix(PARSE_HCV_GLUE_COVERAGE.out.tsv)
            .mix(PLOT_HCV_SUMMARY.out.feature)
            .mix(RENDER_HCV_REPORT.out.pdf)
    }

    BUILD_RUN_SUMMARY (
        ch_summary_triggers.collect(),
        file(params.outdir).toString(),
        workflow.manifest.version ?: ""
    )

    BUILD_MULTIQC_SECTIONS (
        ch_summary_triggers.collect(),
        file(params.outdir).toString()
    )

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

    ch_multiqc_files = ch_multiqc_files.mix(ch_long_raw_stats.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(ch_long_clean_stats.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(ch_short_raw_stats.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(ch_short_clean_stats.collect { it[1] })
    ch_multiqc_files = ch_multiqc_files.mix(BUILD_RUN_SUMMARY.out.mqc)
    ch_multiqc_files = ch_multiqc_files.mix(BUILD_MULTIQC_SECTIONS.out.coverage)
    ch_multiqc_files = ch_multiqc_files.mix(BUILD_MULTIQC_SECTIONS.out.read_stats)
    ch_multiqc_files = ch_multiqc_files.mix(BUILD_MULTIQC_SECTIONS.out.variants)
    if (params.virus_preset == "hcv" && params.run_hcv_glue) {
        ch_multiqc_files = ch_multiqc_files.mix(BUILD_MULTIQC_SECTIONS.out.region_coverage)
    }
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
    run_summary    = BUILD_RUN_SUMMARY.out.csv
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
