include { BWA_INDEX                               } from '../../modules/nf-core/bwa/index/main'
include { BWA_MEM                                 } from '../../modules/nf-core/bwa/mem/main'
include { PREPARE_CONSENSUS_BAM                   } from './prepare_consensus_bam'

workflow SHORT_READ_MAPPING {
    take:
        ch_best_ref_fasta
        ch_clean_reads

    main:
        ch_versions = Channel.empty()

    // format the reads for combining
    ch_clean_reads_for_alignment = ch_clean_reads
        .map{ meta, reads ->
            [meta.id, meta, reads]
        }

    // index the fasta with BWA
    ch_bwa_index_input = ch_best_ref_fasta
        .map{ _id, meta, fasta ->
            [meta, fasta]
        }

    BWA_INDEX (
        ch_bwa_index_input
    )

    ch_versions = ch_versions.mix(BWA_INDEX.out.versions)

    ch_index_for_combining = BWA_INDEX.out.fasta_and_index
        .map{ meta, fasta, index ->
            [meta.id, meta, fasta, index]
        }

    ch_alignment_input = ch_clean_reads_for_alignment
        .combine(ch_index_for_combining, by: 0)
        .map{ _id, meta1, reads, meta2, fasta, index ->
            [meta1, reads, meta2, fasta, index]
        }

    BWA_MEM (
        ch_alignment_input,
        true
    )

    ch_mapped_bam = BWA_MEM.out.bam
    _ch_mapped_bam_idx = BWA_MEM.out.csi
    ch_versions = ch_versions.mix(BWA_MEM.out.versions)

    PREPARE_CONSENSUS_BAM (
        ch_best_ref_fasta,
        ch_mapped_bam
    )

    ch_consensus_bam = PREPARE_CONSENSUS_BAM.out.consensus_bam
    ch_consensus_bam_idx = PREPARE_CONSENSUS_BAM.out.consensus_bam_idx
    ch_versions = ch_versions.mix(PREPARE_CONSENSUS_BAM.out.versions)

    emit:
        consensus_bam = ch_consensus_bam
        consensus_bam_idx = ch_consensus_bam_idx
        versions = ch_versions
}
