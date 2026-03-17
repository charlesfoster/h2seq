include { MINIMAP2_ALIGN                          } from '../../modules/nf-core/minimap2/align/main'
include { CLAIR3                                  } from '../../modules/local/custom/clair3/main'
include { SAMTOOLS_FAIDX                          } from '../../modules/nf-core/samtools/faidx/main'
include { PREPARE_CONSENSUS_BAM                   } from './prepare_consensus_bam'

workflow LONG_READ_MAPPING {
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

    // combine with the best ref
    ch_alignment_input = ch_clean_reads_for_alignment
        .combine(ch_best_ref_fasta, by: 0)
        .map{ _id, meta1, reads, meta2, fasta ->
            [meta1, reads, meta2, fasta]
        }

    // align the reads
    MINIMAP2_ALIGN (
        ch_alignment_input,
        true,
        true,
        "bai",
        false,
        false
    )

    ch_mapped_bam = MINIMAP2_ALIGN.out.bam
    ch_mapped_bam_idx = MINIMAP2_ALIGN.out.index
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)

    PREPARE_CONSENSUS_BAM (
        ch_best_ref_fasta,
        ch_mapped_bam
    )

    ch_consensus_bam = PREPARE_CONSENSUS_BAM.out.consensus_bam
    ch_consensus_bam_idx = PREPARE_CONSENSUS_BAM.out.consensus_bam_idx
    ch_amplicon_bed = PREPARE_CONSENSUS_BAM.out.amplicon_bed
    ch_versions = ch_versions.mix(PREPARE_CONSENSUS_BAM.out.versions)

    // index best_ref_fasta
    SAMTOOLS_FAIDX (
        ch_best_ref_fasta.map { _id, meta, fasta -> [meta, fasta] }
    )

    ch_best_ref_fasta_idx = SAMTOOLS_FAIDX.out.fa_and_idx
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    // run clair3
    ch_clair3_input = ch_consensus_bam
        .map {meta, bam ->
            [meta.id, meta, bam]
        }
        .combine(ch_consensus_bam_idx.map { meta, idx -> [meta.id, idx] }, by: 0)
        .combine(ch_best_ref_fasta_idx.map { meta, fa, fai -> [meta.id, fa, fai] }, by: 0)
        .map {_id, meta, bam, bam_idx, fa, fai ->
            [meta, bam, bam_idx, fa, fai]
        }

    CLAIR3 (
        ch_clair3_input
    )

    ch_clair3_vcf = CLAIR3.out.vcf
    ch_versions = ch_versions.mix(CLAIR3.out.versions)

    emit:
    consensus_bam = ch_consensus_bam
    consensus_bam_idx = ch_consensus_bam_idx
    amplicon_bed = ch_amplicon_bed
    versions = ch_versions
}
