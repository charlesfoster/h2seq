include { BWA_INDEX as BWA_INDEX_PRIMERS } from '../../modules/nf-core/bwa/index/main'
include { BWA_MEM as MAP_PRIMERS } from '../../modules/nf-core/bwa/mem/main'
include { BEDTOOLS_BAMTOBED } from '../../modules/nf-core/bedtools/bamtobed/main'
include { SAMTOOLS_AMPLICONCLIP } from '../../modules/nf-core/samtools/ampliconclip/main'
include { SAMTOOLS_SORT } from '../../modules/nf-core/samtools/sort/main'

workflow PREPARE_CONSENSUS_BAM {
    take:
    ch_best_ref_fasta
    ch_mapped_bam

    main:
    ch_versions = Channel.empty()
    ch_fasta_for_samtools_sort = ch_best_ref_fasta
    ch_amplicon_bed = Channel.empty()

    if (!params.skip_primer_trimming) {
        if (!params.primer_bed) {
            ch_bwa_input_fasta = ch_best_ref_fasta
                .map { _id, meta, fasta ->
                    [meta, fasta]
                }

            BWA_INDEX_PRIMERS (
                ch_bwa_input_fasta
            )
            ch_versions = ch_versions.mix(BWA_INDEX_PRIMERS.out.versions)

            ch_primer_fasta = Channel.fromPath(params.primer_fasta)
                .map { fasta ->
                    def meta = [id: 'primers']
                    [meta, fasta]
                }

            ch_reference_with_index = ch_best_ref_fasta
                .map { _id, meta, fasta ->
                    [meta.id, meta, fasta]
                }
                .combine(
                    BWA_INDEX_PRIMERS.out.index
                        .map { meta, index ->
                            [meta.id, meta, index]
                        },
                    by: 0
                )
                .map { _id, meta, fasta, _meta2, index ->
                    [meta, fasta, index]
                }

            ch_bwa_mem_input = ch_primer_fasta
                .combine(ch_reference_with_index)
                .map { primer_meta, primer_fasta, ref_meta, ref_fasta, ref_index ->
                    [primer_meta, primer_fasta, ref_meta, ref_fasta, ref_index]
                }

            MAP_PRIMERS (
                ch_bwa_mem_input,
                true
            )
            ch_versions = ch_versions.mix(MAP_PRIMERS.out.versions)

            BEDTOOLS_BAMTOBED (
                MAP_PRIMERS.out.bam
            )
            ch_versions = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions)

            ch_primer_bed = BEDTOOLS_BAMTOBED.out.bed
                .map { meta, bed ->
                    [meta.id, meta, bed]
                }

            ch_amplicon_bed = ch_primer_bed
                .map { _id, meta, bed ->
                    [meta, bed]
                }

            ch_samtools_ampliconclip_input = ch_mapped_bam
                .map { meta, bam ->
                    [meta.id, meta, bam]
                }
                .combine(ch_primer_bed, by: 0)
                .map { _id, meta, bam, _meta2, bed ->
                    [meta, bam, bed]
                }
        } else {
            ch_primer_bed = Channel.fromPath(params.primer_bed, checkIfExists: true)

            ch_amplicon_bed = ch_mapped_bam
                .map { meta, bam ->
                    [meta, bam]
                }
                .combine(ch_primer_bed)
                .map { meta, _bam, bed ->
                    [meta, bed]
                }

            ch_samtools_ampliconclip_input = ch_mapped_bam
                .combine(ch_primer_bed)
        }

        SAMTOOLS_AMPLICONCLIP (
            ch_samtools_ampliconclip_input,
            true,
            false
        )

        ch_bam_for_sort = SAMTOOLS_AMPLICONCLIP.out.bam
        ch_versions = ch_versions.mix(SAMTOOLS_AMPLICONCLIP.out.versions)
    } else {
        ch_bam_for_sort = ch_mapped_bam
        ch_amplicon_bed = ch_mapped_bam
            .map { meta, bam ->
                [meta, file("${projectDir}/assets/reference_data/empty.bed")]
            }
    }

    ch_samtools_sort_input = ch_bam_for_sort
        .map { meta, bam ->
            [meta.id, meta, bam]
        }
        .combine(ch_fasta_for_samtools_sort, by: 0)
        .map { _id, meta, bam, _meta2, fasta ->
            [meta, bam, fasta]
        }

    SAMTOOLS_SORT (
        ch_samtools_sort_input
    )

    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    emit:
    consensus_bam = SAMTOOLS_SORT.out.bam
    consensus_bam_idx = SAMTOOLS_SORT.out.csi
    amplicon_bed = ch_amplicon_bed
    versions = ch_versions
}
