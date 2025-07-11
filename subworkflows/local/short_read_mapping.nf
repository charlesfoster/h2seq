include { BWA_INDEX as BWA_INDEX_PRIMERS          } from '../../modules/nf-core/bwa/index/main'
include { BWA_INDEX                               } from '../../modules/nf-core/bwa/index/main'
include { BWA_MEM as MAP_PRIMERS                  } from '../../modules/nf-core/bwa/mem/main'
include { BWA_MEM                                 } from '../../modules/nf-core/bwa/mem/main'
include { BEDTOOLS_BAMTOBED                       } from '../../modules/nf-core/bedtools/bamtobed/main'
include { SAMTOOLS_AMPLICONCLIP                   } from '../../modules/nf-core/samtools/ampliconclip/main'
include { SAMTOOLS_SORT                           } from '../../modules/nf-core/samtools/sort/main'

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

    // set up the sorting fasta
    ch_fasta_for_samtools_sort = ch_best_ref_fasta

    // index the fasta with BWA
    ch_bwa_index_input = ch_best_ref_fasta
        .map{ id, meta, fasta ->
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
        .combine(ch_index_for_combining, by:0)
        .map{ id, meta1, reads, meta2, fasta, index ->
            [meta1, reads, meta2, fasta, index]
        }

    BWA_MEM (
        ch_alignment_input,
        true
    )

    ch_mapped_bam = BWA_MEM.out.bam
    ch_mapped_bam_idx = BWA_MEM.out.csi
    ch_versions = ch_versions.mix(BWA_MEM.out.versions)

    // clip primers if needed
    if (!params.skip_primer_trimming){
        if (!params.primer_bed){
            // reshape the channel
            ch_bwa_input_fasta = ch_best_ref_fasta
                .map{id,meta,fasta ->
                    [meta, fasta]
                }

            // index the fasta
            BWA_INDEX_PRIMERS (
                ch_bwa_input_fasta
            )
            ch_versions = ch_versions.mix(BWA_INDEX.out.versions)

            ch_primer_fasta = Channel.fromPath(params.primer_fasta)
                .map { fasta ->
                    def meta = [id: 'primers']
                    return [meta, fasta]
                }

            // no need to match on a key here --> same primer fasta is recycled for all samples
            ch_bwa_mem_input = ch_primer_fasta
                .combine(BWA_INDEX.out.fasta_and_index)

            // map the primer fasta
            MAP_PRIMERS (
                ch_bwa_mem_input,
                true
            )
            ch_primer_bam = MAP_PRIMERS.out.bam
            ch_versions = ch_versions.mix(MAP_PRIMERS.out.versions)

            // convert to bed format
            BEDTOOLS_BAMTOBED (
                ch_primer_bam
            )

            ch_primer_bed = BEDTOOLS_BAMTOBED.out.bed
                .map { meta, bed ->
                    [ meta.id, meta, bed ]
                }

            ch_versions = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions)

            // prepare input for ampliconclip
            ch_samtools_ampliconclip_input = ch_mapped_bam
                .map { meta, bam ->
                    [meta.id, meta, bam]
                }
                .combine(ch_primer_bed, by:0)
                .map { id, meta, bam, meta2, bed ->
                    [meta, bam, bed]
                }
        } else {
            ch_primer_bed = Channel.fromPath(params.primer_bed, checkIfExists: true)

            ch_samtools_ampliconclip_input = ch_mapped_bam
                .combine(ch_primer_bed)
        }

        // clip the primers
        SAMTOOLS_AMPLICONCLIP (
            ch_samtools_ampliconclip_input,
            true,
            false
        )

        ch_clipped_bam = SAMTOOLS_AMPLICONCLIP.out.bam
        ch_versions = ch_versions.mix(SAMTOOLS_AMPLICONCLIP.out.versions)

        // sort the clipped bam
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
        ch_consensus_bam = ch_mapped_bam
    }

    emit:
        consensus_bam = ch_consensus_bam
        versions = ch_versions
}
