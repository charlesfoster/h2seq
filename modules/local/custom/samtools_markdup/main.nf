process SAMTOOLS_MARKDUP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.markdup.bam"), path("*.markdup.bam.csi"), emit: bam
    tuple val(meta), path("*.markdup_stats.txt"), emit: stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sortThreads = task.cpus > 1 ? task.cpus - 1 : 1
    """
    # Duplicates are only flagged (0x400), never removed: no -r anywhere.
    # markdup needs mate scores (fixmate -m) on a coordinate-sorted BAM.
    samtools sort -@ ${sortThreads} -n -u -T ${prefix}.namesort ${bam} \\
        | samtools fixmate -@ ${sortThreads} -m -u - - \\
        | samtools sort -@ ${sortThreads} -u -T ${prefix}.possort - \\
        | samtools markdup -@ ${sortThreads} -s -f ${prefix}.markdup_stats.txt - ${prefix}.markdup.bam

    samtools index -@ ${sortThreads} -c ${prefix}.markdup.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.markdup.bam
    touch ${prefix}.markdup.bam.csi
    touch ${prefix}.markdup_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: stub
    END_VERSIONS
    """
}
