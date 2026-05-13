process MINIMAP2_COMPETITIVE_MAP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' :
        'biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' }"

    input:
    tuple val(meta), path(reads), val(reference_meta), path(reference_fasta), path(reference_index)

    output:
    tuple val(meta), path("*.competitive.bam"), path("*.competitive.bam.bai"), path("*.competitive.primary.sam.gz"), val(reference_meta), path(reference_fasta), emit: bam
    path "versions.yml"                                                                                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sort_threads = task.cpus > 1 ? task.cpus - 1 : 1
    """
    minimap2 \\
        $args \\
        -a \\
        -t $task.cpus \\
        ${reference_index} \\
        ${reads} \\
        | samtools sort -@ ${sort_threads} -o ${prefix}.competitive.bam -

    samtools index -@ ${sort_threads} ${prefix}.competitive.bam
    samtools view -@ ${sort_threads} -h -F 2308 ${prefix}.competitive.bam | gzip -c > ${prefix}.competitive.primary.sam.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.competitive.bam
    touch ${prefix}.competitive.bam.bai
    gzip -c /dev/null > ${prefix}.competitive.primary.sam.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: stub
        samtools: stub
    END_VERSIONS
    """
}
