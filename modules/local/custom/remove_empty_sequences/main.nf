process REMOVE_EMPTY_SEQUENCES {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.70--np112py36_0' :
        'biocontainers/biopython:1.70--np112py36_0'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.reporting.fasta") , emit: fasta, optional: true
    tuple val(meta), path("*.empty.fasta") , emit: empty_fasta, optional: true
    path "versions.yml"                                                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    filter_fasta.py $fasta ${prefix}

    [ -s ${prefix}.reporting.fasta ] || rm ${prefix}.reporting.fasta
    [ -s ${prefix}.empty.fasta ] || rm ${prefix}.empty.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biopython: 1.70
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_filtered"
    """
    touch sample.reporting.fasta
    touch sample.reporting.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biopython: 1.70
    END_VERSIONS
    """
}
