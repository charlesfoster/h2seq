process CALCULATE_READ_STATS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.70--np112py36_0' :
        'biocontainers/biopython:1.70--np112py36_0'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path(reads), (env mean_length), (env std_dev) , emit: reads_and_stats
    path "versions.yml"                                                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    calculate_read_stats.py $reads
    mean_length=\$(head -n 1 mean_length.txt)
    std_dev=\$(head -n 1 std_dev.txt)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biopython: 1.70
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_filtered"
    """
    touch mean_length.txt
    touch std_dev.txt
    mean_length=123
    std_dev=456

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biopython: 1.70
    END_VERSIONS
    """
}
