process SPLIT_CONSENSUS_GENOMES {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.70--np112py36_0' :
        'biocontainers/biopython:1.70--np112py36_0'}"

    input:
    tuple val(meta), path(fasta), path(pattern)

    output:
    tuple val(meta), path("*.consensus_*.fa")  , emit: fastas

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    split_consensus_sequences.py --prefix ${prefix} $pattern $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biopython: 1.70
    END_VERSIONS
    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.consensus_reference_main.fa
    touch ${prefix}.consensus_reference_alt1.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biopython: 1.70
    END_VERSIONS
    """
}

