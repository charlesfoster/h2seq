process CREATE_PATTERN_FILE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(reference_fasta)

    output:
    tuple val(meta), path("*.pattern.txt"), emit: pattern

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Extract the first sequence header from the reference FASTA (without the >)
    head -n 1 ${reference_fasta} | sed 's/^>//' > ${prefix}.pattern.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n 1 | cut -d' ' -f4)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "dummy_pattern" > ${prefix}.pattern.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n 1 | cut -d' ' -f4)
    END_VERSIONS
    """
}
