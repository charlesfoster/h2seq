process GENERATE_WHOLE_GENOME_BED {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(reference_fai)

    output:
    tuple val(meta), path("*.whole_genome.bed"), emit: bed
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    awk 'NR==1 { print \$1 "\\t0\\t" \$2 "\\twhole_genome" }' ${reference_fai} > ${prefix}.whole_genome.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n 1 | cut -d' ' -f4)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat <<-EOF > ${prefix}.whole_genome.bed
    ref\t0\t100\twhole_genome
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n 1 | cut -d' ' -f4)
    END_VERSIONS
    """
}
