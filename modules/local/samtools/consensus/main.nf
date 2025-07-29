process SAMTOOLS_CONSENSUS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_0' :
        'biocontainers/samtools:1.20--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.consensus.fa")  , emit: fasta
    tuple val(meta), path("*.consensus.bam"), path("*.csi")  , emit: bam
    tuple val(meta), path("*.simple_consensus.fa")  , emit: simple_fasta
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ln -s $bam ./${prefix}.consensus.bam
    samtools index --csi ./${prefix}.consensus.bam

    samtools \\
        consensus \\
        --threads ${task.cpus-1} \\
        $args \\
        -o tmp.fa \\
        $bam

    samtools \\
        consensus \\
        --threads ${task.cpus-1} \\
        $args2 \\
        -o tmp2.fa \\
        $bam

    sed -e "/^>/s/>/>${prefix} /g" -e "s/\\*/-/g" tmp.fa > ${prefix}.consensus.fa
    sed -e "/^>/s/>/>${prefix} /g" -e "s/\\*/-/g" tmp2.fa > ${prefix}.simple_consensus.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.consensus.fa
    touch ${prefix}.simple_consensus.fa
    touch ${prefix}.consensus.bam
    touch ${prefix}.consensus.bam.csi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
