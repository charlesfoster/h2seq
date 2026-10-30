process SAMTOOLS_CONSENSUS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.24--h9dcdb79_1' :
        'biocontainers/samtools:1.24--h9dcdb79_1' }"

    input:
    tuple val(meta), path(bam), path(reference)

    output:
    tuple val(meta), path("*.draft_consensus.fa"), emit: fasta
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // -T/--ref-qual 1: fill zero-coverage positions from the mapping reference instead of N,
    // so round 2 can map through divergent regions (e.g. HCV HVR1) rather than an N run.
    // -c 0/--call-fract 0: always call the most common base. The default (0.75) writes N wherever
    // the top base is under 75% of the depth, and that N survives into the final consensus because
    // round 2 only applies variants onto the draft. Minor alleles are instead resolved in round 2
    // by variant calling and the --consensus_call_af / IUPAC thresholds.
    """
    samtools consensus \\
        -m simple \\
        -f fasta \\
        -a \\
        -T $reference \\
        --ref-qual 1 \\
        -c 0 \\
        $args \\
        $bam > ${prefix}.draft_consensus.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.draft_consensus.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
