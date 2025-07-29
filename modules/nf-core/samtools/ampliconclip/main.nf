process SAMTOOLS_AMPLICONCLIP {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.20--h50ea8bc_0' :
        'biocontainers/samtools:1.20--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam), path(bed)
    val save_cliprejects
    val save_clipstats

    output:
    tuple val(meta), path("*.clipallowed.bam")  , emit: bam
    tuple val(meta), path("*.clipstats.txt")    , optional:true, emit: stats
    tuple val(meta), path("*.cliprejects.bam")  , optional:true, emit: rejects_bam
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def rejects = save_cliprejects ? "--rejects-file ${prefix}.cliprejects.bam" : ""
    def stats   = save_clipstats   ? "-f ${prefix}.clipstats.txt"               : ""
    if ("$bam" == "${prefix}.clipallowed.bam") error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    samtools \\
        ampliconclip \\
        --threads ${task.cpus-1} \\
        $args \\
        $rejects \\
        $stats \\
        -b $bed \\
        -o ${prefix}.clipallowed.bam \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def rejects_file = save_cliprejects ? "${prefix}.cliprejects.bam" : ""
    def stats_file = save_clipstats ? "${prefix}.clipstats.txt" : ""
    """
    touch ${prefix}.clipallowed.bam
    if [[ "$rejects_file" != "" ]]; then
        touch $rejects_file
    fi
    if [[ "$stats_file" != "" ]]; then
        touch $stats_file
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
