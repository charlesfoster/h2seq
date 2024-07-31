process KALLISTO_QUANT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kallisto:0.48.0--h15996b6_2':
        'biocontainers/kallisto:0.48.0--h15996b6_2' }"

    input:
    tuple val(meta), path(reads), val(fragment_length), val(fragment_length_sd), val(meta2), path(index)

    output:
    tuple val(meta), path("${prefix}.abundance.tsv")        , emit: tsv
    tuple val(meta), path("*.run_info.json")  , emit: json_info
    tuple val(meta), path("*.kallisto_quant.log")            , emit: log
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def single_end_params = ''
    if (meta.single_end) {
        if (!(fragment_length =~ /^\d+$/)) {
            error "fragment_length must be set and numeric for single-end data"
        }
        if (!(fragment_length_sd =~ /^\d+$/)) {
            error "fragment_length_sd must be set and numeric for single-end data"
        }
        single_end_params = "--single --single-overhang --fragment-length=${fragment_length} --sd=${fragment_length_sd}"
    }

    """
    mkdir -p $prefix && kallisto quant \\
            --threads ${task.cpus} \\
            --index ${index} \\
            ${single_end_params} \\
            ${args} \\
            -o $prefix \\
            ${reads} 2> >(tee -a ${prefix}.kallisto_quant.log >&2)

    mv ${prefix}/run_info.json ${prefix}.run_info.json
    mv ${prefix}/abundance.tsv ${prefix}.abundance.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallisto: \$(echo \$(kallisto version) | sed "s/kallisto, version //g" )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p $prefix
    touch ${prefix}.abundance.tsv
    touch ${prefix}.kallisto_quant.log
    touch ${prefix}.run_info.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallisto: \$(echo \$(kallisto version) | sed "s/kallisto, version //g" )
    END_VERSIONS
    """
}
