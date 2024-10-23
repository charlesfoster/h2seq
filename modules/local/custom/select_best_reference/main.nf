process SELECT_BEST_REFERENCE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'biocontainers/pandas:2.2.1'}"

    input:
    tuple val(meta), path(abundances)

    output:
    tuple val(meta), path("*.best_reference.tsv"), emit: best_ref_tsv
    tuple val(meta), path("*.best_reference.txt"), emit: best_ref_txt
    tuple val(meta), path("*.alternate_subtypes.txt"), emit: alt_ref_txt, optional: true
    path "versions.yml"                                                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    select_best_reference.py \\
        --input $abundances \\
        --sample_name $prefix \\
        --output ${prefix}.best_reference.tsv \\
        --best_ref_txt ${prefix}.best_reference.txt \\
        --alternate_subtype_txt ${prefix}.alternate_subtypes.txt \\
        $args 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pandas: 2.2.1
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.best_reference.tsv
    touch ${prefix}.best_reference.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pandas: 2.2.1
    END_VERSIONS
    """
}
