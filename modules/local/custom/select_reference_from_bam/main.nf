process SELECT_REFERENCE_FROM_BAM {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'biocontainers/pandas:2.2.1'}"

    input:
    tuple val(meta), path(competitive_bam), path(competitive_bai), path(primary_sam), val(reference_meta), path(reference_fasta)

    output:
    tuple val(meta), path("*.best_reference.tsv")              , emit: best_ref_tsv
    tuple val(meta), path("*.best_reference.txt")              , emit: best_ref_txt
    tuple val(meta), path("*.alternate_subtypes.txt")          , emit: alt_ref_txt
    tuple val(meta), path("*.reference_selection_ranking.tsv") , emit: ranking_tsv
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def read_type = meta.long_reads ? "long" : "short"
    """
    python3 ${projectDir}/bin/select_reference_from_bam.py \\
        --primary-sam ${primary_sam} \\
        --bam ${competitive_bam} \\
        --panel-fasta ${reference_fasta} \\
        --sample-name ${prefix} \\
        --read-type ${read_type} \\
        --output ${prefix}.best_reference.tsv \\
        --best-ref-txt ${prefix}.best_reference.txt \\
        --alternate-subtype-txt ${prefix}.alternate_subtypes.txt \\
        --ranking-output ${prefix}.reference_selection_ranking.tsv \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.best_reference.tsv
    touch ${prefix}.best_reference.txt
    touch ${prefix}.alternate_subtypes.txt
    touch ${prefix}.reference_selection_ranking.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
