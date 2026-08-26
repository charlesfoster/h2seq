process PARSE_HCV_GLUE_COVERAGE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta), path(report_html)

    output:
    tuple val(meta), path("*.hcv_glue_coverage.tsv"), emit: tsv
    path "versions.yml"                             , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.reference_name}.${meta.component_role}"
    def readType = meta.long_reads ? 'long' : 'short'
    """
    python3 ${projectDir}/bin/parse_hcv_glue_coverage.py \\
        --sample-id ${meta.id} \\
        --read-type ${readType} \\
        --reference-name "${meta.reference_name}" \\
        --input ${report_html} \\
        --output ${prefix}.hcv_glue_coverage.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.reference_name}.${meta.component_role}"
    def readType = meta.long_reads ? 'long' : 'short'
    """
    cat <<-EOF > ${prefix}.hcv_glue_coverage.tsv
    sample_id\tread_type\treference_name\tfeature\tcoverage_pct
    ${meta.id}\t${readType}\t${meta.reference_name}\tPolyprotein\t99.0
    ${meta.id}\t${readType}\t${meta.reference_name}\tCore\t100.0
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
