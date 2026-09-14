process RENDER_HCV_REPORT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/longqc:1.2.0c--hdfd78af_0' :
        'quay.io/biocontainers/longqc:1.2.0c--hdfd78af_0' }"

    input:
    tuple val(meta), path(best_reference_tsv), path(coverage_summary), path(depth_plot), path(assignment_summary), path(feature_plot), path(hcv_coverage)
    path logo
    val pipeline_version

    output:
    tuple val(meta), path("*.summary_report.pdf"), emit: pdf
    path "versions.yml"                              , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.reference_name}.${meta.component_role}"
    def readType = meta.long_reads ? 'long' : 'short'
    def snvMinAf = meta.long_reads ? params.ont_min_snv_af : params.illumina_min_snv_af
    def indelMinAf = meta.long_reads ? params.ont_min_indel_af : params.illumina_min_indel_af
    def renderScriptMtime = file("${projectDir}/bin/render_hcv_report.py").lastModified()
    """
    # render_hcv_report_py_mtime=${renderScriptMtime}
    export MPLCONFIGDIR=\$PWD/.mplconfig

    python3 ${projectDir}/bin/render_hcv_report.py \\
        --sample-id ${meta.id} \\
        --read-type ${readType} \\
        --best-reference-tsv ${best_reference_tsv} \\
        --reference-name "${meta.reference_name}" \\
        --component-role "${meta.component_role}" \\
        --coverage-summary ${coverage_summary} \\
        --assignment-summary ${assignment_summary} \\
        --depth-plot ${depth_plot} \\
        --feature-plot ${feature_plot} \\
        --hcv-coverage ${hcv_coverage} \\
        --include-feature-plot \\
        --logo ${logo} \\
        --pipeline-version "${pipeline_version}" \\
        --reference-selection-tool "${params.reference_selection_tool}" \\
        --consensus-min-depth "${params.consensus_min_depth}" \\
        --snv-min-af "${snvMinAf}" \\
        --indel-min-af "${indelMinAf}" \\
        --long-reads-min-len "${params.long_reads_min_len}" \\
        --long-reads-max-len "${params.long_reads_max_len}" \\
        --short-reads-min-len "${params.short_reads_min_len}" \\
        --output ${prefix}.summary_report.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.reference_name}.${meta.component_role}"
    def renderScriptMtime = file("${projectDir}/bin/render_hcv_report.py").lastModified()
    """
    # render_hcv_report_py_mtime=${renderScriptMtime}
    touch ${prefix}.summary_report.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
