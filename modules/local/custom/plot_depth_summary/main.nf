process PLOT_DEPTH_SUMMARY {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/longqc:1.2.0c--hdfd78af_0' :
        'quay.io/biocontainers/longqc:1.2.0c--hdfd78af_0' }"

    input:
    tuple val(meta), path(per_base_bed_gz), path(coverage_summary)

    output:
    tuple val(meta), path("*.depth_plot.png"), emit: depth
    path "versions.yml"                            , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.reference_name}.${meta.component_role}"
    def readType = meta.long_reads ? 'long' : 'short'
    def plotScriptMtime = file("${projectDir}/bin/plot_depth_summary.py").lastModified()
    """
    # plot_depth_summary_py_mtime=${plotScriptMtime}
    export MPLCONFIGDIR=\$PWD/.mplconfig

    python3 ${projectDir}/bin/plot_depth_summary.py \\
        --sample-id ${meta.id} \\
        --read-type ${readType} \\
        --reference-name "${meta.reference_name ?: ''}" \\
        --per-base-bed-gz ${per_base_bed_gz} \\
        --coverage-summary ${coverage_summary} \\
        --min-depth ${params.consensus_min_depth} \\
        --depth-plot-output ${prefix}.depth_plot.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.${meta.reference_name}.${meta.component_role}"
    def plotScriptMtime = file("${projectDir}/bin/plot_depth_summary.py").lastModified()
    """
    # plot_depth_summary_py_mtime=${plotScriptMtime}
    python3 - <<-'PY'
    import base64
    png = base64.b64decode("iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAusB9pM8S2QAAAAASUVORK5CYII=")
    open("${prefix}.depth_plot.png", "wb").write(png)
    PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
