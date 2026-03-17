process PLOT_HCV_SUMMARY {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/longqc:1.2.0c--0' :
        'quay.io/biocontainers/longqc:1.2.0c--hdfd78af_0' }"

    input:
    tuple val(meta), path(per_base_bed_gz), path(coverage_summary), path(hcv_coverage_tsv)

    output:
    tuple val(meta), path("*.depth_plot.png"), emit: depth
    tuple val(meta), path("*.feature_coverage.png"), emit: feature
    path "versions.yml"                               , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def readType = meta.long_reads ? 'long' : 'short'
    def plotScriptMtime = file("${projectDir}/bin/plot_hcv_summary.py").lastModified()
    """
    # plot_hcv_summary_py_mtime=${plotScriptMtime}
    export MPLCONFIGDIR=\$PWD/.mplconfig

    python3 ${projectDir}/bin/plot_hcv_summary.py \\
        --sample-id ${meta.id} \\
        --read-type ${readType} \\
        --reference-name "${meta.reference_name ?: ''}" \\
        --per-base-bed-gz ${per_base_bed_gz} \\
        --coverage-summary ${coverage_summary} \\
        --hcv-coverage ${hcv_coverage_tsv} \\
        --min-depth ${params.consensus_min_depth} \\
        --depth-plot-output ${prefix}.depth_plot.png \\
        --feature-plot-output ${prefix}.feature_coverage.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def plotScriptMtime = file("${projectDir}/bin/plot_hcv_summary.py").lastModified()
    """
    # plot_hcv_summary_py_mtime=${plotScriptMtime}
    python3 - <<-'PY'
    import base64
    png = base64.b64decode("iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAusB9pM8S2QAAAAASUVORK5CYII=")
    open("${prefix}.depth_plot.png", "wb").write(png)
    open("${prefix}.feature_coverage.png", "wb").write(png)
    PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
