process BUILD_RUN_SUMMARY {
    tag "run_summary"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'biocontainers/pandas:2.2.1'}"

    input:
    val trigger_files
    val outdir
    val pipeline_version

    output:
    path "combined_results_summary.csv", emit: csv
    path "combined_results_summary_mqc.json", emit: mqc

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 ${projectDir}/bin/build_run_summary.py \\
        --outdir "${outdir}" \\
        --pipeline-version "${pipeline_version}" \\
        --output combined_results_summary.csv \\
        --multiqc-output combined_results_summary_mqc.json
    """

    stub:
    """
    cat <<-EOF > combined_results_summary.csv
    sample_id,read_type,pipeline_version
    stub,stub,stub
    EOF

    cat <<-EOF > combined_results_summary_mqc.json
    {}
    EOF
    """
}
