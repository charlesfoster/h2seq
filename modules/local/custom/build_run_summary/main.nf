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
    path "reference_component_summary.csv", emit: components
    path "combined_results_summary_mqc.json", emit: mqc

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 ${projectDir}/bin/build_run_summary.py \\
        --outdir "${outdir}" \\
        --pipeline-version "${pipeline_version}" \\
        --min-reference-coverage-pct ${params.qc_min_ref_coverage_pct} \\
        --output combined_results_summary.csv \\
        --component-output reference_component_summary.csv \\
        --multiqc-output combined_results_summary_mqc.json
    """

    stub:
    """
    cat <<-EOF > combined_results_summary.csv
    sample_id,read_type,pipeline_version
    stub,stub,stub
    EOF

    cat <<-EOF > reference_component_summary.csv
    sample_id,read_type,component_role,reference_name,assigned_fraction_of_assigned
    stub,stub,main,ref,1.0
    EOF

    cat <<-EOF > combined_results_summary_mqc.json
    {}
    EOF
    """
}
