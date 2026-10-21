process BUILD_RUN_SUMMARY {
    tag "run_summary"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'biocontainers/pandas:2.2.1'}"

    input:
    tuple val(rels), path(files, stageAs: 'staged_inputs/?/*')
    val pipeline_version

    output:
    path "combined_results_summary.csv", emit: csv
    path "reference_component_summary.csv", emit: components
    path "combined_results_summary_mqc.json", emit: mqc

    when:
    task.ext.when == null || task.ext.when

    script:
    def rel_words = rels.collect { rel -> "'" + rel.toString().replace("'", "'\\''") + "'" }.join(' ')
    """
    # Rebuild <sample>/<read_type>/<file> from the staged inputs (staged_inputs/N is the Nth entry of rels)
    mkdir -p staged_outdir
    n=0
    for rel in ${rel_words}; do
        n=\$((n + 1))
        mkdir -p "staged_outdir/\$(dirname "\$rel")"
        ln -sf "\$PWD/staged_inputs/\$n/\$(basename "\$rel")" "staged_outdir/\$rel"
    done

    python3 ${projectDir}/bin/build_run_summary.py \\
        --outdir staged_outdir \\
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
