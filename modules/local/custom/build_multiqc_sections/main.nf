process BUILD_MULTIQC_SECTIONS {
    tag "multiqc_sections"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'biocontainers/pandas:2.2.1'}"

    input:
    val trigger_files
    val outdir

    output:
    path "coverage_statistics_mqc.json", emit: coverage
    path "read_statistics_mqc.json", emit: read_stats
    path "variant_calling_mqc.json", emit: variants

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 ${projectDir}/bin/build_multiqc_sections.py \\
        --outdir "${outdir}" \\
        --coverage-output coverage_statistics_mqc.json \\
        --read-output read_statistics_mqc.json \\
        --variant-output variant_calling_mqc.json
    """

    stub:
    """
    echo '{}' > coverage_statistics_mqc.json
    echo '{}' > read_statistics_mqc.json
    echo '{}' > variant_calling_mqc.json
    """
}
