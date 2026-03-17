process BUILD_MULTIQC_SECTIONS {
    tag "multiqc_sections"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/longqc:1.2.0c--0' :
        'quay.io/biocontainers/longqc:1.2.0c--hdfd78af_0' }"

    input:
    val trigger_files
    val outdir

    output:
    path "coverage_statistics_mqc.json", emit: coverage
    path "region_coverage_mqc.png", emit: region_coverage
    path "coverage_per_genomic_region.pdf", emit: region_coverage_pdf
    path "read_statistics_mqc.json", emit: read_stats
    path "variant_calling_mqc.json", emit: variants

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export MPLCONFIGDIR=\$PWD/.mplconfig
    export XDG_CACHE_HOME=\$PWD/.cache

    python3 ${projectDir}/bin/build_multiqc_sections.py \\
        --outdir "${outdir}" \\
        --coverage-output coverage_statistics_mqc.json \\
        --region-coverage-output region_coverage_mqc.png \\
        --region-coverage-pdf-output coverage_per_genomic_region.pdf \\
        --read-output read_statistics_mqc.json \\
        --variant-output variant_calling_mqc.json
    """

    stub:
    """
    echo '{}' > coverage_statistics_mqc.json
    python3 - <<-'PY'
    import base64
    png = base64.b64decode("iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAusB9pM8S2QAAAAASUVORK5CYII=")
    open("region_coverage_mqc.png", "wb").write(png)
    open("coverage_per_genomic_region.pdf", "wb").write(b"%PDF-1.1\\n%%EOF\\n")
    PY
    echo '{}' > read_statistics_mqc.json
    echo '{}' > variant_calling_mqc.json
    """
}
