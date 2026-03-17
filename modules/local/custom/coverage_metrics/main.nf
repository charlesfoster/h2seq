process COVERAGE_METRICS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1' :
        'biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta), path(genome_regions_bed_gz), path(per_base_bed_gz)

    output:
    tuple val(meta), path("*.coverage_summary.tsv"), emit: summary
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def readType = meta.long_reads ? 'long' : 'short'
    """
    python3 ${projectDir}/bin/summarise_mosdepth_coverage.py \\
        --sample-id ${meta.id} \\
        --read-type ${readType} \\
        --reference-name ${meta.reference_name} \\
        --genome-bed-gz ${genome_regions_bed_gz} \\
        --per-base-bed-gz ${per_base_bed_gz} \\
        --min-depth ${params.consensus_min_depth} \\
        --summary-output ${prefix}.coverage_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def readType = meta.long_reads ? 'long' : 'short'
    """
    cat <<-EOF > ${prefix}.coverage_summary.tsv
    sample_id\tread_type\treference_name\treference_length\tpositions_covered\tgenome_coverage_pct\tmean_depth
    ${meta.id}\t${readType}\t${meta.reference_name ?: 'ref'}\t100\t90\t90\t12
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
