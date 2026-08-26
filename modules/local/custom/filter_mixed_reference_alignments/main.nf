process FILTER_MIXED_REFERENCE_ALIGNMENTS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.mixed_filtered.bam"), path("*.mixed_filtered.bam.csi"), emit: bam
    tuple val(meta), path("*.mixed_assignment.tsv"), emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def readType = meta.long_reads ? 'long' : 'short'
    def sortThreads = task.cpus > 1 ? task.cpus - 1 : 1
    """
    samtools collate -@ ${sortThreads} -O -u ${bam} \\
        | samtools view -h - \\
        | awk \\
            -v sample_id="${meta.id}" \\
            -v read_type="${readType}" \\
            -v min_mapq="${params.mixed_assignment_min_mapq}" \\
            -v summary_output="${prefix}.mixed_assignment.tsv" \\
            -f ${projectDir}/bin/filter_mixed_reference_alignments.awk \\
        | samtools view -b - \\
        | samtools sort -@ ${sortThreads} -o ${prefix}.mixed_filtered.bam -

    samtools index -@ ${sortThreads} -c ${prefix}.mixed_filtered.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def readType = meta.long_reads ? 'long' : 'short'
    """
    touch ${prefix}.mixed_filtered.bam
    touch ${prefix}.mixed_filtered.bam.csi
    cat <<-EOF > ${prefix}.mixed_assignment.tsv
    sample_id\tread_type\tassignment\treference_name\tfragments\tprimary_alignment_records\tmin_mapq\treference_count
    ${meta.id}\t${readType}\tassigned\tref\t1\t1\t${params.mixed_assignment_min_mapq}\t1
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: stub
    END_VERSIONS
    """
}
