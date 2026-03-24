process LOFREQ_CALL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container {
        if (workflow.stubRun) {
            null
        } else if (workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container) {
            'https://depot.galaxyproject.org/singularity/lofreq:2.1.5--py39h9f2253c_15'
        } else {
            'quay.io/biocontainers/lofreq:2.1.5--py39h9f2253c_15'
        }
    }

    input:
    tuple val(meta), path(bam), path(bai), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.lofreq.raw.vcf.gz"), path("*.lofreq.raw.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                                             , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    lofreq call-parallel \\
        --no-baq \\
        --call-indels \\
        --pp-threads ${task.cpus} \\
        -f $fasta \\
        -o ${prefix}.lofreq.raw.vcf \\
        $bam

    bgzip -c ${prefix}.lofreq.raw.vcf > ${prefix}.lofreq.raw.vcf.gz
    tabix -f -p vcf ${prefix}.lofreq.raw.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq: \$(lofreq version 2>&1 | head -n 1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.lofreq.raw.vcf.gz
    touch ${prefix}.lofreq.raw.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq: stub
    END_VERSIONS
    """
}
