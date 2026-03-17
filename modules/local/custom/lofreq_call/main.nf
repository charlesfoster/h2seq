process LOFREQ_CALL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://quay.io/biocontainers/lofreq:2.1.5--py39h917a906_8' :
        'quay.io/biocontainers/lofreq:2.1.5--py39h917a906_8' }"

    input:
    tuple val(meta), path(bam), path(bai), path(fasta)

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
