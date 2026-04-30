process LOFREQ_INDELQUAL {
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
    tuple val(meta), path(bam), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.indelqual.bam"), path("*.indelqual.bam.bai"), emit: bam
    path "versions.yml"                                                     , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    lofreq indelqual --dindel $bam -f $fasta -o ${prefix}.indelqual.bam
    samtools index ${prefix}.indelqual.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq: \$(lofreq version 2>&1 | head -n 1 | awk '{print \$2}')
        samtools: \$(samtools --version 2>&1 | head -n 1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.indelqual.bam
    touch ${prefix}.indelqual.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lofreq: stub
        samtools: stub
    END_VERSIONS
    """
}
