process PREPARE_LOFREQ_VCF {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0' :
        'quay.io/biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(vcf_idx)

    output:
    tuple val(meta), path("*.prepared.vcf.gz"), path("*.prepared.vcf.gz.csi"), emit: vcf
    path "versions.yml"                                                       , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    python3 ${projectDir}/bin/prepare_lofreq_vcf.py \\
        --input $vcf \\
        --output ${prefix}.prepared.vcf \\
        --sample ${meta.id}

    bgzip -f ${prefix}.prepared.vcf
    bcftools index -f ${prefix}.prepared.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n 1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.prepared.vcf.gz
    touch ${prefix}.prepared.vcf.gz.csi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: stub
    END_VERSIONS
    """
}
