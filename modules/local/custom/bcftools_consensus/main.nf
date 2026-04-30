process BCFTOOLS_CONSENSUS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0' :
        'quay.io/biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(vcf_idx), path(fasta), path(mask_bed)

    output:
    tuple val(meta), path("*.consensus.fa"), emit: fasta
    path "versions.yml"                     , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def haplotypeArg = params.majority_allele_consensus ? '' : "        -H I \\\n"
    def headerLabel = meta.reference_name ? "${prefix} ${meta.reference_name}" : prefix
    """
    bcftools consensus \\
        -f $fasta \\
        -m $mask_bed \\
        --mark-del '-' \\
${haplotypeArg}        \
        $vcf | sed "/^>/s/.*/>${headerLabel}/" > ${prefix}.consensus.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n 1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf ">${prefix}\\nNNNN\\n" > ${prefix}.consensus.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: stub
    END_VERSIONS
    """
}
