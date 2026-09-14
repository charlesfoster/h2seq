process BCFTOOLS_CONSENSUS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0' :
        'quay.io/biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(vcf_idx), path(fasta), path(mask_bed), path(simple_vcf), path(simple_vcf_idx)

    output:
    tuple val(meta), path("*.consensus.fa")              , emit: fasta
    tuple val(meta), path("*.consensus_simple.fa")       , emit: simple_fasta
    path "versions.yml"                                   , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bcftools consensus \\
        -f $fasta \\
        -m $mask_bed \\
        --mark-del '-' \\
        -H I \\
        $vcf | awk -v p="${prefix}" '/^>/{print ">" p " " substr(\$0,2); next}{print}' > ${prefix}.consensus.fa

    bcftools consensus \\
        -f $fasta \\
        -m $mask_bed \\
        --mark-del '-' \\
        -H 1 \\
        $simple_vcf | awk -v p="${prefix}" '/^>/{print ">" p " " substr(\$0,2); next}{print}' > ${prefix}.consensus_simple.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n 1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf ">${prefix}\\nNNNN\\n" > ${prefix}.consensus.fa
    printf ">${prefix}\\nNNNN\\n" > ${prefix}.consensus_simple.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: stub
    END_VERSIONS
    """
}
