process ANNOTATE_VARIANTS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0' :
        'quay.io/biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(vcf_idx)

    output:
    tuple val(meta), path("*.annotated.tsv"), emit: tsv
    path "versions.yml"                           , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf "sample_id\\tread_type\\tchrom\\tpos\\tref\\talt\\tqual\\tdepth\\taf\\ttype\\tgt\\n" > ${prefix}.annotated.tsv
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%QUAL\\t%INFO/DP\\t%INFO/AF\\t%INFO/TYPE[\\t%GT]\\n' $vcf | \\
    awk -v OFS="\\t" -v sample="${meta.id}" -v read_type="${meta.long_reads ? 'long' : 'short'}" '{print sample, read_type, \$0}' >> ${prefix}.annotated.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n 1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    printf "sample_id\\tread_type\\tchrom\\tpos\\tref\\talt\\tqual\\tdepth\\taf\\ttype\\tgt\\n" > ${prefix}.annotated.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: stub
    END_VERSIONS
    """
}
