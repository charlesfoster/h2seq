process FILTER_VARIANTS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0' :
        'quay.io/biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(vcf_idx)

    output:
    tuple val(meta), path("*.final.vcf.gz"), path("*.final.vcf.gz.csi")                 , emit: vcf
    tuple val(meta), path("*.final_consensus.vcf.gz"), path("*.final_consensus.vcf.gz.csi"), emit: consensus_vcf
    tuple val(meta), path("*.final_simple.vcf.gz"), path("*.final_simple.vcf.gz.csi")   , emit: simple_vcf
    path "versions.yml"                                                                 , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def minQual = meta.long_reads ? params.ont_variant_min_qual : params.illumina_variant_min_qual
    def minSnvAf = meta.long_reads ? params.ont_min_snv_af : params.illumina_min_snv_af
    def minIndelAf = meta.long_reads ? params.ont_min_indel_af : params.illumina_min_indel_af
    def iupacMinAf = meta.long_reads ? params.ont_iupac_min_af : params.illumina_iupac_min_af
    def iupacMinAltReads = meta.long_reads ? params.ont_iupac_min_alt_reads : params.illumina_iupac_min_alt_reads
    // Mixed (GT 0/1) SNP sites that are too weakly supported are kept in final.vcf.gz (reported) but are not
    // IUPAC-encoded in the consensus. Alt reads are estimated as DP*AF; the 0.5 offset guards against float error.
    def iupacMinAltReadsNum = iupacMinAltReads as Double
    def weakMixedCond = "INFO/AF < ${iupacMinAf}"
    if (iupacMinAltReadsNum > 0) {
        weakMixedCond += " || INFO/DP*INFO/AF < ${iupacMinAltReadsNum - 0.5}"
    }
    """
    bcftools index -f $vcf
    bcftools +fill-tags $vcf -Ou -- -t TYPE | \\
        bcftools norm -Ou -a -m - | \\
        bcftools view -f 'PASS,dn,dp,.' -i '((TYPE="snp" && INFO/AF >= ${minSnvAf}) || (TYPE="indel" && INFO/AF >= ${minIndelAf})) && INFO/DP >= ${params.consensus_min_depth} && QUAL >= ${minQual}' -Oz -o ${prefix}.filtered.vcf.gz
    bcftools index -f ${prefix}.filtered.vcf.gz

    bcftools +setGT ${prefix}.filtered.vcf.gz -Ou -- -t a -n 'c:1/1' | \\
        bcftools +setGT -Ou -- -t q -i 'GT="1/1" && TYPE="snp" && INFO/AF < ${params.consensus_call_af}' -n 'c:0/1' | \\
        bcftools +setGT -Ou -- -t q -i 'TYPE="indel" && INFO/AF < ${minIndelAf}' -n . | \\
        bcftools +setGT -Oz -o ${prefix}.final.vcf.gz -- -t q -i 'GT="1/1" && INFO/AF >= ${params.consensus_call_af}' -n 'c:1/1'
    bcftools index -f ${prefix}.final.vcf.gz

    # VCF used only for the IUPAC consensus: drop weakly supported mixed SNP sites (fixed sites and indels are untouched)
    bcftools view -e 'GT="het" && TYPE="snp" && (${weakMixedCond})' -Oz -o ${prefix}.final_consensus.vcf.gz ${prefix}.final.vcf.gz
    bcftools index -f ${prefix}.final_consensus.vcf.gz

    bcftools +fill-tags $vcf -Ou -- -t TYPE | \\
        bcftools norm -Ou -a -m - | \\
        bcftools view -f 'PASS,dn,dp,.' -i '((TYPE="snp" && INFO/AF >= 0.5) || (TYPE="indel" && INFO/AF >= 0.5)) && INFO/DP >= ${params.consensus_min_depth} && QUAL >= ${minQual}' -Oz -o ${prefix}.filtered_simple.vcf.gz
    bcftools index -f ${prefix}.filtered_simple.vcf.gz

    bcftools +setGT ${prefix}.filtered_simple.vcf.gz -Oz -o ${prefix}.final_simple.vcf.gz -- -t a -n 'c:1/1'
    bcftools index -f ${prefix}.final_simple.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n 1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.final.vcf.gz
    touch ${prefix}.final.vcf.gz.csi
    touch ${prefix}.final_consensus.vcf.gz
    touch ${prefix}.final_consensus.vcf.gz.csi
    touch ${prefix}.final_simple.vcf.gz
    touch ${prefix}.final_simple.vcf.gz.csi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: stub
    END_VERSIONS
    """
}
