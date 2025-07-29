process CLAIR3 {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularityclair3:1.1.1--py310h779eee5_0' :
        'docker.io/hkubal/clair3:v1.1.1' }"

    input:
    tuple val(meta), path(bam), path(bam_index), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.vcf.gz"), path("*.tbi") , emit: vcf
    path "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    /opt/bin/run_clair3.sh  \\
        --threads=${task.cpus} \\
        --bam_fn=$bam \\
        --sample_name=$prefix \\
        --ref_fn=$fasta \\
        --output=clair3_output \\
        $args

    cp clair3_output/merge_output.vcf.gz ./${prefix}.vcf.gz
    cp clair3_output/merge_output.vcf.gz.tbi ./${prefix}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(echo \$(/opt/bin/run_clair3.sh --version 2>&1) | sed 's/^.*v//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(echo \$(/opt/bin/run_clair3.sh --version 2>&1) | sed 's/^.*v//')
    END_VERSIONS
    """
}

