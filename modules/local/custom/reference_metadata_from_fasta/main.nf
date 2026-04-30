process REFERENCE_METADATA_FROM_FASTA {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(sample_id), val(meta), path(reference_fasta)

    output:
    tuple val(meta), path("*.best_reference.tsv"), emit: best_ref_tsv
    tuple val(meta), path("*.best_reference.txt"), emit: best_ref_txt
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ref_name=\$(awk '/^>/{sub(/^>/,""); split(\$0,a,/ /); print a[1]; exit}' ${reference_fasta})
    genotype=""
    subtype=""
    ref_prefix="\${ref_name%%_*}"

    if [[ "\$ref_prefix" =~ ^([0-9]+)([A-Za-z].*)\$ ]]; then
        genotype="\${BASH_REMATCH[1]}"
        subtype="\$ref_prefix"
    fi

    printf "%s\n" "\$ref_name" > ${prefix}.best_reference.txt
    {
        printf "sample_id\tgenotype\tsubtype\tbest_ref\tclose_hits\tother_potential_subtypes\n"
        printf "%s\t%s\t%s\t%s\t\t\n" "${meta.id}" "\$genotype" "\$subtype" "\$ref_name"
    } > ${prefix}.best_reference.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n 1 | cut -d' ' -f4)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat <<-EOF > ${prefix}.best_reference.txt
    stub_reference
    EOF

    cat <<-EOF > ${prefix}.best_reference.tsv
    sample_id\tgenotype\tsubtype\tbest_ref\tclose_hits\tother_potential_subtypes
    ${meta.id}\t1\t1a\tstub_reference\t\t
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n 1 | cut -d' ' -f4)
    END_VERSIONS
    """
}
