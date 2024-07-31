process HCV_GLUE {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.html") , emit: report

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cp ${fasta} ${prefix}.reporting.fasta
    docker run --rm --name ${prefix} -v \$PWD:\$PWD -w \$PWD --link gluetools-mysql cvrbioinformatics/gluetools:latest gluetools.sh -i project hcv  module phdrReportingController invoke-function reportFastaAsHtml ${prefix}.reporting.fasta ${prefix}.html
    """
}
