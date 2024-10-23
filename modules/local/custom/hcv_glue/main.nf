process HCV_GLUE {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.html") , emit: report, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def randomName = UUID.randomUUID().toString().take(8)
    """
    BNAME=\$(basename "$fasta" .reporting.fasta)
    FASTA_DIR=\$(dirname \$(readlink -f ${fasta}))
    docker run --rm --name ${prefix}_${randomName} -v \$PWD:\$PWD -v \$FASTA_DIR:\$FASTA_DIR -w \$PWD --link gluetools-mysql cvrbioinformatics/gluetools:latest gluetools.sh -i project hcv  module phdrReportingController invoke-function reportFastaAsHtml ${fasta} \${BNAME}.html
    """
}
