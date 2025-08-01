//
// Subworkflow with functionality specific to the charlesfoster/h2seq pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { dashedLine                } from '../../nf-core/utils_nfcore_pipeline'
include { nfCoreLogo                } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { workflowCitation          } from '../../nf-core/utils_nfcore_pipeline'

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE
========================================================================================
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    help              // boolean: Display help text
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    pre_help_text = nfCoreLogo(monochrome_logs)
    post_help_text = '\n' + workflowCitation() + '\n' + dashedLine(monochrome_logs)
    def String workflow_command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        "${projectDir}/nextflow_schema.json"
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    // Validate FASTQ input
    ch_samplesheet = Channel
        .fromList(samplesheetToList(params.input,"${projectDir}/assets/schema_input.json"))
        .map {
            validateInputSamplesheet(it[0], it[1], it[2], it[3])
        }

    //
    // Create fastq channels by separating long and short reads
    //
    ch_raw_long_reads = ch_samplesheet
        .map { meta, lr, sr1, sr2 ->
            if (lr) {
                meta.single_end = true
                return [ meta, lr ]
            }
        }

    ch_raw_short_reads = ch_samplesheet
        .map { meta, lr, sr1, sr2 ->
            if (sr1 && sr2) {
                meta.single_end = false
                return [ meta, [ sr1, sr2 ] ]
            } else if (sr1) {
                meta.single_end = true
                return [ meta, [ sr1 ] ]
            } else {
                return null // This ensures the channel remains empty if both sr1 and sr2 are missing
            }
        }
        .filter { it != null } // Remove any null values from the channel

    ch_raw_long_reads = ch_raw_long_reads
        .map { meta, lr ->
            def meta_new = meta + [long_reads: true]
            return [ meta_new, lr ]
        }

    ch_raw_short_reads = ch_raw_short_reads
        .map { meta, reads ->
            def meta_new = meta + [long_reads: false]
            return [ meta_new, reads ]
        }

    emit:
    raw_long_reads   = ch_raw_long_reads
    raw_short_reads  = ch_raw_short_reads
    versions    = ch_versions
}

/*
========================================================================================
    SUBWORKFLOW FOR PIPELINE COMPLETION
========================================================================================
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:

    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(summary_params, email, email_on_fail, plaintext_email, outdir, monochrome_logs, multiqc_report.toList())
        }

        completionSummary(monochrome_logs)

        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/
//
// Check and validate pipeline parameters
//
// TODO nf-core: add in other checks here based on all params
def validateInputParameters() {
    genomeExistsError()

    // check for presets
    // Check if proper files provided for reference selection and alignment

    if ( !params.skip_reference_selection && params.virus_preset ) {
        if ( params.virus_preset != "hcv" ) {
            error("[charlesfoster/h2seq] ERROR: Invalid option for parameter '--virus_preset'.")
        }
    } else if (!params.skip_reference_selection && !params.virus_preset && !params.possible_references) {
        error("[charlesfoster/h2seq] ERROR: Invalid combination of parameter '--skip_reference_selection' and parameter '--possible_references'. If '--skip_reference_selection' is not specified, a (multi)fasta file must be specified with '--possible_references'.")
    }

    if (params.skip_reference_selection && !params.reference_fasta) {
        error("[charlesfoster/h2seq] ERROR: Invalid combination of parameter '--skip_reference_selection' and parameter '--reference_fasta'. If '--skip_reference_selection' is NOT specified, a fasta file MUST be specified with '--reference_fasta'.")
    }

    // Check if proper files provided for primer clipping
    if (!params.skip_primer_trimming && !params.primer_fasta && !params.primer_bed) {
        error("[charlesfoster/h2seq] ERROR: Invalid combination of parameters '--skip_primer_trimming', '--primer_fasta' and '--primer_bed'. The '--primer_fasta' parameter is *mandatory* if `--skip_primer_trimming` is not specified and `--primer_bed' is not specified. The '--primer_bed' parameter is *mandatory* if `--skip_primer_trimming` is not specified and `--primer_fasta' is not specified.")
    }

    if (params.filter_host_reads) {
        log.warn("[charlesfoster/h2seq] WARNING: host filtration not yet implemented so nothing will happen. Watch this space...")
        // Check if parameters for host contamination removal are valid
        if ( params.host_fasta && params.host_genome) {
            error('[charlesfoster/h2seq] ERROR: Both host fasta reference and iGenomes genome are specified to remove host contamination! Invalid combination, please specify either --host_fasta or --host_genome.')
        }

        // below taken from nf-core/MAG
        if ( params.host_genome ) {
            if (!params.genomes) {
                error('[charlesfoster/h2seq] ERROR: No config file containing genomes provided!')
            }
            // Check if host genome exists in the config file
            if (!params.genomes.containsKey(params.host_genome)) {
                error('=============================================================================\n' +
                        "  Host genome '${params.host_genome}' not found in any config files provided to the pipeline.\n" +
                        '  Currently, the available genome keys are:\n' +
                        "  ${params.genomes.keySet().join(', ')}\n" +
                        '===================================================================================')
            }
            if ( !params.genomes[params.host_genome].fasta ) {
                error("[charlesfoster/h2seq] ERROR: No fasta file specified for the host genome ${params.host_genome}!")
            }
            if ( !params.genomes[params.host_genome].bowtie2 ) {
                error("[charlesfoster/h2seq] ERROR: No Bowtie 2 index file specified for the host genome ${params.host_genome}!")
            }
        }
    }
}

//
// Validate channels from input samplesheet - taken from MAG
//
def validateInputSamplesheet(meta, lr, sr1, sr2 ) {

        if ( sr1 && !sr2 && !params.single_end ) { error("[charlesfoster/h2seq] ERROR: Single-end data must be executed with `--single_end`. Note that it is not possible to mix single- and paired-end data (for short reads) in one run! Check input TSV for sample: ${meta.id}") }
        if ( sr2 && params.single_end ) { error("[charlesfoster/h2seq] ERROR: Paired-end data must be executed without `--single_end`. Note that it is not possible to mix single- and paired-end data (for short reads) in one run! Check input TSV for sample: ${meta.id}") }

    return [meta, lr, sr1, sr2]
}

//
// Validate channels from input samplesheet - default, not using
//
// def validateInputSamplesheet(input) {
//     def (metas, fastqs) = input[1..2]

//     // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
//     def endedness_ok = metas.collect{ it.single_end }.unique().size == 1
//     if (!endedness_ok) {
//         error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
//     }

//     return [ metas[0], fastqs ]
// }

//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        String[] manifest_doi = meta.manifest_map.doi.tokenize(",")
        for (String doi_ref: manifest_doi) temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
