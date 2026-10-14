#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Julio92-C/pneuma
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/Julio92-C/pneuma
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PNEUMA  } from './workflows/pneuma'
include { PREPARE_DATABASES       } from './subworkflows/local/prepare_databases/main'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_pneuma_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_pneuma_pipeline'
include { resolveDatabases        } from './subworkflows/local/utils_nfcore_pneuma_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow JULIO92C_PNEUMA {

    take:
    samplesheet // channel: samplesheet read in from --input
    databases   // map: kraken2_db, taxdump, host_index, abricate_db, database_manifest

    main:

    //
    // WORKFLOW: Run pipeline
    //
    PNEUMA (
        samplesheet,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir,
        databases,
    )
    emit:
    multiqc_report = PNEUMA.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.build_databases ? params.db_outdir : params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    def ch_multiqc_report = channel.empty()
    if (params.build_databases) {
        //
        // WORKFLOW: Database mode - download and prepare the databases, no analysis
        //
        PREPARE_DATABASES (
            params.kraken_release,
            params.taxdump_month,
            file(params.db_outdir).toAbsolutePath().toString()
        )
    } else {
        //
        // WORKFLOW: Run main workflow (databases resolved from explicit params and/or --db_dir)
        //
        JULIO92C_PNEUMA (
            PIPELINE_INITIALISATION.out.samplesheet,
            resolveDatabases(params)
        )
        ch_multiqc_report = JULIO92C_PNEUMA.out.multiqc_report
    }
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.build_databases ? params.db_outdir : params.outdir,
        params.monochrome_logs,
        ch_multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
