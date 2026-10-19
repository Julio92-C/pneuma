/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_pneuma_pipeline'
include { PREPROCESS             } from '../subworkflows/local/preprocess/main'
include { TAXONOMY               } from '../subworkflows/local/taxonomy/main'
include { CONTAMINATION          } from '../subworkflows/local/contamination/main'
include { FUNCTIONAL_READS       } from '../subworkflows/local/functional_reads/main'
include { MERIDIAN               } from '../subworkflows/local/meridian/main'
include { MULTIQC_TABLES         } from '../modules/local/multiqc_tables/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PNEUMA {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir
    databases // map: kraken2_db, taxdump, host_index, abricate_db, database_manifest (see resolveDatabases)

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    PREPROCESS(ch_samplesheet, file(databases.host_index, checkIfExists: true))
    TAXONOMY(PREPROCESS.out.reads, databases.kraken2_db)
    CONTAMINATION(TAXONOMY.out.assignments, databases.taxdump)
    FUNCTIONAL_READS(TAXONOMY.out.classified, databases.abricate_db)

    ch_versions = ch_versions.mix(PREPROCESS.out.versions, TAXONOMY.out.versions,
                                  CONTAMINATION.out.versions, FUNCTIONAL_READS.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(PREPROCESS.out.multiqc)

    // Hostile, Kraken2 classification, ABRicate and Recentrifuge have no (usable) MultiQC module:
    // summarise them as custom-content tables. Kraken2/Bracken reports go to MultiQC's kraken module as is.
    MULTIQC_TABLES(
        PREPROCESS.out.hostile_log.collect(),
        TAXONOMY.out.reports.map { _meta, f -> f }.collect(),
        FUNCTIONAL_READS.out.per_sample.map { _meta, f -> f }.collect(),
        CONTAMINATION.out.rcf_stat
    )
    ch_versions = ch_versions.mix(MULTIQC_TABLES.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        TAXONOMY.out.reports.map { _meta, f -> f },
        TAXONOMY.out.bracken_reports.map { _meta, f -> f },
        MULTIQC_TABLES.out.tables
    )

    if (!params.skip_meridian) {
        def ch_controls = ch_samplesheet
            .filter { meta, _fastqs -> meta.is_control }
            .map { meta, _fastqs -> meta.id }
            .collect()
            .ifEmpty([])
        def ch_samples = ch_samplesheet.map { meta, _fastqs -> meta.id }.collect().map { it.sort() }
        MERIDIAN(ch_controls, ch_samples, FUNCTIONAL_READS.out.summary, TAXONOMY.out.kraken_combined,
                 TAXONOMY.out.bracken_arranged, CONTAMINATION.out.rcf_csv)
        ch_versions = ch_versions.mix(MERIDIAN.out.versions)
    }

    if (databases.database_manifest) {
        channel.fromPath(databases.database_manifest, checkIfExists: true)
            .collectFile(name: 'database_manifest.tsv', storeDir: "${outdir}/pipeline_info")
    }

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'pneuma_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'pneuma'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
