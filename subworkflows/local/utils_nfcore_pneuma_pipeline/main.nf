//
// Subworkflow with functionality specific to the Julio92-C/pneuma pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

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

    def before_text = ""
    def after_text = ""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
        false
    )

    //
    // Parameter guards that the JSON schema cannot express (nf-schema has no conditional "required")
    //
    if (params.build_databases) {
        checkBuildDatabasesParams(params)
    } else {
        if (!input || !params.outdir) {
            error("--input and --outdir are required (or use --build_databases with --db_outdir to only download the databases)")
        }
        def dbs = resolveDatabases(params)
        def keys = ['kraken2_db', 'taxdump', 'host_index', 'abricate_db', 'database_manifest'] +
                   (params.pathogens ? ['checkm2_db', 'bakta_db', 'amrfinderplus_db', 'mobsuite_db', 'medaka_model_path', 'panel_references_dir'] : [])
        keys.each { k ->
            def source = params[k] ? 'explicit parameter' : (dbs[k] ? "--db_dir (${params.db_dir})" : 'not set')
            log.info "database ${k.padRight(20)} <- ${source}${dbs[k] ? ': ' + dbs[k] : ''}"
        }
        if (!params.skip_meridian && !(params.meridian_metadata && params.meridian_group)) {
            error("--meridian_metadata and --meridian_group are required unless --skip_meridian is set")
        }
        if (!params.skip_meridian && !params.meridian_group.toString().split(',').any { it.trim() }) {
            error("--meridian_group lists no metadata column")
        }
        if (!params.skip_meridian && params.meridian_random_effect_tests && !params.meridian_random_effect) {
            error("--meridian_random_effect_tests needs --meridian_random_effect (the patient column to permute within)")
        }
    }

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Create channel from input file provided through params.input (analysis mode only)
    //
    def ch_samplesheet = channel.empty()
    if (!params.build_databases) {
        def rows = samplesheetToList(input, "${projectDir}/assets/schema_input.json")
        warnRunsWithoutControls(rows)
        notePatients(rows)

        ch_samplesheet = channel
            .fromList(rows)
            .map { meta, fastq -> [ meta.id, meta, fastq ] }
            .groupTuple()
            .map { id, metas, fastqs -> groupRuns(id, metas, fastqs) }
    }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Guards for --build_databases mode (no samplesheet, no analysis).
//
def checkBuildDatabasesParams(p) {
    def problems = []
    if (!p.db_outdir) {
        problems.add("--db_outdir is required with --build_databases")
    }
    if (!p.skip_kraken_download && !(p.kraken_release ==~ /\d{8}/)) {
        problems.add("--kraken_release YYYYMMDD is required (see https://benlangmead.github.io/aws-indexes/k2), or pass --skip_kraken_download")
    }
    if (!p.skip_taxdump_download && !(p.taxdump_month ==~ /\d{4}-(0[1-9]|1[0-2])/)) {
        problems.add("--taxdump_month YYYY-MM is required (match the Kraken2 release month), or pass --skip_taxdump_download")
    }
    if (p.skip_kraken_download && p.skip_taxdump_download && p.skip_hostile_download && p.skip_abricate_download && p.skip_pathogen_download) {
        problems.add("every download is skipped; nothing to build")
    }
    if (!p.skip_pathogen_download && !(p.bakta_db_type in ['light', 'full'])) {
        problems.add("--bakta_db_type must be 'light' or 'full'")
    }
    if (problems) {
        error("Invalid --build_databases options:\n  - " + problems.join("\n  - "))
    }
    if (p.input) {
        log.warn "--input is ignored with --build_databases (no analysis is run)"
    }
}

//
// Read database_manifest.tsv (name, version, path, sha256) into [name: [version:, path:, sha256:]].
//
def readDatabaseManifest(manifest) {
    def rows = [:]
    manifest.readLines().drop(1).findAll { it.trim() }.each { line ->
        def c = line.split('\t', -1)
        if (c.size() >= 4) {
            rows[c[0]] = [version: c[1], path: c[2], sha256: c[3]]
        } else {
            log.warn "Ignoring malformed row in ${manifest} (expected 4 tab-separated columns, found ${c.size()}): '${line}'"
        }
    }
    return rows
}

//
// Sub-directories of root named <prefix>*, narrowed to <prefix><version> when that exists.
//
def findVersionedDirs(root, prefix, version) {
    def all = (root.listFiles() ?: []).findAll { it.isDirectory() && it.name.startsWith(prefix) }.sort { it.name }
    def exact = version ? all.findAll { it.name == prefix + version } : []
    return exact ?: all
}

//
// Hostile minimap2 indexes (human-t2t-hla*.mmi) anywhere below dir.
//
def findHostileIndexes(dir) {
    def entries = dir.listFiles() ?: []
    def here = entries.findAll { it.isFile() && it.name.startsWith('human-t2t-hla') && it.name.endsWith('.mmi') }
    return (here + entries.findAll { it.isDirectory() }.collectMany { findHostileIndexes(it) }).sort { it.toString() }
}

//
// Pick the only candidate; add a problem when there are none or several.
//
def pickOne(candidates, what, problems, hint) {
    if (candidates.size() == 1) {
        return candidates[0]
    }
    if (candidates.isEmpty()) {
        problems.add("${what}: nothing found. ${hint}")
    } else {
        problems.add("${what}: several candidates (${candidates.collect { it.name }.join(', ')}). ${hint}")
    }
    return null
}

def missingFiles(dir, names) {
    return names.findAll { !dir.resolve(it).exists() }
}

//
// The panel files default to the ones bundled in assets/ (params are null so that the schema
// carries no machine-specific default). Every consumer goes through these two helpers.
//
def pathogenPanelPath(p) {
    return p.pathogen_panel ?: "${workflow.projectDir}/assets/pathogen_panel.csv"
}

def panelReferencesPath(p) {
    return p.panel_references ?: "${workflow.projectDir}/assets/panel_references.tsv"
}

//
// Resolve the database parameters for analysis mode. Explicit parameters always win; the rest
// is filled from --db_dir (database_manifest.tsv + the layout written by --build_databases).
// Returns [kraken2_db:, taxdump:, host_index:, abricate_db:, database_manifest:, checkm2_db:, bakta_db:,
// amrfinderplus_db:, mobsuite_db:, medaka_model_path:, panel_references_dir:] as strings or null.
// Calls error() listing every problem found. kraken2_db, taxdump and host_index are mandatory; with
// --pathogens, so are checkm2_db, bakta_db, amrfinderplus_db, medaka_model_path and
// panel_references_dir (mobsuite_db only with --run_mobsuite).
//
def resolveDatabases(p) {
    def res = [
        kraken2_db        : p.kraken2_db ?: null,
        taxdump           : p.taxdump ?: null,
        host_index        : p.host_index ?: null,
        abricate_db       : p.abricate_db ?: null,
        database_manifest : p.database_manifest ?: null,
        checkm2_db           : p.checkm2_db ?: null,
        bakta_db             : p.bakta_db ?: null,
        amrfinderplus_db     : p.amrfinderplus_db ?: null,
        mobsuite_db          : p.mobsuite_db ?: null,
        medaka_model_path    : p.medaka_model_path ?: null,
        panel_references_dir : p.panel_references_dir ?: null,
    ]
    def problems = []

    if (p.db_dir) {
        def root = file(p.db_dir)
        if (!root.isDirectory()) {
            error("--db_dir '${p.db_dir}' is not a directory")
        }
        def manifestFile = root.resolve('database_manifest.tsv')
        def manifest = manifestFile.exists() ? readDatabaseManifest(manifestFile) : [:]
        if (!manifestFile.exists()) {
            log.warn "No database_manifest.tsv in --db_dir '${p.db_dir}': database versions cannot be recorded and the newest matching directories are picked by name"
        }
        def hint = "Pass the parameter explicitly or rebuild with --build_databases."

        if (!res.database_manifest && manifestFile.exists()) {
            res.database_manifest = manifestFile.toString()
        }

        if (!res.kraken2_db) {
            def found = pickOne(findVersionedDirs(root, 'kraken2_pluspf_', manifest.kraken2_pluspf?.version),
                                "--kraken2_db (kraken2_pluspf_* in ${p.db_dir})", problems, hint)
            if (found) {
                def missing = missingFiles(found, ['hash.k2d', 'opts.k2d', 'taxo.k2d'])
                if (missing) {
                    problems.add("--kraken2_db: ${found} is missing ${missing.join(', ')}")
                } else {
                    res.kraken2_db = found.toString()
                }
            }
        }
        if (!res.taxdump) {
            def found = pickOne(findVersionedDirs(root, 'taxdump_', manifest.ncbi_taxdump?.version),
                                "--taxdump (taxdump_* in ${p.db_dir})", problems, hint)
            if (found) {
                def missing = missingFiles(found, ['nodes.dmp', 'names.dmp'])
                if (missing) {
                    problems.add("--taxdump: ${found} is missing ${missing.join(', ')}")
                } else {
                    res.taxdump = found.toString()
                }
            }
        }
        if (!res.host_index) {
            def indexes = root.resolve('hostile').isDirectory() ? findHostileIndexes(root.resolve('hostile')) : []
            def wanted = manifest.hostile_index ? file(manifest.hostile_index.path).name : null
            def named = wanted ? indexes.findAll { it.name == wanted } : []
            def found = pickOne(named ?: indexes, "--host_index (hostile/human-t2t-hla*.mmi in ${p.db_dir})", problems, hint)
            if (found) {
                res.host_index = found.toString()
            }
        }
        if (!res.abricate_db && root.resolve('abricate').isDirectory()) {
            def adir = root.resolve('abricate')
            def missing = p.abricate_dbs.tokenize(',').findAll { !adir.resolve("${it}/sequences").exists() }
            if (missing) {
                problems.add("--abricate_db: ${adir} has no sequences file for ABRicate database(s) ${missing.join(', ')} (see --abricate_dbs)")
            } else {
                res.abricate_db = adir.toString()
            }
        }
        if (!res.abricate_db && !root.resolve('abricate').isDirectory()) {
            log.warn "No abricate/ folder in --db_dir '${p.db_dir}' and no --abricate_db: ABRicate will use the databases bundled in the container (2021-Mar-27)"
        }
        if (p.pathogens) {
            resolvePathogenDatabases(p, root, res, problems, hint)
        }
    }

    ['host_index', 'kraken2_db', 'taxdump'].each { k ->
        if (!res[k] && !problems.any { it.contains("--${k}") }) {
            problems.add("--${k} is required: pass it or provide --db_dir containing it")
        }
    }
    if (p.pathogens) {
        def needed = ['checkm2_db', 'bakta_db', 'amrfinderplus_db', 'medaka_model_path', 'panel_references_dir'] +
                     (p.run_mobsuite ? ['mobsuite_db'] : [])
        needed.each { k ->
            if (!res[k] && !problems.any { pr -> pr.contains("--${k}") }) {
                problems.add("--${k} is required for pathogen reconstruction: pass it, provide --db_dir containing it, or do not pass --pathogens")
            }
        }
        if (res.panel_references_dir) {
            checkPanelReferences(file(res.panel_references_dir), file(panelReferencesPath(p)), problems)
        }
    }
    if (problems) {
        error("Cannot resolve the databases:\n  - " + problems.join("\n  - "))
    }
    return res
}

//
// Fill the pathogen database paths from a --db_dir built by --build_databases.
//
def resolvePathogenDatabases(p, root, res, problems, hint) {
    if (!res.checkm2_db && root.resolve('checkm2').isDirectory()) {
        def dmnd = (root.resolve('checkm2').listFiles() ?: []).findAll { f -> f.name.endsWith('.dmnd') }.sort { f -> f.name }
        def found = pickOne(dmnd, "--checkm2_db (checkm2/*.dmnd in ${p.db_dir})", problems, hint)
        if (found) { res.checkm2_db = found.toString() }
    }
    if (!res.bakta_db && root.resolve('bakta').isDirectory()) {
        def found = pickOne(findVersionedDirs(root.resolve('bakta'), 'db', null),
                            "--bakta_db (bakta/db* in ${p.db_dir})", problems, hint)
        if (found) {
            if (found.resolve('version.json').exists()) {
                res.bakta_db = found.toString()
            } else {
                problems.add("--bakta_db: ${found} has no version.json")
            }
        }
    }
    def amr = root.resolve('amrfinderplus/amrfinderdb.tar.gz')
    if (!res.amrfinderplus_db && amr.exists()) { res.amrfinderplus_db = amr.toString() }
    def mob = root.resolve('mobsuite')
    if (!res.mobsuite_db && mob.isDirectory()) { res.mobsuite_db = mob.toString() }
    def model = root.resolve("medaka/${p.medaka_model}_model_pt.tar.gz")
    if (!res.medaka_model_path && model.exists()) { res.medaka_model_path = model.toString() }
    def refs = root.resolve('panel_references')
    if (!res.panel_references_dir && refs.isDirectory()) { res.panel_references_dir = refs.toString() }
}

//
// Every accession of --panel_references needs <dir>/<accession>.fna.gz.
//
def checkPanelReferences(dir, references, problems) {
    def accessions = references.readLines().drop(1).findAll { line -> line.trim() }.collect { line -> line.split('\t')[1].trim() }
    def missing = accessions.findAll { acc -> !dir.resolve("${acc}.fna.gz").exists() }
    if (missing) {
        problems.add("--panel_references_dir: ${dir} lacks ${missing.size()} reference genome(s): ${missing.take(10).join(', ')}${missing.size() > 10 ? ', ...' : ''}. Rebuild with --build_databases (or --panel_references_dir).")
    }
}

//
// Optional samplesheet values arrive as null or [] from nf-schema; normalise to a String.
//
def metaText(value) {
    return value in [null, []] ? '' : value.toString()
}

//
// Merge all runs of one sample into a single meta + ordered fastq list.
// batch / site / is_control / patient must agree across a sample's runs.
//
def groupRuns(id, metas, fastqs) {
    ['batch', 'site', 'is_control', 'patient'].each { key ->
        def values = metas.collect { m -> key == 'is_control' ? m[key] : metaText(m[key]) }.unique()
        if (values.size() > 1) {
            error("Sample '${id}' has conflicting '${key}' values across runs: ${values}")
        }
    }
    def order = (0..<metas.size()).toList().sort { i -> metas[i].run_id }
    def meta = [
        id         : id,
        batch      : metas[0].batch,
        site       : metaText(metas[0].site),
        patient    : metaText(metas[0].patient),
        is_control : metas[0].is_control as Boolean,
        run_ids    : order.collect { i -> metas[i].run_id },
        single_end : true
    ]
    return [ meta, order.collect { i -> fastqs[i] } ]
}

//
// Patient-level pooling of pathogen reads needs `patient`; say once when there is none.
//
def notePatients(rows) {
    def withPatient = rows.findAll { row -> metaText(row[0].patient) }
    if (!withPatient) {
        log.info "No 'patient' values in the samplesheet: patient-level pooling of pathogen reads is skipped (each sample is reconstructed on its own)."
    }
    def controls = withPatient.findAll { row -> row[0].is_control }.collect { row -> row[0].id }.unique()
    if (controls) {
        log.warn "Negative control(s) with a patient value are never pooled: ${controls.join(', ')}"
    }
}

//
// Warn (not fail) for every run_id that has no negative control.
//
def warnRunsWithoutControls(rows) {
    def runs = rows.collect { it[0].run_id }.unique().sort()
    def withCtrl = rows.findAll { it[0].is_control }.collect { it[0].run_id }.unique()
    def missing = runs - withCtrl
    if (missing) {
        log.warn "No negative control declared for run(s): ${missing.join(', ')}. " +
                 "Contamination assessment for these runs relies on controls from other runs."
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
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
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
