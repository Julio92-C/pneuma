include { BUILD_PATHOGEN_REPORT } from '../../../modules/local/pathogen_report/main'
include { pathogenPanelPath     } from '../utils_nfcore_pneuma_pipeline/main'

// All files of one kind from the [ kind, file ] channel, as one list ([] when there are none).
def filesOf(ch, String kind) {
    return ch.filter { k, _f -> k == kind }.map { _k, f -> f }.collect().ifEmpty([])
}

workflow PATHOGEN_REPORT {
    take:
    panel_resolved   // path: panel_resolved.tsv
    candidates       // path: pathogen_candidates.tsv
    samples_tsv      // path: pathogen_samples.tsv
    units_tsv        // path: units.tsv
    report_files     // channel: [ kind, file ]
    read_abricate    // value: list of <sample>_abricate.csv (FUNCTIONAL_READS)
    manifest         // path or []

    main:
    BUILD_PATHOGEN_REPORT(
        panel_resolved, candidates, samples_tsv, units_tsv,
        filesOf(report_files, 'gate'), filesOf(report_files, 'species'), filesOf(report_files, 'assemblies'),
        filesOf(report_files, 'mlst'), filesOf(report_files, 'amrfinderplus'), filesOf(report_files, 'abricate'),
        filesOf(report_files, 'mobsuite'), filesOf(report_files, 'mixed'), filesOf(report_files, 'checkm2'),
        filesOf(report_files, 'quast'), read_abricate,
        file(pathogenPanelPath(params), checkIfExists: true), manifest
    )

    emit:
    cohort = BUILD_PATHOGEN_REPORT.out.cohort   // path: pathogen_cohort.tsv
    html   = BUILD_PATHOGEN_REPORT.out.html     // path(s): html/<sample>.pathogens.html
}
