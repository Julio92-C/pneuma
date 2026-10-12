include { RECENTRIFUGE } from '../../../modules/local/recentrifuge/main'

workflow CONTAMINATION {
    take:
    ch_assignments   // [ meta, classifiedreads.txt ]
    taxdump          // path: directory with nodes.dmp + names.dmp

    main:
    ch_controls = ch_assignments.filter { meta, _f -> meta.is_control }.map { meta, _f -> meta.id }.collect().ifEmpty([])
    ch_files    = ch_assignments.map { _meta, f -> f }.collect()

    RECENTRIFUGE(ch_controls, ch_files, channel.value(file(taxdump, checkIfExists: true)))

    emit:
    rcf_csv  = RECENTRIFUGE.out.csv
    html     = RECENTRIFUGE.out.html.mix(RECENTRIFUGE.out.controls_html)
    versions = RECENTRIFUGE.out.versions
}
