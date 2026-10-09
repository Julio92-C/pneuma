include { UNTAR                     } from '../../../modules/nf-core/untar/main'
include { KRAKEN2_KRAKEN2           } from '../../../modules/nf-core/kraken2/kraken2/main'
include { BRACKEN_BRACKEN           } from '../../../modules/nf-core/bracken/bracken/main'
include { KRAKENTOOLS_KREPORT2KRONA } from '../../../modules/nf-core/krakentools/kreport2krona/main'
include { KRONA_KTIMPORTTEXT        } from '../../../modules/nf-core/krona/ktimporttext/main'
include { COMBINE_KREPORTS as COMBINE_KRAKEN  } from '../../../modules/local/combine_kreports/main'
include { COMBINE_KREPORTS as COMBINE_BRACKEN } from '../../../modules/local/combine_kreports/main'
include { ARRANGE_REPORTS           } from '../../../modules/local/arrange_reports/main'

// Collect [meta, file] items into a single [ [ids], [files] ] tuple, sorted by sample id
def collectSorted(ch) {
    return ch.toSortedList { a, b -> a[0].id <=> b[0].id }
        .map { items -> [ items.collect { item -> item[0].id }, items.collect { item -> item[1] } ] }
}

workflow TAXONOMY {
    take:
    ch_reads     // [ meta, fastq.gz ]
    kraken2_db   // path: directory or .tar.gz

    main:
    ch_versions = channel.empty()

    if (kraken2_db.toString().endsWith('.tar.gz')) {
        UNTAR(channel.value([ [id:'kraken2_db'], file(kraken2_db, checkIfExists: true) ]))
        ch_db = UNTAR.out.untar.map { _meta, dir -> dir }.first()
    } else {
        ch_db = channel.value(file(kraken2_db, checkIfExists: true))
    }

    KRAKEN2_KRAKEN2(ch_reads, ch_db, true, true)

    // Bracken is skipped for LOW_DEPTH samples (no reliable re-estimation)
    BRACKEN_BRACKEN(KRAKEN2_KRAKEN2.out.report.filter { meta, _r -> !meta.low_depth }, ch_db)

    def k = collectSorted(KRAKEN2_KRAKEN2.out.report).multiMap { names, files -> names: names; files: files }
    def b = collectSorted(BRACKEN_BRACKEN.out.txt).multiMap { names, files -> names: names; files: files }

    COMBINE_KRAKEN('kraken2_combined', k.names, k.files)
    COMBINE_BRACKEN('bracken_combined', b.names, b.files)
    ch_versions = ch_versions.mix(COMBINE_KRAKEN.out.versions)

    ARRANGE_REPORTS(COMBINE_KRAKEN.out.txt, COMBINE_BRACKEN.out.txt)
    ch_versions = ch_versions.mix(ARRANGE_REPORTS.out.versions)

    KRAKENTOOLS_KREPORT2KRONA(KRAKEN2_KRAKEN2.out.report)
    KRONA_KTIMPORTTEXT(KRAKENTOOLS_KREPORT2KRONA.out.txt)

    emit:
    classified       = KRAKEN2_KRAKEN2.out.classified_reads_fastq
    assignments      = KRAKEN2_KRAKEN2.out.classified_reads_assignment
    reports          = KRAKEN2_KRAKEN2.out.report
    kraken_combined  = ARRANGE_REPORTS.out.kraken
    bracken_arranged = ARRANGE_REPORTS.out.bracken
    versions         = ch_versions
}
