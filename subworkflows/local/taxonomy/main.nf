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

// True if a Kraken2 report has any non-unclassified taxon line with a clade count > 0
def hasClassifiedReads(report) {
    return file(report).readLines().any { line ->
        def f = line.split('	')
        f.size() >= 6 && f[3].trim() != 'U' && f[1].trim().isLong() && f[1].trim().toLong() > 0
    }
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

    // Reports with no classified reads (e.g. negative controls, near-empty samples) make
    // KrakenTools combine_kreports.py crash (no root node). Keep them out of the combined
    // tables and Bracken; they still flow through assignments/classified/Krona below.
    ch_report_split = KRAKEN2_KRAKEN2.out.report.branch { _meta, report ->
        classified: hasClassifiedReads(report)
        empty: true
    }
    ch_report_split.empty.subscribe { meta, _report ->
        log.warn "TAXONOMY: sample '${meta.id}' has no classified reads; excluded from the Kraken2/Bracken combined tables"
    }
    ch_kraken_ok = ch_report_split.classified

    // Bracken is skipped for LOW_DEPTH samples (no reliable re-estimation)
    BRACKEN_BRACKEN(ch_kraken_ok.filter { meta, _r -> !meta.low_depth }, ch_db)

    // Bracken's per-sample errorStrategy is 'ignore' (exit 1 = no taxa above threshold), so a
    // sample can silently lack a Bracken output. Compare the two sets here: warn per missing
    // sample, and fail loudly if there is no Bracken output at all (likely misconfiguration).
    ch_checked = collectSorted(ch_kraken_ok)
        .combine(collectSorted(BRACKEN_BRACKEN.out.txt))
        .map { k_names, k_files, b_names, b_files ->
            if (!k_names) {
                error "TAXONOMY: no sample has classified Kraken2 reads; cannot build combined tables"
            }
            if (!b_names) {
                error "TAXONOMY: Bracken produced no output for any sample. Likely causes: params.bracken_read_length (${params.bracken_read_length}) has no matching database<N>mers.kmer_distrib in the Kraken2 DB, the DB lacks Bracken files, or params.bracken_level/bracken_threshold are invalid. Check the BRACKEN_BRACKEN task logs."
            }
            def missing = k_names - b_names
            if (missing) {
                log.warn "TAXONOMY: no Bracken output for ${missing.size()} sample(s) with Kraken2 reports (low depth, or no taxa above threshold): ${missing.join(', ')}; absent from bracken_arranged.csv"
            }
            [ k_names, k_files, b_names, b_files ]
        }
        .multiMap { k_names, k_files, b_names, b_files ->
            k_names: k_names
            k_files: k_files
            b_names: b_names
            b_files: b_files
        }

    COMBINE_KRAKEN('kraken2_combined', ch_checked.k_names, ch_checked.k_files)
    COMBINE_BRACKEN('bracken_combined', ch_checked.b_names, ch_checked.b_files)
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
