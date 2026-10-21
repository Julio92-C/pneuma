include { MERIDIAN_CONFIG } from '../../../modules/local/meridian_config/main'
include { RCF_FOR_MERIDIAN } from '../../../modules/local/rcf_for_meridian/main'
include { MERIDIAN_RUN    } from '../../../modules/local/meridian_run/main'

workflow MERIDIAN {
    take:
    ch_controls        // value channel: List<String>
    ch_samples         // value channel: List<String>, all sample ids (controls included)
    abricate_summary
    kraken_combined
    bracken_arranged
    rcf_csv

    main:
    def metadata = file(params.meridian_metadata, checkIfExists: true)
    // --meridian_group may list several metadata columns: one MERIDIAN comparison (config + run) each.
    // MERIDIAN_CONFIG needs only the metadata and the samplesheet, so it also validates every column
    // (exists, >= 2 non-NA levels among non-controls) at the start of the run, before heavy work.
    def groups = params.meridian_group.toString().split(',').collect { it.trim() }.findAll { it }.unique()
    MERIDIAN_CONFIG(file("${projectDir}/assets/meridian_config.template.yaml"), metadata, ch_controls, channel.fromList(groups))
    // MERIDIAN only understands column 1 + one cnt/una/sco triplet per sample; the raw
    // Recentrifuge CSV (published unchanged under contamination/) has derived columns too.
    RCF_FOR_MERIDIAN(rcf_csv, ch_samples)
    // Shared inputs must be VALUE channels: a queue channel with one item is consumed by the first
    // (group, config) tuple and MERIDIAN_RUN would run for one comparison only.
    MERIDIAN_RUN(MERIDIAN_CONFIG.out.config, metadata, abricate_summary.first(), kraken_combined.first(),
                 bracken_arranged.first(), RCF_FOR_MERIDIAN.out.csv.first())

    // MERIDIAN's stage safety net keeps the run (and the report) alive when one analysis stage fails;
    // surface those failures in the Nextflow log so they are not only visible inside the HTML report.
    MERIDIAN_RUN.out.stage_errors.subscribe { group, csv ->
        def failed = csv.readLines().drop(1).collect { line -> line.tokenize(',')[0].replace('"', '') }
        log.warn "MERIDIAN (${group}): ${failed.size()} stage(s) failed and were skipped: ${failed.join(', ')}. " +
                 "See ${params.outdir}/meridian/${group}/Reports/stage_errors.csv and the red box in that report."
    }

    emit:
    reports  = MERIDIAN_RUN.out.reports.collect()
    versions = MERIDIAN_CONFIG.out.versions.mix(RCF_FOR_MERIDIAN.out.versions, MERIDIAN_RUN.out.versions)
}
