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

    emit:
    reports  = MERIDIAN_RUN.out.reports.collect()
    versions = MERIDIAN_CONFIG.out.versions.mix(RCF_FOR_MERIDIAN.out.versions, MERIDIAN_RUN.out.versions)
}
