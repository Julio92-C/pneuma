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
    MERIDIAN_CONFIG(file("${projectDir}/assets/meridian_config.template.yaml"), metadata, ch_controls)
    // MERIDIAN only understands column 1 + one cnt/una/sco triplet per sample; the raw
    // Recentrifuge CSV (published unchanged under contamination/) has derived columns too.
    RCF_FOR_MERIDIAN(rcf_csv, ch_samples)
    MERIDIAN_RUN(MERIDIAN_CONFIG.out.config, metadata, abricate_summary, kraken_combined, bracken_arranged, RCF_FOR_MERIDIAN.out.csv)

    emit:
    reports  = MERIDIAN_RUN.out.reports
    versions = MERIDIAN_CONFIG.out.versions.mix(RCF_FOR_MERIDIAN.out.versions, MERIDIAN_RUN.out.versions)
}
