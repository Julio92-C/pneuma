include { MERIDIAN_CONFIG } from '../../../modules/local/meridian_config/main'
include { MERIDIAN_RUN    } from '../../../modules/local/meridian_run/main'

workflow MERIDIAN {
    take:
    ch_controls        // value channel: List<String>
    abricate_summary
    kraken_combined
    bracken_arranged
    rcf_csv

    main:
    def metadata = file(params.meridian_metadata, checkIfExists: true)
    MERIDIAN_CONFIG(file("${projectDir}/assets/meridian_config.template.yaml"), metadata, ch_controls)
    MERIDIAN_RUN(MERIDIAN_CONFIG.out.config, metadata, abricate_summary, kraken_combined, bracken_arranged, rcf_csv)

    emit:
    reports  = MERIDIAN_RUN.out.reports
    versions = MERIDIAN_CONFIG.out.versions.mix(MERIDIAN_RUN.out.versions)
}
