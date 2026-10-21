include { PREPARE_KRAKEN2_DB    } from '../../../modules/local/prepare_kraken2_db/main'
include { PREPARE_TAXDUMP       } from '../../../modules/local/prepare_taxdump/main'
include { PREPARE_HOSTILE_INDEX } from '../../../modules/local/prepare_hostile_index/main'
include { PREPARE_ABRICATE_DB   } from '../../../modules/local/prepare_abricate_db/main'
include { DATABASE_MANIFEST     } from '../../../modules/local/database_manifest/main'

//
// Download and prepare the reference databases (--build_databases). Each database is published to
// --db_outdir by its own process (see conf/modules.config); database_manifest.tsv is written last.
// --skip_*_download leaves a database out entirely (no download, no manifest row).
//
workflow PREPARE_DATABASES {
    take:
    kraken_release   // string: YYYYMMDD (unused when params.skip_kraken_download)
    taxdump_month    // string: YYYY-MM  (unused when params.skip_taxdump_download)
    db_root          // string: absolute --db_outdir, written into the manifest paths

    main:
    def ch_rows     = channel.empty()
    def ch_versions = channel.empty()

    if (!params.skip_kraken_download) {
        PREPARE_KRAKEN2_DB(channel.value(kraken_release))
        ch_rows     = ch_rows.mix(PREPARE_KRAKEN2_DB.out.row)
        ch_versions = ch_versions.mix(PREPARE_KRAKEN2_DB.out.versions)
    }
    if (!params.skip_taxdump_download) {
        PREPARE_TAXDUMP(channel.value(taxdump_month))
        ch_rows     = ch_rows.mix(PREPARE_TAXDUMP.out.row)
        ch_versions = ch_versions.mix(PREPARE_TAXDUMP.out.versions)
    }
    if (!params.skip_hostile_download) {
        PREPARE_HOSTILE_INDEX()
        ch_rows     = ch_rows.mix(PREPARE_HOSTILE_INDEX.out.row)
        ch_versions = ch_versions.mix(PREPARE_HOSTILE_INDEX.out.versions)
    }
    if (!params.skip_abricate_download) {
        PREPARE_ABRICATE_DB()
        ch_rows     = ch_rows.mix(PREPARE_ABRICATE_DB.out.row)
        ch_versions = ch_versions.mix(PREPARE_ABRICATE_DB.out.versions)
    }

    DATABASE_MANIFEST(ch_rows.collect(), channel.value(db_root))

    emit:
    manifest = DATABASE_MANIFEST.out.manifest
    versions = ch_versions
}
