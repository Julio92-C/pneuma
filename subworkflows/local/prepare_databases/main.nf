include { PREPARE_KRAKEN2_DB        } from '../../../modules/local/prepare_kraken2_db/main'
include { PREPARE_TAXDUMP           } from '../../../modules/local/prepare_taxdump/main'
include { PREPARE_HOSTILE_INDEX     } from '../../../modules/local/prepare_hostile_index/main'
include { PREPARE_ABRICATE_DB       } from '../../../modules/local/prepare_abricate_db/main'
include { CHECKM2_DATABASEDOWNLOAD  } from '../../../modules/nf-core/checkm2/databasedownload/main'
include { BAKTA_BAKTADBDOWNLOAD     } from '../../../modules/nf-core/bakta/baktadbdownload/main'
include { AMRFINDERPLUS_UPDATE      } from '../../../modules/nf-core/amrfinderplus/update/main'
include { DATABASE_ROW              } from '../../../modules/local/database_row/main'
include { PREPARE_MOBSUITE_DB       } from '../../../modules/local/prepare_mobsuite_db/main'
include { PREPARE_MEDAKA_MODEL      } from '../../../modules/local/prepare_medaka_model/main'
include { PREPARE_MLST_SCHEMES      } from '../../../modules/local/prepare_mlst_schemes/main'
include { PREPARE_PANEL_REFERENCES  } from '../../../modules/local/prepare_panel_references/main'
include { DATABASE_MANIFEST         } from '../../../modules/local/database_manifest/main'

//
// Download and prepare the reference databases (--build_databases). Each database is published to
// --db_outdir by its own process (see conf/modules.config); database_manifest.tsv is written last
// and keeps the rows of an existing manifest for databases not rebuilt in this run.
// --skip_*_download leaves a database out (no download, no new manifest row).
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
    if (!params.skip_pathogen_download) {
        CHECKM2_DATABASEDOWNLOAD(channel.value(params.checkm2_zenodo_id ?: ''))
        BAKTA_BAKTADBDOWNLOAD()
        AMRFINDERPLUS_UPDATE()
        def ch_downloaded = CHECKM2_DATABASEDOWNLOAD.out.database
            .map { meta, dmnd -> [ 'checkm2_db', '07', meta.version.toString(), "checkm2/${dmnd.name}".toString(), dmnd ] }
            .mix(BAKTA_BAKTADBDOWNLOAD.out.db.map { dir -> [ 'bakta_db', '08', '', "bakta/${dir.name}".toString(), dir ] })
            .mix(AMRFINDERPLUS_UPDATE.out.db.map { tgz -> [ 'amrfinderplus_db', '09', '', "amrfinderplus/${tgz.name}".toString(), tgz ] })
        DATABASE_ROW(ch_downloaded)
        PREPARE_MEDAKA_MODEL(channel.value(params.medaka_model))
        PREPARE_MLST_SCHEMES()
        PREPARE_PANEL_REFERENCES(file(params.panel_references, checkIfExists: true))
        ch_rows = ch_rows.mix(DATABASE_ROW.out.row, PREPARE_MEDAKA_MODEL.out.row,
                              PREPARE_MLST_SCHEMES.out.row, PREPARE_PANEL_REFERENCES.out.row)
        if (!params.skip_mobsuite_download) {
            PREPARE_MOBSUITE_DB()
            ch_rows = ch_rows.mix(PREPARE_MOBSUITE_DB.out.row)
        }
    }

    def existing = file("${db_root}/database_manifest.tsv")
    DATABASE_MANIFEST(ch_rows.collect(), channel.value(db_root), existing.exists() ? existing : [])

    emit:
    manifest = DATABASE_MANIFEST.out.manifest
    versions = ch_versions
}
