include { MLST                                     } from '../../../modules/nf-core/mlst/main'
include { AMRFINDERPLUS_RUN                        } from '../../../modules/nf-core/amrfinderplus/run/main'
include { ABRICATE_RUN as ABRICATE_VFDB            } from '../../../modules/nf-core/abricate/run/main'
include { ABRICATE_RUN as ABRICATE_PLASMIDFINDER   } from '../../../modules/nf-core/abricate/run/main'
include { BAKTA_BAKTA                              } from '../../../modules/nf-core/bakta/bakta/main'
include { MOBSUITE_RECON                           } from '../../../modules/local/mobsuite_recon/main'

//
// Typing and annotation of confirmed assemblies (and, at patient level, single-strain ones):
// mlst (PubMLST schemes bundled in the container), AMRFinderPlus (--organism when the panel maps
// one), ABRicate VFDB + PlasmidFinder on contigs, Bakta, optional MOB-suite.
//
workflow PATHOGEN_TYPING {
    take:
    ch_assembly   // [ meta, <id>.fasta ]
    databases     // map from resolveDatabases()

    main:
    MLST(ch_assembly)
    // nf-core AMRFINDERPLUS_RUN adds --organism when meta has an 'organism' key: set it only when mapped
    AMRFINDERPLUS_RUN(
        ch_assembly.map { meta, fasta -> [ meta.amr_organism ? meta + [ organism: meta.amr_organism ] : meta, fasta ] },
        file(databases.amrfinderplus_db, checkIfExists: true)
    )
    def abricate_db = databases.abricate_db ? file(databases.abricate_db, checkIfExists: true, type: 'dir') : []
    ABRICATE_VFDB(ch_assembly, abricate_db)
    ABRICATE_PLASMIDFINDER(ch_assembly, abricate_db)
    BAKTA_BAKTA(ch_assembly, file(databases.bakta_db, checkIfExists: true, type: 'dir'), [], [], [], [])

    def ch_mob = channel.empty()
    if (params.run_mobsuite) {
        MOBSUITE_RECON(ch_assembly, file(databases.mobsuite_db, checkIfExists: true, type: 'dir'))
        ch_mob = MOBSUITE_RECON.out.mobtyper
    }

    def ch_files = MLST.out.tsv.map { _meta, f -> [ 'mlst', f ] }
        .mix(AMRFINDERPLUS_RUN.out.report.map { _meta, f -> [ 'amrfinderplus', f ] })
        .mix(ABRICATE_VFDB.out.report.map { _meta, f -> [ 'abricate', f ] })
        .mix(ABRICATE_PLASMIDFINDER.out.report.map { _meta, f -> [ 'abricate', f ] })
        .mix(ch_mob.map { _meta, f -> [ 'mobsuite', f ] })

    emit:
    files = ch_files            // [ kind, file ] for the pathogen report
    gff   = BAKTA_BAKTA.out.gff // [ meta, <id>.gff3 ]
}
