include { FLYE                  } from '../../../modules/nf-core/flye/main'
include { GUNZIP as GUNZIP_FLYE } from '../../../modules/nf-core/gunzip/main'
include { MEDAKA_CONSENSUS      } from '../../../modules/local/medaka_consensus/main'
include { CHECKM2_PREDICT       } from '../../../modules/nf-core/checkm2/predict/main'
include { QUAST                 } from '../../../modules/nf-core/quast/main'

//
// Flye -> Medaka per unit that passed the coverage gate, then ONE CheckM2 and ONE QUAST run over
// all assemblies (one database load). A unit whose Flye/Medaka task failed (errorStrategy: memory
// retries, then 'ignore' for these two processes only) has no assembly -> ASSEMBLY_FAILED.
//
workflow PATHOGEN_ASSEMBLY {
    take:
    ch_reads       // [ meta, reads.fastq.gz ]  units with meta.gate_pass
    medaka_model   // path: <model>_model_pt.tar.gz
    checkm2_db     // path: CheckM2 .dmnd

    main:
    FLYE(ch_reads, params.flye_mode)
    GUNZIP_FLYE(FLYE.out.fasta)
    def ch_polish = ch_reads.map { meta, reads -> [ meta.id, meta, reads ] }
        .join(GUNZIP_FLYE.out.gunzip.map { meta, draft -> [ meta.id, draft ] })
        .map { _id, meta, reads, draft -> [ meta, reads, draft ] }
    MEDAKA_CONSENSUS(ch_polish, medaka_model)

    def ch_all = MEDAKA_CONSENSUS.out.assembly.map { _meta, fasta -> fasta }.collect()
    CHECKM2_PREDICT(ch_all.map { fastas -> [ [ id: 'pathogens' ], fastas ] }, [ [ id: 'checkm2_db' ], checkm2_db ])
    QUAST(ch_all.map { fastas -> [ [ id: 'pathogens_quast' ], fastas ] }, [ [:], [] ], [ [:], [] ])

    emit:
    assembly    = MEDAKA_CONSENSUS.out.assembly    // [ meta, <id>.fasta ]
    checkm2_tsv = CHECKM2_PREDICT.out.checkm2_tsv  // [ [id:'pathogens'], pathogens_checkm2_report.tsv ]
    quast_tsv   = QUAST.out.tsv                    // [ [id:'pathogens_quast'], pathogens_quast.tsv ]
    flye_log    = FLYE.out.log
}
