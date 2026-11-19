include { CIRCULAR_PLASMIDS                   } from '../../../modules/local/circular_plasmids/main'
include { MEDAKA_CONSENSUS as MEDAKA_PLASMIDS } from '../../../modules/local/medaka_consensus/main'
include { POLISHED_CIRCLES                    } from '../../../modules/local/polished_circles/main'
include { ISESCAN                             } from '../../../modules/local/isescan/main'
include { IS_FLANKS                           } from '../../../modules/local/is_flanks/main'
include { PLASMID_UNITS                       } from '../../../modules/local/plasmid_units/main'
include { warnFailedSamples                   } from '../utils_context/main'

//
// Plasmids in their own right: circular local-assembly contigs carrying a hit become plasmid drafts, polished with
// Medaka against all of the sample's contigs (reads of a chromosomal copy of an IS element or gene then align to the
// chromosome, not to the plasmid) and cut out again; ISEScan on the hit contigs gives the IS elements next to each ARG
// (composite transposons), and the drafts of all samples are grouped into plasmid units.
//
workflow MOBILOME {
    take:
    ch_info          // [ meta, <id>.assembly_info.tsv ]
    contig_hits      // contig_hits.tsv
    contig_fasta     // all_contigs.fa
    hit_contigs      // hit_contigs.fa
    genomad_scores   // hit_contigs_aggregated_classification.tsv
    genomad_plasmids // hit_contigs_plasmid_summary.tsv
    mobtyper         // hit_contigs.mobtyper.txt or []
    ch_pafs          // [ <id>.paf, ... ] or []
    hit_taxa         // hit_reads.taxa.tsv
    ch_recruited     // [ meta, <id>.r2.fastq.gz ]
    ch_recruit_counts// [ meta, <id>.r<n>.recruit.tsv ]
    medaka_model     // string or null
    taxdump          // string or null: directory with names.dmp (plasmid-unit read taxa by name)
    ch_contigs       // [ meta, <id>.contigs.fa ]   (the sample's local-assembly contigs: the Medaka draft)

    main:
    ISESCAN(hit_contigs)
    CIRCULAR_PLASMIDS(ch_info.map { _meta, f -> f }.collect().ifEmpty([]), contig_hits, contig_fasta, genomad_scores,
                      genomad_plasmids, mobtyper, ch_pafs, hit_taxa.ifEmpty([]))
    IS_FLANKS(ISESCAN.out.tsv, contig_hits)
    PLASMID_UNITS(CIRCULAR_PLASMIDS.out.tsv, CIRCULAR_PLASMIDS.out.fasta, ISESCAN.out.tsv,
                  ch_recruit_counts.map { _meta, f -> f }.collect().ifEmpty([]),
                  taxdump ? file(taxdump, checkIfExists: true, type: 'dir') : [])

    def ch_polished = channel.empty()
    if (params.polish_plasmids && medaka_model) {
        def ch_circles = CIRCULAR_PLASMIDS.out.per_sample.flatten().map { f -> [ f.name.replaceAll(/\.circles\.fa$/, ''), f ] }
        def ch_polish = ch_recruited.map { meta, reads -> [ meta.id, meta, reads ] }
            .join(ch_contigs.map { meta, contigs -> [ meta.id, contigs ] })
            .join(ch_circles)
            .map { _id, meta, reads, contigs, circles -> [ meta, reads, contigs, circles ] }
        MEDAKA_PLASMIDS(ch_polish.map { meta, reads, contigs, _circles -> [ meta, reads, contigs ] },
                        file(medaka_model, checkIfExists: true))
        POLISHED_CIRCLES(MEDAKA_PLASMIDS.out.assembly.map { meta, f -> [ meta.id, meta, f ] }
            .join(ch_polish.map { meta, _reads, _contigs, circles -> [ meta.id, circles ] })
            .map { _id, meta, polished, circles -> [ meta, polished, circles ] })
        warnFailedSamples(ch_polish.map { meta, _reads, _contigs, circles -> [ meta, circles ] }, MEDAKA_PLASMIDS.out.assembly,
                          'plasmid polishing (Medaka)', 'their unpolished plasmid drafts were used',
                          'Check --medaka_model_path (R9.4.1 reads need an R9.4.1 --medaka_model) and the memory in the task logs.')
        ch_polished = POLISHED_CIRCLES.out.fasta
    }

    emit:
    candidates = CIRCULAR_PLASMIDS.out.tsv
    units      = PLASMID_UNITS.out.units
    members    = PLASMID_UNITS.out.members
    is_flanks  = IS_FLANKS.out.tsv
    polished   = ch_polished          // [ meta, <id>.circles.polished.fasta ]
}
