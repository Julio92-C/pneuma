include { CIRCULAR_PLASMIDS                   } from '../../../modules/local/circular_plasmids/main'
include { MEDAKA_CONSENSUS as MEDAKA_PLASMIDS } from '../../../modules/local/medaka_consensus/main'
include { ISESCAN                             } from '../../../modules/local/isescan/main'
include { IS_FLANKS                           } from '../../../modules/local/is_flanks/main'
include { PLASMID_UNITS                       } from '../../../modules/local/plasmid_units/main'

//
// Plasmids in their own right: circular local-assembly contigs carrying a hit become plasmid drafts (Medaka-
// polished with the sample's recruited reads), ISEScan on the hit contigs gives the IS elements next to each ARG
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

    main:
    ISESCAN(hit_contigs)
    CIRCULAR_PLASMIDS(ch_info.map { _meta, f -> f }.collect().ifEmpty([]), contig_hits, contig_fasta, genomad_scores,
                      genomad_plasmids, mobtyper, ch_pafs, hit_taxa.ifEmpty([]))
    IS_FLANKS(ISESCAN.out.tsv, contig_hits)
    PLASMID_UNITS(CIRCULAR_PLASMIDS.out.tsv, CIRCULAR_PLASMIDS.out.fasta, ISESCAN.out.tsv,
                  ch_recruit_counts.map { _meta, f -> f }.collect().ifEmpty([]))

    if (params.polish_plasmids && medaka_model) {
        def ch_circles = CIRCULAR_PLASMIDS.out.per_sample.flatten().map { f -> [ f.name.replaceAll(/\.circles\.fa$/, ''), f ] }
        def ch_polish = ch_recruited.map { meta, reads -> [ meta.id, meta, reads ] }
            .join(ch_circles)
            .map { _id, meta, reads, circles -> [ meta, reads, circles ] }
        MEDAKA_PLASMIDS(ch_polish, file(medaka_model, checkIfExists: true))
    }

    emit:
    candidates = CIRCULAR_PLASMIDS.out.tsv
    units      = PLASMID_UNITS.out.units
    members    = PLASMID_UNITS.out.members
    is_flanks  = IS_FLANKS.out.tsv
}
