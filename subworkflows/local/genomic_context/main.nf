include { HIT_READS        } from '../../../modules/local/hit_reads/main'
include { GENOMAD_ENDTOEND } from '../../../modules/local/genomad_endtoend/main'
include { ARG_CONTEXT      } from '../../../modules/local/arg_context/main'

//
// Layer 1 of the genomic-context analysis: every ABRicate hit on a Kraken2-classified read is placed as plasmid /
// chromosomal_mge / chromosome / unresolved from geNomad on the pooled hit reads, with the plasmid-evidence rule
// (a hallmark or conjugation gene in geNomad's genes table, or a PlasmidFinder / MOB-suite replicon on the read).
//
workflow GENOMIC_CONTEXT {
    take:
    ch_tagged       // [ meta, <sample>_kraken2ID_PFP.fasta ]   (FUNCTIONAL_READS.out.tagged)
    ch_abricate     // [ meta, <sample>_abricate.csv ]          (FUNCTIONAL_READS.out.per_sample)
    genomad_db      // string: geNomad database directory

    main:
    HIT_READS(ch_tagged.join(ch_abricate))
    def ch_fa      = HIT_READS.out.fasta.map { _meta, f -> f }.collectFile(name: 'hit_reads.fa', sort: true)
    def ch_hits    = HIT_READS.out.hits.map { _meta, f -> f }.collectFile(name: 'hits.tsv', sort: true)
    def ch_lengths = HIT_READS.out.lengths.map { _meta, f -> f }.collectFile(name: 'hit_reads.lengths.tsv', sort: true)
    def ch_taxa    = HIT_READS.out.taxa.map { _meta, f -> f }.collectFile(name: 'hit_reads.taxa.tsv', sort: true)

    GENOMAD_ENDTOEND(ch_fa.map { f -> [ [id: 'hit_reads'], f ] }, file(genomad_db, checkIfExists: true))

    def ch_inputs = ch_hits
        .combine(ch_lengths)
        .combine(ch_fa)
        .combine(GENOMAD_ENDTOEND.out.scores.map { _meta, f -> f })
        .combine(GENOMAD_ENDTOEND.out.virus.map { _meta, f -> f })
        .combine(GENOMAD_ENDTOEND.out.genes.map { _meta, f -> f })
    def ch_dbs = channel.fromList(params.context_databases.tokenize(','))
    ARG_CONTEXT(ch_dbs.combine(ch_inputs).map { db, hits, lengths, fa, scores, virus, genes ->
        [ 'reads', db, hits, lengths, fa, scores, virus, genes, [] ] })

    emit:
    loci       = ARG_CONTEXT.out.loci            // [ 'reads', db, loci_reads_<db>.tsv ]
    score_only = ARG_CONTEXT.out.score_only
    hit_reads  = ch_fa                           // pooled hit reads (Plan B: local assembly seeds)
    taxa       = ch_taxa                         // read id -> Kraken2 taxid
    genomad    = GENOMAD_ENDTOEND.out.scores.map { _m, f -> f }
                     .combine(GENOMAD_ENDTOEND.out.virus.map { _m, f -> f })
                     .combine(GENOMAD_ENDTOEND.out.plasmids.map { _m, f -> f })
                     .combine(GENOMAD_ENDTOEND.out.genes.map { _m, f -> f })
}
