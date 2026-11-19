include { MINIMAP2_ALIGN as RECRUIT_R1        } from '../../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as RECRUIT_R2        } from '../../../modules/nf-core/minimap2/align/main'
include { MINIMAP2_ALIGN as MAP_TO_CONTIGS    } from '../../../modules/nf-core/minimap2/align/main'
include { SELECT_READS as SELECT_R1           } from '../../../modules/local/select_reads/main'
include { SELECT_READS as SELECT_R2           } from '../../../modules/local/select_reads/main'
include { FLYE as FLYE_CONTEXT                } from '../../../modules/nf-core/flye/main'
include { CONTEXT_CONTIGS                     } from '../../../modules/local/context_contigs/main'
include { ABRICATE_CONTIGS                    } from '../../../modules/local/abricate_contigs/main'
include { GENOMAD_ENDTOEND as GENOMAD_CONTIGS } from '../../../modules/local/genomad_endtoend/main'
include { MOBSUITE_TYPER                      } from '../../../modules/local/mobsuite_typer/main'
include { ARG_CONTEXT as ARG_CONTEXT_CONTIGS  } from '../../../modules/local/arg_context/main'
include { EXTEND_CONTEXT                      } from '../../../modules/local/extend_context/main'

//
// [ meta, a ] and [ meta, b ] joined on meta.id -> [ meta, a, b ]
//
def joinById(ch_a, ch_b) {
    return ch_a.map { meta, a -> [ meta.id, meta, a ] }
        .join(ch_b.map { meta, b -> [ meta.id, b ] })
        .map { _id, meta, a, b -> [ meta, a, b ] }
}

//
// Warn about samples whose process failed (errorStrategy 'ignore'): attempted ids vs ids that produced output.
//
def warnFailedSamples(ch_attempted, ch_done, label, fallback) {
    def ch_a = ch_attempted.map { meta, _f -> meta.id }.collect().ifEmpty([]).map { ids -> ids.sort().join(',') }
    def ch_d = ch_done.map { meta, _f -> meta.id }.collect().ifEmpty([]).map { ids -> ids.sort().join(',') }
    ch_a.combine(ch_d).subscribe { attempted, done ->
        def att  = attempted ? attempted.tokenize(',') : []
        def ok   = done ? done.tokenize(',') : []
        def lost = att.findAll { id -> !(id in ok) }
        if (att && !ok) {
            log.warn "PNEUMA: ${label} failed for every attempted sample (${lost.join(', ')}); ${fallback}. A configuration problem (--flye_mode, --medaka_model_path, memory) is more likely than biology: check the task logs."
        } else if (lost) {
            log.warn "PNEUMA: ${label} failed for ${lost.size()} of ${att.size()} sample(s): ${lost.join(', ')}; ${fallback}."
        }
    }
}

//
// Layer 2 of the genomic-context analysis: per sample, reads overlapping a hit read (round 1) and reads overlapping
// those (round 2, walks out along a plasmid backbone) are assembled with Flye --meta; ABRicate, geNomad and
// MOB-suite on the pooled contigs give contig calls, which EXTEND_CONTEXT carries back to each read hit.
// Samples without hit reads are skipped; every pooled step runs on empty input (header-only tables).
//
workflow CONTEXT_ASSEMBLY {
    take:
    ch_reads       // [ meta, clean reads.fastq.gz ]
    ch_hit_reads   // [ meta, <sample>.seed_reads.fa ]   (HIT_READS seeds: no CARD-variant-only reads)
    ch_read_loci   // [ 'reads', db, loci_reads_<db>.tsv ]
    abricate_db    // string or null (bundled ABRicate databases)
    genomad_db     // string
    mobsuite_db    // string or null (no MOB-suite typing of the contigs)

    main:
    def ch_r1 = joinById(ch_reads, ch_hit_reads.filter { _meta, fa -> fa.size() > 0 })
        .multiMap { meta, reads, fa -> reads: [ meta, reads ]; target: [ meta, fa ] }
    RECRUIT_R1(ch_r1.reads, ch_r1.target, false, '', false, false)
    SELECT_R1(joinById(ch_reads, RECRUIT_R1.out.paf), 1)

    def ch_r2 = joinById(ch_reads, SELECT_R1.out.reads)
        .multiMap { meta, reads, r1 -> reads: [ meta, reads ]; target: [ meta, r1 ] }
    RECRUIT_R2(ch_r2.reads, ch_r2.target, false, '', false, false)
    SELECT_R2(joinById(ch_reads, RECRUIT_R2.out.paf), 2)

    FLYE_CONTEXT(SELECT_R2.out.reads, params.flye_mode)
    warnFailedSamples(SELECT_R2.out.reads, FLYE_CONTEXT.out.fasta, 'local assembly (Flye)', 'their read-level context calls were used')
    CONTEXT_CONTIGS(joinById(FLYE_CONTEXT.out.fasta, FLYE_CONTEXT.out.txt))

    def ch_abricate_db = abricate_db ? channel.value(file(abricate_db, checkIfExists: true, type: 'dir')) : channel.value([])
    ABRICATE_CONTIGS(CONTEXT_CONTIGS.out.fasta.map { _meta, f -> f }.collect().ifEmpty([]), ch_abricate_db)
    GENOMAD_CONTIGS(ABRICATE_CONTIGS.out.hit_fasta.map { f -> [ [id: 'hit_contigs'], f ] }, file(genomad_db, checkIfExists: true))

    def ch_contig_in = ABRICATE_CONTIGS.out.hits
        .combine(ABRICATE_CONTIGS.out.lengths)
        .combine(ABRICATE_CONTIGS.out.hit_fasta)
        .combine(GENOMAD_CONTIGS.out.scores.map { _m, f -> f })
        .combine(GENOMAD_CONTIGS.out.virus.map { _m, f -> f })
        .combine(GENOMAD_CONTIGS.out.genes.map { _m, f -> f })
    def ch_mob = channel.value([])
    if (mobsuite_db) {
        MOBSUITE_TYPER(ABRICATE_CONTIGS.out.hit_fasta, file(mobsuite_db, checkIfExists: true))
        ch_mob = MOBSUITE_TYPER.out.results
        ch_contig_in = ch_contig_in.combine(ch_mob)
    } else {
        ch_contig_in = ch_contig_in.map { t -> t + [ [] ] }
    }
    def ch_dbs = channel.fromList(params.context_databases.tokenize(','))
    ARG_CONTEXT_CONTIGS(ch_dbs.combine(ch_contig_in).map { db, hits, lengths, fa, scores, virus, genes, mob ->
        [ 'contigs', db, hits, lengths, fa, scores, virus, genes, mob ] })

    def ch_map = joinById(SELECT_R2.out.reads, CONTEXT_CONTIGS.out.fasta.filter { _meta, f -> f.size() > 0 })
        .multiMap { meta, reads, contigs -> reads: [ meta, reads ]; ref: [ meta, contigs ] }
    MAP_TO_CONTIGS(ch_map.reads, ch_map.ref, false, '', false, false)
    def ch_pafs = MAP_TO_CONTIGS.out.paf.map { _meta, f -> f }.collect().ifEmpty([])

    EXTEND_CONTEXT(
        ch_read_loci.map { _route, db, f -> [ db, f ] }
            .join(ARG_CONTEXT_CONTIGS.out.loci.map { _route, db, f -> [ db, f ] }),
        ch_pafs)

    emit:
    contig_loci     = ARG_CONTEXT_CONTIGS.out.loci         // [ 'contigs', db, loci_contigs_<db>.tsv ]
    hit_context     = EXTEND_CONTEXT.out.context           // [ db, hit_context_<db>.tsv ]
    recruited       = SELECT_R2.out.reads                  // [ meta, <id>.r2.fastq.gz ]
    recruit_counts  = SELECT_R1.out.counts.mix(SELECT_R2.out.counts)
    info            = CONTEXT_CONTIGS.out.info             // [ meta, <id>.assembly_info.tsv ]
    contig_fasta    = ABRICATE_CONTIGS.out.fasta
    contig_hits     = ABRICATE_CONTIGS.out.hits
    hit_contigs     = ABRICATE_CONTIGS.out.hit_fasta
    genomad_scores  = GENOMAD_CONTIGS.out.scores.map { _m, f -> f }
    genomad_plasmid = GENOMAD_CONTIGS.out.plasmids.map { _m, f -> f }
    mobtyper        = ch_mob                               // hit_contigs.mobtyper.txt or []
    pafs            = ch_pafs                              // [ <id>.paf, ... ] or []
}
