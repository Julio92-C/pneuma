include { HIT_READS        } from '../../../modules/local/hit_reads/main'
include { GENOMAD_ENDTOEND } from '../../../modules/local/genomad_endtoend/main'
include { ARG_CONTEXT      } from '../../../modules/local/arg_context/main'
include { CONTEXT_ASSEMBLY } from '../context_assembly/main'
include { MOBILOME         } from '../mobilome/main'
include { CONTEXT_SUMMARY  } from '../../../modules/local/context_summary/main'

//
// Layer 1 of the genomic-context analysis: every ABRicate hit on a Kraken2-classified read is placed as plasmid /
// chromosomal_mge / chromosome / unresolved from geNomad on the pooled hit reads, with the plasmid-evidence rule
// (a hallmark or conjugation gene in geNomad's genes table, or a PlasmidFinder / MOB-suite replicon on the read).
// Layer 2 (`--context_local_assembly`) and the mobilome branch (`--mobilome`) follow in CONTEXT_ASSEMBLY / MOBILOME.
// CONTEXT_SUMMARY then combines the calls of all databases into context/hit_context.tsv (with the read's Kraken2 taxon) and the summaries.
//
workflow GENOMIC_CONTEXT {
    take:
    ch_tagged       // [ meta, <sample>_kraken2ID_PFP.fasta ]   (FUNCTIONAL_READS.out.tagged)
    ch_abricate     // [ meta, <sample>_abricate.csv ]          (FUNCTIONAL_READS.out.per_sample)
    ch_reads        // [ meta, clean reads.fastq.gz ]           (PREPROCESS.out.reads; local assembly)
    genomad_db      // string: geNomad database directory
    abricate_db     // string or null: ABRicate databases (null = bundled)
    mobsuite_db     // string or null: MOB-suite database (required by the mobilome)
    medaka_model    // string or null: Medaka model (plasmid draft polishing)
    taxdump         // string: NCBI taxdump directory (taxon names)

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

    def ch_contig_loci = channel.empty()
    def ch_hit_context = channel.empty()
    def ch_candidates  = channel.empty()
    def ch_units       = channel.empty()
    def ch_is_flanks   = channel.empty()
    if (params.context_local_assembly) {
        CONTEXT_ASSEMBLY(ch_reads, HIT_READS.out.seeds, ARG_CONTEXT.out.loci, abricate_db, genomad_db, mobsuite_db)
        ch_contig_loci = CONTEXT_ASSEMBLY.out.contig_loci
        ch_hit_context = CONTEXT_ASSEMBLY.out.hit_context
        if (params.mobilome) {
            MOBILOME(CONTEXT_ASSEMBLY.out.info, CONTEXT_ASSEMBLY.out.contig_hits, CONTEXT_ASSEMBLY.out.contig_fasta,
                     CONTEXT_ASSEMBLY.out.hit_contigs, CONTEXT_ASSEMBLY.out.genomad_scores,
                     CONTEXT_ASSEMBLY.out.genomad_plasmid, CONTEXT_ASSEMBLY.out.mobtyper, CONTEXT_ASSEMBLY.out.pafs,
                     ch_taxa, CONTEXT_ASSEMBLY.out.recruited, CONTEXT_ASSEMBLY.out.recruit_counts, medaka_model, taxdump)
            ch_candidates = MOBILOME.out.candidates
            ch_units      = MOBILOME.out.units
            ch_is_flanks  = MOBILOME.out.is_flanks
        }
    }

    def ch_tables = params.context_local_assembly
        ? ch_hit_context.map { _db, f -> f }
        : ARG_CONTEXT.out.loci.map { _route, _db, f -> f }
    CONTEXT_SUMMARY(ch_tables.collect().ifEmpty([]), ch_taxa.ifEmpty([]), file(taxdump, checkIfExists: true, type: 'dir'),
                    params.meridian_metadata ? file(params.meridian_metadata, checkIfExists: true) : [])

    emit:
    loci       = ARG_CONTEXT.out.loci            // [ 'reads', db, loci_reads_<db>.tsv ]
    score_only = ARG_CONTEXT.out.score_only
    hit_reads  = ch_fa                           // pooled hit reads (geNomad input; the local assembly is seeded by the per-sample HIT_READS.out.seeds)
    taxa       = ch_taxa                         // read id -> Kraken2 taxid
    contig_loci        = ch_contig_loci    // [ 'contigs', db, loci_contigs_<db>.tsv ]
    hit_context        = ch_hit_context    // [ db, hit_context_<db>.tsv ]
    plasmid_candidates = ch_candidates
    plasmid_units      = ch_units
    is_flanks          = ch_is_flanks
    hit_table          = CONTEXT_SUMMARY.out.hit_context   // context/hit_context.tsv
    summary            = CONTEXT_SUMMARY.out.summary       // [ summary/*.tsv ]
    genomad    = GENOMAD_ENDTOEND.out.scores.map { _m, f -> f }
                     .combine(GENOMAD_ENDTOEND.out.virus.map { _m, f -> f })
                     .combine(GENOMAD_ENDTOEND.out.plasmids.map { _m, f -> f })
                     .combine(GENOMAD_ENDTOEND.out.genes.map { _m, f -> f })
}
