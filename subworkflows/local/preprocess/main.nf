include { CAT_FASTQ     } from '../../../modules/nf-core/cat/fastq/main'
include { PORECHOP_ABI  } from '../../../modules/nf-core/porechop/abi/main'
include { NANOQ         } from '../../../modules/nf-core/nanoq/main'
include { HOSTILE_CLEAN } from '../../../modules/local/hostile_clean/main'
include { NANOPLOT      } from '../../../modules/nf-core/nanoplot/main'
include { SEQKIT_STATS  } from '../../../modules/nf-core/seqkit/stats/main'

workflow PREPROCESS {
    take:
    ch_reads      // [ meta, [fastq...] ]
    host_index    // path: fasta or .mmi

    main:
    ch_versions = channel.empty()
    ch_multiqc  = channel.empty()

    ch_reads
        .branch { meta, fastqs ->
            multi : fastqs.size() > 1
            single: true
                return [ meta, fastqs[0] ]
        }
        .set { ch_split }

    CAT_FASTQ(ch_split.multi)
    ch_merged   = CAT_FASTQ.out.reads.mix(ch_split.single)

    if (params.run_porechop) {
        PORECHOP_ABI(ch_merged, [])
        ch_merged   = PORECHOP_ABI.out.reads
        ch_multiqc  = ch_multiqc.mix(PORECHOP_ABI.out.log.map { _meta, f -> f })
    }

    NANOQ(ch_merged, 'fastq.gz')

    HOSTILE_CLEAN(NANOQ.out.reads, host_index)
    ch_versions = ch_versions.mix(HOSTILE_CLEAN.out.versions.first())

    NANOPLOT(HOSTILE_CLEAN.out.reads)
    SEQKIT_STATS(HOSTILE_CLEAN.out.reads)
    ch_multiqc  = ch_multiqc.mix(NANOPLOT.out.txt.collect { _meta, f -> f }, NANOQ.out.stats.collect { _meta, f -> f },
                                 SEQKIT_STATS.out.stats.collect { _meta, f -> f })

    // seqkit stats -T: header + one row; column 3 (0-based) = num_seqs
    ch_counts = SEQKIT_STATS.out.stats.map { meta, tsv ->
        [ meta.id, tsv.readLines()[1].split('\t')[3].replace(',', '') as Long ]
    }

    HOSTILE_CLEAN.out.reads
        .map { meta, fq -> [ meta.id, meta, fq ] }
        .join(ch_counts)
        .map { id, meta, fq, n ->
            if (n < params.min_reads) {
                log.warn "Sample ${id}: ${n} reads after host removal (< --min_reads ${params.min_reads}); flagged LOW_DEPTH"
            }
            [ meta + [ n_reads: n, low_depth: n < params.min_reads ], fq ]
        }
        .set { ch_clean }

    emit:
    reads    = ch_clean
    multiqc  = ch_multiqc
    hostile_log = HOSTILE_CLEAN.out.log.map { _meta, f -> f }   // Hostile JSON logs, for the MultiQC tables
    versions = ch_versions
}
