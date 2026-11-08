include { MINIMAP2_ALIGN as MINIMAP2_GATE } from '../../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_COVERAGE               } from '../../../modules/nf-core/samtools/coverage/main'
include { COVERAGE_GATE                   } from '../../../modules/local/coverage_gate/main'
include { tsvField                        } from '../utils_pathogens/main'

//
// Map a unit's reads to its panel reference and apply the coverage gate (depth and breadth).
// Called twice by PATHOGEN_RECONSTRUCTION: GATE_SAMPLE (per sample) and GATE_PATIENT (pooled).
//
workflow PATHOGEN_GATE {
    take:
    ch_reads   // [ meta (unit), reads.fastq.gz ]
    refdir     // path: directory with <accession>.fna.gz

    main:
    def ch_in = ch_reads.multiMap { meta, reads ->
        reads: [ meta, reads ]
        ref:   [ meta, refdir.resolve("${meta.reference}.fna.gz") ]
    }
    MINIMAP2_GATE(ch_in.reads, ch_in.ref, true, 'bai', false, false)
    SAMTOOLS_COVERAGE(MINIMAP2_GATE.out.bam.join(MINIMAP2_GATE.out.index), [ [:], [], [] ])
    COVERAGE_GATE(SAMTOOLS_COVERAGE.out.coverage)

    def ch_gated = COVERAGE_GATE.out.tsv
        .map { meta, tsv ->
            [ meta.id, meta + [
                reads_mapped : tsvField(tsv, 'reads_mapped') as Long,
                depth        : tsvField(tsv, 'mean_depth') as Double,
                breadth      : tsvField(tsv, 'breadth') as Double,
                gate_pass    : tsvField(tsv, 'pass') == 'true' ] ]
        }
        .join(ch_reads.map { meta, reads -> [ meta.id, reads ] })
        .map { _id, meta, reads -> [ meta, reads ] }

    emit:
    gated = ch_gated               // [ meta + gate fields, reads ]
    tsv   = COVERAGE_GATE.out.tsv  // [ meta, <id>.gate.tsv ]
}
