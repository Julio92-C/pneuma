include { TAG_READS      } from '../../../modules/local/tag_reads/main'
include { ABRICATE_READS } from '../../../modules/local/abricate_reads/main'
include { ABRICATE_MERGE } from '../../../modules/local/abricate_merge/main'

workflow FUNCTIONAL_READS {
    take:
    ch_classified   // [ meta, classified.fastq.gz ]

    main:
    TAG_READS(ch_classified)
    ABRICATE_READS(TAG_READS.out.fasta)
    ABRICATE_MERGE(ABRICATE_READS.out.csv.map { _meta, csv -> csv }.collect())

    emit:
    summary    = ABRICATE_MERGE.out.summary
    per_sample = ABRICATE_READS.out.csv
    versions   = TAG_READS.out.versions.first().mix(ABRICATE_READS.out.versions.first(), ABRICATE_MERGE.out.versions)
}
