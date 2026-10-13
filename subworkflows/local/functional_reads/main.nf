include { TAG_READS      } from '../../../modules/local/tag_reads/main'
include { ABRICATE_READS } from '../../../modules/local/abricate_reads/main'
include { ABRICATE_MERGE } from '../../../modules/local/abricate_merge/main'

workflow FUNCTIONAL_READS {
    take:
    ch_classified   // [ meta, classified.fastq.gz ]
    abricate_db     // string: directory of current ABRicate databases, or null = databases bundled in the container

    main:
    TAG_READS(ch_classified)
    ch_abricate_db = abricate_db ? channel.value(file(abricate_db, checkIfExists: true, type: 'dir')) : channel.value([])
    ABRICATE_READS(TAG_READS.out.fasta, ch_abricate_db)
    ABRICATE_MERGE(ABRICATE_READS.out.csv.map { _meta, csv -> csv }.collect())

    emit:
    summary    = ABRICATE_MERGE.out.summary
    per_sample = ABRICATE_READS.out.csv
    versions   = TAG_READS.out.versions.first().mix(ABRICATE_READS.out.versions.first(), ABRICATE_MERGE.out.versions)
}
