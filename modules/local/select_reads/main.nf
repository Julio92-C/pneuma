//
// Reads recruited for the local assembly: ids of reads aligned over >= --context_recruit_min_aln bp in a minimap2
// PAF (column 11 = alignment block length), extracted from the sample's reads. No read recruited -> no FASTQ.
//
process SELECT_READS {
    tag "${meta.id}:r${round}"
    label 'process_low'
    container 'quay.io/biocontainers/seqkit:2.13.0--he881be0_0'

    input:
    tuple val(meta), path(reads), path(paf)
    val round

    output:
    tuple val(meta), path("${meta.id}.r${round}.fastq.gz")  , emit: reads, optional: true
    tuple val(meta), path("${meta.id}.r${round}.recruit.tsv"), emit: counts
    tuple val("${task.process}"), val('seqkit'), eval("seqkit version | sed 's/seqkit v//'"), topic: versions, emit: versions_seqkit

    script:
    def out = "${meta.id}.r${round}"
    def total = round == 1 ? "\$(seqkit stats -T ${reads} | awk 'NR == 2 {print \$4}')" : 'NA'
    """
    awk -v m=${params.context_recruit_min_aln} '\$11 >= m {print \$1}' ${paf} | sort -u > ids.txt
    if [ -s ids.txt ]; then
        seqkit grep -j ${task.cpus} -f ids.txt ${reads} -o ${out}.fastq.gz
    fi
    printf 'sample\tround\trecruited\ttotal_reads\n%s\t%s\t%s\t%s\n' ${meta.id} ${round} \$(wc -l < ids.txt) ${total} > ${out}.recruit.tsv
    """

    stub:
    def out = "${meta.id}.r${round}"
    """
    echo '@stub' | gzip > ${out}.fastq.gz
    printf 'sample\tround\trecruited\ttotal_reads\n${meta.id}\t${round}\t1\t1\n' > ${out}.recruit.tsv
    """
}
