//
// Functions shared by the pathogen subworkflows (no processes here).
//

//
// Meta of one reconstruction unit (a sample x species, or a pooled patient x species).
// panelRow: one row of panel_resolved.tsv (or panelOf(unitMeta)).
//
def unitMeta(Map panelRow, String level, String sample, String patient, List samples, Boolean isControl) {
    def id = level == 'sample' ? "${sample}.${panelRow.slug}" : "patient-${patient}.${panelRow.slug}"
    return [
        id                    : id.toString(),
        level                 : level,
        sample                : sample,
        patient               : patient,
        samples               : samples,
        is_control            : isControl,
        species               : panelRow.species,
        slug                  : panelRow.slug,
        taxid                 : panelRow.taxid,
        reference             : panelRow.reference_accession,
        confusable_accessions : panelRow.confusable_accessions ? panelRow.confusable_accessions.tokenize(';') : [],
        ani_min               : panelRow.ani_min,
        amr_organism          : panelRow.amrfinder_organism ?: '',
        single_end            : true,
    ]
}

//
// The panel fields of a unit meta, in panel_resolved.tsv column names (for unitMeta).
//
def panelOf(Map meta) {
    return [
        species               : meta.species,
        slug                  : meta.slug,
        taxid                 : meta.taxid,
        reference_accession   : meta.reference,
        confusable_accessions : meta.confusable_accessions.join(';'),
        ani_min               : meta.ani_min,
        amrfinder_organism    : meta.amr_organism,
    ]
}

//
// Patient-level fallback: given the gated sample-level metas of one patient x species, return the
// contributing samples (sorted) when NO sample passed the gate and at least two non-control
// samples have extracted reads; otherwise [] (no pooling).
//
def poolingDecision(String patient, List metas) {
    if (!patient) {
        return []
    }
    if (metas.any { m -> m.gate_pass }) {
        return []
    }
    def contributing = metas.findAll { m -> !m.is_control }.collect { m -> m.sample }.unique().sort()
    return contributing.size() >= 2 ? contributing : []
}

//
// Value of `column` in the first data row of a 2-line TSV written by a PNEUMA helper.
//
def tsvField(tsv, String column) {
    def lines = tsv.readLines()
    def header = lines[0].split('\t') as List
    def idx = header.indexOf(column)
    if (idx < 0 || lines.size() < 2) {
        error "PATHOGENS: ${tsv.name} has no '${column}' value"
    }
    return lines[1].split('\t', -1)[idx]
}

//
// skani comparison set of a unit: its panel reference plus the confusable references.
//
def referenceFiles(Map meta, refdir) {
    return ([ meta.reference ] + meta.confusable_accessions).unique().collect { acc -> refdir.resolve("${acc}.fna.gz") }
}

//
// Placeholder so `nf-core pipelines lint` (which expects a workflow block in every local
// subworkflow) does not crash on this functions-only file. It is never called.
//
workflow UTILS_PATHOGENS {
}
