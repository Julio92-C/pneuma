//
// Functions shared by the genomic-context subworkflows (no processes here).
//

//
// Warning text for the samples whose process (errorStrategy 'ignore') produced no output, or null when every attempted
// sample produced one. attempted / done: sample ids. hint: the likely causes for this step (small read sets are the
// usual reason a local assembly fails, so it is not presented as a configuration problem).
//
def failedSamplesMessage(List attempted, List done, String label, String fallback, String hint) {
    def tried = attempted.unique(false).sort(false)
    def lost = tried.findAll { id -> !(id in done) }
    if (!lost) {
        return null
    }
    def which = lost.size() == tried.size() ? "every attempted sample (${lost.join(', ')})" :
                "${lost.size()} of ${tried.size()} sample(s): ${lost.join(', ')}"
    return "PNEUMA: ${label} produced no output for ${which}; ${fallback}. ${hint}".toString()
}

//
// Log failedSamplesMessage once both channels are complete. ch_attempted / ch_done: [ meta, file ]. Sample ids cannot
// contain ',' (samplesheet schema), so the id lists travel as joined strings (a combined list would be flattened).
//
def warnFailedSamples(ch_attempted, ch_done, String label, String fallback, String hint) {
    def ch_a = ch_attempted.map { meta, _f -> meta.id }.collect().ifEmpty([]).map { ids -> ids.sort().join(',') }
    def ch_d = ch_done.map { meta, _f -> meta.id }.collect().ifEmpty([]).map { ids -> ids.sort().join(',') }
    ch_a.combine(ch_d).subscribe { attempted, done ->
        def msg = failedSamplesMessage(attempted ? attempted.tokenize(',') : [], done ? done.tokenize(',') : [],
                                       label, fallback, hint)
        if (msg) {
            log.warn msg
        }
    }
}

//
// Placeholder so `nf-core pipelines lint` (which expects a workflow block in every local subworkflow) does not crash
// on this functions-only file. It is never called.
//
workflow UTILS_CONTEXT {
}
