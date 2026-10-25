include { RESOLVE_PANEL                     } from '../../../modules/local/resolve_panel/main'
include { PATHOGEN_CANDIDATES               } from '../../../modules/local/pathogen_candidates/main'
include { KRAKENTOOLS_EXTRACTKRAKENREADS    } from '../../../modules/nf-core/krakentools/extractkrakenreads/main'
include { CAT_FASTQ as CAT_POOLED           } from '../../../modules/nf-core/cat/fastq/main'
include { SKANI_DIST                        } from '../../../modules/nf-core/skani/dist/main'
include { SPECIES_CONFIRM                   } from '../../../modules/local/species_confirm/main'
include { MINIMAP2_ALIGN as MINIMAP2_MIXED  } from '../../../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_MPILEUP                  } from '../../../modules/nf-core/samtools/mpileup/main'
include { MIXED_STRAIN_CHECK                } from '../../../modules/local/mixed_strain_check/main'
include { PATHOGEN_GATE as GATE_SAMPLE      } from '../pathogen_gate/main'
include { PATHOGEN_GATE as GATE_PATIENT     } from '../pathogen_gate/main'
include { PATHOGEN_ASSEMBLY                 } from '../pathogen_assembly/main'
include { PATHOGEN_TYPING                   } from '../pathogen_typing/main'
include { pathogenPanelPath; panelReferencesPath } from '../utils_nfcore_pneuma_pipeline/main'
include { unitMeta; panelOf; poolingDecision; tsvField; referenceFiles } from '../utils_pathogens/main'

//
// Spec 5.5 + patient-level fallback. Per sample x reconstruct-species: Kraken2 read extraction
// (--include-children) -> coverage gate -> Flye/Medaka -> CheckM2/QUAST -> skani species
// confirmation -> mlst / AMRFinderPlus / ABRicate / Bakta (/ MOB-suite). Patient x species where
// no sample passed the gate: pooled reads -> gate -> assembly -> species -> mixed-strain check -> typing.
//
workflow PATHOGEN_RECONSTRUCTION {
    take:
    ch_samples      // [ meta ] one per sample: id, patient, is_control, low_depth, n_reads
    ch_assignments  // [ meta, <id>.kraken2.classifiedreads.txt ]
    ch_classified   // [ meta, <id>.classified.fastq.gz ]
    ch_reports      // [ meta, <id>.kraken2.report.txt ]
    ch_bracken      // [ meta, <id>.tsv ] (absent for LOW_DEPTH / no-taxa samples)
    databases       // map from resolveDatabases()

    main:
    def refdir = file(databases.panel_references_dir, checkIfExists: true, type: 'dir')

    RESOLVE_PANEL(
        file(pathogenPanelPath(params), checkIfExists: true),
        file(panelReferencesPath(params), checkIfExists: true),
        file(databases.taxdump, checkIfExists: true, type: 'dir')
    )

    def ch_samples_tsv = channel.of('sample\tpatient\tis_control\tlow_depth\tn_reads')
        .concat(ch_samples.map { meta -> [ meta.id, meta.patient ?: '', meta.is_control, meta.low_depth, meta.n_reads ].join('\t') })
        .collectFile(name: 'samples.tsv', newLine: true, sort: false)

    PATHOGEN_CANDIDATES(
        RESOLVE_PANEL.out.tsv,
        ch_samples_tsv,
        ch_reports.map { _meta, f -> f }.collect(),
        ch_bracken.map { _meta, f -> f }.collect().ifEmpty([])
    )

    def ch_panel = RESOLVE_PANEL.out.tsv
        .splitCsv(sep: '\t', header: true)
        .toList()
        .map { rows -> rows.collectEntries { row -> [ (row.species): row ] } }

    def ch_kraken = ch_assignments.map { meta, f -> [ meta.id, meta, f ] }
        .join(ch_classified.map { meta, f -> [ meta.id, f ] })
        .join(ch_reports.map { meta, f -> [ meta.id, f ] })

    // ---- sample level: extraction units (decision = extract)
    def ch_extract = PATHOGEN_CANDIDATES.out.candidates
        .splitCsv(sep: '\t', header: true)
        .filter { row -> row.decision == 'extract' }
        .map { row -> [ row.sample, row ] }
        .combine(ch_kraken, by: 0)
        .combine(ch_panel)
        .map { _id, row, smeta, txt, fq, report, panel ->
            def meta = unitMeta(panel[row.species], 'sample', smeta.id, smeta.patient ?: '', [ smeta.id ], smeta.is_control as Boolean) +
                       [ clade_reads: row.clade_reads as Long ]
            [ meta, txt, fq, report ]
        }
        .multiMap { meta, txt, fq, report ->
            taxid:      meta.taxid
            assignment: [ meta, txt ]
            reads:      [ meta, fq ]
            report:     [ meta, report ]
        }
    KRAKENTOOLS_EXTRACTKRAKENREADS(ch_extract.taxid, ch_extract.assignment, ch_extract.reads, ch_extract.report)
    GATE_SAMPLE(KRAKENTOOLS_EXTRACTKRAKENREADS.out.extracted_kraken2_reads, refdir)

    // ---- patient level: pool a species' reads when no sample of the patient passed the gate
    def ch_pool = GATE_SAMPLE.out.gated
        .filter { meta, _reads -> meta.patient && !meta.is_control }
        .map { meta, reads -> [ "${meta.patient}\t${meta.species}".toString(), meta, reads ] }
        .groupTuple()
        .flatMap { _key, metas, reads ->
            def contributing = poolingDecision(metas[0].patient, metas)
            if (!contributing) {
                return []
            }
            def idx = (0..<metas.size()).toList().findAll { i -> metas[i].sample in contributing }.sort { i -> metas[i].sample }
            def meta = unitMeta(panelOf(metas[0]), 'patient', '', metas[0].patient, contributing, false) +
                       [ clade_reads: idx.collect { i -> metas[i].clade_reads }.sum() ]
            log.info "PATHOGENS: no sample of patient ${metas[0].patient} passed the gate for ${metas[0].species}; pooling ${contributing.join(', ')}"
            return [ [ meta, idx.collect { i -> reads[i] } ] ]
        }
    CAT_POOLED(ch_pool)
    GATE_PATIENT(CAT_POOLED.out.reads, refdir)

    // ---- assembly + QC (both levels in one call)
    def ch_gated = GATE_SAMPLE.out.gated.mix(GATE_PATIENT.out.gated)
    def ch_pass  = ch_gated.filter { meta, _reads -> meta.gate_pass }
    PATHOGEN_ASSEMBLY(ch_pass, file(databases.medaka_model_path, checkIfExists: true), file(databases.checkm2_db, checkIfExists: true))
    def ch_asm = PATHOGEN_ASSEMBLY.out.assembly

    ch_pass.count()
        .combine(ch_asm.count())
        .subscribe { attempted, assembled ->
            if (attempted > 0 && assembled == 0) {
                log.warn "PATHOGENS: all ${attempted} assembly attempt(s) failed (ASSEMBLY_FAILED). A configuration problem (--medaka_model_path, --flye_mode, memory) is more likely than biology: check the FLYE / MEDAKA_CONSENSUS task logs."
            }
        }

    // ---- species confirmation (skani vs panel reference + confusables)
    def ch_skani = ch_asm.multiMap { meta, fasta ->
        query: [ meta, fasta ]
        refs:  [ [ id: 'refs' ], referenceFiles(meta, refdir) ]
    }
    SKANI_DIST(ch_skani.query, ch_skani.refs)
    SPECIES_CONFIRM(SKANI_DIST.out.dist, file(panelReferencesPath(params), checkIfExists: true))
    def ch_confirmed = SPECIES_CONFIRM.out.tsv
        .map { meta, tsv -> [ meta.id, tsvField(tsv, 'confirmed') == 'true' ] }
        .join(ch_asm.map { meta, fasta -> [ meta.id, meta, fasta ] })
        .filter { _id, confirmed, _meta, _fasta -> confirmed }
        .map { _id, _confirmed, meta, fasta -> [ meta, fasta ] }

    // ---- mixed-strain check of confirmed patient-level assemblies (pooled reads -> pooled assembly)
    def ch_patient_asm = ch_confirmed.filter { meta, _fasta -> meta.level == 'patient' }
    def ch_mix = ch_patient_asm.map { meta, fasta -> [ meta.id, meta, fasta ] }
        .join(GATE_PATIENT.out.gated.map { meta, reads -> [ meta.id, reads ] })
        .multiMap { _id, meta, fasta, reads ->
            reads: [ meta, reads ]
            ref:   [ meta, fasta ]
        }
    MINIMAP2_MIXED(ch_mix.reads, ch_mix.ref, true, 'bai', false, false)
    def ch_pile = MINIMAP2_MIXED.out.bam.join(MINIMAP2_MIXED.out.index)
        .map { meta, bam, bai -> [ meta.id, meta, bam, bai ] }
        .join(ch_patient_asm.map { meta, fasta -> [ meta.id, fasta ] })
        .multiMap { _id, meta, bam, bai, fasta ->
            bam:   [ meta, bam, bai, [] ]
            fasta: [ meta, fasta, [] ]
        }
    SAMTOOLS_MPILEUP(ch_pile.bam, ch_pile.fasta)
    MIXED_STRAIN_CHECK(SAMTOOLS_MPILEUP.out.mpileup)
    def ch_single_strain = MIXED_STRAIN_CHECK.out.tsv
        .map { meta, tsv -> [ meta.id, tsvField(tsv, 'verdict') != 'mixed' ] }
        .join(ch_patient_asm.map { meta, fasta -> [ meta.id, meta, fasta ] })
        .filter { _id, single, _meta, _fasta -> single }
        .map { _id, _single, meta, fasta -> [ meta, fasta ] }

    // ---- typing / annotation: confirmed sample-level + confirmed single-strain patient-level
    PATHOGEN_TYPING(ch_confirmed.filter { meta, _fasta -> meta.level == 'sample' }.mix(ch_single_strain), databases)

    // ---- report inputs
    // One row per unit that was extracted (or pooled), taken BEFORE the gate: a unit that is lost
    // between extraction and the coverage gate (e.g. no mapping output) still reaches the report,
    // which reports it as LOW_COVERAGE (no gate result), never as a silent omission.
    def ch_units_tsv = channel.of('unit_id	level	sample	patient	samples	species	taxid	reference')
        .concat(KRAKENTOOLS_EXTRACTKRAKENREADS.out.extracted_kraken2_reads.mix(CAT_POOLED.out.reads).map { meta, _reads ->
            [ meta.id, meta.level, meta.sample, meta.patient, meta.samples.join(';'), meta.species, meta.taxid, meta.reference ].join('	') })
        .collectFile(name: 'units.tsv', newLine: true, sort: false)

    def ch_report_files = GATE_SAMPLE.out.tsv.mix(GATE_PATIENT.out.tsv).map { _meta, f -> [ 'gate', f ] }
        .mix(ch_asm.map { _meta, f -> [ 'assemblies', f ] })
        .mix(SPECIES_CONFIRM.out.tsv.map { _meta, f -> [ 'species', f ] })
        .mix(MIXED_STRAIN_CHECK.out.tsv.map { _meta, f -> [ 'mixed', f ] })
        .mix(PATHOGEN_ASSEMBLY.out.checkm2_tsv.map { _meta, f -> [ 'checkm2', f ] })
        .mix(PATHOGEN_ASSEMBLY.out.quast_tsv.map { _meta, f -> [ 'quast', f ] })
        .mix(PATHOGEN_TYPING.out.files)

    emit:
    panel_resolved = RESOLVE_PANEL.out.tsv
    candidates     = PATHOGEN_CANDIDATES.out.candidates
    samples_tsv    = PATHOGEN_CANDIDATES.out.samples
    units_tsv      = ch_units_tsv
    report_files   = ch_report_files
    multiqc        = PATHOGEN_ASSEMBLY.out.quast_tsv.map { _meta, f -> f }.mix(PATHOGEN_ASSEMBLY.out.checkm2_tsv.map { _meta, f -> f })
}
