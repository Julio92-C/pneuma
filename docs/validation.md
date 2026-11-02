# Julio92-C/pneuma: Validation

## v0.3.0 on a clinical ONT cohort — 2 November 2026

First end-to-end run of `v0.3.0` (`--pathogens`) on real clinical data, reviewed
stage by stage against the previous nf-core/taxprofiler processing of the same
libraries. This page records what the run demonstrated and the eight issues it
surfaced, so they can be fixed or documented rather than rediscovered.

**Run:** 143 Oxford Nanopore libraries — 137 airway samples from 32 patients
across two sequencing batches and seven runs, plus 6 negative extraction
controls. Nextflow 26.04.6, pipeline `v0.3.0-g0741045`, MERIDIAN v1.3.0,
Kraken2 2.1.6 (PlusPF-8), hostile 2.0.2, nanoq 0.10.0, Flye 2.9.5-b1801,
Medaka 2.2.2, CheckM2 1.1.0, skani 0.2.2, AMRFinderPlus 4.2.7 (DB 2026-08-07.1).
Wall time 31 m 37 s per MERIDIAN config. Run settings as recorded in
`pipeline_info/params_2026-11-01_14-04-06.json`: `--pathogens --run_mobsuite`,
Kraken2 PlusPF-8, host index GRCh38.p14, `nanoq_min_length 100`,
`nanoq_min_qscore 7`, `min_reads 100`, `bracken_level S`, `bracken_threshold 10`,
ABRicate CARD/VFDB/PlasmidFinder at 80% identity and coverage, coverage gate 10x
depth over 50% breadth, Flye `--nano-hq --meta`, `run_porechop false`.

### Stage verdict

| Stage | Result |
|---|---|
| nanoq, hostile, Kraken2, Krona, ABRicate, Recentrifuge | 143/143 libraries ✅ |
| Bracken | 141/143 — see **V6** |
| MultiQC, custom tables, software manifest | complete ✅ |
| Pathogen candidates → coverage gate → Flye/Medaka | 14 assemblies produced ✅ |
| Species confirmation | 1 of 14 confirmed — see **V2**, **V3**, **V4** |
| Typing / annotation (MLST, AMRFinder, MOB-suite, Bakta, ABRicate) | ran for the 1 confirmed unit — correct gating, but see **V2** |
| MERIDIAN | 16 of 18 stages — see **V5**; statistics see **V1** |

### Confirmed working: the QC path

The read-handling rewrite is the clearest win and is corroborated from three
directions. Against taxprofiler on the same libraries:

| | taxprofiler | pneuma v0.3.0 |
|---|---|---|
| Reads surviving QC | 4.4% (median) | **98.4%** |
| Non-host reads | baseline | **2.6× more** (median; up to 126×) |
| Host removed | 58% / 91% by sample type | 94.2% / 98.4% |
| Classified fraction | 73.2% | 73.0% / 86.9% by sample type |

The classified fraction reproducing the earlier value is the cross-check: the
extra reads are real sequence, not noise. The host percentage rises because the
short reads `nanoq` now retains are predominantly human — the previous QC
discarded them *before* host depletion. Libraries that the old processing left
with 18–650 non-host reads came through this run with 530–15,300. One sequencing
run that under-delivered by ~40× at the instrument remains under-powered, as
expected: reprocessing cannot recover reads that were never generated.

### Issues

| | Severity | Finding | Suggested change |
|---|---|---|---|
| **V1** | **High** | `--meridian_random_effect_tests` defaults to `false`, so a config that sets `random_effect` still gets `Permutation: free` in PERMANOVA. On repeated-measures data (several libraries per patient) this is anticonservative: a variable whose R² sat at the cohort's no-effect level returned p = 0.006. | When `random_effect` is set and `random_effect_tests` is false, warn at config load. Consider making the patient-aware default the moment `random_effect` is present. |
| **V2** | **High** | 13 of 14 assemblies failed `SPECIES_CONFIRM`, so typing and annotation silently ran for one unit. The run exits 0 and nothing in MultiQC states how many units were rejected or why. | Add a gate-summary row to MultiQC and a `confirmed`/`reason` tally to `pathogen_cohort.tsv`; distinguish `ani_below_min` from `best_hit_other_species` in the report. |
| **V3** | **High** | Nothing gates assembly purity before species confirmation. CheckM2 contamination reached 44%, assembly size exceeded the panel reference by 40–70%, and skani aligned fraction fell to 37% — all consistent with mixed-species bins from read recruitment of close congeners. Those bins then fail ANI, which reads as a threshold problem rather than an assembly problem. Typing the failed assemblies as a diagnostic (below) corroborates this from a fourth direction. | Gate on CheckM2 contamination and on assembly size versus the panel's `genome_mb` *before* skani, and surface a MIMAG grade per unit (this run: 2 high, 2 medium, 10 low/failed). |
| **V4** | **Medium** | For species whose panel row carries `ani_discriminates=no`, ANI still produces a `best_hit_other_species` call against a sister taxon. A genuine detection was reported as its ~98%-ANI relative. | Suppress the species call for `ani_discriminates=no` rows and emit an explicit "requires marker genes / MLST" status instead. |
| **V5** | **Medium** | Two MERIDIAN stages fail deterministically in every config: `shared_fraction` (`factor level [N] is duplicated`) and `network` (circlize `gap.degree` — no space to allocate sectors). The pipeline still reports success. | Fix both upstream; until then, propagate a non-zero signal or a MultiQC warning when `stage_errors.csv` is non-empty. |
| **V6** | **Low** | Bracken emits no file for near-empty libraries (2 of 143 here), so `taxonomy/bracken/` is silently short and any downstream join loses those samples without notice. | Write an empty-but-present table, or record the skip in the manifest. |
| **V7** | **Low** | Negative controls are carried into pathogen candidate generation and the coverage gate, producing per-control coverage files for panel species. | Exclude `is_control` samples from the pathogen subworkflow, or make it an opt-in flag. |
| **V8** | **Info** | `database_manifest.tsv` records a SHA-256 for the CheckM2, Bakta, AMRFinderPlus, MOB-suite, Medaka, MLST and panel-reference databases, but not for the Kraken2 database or the host index — the two that determine every taxonomic result. | Record both in the manifest (path plus checksum), since they are the provenance a reader of the taxonomy needs. |

### The species-confirmation gate was worth having

Typing was run separately over the 13 assemblies the gate rejected, to test whether holding
their results back was over-cautious. It was not, and the exercise is worth recording because
it is the strongest available evidence that `SPECIES_CONFIRM` is doing real work:

- **No assembly yielded a sequence type.** All 13 returned a blank ST; two returned no scheme
  at all, and those two were the highest-value targets in the panel.
- **`mlst` selected a relative's scheme** where the species has none of its own: a
  *Haemophilus* commensal was typed against the *H. influenzae* scheme. (The *Neisseria*
  commensals were typed with the PubMLST *Neisseria* scheme, which covers the whole genus.)
  Most blank STs are new alleles rather than missing loci: commensal *Neisseria* and mitis-group
  streptococci are highly diverse and poorly represented in PubMLST [7, 9], and a single
  residual indel in a Nanopore consensus also turns an exact allele into a near match [4].
- **Resistance genes may not belong to the species named.** A *Neisseria* bin carried *bla*TEM
  with two aminoglycoside phosphotransferases on one contig, and another carried *msr*(D) and
  *tet*(M). Both are plausible in commensal *Neisseria*, where *tet*(M) (often on the conjugative
  plasmid pConj) and macrolide resistance are common [7, 8, 9]; but in pooled or contaminated
  bins they may equally come from another organism. Reported under the panel species' name
  without that caveat, a consumer of `pathogen_cohort.tsv` could have acted on them. The read-back
  check and the gene table (contig, chromosome or plasmid) now give the evidence to decide.

Two consequences for the pipeline, both now implemented. First, the panel's `mlst_scheme`
column pins the scheme per species and `none` marks species without one. Second, a failed
species confirmation is surfaced instead of silently withholding results: every assembly is
typed, but an unconfirmed one keeps its status, is reported under its best match (`typed_as`)
and carries its purity grade, so the evidence is visible and labelled rather than hidden.
Typing directly from reads (Krocus-style, per-locus consensus) recovers STs that Nanopore
assemblies miss [4, 5] and would also show minor alleles in pooled samples; it is planned as a
prototype.

Run parameters are recorded: `pipeline_info/` carries the params JSON, execution
report, timeline, trace and DAG for each attempt, and `database_manifest.tsv` the
database versions. The findings above were checked against
`params_2026-11-01_14-04-06.json` rather than inferred from the outputs.

### Follow-up (development branch, 2 November 2026)

| | Status |
|---|---|
| **V1** | Fixed: `--meridian_random_effect_tests` defaults to `true` and applies whenever `--meridian_random_effect` is set; turning it off logs a warning. |
| **V2** | Fixed: `confirm_reason` per unit, a "Pathogen reconstruction funnel" MultiQC section, and every assembly typed (unconfirmed ones under their best match, `typed_as`). |
| **V3** | Fixed as flags rather than a gate (author's decision): `mimag_grade` (CheckM2) and `size_ratio` (assembly / panel reference) per unit; the read-back check (own reads mapped to the assembly) runs for every assembly. |
| **V4** | Fixed with the panel's existing `confusable_with` list (there is no `ani_discriminates` column): a confusable best hit within `--pathogen_ani_confusable_delta` (0.5) ANI points of the panel species gives `NOT_RESOLVED_BY_ANI`. |
| **V5** | Fixed upstream in MERIDIAN v1.3.1 (both failures reproduced and fixed on this cohort); failed stages also appear in MultiQC. |
| **V6** | Fixed: `taxonomy/bracken/bracken_skipped.tsv` lists each skipped sample and why. |
| **V7** | Fixed: negative controls are detection only unless `--pathogen_gate_controls`. |
| **V8** | Fixed: databases passed by path at run time (Kraken2, taxdump, host index) get a row with path and sha256 in `pipeline_info/database_manifest.tsv`, alongside the `--db_dir` manifest. |

### Not covered

`--build_databases`, resumption behaviour, and the detect-only panel tiers. One
unexplained observation is left open rather than attributed to the pipeline: two
high-input negative controls returned large non-host read counts that were
87–93% unclassified, and the two pipelines disagree about their composition.

### References

1. Charalampous T, et al. Nanopore metagenomics enables rapid clinical diagnosis of bacterial lower respiratory infection. *Nat Biotechnol.* 2019.
2. Behruznia M, et al. Metagenomic sequencing of sputum enables accurate, non-invasive diagnosis of community-acquired pneumonia. *JAC-Antimicrob Resist.* 2025.
3. Wagner GE, et al. Real-time nanopore Q20+ sequencing enables extremely fast and accurate core genome MLST typing. *J Clin Microbiol.* 2023.
4. Bejaoui S, et al. Comparison of Illumina and Oxford Nanopore sequencing data quality for *Clostridioides difficile* genome analysis. *BMC Genomics.* 2025.
5. García-Pérez I, et al. NanoMLST: a high-throughput MLST workflow using Oxford Nanopore sequencing for ESKAPE+E pathogens. *MicrobiologyOpen.* 2025.
6. Mac Aogáin M, et al. Metagenomics reveals a core macrolide resistome related to microbiota in chronic respiratory disease. *Am J Respir Crit Care Med.* 2020.
7. Regan MR, et al. Surveillance of antibiotic resistance within the commensal *Neisseria*. *Microbiol Spectr.* 2026.
8. Robinson LR, et al. Doxycycline post-exposure prophylaxis and the selection of doxycycline resistance in *Neisseria* commensals. *Sci Rep.* 2025.
9. Miari VF, et al. Carriage and antimicrobial susceptibility of commensal *Neisseria* species from the human oropharynx. *Sci Rep.* 2024.
10. Dmitrijeva M, et al. Strain-resolved dynamics of the lung microbiome in patients with cystic fibrosis. *mBio.* 2020.
11. Fan Y, et al. Long-read metagenomics for strain tracking after fecal microbiota transplant. *Nat Microbiol.* 2025.
12. Podlesny D, et al. Metagenomic strain detection with SameStr. *Microbiome.* 2022.
