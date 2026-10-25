"""Synthetic, contract-shaped inputs for bin/pathogen_report.py: one case per status.

Used by tests/python/test_pathogen_report.py, and materialised for the BUILD_PATHOGEN_REPORT nf-test:
    python3 tests/python/report_fixture.py tests/data/pathogen_report
Layouts follow the real tool outputs captured in tests/data/pathogen_contracts/. IDs are invented.
"""
import sys
from pathlib import Path

PANEL_COLUMNS = ["species", "slug", "tier", "who_category", "action", "taxid", "scientific_name", "matched_name_class",
                 "taxids", "reference_accession", "confusable_with", "confusable_accessions", "ani_min", "amrfinder_organism"]
PANEL = [
    ["Streptococcus pneumoniae", "Streptococcus_pneumoniae", "1", "BPPL 2024 medium", "reconstruct", "1313",
     "Streptococcus pneumoniae", "scientific name", "1313;170187", "GCF_001457635.1", "Streptococcus mitis",
     "GCF_018603625.2", "95.0", "Streptococcus_pneumoniae"],
    ["Haemophilus influenzae", "Haemophilus_influenzae", "1", "BPPL 2024 medium", "reconstruct", "727",
     "Haemophilus influenzae", "scientific name", "727", "GCF_020736045.1", "", "", "95.0", "Haemophilus_influenzae"],
    ["Mycobacterium tuberculosis", "Mycobacterium_tuberculosis", "D", "BPPL 2024 critical", "detect", "1773",
     "Mycobacterium tuberculosis", "scientific name", "1773", "", "", "", "", ""],
    ["Streptococcus mitis", "Streptococcus_mitis", "R", "", "reconstruct", "28037", "Streptococcus mitis",
     "scientific name", "28037", "GCF_018603625.2", "Streptococcus pneumoniae", "GCF_001457635.1", "95.0", ""],
]
SPN, HIN, MTB, SMI = (row[0] for row in PANEL)
SAMPLES = [("S1", "P1", "false", "false", 20000), ("S2", "P2", "false", "false", 15000),
           ("S3", "P3", "false", "false", 18000), ("S4", "P4", "false", "false", 9000),
           ("S5", "P4", "false", "false", 8000), ("S6", "P5", "false", "false", 9000),
           ("S7", "P5", "false", "false", 9500), ("S8", "P6", "false", "true", 50),
           ("CTRL", "", "true", "false", 800)]
CANDIDATES = {  # (sample, species): (clade_reads, decision); default (0, not_detected | detect | low_depth)
    ("S1", SPN): (3000, "extract"), ("S1", MTB): (5, "detect"),
    ("S2", SPN): (6, "below_min_reads"), ("S2", SMI): (900, "extract"),
    ("S3", SPN): (1500, "extract"), ("S3", HIN): (1200, "extract"), ("S3", SMI): (1100, "extract"),
    ("S4", HIN): (400, "extract"), ("S5", HIN): (380, "extract"),
    ("S6", SPN): (350, "extract"), ("S7", SPN): (360, "extract"),
    ("S8", SPN): (12, "low_depth"),
}
UNITS = [  # unit, level, sample, patient, samples, species, gate pass, depth, breadth
    ("S1.Streptococcus_pneumoniae", "sample", "S1", "P1", "S1", SPN, "true", "25.1", "0.98"),
    ("S2.Streptococcus_mitis", "sample", "S2", "P2", "S2", SMI, "true", "14.0", "0.91"),
    ("S3.Streptococcus_pneumoniae", "sample", "S3", "P3", "S3", SPN, "true", "12.2", "0.62"),
    ("S3.Haemophilus_influenzae", "sample", "S3", "P3", "S3", HIN, "true", "11.0", "0.95"),
    ("S3.Streptococcus_mitis", "sample", "S3", "P3", "S3", SMI, "true", "15.3", "0.93"),
    ("S4.Haemophilus_influenzae", "sample", "S4", "P4", "S4", HIN, "false", "6.1", "0.88"),
    ("S5.Haemophilus_influenzae", "sample", "S5", "P4", "S5", HIN, "false", "5.9", "0.87"),
    ("patient-P4.Haemophilus_influenzae", "patient", "", "P4", "S4;S5", HIN, "true", "12.0", "0.97"),
    ("S6.Streptococcus_pneumoniae", "sample", "S6", "P5", "S6", SPN, "false", "6.0", "0.90"),
    ("S7.Streptococcus_pneumoniae", "sample", "S7", "P5", "S7", SPN, "false", "6.2", "0.90"),
    ("patient-P5.Streptococcus_pneumoniae", "patient", "", "P5", "S6;S7", SPN, "true", "12.2", "0.97"),
]
ASSEMBLED = ["S1.Streptococcus_pneumoniae", "S3.Streptococcus_pneumoniae", "S3.Haemophilus_influenzae",
             "S3.Streptococcus_mitis", "patient-P4.Haemophilus_influenzae", "patient-P5.Streptococcus_pneumoniae"]
SPECIES = {  # unit: (best species, accession, ani, confirmed, reason)
    "S1.Streptococcus_pneumoniae": (SPN, "GCF_001457635.1", "98.9", "true", "confirmed"),
    "S3.Streptococcus_pneumoniae": (SMI, "GCF_018603625.2", "96.8", "false", "best_hit_other_species"),
    "S3.Haemophilus_influenzae": (HIN, "GCF_020736045.1", "97.5", "true", "confirmed"),
    "S3.Streptococcus_mitis": (SMI, "GCF_018603625.2", "95.6", "true", "confirmed"),
    "patient-P4.Haemophilus_influenzae": (HIN, "GCF_020736045.1", "98.1", "true", "confirmed"),
    "patient-P5.Streptococcus_pneumoniae": (SPN, "GCF_001457635.1", "98.7", "true", "confirmed"),
}
MLST = {
    "S1.Streptococcus_pneumoniae": "spneumoniae\t1\taroE(1)\tgdh(1)\tgki(1)\trecP(1)\tspi(1)\txpt(1)\tddl(1)",
    "S3.Haemophilus_influenzae": "hinfluenzae\t-\tadk(1)\tatpG(~3)\tfrdB(-)\tfucK(1)\tmdh(1)\tpgi(1)\trecA(1)",
    "S3.Streptococcus_mitis": "-\t-",
    "patient-P4.Haemophilus_influenzae": "hinfluenzae\t3\tadk(1)\tatpG(1)\tfrdB(1)\tfucK(1)\tmdh(1)\tpgi(1)\trecA(1)",
}
MIXED = {"patient-P4.Haemophilus_influenzae": ("40000", "2", "0.050", "single"),
         "patient-P5.Streptococcus_pneumoniae": ("38000", "190", "5.000", "mixed")}
AMR_HEADER = ("Protein id\tContig id\tStart\tStop\tStrand\tElement symbol\tElement name\tScope\tType\tSubtype\tClass\t"
              "Subclass\tMethod\tTarget length\tReference sequence length\t% Coverage of reference\t"
              "% Identity to reference\tAlignment length\tClosest reference accession\tClosest reference name\t"
              "HMM accession\tHMM description")
ABRICATE_TAB = ("#FILE\tSEQUENCE\tSTART\tEND\tSTRAND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS\t%COVERAGE\t%IDENTITY\t"
                "DATABASE\tACCESSION\tPRODUCT\tRESISTANCE")
ABRICATE_CSV = ("#FILE,SEQUENCE,START,END,STRAND,GENE,COVERAGE,COVERAGE_MAP,GAPS,%COVERAGE,%IDENTITY,DATABASE,"
                "ACCESSION,PRODUCT,RESISTANCE")


def put(path, text):
    """Write with LF endings on every OS so the committed nf-test fixture is byte-identical everywhere."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="\n") as fh:
        fh.write(text)


def tsv(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    put(path, "\t".join(header) + "\n" + "".join("\t".join(str(c) for c in r) + "\n" for r in rows))


def write(root):
    root = Path(root)
    tsv(root / "panel_resolved.tsv", PANEL_COLUMNS, PANEL)
    tsv(root / "pathogen_samples.tsv",
        ["sample", "patient", "is_control", "low_depth", "n_reads", "classified_reads", "unclassified_reads", "unclassified_pct"],
        [(s, p, c, ld, n, n - 100, 100, "%.2f" % (10000.0 / n)) for s, p, c, ld, n in SAMPLES])
    cand = []
    for s, p, c, ld, _n in SAMPLES:
        for row in PANEL:
            sp, slug, tier, action, taxid = row[0], row[1], row[2], row[4], row[5]
            default = "low_depth" if ld == "true" else ("detect" if action == "detect" else "not_detected")
            clade, dec = CANDIDATES.get((s, sp), (0, default))
            cand.append((s, p, c, ld, sp, slug, tier, action, taxid, clade, clade, "0.01" if clade else "0", dec))
    tsv(root / "pathogen_candidates.tsv",
        ["sample", "patient", "is_control", "low_depth", "species", "slug", "tier", "action", "taxid", "clade_reads",
         "bracken_reads", "bracken_fraction", "decision"], cand)
    taxid = {row[0]: row[5] for row in PANEL}
    ref = {row[0]: row[9] for row in PANEL}
    tsv(root / "units.tsv", ["unit_id", "level", "sample", "patient", "samples", "species", "taxid", "reference"],
        [(u, lv, s, p, ss, sp, taxid[sp], ref[sp]) for u, lv, s, p, ss, sp, *_g in UNITS])
    for u, _lv, _s, _p, _ss, _sp, ok, depth, breadth in UNITS:
        tsv(root / "gate" / ("%s.gate.tsv" % u),
            ["unit_id", "reads_mapped", "mean_depth", "breadth", "gate_min_depth", "gate_min_breadth", "pass"],
            [(u, 1000, depth, breadth, 10.0, 0.5, ok)])
    for u in ASSEMBLED:
        (root / "assemblies").mkdir(parents=True, exist_ok=True)
        put((root / "assemblies" / ("%s.fasta" % u)), ">contig_1\n" + "ACGT" * 50 + "\n")
    for u, (best, acc, ani, ok, reason) in SPECIES.items():
        panel_sp = [x[5] for x in UNITS if x[0] == u][0]
        tsv(root / "species" / ("%s.species.tsv" % u),
            ["unit_id", "panel_species", "best_species", "best_accession", "ani", "af_query", "confirmed", "reason"],
            [(u, panel_sp, best, acc, ani, "88.0", ok, reason)])
    for u, line in MLST.items():
        (root / "mlst").mkdir(parents=True, exist_ok=True)
        put((root / "mlst" / ("%s.tsv" % u)), "%s.fasta\t%s\n" % (u, line))
    tsv(root / "amrfinderplus" / "S1.Streptococcus_pneumoniae.tsv", AMR_HEADER.split("\t"),
        [("NA", "contig_1", 1000, 2920, "+", "tet(M)", "Tet(M)", "core", "AMR", "AMR", "TETRACYCLINE", "TETRACYCLINE",
          "EXACTX", 639, 639, "100.00", "100.00", 639, "WP_1", "Tet(M)", "NA", "NA"),
         ("NA", "contig_1", 5000, 6000, "+", "ply", "pneumolysin", "plus", "VIRULENCE", "VIRULENCE", "", "",
          "EXACTX", 471, 471, "100.00", "100.00", 471, "WP_2", "pneumolysin", "NA", "NA")])
    tsv(root / "abricate" / "S1.Streptococcus_pneumoniae.vfdb.txt", ABRICATE_TAB.split("\t"),
        [("S1.Streptococcus_pneumoniae.vfdb.fasta", "contig_1", 1, 1416, "+", "ply", "1-1416/1416", "===", "0/0",
          "100.00", "99.93", "vfdb", "NP_1", "pneumolysin", "")])
    tsv(root / "abricate" / "S1.Streptococcus_pneumoniae.plasmidfinder.txt", ABRICATE_TAB.split("\t"), [])
    for u, (ev, het, per_kb, verdict) in MIXED.items():
        tsv(root / "mixed" / ("%s.mixed.tsv" % u),
            ["unit_id", "evaluable_sites", "het_sites", "het_per_kb", "verdict", "min_depth", "min_maf", "min_alt", "max_het_per_kb"],
            [(u, ev, het, per_kb, verdict, 8, 0.25, 3, 1.0)])
    tsv(root / "checkm2" / "pathogens_checkm2_report.tsv",
        ["Name", "Completeness", "Contamination", "Completeness_Model_Used", "Translation_Table_Used"],
        [(u, "97.5", "1.2", "Gradient Boost (General Model)", 11) for u in ASSEMBLED])
    tsv(root / "quast" / "pathogens_quast.tsv", ["Assembly"] + ASSEMBLED,
        [["# contigs"] + ["3"] * len(ASSEMBLED), ["Total length"] + ["2100000"] * len(ASSEMBLED),
         ["N50"] + ["900000"] * len(ASSEMBLED)])
    (root / "functional").mkdir(parents=True, exist_ok=True)
    put((root / "functional" / "S2_abricate.csv"),
        ABRICATE_CSV + "\n"
        "S2_kraken2ID_PFP.fasta,read7_kraken:taxid|1313,1,1900,+,tet(M),1-1920/1920,===,0/0,99.00,99.10,card,AM1,Tet(M),tetracycline\n"
        "S2_kraken2ID_PFP.fasta,read9_kraken:taxid|562,1,800,+,blaTEM-1,1-861/861,===,0/0,93.00,99.00,card,AM2,TEM-1,penam\n")
    put((root / "functional" / "S1_abricate.csv"), ABRICATE_CSV + "\n")
    put((root / "pathogen_panel.csv"), "species,tier\n" + "".join("%s,%s\n" % (r[0], r[2]) for r in PANEL))
    tsv(root / "database_manifest.tsv", ["name", "version", "path", "sha256"],
        [("checkm2_db", "3", "/db/checkm2/checkm2_db_v3.dmnd", "a" * 64)])
    return root


if __name__ == "__main__":
    write(sys.argv[1] if len(sys.argv) > 1 else "tests/data/pathogen_report")
