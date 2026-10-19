#!/usr/bin/env python3
"""Generate the SYNTHETIC MERIDIAN integration fixture (stdlib only, deterministic).

Writes tests/data/meridian_fixture/ with the files PNEUMA hands to MERIDIAN, in the
exact layouts PNEUMA produces (no real or patient data; IDs are anonymised):

  kraken2_combined_reports.txt   KrakenTools combine_kreports, preamble stripped by ARRANGE_REPORTS
                                 (starts at '#perc', numbered <i>_all/<i>_lvl columns)
  bracken_arranged.csv           #perc,tot_all,<samples>,taxid,name
  pneuma.rcf.data.csv            raw Recentrifuge 2.1.2 CSV after the RECENTRIFUGE module's
                                 sed (row1 'Samples', row2 'Stats'); derived EXCLUSIVE/SHARED
                                 triplets and Details/Details/Name/Rank columns included
  summary_reportClean.csv        merged ABRicate (CARD / VFDB / PlasmidFinder) hits
  metadata.csv                   sample,group,ID
  raw/kraken_combined_raw.txt    combine_kreports output with preamble (input of ARRANGE_REPORTS)
  raw/bracken_combined_raw.txt   same for Bracken (totals-only unclassified row)

Samples: S01-S04 (group A), S05-S08 (group B), CTRL1 (negative control). S07 has zero
ABRicate hits. Columns are ordered CTRL1, S01..S08 (PNEUMA sorts by sample id).

Usage: python tests/make_meridian_fixture.py [outdir]
"""
import math
import random
import sys
import uuid
from pathlib import Path

SEED = 20261019
SAMPLES = ["CTRL1"] + [f"S{i:02d}" for i in range(1, 9)]
GROUP = {f"S{i:02d}": ("A" if i <= 4 else "B") for i in range(1, 9)}
ZERO_HIT_SAMPLE = "S07"

CELL = [("R1", 131567, "cellular organisms")]
BAC = CELL + [("D", 2, "Bacteria")]
EUK = CELL + [("D", 2759, "Eukaryota")]
FIRM = BAC + [("P", 1239, "Bacillota")]
ACTI = BAC + [("P", 201174, "Actinomycetota")]
BACT = BAC + [("P", 976, "Bacteroidota")]
FUSO = BAC + [("P", 32066, "Fusobacteriota")]
PSEU = BAC + [("P", 1224, "Pseudomonadota")]
STREP = FIRM + [("C", 91061, "Bacilli"), ("O", 186826, "Lactobacillales"),
                ("F", 1300, "Streptococcaceae"), ("G", 1301, "Streptococcus")]
STAPH = FIRM + [("C", 91061, "Bacilli"), ("O", 1385, "Bacillales"),
                ("F", 90964, "Staphylococcaceae"), ("G", 1279, "Staphylococcus")]
VEIL = FIRM + [("C", 909932, "Negativicutes"), ("O", 909929, "Veillonellales"),
               ("F", 31977, "Veillonellaceae"), ("G", 29465, "Veillonella")]
ROTH = ACTI + [("C", 1760, "Actinomycetes"), ("O", 85006, "Micrococcales"),
               ("F", 1268, "Micrococcaceae"), ("G", 32207, "Rothia")]
CORY = ACTI + [("C", 1760, "Actinomycetes"), ("O", 85007, "Corynebacteriales"),
               ("F", 1653, "Corynebacteriaceae"), ("G", 1716, "Corynebacterium")]
CUTI = ACTI + [("C", 1760, "Actinomycetes"), ("O", 85009, "Propionibacteriales"),
               ("F", 31957, "Propionibacteriaceae"), ("G", 1912216, "Cutibacterium")]
PREV = BACT + [("C", 200643, "Bacteroidia"), ("O", 171549, "Bacteroidales"),
               ("F", 171552, "Prevotellaceae"), ("G", 838, "Prevotella")]
FUSOB = FUSO + [("C", 203490, "Fusobacteriia"), ("O", 203491, "Fusobacteriales"),
                ("F", 203492, "Fusobacteriaceae"), ("G", 848, "Fusobacterium")]
GAMMA = PSEU + [("C", 1236, "Gammaproteobacteria")]
HAEM = GAMMA + [("O", 135625, "Pasteurellales"), ("F", 712, "Pasteurellaceae"), ("G", 724, "Haemophilus")]
PSEUD = GAMMA + [("O", 72274, "Pseudomonadales"), ("F", 135621, "Pseudomonadaceae"), ("G", 286, "Pseudomonas")]
MORAX = GAMMA + [("O", 72274, "Pseudomonadales"), ("F", 468, "Moraxellaceae"), ("G", 475, "Moraxella")]
ACIN = GAMMA + [("O", 72274, "Pseudomonadales"), ("F", 468, "Moraxellaceae"), ("G", 469, "Acinetobacter")]
KLEB = GAMMA + [("O", 91347, "Enterobacterales"), ("F", 543, "Enterobacteriaceae"), ("G", 570, "Klebsiella")]
NEIS = PSEU + [("C", 28216, "Betaproteobacteria"), ("O", 206351, "Neisseriales"),
               ("F", 481, "Neisseriaceae"), ("G", 482, "Neisseria")]
RALS = PSEU + [("C", 28216, "Betaproteobacteria"), ("O", 80840, "Burkholderiales"),
               ("F", 119060, "Burkholderiaceae"), ("G", 48736, "Ralstonia")]
BRAD = PSEU + [("C", 28211, "Alphaproteobacteria"), ("O", 356, "Hyphomicrobiales"),
               ("F", 41294, "Nitrobacteraceae"), ("G", 374, "Bradyrhizobium")]
HOMO = EUK + [("P", 7711, "Chordata"), ("C", 40674, "Mammalia"), ("O", 9443, "Primates"),
              ("F", 9604, "Hominidae"), ("G", 9605, "Homo")]

# (lineage, taxid, name, mean_A, mean_B, mean_CTRL)
SPECIES = [
    (STREP, 28037, "Streptococcus mitis", 6000, 3000, 0),
    (STREP, 1303, "Streptococcus oralis", 3500, 1800, 0),
    (STREP, 1304, "Streptococcus salivarius", 2500, 1500, 0),
    (STREP, 1318, "Streptococcus parasanguinis", 1500, 900, 0),
    (STREP, 1314, "Streptococcus pyogenes", 40, 200, 0),
    (STREP, 1313, "Streptococcus pneumoniae", 600, 2500, 0),
    (STAPH, 1280, "Staphylococcus aureus", 150, 4000, 0),
    (STAPH, 1282, "Staphylococcus epidermidis", 400, 500, 60),
    (VEIL, 29466, "Veillonella parvula", 5000, 1800, 0),
    (VEIL, 39777, "Veillonella atypica", 2500, 900, 0),
    (VEIL, 39778, "Veillonella dispar", 3000, 1200, 0),
    (ROTH, 43675, "Rothia mucilaginosa", 4000, 1500, 0),
    (ROTH, 2047, "Rothia dentocariosa", 800, 300, 0),
    (CORY, 43770, "Corynebacterium striatum", 100, 700, 0),
    (CUTI, 1747, "Cutibacterium acnes", 300, 300, 900),
    (PREV, 28132, "Prevotella melaninogenica", 7000, 2500, 0),
    (PREV, 28131, "Prevotella intermedia", 1500, 600, 0),
    (PREV, 1177574, "Prevotella jejuni", 2000, 700, 0),
    (PREV, 28129, "Prevotella denticola", 900, 300, 0),
    (FUSOB, 851, "Fusobacterium nucleatum", 1200, 1500, 0),
    (HAEM, 727, "Haemophilus influenzae", 500, 4500, 0),
    (HAEM, 729, "Haemophilus parainfluenzae", 1500, 900, 0),
    (PSEUD, 287, "Pseudomonas aeruginosa", 60, 6000, 0),
    (MORAX, 480, "Moraxella catarrhalis", 80, 900, 0),
    (ACIN, 470, "Acinetobacter baumannii", 20, 500, 0),
    (KLEB, 573, "Klebsiella pneumoniae", 30, 2500, 0),
    (NEIS, 28449, "Neisseria subflava", 3000, 800, 0),
    (NEIS, 488, "Neisseria mucosa", 1200, 400, 0),
    (NEIS, 487, "Neisseria meningitidis", 100, 60, 0),
    (RALS, 329, "Ralstonia pickettii", 80, 90, 2500),   # kit taxa: control-only baseline
    (RALS, 190721, "Ralstonia insidiosa", 60, 70, 1500),
    (BRAD, 375, "Bradyrhizobium japonicum", 40, 50, 1200),
    (HOMO, 9606, "Homo sapiens", 60000, 45000, 300),
]
SPECIES_PROFILE = {s[1]: s for s in SPECIES}

# Fraction of a node's clade reads assigned directly to it (not to a child), by rank code.
INTERNAL_LVL_FRAC = {"G": 0.05, "F": 0.012, "O": 0.005, "C": 0.003, "P": 0.0015,
                     "D": 0.0006, "R1": 0.0003}

# ABRicate gene catalogue per species: (database, gene, accession, product, resistance)
CARD, VFDB, PF = "card", "vfdb", "plasmidfinder"
GENES = {
    1313: [(CARD, "tet(M)", "NG_048200.1", "tetracycline-resistant ribosomal protection protein Tet(M)", "tetracycline"),
           (CARD, "erm(B)", "NG_047831.1", "23S rRNA methyltransferase Erm(B)", "macrolide;lincosamide;streptogramin"),
           (CARD, "mef(A)", "NG_048003.1", "macrolide efflux protein Mef(A)", "macrolide"),
           (CARD, "msr(D)", "NG_048004.1", "ABC-F ribosomal protection protein Msr(D)", "macrolide"),
           (VFDB, "ply", "VFG000100", "pneumolysin", ""),
           (VFDB, "nanA", "VFG000101", "neuraminidase A", ""),
           (VFDB, "pspA", "VFG000102", "pneumococcal surface protein A", "")],
    28037: [(CARD, "tet(M)", "NG_048200.1", "tetracycline-resistant ribosomal protection protein Tet(M)", "tetracycline"),
            (CARD, "mef(A)", "NG_048003.1", "macrolide efflux protein Mef(A)", "macrolide"),
            (VFDB, "nanA", "VFG000101", "neuraminidase A", "")],
    1303: [(CARD, "tet(M)", "NG_048200.1", "tetracycline-resistant ribosomal protection protein Tet(M)", "tetracycline"),
           (CARD, "erm(B)", "NG_047831.1", "23S rRNA methyltransferase Erm(B)", "macrolide;lincosamide;streptogramin")],
    1280: [(CARD, "mecA", "NG_047945.1", "penicillin-binding protein PBP2a", "penam;cephalosporin"),
           (CARD, "blaZ", "NG_047886.1", "class A beta-lactamase BlaZ", "penam"),
           (CARD, "erm(C)", "NG_047832.1", "23S rRNA methyltransferase Erm(C)", "macrolide;lincosamide;streptogramin"),
           (CARD, "norA", "NG_048001.1", "multidrug efflux pump NorA", "fluoroquinolone"),
           (VFDB, "hla", "VFG001200", "alpha-hemolysin", ""),
           (VFDB, "hlgB", "VFG001201", "gamma-hemolysin component B", ""),
           (VFDB, "sak", "VFG001202", "staphylokinase", ""),
           (VFDB, "scn", "VFG001203", "staphylococcal complement inhibitor", ""),
           (PF, "rep7a_1_repC(Cassette)", "CP000000", "rep7a_1_repC(Cassette)", "")],
    1282: [(CARD, "mecA", "NG_047945.1", "penicillin-binding protein PBP2a", "penam;cephalosporin"),
           (VFDB, "icaA", "VFG001210", "intercellular adhesion protein A", "")],
    29466: [(CARD, "tet(Q)", "NG_048210.1", "tetracycline-resistant ribosomal protection protein Tet(Q)", "tetracycline"),
            (CARD, "tet(M)", "NG_048200.1", "tetracycline-resistant ribosomal protection protein Tet(M)", "tetracycline")],
    39777: [(CARD, "tet(Q)", "NG_048210.1", "tetracycline-resistant ribosomal protection protein Tet(Q)", "tetracycline")],
    43675: [(CARD, "tet(M)", "NG_048200.1", "tetracycline-resistant ribosomal protection protein Tet(M)", "tetracycline"),
            (CARD, "erm(X)", "NG_047840.1", "23S rRNA methyltransferase Erm(X)", "macrolide;lincosamide;streptogramin")],
    43770: [(CARD, "erm(X)", "NG_047840.1", "23S rRNA methyltransferase Erm(X)", "macrolide;lincosamide;streptogramin"),
            (CARD, "aph(3')-Ia", "NG_047400.1", "aminoglycoside phosphotransferase APH(3')-Ia", "aminoglycoside")],
    28132: [(CARD, "CfxA6", "NG_047922.1", "class A beta-lactamase CfxA6", "cephalosporin;penam"),
            (CARD, "tet(Q)", "NG_048210.1", "tetracycline-resistant ribosomal protection protein Tet(Q)", "tetracycline"),
            (CARD, "erm(F)", "NG_047835.1", "23S rRNA methyltransferase Erm(F)", "macrolide;lincosamide;streptogramin")],
    28131: [(CARD, "tet(Q)", "NG_048210.1", "tetracycline-resistant ribosomal protection protein Tet(Q)", "tetracycline")],
    1177574: [(CARD, "CfxA6", "NG_047922.1", "class A beta-lactamase CfxA6", "cephalosporin;penam")],
    851: [(VFDB, "fadA", "VFG002000", "adhesin FadA", ""),
          (VFDB, "fap2", "VFG002001", "autotransporter adhesin Fap2", "")],
    727: [(CARD, "TEM-1", "NG_050145.1", "class A beta-lactamase TEM-1", "penam;cephalosporin"),
          (VFDB, "iga", "VFG002100", "IgA1 protease", ""),
          (VFDB, "hxuA", "VFG002101", "heme-hemopexin utilization protein A", ""),
          (PF, "Col(BS512)_1", "CP000001", "Col(BS512)_1", "")],
    729: [(CARD, "TEM-1", "NG_050145.1", "class A beta-lactamase TEM-1", "penam;cephalosporin")],
    287: [(CARD, "OXA-50", "NG_049900.1", "class D beta-lactamase OXA-50", "penam;carbapenem"),
          (CARD, "PDC-3", "NG_049901.1", "Pseudomonas-derived cephalosporinase PDC-3", "cephalosporin"),
          (CARD, "mexA", "NG_048020.1", "MexAB-OprM efflux membrane fusion protein MexA", "fluoroquinolone;tetracycline;penam"),
          (CARD, "mexB", "NG_048021.1", "MexAB-OprM efflux RND transporter MexB", "fluoroquinolone;tetracycline;penam"),
          (CARD, "oprM", "NG_048022.1", "MexAB-OprM outer membrane channel OprM", "fluoroquinolone;tetracycline;penam"),
          (CARD, "aph(3')-IIb", "NG_047410.1", "aminoglycoside phosphotransferase APH(3')-IIb", "aminoglycoside"),
          (CARD, "catB7", "NG_047890.1", "chloramphenicol acetyltransferase CatB7", "phenicol"),
          (VFDB, "lasB", "VFG003000", "elastase LasB", ""),
          (VFDB, "plcH", "VFG003001", "hemolytic phospholipase C", ""),
          (VFDB, "toxA", "VFG003002", "exotoxin A", ""),
          (VFDB, "exoS", "VFG003003", "exoenzyme S", "")],
    480: [(CARD, "BRO-1", "NG_049100.1", "class A beta-lactamase BRO-1", "penam"),
          (VFDB, "uspA1", "VFG004000", "ubiquitous surface protein A1", "")],
    470: [(CARD, "OXA-66", "NG_049910.1", "class D beta-lactamase OXA-66", "carbapenem;penam"),
          (CARD, "adeB", "NG_048030.1", "AdeABC efflux RND transporter AdeB", "aminoglycoside;tetracycline")],
    573: [(CARD, "SHV-1", "NG_050190.1", "class A beta-lactamase SHV-1", "penam"),
          (CARD, "oqxB", "NG_048040.1", "OqxAB efflux RND transporter OqxB", "fluoroquinolone"),
          (CARD, "fosA", "NG_047900.1", "fosfomycin thiol transferase FosA", "fosfomycin"),
          (VFDB, "fimH", "VFG005000", "type 1 fimbrin adhesin FimH", ""),
          (VFDB, "mrkD", "VFG005001", "type 3 fimbrial adhesin MrkD", ""),
          (VFDB, "ybtS", "VFG005002", "yersiniabactin biosynthesis salicylate synthase", ""),
          (PF, "IncFIB(K)_1_Kpn3", "FN000001", "IncFIB(K)_1_Kpn3", "")],
    28449: [(VFDB, "pilT", "VFG006000", "twitching motility protein PilT", "")],
    488: [(VFDB, "pilT", "VFG006000", "twitching motility protein PilT", "")],
    487: [(VFDB, "pilE", "VFG006001", "type IV pilin PilE", "")],
    329: [(CARD, "OXA-22", "NG_049950.1", "class D beta-lactamase OXA-22", "penam"),
          (CARD, "mexA", "NG_048020.1", "MexAB-OprM efflux membrane fusion protein MexA", "fluoroquinolone;tetracycline;penam")],
    190721: [(CARD, "OXA-22", "NG_049950.1", "class D beta-lactamase OXA-22", "penam")],
    1747: [(CARD, "erm(X)", "NG_047840.1", "23S rRNA methyltransferase Erm(X)", "macrolide;lincosamide;streptogramin")],
}
# VFDB function classes (ABRicate PRODUCT carries "- <Function> (VFCxxxx)")
VF_FUNC = {
    "ply": ("Exotoxin", "VFC0235"), "nanA": ("Immune modulation", "VFC0258"),
    "pspA": ("Immune modulation", "VFC0258"), "hla": ("Exotoxin", "VFC0235"),
    "hlgB": ("Exotoxin", "VFC0235"), "sak": ("Immune modulation", "VFC0258"),
    "scn": ("Immune modulation", "VFC0258"), "icaA": ("Biofilm", "VFC0271"),
    "fadA": ("Adherence", "VFC0001"), "fap2": ("Adherence", "VFC0001"),
    "iga": ("Immune modulation", "VFC0258"), "hxuA": ("Nutritional/Metabolic factor", "VFC0272"),
    "lasB": ("Exoenzyme", "VFC0301"), "plcH": ("Exoenzyme", "VFC0301"),
    "toxA": ("Exotoxin", "VFC0235"), "exoS": ("Effector delivery system", "VFC0086"),
    "uspA1": ("Adherence", "VFC0001"), "fimH": ("Adherence", "VFC0001"),
    "mrkD": ("Adherence", "VFC0001"), "ybtS": ("Nutritional/Metabolic factor", "VFC0272"),
    "pilT": ("Adherence", "VFC0001"), "pilE": ("Adherence", "VFC0001"),
}
HIT_RATE = 0.0006  # expected hits per read per gene


# --- taxonomy tree ---------------------------------------------------------------------------
class Node:
    def __init__(self, rank, taxid, name):
        self.rank, self.taxid, self.name = rank, taxid, name
        self.children = {}
        self.lvl = {}   # sample -> directly assigned reads
        self.all = {}   # sample -> clade reads


def build_tree():
    root = Node("R", 1, "root")
    for lineage, taxid, name, *_ in SPECIES:
        cur = root
        for rank, tid, nm in lineage:
            cur = cur.children.setdefault(tid, Node(rank, tid, nm))
        cur.children[taxid] = Node("S", taxid, name)
    return root


def walk(node):
    yield node
    for ch in node.children.values():
        yield from walk(ch)


def sample_counts(rng):
    """Species-level kraken counts per sample (deterministic)."""
    depth = {s: math.exp(rng.gauss(0, 0.25)) for s in SAMPLES}
    depth["CTRL1"] = 0.4
    counts = {}
    for _l, taxid, _n, ma, mb, mc in SPECIES:
        for s in SAMPLES:
            mean = mc if s == "CTRL1" else (ma if GROUP[s] == "A" else mb)
            if mean == 0:
                counts[taxid, s] = 0
                continue
            # Kit taxa bleed into real samples below the control baseline
            val = mean * depth[s] * math.exp(rng.gauss(0, 0.35))
            counts[taxid, s] = max(0, int(round(val)))
    return counts, depth


def fill_kraken(root, species_counts, rng):
    for node in walk(root):
        for s in SAMPLES:
            node.lvl[s] = species_counts.get((node.taxid, s), 0) if node.rank == "S" else 0

    def post(node, s):
        child_sum = sum(post(c, s) for c in node.children.values())
        if node.rank != "S" and node.children and node.rank != "R":
            frac = INTERNAL_LVL_FRAC.get(node.rank, 0.001)
            node.lvl[s] = int(round(child_sum * frac * rng.uniform(0.5, 1.5)))
        node.all[s] = node.lvl[s] + child_sum
        return node.all[s]

    for s in SAMPLES:
        post(root, s)


def bracken_counts(root):
    """Bracken: all non-species reads are pushed down to species proportionally."""
    out = {}  # (taxid, s) -> species count after redistribution
    for s in SAMPLES:
        def push(node, extra):
            if node.rank == "S":
                out[node.taxid, s] = node.lvl[s] + int(round(extra))
                return
            total = node.lvl[s] + extra
            kids = list(node.children.values())
            weight = sum(k.all[s] for k in kids)
            for k in kids:
                share = total * (k.all[s] / weight) if weight else 0
                push(k, share)
        # root.lvl is 0 by construction
        push(root, 0)
    return out


def dfs_sorted(node, key):
    yield node
    for ch in sorted(node.children.values(), key=lambda n: (-key(n), n.taxid)):
        yield from dfs_sorted(ch, key)


def depth_of(root):
    d = {}

    def rec(n, k):
        d[n.taxid] = k
        for c in n.children.values():
            rec(c, k + 1)
    rec(root, 0)
    return d


# --- writers ---------------------------------------------------------------------------------
def _w(path, text):
    with open(path, "w", newline="\n", encoding="utf-8") as fh:
        fh.write(text)


def pct(x, total):
    return f"{100 * x / total:.2f}"


def write_kraken(root, unclass, path_raw, path_arranged):
    total = {s: root.all[s] + unclass[s] for s in SAMPLES}
    grand = sum(total.values())
    depth = depth_of(root)
    order = list(dfs_sorted(root, lambda n: sum(n.all.values())))
    n = len(SAMPLES)
    cols = []
    for i in range(1, n + 1):
        cols += [f"{i}_all", f"{i}_lvl"]
    header = "#perc\ttot_all\ttot_lvl\t" + "\t".join(cols) + "\tlvl_type\ttaxid\tname"
    lines = []
    u_all = sum(unclass.values())
    lines.append("\t".join([pct(u_all, grand), str(u_all), str(u_all)]
                           + sum(([str(unclass[s]), str(unclass[s])] for s in SAMPLES), [])
                           + ["U", "0", "unclassified"]))
    for node in order:
        tall = sum(node.all.values())
        tlvl = sum(node.lvl.values())
        if node.rank != "R" and tall == 0:
            continue
        lines.append("\t".join([pct(tall, grand), str(tall), str(tlvl)]
                               + sum(([str(node.all[s]), str(node.lvl[s])] for s in SAMPLES), [])
                               + [node.rank, str(node.taxid), "  " * depth[node.taxid] + node.name]))
    preamble = [f"#Number of Samples: {n}", f"#Total Number of Reads: {grand}"]
    preamble += [f"#{s}\treports/{s}.kraken2.report.txt" for s in SAMPLES]
    _w(path_raw, "\n".join(preamble + [header] + lines) + "\n")
    _w(path_arranged, "\n".join([header] + lines) + "\n")
    return order


def fmt_number(value):  # same as bin/arrange_reports.py
    value = round(value, 4)
    if float(value).is_integer():
        return str(int(value))
    return f"{value:.4f}".rstrip("0").rstrip(".")


def write_bracken(root, order, brk, path_raw, path_arranged):
    """Bracken kreport tree: species lvl = redistributed count, ancestors are sums."""
    bnode = {}
    for node in walk(root):
        bnode[node.taxid] = {s: 0 for s in SAMPLES}
    for node in walk(root):
        if node.rank == "S":
            for s in SAMPLES:
                bnode[node.taxid][s] = brk[node.taxid, s]
    parent = {}
    for node in walk(root):
        for c in node.children.values():
            parent[c.taxid] = node
    blvl = {t: dict(v) for t, v in bnode.items()}
    ball = {t: {s: 0 for s in SAMPLES} for t in bnode}
    for node in walk(root):
        if node.rank == "S":
            for s in SAMPLES:
                cur = node
                while True:
                    ball[cur.taxid][s] += blvl[node.taxid][s]
                    if cur.taxid not in parent:
                        break
                    cur = parent[cur.taxid]
    depth = depth_of(root)
    grand = sum(sum(ball[1].values()) for _ in [0])
    n = len(SAMPLES)
    cols = []
    for i in range(1, n + 1):
        cols += [f"{i}_all", f"{i}_lvl"]
    header = "#perc\ttot_all\ttot_lvl\t" + "\t".join(cols) + "\tlvl_type\ttaxid\tname"
    raw_lines = ["\t".join(["0.00", "0", "0", "U", "0", "unclassified"])]  # totals-only row
    rows = [["#perc", "tot_all"] + SAMPLES + ["taxid", "name"]]
    rows.append([fmt_number(0), "0"] + ["0"] * n + ["0", "unclassified"])
    for node in order:
        t = node.taxid
        tall = sum(ball[t].values())
        if node.rank != "R" and tall == 0:
            continue
        tlvl = sum(blvl[t].values())
        raw_lines.append("\t".join([pct(tall, grand), str(tall), str(tlvl)]
                                   + sum(([str(ball[t][s]), str(blvl[t][s])] for s in SAMPLES), [])
                                   + [node.rank, str(t), "  " * depth[t] + node.name]))
        rows.append([fmt_number(100 * tall / grand), fmt_number(tall)]
                    + [fmt_number(ball[t][s]) for s in SAMPLES] + [str(t), node.name])
    preamble = [f"#Number of Samples: {n}", f"#Total Number of Reads: {grand}"]
    preamble += [f"#{s}\treports/{s}_bracken_species.kreport" for s in SAMPLES]
    _w(path_raw, "\n".join(preamble + [header] + raw_lines) + "\n")
    _w(path_arranged, "".join(",".join(r) + "\n" for r in rows))


def write_rcf(root, order, rng, path):
    """Raw Recentrifuge 2.1.2 data CSV (post RECENTRIFUGE-module header normalisation)."""
    def trip(prefix_cell, stats=False):
        return [prefix_cell, "", ""]

    h, st = ["Samples"], ["Stats"]
    for s in SAMPLES:
        h += [f"kraken/{s}.kraken2.classifiedreads", "", ""]
        st += ["cnt", "una", "sco"]
    derived = []
    for rank in ("genus", "species"):
        for s in SAMPLES:
            derived.append(f"kraken/{s}.kraken2.classifiedreads_EXCLUSIVE_{rank}")
    for rank in ("genus", "species"):
        for s in SAMPLES:
            derived.append(f"kraken/{s}.kraken2.classifiedreads_SHARED_{rank}")
    for d in derived:
        h += [d, "", ""]
        st += ["cnt", "una", "sco"]
    h += ["SHARED_SUMMARY", "", "", "Details", "Details", "Name", "Rank"]
    st += ["cnt", "una", "sco", "Details", "Details", "Name", "Rank"]
    idrow = ["Id"] + [""] * (len(h) - 1)
    rows = [h, st, idrow]
    rank_name = {"G": "genus", "S": "species", "F": "family", "O": "order", "C": "class",
                 "P": "phylum", "D": "superkingdom", "R": "no rank", "R1": "no rank"}
    for node in order:
        if node.rank == "R":
            continue
        # Recentrifuge reports reads assigned at the node (cnt) with score (sco)
        cells = [str(node.taxid)]
        present = {s: node.lvl[s] > 0 for s in SAMPLES}
        if not any(present.values()):
            continue
        for s in SAMPLES:
            c = node.lvl[s]
            cells += [str(c), str(c), f"{rng.uniform(18, 32):.1f}" if c else "0.0"]
        for rank in ("genus", "species"):
            for s in SAMPLES:
                only = present[s] and sum(present.values()) == 1 and node.rank == rank[0].upper()
                c = node.lvl[s] if only else 0
                cells += [str(c), str(c), "0.0"]
        for rank in ("genus", "species"):
            for s in SAMPLES:
                shared = present[s] and sum(present.values()) > 1 and node.rank == rank[0].upper()
                c = node.lvl[s] if shared else 0
                cells += [str(c), str(c), "0.0"]
        tot = sum(node.lvl.values())
        cells += [str(tot), str(tot), "0.0", "", "", node.name, rank_name.get(node.rank, "no rank")]
        rows.append(cells)
    _w(path, "".join(",".join(f'"{c}"' if "," in c else c for c in r) + "\n" for r in rows))


def hash_gene(gene):
    """Stable small hash (python's hash() is salted per process)."""
    return sum((i + 1) * ord(c) for i, c in enumerate(gene))


def write_abricate(species_counts, rng, path):
    header = ("#FILE,SEQUENCE,START,END,STRAND,GENE,COVERAGE,COVERAGE_MAP,GAPS,"
              "%COVERAGE,%IDENTITY,DATABASE,ACCESSION,PRODUCT,RESISTANCE")
    lines = [header]
    for s in SAMPLES:
        if s == ZERO_HIT_SAMPLE:
            continue
        for taxid, genes in GENES.items():
            reads = species_counts.get((taxid, s), 0)
            for db, gene, acc, product, resist in genes:
                mean = reads * HIT_RATE
                n = max(0, int(round(rng.gauss(mean, math.sqrt(mean))))) if mean > 0 else 0
                if db == VFDB:
                    fn, vfc = VF_FUNC[gene]
                    prod = (f"({gene}) {product} [{gene} (VF{abs(hash_gene(gene)) % 9000 + 1000})] - "
                            f"{fn} ({vfc}) [{SPECIES_PROFILE[taxid][2]}]")
                else:
                    prod = product
                for _ in range(n):
                    rid = uuid.UUID(int=rng.getrandbits(128), version=4)
                    length = rng.randint(180, 1800)
                    start = rng.randint(1, 200)
                    end = start + length - 1
                    strand = rng.choice("+-")
                    cov = rng.uniform(88, 100)
                    ident = rng.uniform(86, 100)
                    lines.append(",".join([
                        f"{s}_kraken2ID_PFP.fasta", f"{rid}_kraken:taxid|{taxid}",
                        str(start), str(end), strand, gene, f"1-{length}/{length}",
                        "=" * 12, "0/0", f"{cov:.2f}", f"{ident:.2f}", db, acc, prod, resist]))
    _w(path, "\n".join(lines) + "\n")
    return len(lines) - 1


def write_metadata(path):
    rows = ["sample,group,ID"]
    for i in range(1, 9):
        s = f"S{i:02d}"
        rows.append(f"{s},{GROUP[s]},P{i:02d}")
    rows.append("CTRL1,Control,CTRL")
    _w(path, "\n".join(rows) + "\n")


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    out = Path(argv[0]) if argv else Path(__file__).resolve().parent / "data" / "meridian_fixture"
    (out / "raw").mkdir(parents=True, exist_ok=True)
    rng = random.Random(SEED)
    counts, _depth = sample_counts(rng)
    root = build_tree()
    fill_kraken(root, counts, rng)
    # classified totals are the root sums; unclassified ~ 20% of reads
    unclass = {s: int(round(root.all[s] * rng.uniform(0.15, 0.3))) for s in SAMPLES}
    order = write_kraken(root, unclass, out / "raw" / "kraken_combined_raw.txt",
                         out / "kraken2_combined_reports.txt")
    brk = bracken_counts(root)
    write_bracken(root, order, brk, out / "raw" / "bracken_combined_raw.txt",
                  out / "bracken_arranged.csv")
    write_rcf(root, order, rng, out / "pneuma.rcf.data.csv")
    n_hits = write_abricate(counts, rng, out / "summary_reportClean.csv")
    write_metadata(out / "metadata.csv")
    print(f"wrote fixture to {out} ({n_hits} ABRicate hits)")


if __name__ == "__main__":
    main()
