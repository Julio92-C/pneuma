#!/usr/bin/env python3
"""Generate one SVG "module card" per PNEUMA subworkflow into docs/images/cards/.

A card reads left to right: the inputs (data objects drawn by kind: reads as pills, sequences as parallelograms,
tables as boxes with a header band, databases as cylinders), the tool steps from top to bottom (tool, the question it
answers, its defaults above it, its blind spot underneath; optional steps dashed, default steps solid), and the
outputs. Same palette and font as make_metro_map.py. Standard library only; output is deterministic, so re-running
must leave the committed SVGs unchanged:

    python docs/images/make_module_cards.py
"""
import argparse
import os
from pathlib import Path
from xml.sax.saxutils import escape

OUT_DIR = Path(__file__).with_name("cards")
BG, INK, MUTED, FONT = "#f6f8fa", "#1f2328", "#57606a", "Helvetica, Arial, sans-serif"
C_QC, C_AB, C_CO, C_FN = "#0072B2", "#009E73", "#D55E00", "#9B3FA8"      # line colours of make_metro_map.py
C_DS, C_PA, C_GC = "#D4A017", "#CC79A7", "#56B4E9"
W, TOP, STEP_H, BOX_H, OBJ_H, OBJ_GAP = 1100, 110, 100, 46, 40, 16
IN_X, IN_W, ST_X, ST_W, OUT_X, OUT_W = 24, 220, 300, 500, 856, 220


def S(tool, question, default="", blind="", optional=False):
    return {"tool": tool, "question": question, "default": default, "blind": blind, "optional": optional}


CARDS = [
    {"slug": "preprocess", "title": "Preprocess: QC and host removal", "colour": C_QC,
     "inputs": [("reads", "ONT reads (FASTQ, per run)"), ("db", "Host index (Hostile)")],
     "steps": [S("cat", "one FASTQ per sample: its runs merged", "", "a run's own quality is only seen before merging"),
               S("Porechop_ABI", "are adapters left on the reads?", "--run_porechop",
                 "adapters only; no barcode demultiplexing", optional=True),
               S("nanoq", "is each read long and accurate enough?", "--nanoq_min_length 100 · --nanoq_min_qscore 7",
                 "mean read quality; local errors stay"),
               S("Hostile", "which reads are human?", "human-t2t-hla index",
                 "only the host in the index (--host_index for others)"),
               S("NanoPlot · SeqKit stats", "what is left after QC and host removal?", "LOW_DEPTH below --min_reads 100",
                 "summaries only; nothing is filtered here")],
     "outputs": [("reads", "Clean reads"), ("table", "QC tables (MultiQC)")]},
    {"slug": "taxonomy", "title": "Taxonomy and abundance", "colour": C_AB,
     "inputs": [("reads", "Clean reads"), ("db", "Kraken2 PlusPF index"), ("db", "NCBI taxdump")],
     "steps": [S("Kraken2", "which taxon is each read?", "--kraken2_confidence 0",
                 "taxa missing from the index go to a relative or unclassified"),
               S("Bracken", "how abundant is each species?", "--bracken_level S · --bracken_threshold 10",
                 "re-estimated counts, not proof of presence"),
               S("KrakenTools · Krona", "the cohort tables and the interactive plot")],
     "outputs": [("table", "Kraken2 / Bracken tables"), ("reads", "Kraken2-classified reads"), ("table", "Krona HTML")]},
    {"slug": "contamination", "title": "Contamination", "colour": C_CO,
     "inputs": [("table", "Kraken2 per-read output"), ("db", "NCBI taxdump"), ("table", "Negative controls")],
     "steps": [S("Recentrifuge", "which taxa are also in the negative controls?", "is_control = true in the samplesheet",
                 "needs controls; a run without one only gets a warning")],
     "outputs": [("table", "pneuma.rcf.data.csv"), ("table", "pneuma_controls.rcf.*"), ("table", "Recentrifuge HTML")]},
    {"slug": "functional", "title": "Functional profiling", "colour": C_FN,
     "inputs": [("reads", "Kraken2-classified reads"), ("db", "CARD · VFDB · PlasmidFinder")],
     "steps": [S("Tag reads", "which Kraken2 taxon does each read carry?", "taxid written into the read id",
                 "one taxon per read, also for a read that carries a plasmid"),
               S("ABRicate", "which ARG, VF and replicon genes are on each read?",
                 "--abricate_minid 80 · --abricate_mincov 80",
                 "known genes only; CARD variant models: mutation not checked")],
     "outputs": [("table", "<sample>_abricate.csv"), ("table", "summary_reportClean.csv")]},
    {"slug": "genomic_context", "title": "Genomic context", "colour": C_GC,
     "inputs": [("reads", "Kraken2-tagged hit reads"), ("table", "ABRicate hits"), ("reads", "Clean reads"),
                ("db", "geNomad database")],
     "steps": [S("Hit reads", "which reads carry a hit, and which are seeds?",
                 "--context_databases card,vfdb,plasmidfinder", "CARD variant models are not recruitment seeds"),
               S("geNomad", "does each hit read look like plasmid, virus or chromosome?",
                 "--context_plasmid_min 0.7 · --context_chrom_min 0.7",
                 "scores alone call chromosomal operons plasmid: evidence rule"),
               S("ARG context (reads)", "where is each gene: plasmid, MGE or chromosome?",
                 "--context_min_flank 1000 · --context_mge_window 10000",
                 "most single reads end too close to the gene: unresolved"),
               S("Read recruitment (minimap2)", "which reads overlap the seeds? (two rounds)",
                 "--context_recruit_min_aln 1000", "a region with few reads cannot be extended"),
               S("Flye --meta", "can the region around the hit be assembled?", "--flye_mode --nano-hq",
                 "low-coverage samples fail and keep their read calls"),
               S("Contig typing", "ABRicate, geNomad and MOB-suite on the contigs"),
               S("Extend context", "which contig call does each read hit take?", "identity >= 0.8 · CIGAR projection",
                 "an unresolved contig never overrides a read call"),
               S("Context summary", "which taxon, and how strong is the host vote?", "weak_vote: share < 0.5 or < 2 reads",
                 "Kraken2 names the host of a plasmid gene weakly")],
     "outputs": [("table", "context/hit_context.tsv"), ("table", "context/summary/"), ("table", "loci_* tables"),
                 ("table", "MultiQC: Genomic context")]},
    {"slug": "mobilome", "title": "Mobilome", "colour": C_GC,
     "inputs": [("seqs", "Local-assembly contigs"), ("reads", "Recruited reads"), ("db", "MOB-suite database"),
                ("db", "Medaka model")],
     "steps": [S("Circular plasmids", "which circular contigs carry a hit?", "Flye circ. = Y",
                 "linear plasmid pieces are not drafts; repeat contigs flagged"),
               S("Medaka", "polish the drafts against all of the sample's contigs", "--polish_plasmids true",
                 "polishing does not join or close a draft"),
               S("MOB-suite typer", "replicon, relaxase, MPF and mobility of each contig", "", "typing only: no host"),
               S("ISEScan · IS flanks", "which IS elements flank each ARG? composite transposon?",
                 "--context_is_window 500", "IS elements beyond the contig end are missed"),
               S("Plasmid units", "is it the same plasmid in other samples?", "--plasmid_unit_min_identity 0.95",
                 "drafts under 80% of each other's length never cluster")],
     "outputs": [("table", "mobilome/plasmid_units.tsv"), ("table", "mobilome/is_flanks.tsv"),
                 ("seqs", "Polished plasmid drafts"), ("table", "MultiQC: plasmid units")]},
    {"slug": "pathogens", "title": "Pathogen reconstruction (opt-in)", "colour": C_PA,
     "inputs": [("table", "Kraken2 assignments"), ("reads", "Kraken2-classified reads"),
                ("table", "Pathogen panel (WHO priority)"), ("db", "CheckM2 · Bakta · AMRFinderPlus")],
     "steps": [S("Panel candidates", "which panel species have enough reads?",
                 "--pathogen_min_reads 10 · --pathogen_min_est_depth 2.5", "species outside the panel are not rebuilt"),
               S("KrakenTools extract", "the reads of each candidate clade", "--include-children",
                 "reads Kraken2 put elsewhere are left behind"),
               S("Coverage gate (minimap2)", "is there enough depth and breadth to assemble?",
                 "--gate_min_depth 10 · --gate_min_breadth 0.5", "patient-level pooling only when no sample passes"),
               S("Flye · Medaka", "a draft genome per sample x species", "--flye_mode --nano-hq",
                 "mixed strains collapse into one consensus"),
               S("CheckM2 · QUAST", "how complete and how contaminated is the draft?"),
               S("skani", "is the draft the panel species or a close relative?", "ani_min per species",
                 "confusables within 0.5 ANI points stay unresolved"),
               S("MLST · AMRFinderPlus · ABRicate · Bakta", "sequence type, resistance, virulence, annotation"),
               S("MOB-suite", "the plasmids of the draft genome", "--run_mobsuite", "", optional=True)],
     "outputs": [("table", "pathogen_cohort.tsv"), ("table", "Per-sample HTML report"), ("seqs", "Pathogen assemblies")]},
    {"slug": "meridian", "title": "MERIDIAN hand-off", "colour": C_DS,
     "inputs": [("table", "summary_reportClean.csv"), ("table", "Kraken2 / Bracken tables"),
                ("table", "pneuma.rcf.data.csv"), ("table", "Metadata (--meridian_metadata)"),
                ("table", "context/hit_context.tsv")],
     "steps": [S("MERIDIAN config", "is each --meridian_group column usable?", "one comparison per column",
                 "checks the metadata columns, not the study design"),
               S("MERIDIAN (container)", "diversity, differential abundance, resistome, networks",
                 "--meridian_container v1.3.1", "does not read hit_context.tsv yet; one TPM denominator for all")],
     "outputs": [("table", "meridian/<group>/Reports/"), ("table", "MERIDIAN HTML report")]},
]


def text(x, y, s, size=12, weight="400", fill=INK, anchor="start", italic=False):
    st = ' font-style="italic"' if italic else ""
    return ('<text x="%d" y="%d" font-size="%s" font-weight="%s" fill="%s" text-anchor="%s"%s>%s</text>'
            % (x, y, size, weight, fill, anchor, st, escape(s)))


def data_object(kind, x, y, colour, label, w=IN_W, h=OBJ_H):
    """An input / output: reads = pill, seqs = parallelogram, table = box with a header band, db = cylinder."""
    if kind == "reads":
        parts = ['<rect x="%d" y="%d" width="%d" height="%d" rx="%d" fill="#ffffff" stroke="%s" stroke-width="2"/>'
                 % (x, y, w, h, h // 2, colour)]
    elif kind == "seqs":
        parts = ['<polygon points="%d,%d %d,%d %d,%d %d,%d" fill="#ffffff" stroke="%s" stroke-width="2"/>'
                 % (x + 12, y, x + w, y, x + w - 12, y + h, x, y + h, colour)]
    elif kind == "table":
        parts = ['<rect x="%d" y="%d" width="%d" height="%d" fill="#ffffff" stroke="%s" stroke-width="2"/>'
                 % (x, y, w, h, colour),
                 '<rect x="%d" y="%d" width="%d" height="8" fill="%s"/>' % (x, y, w, colour)]
    elif kind == "db":
        e = 6
        parts = ['<path d="M%d,%d v%d a%d,%d 0 0 0 %d,0 v-%d" fill="#ffffff" stroke="%s" stroke-width="2"/>'
                 % (x, y + e, h - 2 * e, w // 2, e, w, h - 2 * e, colour),
                 '<ellipse cx="%d" cy="%d" rx="%d" ry="%d" fill="#ffffff" stroke="%s" stroke-width="2"/>'
                 % (x + w // 2, y + e, w // 2, e, colour)]
    else:
        raise ValueError("unknown data object kind %r" % kind)
    parts.append(text(x + w // 2, y + h // 2 + (7 if kind == "table" else 4), label, anchor="middle"))
    return parts


def arrow(x1, y1, x2, y2):
    return ('<line x1="%d" y1="%d" x2="%d" y2="%d" stroke="%s" stroke-width="1.6" marker-end="url(#arrow)"/>'
            % (x1, y1, x2, y2, MUTED))


def stack_top(n, centre):
    """Top y of the first of n stacked data objects centred on `centre`, never above TOP."""
    return max(TOP, centre - (n * (OBJ_H + OBJ_GAP) - OBJ_GAP) // 2)


def render(card):
    steps, colour = card["steps"], card["colour"]
    n_in, n_out = len(card["inputs"]), len(card["outputs"])
    last_y = TOP + (len(steps) - 1) * STEP_H
    out_top = stack_top(n_out, last_y + BOX_H // 2)
    body = max(last_y + STEP_H, TOP + n_in * (OBJ_H + OBJ_GAP), out_top + n_out * (OBJ_H + OBJ_GAP))
    h = body + 50
    out = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 %d %d" width="%d" height="%d" role="img" '
           'font-family="%s">' % (W, h, W, h, FONT),
           "<title>%s</title>" % escape("PNEUMA module card: " + card["title"]),
           '<defs><marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" '
           'orient="auto-start-reverse"><path d="M0,0 L10,5 L0,10 z" fill="%s"/></marker></defs>' % MUTED,
           '<rect x="1" y="1" width="%d" height="%d" rx="16" fill="%s" stroke="#d0d7de" stroke-width="2"/>'
           % (W - 2, h - 2, BG),
           '<rect x="24" y="22" width="10" height="28" rx="3" fill="%s"/>' % colour,
           text(44, 43, card["title"], size=18, weight="700")]
    for x, w, label in ((IN_X, IN_W, "INPUTS"), (ST_X, ST_W, "STEPS"), (OUT_X, OUT_W, "OUTPUTS")):
        out.append(text(x + w // 2, TOP - 30, label, size=10, weight="700", fill=MUTED, anchor="middle"))
    for j, (kind, label) in enumerate(card["inputs"]):
        y = TOP + j * (OBJ_H + OBJ_GAP)
        out += data_object(kind, IN_X, y, colour, label)
        out.append(arrow(IN_X + IN_W, y + OBJ_H // 2, ST_X - 4, TOP + BOX_H // 2))
    for i, s in enumerate(steps):
        y = TOP + i * STEP_H
        if i:
            out.append(arrow(ST_X + 40, y - STEP_H + BOX_H + 24, ST_X + 40, y - 3))
        if s["default"]:
            out.append(text(ST_X + ST_W, y - 7, s["default"], size=11, fill=colour, anchor="end", italic=True))
        dash = ' stroke-dasharray="7 5"' if s["optional"] else ""
        out.append('<rect x="%d" y="%d" width="%d" height="%d" rx="8" fill="#ffffff" stroke="%s" stroke-width="2.5"%s/>'
                   % (ST_X, y, ST_W, BOX_H, colour, dash))
        out.append(text(ST_X + 12, y + 19, s["tool"], size=13, weight="700"))
        out.append(text(ST_X + 12, y + 37, s["question"]))
        if s["blind"]:
            out.append(text(ST_X + 12, y + BOX_H + 16, "blind spot: " + s["blind"], size=11, fill=MUTED, italic=True))
    for j, (kind, label) in enumerate(card["outputs"]):
        y = out_top + j * (OBJ_H + OBJ_GAP)
        out.append(arrow(ST_X + ST_W, last_y + BOX_H // 2, OUT_X - 4, y + OBJ_H // 2))
        out += data_object(kind, OUT_X, y, colour, label, w=OUT_W)
    out.append(text(24, h - 18, "Inputs and outputs: reads = pill, sequences = parallelogram, table = banded box, "
                    "database = cylinder. Dashed step = optional (off by default).", size=11, fill=MUTED))
    out.append("</svg>")
    return "\n".join(out) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--outdir", default=str(OUT_DIR))
    a = ap.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    for card in CARDS:
        with open(os.path.join(a.outdir, card["slug"] + ".svg"), "w", encoding="utf-8", newline="\n") as fh:
            fh.write(render(card))
    print("wrote %d cards to %s" % (len(CARDS), a.outdir))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
