#!/usr/bin/env python3
"""Generate docs/images/pneuma_metro_map.svg (nf-core "tube map" style).

Stage columns (faint alternating bands, small uppercase labels), a blue QC
trunk into the Kraken2 interchange, three coloured branch tracks (abundance,
contamination, functional) that converge on the MERIDIAN interchange and the
Report terminus, a planned (grey, dashed) pathogen line, and a thin MultiQC
bus along the bottom fed by short drop lines. Standard library only; output
is deterministic, so re-running must leave the committed SVG unchanged:

    python docs/images/make_metro_map.py
"""
from pathlib import Path

OUT = Path(__file__).with_name("pneuma_metro_map.svg")

BG = "#f6f8fa"       # light backing panel (readable on light and dark GitHub themes)
BAND = "#eceff3"     # alternating column band
GUIDE = "#d0d7de"
INK = "#1f2328"
MUTED = "#57606a"
GREY = "#8c959f"
FONT = "Helvetica, Arial, sans-serif"

C_QC, C_AB, C_CO = "#0072B2", "#009E73", "#D55E00"
C_FN, C_DS = "#9B3FA8", "#D4A017"

W, H = 1400, 610
LW = 6                      # line width
SZ, SZS, SZN = 13, 11, 12   # station / stage / note font sizes

# track rows
Y0 = 205                    # QC trunk + abundance (green)
YC, YF, YP = 270, 335, 400  # contamination, functional, planned
BUS = 468                   # MultiQC bus

out = []


def add(s):
    out.append(s)


def rpath(pts, r=14):
    """Polyline through pts with rounded (quadratic) corners."""
    d = f"M{pts[0][0]} {pts[0][1]}"
    for i in range(1, len(pts) - 1):
        (x0, y0), (x1, y1), (x2, y2) = pts[i - 1], pts[i], pts[i + 1]

        def toward(ax, ay, bx, by):
            dx, dy = bx - ax, by - ay
            n = (dx * dx + dy * dy) ** 0.5
            return ax + dx / n * r, ay + dy / n * r

        sx, sy = toward(x1, y1, x0, y0)
        ex, ey = toward(x1, y1, x2, y2)
        d += f" L{sx:.1f} {sy:.1f} Q{x1} {y1} {ex:.1f} {ey:.1f}"
    d += f" L{pts[-1][0]} {pts[-1][1]}"
    return d


def line(pts, col, dashed=False, lw=LW):
    dash = ' stroke-dasharray="11 8"' if dashed else ""
    add(f'<path d="{rpath(pts)}" fill="none" stroke="{col}" stroke-width="{lw}" '
        f'stroke-linecap="round" stroke-linejoin="round"{dash}/>')


def text(x, y, lines, anchor="middle", size=SZ, weight="400", fill=INK,
         style="", lh=15, spacing=None):
    tsp = "".join(f'<tspan x="{x}" dy="{0 if i == 0 else lh}">{t}</tspan>'
                  for i, t in enumerate(lines))
    ls = f' letter-spacing="{spacing}"' if spacing else ""
    st = f' font-style="{style}"' if style else ""
    add(f'<text x="{x}" y="{y}" text-anchor="{anchor}" font-size="{size}" '
        f'font-weight="{weight}" fill="{fill}"{ls}{st} stroke="{BG}" '
        f'stroke-width="3.5" paint-order="stroke" stroke-linejoin="round">{tsp}</text>')


def label(x, y, lines, side=-1, **kw):
    """Station label above (side=-1) or below (side=+1) a station at (x, y)."""
    kw.setdefault("weight", "700")
    if side < 0:
        text(x, y - 16 - 15 * (len(lines) - 1), lines, **kw)
    else:
        text(x, y + 26, lines, **kw)


def stage(x, y, lab, anchor="middle"):
    text(x, y, [lab], anchor=anchor, size=SZS, weight="700", fill=MUTED, spacing="1.2")


def station(x, y, kind="", r=6.5):
    """kind: '' process | 'opt' optional (dotted) | 'plan' planned (hollow dotted)."""
    if kind == "opt":
        add(f'<circle cx="{x}" cy="{y}" r="{r}" fill="#ffffff" stroke="{INK}" '
            f'stroke-width="2.5" stroke-dasharray="1.5 2.5" stroke-linecap="round"/>')
    elif kind == "plan":
        add(f'<circle cx="{x}" cy="{y}" r="{r}" fill="{BG}" stroke="{GREY}" '
            f'stroke-width="2.5" stroke-dasharray="1.5 2.5" stroke-linecap="round"/>')
    else:
        add(f'<circle cx="{x}" cy="{y}" r="{r}" fill="#ffffff" stroke="{INK}" '
            f'stroke-width="2.5"/>')


def interchange(cx, y, w=100, h=24):
    add(f'<rect x="{cx - w / 2:.0f}" y="{y - h / 2:.0f}" width="{w}" height="{h}" '
        f'rx="{h / 2:.0f}" fill="#ffffff" stroke="{INK}" stroke-width="3"/>')


def terminus(x, y, col=INK):
    add(f'<circle cx="{x}" cy="{y}" r="9" fill="#ffffff" stroke="{col}" stroke-width="3"/>')
    add(f'<circle cx="{x}" cy="{y}" r="3.5" fill="{col}"/>')


def drop(x, y):
    """Thin grey drop line from a station down to the MultiQC bus."""
    add(f'<line x1="{x}" y1="{y + 8}" x2="{x}" y2="{BUS}" stroke="{GREY}" '
        f'stroke-width="2.5" stroke-linecap="round"/>')
    add(f'<circle cx="{x}" cy="{BUS}" r="3.5" fill="{GREY}"/>')


# ---------------------------------------------------------------- geometry
TX = 45                               # input terminus
KX = 625                              # Kraken2 interchange centre
MX = 1290                             # MERIDIAN interchange centre
RX = 1362                             # Report terminus
MQX = 1290                            # MultiQC terminus
PTX = 1330                            # planned Pathogen report terminus
BX0 = 50                              # band / guide extents
BT, BB = 70, 515
COLS = [(15, 560, "QC &amp; HOST REMOVAL"), (560, 700, "TAXONOMY"),
        (700, 965, "ABUNDANCE"), (965, 1090, "CONTAMINATION"),
        (1090, 1255, "FUNCTIONAL"), (1255, 1385, "DOWNSTREAM")]

# ---------------------------------------------------------------- header
add(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" '
    f'height="{H}" role="img" font-family="{FONT}">')
add("<title>PNEUMA pipeline metro map</title>")
add("<desc>Tube-map overview of the PNEUMA pipeline in six stage columns. A blue QC "
    "and host removal trunk (cat, Porechop_ABI, nanoq, NanoPlot, SeqKit stats, "
    "Hostile) reaches the Kraken2 interchange, which branches into Abundance "
    "(Bracken, KrakenTools combine, Krona), Contamination (Recentrifuge, negative-"
    "control aware) and Functional profiling (Tag reads, ABRicate with CARD, VFDB "
    "and PlasmidFinder). The three tracks converge on the MERIDIAN interchange "
    "and the Report. A grey dashed planned line (assembly, CheckM2, Bakta, MLST, "
    "Pathogen report) also leaves Kraken2. A thin MultiQC bus along the bottom "
    "collects QC, Kraken2/Bracken, Recentrifuge and ABRicate outputs into the QC "
    "report. Reference databases are built with --build_databases or supplied "
    "with --db_dir.</desc>")
add(f'<rect x="1" y="1" width="{W - 2}" height="{H - 2}" rx="18" fill="{BG}" '
    f'stroke="#d0d7de" stroke-width="2"/>')
text(30, 42, ["PNEUMA: ONT metagenomics for low-biomass respiratory samples"],
     anchor="start", size=20, weight="700")

# stage columns: faint alternating bands + guide lines + labels
for i, (x0, x1, lab) in enumerate(COLS):
    if i % 2 == 0:
        add(f'<rect x="{x0}" y="{BT}" width="{x1 - x0}" height="{BB - BT}" '
            f'fill="{BAND}" opacity="0.55"/>')
    if i > 0:
        add(f'<line x1="{x0}" y1="{BT}" x2="{x0}" y2="{BB}" stroke="{GUIDE}" '
            f'stroke-width="1" opacity="0.7"/>')
for x0, x1, lab in COLS:
    stage((x0 + x1) / 2, BT + 20, lab)

# ---------------------------------------------------------------- lines
PX = MX - 14                          # left edge of the tall MERIDIAN pill
D = 30                                # bend start offset right of Kraken2 centre
line([(TX, Y0), (KX, Y0)], C_QC)
line([(KX, Y0), (PX, Y0)], C_AB)                                        # abundance
line([(KX, Y0), (KX + D, Y0), (KX + D + YP - Y0, YP), (PTX, YP)], GREY,  # planned
     dashed=True)
line([(KX, Y0), (KX + D, Y0), (KX + D + YF - Y0, YF), (PX, YF)], C_FN)   # functional
line([(KX, Y0), (KX + D, Y0), (KX + D + YC - Y0, YC), (PX, YC)], C_CO)   # contamination
line([(MX, Y0), (RX, Y0)], C_DS)                                        # downstream
line([(280, BUS), (MQX, BUS)], GREY, lw=2.5)                            # MultiQC bus

# ---------------------------------------------------------------- QC trunk
terminus(TX, Y0)
text(20, Y0 + 28, ["samplesheet +", "ONT FASTQ"], anchor="start", size=SZN, fill=MUTED)
qc = [("cat", 125, "", ["cat"]), ("Porechop_ABI", 205, "opt", ["Porechop_ABI"]),
      ("nanoq", 285, "", ["nanoq"]), ("NanoPlot", 365, "", ["NanoPlot"]),
      ("SeqKit stats", 445, "", ["SeqKit", "stats"]), ("Hostile", 525, "", ["Hostile"])]
for name, x, kind, lines in qc:
    station(x, Y0, kind)
    label(x, Y0, lines)
    if x >= 285:
        drop(x, Y0)
text(525, Y0 - 34, ["LOW_DEPTH flag"], size=SZN, fill=MUTED, style="italic")

interchange(KX, Y0)
text(KX, Y0 - 40, ["Kraken2", "(PlusPF)"], weight="700")
drop(KX, Y0)

# ---------------------------------------------------------------- abundance
for lines, x in [(["Bracken"], 760), (["KrakenTools", "combine"], 850), (["Krona"], 940)]:
    station(x, Y0)
    label(x, Y0, lines)

# ---------------------------------------------------------------- contamination
RCX = 1025
station(RCX, YC)
label(RCX, YC, ["Recentrifuge"])
text(RCX, YC - 34, ["negative-control aware"], size=SZN, fill=MUTED, style="italic")
drop(RCX, YC)

# ---------------------------------------------------------------- functional
for name, x in [("Tag reads", 1135), ("ABRicate", 1205)]:
    station(x, YF)
    label(x, YF, [name])
text(1170, YF - 36, ["CARD · VFDB · PlasmidFinder"], size=SZN, fill=MUTED, style="italic")
drop(1205, YF)

# ---------------------------------------------------------------- downstream
add(f'<rect x="{MX - 14}" y="{Y0 - 14}" width="28" height="{YF - Y0 + 28}" rx="14" '
    f'fill="#ffffff" stroke="{INK}" stroke-width="3"/>')
text(MX, Y0 - 24, ["MERIDIAN"], size=SZ, weight="700")
terminus(RX, Y0, C_DS)
text(RX, Y0 - 24, ["Report"], size=SZ, weight="700")

# ---------------------------------------------------------------- planned line
for lines, x in [(["binning /", "assembly"], 905), (["CheckM2"], 975), (["Bakta"], 1085),
                 (["MLST"], 1150)]:
    station(x, YP, "plan")
    label(x, YP, lines, weight="400", fill=MUTED)
terminus(PTX, YP, GREY)
text(PTX, YP - 24, ["Pathogen report"], size=SZ, fill=MUTED)
text(1085, YP + 26, ["planned"], size=SZN, fill=MUTED, style="italic")

# ---------------------------------------------------------------- MultiQC terminus
terminus(MQX, BUS, GREY)
text(MQX, BUS + 28, ["MultiQC"], size=SZ, weight="700")
text(MQX, BUS + 44, ["QC report (all stages)"], size=SZN, fill=MUTED, style="italic")

# ---------------------------------------------------------------- legend
LY = 552
x = 30
leg = [("QC &amp; taxonomy", C_QC, "solid"), ("Abundance", C_AB, "solid"),
       ("Contamination", C_CO, "solid"), ("Functional (AMR / VF / MGE)", C_FN, "solid"),
       ("Downstream", C_DS, "solid"), ("Planned", GREY, "dash"),
       ("MultiQC bus", GREY, "thin")]
for lab, col, kind in leg:
    dash = ' stroke-dasharray="8 6"' if kind == "dash" else ""
    lw = 2.5 if kind == "thin" else LW
    add(f'<line x1="{x}" y1="{LY}" x2="{x + 34}" y2="{LY}" stroke="{col}" '
        f'stroke-width="{lw}" stroke-linecap="round"{dash}/>')
    text(x + 44, LY + 5, [lab], anchor="start", size=SZ)
    x += 44 + int(len(lab.replace("&amp;", "&")) * 7.3) + 30
SY = LY + 32
x = 30
for lab, kind in [("Process", ""), ("Optional step", "opt"), ("Planned step", "plan")]:
    station(x + 8, SY, kind)
    text(x + 24, SY + 5, [lab], anchor="start", size=SZ)
    x += 24 + int(len(lab) * 7.3) + 34
interchange(x + 24, SY, w=44)
text(x + 56, SY + 5, ["Interchange"], anchor="start", size=SZ)
x += 56 + int(len("Interchange") * 7.3) + 40
text(x, SY + 5, ["Databases: --build_databases to download, or --db_dir to reuse"],
     anchor="start", size=SZN, fill=MUTED, style="italic")

add("</svg>")

OUT.write_text("\n".join(out) + "\n", encoding="utf-8", newline="\n")
print(f"wrote {OUT}")
