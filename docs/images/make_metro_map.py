#!/usr/bin/env python3
"""Generate docs/images/pneuma_metro_map.svg (nf-core "tube map" style).

Coloured transit lines with stations, 45-degree branch / merge bends and a
wider interchange at Kraken2 and MERIDIAN; small uppercase stage labels
replace the old boxed blocks. Standard library only; output is
deterministic, so re-running must leave the committed SVG unchanged:

    python docs/images/make_metro_map.py
"""
from pathlib import Path

OUT = Path(__file__).with_name("pneuma_metro_map.svg")

BG = "#f6f8fa"       # light backing panel (readable on light and dark GitHub themes)
INK = "#1f2328"
MUTED = "#57606a"
GREY = "#8c959f"
FONT = "Helvetica, Arial, sans-serif"

C_QC, C_AB, C_CO = "#0072B2", "#009E73", "#D55E00"
C_FN, C_DS = "#9B3FA8", "#D4A017"
C_MQ = "#57606a"

W, H = 1330, 486
LW = 6                      # line width

# lane geometry
Y0 = 240                    # trunk / functional / downstream line
YA, YC, YR = 96, 150, 350  # abundance, contamination, reconstruction lanes

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


def line(pts, col, dashed=False):
    dash = ' stroke-dasharray="11 8"' if dashed else ""
    add(f'<path d="{rpath(pts)}" fill="none" stroke="{col}" stroke-width="{LW}" '
        f'stroke-linecap="round" stroke-linejoin="round"{dash}/>')


def text(x, y, lines, anchor="middle", size=12, weight="400", fill=INK,
         style="", lh=13, spacing=None):
    tsp = "".join(f'<tspan x="{x}" dy="{0 if i == 0 else lh}">{t}</tspan>'
                  for i, t in enumerate(lines))
    ls = f' letter-spacing="{spacing}"' if spacing else ""
    st = f' font-style="{style}"' if style else ""
    add(f'<text x="{x}" y="{y}" text-anchor="{anchor}" font-size="{size}" '
        f'font-weight="{weight}" fill="{fill}"{ls}{st} stroke="{BG}" '
        f'stroke-width="3.5" paint-order="stroke" stroke-linejoin="round">{tsp}</text>')


def stage(x, y, label, anchor="start"):
    text(x, y, [label], anchor=anchor, size=9.5, weight="700", fill=MUTED, spacing="1.2")


def station(x, y, kind="", r=6.5, fill="#ffffff"):
    """kind: '' process | 'opt' optional (dotted) | 'plan' planned (hollow dotted)."""
    if kind == "opt":
        add(f'<circle cx="{x}" cy="{y}" r="{r}" fill="#ffffff" stroke="{INK}" '
            f'stroke-width="2.5" stroke-dasharray="1.5 2.5" stroke-linecap="round"/>')
    elif kind == "plan":
        add(f'<circle cx="{x}" cy="{y}" r="{r}" fill="{BG}" stroke="{GREY}" '
            f'stroke-width="2.5" stroke-dasharray="1.5 2.5" stroke-linecap="round"/>')
    else:
        add(f'<circle cx="{x}" cy="{y}" r="{r}" fill="{fill}" stroke="{INK}" '
            f'stroke-width="2.5"/>')


def interchange(cx, y, w=84, fill="#ffffff"):
    add(f'<rect x="{cx - w / 2:.0f}" y="{y - 12}" width="{w}" height="24" rx="12" '
        f'fill="{fill}" stroke="{INK}" stroke-width="3"/>')


def terminus(x, y, col=INK):
    add(f'<circle cx="{x}" cy="{y}" r="9" fill="#ffffff" stroke="{col}" stroke-width="3"/>')
    add(f'<circle cx="{x}" cy="{y}" r="3.5" fill="{col}"/>')


# ---------------------------------------------------------------- header
add(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" '
    f'height="{H}" role="img" font-family="{FONT}">')
add("<title>PNEUMA pipeline metro map</title>")
add("<desc>Tube-map overview of the PNEUMA pipeline. A blue QC and host removal "
    "trunk (cat, Porechop_ABI, nanoq, NanoPlot, SeqKit stats, Hostile) reaches the "
    "Kraken2 interchange which branches into Abundance (Bracken, KrakenTools "
    "combine, Krona), Contamination (Recentrifuge), Functional profiling (Tag "
    "reads, ABRicate, CARD, VFDB, PlasmidFinder) and a planned Metagenome "
    "reconstruction line. The first three merge at MERIDIAN for downstream "
    "analysis and the Report. A short spur just before MERIDIAN ends at MultiQC (report of all steps); "
    "filled stops mark the stations whose output it reports.</desc>")
add(f'<rect x="1" y="1" width="{W - 2}" height="{H - 2}" rx="18" fill="{BG}" '
    f'stroke="#d0d7de" stroke-width="2"/>')
text(30, 42, ["PNEUMA: ONT metagenomics for low-biomass respiratory samples"],
     anchor="start", size=18, weight="700")

# ---------------------------------------------------------------- geometry
KX = 505                              # Kraken2 interchange centre
MX = 1195                             # MERIDIAN interchange centre
M0 = MX - 42                          # its left edge
RX = 1296                             # Report terminus
TX = 60                               # input terminus
MQD = 36                              # MultiQC spur drop
MQX, MQY = 1268, Y0 + MQD             # MultiQC terminus

# lines (drawn first, stations on top)
line([(TX, Y0), (KX - 20, Y0)], C_QC)
# abundance (green): fan out, run, merge
xa = 570
line([(KX, Y0), (xa, Y0), (xa + Y0 - YA, YA), (916, YA),
      (916 + Y0 - YA, Y0), (M0, Y0)], C_AB)
# contamination (orange)
xc = 610
line([(KX, Y0), (xc, Y0), (xc + Y0 - YC, YC), (985, YC),
      (985 + Y0 - YC, Y0), (M0, Y0)], C_CO)
# functional (purple): straight lane with three spurs
FX0, FX1 = 860, 1040                  # fork / merge start
SP = 44
line([(KX, Y0), (FX0, Y0)], C_FN)
line([(FX0, Y0), (FX0 + SP, Y0 - SP), (FX1, Y0 - SP), (FX1 + SP, Y0)], C_FN)
line([(FX0, Y0), (FX0 + SP, Y0 + SP), (FX1, Y0 + SP), (FX1 + SP, Y0)], C_FN)
line([(FX0, Y0), (FX1 + SP, Y0), (M0, Y0)], C_FN)
# MultiQC spur (neutral grey): leaves the merged bundle after every line has joined
MQS = 1124
line([(MQS, Y0), (MQS + MQD, Y0 + MQD), (MQX, Y0 + MQD)], C_MQ)
# downstream (gold)
line([(M0, Y0), (RX, Y0)], C_DS)
# reconstruction (grey, dashed)
xr = 570
line([(KX, Y0), (xr, Y0), (xr + YR - Y0, YR), (1060, YR)], GREY, dashed=True)

# ---------------------------------------------------------------- trunk stations
qc = [("cat", 116, "", 1), ("Porechop_ABI", 178, "opt", -1), ("nanoq", 240, "", 1),
      ("NanoPlot", 302, "", -1), ("SeqKit stats", 364, "", 1), ("Hostile", 426, "", -1)]
for name, x, kind, side in qc:
    station(x, Y0, kind, fill=C_QC if name in ("nanoq", "NanoPlot", "SeqKit stats", "Hostile") else "#ffffff")
    ly = Y0 - 14 if side > 0 else Y0 + 22
    text(x, ly, [name], size=12, weight="700")
text(426, Y0 + 35, ["LOW_DEPTH flag"], size=10, fill=MUTED, style="italic")
stage(100, Y0 - 44, "QC &amp; HOST REMOVAL")

# MultiQC terminus (second output, grey)
terminus(MQX, MQY, C_MQ)
text(MQX, MQY + 28, ["MultiQC"], size=12, weight="700")
text(MQX, MQY + 42, ["report of all steps"], size=10, fill=MUTED, style="italic")

# input terminus
terminus(TX, Y0)
text(24, Y0 + 28, ["samplesheet +", "ONT FASTQ", "(multi-run pooling)"],
     anchor="start", size=11, fill=MUTED)

# Kraken2 interchange
interchange(KX, Y0, fill="#bfdcec")
text(KX, Y0 - 20, ["Kraken2 (PlusPF)"], size=12, weight="700")
stage(KX, Y0 - 44, "TAXONOMY", anchor="middle")

# ---------------------------------------------------------------- abundance
for name, x in [("Bracken", 750), ("KrakenTools|combine", 820), ("Krona", 890)]:
    station(x, YA, fill=C_AB if name == "Bracken" else "#ffffff")
    text(x, YA - 14 - (13 if "|" in name else 0), name.split("|"), size=12, weight="700")
stage(722, YA + 20, "ABUNDANCE")

# ---------------------------------------------------------------- contamination
station(770, YC, fill=C_CO)
text(770, YC - 14, ["Recentrifuge"], size=12, weight="700")
text(826, YC - 12, ["negative-control aware"], anchor="start", size=10, fill=MUTED,
     style="italic")
stage(725, YC + 20, "CONTAMINATION")

# ---------------------------------------------------------------- functional
for name, x in [("Tag reads", 690), ("ABRicate", 775)]:
    station(x, Y0, fill=C_FN if name == "ABRicate" else "#ffffff")
    text(x, Y0 - 14, [name], size=12, weight="700")
SX = (FX0 + SP + FX1) // 2
for name, y, dy in [("CARD", Y0 - SP, -12), ("VFDB", Y0, -12), ("PlasmidFinder", Y0 + SP, 20)]:
    station(SX, y, r=5.5)
    text(SX, y + dy, [name], size=12, weight="700")
stage(660, Y0 + 28, "FUNCTIONAL")

# ---------------------------------------------------------------- downstream
interchange(MX, Y0)
text(MX, Y0 - 20, ["MERIDIAN"], size=12, weight="700")
terminus(RX, Y0, C_DS)
text(RX, Y0 - 20, ["Report"], size=12, weight="700")
stage(MX, Y0 - 44, "DOWNSTREAM", anchor="middle")

# ---------------------------------------------------------------- reconstruction (planned)
for name, x in [("binning /|assembly", 770), ("CheckM2", 850), ("Bakta", 920),
                ("MLST", 985)]:
    station(x, YR, "plan")
    parts = name.split("|")
    text(x, YR - 14 - 13 * (len(parts) - 1), parts, size=12, weight="400", fill=MUTED)
terminus(1060, YR, GREY)
text(1060, YR - 20, ["Pathogen report"], size=12, weight="400", fill=MUTED)
stage(725, YR + 28, "METAGENOME RECONSTRUCTION (PLANNED)")

# ---------------------------------------------------------------- legend
LY = 424
stage(30, LY - 18, "LEGEND")
lines_leg = [("QC &amp; taxonomy", C_QC, False), ("Abundance", C_AB, False),
             ("Contamination", C_CO, False), ("Functional (AMR / VF / MGE)", C_FN, False),
             ("Downstream", C_DS, False), ("Planned", GREY, True)]
xs = [30, 190, 320, 470, 700, 840]
for (lab, col, dashed), x in zip(lines_leg, xs):
    dash = ' stroke-dasharray="8 6"' if dashed else ""
    add(f'<line x1="{x}" y1="{LY}" x2="{x + 34}" y2="{LY}" stroke="{col}" '
        f'stroke-width="{LW}" stroke-linecap="round"{dash}/>')
    text(x + 44, LY + 4, [lab], anchor="start", size=12)
SY = LY + 32
sym = [("Process", "", 30), ("Optional step", "opt", 150), ("Planned step", "plan", 300)]
for lab, kind, x in sym:
    station(x + 8, SY, kind)
    text(x + 24, SY + 4, [lab], anchor="start", size=12)
interchange(470 + 24, SY, w=44)
text(470 + 56, SY + 4, ["Interchange"], anchor="start", size=12)
station(640, SY, fill=C_QC)
text(654, SY + 4, ["Reported in MultiQC"], anchor="start", size=12)
text(830, SY + 4, ["Databases: --build_databases to download, or --db_dir to reuse"],
     anchor="start", size=10.5, fill=MUTED, style="italic")

add("</svg>")

OUT.write_text("\n".join(out) + "\n", encoding="utf-8", newline="\n")
print(f"wrote {OUT}")
