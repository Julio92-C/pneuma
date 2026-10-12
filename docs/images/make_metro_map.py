#!/usr/bin/env python3
"""Generate docs/images/pneuma_metro_map.svg (nf-core "tube map" style).

Eight labelled stage blocks (two rows of four), one coloured line per block.
Standard library only; output is deterministic, so re-running must leave
the committed SVG unchanged:

    python docs/images/make_metro_map.py
"""
from pathlib import Path

OUT = Path(__file__).with_name("pneuma_metro_map.svg")

BG = "#f6f8fa"       # light backing panel (readable on light and dark GitHub themes)
BOX = "#ffffff"      # block fill
INK = "#1f2328"
MUTED = "#57606a"
GREY = "#8c959f"
FONT = "Helvetica, Arial, sans-serif"

BW, GAP, MX = 270, 30, 30          # block width, column gap, outer margin
W = MX * 2 + 4 * BW + 3 * GAP      # 1230
TOP1 = 76                          # top of row 1
ROW_GAP = 46                       # vertical gap between the two rows

# key, colour, title, subtitle, planned, items
# item = (name, descriptor, kind)   kind: "" station | "opt" optional | "sub" sub-station
BLOCKS = [
    ("inp", "#56B4E9", "Input", "samplesheet + ONT FASTQ", False, [
        ("samplesheet.csv", "one row per sample x run", ""),
        ("ONT FASTQ", "multi-run / multi-batch pooling", ""),
    ]),
    ("qc", "#0072B2", "QC &amp; host removal", "low-biomass safeguards", False, [
        ("cat", "pool runs per sample", ""),
        ("Porechop_ABI", "adapter trimming (optional)", "opt"),
        ("nanoq", "length / quality filter", ""),
        ("NanoPlot", "", ""),
        ("SeqKit stats", "", ""),
        ("Hostile", "host read removal", ""),
        ("LOW_DEPTH flag", "fewer than --min_reads reads", ""),
    ]),
    ("tax", "#009E73", "Taxonomic classification", "", False, [
        ("Kraken2", "PlusPF index", ""),
    ]),
    ("abu", "#7A9A01", "Abundance estimation", "", False, [
        ("Bracken", "", ""),
        ("KrakenTools combine", "", ""),
        ("Krona", "", ""),
    ]),
    ("con", "#D55E00", "Contamination assessment", "", False, [
        ("Recentrifuge", "negative-control subtraction", ""),
    ]),
    ("fun", "#9B3FA8", "Functional profiling", "ABRicate on taxon-tagged reads", False, [
        ("Tag reads", "by taxid", ""),
        ("ABRicate", "80 % identity / 80 % coverage", ""),
        ("AMR", "CARD", "sub"),
        ("Virulence", "VFDB", "sub"),
        ("MGE", "PlasmidFinder", "sub"),
    ]),
    ("mer", "#B07D00", "Downstream analysis", "MERIDIAN", False, [
        ("Diversity", "", ""),
        ("Differential abundance", "", ""),
        ("Resistome / virulome / mobilome", "", ""),
        ("Report", "", ""),
    ]),
    ("rec", GREY, "Metagenome reconstruction", "planned, Phase 2", True, [
        ("Priority-pathogen binning", "and assembly", ""),
        ("CheckM2", "", ""),
        ("Bakta", "", ""),
        ("MLST", "", ""),
        ("Pathogen report", "", ""),
    ]),
]

out = []


def add(s):
    out.append(s)


def rpath(pts, r=22):
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


def text(x, y, lines, anchor="middle", size=12.5, weight="400", fill=INK,
         halo=BOX, style=""):
    tsp = "".join(f'<tspan x="{x}" dy="{0 if i == 0 else 15}">{t}</tspan>'
                  for i, t in enumerate(lines))
    add(f'<text x="{x}" y="{y}" text-anchor="{anchor}" font-size="{size}" '
        f'font-weight="{weight}" fill="{fill}" stroke="{halo}" stroke-width="3.5" '
        f'paint-order="stroke" stroke-linejoin="round"{style}>{tsp}</text>')


def station(x, y, dashed=False, stroke=INK, r=7.5):
    dash = ' stroke-dasharray="3 2"' if dashed else ""
    add(f'<circle cx="{x}" cy="{y}" r="{r}" fill="#ffffff" stroke="{stroke}" '
        f'stroke-width="2.5"{dash}/>')


def arrow(pts, dashed=False):
    dash = ' stroke-dasharray="7 5"' if dashed else ""
    add(f'<path d="{rpath(pts, 14)}" fill="none" stroke="{MUTED}" stroke-width="3" '
        f'stroke-linejoin="round"{dash}/>')
    (x0, y0), (x1, y1) = pts[-2], pts[-1]
    if abs(x1 - x0) >= abs(y1 - y0):
        s = 1 if x1 > x0 else -1
        tri = f"{x1} {y1} {x1 - 10 * s} {y1 - 6} {x1 - 10 * s} {y1 + 6}"
    else:
        s = 1 if y1 > y0 else -1
        tri = f"{x1} {y1} {x1 - 6} {y1 - 10 * s} {x1 + 6} {y1 - 10 * s}"
    add(f'<polygon points="{tri}" fill="{MUTED}"/>')


def item_height(it):
    return 40 if it[1] and it[2] != "sub" else 30


def block_height(items):
    return 62 + sum(item_height(i) for i in items) + 10


H1 = max(block_height(b[5]) for b in BLOCKS[:4])
H2 = max(block_height(b[5]) for b in BLOCKS[4:])
TOP2 = TOP1 + H1 + ROW_GAP
LEG_Y = TOP2 + H2 + 24
LEG_H = 122
H = LEG_Y + LEG_H + 20

# ---------------------------------------------------------------- header
add(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" '
    f'height="{H}" role="img" font-family="{FONT}">')
add("<title>PNEUMA pipeline metro map</title>")
add("<desc>Metro-map overview of the PNEUMA pipeline as eight labelled stage "
    "blocks: Input, QC and host removal, Taxonomic classification, Abundance "
    "estimation, Contamination assessment, Functional profiling, Downstream "
    "analysis (MERIDIAN) and the planned Phase 2 Metagenome reconstruction.</desc>")
add(f'<rect x="1" y="1" width="{W - 2}" height="{H - 2}" rx="18" fill="{BG}" '
    f'stroke="#d0d7de" stroke-width="2"/>')
text(MX, 40, ["PNEUMA: ONT metagenomics for low-biomass respiratory samples"],
     anchor="start", size=18, weight="700", halo=BG)

# ---------------------------------------------------------------- blocks
geom = []
for n, (key, col, title, sub, planned, items) in enumerate(BLOCKS):
    row, c = divmod(n, 4)
    bx = MX + c * (BW + GAP)
    by = TOP1 if row == 0 else TOP2
    bh = H1 if row == 0 else H2
    geom.append((bx, by, bh))
    dash = ' stroke-dasharray="9 6"' if planned else ""
    add(f'<rect x="{bx}" y="{by}" width="{BW}" height="{bh}" rx="14" fill="{BOX}" '
        f'stroke="{col}" stroke-width="2.5"{dash}/>')
    # numbered header
    add(f'<circle cx="{bx + 24}" cy="{by + 26}" r="12" fill="{col}"/>')
    add(f'<text x="{bx + 24}" y="{by + 31}" text-anchor="middle" font-size="13" '
        f'font-weight="700" fill="#ffffff">{n + 1}</text>')
    text(bx + 44, by + 26, [title], anchor="start", size=14, weight="700")
    if sub:
        text(bx + 44, by + 43, [sub], anchor="start", size=11, fill=MUTED)
    add(f'<line x1="{bx + 12}" y1="{by + 54}" x2="{bx + BW - 12}" y2="{by + 54}" '
        f'stroke="{col}" stroke-width="1.5" stroke-opacity="0.5"/>')
    # station track
    lx = bx + 28
    ys, y = [], by + 62
    for it in items:
        ys.append(y + 15)
        y += item_height(it)
    if len(ys) > 1:
        dl = ' stroke-dasharray="12 8"' if planned else ""
        add(f'<line x1="{lx}" y1="{ys[0]}" x2="{lx}" y2="{ys[-1]}" stroke="{col}" '
            f'stroke-width="6" stroke-linecap="round"{dl}/>')
    for (name, desc, kind), sy in zip(items, ys):
        sub_st = kind == "sub"
        station(lx, sy, dashed=(kind == "opt" or planned),
                stroke=MUTED if planned else INK, r=5.5 if sub_st else 7.5)
        if sub_st:
            add(f'<text x="{lx + 22}" y="{sy + 4}" font-size="12.5" fill="{MUTED}" '
                f'stroke="{BOX}" stroke-width="3.5" paint-order="stroke" '
                f'stroke-linejoin="round">{name}: <tspan font-weight="700" '
                f'fill="{INK}">{desc}</tspan></text>')
        else:
            text(lx + 22, sy + 4, [name], anchor="start", size=12.5,
                 weight="400" if planned else "700",
                 fill=MUTED if planned else INK)
            if desc:
                text(lx + 22, sy + 19, [desc], anchor="start", size=11, fill=MUTED)

# ---------------------------------------------------------------- flow arrows
for n in range(7):
    bx, by, bh = geom[n]
    nx, ny, nbh = geom[n + 1]
    planned = BLOCKS[n + 1][4]
    if n == 3:  # wrap: row 1 -> row 2
        x0, x1 = bx + BW // 2, nx + BW // 2
        mid = TOP1 + H1 + ROW_GAP // 2
        arrow([(x0, by + bh), (x0, mid), (x1, mid), (x1, ny)])
    else:
        arrow([(bx + BW, by + 26), (nx, ny + 26)], dashed=planned)

# ---------------------------------------------------------------- legend
add(f'<rect x="{MX}" y="{LEG_Y}" width="{W - 2 * MX}" height="{LEG_H}" rx="12" '
    f'fill="{BOX}" fill-opacity="0.6" stroke="#d0d7de" stroke-width="1.5"/>')
text(MX + 16, LEG_Y + 24, ["LEGEND"], anchor="start", size=12, weight="700", halo=BG)
colw = (W - 2 * MX - 32) / 4
for n, (key, col, title, sub, planned, _) in enumerate(BLOCKS):
    cx = MX + 16 + (n % 4) * colw
    cy = LEG_Y + 46 + (n // 4) * 24
    dash = ' stroke-dasharray="8 5"' if planned else ""
    add(f'<line x1="{cx:.0f}" y1="{cy}" x2="{cx + 40:.0f}" y2="{cy}" stroke="{col}" '
        f'stroke-width="6"{dash}/>')
    lab = f"{n + 1}  {title}" + (" (planned)" if planned else "")
    text(cx + 52, cy + 4, [lab], anchor="start", size=12.5, halo=BG)
sy = LEG_Y + 104
sx = MX + 16
station(sx + 8, sy)
text(sx + 26, sy + 4, ["Process / tool"], anchor="start", size=12.5, halo=BG)
station(sx + colw + 8, sy, dashed=True)
text(sx + colw + 26, sy + 4, ["Optional step"], anchor="start", size=12.5, halo=BG)
station(sx + 2 * colw + 8, sy, dashed=True, stroke=MUTED)
text(sx + 2 * colw + 26, sy + 4, ["Planned (Phase 2)"], anchor="start", size=12.5, halo=BG)
arrow([(sx + 3 * colw - 2, sy), (sx + 3 * colw + 32, sy)])
text(sx + 3 * colw + 46, sy + 4, ["Data flow"], anchor="start", size=12.5, halo=BG)

add("</svg>")

OUT.write_text("\n".join(out) + "\n", encoding="utf-8", newline="\n")
print(f"wrote {OUT}")
