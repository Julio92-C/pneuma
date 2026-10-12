#!/usr/bin/env python3
"""Generate docs/images/pneuma_metro_map.svg (nf-core "tube map" style).

Standard library only; output is deterministic, so re-running must leave
the committed SVG unchanged:

    python docs/images/make_metro_map.py
"""
from pathlib import Path

OUT = Path(__file__).with_name("pneuma_metro_map.svg")

W, H = 1520, 790
BG = "#f6f8fa"       # light backing panel (readable on light and dark GitHub themes)
INK = "#1f2328"
MUTED = "#57606a"
GREY = "#8c959f"
FONT = "Helvetica, Arial, sans-serif"

COL = {
    "pre": "#0072B2",   # Preprocess
    "tax": "#009E73",   # Taxonomy
    "con": "#D55E00",   # Contamination
    "fun": "#9B3FA8",   # Functional
    "mer": "#B07D00",   # MERIDIAN
    "pat": GREY,        # Phase 2 pathogen (dashed)
}
NAME = {"pre": "PREPROCESS", "tax": "TAXONOMY", "con": "CONTAMINATION",
        "fun": "FUNCTIONAL", "mer": "MERIDIAN", "pat": "PATHOGEN (PHASE 2)"}

Y = {"pre": 110, "tax": 220, "con": 330, "fun": 440, "pat": 560, "mer": 680}
END = 1330  # x of line termini

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


def line(key, pts, dashed=False, branch=False):
    d = rpath(pts)
    # halo underlay; for branches it starts 14 px below the parent line so the
    # parent is not cut where the branch leaves it
    hp = [(pts[0][0], pts[0][1] + 14)] + pts[1:] if branch else pts
    add(f'<path d="{rpath(hp)}" fill="none" stroke="{BG}" stroke-width="13" stroke-linejoin="round"/>')
    dash = ' stroke-dasharray="12 8"' if dashed else ""
    add(f'<path d="{d}" fill="none" stroke="{COL[key]}" stroke-width="6" '
        f'stroke-linejoin="round"{dash}/>')
    x, y = pts[-1]
    add(f'<rect x="{x}" y="{y - 11}" width="5" height="22" rx="1.5" fill="{COL[key]}"/>')


def text(x, y, lines, anchor="middle", size=12.5, weight="400", fill=INK, style=""):
    """Multi-line text with a light halo; y is the first baseline."""
    tsp = "".join(f'<tspan x="{x}" dy="{0 if i == 0 else 15}">{t}</tspan>'
                  for i, t in enumerate(lines))
    add(f'<text x="{x}" y="{y}" text-anchor="{anchor}" font-size="{size}" '
        f'font-weight="{weight}" fill="{fill}" stroke="{BG}" stroke-width="3.5" '
        f'paint-order="stroke" stroke-linejoin="round"{style}>{tsp}</text>')


def station(x, y, dashed=False, stroke=INK):
    dash = ' stroke-dasharray="3 2"' if dashed else ""
    add(f'<circle cx="{x}" cy="{y}" r="7.5" fill="#ffffff" stroke="{stroke}" '
        f'stroke-width="2.5"{dash}/>')


def label(x, y, lines, pos="above"):
    if pos == "above":
        text(x, y - 15 - 15 * (len(lines) - 1), lines)
    else:
        text(x, y + 26, lines)


def pill(x, y, key, txt=None):
    t = txt or NAME[key]
    w = 7.4 * len(t) + 16
    add(f'<rect x="{x}" y="{y}" width="{w:.0f}" height="20" rx="10" fill="{COL[key]}"/>')
    add(f'<text x="{x + w / 2:.1f}" y="{y + 14}" text-anchor="middle" font-size="11" '
        f'font-weight="700" fill="#ffffff" letter-spacing="0.5">{t}</text>')


# ---------------------------------------------------------------- header
add(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" '
    f'height="{H}" role="img" font-family="{FONT}">')
add("<title>PNEUMA pipeline metro map</title>")
add("<desc>Metro-map overview of the PNEUMA pipeline: Preprocess, Taxonomy, "
    "Contamination, Functional and MERIDIAN lines (Phase 1) and the planned "
    "Phase 2 pathogen line.</desc>")
add(f'<rect x="1" y="1" width="{W - 2}" height="{H - 2}" rx="18" fill="{BG}" '
    f'stroke="#d0d7de" stroke-width="2"/>')
text(40, 38, ["PNEUMA: ONT respiratory metagenomics"], anchor="start",
     size=18, weight="700")

# ---------------------------------------------------------------- input box
add(f'<rect x="30" y="66" width="140" height="88" rx="12" fill="#ffffff" '
    f'stroke="{INK}" stroke-width="2"/>')
add(f'<text x="100" y="94" text-anchor="middle" font-size="13" font-weight="700" fill="{INK}">INPUT</text>')
add(f'<text x="100" y="114" text-anchor="middle" font-size="12.5" fill="{INK}">FASTQ (ONT)</text>')
add(f'<text x="100" y="132" text-anchor="middle" font-size="12.5" fill="{INK}">+ samplesheet.csv</text>')

# ---------------------------------------------------------------- lines
PRE_X = [245, 325, 405, 485, 565, 635, 705]
BR = 705  # taxonomy branches off after LOW_DEPTH
line("pre", [(170, Y["pre"]), (END, Y["pre"])])
line("tax", [(BR, Y["pre"]), (BR, Y["tax"]), (END, Y["tax"])], branch=True)
line("con", [(880, Y["tax"]), (880, Y["con"]), (END, Y["con"])], branch=True)
line("fun", [(850, Y["tax"]), (850, Y["fun"]), (END, Y["fun"])], branch=True)
line("pat", [(820, Y["tax"]), (820, Y["pat"]), (END, Y["pat"])], dashed=True, branch=True)
line("mer", [(1100, Y["tax"]), (1100, Y["mer"]), (END, Y["mer"])], branch=True)

# ---------------------------------------------------------------- stations
pre = [("cat", ["cat"], False), ("porechop_abi", ["porechop-ABI", "(optional)"], True),
       ("nanoq", ["nanoq"], False), ("NanoPlot", ["NanoPlot"], False),
       ("seqkit", ["seqkit stats"], False), ("Hostile", ["Hostile"], False),
       ("LOW_DEPTH", ["LOW_DEPTH", "flag"], False)]
for x, (_, lab, opt) in zip(PRE_X, pre):
    station(x, Y["pre"], dashed=opt)
    label(x, Y["pre"], lab, "above")

tax = [(780, ["Kraken2", "(PlusPF)"]), (940, ["Bracken"]),
       (1100, ["KrakenTools", "combine"]), (1250, ["Krona"])]
for x, lab in tax:
    station(x, Y["tax"])
    # the MERIDIAN branch leaves from KrakenTools combine, so keep its label clear of it
    label(x, Y["tax"], lab, "above" if x == 1100 else "below")

for x, lab in [(945, ["Recentrifuge", "(plain)"]), (1040, ["Recentrifuge", "(control)"])]:
    station(x, Y["con"])
    label(x, Y["con"], lab, "below")

for x, lab in [(950, ["Tag reads", "(taxid)"]),
               (1040, ["ABRicate", "card, vfdb,", "plasmidfinder", "(80/80)"])]:
    station(x, Y["fun"])
    label(x, Y["fun"], lab, "below")

pat = [(890, ["Assembly"]), (960, ["CheckM2"]), (1030, ["Bakta"]),
       (1170, ["MLST"]), (1260, ["Pathogen report", "AMR / VF / MGE"])]
for x, lab in pat:
    station(x, Y["pat"], dashed=True, stroke=MUTED)
    text(x, Y["pat"] + 26, lab, fill=MUTED)

for x, lab in [(1180, ["MERIDIAN", "config"]), (1260, ["MERIDIAN run", "(container v1.0.0)"])]:
    station(x, Y["mer"])
    label(x, Y["mer"], lab, "below")

# interchange stations (branch points)
station(BR, Y["pre"])
station(1100, Y["tax"])

# ---------------------------------------------------------------- pills
pill(185, Y["pre"] + 12, "pre")
pill(BR + 16, Y["tax"] - 32, "tax")
pill(896, Y["con"] - 32, "con")
pill(866, Y["fun"] - 32, "fun")
pill(836, Y["pat"] - 32, "pat")
pill(1122, Y["mer"] - 32, "mer")
text(836, Y["pat"] + 80, ["WHO-priority respiratory panel: tiers 1, 2, 3, D",
                          "+ reservoir commensals list (planned)"],
     anchor="start", size=12, fill=MUTED, style=' font-style="italic"')

# ---------------------------------------------------------------- outputs
add(f'<rect x="1360" y="66" width="140" height="690" rx="12" fill="#ffffff" '
    f'stroke="{INK}" stroke-width="2"/>')
add(f'<text x="1430" y="94" text-anchor="middle" font-size="13" font-weight="700" fill="{INK}">OUTPUTS</text>')
for key, lab in [("pre", "preprocess/"), ("tax", "taxonomy/"), ("con", "contamination/"),
                 ("fun", "functional/"), ("pat", "pathogen/ (Ph. 2)"), ("mer", "meridian/")]:
    y = Y[key]
    add(f'<rect x="1372" y="{y - 6}" width="12" height="12" rx="3" fill="{COL[key]}"/>')
    add(f'<text x="1392" y="{y + 4}" font-size="11.5" fill="{INK}">{lab}</text>')
add(f'<text x="1430" y="732" text-anchor="middle" font-size="11.5" fill="{MUTED}">multiqc/ · pipeline_info/</text>')
add(f'<text x="1430" y="748" text-anchor="middle" font-size="11.5" fill="{MUTED}">HTML report</text>')

# ---------------------------------------------------------------- legend
LX, LY = 40, 560
add(f'<rect x="{LX - 10}" y="{LY - 10}" width="660" height="200" rx="12" fill="#ffffff" '
    f'fill-opacity="0.6" stroke="#d0d7de" stroke-width="1.5"/>')
text(LX + 4, LY + 14, ["LEGEND"], anchor="start", size=12, weight="700")
rows = [("pre", "Preprocess"), ("tax", "Taxonomy"), ("con", "Contamination"),
        ("fun", "Functional"), ("mer", "MERIDIAN"), ("pat", "Phase 2 pathogen (planned)")]
for i, (key, lab) in enumerate(rows):
    cx = LX + 4 + (i // 3) * 300
    cy = LY + 40 + (i % 3) * 28
    dash = ' stroke-dasharray="8 5"' if key == "pat" else ""
    add(f'<line x1="{cx}" y1="{cy}" x2="{cx + 44}" y2="{cy}" stroke="{COL[key]}" '
        f'stroke-width="6"{dash}/>')
    text(cx + 56, cy + 4, [lab], anchor="start", size=12.5)
sy = LY + 130
station(LX + 26, sy)
text(LX + 56, sy + 4, ["Process / tool"], anchor="start", size=12.5)
station(LX + 326, sy, dashed=True, stroke=INK)
text(LX + 356, sy + 4, ["Optional step"], anchor="start", size=12.5)
station(LX + 26, sy + 28)
add(f'<circle cx="{LX + 26}" cy="{sy + 28}" r="11.5" fill="none" stroke="{INK}" stroke-width="1.5"/>')
text(LX + 56, sy + 32, ["Interchange (lines branch)"], anchor="start", size=12.5)
station(LX + 326, sy + 28, dashed=True, stroke=MUTED)
text(LX + 356, sy + 32, ["Planned (Phase 2)"], anchor="start", size=12.5)

add("</svg>")

OUT.write_text("\n".join(out) + "\n", encoding="utf-8", newline="\n")
print(f"wrote {OUT}")
