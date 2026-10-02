# makes the research charts in index.html from the genie3_brca_project results (run: python python/make_figures.py)

import csv
import json
import math
import re
import sys
from pathlib import Path

SITE = Path(__file__).resolve().parent.parent
INDEX = SITE / "index.html"
SUMMARY = SITE / "data" / "research-summary.json"
FREE_THROWS = SITE / "data" / "free-throws.csv"
RESULTS = SITE.parent / "genie3_brca_project" / "results"

SUBTYPES = ["Basal", "LumA", "LumB"]
NAMES = {"Basal": "Basal-like", "LumA": "Luminal A", "LumB": "Luminal B"}
MODULES = ["LumA_M1", "LumB_M1", "Basal_M1"]


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def cox(row):
    return {
        "hr": round(float(row["hr_uni"]), 3),
        "low": round(float(row["ci_lower_uni"]), 3),
        "high": round(float(row["ci_upper_uni"]), 3),
        "p": round(float(row["p_uni"]), 4),
        "patients": int(row["n"]),
        "deaths": int(row["events"]),
    }


def load_results(results):
    data = {"overlap": [], "hubs": {}}

    for row in read_tsv(results / "rewiring" / "pairwise_overlap.tsv"):
        a, b = row["group_a"], row["group_b"]
        sig = read_tsv(results / "rewiring" / f"{a}_vs_{b}_significant_edges.tsv")
        data["overlap"].append({"a": a, "b": b, "top": int(row["n_a"]), "shared": int(row["shared"]), "differ": len(sig)})

    for s in SUBTYPES:
        rows = read_tsv(results / "hubs" / f"{s}_hub_frequency.tsv")
        rows.sort(key=lambda r: -float(r["bootstrap_frequency"]))
        data["hubs"][s] = []
        for r in rows[:9]:
            out_of = int(r["n_bootstrap"])
            times = round(float(r["bootstrap_frequency"]) * out_of)
            data["hubs"][s].append({"gene": r["gene"], "times": times, "out_of": out_of, "tf": r["is_tf"] == "True"})

    tcga = {r["module"]: r for r in read_tsv(results / "survival" / "tcga_module_cox.tsv")}
    metabric = {r["module"]: r for r in read_tsv(results / "validation" / "metabric_module_cox.tsv")}
    sizes = {r["subtype"] + "_" + r["module"]: int(r["size"]) for r in read_tsv(results / "modules" / "module_coherence.tsv")}
    data["survival"] = []
    for m in MODULES:
        data["survival"].append({
            "module": m,
            "subtype": m.split("_")[0],
            "genes": sizes[m],
            "tcga": cox(tcga[m]),
            "metabric": cox(metabric[m]),
        })

    edges = read_tsv(results / "networks" / "Basal_top_edges.tsv")
    strength = {}
    for e in edges:
        if "SOX10" in (e["regulator"], e["target"]):
            other = e["target"] if e["regulator"] == "SOX10" else e["regulator"]
            strength[other] = max(strength.get(other, 0), float(e["importance"]))
    partners = sorted(strength, key=strength.get, reverse=True)[:12]

    links = {}
    for e in edges:
        if e["regulator"] in partners and e["target"] in partners:
            pair = tuple(sorted([e["regulator"], e["target"]]))
            links[pair] = max(links.get(pair, 0), float(e["importance"]))

    data["network"] = {
        "hub": "SOX10",
        "partners": [{"gene": g, "weight": round(strength[g], 3)} for g in partners],
        "links": [{"a": a, "b": b, "weight": round(w, 3)} for (a, b), w in sorted(links.items())],
    }
    return data


def network_svg(net):
    cx, cy, rx, ry = 200, 165, 122, 112
    partners = net["partners"]
    biggest = max(p["weight"] for p in partners)

    spots = {}
    for i, p in enumerate(partners):
        angle = -math.pi / 2 + 2 * math.pi * i / len(partners)
        spots[p["gene"]] = (angle, cx + rx * math.cos(angle), cy + ry * math.sin(angle))

    out = [
        '<svg class="net" viewBox="0 0 400 330" role="img" aria-labelledby="net-title">',
        f'  <title id="net-title">{net["hub"]} and its {len(partners)} strongest partner genes in the Basal-like network</title>',
    ]
    for link in net["links"]:
        _, x1, y1 = spots[link["a"]]
        _, x2, y2 = spots[link["b"]]
        qx = 0.3 * (x1 + x2) + 0.4 * cx
        qy = 0.3 * (y1 + y2) + 0.4 * cy
        out.append(f'  <path class="net-link" fill="none" d="M{x1:.1f} {y1:.1f} Q{qx:.1f} {qy:.1f} {x2:.1f} {y2:.1f}"/>')

    for p in partners:
        _, x, y = spots[p["gene"]]
        width = 0.8 + 3.2 * p["weight"] / biggest
        out.append(f'  <line class="net-spoke" x1="{cx}" y1="{cy}" x2="{x:.1f}" y2="{y:.1f}" stroke-width="{width:.2f}"/>')

    out.append(f'  <circle class="net-halo" cx="{cx}" cy="{cy}" r="34"/>')

    for p in partners:
        angle, x, y = spots[p["gene"]]
        r = 3.5 + 3.5 * p["weight"] / biggest
        if math.cos(angle) > 0.35:
            anchor, lx, ly = "start", x + r + 6, y + 4
        elif math.cos(angle) < -0.35:
            anchor, lx, ly = "end", x - r - 6, y + 4
        elif math.sin(angle) < 0:
            anchor, lx, ly = "middle", x, y - r - 7
        else:
            anchor, lx, ly = "middle", x, y + r + 15
        out.append(f'  <circle class="net-node" cx="{x:.1f}" cy="{y:.1f}" r="{r:.1f}"/>')
        out.append(f'  <text class="net-label" x="{lx:.1f}" y="{ly:.1f}" text-anchor="{anchor}">{p["gene"]}</text>')

    out.append(f'  <circle class="net-hub" cx="{cx}" cy="{cy}" r="24"/>')
    out.append(f'  <text class="net-hub-label" x="{cx}" y="{cy + 4}" text-anchor="middle">{net["hub"]}</text>')
    out.append('</svg>')
    return out


def overlap_table(overlap):
    out = [
        '<div class="table-scroll">',
        '  <table class="overlap-table">',
        '    <thead>',
        '      <tr><th scope="col">Networks compared</th><th scope="col">Shared</th><th scope="col">Really differ</th></tr>',
        '    </thead>',
        '    <tbody>',
    ]
    for row in overlap:
        share = 100 * row["shared"] / row["top"]
        zero = " is-zero" if row["differ"] == 0 else ""
        out += [
            '      <tr>',
            f'        <th scope="row">{NAMES[row["a"]]} <span class="vs">vs</span> {NAMES[row["b"]]}'
            f'<svg class="meter" aria-hidden="true"><rect class="meter-track" width="100%" height="100%" rx="4"/>'
            f'<rect class="meter-fill" width="{share:.1f}%" height="100%" rx="4"/></svg></th>',
            f'        <td class="num">{row["shared"]:,}<span class="of"> / {row["top"]:,}</span></td>',
            f'        <td class="num num-big{zero}">{row["differ"]}</td>',
            '      </tr>',
        ]
    out += ['    </tbody>', '  </table>', '</div>']
    return out


def hub_bars(hubs):
    out = ['<div class="hubs">']
    for s in SUBTYPES:
        out.append('  <div>')
        out.append(f'    <h4 class="hub-title">{NAMES[s]}</h4>')
        out.append('    <ol class="hub-list">')
        for h in hubs[s]:
            share = 100 * h["times"] / h["out_of"]
            tf = ' class="is-tf"' if h["tf"] else ""
            note = '<span class="visually-hidden"> (transcription factor)</span>' if h["tf"] else ""
            out.append(
                f'      <li{tf}><span class="gene">{h["gene"]}</span>{note}'
                f'<svg class="hub-bar" aria-hidden="true"><rect class="bar-track" width="100%" height="100%" rx="3"/>'
                f'<rect class="bar-fill" width="{share:.1f}%" height="100%" rx="3"/>'
                f'<line class="bar-cut" x1="70%" y1="-25%" x2="70%" y2="125%"/></svg>'
                f'<span class="hub-val">{h["times"]}/{h["out_of"]}</span></li>'
            )
        out.append('    </ol>')
        out.append('  </div>')
    out.append('</div>')
    return out


def forest_plot(survival):
    lo, hi = 0.6, 1.1
    ticks = [0.6, 0.7, 0.8, 0.9, 1.0, 1.1]

    def x(v):
        return 100 * (math.log(v) - math.log(lo)) / (math.log(hi) - math.log(lo))

    def interval(cls, c, y):
        a, b, m = x(c["low"]), x(c["high"]), x(c["hr"])
        return (
            f'<line class="ci {cls}" x1="{a:.1f}%" y1="{y}%" x2="{b:.1f}%" y2="{y}%"/>'
            f'<line class="ci {cls}" x1="{a:.1f}%" y1="{y - 10}%" x2="{a:.1f}%" y2="{y + 10}%"/>'
            f'<line class="ci {cls}" x1="{b:.1f}%" y1="{y - 10}%" x2="{b:.1f}%" y2="{y + 10}%"/>'
            f'<circle class="dot {cls}" cx="{m:.1f}%" cy="{y}%" r="6.5"/>'
        )

    grid = ""
    for t in ticks:
        one = " is-one" if t == 1.0 else ""
        grid += f'<line class="forest-grid{one}" x1="{x(t):.1f}%" y1="0" x2="{x(t):.1f}%" y2="100%"/>'

    out = ['<div class="forest">']
    for m in survival:
        t, v = m["tcga"], m["metabric"]
        out += [
            '  <div class="forest-row">',
            f'    <p class="forest-name">{NAMES[m["subtype"]]} module<span>{m["genes"]} genes</span></p>',
            f'    <svg class="forest-plot" aria-hidden="true">{grid}{interval("is-tcga", t, 34)}{interval("is-metabric", v, 68)}</svg>',
            f'    <p class="forest-vals" aria-hidden="true"><span class="is-tcga">{t["hr"]:.2f}</span> vs <span class="is-metabric">{v["hr"]:.2f}</span></p>',
            f'    <p class="visually-hidden">Hazard ratio {t["hr"]:.2f} in TCGA ({t["low"]:.2f} to {t["high"]:.2f}) and {v["hr"]:.2f} in METABRIC ({v["low"]:.2f} to {v["high"]:.2f}).</p>',
            '  </div>',
        ]

    labels = "".join(f'<text x="{x(t):.1f}%" y="14" text-anchor="middle">{t:.1f}</text>' for t in ticks)
    out += [
        '  <div class="forest-row" aria-hidden="true">',
        f'    <svg class="forest-axis">{labels}</svg>',
        '  </div>',
        '</div>',
    ]
    return out


def sparkline(weeks):
    pct = [100 * wk["made"] / wk["attempted"] for wk in weeks]
    lo, hi = 70, 85
    width, height, pad = 300, 100, 8

    def y(p):
        return pad + (hi - p) * (height - 2 * pad) / (hi - lo)

    xs = [pad + i * (width - 2 * pad) / (len(pct) - 1) for i in range(len(pct))]
    ys = [y(p) for p in pct]
    points = " ".join(f"{a:.1f},{b:.1f}" for a, b in zip(xs, ys))
    area = f"M{xs[0]:.1f},{height - pad} L" + " L".join(points.split(" ")) + f" L{xs[-1]:.1f},{height - pad} Z"
    lines = "".join(f'<line class="spark-grid" x1="{pad}" y1="{y(v):.1f}" x2="{width - pad}" y2="{y(v):.1f}"/>' for v in range(lo, hi + 1, 5))
    made = sum(wk["made"] for wk in weeks)
    tried = sum(wk["attempted"] for wk in weeks)

    return [
        f'<svg class="spark-svg" viewBox="0 0 {width} {height}" role="img" aria-label="Free throw percentage by week, from {pct[0]:.0f}% to {pct[-1]:.0f}%">',
        f'  {lines}',
        f'  <path class="spark-area" d="{area}"/>',
        f'  <polyline class="spark-line" fill="none" points="{points}"/>',
        f'  <circle class="spark-dot" cx="{xs[-1]:.1f}" cy="{ys[-1]:.1f}" r="5"/>',
        '</svg>',
        f'<figcaption class="spark-caption"><span>Week 1: {pct[0]:.0f}%</span><span>Week {len(weeks)}: {pct[-1]:.0f}%</span>'
        f'<span class="spark-total">{made:,} of {tried:,} made overall</span></figcaption>',
    ]


def put_in_page(html, name, lines):
    m = re.search(rf"([ \t]*)(<!-- figure:{name} -->).*?(<!-- /figure:{name} -->)", html, re.S)
    if not m:
        sys.exit(f"couldn't find <!-- figure:{name} --> in index.html")
    pad = m.group(1)
    new = pad + m.group(2) + "\n" + "\n".join(pad + line for line in lines) + "\n" + pad + m.group(3)
    return html[:m.start()] + new + html[m.end():]


results = Path(sys.argv[1]) if len(sys.argv) > 1 else RESULTS
if results.exists():
    data = load_results(results)
    SUMMARY.write_text(json.dumps(data, indent=2), encoding="utf-8")
    print("read results from", results)
else:
    data = json.loads(SUMMARY.read_text(encoding="utf-8"))
    print("no results folder, using", SUMMARY.name)

with open(FREE_THROWS, encoding="utf-8") as f:
    weeks = [{"made": int(r["made"]), "attempted": int(r["attempted"])} for r in csv.DictReader(f)]

html = INDEX.read_text(encoding="utf-8")
html = put_in_page(html, "network", network_svg(data["network"]))
html = put_in_page(html, "overlap", overlap_table(data["overlap"]))
html = put_in_page(html, "hubs", hub_bars(data["hubs"]))
html = put_in_page(html, "survival", forest_plot(data["survival"]))
html = put_in_page(html, "freethrows", sparkline(weeks))

with open(INDEX, "w", encoding="utf-8", newline="\n") as f:
    f.write(html)
print("updated index.html")
