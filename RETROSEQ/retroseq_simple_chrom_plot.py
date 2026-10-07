# simplifies chromosome plots of TE insertion mutations from Retroseq

#!/usr/bin/env python3
"""
Non-reference TE insertions (RetroSeq) - single mirrored chromosome map, Supp. Fig. 7 A/B style.

  28°C (control) insertions sit to the LEFT of each chromosome, 34°C (heat) to the RIGHT.
  Points are coloured by TE class (rocket palette, as in the original retroseq plot);
  use --colour group to colour by temperature instead (blue 28°C / red 34°C, as in the SNP map).
  Open chromatin (DANIO-CODE DOPEs track) is drawn as dark-red ticks across the chromosome.
  Only insertions UNIQUE to one temperature are plotted and counted: an insertion is treated
  as shared (and dropped from both sets) if the other temperature has an insertion on the same
  chromosome within --window bp (default 100 bp, any TE class).
  Numbers above each chromosome: unique insertions at 28°C (blue) | 34°C (red).
TE class = first element of MEINFO NAME (e.g. "LTR-DNA-hybrid" -> LTR).

Usage:
  Just run it: by default it reads, from the same folder as this script,
    ctrl_female_exp_unique_8_win100_gq750_fl8.vcf          (28°C)
    female_exp_unique_8_win100_filterpy_gq750_fl8.vcf      (34°C)
    daniocode_hub_280355_dopes_all.txt                     (open chromatin; skipped if absent)
  and writes TE_map_mirrored_ovaries.pdf/.png/.svg there.
  Override with --ctrl / --heat / --dopes / --tissue / --colour / --out.
"""
import argparse, re
from pathlib import Path
from collections import Counter
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "FreeSans", "DejaVu Sans"],
    "font.weight": "bold", "axes.labelweight": "bold", "axes.titleweight": "bold",
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "svg.fonttype": "none",
})
HEAT, CTRL, GREY, OPEN = "#C0392B", "#2E6DA4", "#D9D9D9", "darkred"
TE_ORDER = ["DNA", "SINE", "LINE", "LTR", "RC", "Satellite", "Unknown"]
TE_COL = {"DNA": "#35193e", "SINE": "#701f57", "LINE": "#ad1759", "LTR": "#e13342",
          "RC": "#f37651", "Satellite": "#f6b48f", "Unknown": "#9E9E9E"}
GRCZ11 = {"1": 59578282, "2": 59640629, "3": 62628489, "4": 78093715, "5": 72500376,
          "6": 60270059, "7": 74282399, "8": 54304671, "9": 56459846, "10": 45420867,
          "11": 45484837, "12": 49182954, "13": 52186027, "14": 52660232, "15": 48040578,
          "16": 55266484, "17": 53461100, "18": 51023478, "19": 48449771, "20": 55201332,
          "21": 45934066, "22": 39133080, "23": 46223584, "24": 42172926, "25": 37502051}
OFF = 0.22   # horizontal offset of points from the chromosome centre


def te_class(name):
    return {"DNA": "DNA", "LINE": "LINE", "LTR": "LTR", "SINE": "SINE", "RC": "RC",
            "SATELLITE": "Satellite"}.get(name.split("-")[0].upper(), "Unknown")


def read_vcf(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"): continue
            c = line.rstrip("\n").split("\t")
            info = dict(kv.split("=", 1) for kv in c[7].split(";") if "=" in kv)
            name = info.get("MEINFO", "").split(",")[0]
            fmt = dict(zip(c[8].split(":"), c[9].split(":"))) if len(c) > 9 else {}
            rows.append(dict(chr=c[0], pos=int(c[1]), te=name, te_class=te_class(name),
                             GQ=fmt.get("GQ"), SP=fmt.get("SP")))
    d = pd.DataFrame(rows)
    n_all = len(d)
    d = d[d.chr.str.fullmatch(r"\d+")].copy()
    d["chrn"] = d.chr.astype(int)
    return d[d.chrn.between(1, 25)], n_all


def mark_shared(a, b, window):
    """True for rows of a with an insertion in b on the same chromosome within window bp."""
    out = pd.Series(False, index=a.index)
    for ch, g in a.groupby("chrn"):
        q = np.sort(b.loc[b.chrn == ch, "pos"].to_numpy())
        if len(q) == 0: continue
        i = np.searchsorted(q, g.pos.to_numpy())
        lo = np.abs(g.pos.to_numpy() - q[np.clip(i - 1, 0, len(q) - 1)])
        hi = np.abs(q[np.clip(i, 0, len(q) - 1)] - g.pos.to_numpy())
        out.loc[g.index] = np.minimum(lo, hi) <= window
    return out


def read_dopes(path):
    """Same format as the original script: header line, then chrom start end name score strand."""
    out = []
    with open(path) as f:
        next(f)
        for line in f:
            p = line.split()
            if len(p) < 3: continue
            ch = p[0].replace("chr", "")
            if ch.isdigit() and 1 <= int(ch) <= 25:
                out.append((int(ch), int(p[1]), int(p[2])))
    return out


def main():
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--ctrl", default=str(here / "ctrl_male_exp_unique_8_win100_gq750_fl8.vcf"))
    ap.add_argument("--heat", default=str(here / "male_exp_unique_8_win100_filterpy_gq750_fl8.vcf"))
    ap.add_argument("--dopes", default=str(here / "daniocode_hub_280355_dopes_all.txt"))
    ap.add_argument("--tissue", default="Testes")
    ap.add_argument("--window", type=int, default=100,
                    help="bp distance within which a 28C and a 34C insertion count as the same")
    ap.add_argument("--colour", choices=["class", "group"], default="class")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    out = a.out or str(here / f"TE_map_mirrored_{a.tissue.lower()}")
    for f in (a.ctrl, a.heat):
        if not Path(f).exists(): raise SystemExit(f"file not found: {f}")
    ctrl, n_c = read_vcf(a.ctrl); heat, n_h = read_vcf(a.heat)
    sh_c, sh_h = mark_shared(ctrl, heat, a.window), mark_shared(heat, ctrl, a.window)
    for nm, d, sh in (("28°C", ctrl, sh_c), ("34°C", heat, sh_h)):
        print(f"{nm}: {len(d)} insertions on chr 1-25, {sh.sum()} shared with the other "
              f"temperature (within {a.window} bp), {(~sh).sum()} unique")
    ctrl, heat = ctrl[~sh_c].copy(), heat[~sh_h].copy()
    for nm, d in (("28°C unique", ctrl), ("34°C unique", heat)):
        cc = Counter(d.te_class)
        print(f"{nm}: {len(d)}; "
              + ", ".join(f"{k} {cc[k]}" for k in TE_ORDER if cc[k]))
    dopes = read_dopes(a.dopes) if Path(a.dopes).exists() else []
    print(f"open chromatin: {len(dopes)} regions" if dopes else
          f"open chromatin file not found ({a.dopes}) - plotting without it")

    chroms = list(range(1, 26)); lengths = GRCZ11
    maxlen = max(lengths[str(c)] for c in chroms) / 1e6
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    for c in chroms:
        ax.add_patch(mpl.patches.FancyBboxPatch((c - 0.08, 0), 0.16, lengths[str(c)] / 1e6,
                     boxstyle="round,pad=0,rounding_size=0.08", fc=GREY, ec="none", zorder=1))
    if dopes:
        segs = [[(c - 0.08, (s + e) / 2e6), (c + 0.08, (s + e) / 2e6)] for c, s, e in dopes]
        ax.add_collection(LineCollection(segs, colors=OPEN, linewidths=0.3, zorder=2))
    for d, off, gcol in ((ctrl, -OFF, CTRL), (heat, OFF, HEAT)):
        for cls in TE_ORDER:   # common classes first, rarer ones on top
            s = d[d.te_class == cls]
            col = TE_COL[cls] if a.colour == "class" else gcol
            ax.scatter(s.chrn + off, s.pos / 1e6, s=9, c=col, alpha=0.8,
                       lw=0.25, ec="white", zorder=3 + TE_ORDER.index(cls))
    nc, nh = ctrl.chrn.value_counts(), heat.chrn.value_counts()
    for c in chroms:
        ax.text(c - 0.04, -2.2, str(nc.get(c, 0)), ha="right", va="bottom", fontsize=5.5, color=CTRL)
        ax.text(c + 0.04, -2.2, str(nh.get(c, 0)), ha="left", va="bottom", fontsize=5.5, color=HEAT)
    ax.set_xlim(0.4, 25.6); ax.set_ylim(maxlen * 1.04, -5)
    ax.set_xticks(chroms); ax.set_xlabel("Chromosome"); ax.set_ylabel("Position (Mb)")
    ax.spines["bottom"].set_visible(False); ax.tick_params(axis="x", length=0)

    if a.colour == "class":
        present = [c for c in TE_ORDER if (pd.concat([ctrl, heat]).te_class == c).any()]
        h = [Line2D([], [], ls="", marker="o", ms=5, mfc=TE_COL[c], mec="none", label=c) for c in present]
    else:
        h = [Line2D([], [], ls="", marker="o", ms=5, mfc=CTRL, mec="none", label="28°C"),
             Line2D([], [], ls="", marker="o", ms=5, mfc=HEAT, mec="none", label="34°C")]
    if dopes:
        h.append(Line2D([], [], color=OPEN, lw=1.5, label="open chromatin"))
    ax.legend(handles=h, loc="lower right", ncol=len(h), frameon=False, handletextpad=0.2,
              columnspacing=0.9, bbox_to_anchor=(1.0, -0.02))
    sym = "♀" if a.tissue.lower().startswith("ov") else "♂"
    ax.set_title(f"{a.tissue} {sym}   temperature-unique non-reference TE insertions", loc="left", pad=16)
    ax.text(0.0, 1.035, f"left: unique to 28°C (n = {len(ctrl)})", color=CTRL,
            transform=ax.transAxes, fontsize=7)
    ax.text(0.33, 1.035, f"right: unique to 34°C (n = {len(heat)})", color=HEAT,
            transform=ax.transAxes, fontsize=7)

    pd.concat([ctrl.assign(group="28C"), heat.assign(group="34C")]).drop(columns=["chrn"]) \
        .to_csv(f"{out}_insertions.tsv", sep="\t", index=False)
    for ext in ["pdf", "png", "svg"]:
        fig.savefig(f"{out}.{ext}", dpi=600 if ext == "png" else None, bbox_inches="tight")


if __name__ == "__main__":
    main()
