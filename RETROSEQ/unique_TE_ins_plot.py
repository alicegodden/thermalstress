#!/usr/bin/env python3
"""
Non-reference TE insertions (RetroSeq) - mirrored chromosome maps of insertions UNIQUE to one dataset.

An insertion is kept only if it is found in exactly ONE of the four datasets
(ovaries 28°C, ovaries 34°C, testes 28°C, testes 34°C). If any other dataset has a call on the
same chromosome within --window bp (default 100 bp, any TE class), it is dropped everywhere.
Calls at the same locus within one file (the same insertion reported under two TE classes,
e.g. DNA and LTR at one position) are collapsed to a single insertion first.

Note: each RetroSeq VCF here has a single merged genotype column (6 fish pooled per group), so
presence in individual fish cannot be checked from these files - uniqueness is per dataset.

Writes, for each tissue, TE_map_unique_<tissue>.pdf/.png/.svg + _insertions.tsv,
and TE_unique_summary.tsv (the results table). Same plot style as te_chromosome_map_mirrored.py.

Usage (defaults = files in the same folder as this script):
  python te_chromosome_map_unique.py
  python te_chromosome_map_unique.py --ov-ctrl a.vcf --ov-heat b.vcf --te-ctrl c.vcf --te-heat d.vcf
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


def collapse_within(d, window):
    """Merge calls within window bp of each other in one file into one insertion (first class kept)."""
    keep = []
    for ch, g in d.sort_values(["chrn", "pos"]).groupby("chrn"):
        last = -10**12
        for idx, p in zip(g.index, g.pos):
            if p - last > window: keep.append(idx)
            last = p
    return d.loc[keep]


def draw(ax, ctrl, heat, dopes, colour, tissue):
    chroms = list(range(1, 26)); lengths = GRCZ11
    maxlen = max(lengths[str(c)] for c in chroms) / 1e6
    for c in chroms:
        ax.add_patch(mpl.patches.FancyBboxPatch((c - 0.08, 0), 0.16, lengths[str(c)] / 1e6,
                     boxstyle="round,pad=0,rounding_size=0.08", fc=GREY, ec="none", zorder=1))
    if dopes:
        segs = [[(c - 0.08, (s + e) / 2e6), (c + 0.08, (s + e) / 2e6)] for c, s, e in dopes]
        ax.add_collection(LineCollection(segs, colors=OPEN, linewidths=0.3, zorder=2))
    for d, off, gcol in ((ctrl, -OFF, CTRL), (heat, OFF, HEAT)):
        for cls in TE_ORDER:
            s = d[d.te_class == cls]
            col = TE_COL[cls] if colour == "class" else gcol
            ax.scatter(s.chrn + off, s.pos / 1e6, s=9, c=col, alpha=0.8,
                       lw=0.25, ec="white", zorder=3 + TE_ORDER.index(cls))
    nc, nh = ctrl.chrn.value_counts(), heat.chrn.value_counts()
    for c in chroms:
        ax.text(c - 0.04, -2.2, str(nc.get(c, 0)), ha="right", va="bottom", fontsize=5.5, color=CTRL)
        ax.text(c + 0.04, -2.2, str(nh.get(c, 0)), ha="left", va="bottom", fontsize=5.5, color=HEAT)
    ax.set_xlim(0.4, 25.6); ax.set_ylim(maxlen * 1.04, -5)
    ax.set_xticks(chroms); ax.set_xlabel("Chromosome"); ax.set_ylabel("Position (Mb)")
    ax.spines["bottom"].set_visible(False); ax.tick_params(axis="x", length=0)
    if colour == "class":
        present = [c for c in TE_ORDER if (pd.concat([ctrl, heat]).te_class == c).any()]
        h = [Line2D([], [], ls="", marker="o", ms=5, mfc=TE_COL[c], mec="none", label=c) for c in present]
    else:
        h = [Line2D([], [], ls="", marker="o", ms=5, mfc=CTRL, mec="none", label="28°C"),
             Line2D([], [], ls="", marker="o", ms=5, mfc=HEAT, mec="none", label="34°C")]
    if dopes:
        h.append(Line2D([], [], color=OPEN, lw=1.5, label="open chromatin"))
    ax.legend(handles=h, loc="lower right", ncol=len(h), frameon=False, handletextpad=0.2,
              columnspacing=0.9, bbox_to_anchor=(1.0, -0.02))
    sym = "♀" if tissue.lower().startswith("ov") else "♂"
    ax.set_title(f"{tissue} {sym}   dataset-unique non-reference TE insertions", loc="left", pad=16)
    ax.text(0.0, 1.035, f"left: unique to 28°C (n = {len(ctrl)})", color=CTRL,
            transform=ax.transAxes, fontsize=7)
    ax.text(0.33, 1.035, f"right: unique to 34°C (n = {len(heat)})", color=HEAT,
            transform=ax.transAxes, fontsize=7)


def main():
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--ov-ctrl", default=str(here / "ctrl_female_exp_unique_8_win100_gq750_fl8.vcf"))
    ap.add_argument("--ov-heat", default=str(here / "female_exp_unique_8_win100_filterpy_gq750_fl8.vcf"))
    ap.add_argument("--te-ctrl", default=str(here / "ctrl_male_exp_unique_8_win100_gq750_fl8.vcf"))
    ap.add_argument("--te-heat", default=str(here / "male_exp_unique_8_win100_filterpy_gq750_fl8.vcf"))
    ap.add_argument("--dopes", default=str(here / "daniocode_hub_280355_dopes_all.txt"))
    ap.add_argument("--window", type=int, default=100,
                    help="bp distance within which two calls count as the same insertion")
    ap.add_argument("--colour", choices=["class", "group"], default="class")
    ap.add_argument("--outdir", default=str(here))
    a = ap.parse_args()
    sets = {("Ovaries", "28°C"): a.ov_ctrl, ("Ovaries", "34°C"): a.ov_heat,
            ("Testes", "28°C"): a.te_ctrl, ("Testes", "34°C"): a.te_heat}
    missing = [f for f in sets.values() if not Path(f).exists()]
    if missing:
        print("WARNING - not found, uniqueness checked only against the files present:\n  "
              + "\n  ".join(missing))
    data = {k: read_vcf(f)[0] for k, f in sets.items() if Path(f).exists()}
    raw = {k: len(d) for k, d in data.items()}
    data = {k: collapse_within(d, a.window) for k, d in data.items()}

    rows, unique = [], {}
    for k, d in data.items():
        others = [o for j, o in data.items() if j != k]
        hit = pd.Series(False, index=d.index)
        for o in others: hit |= mark_shared(d, o, a.window)
        unique[k] = d[~hit].copy()
        cc = Counter(unique[k].te_class)
        rows.append(dict(tissue=k[0], temperature=k[1], calls_chr1_25=raw[k],
                         loci_after_collapsing=len(d), found_in_another_dataset=int(hit.sum()),
                         unique=int((~hit).sum()), **{c: cc.get(c, 0) for c in TE_ORDER[:-1]}))
    tab = pd.DataFrame(rows)
    print(tab.to_string(index=False))
    tab.to_csv(Path(a.outdir) / "TE_unique_summary.tsv", sep="\t", index=False)

    dopes = read_dopes(a.dopes) if Path(a.dopes).exists() else []
    print(f"open chromatin: {len(dopes)} regions" if dopes else
          f"open chromatin file not found ({a.dopes}) - plotting without it")
    for tissue in ("Ovaries", "Testes"):
        if (tissue, "28°C") not in unique or (tissue, "34°C") not in unique: continue
        ctrl, heat = unique[(tissue, "28°C")], unique[(tissue, "34°C")]
        fig, ax = plt.subplots(figsize=(7.2, 4.2))
        draw(ax, ctrl, heat, dopes, a.colour, tissue)
        out = Path(a.outdir) / f"TE_map_unique_{tissue.lower()}"
        pd.concat([ctrl.assign(group="28C"), heat.assign(group="34C")]).drop(columns=["chrn"]) \
            .to_csv(f"{out}_insertions.tsv", sep="\t", index=False)
        for ext in ["pdf", "png", "svg"]:
            fig.savefig(f"{out}.{ext}", dpi=600 if ext == "png" else None, bbox_inches="tight")
        plt.close(fig)


if __name__ == "__main__":
    main()
