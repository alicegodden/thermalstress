#!/usr/bin/env python3
"""
Supplementary Figure 7 - HIGH-impact case-control variants (SnpEff + SnpSift CaseControl)

Layout: one chromosome map per sex (top), then summary panels:
  (i)   mirrored chromosome map: variants more frequent in controls on the left of each
        chromosome (blue), more frequent under heat on the right (red); size = -log10 p
  (ii)  alt-allele frequency, 28C vs 34C (bubble size = number of variants at that point)
  (iii) variant consequence, split by direction of frequency change

Usage:
  Just run it: by default it reads
    high_filtered_sift_CC_snpeff_female_allchrs_eff.vcf  (ovaries)
    high_filtered_sift_CC_snpeff_male_allchrs_eff.vcf    (testes)
  from the same folder as this script and writes SuppFig7.pdf/.png/.svg there.
  Override with --ovaries / --testes / --p / --out if needed.
Group is read from sample names (..._3MC_.. / ..._4MT_.. ; C = control, T = temperature);
if a VCF has no header, groups are inferred from the SnpSift Cases/Controls counts.
"""
import argparse, re, itertools
from collections import Counter
import numpy as np, pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import binomtest

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "FreeSans", "DejaVu Sans"],
    "font.weight": "bold",          # all text: tick labels, legends, annotations, gene names
    "axes.labelweight": "bold",     # axis labels
    "axes.titleweight": "bold",     # panel titles
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "svg.fonttype": "none",
})
HEAT, CTRL, GREY = "#C0392B", "#2E6DA4", "#D9D9D9"
EFFECT_ORDER = ["frameshift", "stop gained", "splice donor", "splice acceptor",
                "stop lost", "start lost", "other"]
EFFECT_COL = dict(zip(EFFECT_ORDER, ["#4C2A85", "#B5367A", "#E8743B", "#F2B134",
                                     "#3A9E8F", "#6BAED6", "#9E9E9E"]))
# GRCz11 chromosome lengths (bp), used when a VCF has no ##contig header lines
GRCZ11 = {"1": 59578282, "2": 59640629, "3": 62628489, "4": 78093715, "5": 72500376,
          "6": 60270059, "7": 74282399, "8": 54304671, "9": 56459846, "10": 45420867,
          "11": 45484837, "12": 49182954, "13": 52186027, "14": 52660232, "15": 48040578,
          "16": 55266484, "17": 53461100, "18": 51023478, "19": 48449771, "20": 55201332,
          "21": 45934066, "22": 39133080, "23": 46223584, "24": 42172926, "25": 37502051}
UNNAMED = re.compile(r"^(ENSDARG|si:|zgc:|CR\d|BX\d|CU\d|AL\d|FO\d|LT\d|CT\d|CABZ|wu:|im:|LOC)")


def effect_class(e):
    e = e.split("&")[0].replace("_variant", "").replace("_", " ")
    return e if e in EFFECT_ORDER else "other"


def _carriers(gts, idx):
    hom = het = 0
    for i in idx:
        g = gts[i].replace("|", "/")
        if "." in g: continue
        a = [int(x) > 0 for x in g.split("/")]
        hom += all(a); het += any(a) and not all(a)
    return hom, het


def infer_groups(records, n):
    """No header: find the sample split whose genotypes reproduce SnpSift Cases/Controls."""
    best = (-1, None)
    for S in itertools.combinations(range(n), n // 2):
        comp = [i for i in range(n) if i not in S]
        ok = sum(_carriers(g, S) == ca and _carriers(g, comp) == co for g, ca, co in records)
        best = max(best, (ok, S))
    ok, S = best
    print(f"  no sample header: inferred heat-group columns {[i + 1 for i in S]} "
          f"(matches Cases/Controls for {ok}/{len(records)} variants)")
    return ["T" if i in S else "C" for i in range(n)]


def read_vcf(path):
    lengths, raw, grp = {}, [], None
    with open(path) as fh:
        for line in fh:
            if line.startswith("##contig"):
                m = re.search(r"ID=([^,>]+),length=(\d+)", line)
                if m: lengths[m.group(1)] = int(m.group(2))
            elif line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                grp = []
                for s in samples:
                    m = re.search(r"_\d+[MF]([CT])_", s)
                    grp.append(m.group(1) if m else None)
                if None in grp:
                    raise SystemExit(f"Could not assign groups from sample names: {samples}")
            elif not line.startswith("#"):
                raw.append(line.rstrip("\n").split("\t"))
    if grp is None:
        recs = []
        for c in raw:
            info = dict(kv.split("=", 1) for kv in c[7].split(";") if "=" in kv)
            recs.append(([s.split(":")[0] for s in c[9:]],
                         tuple(map(int, info["Cases"].split(",")[:2])),
                         tuple(map(int, info["Controls"].split(",")[:2]))))
        grp = infer_groups(recs, len(raw[0]) - 9)
    rows = []
    for c in raw:
        if True:
            info = dict(kv.split("=", 1) for kv in c[7].split(";") if "=" in kv)
            anns = [a.split("|") for a in info["ANN"].split(",")]
            hi = next((a for a in anns if a[2] == "HIGH"), anns[0])
            # alt-allele frequency per group from genotypes (any non-ref allele counts)
            al = {"C": [], "T": []}
            for s, g in zip(c[9:], grp):
                gt = s.split(":")[0].replace("|", "/")
                if "." in gt: continue
                al[g] += [int(x) > 0 for x in gt.split("/")]
            p = info.get("CC_TREND", "nan")
            rows.append(dict(chr=c[0], pos=int(c[1]), gene=hi[3], effect=effect_class(hi[1]),
                             p=float(p) if p != "nan" else np.nan,
                             af_ctrl=np.mean(al["C"]) if al["C"] else np.nan,
                             af_heat=np.mean(al["T"]) if al["T"] else np.nan))
    d = pd.DataFrame(rows)
    d = d[d.chr.str.fullmatch(r"\d+")].copy()
    d["chrn"] = d.chr.astype(int)
    d["delta"] = d.af_heat - d.af_ctrl
    d["direction"] = np.select([d.delta > 0, d.delta < 0], ["heat", "ctrl"], "none")
    return d, lengths


def panel_map(ax, d, lengths, label_p=0.003):
    chroms = list(range(1, 26))
    maxlen = max(lengths.get(str(c), 0) for c in chroms) / 1e6
    for c in chroms:
        L = lengths.get(str(c), 0) / 1e6
        ax.add_patch(mpl.patches.FancyBboxPatch((c - 0.08, 0), 0.16, L,
                     boxstyle="round,pad=0,rounding_size=0.08", fc=GREY, ec="none", zorder=1))
    size = lambda p: 3 + 7 * (-np.log10(np.clip(p, 1e-6, 1)))
    for direc, col, off in [("ctrl", CTRL, -0.22), ("heat", HEAT, 0.22)]:
        s = d[d.direction == direc]
        ax.scatter(s.chrn + off, s.pos / 1e6, s=size(s.p), c=col, alpha=0.7,
                   lw=0.25, ec="white", zorder=3)
    # label only the strongest named genes
    lab = d[(d.p < label_p) & ~d.gene.str.match(UNNAMED)].sort_values("p").drop_duplicates("gene")
    for _, r in lab.iterrows():
        side = 0.22 if r.direction == "heat" else -0.22
        ax.annotate(r.gene, (r.chrn + side, r.pos / 1e6),
                    xytext=(r.chrn + side * 3.2, r.pos / 1e6 - 4), fontsize=6.5, fontstyle="italic",
                    ha="left" if side > 0 else "right", va="center",
                    bbox=dict(fc="white", ec="none", pad=0.3, alpha=0.8),
                    arrowprops=dict(arrowstyle="-", lw=0.4, color="0.25"), zorder=5)
    ax.set_xlim(0.4, 25.6); ax.set_ylim(maxlen * 1.04, -2)
    ax.set_xticks(chroms); ax.set_xlabel("Chromosome"); ax.set_ylabel("Position (Mb)")
    ax.spines["bottom"].set_visible(False); ax.tick_params(axis="x", length=0)
    h = [ax.scatter([], [], s=25, c=CTRL, label="higher at 28°C"),
         ax.scatter([], [], s=25, c=HEAT, label="higher at 34°C")]
    h += [ax.scatter([], [], s=size(p_), c="0.55", label=f"p = {p_:g}") for p_ in (0.05, 0.005, 0.0005)]
    ax.legend(handles=h, loc="lower right", ncol=5, frameon=False, handletextpad=0.1,
              columnspacing=0.9, bbox_to_anchor=(1.0, -0.02))


def panel_af(ax, d, title=None, size_legend=True):
    cnt = d.groupby(["af_ctrl", "af_heat", "direction"]).size().reset_index(name="n")
    col = cnt.direction.map({"heat": HEAT, "ctrl": CTRL, "none": "0.5"})
    ax.scatter(cnt.af_ctrl, cnt.af_heat, s=6 + 5 * cnt.n, c=col, alpha=0.6, lw=0.3, ec="white")
    ax.plot([0, 1], [0, 1], ls="--", lw=0.6, c="0.4")
    up, down = (d.delta > 0).sum(), (d.delta < 0).sum()
    ax.text(0.02, 1.02, f"higher at 34°C: {up}", color=HEAT, transform=ax.transAxes, fontsize=7)
    ax.text(1.06, 0.02, f"higher at\n28°C: {down}", color=CTRL, transform=ax.transAxes,
            ha="left", va="bottom", fontsize=7)
    if not size_legend: return _af_axes(ax, title)
    for n_ in (1, 5, 20):
        ax.scatter([], [], s=6 + 5 * n_, c="0.6", label=str(n_))
    ax.legend(title="variants", loc="center left", bbox_to_anchor=(1.0, 0.62), frameon=False,
              title_fontsize=6.5, labelspacing=1.0)
    _af_axes(ax, title)


def _af_axes(ax, title):
    ax.set_xlim(-0.05, 1.05); ax.set_ylim(-0.05, 1.05); ax.set_aspect("equal")
    ax.set_xticks([0, .5, 1]); ax.set_yticks([0, .5, 1])
    ax.set_xlabel("Alt allele freq., 28°C"); ax.set_ylabel("Alt allele freq., 34°C")
    if title: ax.set_title(title, loc="left", pad=12)


def effect_table(d):
    t = d[d.direction != "none"].groupby(["effect", "direction"]).size().unstack(fill_value=0)
    for c_ in ("heat", "ctrl"):
        if c_ not in t: t[c_] = 0
    return t


def panel_effect(ax, d, effects, xmax, title=None):
    """Diverging bars for one sex; effects / xmax shared across sexes so panels are comparable."""
    t = effect_table(d).reindex(effects[::-1], fill_value=0)
    y = np.arange(len(t))
    ax.barh(y, -t.ctrl, height=0.7, color=CTRL); ax.barh(y, t.heat, height=0.7, color=HEAT)
    off = xmax * 0.02
    for yi, c_, h_ in zip(y, t.ctrl, t.heat):
        ax.text(-c_ - off, yi, str(c_), ha="right", va="center", fontsize=6)
        ax.text(h_ + off, yi, str(h_), ha="left", va="center", fontsize=6)
    m = xmax * 1.3
    ax.set_yticks(y); ax.set_yticklabels(t.index)
    ax.axvline(0, c="0.3", lw=0.6); ax.set_xlim(-m, m)
    ticks = ax.get_xticks(); ax.set_xticks(ticks)
    ax.set_xticklabels([f"{abs(int(v))}" for v in ticks]); ax.set_xlim(-m, m)
    ax.set_xlabel("Number of variants"); ax.tick_params(axis="y", length=0)
    ax.text(0.25, 1.02, "higher at 28°C", color=CTRL, transform=ax.transAxes, ha="center", fontsize=7)
    ax.text(0.75, 1.02, "higher at 34°C", color=HEAT, transform=ax.transAxes, ha="center", fontsize=7)
    if title: ax.set_title(title, loc="left", pad=14)


def main():
    from pathlib import Path
    here = Path(__file__).resolve().parent          # folder containing this script
    ap = argparse.ArgumentParser()
    ap.add_argument("--ovaries", default=str(here / "high_filtered_sift_CC_snpeff_female_allchrs.eff.vcf"))
    ap.add_argument("--testes",  default=str(here / "high_filtered_sift_CC_snpeff_male_allchrs.eff.vcf"))
    ap.add_argument("--p", type=float, default=0.05, help="CC_TREND threshold")
    ap.add_argument("--out", default=str(here / "SuppFig7"))
    a = ap.parse_args()
    sets = [(n, f) for n, f in [("Ovaries", a.ovaries), ("Testes", a.testes)] if f]
    for name, f in sets:
        if not Path(f).exists():
            raise SystemExit(f"{name} file not found: {f}")
    n = len(sets)
    # full-page width (7.2 in): chromosome map per sex, then allele-frequency row, then consequence row
    fig = plt.figure(figsize=(7.2, 3.0 * n + 5.0))
    gs = fig.add_gridspec(n + 2, 1, height_ratios=[3.0] * n + [2.5, 2.2], hspace=0.5)
    row_af = gs[n].subgridspec(1, n, wspace=0.6)
    row_ef = gs[n + 1].subgridspec(1, n, wspace=0.55)
    letters = iter("ABCDEFGHI"); data = {}
    loaded = {name: read_vcf(f) for name, f in sets}
    shared = {}
    for _, (_, L) in loaded.items(): shared.update(L)
    for i, (name, f) in enumerate(sets):
        d, lengths = loaded[name]
        lengths = lengths or shared or GRCZ11
        n_all = len(d); d = d[d.p < a.p].copy(); data[name] = d
        b = binomtest(int((d.delta > 0).sum()), int((d.delta != 0).sum()))
        print(f"{name}: {n_all} variants in file, {len(d)} with CC_TREND < {a.p}; "
              f"{(d.delta>0).sum()} higher at 34C vs {(d.delta<0).sum()} higher at 28C "
              f"(binomial p = {b.pvalue:.2g}); {d.gene.nunique()} genes")
        d.drop(columns=["chrn"]).to_csv(f"{a.out}_{name.lower()}_variants.tsv", sep="\t", index=False)
        ax = fig.add_subplot(gs[i]); panel_map(ax, d, lengths)
        sym = "♀" if name == "Ovaries" else "♂"
        ax.set_title(f"{name} {sym}   {len(d)} HIGH-impact variants (CC_TREND p < {a.p})",
                     loc="left", fontweight="bold")
        ax.text(-0.07, 1.04, next(letters), transform=ax.transAxes, fontsize=11, fontweight="bold")
    for j, (name, _) in enumerate(sets):
        ax = fig.add_subplot(row_af[j])
        panel_af(ax, data[name], title=name, size_legend=(j == n - 1))
        ax.text(-0.42, 1.12, next(letters), transform=ax.transAxes, fontsize=11, fontweight="bold")
    effects = [e for e in EFFECT_ORDER if any((d.effect == e).any() for d in data.values())]
    xmax = max(effect_table(d).values.max() for d in data.values())
    for j, (name, _) in enumerate(sets):
        ax = fig.add_subplot(row_ef[j])
        panel_effect(ax, data[name], effects, xmax, title=name)
        ax.text(-0.42, 1.12, next(letters), transform=ax.transAxes, fontsize=11, fontweight="bold")
    for ext in ["pdf", "png", "svg"]:
        fig.savefig(f"{a.out}.{ext}", dpi=600 if ext == "png" else None, bbox_inches="tight")


if __name__ == "__main__":
    main()
