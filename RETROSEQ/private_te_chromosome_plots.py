#!/usr/bin/env python3
"""
Non-reference TE insertions (RetroSeq) from PER-FISH VCFs - mirrored chromosome maps of unique
insertions, judged WITHIN each tissue (ovaries and testes handled separately; fish of the other
tissue are never used). Same plot style as te_chromosome_map_mirrored.py:
28°C insertions LEFT of each chromosome, 34°C RIGHT, coloured by TE class.

--mode private (default): an insertion is plotted only if it is called in ONE fish of that tissue.
   If any other fish of the same tissue (either temperature) has a call on the same chromosome
   within --window bp (default 100 bp, any TE class), it is dropped.
--mode temperature: an insertion is plotted if no fish of the OTHER temperature (same tissue) has
   it; it may be shared by several fish at the same temperature (plotted once).
Nearby calls in the same fish are kept as separate calls (merge them with --collapse).
Chromosomes 1-25 only.

Input: per-fish RetroSeq VCFs (one per fish, or multi-sample VCFs with one column per fish).
Fish, tissue and temperature are read from sample/file names such as "..._12FC_..." or
"..._22MT_..." (F/M = ovaries/testes, C/T = 28°C/34°C), or from --samplesheet
(TSV: sample  tissue  temperature  vcf). Pooled/merged samples are skipped with a warning.

Usage:
  python te_map_unique_per_fish.py            # reads --root folder layout (default: your miranda/vcf)
  python te_map_unique_per_fish.py --vcfs "per_fish/*.vcf"
  python te_map_unique_per_fish.py --vcfs "per_fish/*.vcf" --mode temperature
Writes per tissue: TE_map_unique_perfish_<tissue>.pdf/.png/.svg + _insertions.tsv,
plus TE_unique_perfish_summary.tsv (per group) and TE_unique_perfish_counts.tsv (per fish).
"""
import argparse, glob, re
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


ID_RE = re.compile(r"_(\d+)([FM])([CT])_")


def te_class(name):
    return {"DNA": "DNA", "LINE": "LINE", "LTR": "LTR", "SINE": "SINE", "RC": "RC",
            "SATELLITE": "Satellite"}.get(name.split("-")[0].upper(), "Unknown")


def parse_id(name):
    m = ID_RE.findall("_" + name + "_")
    if len(m) != 1:
        return None
    num, sex, grp = m[0]
    return dict(fish=f"{num}{sex}{grp}", tissue="Ovaries" if sex == "F" else "Testes",
                temperature="28°C" if grp == "C" else "34°C")


def _open(path):
    import gzip
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def group_from_path(path):
    """Tissue/temperature from folders like .../ovaries/control/x.vcf or .../testes/treatment/x.vcf."""
    parts = [p.lower() for p in Path(path).parts]
    tissue = ("Ovaries" if any(p.startswith(("ovar", "ovary", "female")) for p in parts) else
              "Testes" if any(p.startswith(("test", "male")) for p in parts) else None)
    temp = ("28°C" if any(p in ("control", "ctrl", "28", "28c") for p in parts) else
            "34°C" if any(p in ("treatment", "treated", "heat", "temp", "34", "34c") for p in parts)
            else None)
    return tissue, temp


def read_root(root, fl=8, min_gq=0.0):
    """Read per-fish VCFs from root/{ovaries,testes}/{control,treatment}/*.vcf[.gz].
    One fish per file (fish ID = file name). Keeps calls with FL == fl (0 = no FL filter),
    GQ >= min_gq and a non-reference genotype."""
    vcfs = sorted(p for p in Path(root).rglob("*") if p.name.endswith((".vcf", ".vcf.gz"))
                  and "private_uniq" not in p.parts            # skip outputs of the interim script
                  and not p.name.endswith((".private.vcf", ".private.vcf.gz")))
    rows, skipped, log = [], [], []
    for path in vcfs:
        tissue, temp = group_from_path(path.relative_to(root))
        if tissue is None or temp is None:
            skipped.append(str(path)); continue
        fish = path.name.replace(".vcf.gz", "").replace(".vcf", "")
        n_all = n_keep = 0
        with _open(path) as fh:
            for line in fh:
                if line.startswith("#"): continue
                c = line.rstrip("\n").split("\t")
                if not re.fullmatch(r"\d+", c[0]) or not 1 <= int(c[0]) <= 25: continue
                n_all += 1
                if len(c) > 9:
                    fmt = dict(zip(c[8].split(":"), c[9].split(":")))
                    if "1" not in fmt.get("GT", "1"): continue
                    if fl and fmt.get("FL") not in (None, ".") and int(float(fmt["FL"])) != fl: continue
                    gq = fmt.get("GQ")
                    if min_gq and gq not in (None, ".") and float(gq) < min_gq: continue
                info = dict(kv.split("=", 1) for kv in c[7].split(";") if "=" in kv)
                name = info.get("MEINFO", "").split(",")[0]
                rows.append(dict(sample=fish, fish=fish, tissue=tissue, temperature=temp,
                                 chrn=int(c[0]), pos=int(c[1]), te_class=te_class(name)))
                n_keep += 1
        log.append((tissue, temp, fish, n_all, n_keep))
    if skipped:
        print("SKIPPED (could not tell tissue/temperature from folder names):\n  " + "\n  ".join(skipped))
    if not rows: raise SystemExit(f"no usable per-fish VCFs found under {root}")
    lg = pd.DataFrame(log, columns=["tissue", "temperature", "fish", "calls_chr1_25", "kept_after_filters"])
    print(f"read {len(lg)} per-fish VCFs from {root}  (filters: FL == {fl or 'any'}, GQ >= {min_gq})")
    print(lg.sort_values(["tissue", "temperature", "fish"]).to_string(index=False))
    return pd.DataFrame(rows)


def read_calls(vcfs, sheet=None):
    """Return one row per (fish, insertion call)."""
    rows = []
    for path in vcfs:
        samples = None
        with open(path) as fh:
            for line in fh:
                if line.startswith("#CHROM"):
                    samples = line.rstrip("\n").split("\t")[9:]
                    if len(samples) == 1 and sheet is None and parse_id(samples[0]) is None:
                        samples = [Path(path).name]          # fall back to file name
                    continue
                if line.startswith("#"): continue
                c = line.rstrip("\n").split("\t")
                if samples is None:                           # no #CHROM header line
                    samples = [Path(path).name] if len(c) == 10 else \
                              [f"{Path(path).name}:col{i}" for i in range(len(c) - 9)]
                if not re.fullmatch(r"\d+", c[0]) or not 1 <= int(c[0]) <= 25: continue
                info = dict(kv.split("=", 1) for kv in c[7].split(";") if "=" in kv)
                name = info.get("MEINFO", "").split(",")[0]
                for s, val in zip(samples, c[9:]):
                    gt = val.split(":")[0]
                    if "1" not in gt: continue                # absent in this fish
                    rows.append(dict(sample=s, chrn=int(c[0]), pos=int(c[1]),
                                     te_class=te_class(name)))
    d = pd.DataFrame(rows)
    if sheet is not None:
        meta = pd.read_csv(sheet, sep="\t", dtype=str)
        meta["temperature"] = meta.temperature.str.replace("°C", "").str.strip() + "°C"
        meta = meta.rename(columns={"sample": "sample"}).set_index("sample")
        d = d.join(meta[["tissue", "temperature"]], on="sample")
        d["fish"] = d["sample"]
    else:
        ids = {s: parse_id(s) for s in d["sample"].unique()}
        bad = [s for s, v in ids.items() if v is None]
        if bad:
            print("SKIPPED - could not read fish/tissue/temperature from these sample names "
                  "(pooled/merged sample or not a per-fish VCF; use --vcfs or --samplesheet):\n  "
                  + "\n  ".join(bad))
            d = d[~d["sample"].isin(bad)].copy()
            if d.empty:
                raise SystemExit("No per-fish samples left. Point --vcfs at the per-fish RetroSeq "
                                 "VCFs, e.g. --vcfs \"per_fish/*.vcf\"")
        for k in ("fish", "tissue", "temperature"):
            d[k] = d["sample"].map(lambda s: ids[s][k])
    return d


def collapse_per_fish(d, window):
    d = d.sort_values(["fish", "chrn", "pos"])
    keep, last = [], {}
    for idx, f_, ch, p in zip(d.index, d.fish, d.chrn, d.pos):
        if (f_, ch) in last and p - last[(f_, ch)] <= window:
            last[(f_, ch)] = p; continue
        keep.append(idx); last[(f_, ch)] = p
    return d.loc[keep].copy()


def flag_unique(d, window, mode):
    """Within each tissue: mark calls as unique (private to one fish, or absent from other temperature)."""
    d = d.copy(); d["unique"] = True
    for _, g in d.groupby(["tissue", "chrn"]):
        o = np.argsort(g.pos.to_numpy(), kind="stable")
        pos = g.pos.to_numpy()[o]; fish = g.fish.to_numpy()[o]
        temp = g.temperature.to_numpy()[o]; idx = g.index.to_numpy()[o]
        bad = np.zeros(len(pos), bool)
        for i in range(len(pos)):
            j = i + 1
            while j < len(pos) and pos[j] - pos[i] <= window:
                clash = (fish[j] != fish[i]) if mode == "private" else (temp[j] != temp[i])
                if clash: bad[i] = bad[j] = True
                j += 1
        d.loc[idx[bad], "unique"] = False
    return d


def one_per_locus(d, window):
    """For plotting/counting at group level: merge the same insertion carried by several fish."""
    keep = []
    for _, g in d.sort_values(["chrn", "pos"]).groupby("chrn"):
        last = -10**12
        for idx, p in zip(g.index, g.pos):
            if p - last > window: keep.append(idx)
            last = p
    return d.loc[keep]


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


def draw(ax, ctrl, heat, dopes, colour, tissue, title, nf=(None, None)):
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
    ax.set_title(f"{tissue} {sym}   {title}", loc="left", pad=16)
    ax.text(0.0, 1.035, (f"left: 28°C, all {nf[0]} fish summed (n = {len(ctrl)})" if nf[0] else f"left: unique to 28°C (n = {len(ctrl)})"), color=CTRL,
            transform=ax.transAxes, fontsize=7)
    ax.text(0.45, 1.035, (f"right: 34°C, all {nf[1]} fish summed (n = {len(heat)})" if nf[1] else f"right: unique to 34°C (n = {len(heat)})"), color=HEAT,
            transform=ax.transAxes, fontsize=7)


def main():
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcfs", required=False, default=None, help='glob of per-fish VCFs, e.g. "per_fish/*.vcf"')
    ap.add_argument("--samplesheet", default=None)
    ap.add_argument("--root", default="/Users/alicegodden/PycharmProjects/heatstressvariants/miranda/vcf",
                    help="folder containing ovaries/{control,treatment} and testes/{control,treatment} "
                         "with one VCF per fish (used unless --vcfs or --samplesheet is given)")
    ap.add_argument("--collapse", action="store_true",
                    help="merge calls within --window bp in the SAME fish into one (off by default: every call kept)")
    ap.add_argument("--fl", type=int, default=8, help="keep RetroSeq calls with this FL flag (0 = no filter)")
    ap.add_argument("--min-gq", type=float, default=0, help="minimum GQ per call (pooled files used 750)")
    ap.add_argument("--mode", choices=["private", "temperature"], default="private")
    ap.add_argument("--dopes", default=str(here / "daniocode_hub_280355_dopes_all.txt"))
    ap.add_argument("--window", type=int, default=100)
    ap.add_argument("--colour", choices=["class", "group"], default="class")
    ap.add_argument("--outdir", default=str(here))
    a = ap.parse_args()
    if a.samplesheet or a.vcfs:
        vcfs = (sorted(pd.read_csv(a.samplesheet, sep="\t", dtype=str).vcf.unique()) if a.samplesheet
                else sorted(glob.glob(a.vcfs)))
        if not vcfs: raise SystemExit("no VCFs found")
        print(f"reading {len(vcfs)} VCF(s)")
        d = read_calls(vcfs, a.samplesheet)
    else:
        d = read_root(a.root, a.fl, a.min_gq)
    if a.collapse: d = collapse_per_fish(d, a.window)
    d = flag_unique(d, a.window, a.mode)
    fish = d[["fish", "tissue", "temperature"]].drop_duplicates().sort_values(["tissue", "temperature", "fish"])
    print("fish per group:\n" + fish.groupby(["tissue", "temperature"]).size().to_string())

    u = d[d.unique]
    per_fish = (u.groupby("fish").size().reindex(fish.fish, fill_value=0).rename("unique").reset_index()
                .merge(fish, on="fish"))
    per_fish = per_fish.merge(d.groupby("fish").size().rename("calls").reset_index(), on="fish")
    per_fish[["fish", "tissue", "temperature", "calls", "unique"]].to_csv(
        Path(a.outdir) / "TE_unique_perfish_counts.tsv", sep="\t", index=False)

    dopes = read_dopes(a.dopes) if Path(a.dopes).exists() else []
    print(f"open chromatin: {len(dopes)} regions" if dopes else
          f"open chromatin file not found ({a.dopes}) - plotting without it")
    title = ("fish-private non-reference TE insertions" if a.mode == "private"
             else "temperature-unique non-reference TE insertions")
    rows = []
    for tissue in ("Ovaries", "Testes"):
        # private mode: every private call is its own insertion (no merging);
        # temperature mode: the same insertion carried by several fish is plotted/counted once
        pick = (lambda x: one_per_locus(x, a.window)) if a.mode == "temperature" else (lambda x: x)
        grp = {t: pick(u[(u.tissue == tissue) & (u.temperature == t)]) for t in ("28°C", "34°C")}
        if all(len(fish[(fish.tissue == tissue) & (fish.temperature == t)]) == 0 for t in grp): continue
        for t, g in grp.items():
            cc = Counter(g.te_class)
            nf = len(fish[(fish.tissue == tissue) & (fish.temperature == t)])
            rows.append(dict(tissue=tissue, temperature=t, n_fish=nf, unique=len(g),
                             mean_per_fish=round(per_fish[(per_fish.tissue == tissue) &
                                                          (per_fish.temperature == t)].unique.mean(), 2),
                             **{c: cc.get(c, 0) for c in TE_ORDER[:-1]}))
        ctrl, heat = grp["28°C"], grp["34°C"]
        fig, ax = plt.subplots(figsize=(7.2, 4.2))
        nfish = tuple(len(fish[(fish.tissue == tissue) & (fish.temperature == t)]) for t in ("28°C", "34°C"))
        draw(ax, ctrl, heat, dopes, a.colour, tissue, title, nfish if a.mode == "private" else (None, None))
        out = Path(a.outdir) / f"TE_map_unique_perfish_{tissue.lower()}"
        pd.concat([ctrl.assign(group="28C"), heat.assign(group="34C")]).drop(columns=["chrn"]) \
            .to_csv(f"{out}_insertions.tsv", sep="\t", index=False)
        for ext in ["pdf", "png", "svg"]:
            fig.savefig(f"{out}.{ext}", dpi=600 if ext == "png" else None, bbox_inches="tight")
        plt.close(fig)
    tab = pd.DataFrame(rows)
    pd.set_option("display.width", 200)
    print(f"\nmode: {a.mode} (within tissue, window {a.window} bp)")
    print(tab.to_string(index=False))
    tab.to_csv(Path(a.outdir) / "TE_unique_perfish_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
