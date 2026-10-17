#!/usr/bin/env python3
"""
Figure 3C (remake) - truly unique (private) non-reference TE insertions per fish,
28°C vs 34°C within ovaries and within testes.

Definition of "unique": an insertion called in ONE fish only, judged within tissue (default
--scope tissue). If any other fish OF THE SAME TISSUE (either temperature) has a call on the same
chromosome within --window bp (default 100 bp, any TE class), the insertion is dropped from
those fish. Fish of the other tissue are ignored. --scope all compares against all 24 fish. Calls at one locus in the same fish reported
under two TE classes are merged only with --collapse (default: every call kept). Chromosomes 1-25 only.

Statistics (per tissue, for all TEs and for each TE class):
  With --reads (recommended): negative binomial GLM  private ~ temperature + offset(log(reads)),
  temperature tested by likelihood ratio; rate ratio 34 vs 28°C with 95% Wald CI. The plot then
  shows private insertions per 100 million reads (depth-adjusted).
  Without --reads: Wilcoxon rank-sum test on raw counts, and the plot shows raw counts.
  A Wilcoxon test on raw counts is always reported as well. Benjamini-Hochberg FDR across all
  tests of each type (both tissues, all classes).

--reads: samtools flagstat output, either the *.flagstat.txt files (glob) or one text file with
  lines like  "<file>.flagstat.txt:346028092 + 0 in total ...". Fish are matched on the tag
  _<number><F|M><C|T>_ (e.g. _12FC_) present in both the VCF and flagstat names.

Input: per-fish RetroSeq VCFs (one sample per file, or multi-sample VCFs with one column per
fish). Tissue and temperature are read from sample/file names such as "..._12FC_..." or
"..._22MT_..." (number, F/M = ovaries/testes, C/T = 28°C/34°C). Alternatively give
--samplesheet: TSV with columns  sample  tissue(Ovaries/Testes)  temperature(28/34)  vcf
(use the VCF sample column name, or the file name if the VCF has one sample).

Usage:
  python fig3c_private_TE_per_fish.py            # reads --root folder layout (default: your miranda/vcf)
  python fig3c_private_TE_per_fish.py --reads "flagstat/*.flagstat.txt"   # depth-adjusted (recommended)
  python fig3c_private_TE_per_fish.py --vcfs "per_fish_vcfs/*.vcf"
  python fig3c_private_TE_per_fish.py --samplesheet samples.tsv
Writes Fig3C_private_TEs.pdf/.png/.svg, Fig3C_private_TE_counts.tsv (per fish) and
Fig3C_private_TE_stats.tsv.
"""
import argparse, glob, re
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "FreeSans", "DejaVu Sans"],
    "font.weight": "bold", "axes.labelweight": "bold", "axes.titleweight": "bold",
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42, "svg.fonttype": "none",
})
HEAT, CTRL = "#C0392B", "#2E6DA4"
TE_ORDER = ["DNA", "LTR", "LINE", "SINE", "RC", "Satellite", "Unknown"]
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


def private_insertions(d, window, scope="tissue", collapse=False):
    """Optionally collapse within-fish duplicates, then keep calls with no call in any OTHER fish within window."""
    d = d.sort_values(["fish", "chrn", "pos"])
    keep, last = ([], {}) if collapse else (list(d.index), None)
    for idx, f, ch, p in (zip(d.index, d.fish, d.chrn, d.pos) if collapse else []):
        if (f, ch) in last and p - last[(f, ch)] <= window:
            last[(f, ch)] = p; continue
        keep.append(idx); last[(f, ch)] = p
    d = d.loc[keep].copy()
    d["private"] = True
    groups = ["chrn"] if scope == "all" else ["tissue", "chrn"]
    for _, g in d.groupby(groups):
        pos, fish = g.pos.to_numpy(), g.fish.to_numpy()
        o = np.argsort(pos); pos, fish, idx = pos[o], fish[o], g.index.to_numpy()[o]
        shared = np.zeros(len(pos), bool)
        for i in range(len(pos)):
            j = i + 1
            while j < len(pos) and pos[j] - pos[i] <= window:
                if fish[j] != fish[i]: shared[i] = shared[j] = True
                j += 1
        d.loc[idx[shared], "private"] = False
    return d


TAG_RE = re.compile(r"_(\d+[FM][CT])_")


def fish_tag(name):
    m = TAG_RE.search("_" + str(name) + "_")
    return m.group(1) if m else str(name)


def read_reads(spec):
    """Total reads per fish tag from samtools flagstat files or a grep-style summary file."""
    out = {}
    for f in sorted(glob.glob(spec)):
        with open(f) as fh:
            for line in fh:
                m = re.search(r"(\d+) \+ (\d+) in total", line)
                if not m: continue
                src = line.split(":")[0] if ":" in line.split(" ")[0] else Path(f).name
                out[fish_tag(src)] = int(m.group(1)) + int(m.group(2))
                if src == Path(f).name: break
    return out


def _nb_nll(par, y, X, off):
    from scipy.special import gammaln
    b = par[:-1]; la = np.clip(par[-1], -12, 5); r = np.exp(-la)
    mu = np.exp(np.clip(X @ b + off, -30, 30))
    return -np.sum(gammaln(y + r) - gammaln(r) - gammaln(y + 1)
                   + r * np.log(r / (r + mu)) + y * np.log(mu / (r + mu) + 1e-300))


def _nb_fit(y, X, off):
    from scipy.optimize import minimize
    p0 = np.r_[np.log(max(y.mean(), 0.5)) - off.mean(), np.zeros(X.shape[1] - 1), np.log(0.3)]
    r = minimize(_nb_nll, p0, args=(y, X, off), method="Nelder-Mead",
                 options=dict(maxiter=20000, xatol=1e-8, fatol=1e-10))
    return minimize(_nb_nll, r.x, args=(y, X, off), method="BFGS")


def nb_test(y, heat, reads):
    """NB GLM y ~ heat + offset(log reads): rate ratio, 95% CI, LRT p, dispersion."""
    from scipy.stats import chi2
    y = np.asarray(y, float); off = np.log(np.asarray(reads, float))
    X1 = np.c_[np.ones(len(y)), heat]; X0 = X1[:, :1]
    f1, f0 = _nb_fit(y, X1, off), _nb_fit(y, X0, off)
    p = chi2.sf(max(2 * (f0.fun - f1.fun), 0), 1)
    # SE of the temperature coefficient from the observed information of the regression
    # coefficients (dispersion held at its estimate; robust when dispersion is near 0)
    x = f1.x.copy(); la = np.clip(x[-1], -12, 5); r_ = np.exp(-la)
    mu = np.exp(X1 @ x[:-1] + off)
    W = mu * r_ * (r_ + y) / (r_ + mu) ** 2
    try: se = float(np.sqrt(np.linalg.inv(X1.T @ (W[:, None] * X1))[1, 1]))
    except Exception: se = np.nan
    if not np.isfinite(x[1]) or abs(x[1]) > 20: se = np.nan        # all-zero group: no finite CI
    b = x[1]
    return np.exp(b), np.exp(b - 1.96 * se), np.exp(b + 1.96 * se), p, float(np.exp(np.clip(x[-1], -12, 5)))


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(q, 1); return out


def boot_ci(x, rng, n=10000):
    x = np.asarray(x, float)
    if len(x) < 2: return (np.nan, np.nan)
    m = rng.choice(x, size=(n, len(x)), replace=True).mean(1)
    return np.percentile(m, [2.5, 97.5])


def stars(q):
    return "***" if q < 0.001 else "**" if q < 0.01 else "*" if q < 0.05 else ""


def draw(ax, counts, tissue, cats, stat, rng, ylabel=True, ylab="Unique insertions per fish"):
    sub = counts[counts.tissue == tissue]
    x = np.arange(len(cats))
    for temp, col, off in (("28°C", CTRL, -0.17), ("34°C", HEAT, 0.17)):
        g = sub[sub.temperature == temp]
        for xi, cat in zip(x, cats):
            v = g[cat].to_numpy()
            ax.scatter(xi + off + rng.uniform(-0.05, 0.05, len(v)), v, s=8, c=col, alpha=0.3,
                       lw=0, zorder=2)
            lo, hi = boot_ci(v, rng)
            ax.errorbar(xi + off, v.mean(), yerr=[[v.mean() - lo], [hi - v.mean()]], fmt="o",
                        ms=4.5, c=col, ecolor=col, elinewidth=1, capsize=2, zorder=3)
    ymax = sub[cats].to_numpy().max() if len(sub) else 1
    for xi, cat in zip(x, cats):
        r = stat[(stat.tissue == tissue) & (stat.te_class == cat)]
        if len(r) and stars(r.padj.iloc[0]):
            top = sub[cat].max()
            ax.text(xi, top + ymax * 0.05, stars(r.padj.iloc[0]), ha="center", va="bottom", fontsize=9)
        elif len(r):
            top = sub[cat].max()
            ax.text(xi, top + ymax * 0.05, "ns", ha="center", va="bottom", fontsize=6, color="0.4")
    ax.set_xticks(x); ax.set_xticklabels(["All TEs" if c == "All" else c for c in cats])
    ax.set_ylim(0, ymax * 1.2 + 1)
    if ylabel: ax.set_ylabel(ylab)


def main():
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcfs", default=None, help="glob of per-fish VCFs (alternative to --root)")
    ap.add_argument("--samplesheet", default=None)
    ap.add_argument("--root", default="/Users/alicegodden/PycharmProjects/heatstressvariants/miranda/vcf",
                    help="folder containing ovaries/{control,treatment} and testes/{control,treatment} "
                         "with one VCF per fish (used unless --vcfs or --samplesheet is given)")
    ap.add_argument("--collapse", action="store_true",
                    help="merge calls within --window bp in the SAME fish into one (off by default: every call kept)")
    ap.add_argument("--fl", type=int, default=8, help="keep RetroSeq calls with this FL flag (0 = no filter)")
    ap.add_argument("--min-gq", type=float, default=0, help="minimum GQ per call (pooled files used 750)")
    ap.add_argument("--window", type=int, default=100)
    ap.add_argument("--scope", choices=["tissue", "all"], default="tissue",
                    help="tissue: unique vs all other fish of the SAME tissue (default); "
                         "all: unique vs all 24 fish")
    ap.add_argument("--min-total", type=int, default=20,
                    help="test/plot a TE class only if it has at least this many insertions overall")
    ap.add_argument("--reads", default=None,
                    help='flagstat files (glob, e.g. "flagstat/*.flagstat.txt") or one grep-style summary file')
    ap.add_argument("--out", default=str(here / "Fig3C_private_TE"))
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)

    if a.samplesheet or a.vcfs:
        vcfs = (sorted(pd.read_csv(a.samplesheet, sep="\t", dtype=str).vcf.unique()) if a.samplesheet
                else sorted(glob.glob(a.vcfs)))
        if not vcfs: raise SystemExit("no VCFs found")
        print(f"reading {len(vcfs)} VCF(s)")
        d = read_calls(vcfs, a.samplesheet)
    else:
        d = read_root(a.root, a.fl, a.min_gq)
    d = private_insertions(d, a.window, a.scope, a.collapse)
    p = d[d.private]

    fish = d[["fish", "tissue", "temperature"]].drop_duplicates()
    print("fish per group:\n", fish.groupby(["tissue", "temperature"]).size().to_string())
    counts = (p.groupby(["fish", "te_class"]).size().unstack(fill_value=0)
              .reindex(fish.fish, fill_value=0))
    counts.columns.name = None
    for c in TE_ORDER:
        if c not in counts: counts[c] = 0
    counts["All"] = counts[TE_ORDER].sum(axis=1)
    counts = fish.set_index("fish").join(counts).reset_index()

    reads = read_reads(a.reads) if a.reads else {}
    counts["tag"] = counts.fish.map(fish_tag)
    if reads:
        counts["reads_M"] = counts.tag.map(reads) / 1e6
        miss = counts[counts.reads_M.isna()].fish.tolist()
        if miss: raise SystemExit("no read count found for: " + ", ".join(miss))
        print("\nreads (millions) per group:\n" + counts.groupby(["tissue", "temperature"]).reads_M
              .agg(["mean", "std"]).round(0).to_string())
    counts.to_csv(f"{a.out}_counts.tsv", sep="\t", index=False)

    cats = ["All"] + [c for c in TE_ORDER if counts[c].sum() >= a.min_total]
    rows = []
    for tissue in ("Ovaries", "Testes"):
        s_ = counts[counts.tissue == tissue]
        if s_.temperature.nunique() < 2: continue
        heat = (s_.temperature == "34°C").astype(int).to_numpy()
        for cat in cats:
            c28, c34 = s_[heat == 0][cat], s_[heat == 1][cat]
            try:
                pw = mannwhitneyu(c34, c28, alternative="two-sided").pvalue
            except ValueError:
                pw = 1.0
            row = dict(tissue=tissue, te_class=cat, n_28C=len(c28), n_34C=len(c34),
                       total_28C=int(c28.sum()), total_34C=int(c34.sum()),
                       mean_28C=round(c28.mean(), 1), sd_28C=round(c28.std(), 1),
                       mean_34C=round(c34.mean(), 1), sd_34C=round(c34.std(), 1),
                       p_wilcoxon_raw=pw)
            if reads:
                r28 = c28 / s_[heat == 0].reads_M * 100; r34 = c34 / s_[heat == 1].reads_M * 100
                rr, lo, hi, pnb, alpha = nb_test(s_[cat], heat, s_.reads_M)
                row.update(per100M_28C=f"{r28.mean():.1f} ± {r28.std():.1f}",
                           per100M_34C=f"{r34.mean():.1f} ± {r34.std():.1f}",
                           rate_ratio=round(rr, 2), CI95=f"{lo:.2f}-{hi:.2f}",
                           p_nb_lrt=pnb, nb_dispersion=round(alpha, 3))
            rows.append(row)
    stat = pd.DataFrame(rows)
    stat["padj_wilcoxon_raw"] = bh(stat.p_wilcoxon_raw)
    pcol = "p_nb_lrt" if reads else "p_wilcoxon_raw"
    if reads: stat["padj_nb"] = bh(stat.p_nb_lrt)
    stat["padj"] = stat["padj_nb"] if reads else stat["padj_wilcoxon_raw"]
    stat["sig"] = stat.padj.map(stars)
    out = stat.copy()
    for c in [c for c in out if c.startswith(("p_", "padj"))]: out[c] = out[c].map(lambda x: float(f"{x:.3g}"))
    out.to_csv(f"{a.out}_stats.tsv", sep="\t", index=False)
    pd.set_option("display.width", 300); pd.set_option("display.max_columns", 40)
    print("\nPrivate insertions in total:", len(p), "of", len(d), "calls")
    print(f"test used for stars: {'negative binomial GLM with log(reads) offset, LRT' if reads else 'Wilcoxon on raw counts'}; BH-adjusted")
    print(out.drop(columns=["padj"]).to_string(index=False))

    plot = counts.copy()
    ylab = "Private insertions per fish"
    if reads:
        for c in cats: plot[c] = plot[c] / plot.reads_M * 100
        ylab = "Private insertions\nper 100 M reads"

    fig = plt.figure(figsize=(7.2, 2.8))
    gs = fig.add_gridspec(1, 2, wspace=0.25)
    for i, tissue in enumerate(("Ovaries", "Testes")):
        sgs = gs[i].subgridspec(1, 2, width_ratios=[1, len(cats) - 1], wspace=0.35)
        ax0 = fig.add_subplot(sgs[0]); ax1 = fig.add_subplot(sgs[1])
        draw(ax0, plot, tissue, ["All"], stat, rng, ylabel=True, ylab=ylab)
        draw(ax1, plot, tissue, cats[1:], stat, rng, ylabel=False)
        sym = "♀" if tissue == "Ovaries" else "♂"
        ax0.set_title(f"{tissue} {sym}", loc="left")
    h = [plt.Line2D([], [], ls="", marker="o", c=CTRL, label="28°C"),
         plt.Line2D([], [], ls="", marker="o", c=HEAT, label="34°C")]
    fig.legend(handles=h, loc="upper right", ncol=2, frameon=False, bbox_to_anchor=(0.99, 1.04))
    for ext in ["pdf", "png", "svg"]:
        fig.savefig(f"{a.out}.{ext}", dpi=600 if ext == "png" else None, bbox_inches="tight")


if __name__ == "__main__":
    main()
