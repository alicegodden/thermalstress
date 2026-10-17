#!/usr/bin/env python3
"""
Are dataset-unique non-reference TE insertions enriched in open (or closed) chromatin?

Uses the SAME unique insertions as te_chromosome_map_unique.py: calls on chr 1-25, duplicate
calls at one locus within a file collapsed, then (default --scope tissue) an insertion is kept
only if the other-temperature dataset of the SAME tissue has no call within --window bp
(default 100). --scope all = unique vs all three other datasets.

Open chromatin = DANIO-CODE DOPEs regions (merged). An insertion counts as "open" if it lies
inside a region, or within --flanks bp of one (several flanks are tested; 0 = strict overlap).
Everything else is "closed".

Statistics (one row per dataset x TE class x flank in the results table):
  1. Exact binomial test: observed fraction of insertions in open chromatin vs the fraction
     of chr 1-25 covered by open chromatin (+flank). Gives fold enrichment (obs/exp);
     >1 = enriched in open chromatin, <1 = enriched in closed chromatin.
  2. Permutation test (default 10,000): insertions are placed at random within the same
     chromosomes (per-chromosome counts kept), so chromosome-level differences in insertion
     number and open-chromatin density are controlled for. Two-sided empirical p.
  3. 34°C vs 28°C within each tissue: Fisher's exact test on open/closed counts (odds ratio),
     i.e. does heat change where unique insertions land?
  Benjamini-Hochberg FDR is applied across all tests of each type.

Writes TE_open_chromatin_enrichment.tsv (tests 1-2) and TE_open_chromatin_heat_vs_ctrl.tsv (test 3).

Input (default): per-fish private insertions written by private_TE_per_fish_interim.py
(--private-root, i.e. miranda/vcf/private_uniq). All private insertions of the six fish of a group
are pooled for tests 1-2; test 3 adds a per-fish Wilcoxon test on % of insertions in open chromatin.
If --private-root does not exist, the four pooled VCFs (--ov-ctrl etc.) are used as before.

Usage (defaults = files in the same folder as this script):
  python te_open_chromatin_enrichment.py
  python te_open_chromatin_enrichment.py --flanks 0,1000,5000 --perm 10000
"""
import argparse, re
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import binomtest, fisher_exact, mannwhitneyu

GRCZ11 = {"1": 59578282, "2": 59640629, "3": 62628489, "4": 78093715, "5": 72500376,
          "6": 60270059, "7": 74282399, "8": 54304671, "9": 56459846, "10": 45420867,
          "11": 45484837, "12": 49182954, "13": 52186027, "14": 52660232, "15": 48040578,
          "16": 55266484, "17": 53461100, "18": 51023478, "19": 48449771, "20": 55201332,
          "21": 45934066, "22": 39133080, "23": 46223584, "24": 42172926, "25": 37502051}
TE_ORDER = ["DNA", "SINE", "LINE", "LTR", "RC", "Satellite", "Unknown"]


# ---------- insertions (identical logic to te_chromosome_map_unique.py) ----------
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
            rows.append(dict(chr=c[0], pos=int(c[1]), te=name, te_class=te_class(name)))
    d = pd.DataFrame(rows)
    d = d[d.chr.str.fullmatch(r"\d+")].copy()
    d["chrn"] = d.chr.astype(int)
    return d[d.chrn.between(1, 25)]


def collapse_within(d, window):
    keep = []
    for ch, g in d.sort_values(["chrn", "pos"]).groupby("chrn"):
        last = -10**12
        for idx, p in zip(g.index, g.pos):
            if p - last > window: keep.append(idx)
            last = p
    return d.loc[keep]


def mark_shared(a, b, window):
    out = pd.Series(False, index=a.index)
    for ch, g in a.groupby("chrn"):
        q = np.sort(b.loc[b.chrn == ch, "pos"].to_numpy())
        if len(q) == 0: continue
        p = g.pos.to_numpy(); i = np.searchsorted(q, p)
        lo = np.abs(p - q[np.clip(i - 1, 0, len(q) - 1)])
        hi = np.abs(q[np.clip(i, 0, len(q) - 1)] - p)
        out.loc[g.index] = np.minimum(lo, hi) <= window
    return out


# ---------- open chromatin ----------
def read_dopes(path):
    out = {}
    with open(path) as f:
        next(f)
        for line in f:
            p = line.split()
            if len(p) < 3: continue
            ch = p[0].replace("chr", "")
            if ch.isdigit() and 1 <= int(ch) <= 25:
                out.setdefault(int(ch), []).append((int(p[1]), int(p[2])))
    return out


def merged(regions, flank):
    """Pad by flank, clip to chromosome, merge overlaps -> {chr: (starts, ends)}, covered bp."""
    res, cov = {}, 0
    for ch, iv in regions.items():
        L = GRCZ11[str(ch)]
        iv = sorted((max(0, s - flank), min(L, e + flank)) for s, e in iv if s < L)
        S, E = [], []
        for s, e in iv:
            if S and s <= E[-1]: E[-1] = max(E[-1], e)
            else: S.append(s); E.append(e)
        res[ch] = (np.array(S), np.array(E)); cov += sum(e - s for s, e in zip(S, E))
    return res, cov


def in_open(chrs, pos, iv):
    """Boolean array: position inside a merged interval (half-open [start, end))."""
    out = np.zeros(len(pos), bool)
    for ch in np.unique(chrs):
        if ch not in iv: continue
        S, E = iv[ch]; m = chrs == ch; p = pos[m]
        i = np.searchsorted(S, p, side="right") - 1
        out[m] = (i >= 0) & (p < E[np.clip(i, 0, None)])
    return out


def perm_test(d, iv, n_perm, rng):
    """Random positions on the same chromosomes; returns null counts of insertions in open chromatin."""
    counts = d.chrn.value_counts()
    null = np.zeros(n_perm, int)
    for ch, n in counts.items():
        if ch not in iv: continue
        S, E = iv[ch]
        p = rng.integers(0, GRCZ11[str(ch)], size=(n_perm, n))
        i = np.searchsorted(S, p, side="right") - 1
        null += ((i >= 0) & (p < E[np.clip(i, 0, None)])).sum(1)
    return null


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(q, 1); return out


def read_private_root(root):
    """Private-only per-fish VCFs from private_TE_per_fish_interim.py:
    <root>/<ovaries|testes>/<control|treatment>/<fish>.private.vcf"""
    out = {}
    for f in sorted(Path(root).rglob("*.private.vcf")):
        parts = [x.lower() for x in f.relative_to(root).parts]
        tissue = "Ovaries" if parts[0].startswith("ov") else "Testes" if parts[0].startswith("test") else None
        temp = "28°C" if "control" in parts else "34°C" if "treatment" in parts else None
        if tissue is None or temp is None: continue
        d = read_vcf(f).assign(fish=f.name.replace(".private.vcf", ""))
        out.setdefault((tissue, temp), []).append(d)
    return {k: pd.concat(v, ignore_index=True) for k, v in out.items()}


def main():
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser()
    ap.add_argument("--ov-ctrl", default=str(here / "ctrl_female_exp_unique_8_win100_gq750_fl8.vcf"))
    ap.add_argument("--ov-heat", default=str(here / "female_exp_unique_8_win100_filterpy_gq750_fl8.vcf"))
    ap.add_argument("--te-ctrl", default=str(here / "ctrl_male_exp_unique_8_win100_gq750_fl8.vcf"))
    ap.add_argument("--te-heat", default=str(here / "male_exp_unique_8_win100_filterpy_gq750_fl8.vcf"))
    ap.add_argument("--private-root",
                    default="/Users/alicegodden/PycharmProjects/heatstressvariants/miranda/vcf/private_uniq",
                    help="folder written by private_TE_per_fish_interim.py; used if it exists "
                         "(per-fish private insertions). Otherwise the four pooled VCFs below are used.")
    ap.add_argument("--dopes", default=str(here / "daniocode_hub_280355_dopes_all.txt"))
    ap.add_argument("--window", type=int, default=100)
    ap.add_argument("--scope", choices=["tissue", "all"], default="tissue",
                    help="tissue: unique vs the other temperature of the SAME tissue (default); "
                         "all: unique vs all other datasets, both tissues")
    ap.add_argument("--flanks", default="0,1000,5000", help="comma-separated bp distances")
    ap.add_argument("--perm", type=int, default=10000)
    ap.add_argument("--min-n", type=int, default=10, help="skip TE classes with fewer insertions")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--outdir", default=str(here))
    a = ap.parse_args()
    if not Path(a.dopes).exists(): raise SystemExit(f"open chromatin file not found: {a.dopes}")
    rng = np.random.default_rng(a.seed)
    flanks = [int(x) for x in a.flanks.split(",")]

    if a.private_root and Path(a.private_root).exists():
        unique = read_private_root(a.private_root)
        print(f"using per-fish private insertions from {a.private_root}")
        for k, d in sorted(unique.items()):
            print(f"  {k[0]} {k[1]}: {d.fish.nunique()} fish, {len(d)} private insertions")
    else:
        sets = {("Ovaries", "28°C"): a.ov_ctrl, ("Ovaries", "34°C"): a.ov_heat,
                ("Testes", "28°C"): a.te_ctrl, ("Testes", "34°C"): a.te_heat}
        missing = [f for f in sets.values() if not Path(f).exists()]
        if missing:
            print("WARNING - not found, uniqueness checked only against files present:\n  " + "\n  ".join(missing))
        data = {k: collapse_within(read_vcf(f), a.window) for k, f in sets.items() if Path(f).exists()}
        unique = {}
        for k, d in data.items():
            hit = pd.Series(False, index=d.index)
            for j, o in data.items():
                if j != k and (a.scope == "all" or j[0] == k[0]): hit |= mark_shared(d, o, a.window)
            unique[k] = d[~hit]

    dopes = read_dopes(a.dopes)
    genome = sum(GRCZ11.values())
    print(f"open chromatin: {sum(map(len, dopes.values()))} regions on chr 1-25")

    rows, cmp_rows = [], []
    for flank in flanks:
        iv, cov = merged(dopes, flank)
        exp = cov / genome
        print(f"flank {flank} bp: open chromatin covers {cov/1e6:.2f} Mb = {100*exp:.2f}% of chr 1-25")
        for (tissue, temp), d in unique.items():
            d = d.assign(open=in_open(d.chrn.to_numpy(), d.pos.to_numpy(), iv))
            for cls in ["All"] + TE_ORDER:
                s = d if cls == "All" else d[d.te_class == cls]
                n, k = len(s), int(s.open.sum())
                if n == 0 or (cls != "All" and n < a.min_n): continue
                b = binomtest(k, n, exp)
                null = perm_test(s, iv, a.perm, rng)
                m = null.mean()
                p_perm = (1 + np.sum(np.abs(null - m) >= abs(k - m))) / (a.perm + 1)
                rows.append(dict(tissue=tissue, temperature=temp, te_class=cls, flank_bp=flank,
                                 n_unique=n, n_open=k, n_closed=n - k,
                                 pct_open=round(100 * k / n, 2), pct_expected=round(100 * exp, 2),
                                 perm_mean_open=round(m, 1),
                                 fold_enrichment=round((k / n) / exp, 3) if exp else np.nan,
                                 direction="open" if k / n > exp else "closed",
                                 p_binom=b.pvalue, p_perm=p_perm))
            # heat vs control within tissue (all classes)
        for tissue in ("Ovaries", "Testes"):
            if (tissue, "28°C") not in unique or (tissue, "34°C") not in unique: continue
            c = [r for r in rows if r["tissue"] == tissue and r["te_class"] == "All" and r["flank_bp"] == flank]
            c28 = next(r for r in c if r["temperature"] == "28°C")
            c34 = next(r for r in c if r["temperature"] == "34°C")
            orr, p = fisher_exact([[c34["n_open"], c34["n_closed"]], [c28["n_open"], c28["n_closed"]]])
            pw = np.nan
            if "fish" in unique[(tissue, "28°C")]:
                fr = {}
                for t in ("28°C", "34°C"):
                    g = unique[(tissue, t)]
                    o = in_open(g.chrn.to_numpy(), g.pos.to_numpy(), iv)
                    fr[t] = pd.Series(o, index=g.index).groupby(g.fish).mean() * 100
                pw = mannwhitneyu(fr["34°C"], fr["28°C"], alternative="two-sided").pvalue
            cmp_rows.append(dict(tissue=tissue, flank_bp=flank,
                                 open_34C=c34["n_open"], closed_34C=c34["n_closed"],
                                 open_28C=c28["n_open"], closed_28C=c28["n_closed"],
                                 pct_open_34C=c34["pct_open"], pct_open_28C=c28["pct_open"],
                                 odds_ratio_34_vs_28=round(orr, 3), p_fisher=p,
                                 p_wilcoxon_per_fish=pw))

    res = pd.DataFrame(rows)
    res["padj_binom"] = bh(res.p_binom); res["padj_perm"] = bh(res.p_perm)
    for c in ("p_binom", "padj_binom", "p_perm", "padj_perm"):
        res[c] = res[c].map(lambda x: float(f"{x:.3g}"))
    res.to_csv(Path(a.outdir) / "TE_open_chromatin_enrichment.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
    print("\n=== Unique insertions in open chromatin vs genome expectation ===")
    print(res.to_string(index=False))
    if cmp_rows:
        cmp = pd.DataFrame(cmp_rows); cmp["padj_fisher"] = bh(cmp.p_fisher)
        if cmp.p_wilcoxon_per_fish.notna().all():
            cmp["padj_wilcoxon_per_fish"] = bh(cmp.p_wilcoxon_per_fish)
        else:
            cmp = cmp.drop(columns="p_wilcoxon_per_fish")
        for c in [c for c in cmp if c.startswith(("p_", "padj"))]:
            cmp[c] = cmp[c].map(lambda x: float(f"{x:.3g}"))
        cmp.to_csv(Path(a.outdir) / "TE_open_chromatin_heat_vs_ctrl.tsv", sep="\t", index=False)
        print("\n=== 34°C vs 28°C: proportion of unique insertions in open chromatin ===")
        print(cmp.to_string(index=False))


if __name__ == "__main__":
    main()
