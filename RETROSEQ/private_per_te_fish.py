#!/usr/bin/env python3
"""
Keep only TE insertions PRIVATE to ONE fish ACROSS ALL 24 FISH (both tissues, both
treatments), write a filtered VCF per fish, and sum by treatment.

An insertion is PRIVATE only if no OTHER fish in the whole experiment -- regardless of
tissue or temperature -- has a call within --window bp. This removes insertions shared
between any two individuals (including a male and a female), which are inherited
(shared-ancestry) polymorphisms rather than events arising during the experiment.
Ovaries and testes are still reported separately in the output tables.

Input layout (default --root):
  <root>/ovaries/control/*.vcf   <root>/ovaries/treatment/*.vcf
  <root>/testes/control/*.vcf    <root>/testes/treatment/*.vcf
  (one RetroSeq VCF per fish; .vcf or .vcf.gz; control = 28 C, treatment = 34 C)

Steps:
  1. Read calls on chr 1-25 with a non-reference genotype (and FL == --fl if FL present).
  2. (Only with --collapse) merge calls within --window bp in the same fish into one.
  3. Private = NO other fish in the ENTIRE dataset (all 24, both tissues, both treatments)
     has a call within --window bp (default 100 bp, any TE class). Everything else removed.
Outputs (in --outdir, default <root>/private_across_all24):
  <tissue>/<control|treatment>/<fish>.private.vcf
  private_TE_per_fish.tsv
  private_TE_by_treatment.tsv

Usage:
  python private_TE_across_all24.py --root /path/to/vcf --window 100
"""
import argparse, gzip, re
from pathlib import Path
import numpy as np
import pandas as pd

TE_ORDER = ["DNA", "LTR", "LINE", "SINE", "RC", "Satellite", "Unknown"]


def te_class(name):
    return {"DNA": "DNA", "LINE": "LINE", "LTR": "LTR", "SINE": "SINE", "RC": "RC",
            "SATELLITE": "Satellite"}.get(name.split("-")[0].upper(), "Unknown")


def opener(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def group_from_path(rel):
    parts = [p.lower() for p in rel.parts]
    tissue = ("Ovaries" if any(p.startswith(("ovar", "female")) for p in parts) else
              "Testes" if any(p.startswith(("test", "male")) for p in parts) else None)
    trt = ("control" if any(p in ("control", "ctrl") for p in parts) else
           "treatment" if any(p in ("treatment", "treated", "heat") for p in parts) else None)
    return tissue, trt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="/Users/alicegodden/PycharmProjects/heatstressvariants/miranda/vcf")
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--window", type=int, default=100)
    ap.add_argument("--collapse", action="store_true",
                    help="merge calls within --window bp in the SAME fish into one (off by default)")
    ap.add_argument("--fl", type=int, default=8, help="keep calls with this FL flag (0 = no filter)")
    a = ap.parse_args()
    root = Path(a.root)
    outdir = Path(a.outdir) if a.outdir else root / "private_across_all24"

    outdir.mkdir(parents=True, exist_ok=True)
    vcfs = sorted(p for p in root.rglob("*") if p.name.endswith((".vcf", ".vcf.gz"))
                  and outdir not in p.parents and "private_across_all24" not in p.parts
                  and "private_uniq" not in p.parts
                  and not p.name.endswith((".private.vcf", ".private.vcf.gz")))
    rows, headers, files = [], {}, []
    for p in vcfs:
        tissue, trt = group_from_path(p.relative_to(root))
        if tissue is None or trt is None:
            print("skipped (folder not recognised):", p); continue
        fish = re.sub(r"\.vcf(\.gz)?$", "", p.name)
        files.append(dict(fish=fish, tissue=tissue, treatment=trt, file=p.name))
        hdr = []
        with opener(p) as fh:
            for line in fh:
                if line.startswith("#"): hdr.append(line); continue
                c = line.rstrip("\n").split("\t")
                if not re.fullmatch(r"\d+", c[0]) or not 1 <= int(c[0]) <= 25: continue
                if len(c) > 9:
                    fmt = dict(zip(c[8].split(":"), c[9].split(":")))
                    if "1" not in fmt.get("GT", "1"): continue
                    if a.fl and fmt.get("FL") not in (None, ".") and int(float(fmt["FL"])) != a.fl: continue
                info = dict(kv.split("=", 1) for kv in c[7].split(";") if "=" in kv)
                rows.append(dict(fish=fish, tissue=tissue, treatment=trt, chrn=int(c[0]),
                                 pos=int(c[1]), te_class=te_class(info.get("MEINFO", "").split(",")[0]),
                                 line=line))
        headers[fish] = hdr
    if not rows: raise SystemExit(f"no calls read under {root}")
    fl = pd.DataFrame(files)
    n_fish_total = fl.fish.nunique()
    print(f"\nRead {n_fish_total} fish: "
          + ", ".join(f"{t} {len(g)}" for t, g in fl.groupby('tissue')))
    if n_fish_total != 24:
        print(f"  NOTE: expected 24 fish, found {n_fish_total} -- check the folder tree.")

    d = pd.DataFrame(rows)
    n_calls = d.groupby("fish").size()

    # 1) optional: collapse within fish (off by default)
    d = d.sort_values(["fish", "chrn", "pos"])
    keep, last = ([], {}) if a.collapse else (list(d.index), None)
    for idx, f, ch, p in (zip(d.index, d.fish, d.chrn, d.pos) if a.collapse else []):
        if (f, ch) in last and p - last[(f, ch)] <= a.window:
            last[(f, ch)] = p; continue
        keep.append(idx); last[(f, ch)] = p
    d = d.loc[keep].copy()
    n_coll = d.groupby("fish").size()

    # 2) PRIVATE ACROSS ALL 24 FISH -- group by chromosome only, so every fish in the
    #    experiment (both tissues, both treatments) is compared. An insertion is dropped
    #    if ANY other fish has a call within --window bp.
    d["private"] = True
    for _, g in d.groupby(["chrn"]):
        o = np.argsort(g.pos.to_numpy(), kind="stable")
        pos, fish, idx = g.pos.to_numpy()[o], g.fish.to_numpy()[o], g.index.to_numpy()[o]
        bad = np.zeros(len(pos), bool)
        for i in range(len(pos)):
            j = i + 1
            while j < len(pos) and pos[j] - pos[i] <= a.window:
                if fish[j] != fish[i]: bad[i] = bad[j] = True
                j += 1
        d.loc[idx[bad], "private"] = False
    priv = d[d.private]

    # 3) write per-fish private VCFs
    for f, g in priv.groupby("fish"):
        meta = fl[fl.fish == f].iloc[0]
        od = outdir / meta.tissue.lower() / meta.treatment
        od.mkdir(parents=True, exist_ok=True)
        with open(od / f"{f}.private.vcf", "w") as out:
            out.writelines(headers[f])
            out.writelines(g.sort_values(["chrn", "pos"]).line)

    # 4) tables
    cls = priv.groupby(["fish", "te_class"]).size().unstack(fill_value=0)
    for c in TE_ORDER:
        if c not in cls: cls[c] = 0
    per = fl[["fish", "tissue", "treatment"]].set_index("fish")
    per["calls"] = n_calls; per["after_collapsing"] = n_coll
    per["private"] = priv.groupby("fish").size()
    per = per.join(cls[TE_ORDER]).fillna(0)
    for c in ["calls", "after_collapsing", "private"] + TE_ORDER: per[c] = per[c].astype(int)
    per = per.reset_index().sort_values(["tissue", "treatment", "fish"])
    per.insert(3, "temperature", per.treatment.map({"control": "28C", "treatment": "34C"}))
    per.to_csv(outdir / "private_TE_per_fish.tsv", sep="\t", index=False)

    summ = per.groupby(["tissue", "treatment", "temperature"]).agg(
        n_fish=("fish", "size"), calls=("calls", "sum"), private_sum=("private", "sum"),
        private_mean=("private", "mean"), private_sd=("private", "std"),
        **{c: (c, "sum") for c in TE_ORDER}).reset_index()
    summ[["private_mean", "private_sd"]] = summ[["private_mean", "private_sd"]].round(1)
    summ = summ.loc[:, (summ != 0).any()]
    summ.to_csv(outdir / "private_TE_by_treatment.tsv", sep="\t", index=False)

    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 30)
    short = per.copy(); short["fish"] = short.fish.str.extract(r"_(\d+[FM][CT])_", expand=False).fillna(short.fish)
    print(f"\nPrivate insertions per fish (unique across ALL {n_fish_total} fish, window {a.window} bp):")
    print(short[["tissue", "temperature", "fish", "calls", "after_collapsing", "private"]].to_string(index=False))
    print("\nSummed by treatment:")
    print(summ.to_string(index=False))
    print(f"\nPrivate-only VCFs written to {outdir}")


if __name__ == "__main__":
    main()
