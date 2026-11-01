#!/usr/bin/env python3
"""
Depth-adjusted stats for PRIVATE (de novo) TE insertions, unique across all 24 fish.

Private counts (from the --scope all run) and total reads per fish (samtools flagstat)
are built in below, so this runs with no input files. For each tissue it fits:
  negative-binomial GLM:  private ~ temperature + offset(log total reads)
tests temperature by likelihood-ratio, and reports the rate ratio (34 vs 28 C) with a
95% Wald CI. A Wilcoxon rank-sum test on raw counts is reported alongside. The same is
done per TE class, with Benjamini-Hochberg FDR across classes within each tissue.

Why the offset: heat-exposed testes were sequenced more deeply than controls
(536 vs 448 M reads, p = 0.004), and depth predicts private-insertion count, so raw
counts are not comparable. Modelling per 100 M reads removes this confounder.
"""
import numpy as np
from scipy.stats import mannwhitneyu, chi2
from scipy.optimize import minimize
from scipy.special import gammaln

# ----------------------------------------------------------------------
# DATA (edit here if counts change). Keys are fish tags <num><F|M><C|T>.
#   F = ovary, M = testis ; C = 28 C (control), T = 34 C (heat)
# ----------------------------------------------------------------------
PRIVATE = {
    '1FC': 52, '5FC': 197, '9FC': 236, '12FC': 20, '16FC': 87, '25FC': 166,   # ovary 28
    '2FT': 154, '6FT': 141, '10FT': 127, '11FT': 155, '15FT': 68, '22FT': 102, # ovary 34
    '3MC': 221, '7MC': 56, '18MC': 65, '23MC': 105, '30MC': 112, '34MC': 139,  # testis 28
    '4MT': 253, '8MT': 247, '24MT': 130, '28MT': 60, '29MT': 190, '41MT': 231, # testis 34
}
READS = {  # total reads from samtools flagstat
    '1FC': 346028092, '2FT': 479072924, '3MC': 512177392, '4MT': 585406812,
    '5FC': 576179461, '6FT': 504249843, '7MC': 443406190, '8MT': 547840763,
    '9FC': 584869558, '10FT': 525245716, '11FT': 577406957, '12FC': 317064310,
    '15FT': 425067471, '16FC': 397619924, '18MC': 371799440, '22FT': 439018570,
    '23MC': 432967597, '24MT': 497902886, '25FC': 539422604, '28MT': 549271960,
    '29MT': 513095886, '30MC': 463250366, '34MC': 466225989, '41MT': 522394969,
}
# Optional: per-fish, per-class private counts for the class-level tests.
# Fill this in from Fig3C_private_TE_counts.tsv if you want per-class stats.
# Format: PRIVATE_BY_CLASS[tag] = {'DNA':..,'LTR':..,'LINE':..,'SINE':..,'RC':..,'Satellite':..}
PRIVATE_BY_CLASS = {}

TE_ORDER = ["DNA", "LTR", "LINE", "SINE", "RC", "Satellite"]


def group(tag):
    sex = 'Ovary' if 'F' in tag else 'Testis'
    temp = '28' if tag.endswith('C') else '34'
    return sex, temp


# ---- negative binomial GLM with log-reads offset --------------------
def _nb_nll(par, y, X, off):
    b = par[:-1]; la = np.clip(par[-1], -12, 5); r = np.exp(-la)
    mu = np.exp(np.clip(X @ b + off, -30, 30))
    return -np.sum(gammaln(y + r) - gammaln(r) - gammaln(y + 1)
                   + r * np.log(r / (r + mu)) + y * np.log(mu / (r + mu) + 1e-300))


def _nb_fit(y, X, off):
    p0 = np.r_[np.log(max(y.mean(), 0.5)) - off.mean(), np.zeros(X.shape[1] - 1), np.log(0.3)]
    r = minimize(_nb_nll, p0, args=(y, X, off), method="Nelder-Mead",
                 options=dict(maxiter=20000, xatol=1e-8, fatol=1e-10))
    return minimize(_nb_nll, r.x, args=(y, X, off), method="BFGS")


def nb_test(y, heat, reads):
    """Rate ratio (34 vs 28), 95% Wald CI, LRT p for the temperature term."""
    y = np.asarray(y, float); off = np.log(np.asarray(reads, float))
    X1 = np.c_[np.ones(len(y)), heat]; X0 = X1[:, :1]
    f1, f0 = _nb_fit(y, X1, off), _nb_fit(y, X0, off)
    p = chi2.sf(max(2 * (f0.fun - f1.fun), 0), 1)
    x = f1.x.copy(); la = np.clip(x[-1], -12, 5); r_ = np.exp(-la)
    mu = np.exp(X1 @ x[:-1] + off); W = mu * r_ * (r_ + y) / (r_ + mu) ** 2
    try: se = float(np.sqrt(np.linalg.inv(X1.T @ (W[:, None] * X1))[1, 1]))
    except Exception: se = np.nan
    b = x[1]
    return np.exp(b), np.exp(b - 1.96 * se), np.exp(b + 1.96 * se), p


def bh(p):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(q, 1); return out


def run(label, counts_by_tag):
    print(f"\n===== {label} =====")
    print(f"{'tissue':8s}{'n28':>4s}{'n34':>4s}{'mean28':>8s}{'mean34':>8s}"
          f"{'/100M_28':>10s}{'/100M_34':>10s}{'RR':>6s}{'95%CI':>14s}{'NB_p':>8s}{'MWU_p':>8s}")
    out_rows = []
    for sex in ['Ovary', 'Testis']:
        tags = [t for t in counts_by_tag if group(t)[0] == sex]
        y = np.array([counts_by_tag[t] for t in tags], float)
        rd = np.array([READS[t] for t in tags], float)
        heat = np.array([1 if group(t)[1] == '34' else 0 for t in tags])
        y28, y34 = y[heat == 0], y[heat == 1]
        p100_28 = (y28 / (rd[heat == 0] / 1e8)).mean()
        p100_34 = (y34 / (rd[heat == 1] / 1e8)).mean()
        rr, lo, hi, pnb = nb_test(y, heat, rd)
        mwu = mannwhitneyu(y34, y28, alternative='two-sided').pvalue
        print(f"{sex:8s}{len(y28):>4d}{len(y34):>4d}{y28.mean():>8.1f}{y34.mean():>8.1f}"
              f"{p100_28:>10.1f}{p100_34:>10.1f}{rr:>6.2f}{f'{lo:.2f}-{hi:.2f}':>14s}"
              f"{pnb:>8.3f}{mwu:>8.3f}")
        out_rows.append((sex, pnb, mwu))
    return out_rows


def depth_check():
    print("\n--- sequencing depth per group (were heat testes deeper?) ---")
    for sex in ['Ovary', 'Testis']:
        tags = [t for t in PRIVATE if group(t)[0] == sex]
        rd28 = [READS[t] / 1e6 for t in tags if group(t)[1] == '28']
        rd34 = [READS[t] / 1e6 for t in tags if group(t)[1] == '34']
        p = mannwhitneyu(rd34, rd28, alternative='two-sided').pvalue
        print(f"{sex}: 28C {np.mean(rd28):.0f}M vs 34C {np.mean(rd34):.0f}M  (MWU p={p:.3f})")


def per_class():
    if not PRIVATE_BY_CLASS:
        print("\n(PRIVATE_BY_CLASS is empty - fill it from Fig3C_private_TE_counts.tsv "
              "for per-class tests.)")
        return
    print("\n===== PER TE CLASS (depth-adjusted NB-GLM, BH-FDR within tissue) =====")
    for sex in ['Ovary', 'Testis']:
        tags = [t for t in PRIVATE_BY_CLASS if group(t)[0] == sex]
        rd = np.array([READS[t] for t in tags], float)
        heat = np.array([1 if group(t)[1] == '34' else 0 for t in tags])
        rows, ps = [], []
        for cls in TE_ORDER:
            y = np.array([PRIVATE_BY_CLASS[t].get(cls, 0) for t in tags], float)
            if y.sum() < 20:  # skip very rare classes
                continue
            rr, lo, hi, pnb = nb_test(y, heat, rd)
            rows.append((cls, y[heat == 0].mean(), y[heat == 1].mean(), rr, lo, hi, pnb)); ps.append(pnb)
        q = bh(ps)
        print(f"\n  {sex}:")
        print(f"  {'class':10s}{'m28':>7s}{'m34':>7s}{'RR':>6s}{'95%CI':>14s}{'p':>8s}{'FDR':>8s}")
        for (cls, m28, m34, rr, lo, hi, pnb), qq in zip(rows, q):
            print(f"  {cls:10s}{m28:>7.1f}{m34:>7.1f}{rr:>6.2f}{f'{lo:.2f}-{hi:.2f}':>14s}"
                  f"{pnb:>8.3f}{qq:>8.3f}")


if __name__ == "__main__":
    run("ALL TEs (private, unique across 24 fish, depth-adjusted)", PRIVATE)
    depth_check()
    per_class()
