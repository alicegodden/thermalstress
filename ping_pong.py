# Title: Ping-pong signature — corrected overlap counting + barplot
# Author: Dr. Alice M. Godden

import pysam
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from scipy.stats import ttest_ind, mannwhitneyu
import glob
import os

# =============================
# INPUT BAMs  (edit paths per sex)
# =============================
control_bams = sorted(glob.glob("maleC/*.bam"))
temp_bams    = sorted(glob.glob("maleT/*.bam"))
groups = {"control": control_bams, "temperature": temp_bams}

# overlap axis: 1..30 nt (1-based; a true ping-pong peaks at 10)
OVERLAP_MIN, OVERLAP_MAX = 1, 30
axis = np.arange(OVERLAP_MIN, OVERLAP_MAX + 1)      # 1..30
PEAK = 10

# =============================
# FUNCTION: ping-pong overlap histogram (1-based)
# =============================
def compute_pingpong(bam_file):
    sense, antisense = {}, {}
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            if read.is_unmapped:
                continue
            length = read.query_length
            if length is None or length < 24 or length > 32:
                continue
            chrom  = read.reference_name
            strand = "-" if read.is_reverse else "+"
            start, end = read.reference_start, read.reference_end
            five_prime = start if strand == "+" else end - 1
            (sense if strand == "+" else antisense).setdefault(chrom, []).append(five_prime)

    overlaps = []
    for chrom in sense:
        if chrom not in antisense:
            continue
        s_positions = sorted(sense[chrom])
        a_positions = sorted(antisense[chrom])
        j = 0
        for s_pos in s_positions:
            while j < len(a_positions) and a_positions[j] < s_pos:
                j += 1
            k = j
            while k < len(a_positions):
                diff = a_positions[k] - s_pos
                if diff > OVERLAP_MAX:          # diff is 0-based gap
                    break
                if diff >= 0:
                    overlaps.append(diff + 1)   # <-- FIX: 1-based overlap length
                k += 1

    overlaps = np.array(overlaps)
    if len(overlaps) == 0:
        return None, None, 0, np.nan
    overlaps = overlaps[(overlaps >= OVERLAP_MIN) & (overlaps <= OVERLAP_MAX)]
    counts = Counter(overlaps)

    y = np.array([counts.get(i, 0) for i in axis])   # counts at overlaps 1..30
    observed_peak = counts.get(PEAK, 0)

    # z-score of the 10-nt peak vs the rest of the distribution
    bg = np.delete(y, np.where(axis == PEAK)[0][0])
    z = (observed_peak - np.mean(bg)) / np.std(bg) if (len(bg) and np.std(bg) > 0) else np.nan

    # ping-pong strength = peak / mean background (for the stats/boxplot)
    strength = observed_peak / np.mean(bg) if np.mean(bg) > 0 else np.nan
    return y, observed_peak, z, strength

# =============================
# RUN
# =============================
results = {"control": [], "temperature": []}
strengths = {"control": [], "temperature": []}

for group_name, bam_list in groups.items():
    print(f"\nProcessing {group_name} samples...\n", flush=True)
    for bam in bam_list:
        sample = os.path.basename(bam).replace(".bam", "")
        y, peak, z, strength = compute_pingpong(bam)
        if y is None:
            print(f"{sample}: NO OVERLAPS", flush=True)
            continue
        print(f"{sample}: overlaps={int(y.sum())}, peak(10nt)={peak}, Z={z:.2f}, strength={strength:.2f}", flush=True)
        results[group_name].append(y)
        strengths[group_name].append(strength)

        # per-sample barplot (1-based)
        plt.figure(figsize=(5, 4))
        plt.bar(axis, y, color="#4c72b0")
        plt.axvline(PEAK, color="red", linestyle="--", label="10 nt")
        plt.xlabel("5′ overlap (nt)")
        plt.ylabel("Count")
        plt.title(sample)
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"{sample}.pingpong.png", dpi=300)
        plt.close()

# =============================
# GROUP BARPLOT with error bars (mean ± SEM across replicates)
# =============================
plt.figure(figsize=(8, 5))
colours = {"control": "#4c72b0", "temperature": "#c44e52"}
width = 0.4
for i, (group_name, ys) in enumerate(results.items()):
    if len(ys) == 0:
        print(f"WARNING: no data in {group_name}", flush=True)
        continue
    ys = np.vstack(ys)
    mean_y = ys.mean(axis=0)
    sem_y  = ys.std(axis=0, ddof=1) / np.sqrt(ys.shape[0]) if ys.shape[0] > 1 else np.zeros_like(mean_y)
    plt.bar(axis + (i - 0.5) * width, mean_y, width=width,
            yerr=sem_y, capsize=2, color=colours.get(group_name), label=group_name, alpha=.9)
plt.axvline(PEAK, color="red", linestyle="--")
plt.xlabel("5′ overlap (nt)")
plt.ylabel("Mean count (± SEM)")
plt.title("Ping-pong signature (group average)")
plt.xticks(np.arange(2, 31, 2))
plt.legend()
plt.tight_layout()
plt.savefig("group_comparison_pingpong.png", dpi=300)
plt.close()

# =============================
# STATS on ping-pong strength (peak / background), per sample
# =============================
ctrl = np.array([s for s in strengths["control"]    if not np.isnan(s)])
temp = np.array([s for s in strengths["temperature"] if not np.isnan(s)])
out_lines = ["=== GROUP SUMMARY ===",
             f"Control n = {len(ctrl)}, mean = {np.mean(ctrl):.4f}" if len(ctrl) else "Control: no data",
             f"Temperature n = {len(temp)}, mean = {np.mean(temp):.4f}" if len(temp) else "Temperature: no data", ""]
if len(ctrl) > 1 and len(temp) > 1:
    t_stat, p_t = ttest_ind(ctrl, temp, equal_var=False)
    u_stat, p_u = mannwhitneyu(ctrl, temp, alternative="two-sided")
    out_lines += ["=== STATISTICAL TESTS (ping-pong strength) ===",
                  f"Welch t-test p = {p_t:.4e}",
                  f"Mann-Whitney U p = {p_u:.4e}"]
with open("pingpong_stats_results.txt", "w") as f:
    f.write("\n".join(out_lines) + "\n")
print("\n".join(out_lines))
