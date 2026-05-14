# Title: Stats analysis of ping pong signature
# Author: Dr. Alice M. Godden

# script 1

import pysam
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import glob
import os

# =============================
# INPUT BAMs
# =============================
control_bams = sorted(glob.glob("maleC/*.bam"))
temp_bams = sorted(glob.glob("maleT/*.bam"))

groups = {
    "control": control_bams,
    "temperature": temp_bams
}

# =============================
# FUNCTION: streaming ping-pong
# =============================
def compute_pingpong_streaming(bam_file):

    sense = {}
    antisense = {}

    print(f"Reading {bam_file}", flush=True)

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            if read.is_unmapped:
                continue

            length = read.query_length
            if length < 24 or length > 32:
                continue

            chrom = read.reference_name
            strand = "-" if read.is_reverse else "+"

            #  ^|^e CORRECT coordinates (fixes 9nt issue)
            start = read.reference_start
            end = read.reference_end

            if strand == "+":
                five_prime = start
            else:
                five_prime = end - 1

            if strand == "+":
                sense.setdefault(chrom, []).append(five_prime)
            else:
                antisense.setdefault(chrom, []).append(five_prime)

    overlaps = []

    # =============================
    # SAFE overlap calculation
    # =============================
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

                if diff > 30:
                    break

                if diff >= 0:
                    overlaps.append(diff)

                k += 1

    overlaps = np.array(overlaps)

    #  ^|^e handle empty case safely
    if len(overlaps) == 0:
        return None, None, 0, np.nan

    overlaps = overlaps[(overlaps >= 0) & (overlaps <= 30)]

    counts = Counter(overlaps)

    x = np.arange(0, 31)
    y = np.array([counts.get(i, 0) for i in x])

    observed_10 = counts.get(10, 0)

    background = np.delete(y, 10)
    if len(background) == 0 or np.std(background) == 0:
        z_score = np.nan
    else:
	z_score = (observed_10 - np.mean(background)) / np.std(background)

    return x, y, observed_10, z_score

# =============================
# RUN
# =============================
results = {"control": [], "temperature": []}

for group_name, bam_list in groups.items():
    print(f"\nProcessing {group_name} samples...\n", flush=True)

    for bam in bam_list:

        sample_name = os.path.basename(bam).replace(".bam", "")

        x, y, peak10, z = compute_pingpong_streaming(bam)

        if x is None:
            print(f"{sample_name}: NO OVERLAPS", flush=True)
            continue

        print(f"{sample_name}: overlaps={np.sum(y)}, peak10={peak10}, Z={z}", flush=True)

        results[group_name].append({
            "sample": sample_name,
            "x": x,
            "y": y,
            "peak10": peak10,
            "z": z
        })

        # =============================
        # SAVE INDIVIDUAL PLOT
        # =============================
        plt.figure(figsize=(5,4))
        plt.bar(x, y)
        plt.axvline(10, color='red', linestyle='--', label='10 nt')
        plt.xlabel("5 ^`  overlap (nt)")
        plt.ylabel("Count")
        plt.title(sample_name)
        plt.legend()
        plt.tight_layout()
        plt.savefig(f"{sample_name}.multi_pingpong.png", dpi=300)
        plt.close()

# =============================
# GROUP PLOT (FIXED)
# =============================
plt.figure(figsize=(7,5))

for group_name, group_data in results.items():

    if len(group_data) == 0:
        print(f"WARNING: No data in {group_name}", flush=True)
        continue

    ys = [d["y"] for d in group_data if d["y"] is not None and len(d["y"]) == 31]

    if len(ys) == 0:
        print(f"WARNING: No valid arrays in {group_name}", flush=True)
        continue

    ys = np.array(ys)
    mean_y = ys.mean(axis=0)

    plt.plot(range(31), mean_y, label=group_name, linewidth=2)

plt.axvline(10, color='red', linestyle='--')
plt.xlabel("5 ^`  overlap (nt)")
plt.ylabel("Mean count")
plt.title("Ping-pong signature (group average)")
plt.legend()
plt.yticks(np.arange(0, 2.25, 0.25) * 1e8)
plt.ylim(0, 2 * 1e8)
plt.tight_layout()
plt.savefig("group_comparison_pingpong.png", dpi=300)
plt.close()

# =============================
# SUMMARY
# =============================
print("\n=== SUMMARY ===", flush=True)

for group_name, group_data in results.items():

    if len(group_data) == 0:
        print(f"{group_name.upper()}: NO DATA", flush=True)
        continue

    z_scores = [d["z"] for d in group_data if not np.isnan(d["z"])]

    if len(z_scores) == 0:
        print(f"{group_name.upper()}: NO VALID Z-SCORES", flush=True)
        continue

    print(f"\n{group_name.upper()}", flush=True)
    print(f"Mean Z-score: {np.mean(z_scores):.2f}", flush=True)
    print(f"n samples: {len(z_scores)}", flush=True)




# script 2
import pysam
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import glob
import os
from scipy.stats import ttest_ind, mannwhitneyu

# =============================
# INPUT
# =============================
control_bams = sorted(glob.glob("femaleC/*.bam"))
temp_bams = sorted(glob.glob("femaleT/*.bam"))

groups = {
    "control": control_bams,
    "temperature": temp_bams
}

# =============================
# FUNCTION
# =============================
def compute_pingpong_strength(bam_file):

    sense = {}
    antisense = {}

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            if read.is_unmapped:
                continue

            length = read.query_length
            if length < 24 or length > 32:
                continue

            chrom = read.reference_name
            strand = "-" if read.is_reverse else "+"

            start = read.reference_start
            end = read.reference_end

            pos5 = start if strand == "+" else end - 1

            if strand == "+":
                sense.setdefault(chrom, []).append(pos5)
            else:
                antisense.setdefault(chrom, []).append(pos5)

    overlaps = []

    for chrom in sense:
        if chrom not in antisense:
            continue

        s = sorted(sense[chrom])
        a = sorted(antisense[chrom])

        j = 0

	for sp in s:
            while j < len(a) and a[j] < sp:
                j += 1

            k = j
            while k < len(a):
                diff = a[k] - sp

                if diff > 30:
                    break

                if diff >= 0:
                    overlaps.append(diff)

                k += 1

    overlaps = np.array(overlaps)

    if len(overlaps) == 0:
        return np.nan

    overlaps = overlaps[(overlaps >= 0) & (overlaps <= 30)]

    counts = Counter(overlaps)

    y = np.array([counts.get(i, 0) for i in range(31)])

    peak = y[9]  # your peak position

    background = np.delete(y, 9)

    if np.mean(background) == 0:
        return np.nan

    enrichment = peak / np.mean(background)

    return enrichment


# =============================
# RUN
# =============================
results = {"control": [], "temperature": []}
per_sample = []

for group, bam_list in groups.items():

    print(f"\nProcessing {group}\n")

    for bam in bam_list:
        name = os.path.basename(bam)

        val = compute_pingpong_strength(bam)

        print(f"{name}: enrichment = {val}")

        if not np.isnan(val):
            results[group].append(val)
            per_sample.append((name, group, val))


# =============================
# STATS
# =============================
ctrl = np.array(results["control"])
temp = np.array(results["temperature"])

t_stat, p_t = ttest_ind(ctrl, temp, equal_var=False)
u_stat, p_u = mannwhitneyu(ctrl, temp, alternative="two-sided")

# =============================
# WRITE OUTPUT FILE
# =============================
with open("pingpong_stats_results_female.txt", "w") as f:

    f.write("=== PER SAMPLE RESULTS ===\n\n")

    for name, group, val in per_sample:
        f.write(f"{name}\t{group}\t{val:.6f}\n")

    f.write("\n=== GROUP SUMMARY ===\n\n")

    f.write(f"Control n = {len(ctrl)}\n")
    f.write(f"Temperature n = {len(temp)}\n\n")

    f.write(f"Control mean = {np.mean(ctrl):.6f}\n")
    f.write(f"Temperature mean = {np.mean(temp):.6f}\n\n")

    f.write("=== STATISTICAL TESTS ===\n\n")

    f.write(f"T-test p-value = {p_t:.6e}\n")
    f.write(f"Mann-Whitney p-value = {p_u:.6e}\n")

print("\nResults written to pingpong_stats_results.txt")

# =============================
# PLOT
# =============================
plt.figure(figsize=(5,5))

plt.boxplot([ctrl, temp], labels=["control", "temperature"])
plt.scatter([1]*len(ctrl), ctrl)
plt.scatter([2]*len(temp), temp)

plt.ylabel("Ping-pong strength (peak / background)")
plt.title("Ping-pong enrichment")

plt.text(1.5, max(np.concatenate([ctrl, temp]))*1.1,
         f"p = {p_u:.3e}", ha='center')

plt.tight_layout()
plt.savefig("pingpong_stats_female.png", dpi=300)
plt.show()






