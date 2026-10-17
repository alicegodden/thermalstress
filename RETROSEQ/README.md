# RETROSEQ – non-reference TE insertions under thermal stress

Scripts for calling, filtering and analysing non-reference transposable element (TE) insertions in zebrafish ovaries and testes after heat exposure, from:

> Godden AM, De Coriolis J, Sullivan R, Drake A, Immler S. *Sex differences in sRNA mediated transposon regulation under thermal stress in zebrafish germ cells.* (in review, *Genome Biology and Evolution*)

**Design.** Adult AB zebrafish, 2 weeks at 28 °C (control) or 34 °C (heat); 6 fish per group × 2 tissues (ovaries, testes) = 24 fish. Whole-genome sequencing (NovaSeq, 150 bp PE), aligned to GRCz11. Raw data: ENA PRJEB72689.

**Main question.** Does heat increase new TE insertions in the germline, and does it change which TEs insert or where they land?

---

## Pipeline overview

```
per-fish BAMs
   │  RetroSeq -discover / -call                        (HPC)
   ▼
per-fish RetroSeq VCFs
   │  filter: FL = 8, GQ ≥ 100, chr 1–25
   │  bedtools window -v -w 100  (remove calls ≤100 bp from annotated GRCz11 TEs)
   ▼
filtered per-fish VCFs   →  vcf/{ovaries,testes}/{control,treatment}/*.vcf
   │
   │  private_TE_per_fish_interim.py
   ▼
private insertions (in ONE fish only, within tissue)   →  vcf/private_uniq/
   │
   ├── te_map_unique_per_fish.py       Fig 3A–B  chromosome maps
   ├── fig3c_private_TE_per_fish.py    Fig 3C    per-fish, per-family stats (depth-adjusted)
   ├── te_open_chromatin_enrichment.py Table     open-chromatin enrichment
   ├── active_TE_class_landscape_private.R       class composition vs genome; reweighted landscape
   └── TE_age_transcribed_vs_mobilized.R         age of transcribed vs mobilized TEs
```

**Definition of a private insertion.** A non-reference insertion called in one fish and absent from every other fish **of the same tissue** (either temperature), i.e. no call within 100 bp. Ovaries and testes are processed independently. Nearby calls within the same fish are kept (merge them with `--collapse`).

---

## 1. Filtering per-fish VCFs (command line)

```bash
TE_BED=/path/to/GRCz11_tes_combined_sorted.bed

for vcf_file in *.vcf; do
    [[ "$vcf_file" == *_win100.vcf ]] && continue
    output_file="${vcf_file%.vcf}_win100.vcf"
    bedtools window -a "$vcf_file" -b "$TE_BED" -v -w 100 -header > "$output_file"
done
```

Apply the same steps to control and heat VCFs.

## 2. Folder layout expected by the scripts

```
vcf/
├── ovaries/
│   ├── control/     # 28 °C, one VCF per fish (e.g. ..._1FC_...vcf)
│   └── treatment/   # 34 °C (e.g. ..._2FT_...vcf)
└── testes/
    ├── control/     # ..._3MC_...vcf
    └── treatment/   # ..._4MT_...vcf
```

Tissue and temperature are taken from the folder names; fish IDs from the `_<n>[F|M][C|T]_` tag in the file name.

## 3. Scripts

| Script | Purpose | Main outputs |
| --- | --- | --- |
| `private_TE_per_fish_interim.py` | Keeps private insertions only; writes a private-only VCF per fish and summary tables | `private_uniq/<tissue>/<group>/<fish>.private.vcf`, `private_TE_per_fish.tsv`, `private_TE_by_treatment.tsv` |
| `te_map_unique_per_fish.py` | Mirrored chromosome maps per tissue: 28 °C left of each chromosome, 34 °C right, coloured by TE class, open-chromatin ticks, per-chromosome counts (6 fish summed) | `TE_map_unique_perfish_<tissue>.png/.pdf/.svg`, `_insertions.tsv` |
| `fig3c_private_TE_per_fish.py` | Private insertions per fish, all TEs and per TE class; 34 °C vs 28 °C within each tissue | `Fig3C_private_TE.png/.pdf/.svg`, `Fig3C_private_TE_stats.tsv`, `Fig3C_private_TE_counts.tsv` |
| `te_open_chromatin_enrichment.py` | Enrichment of private insertions near open chromatin (DANIO-CODE DOPEs) at 0 / 1 / 5 kb; heat vs control | `TE_open_chromatin_enrichment.tsv`, `TE_open_chromatin_heat_vs_ctrl.tsv` |
| `active_TE_class_landscape_private.R` | Class composition of private insertions vs genome-wide TE content (28 °C and 34 °C); heat vs control composition; Kimura landscape reweighted by the 34 °C private set | `fig_private_TE_enrichment`, `fig_private_TE_landscape`, `private_TE_class_*.csv` |
| `TE_age_transcribed_vs_mobilized.R` | Kimura age distribution of DE (transcribed) TEs vs private (mobilized) TEs, by class composition | `fig_TE_transcribed_vs_mobilized_age.png/.pdf`, `TE_transcribed_vs_mobilized_age_stats.csv` |
| `reads.txt` | Total reads per fish (samtools flagstat), used to adjust for sequencing depth | – |

### Running

```bash
python private_TE_per_fish_interim.py --root vcf/
python te_map_unique_per_fish.py      --root vcf/ --dopes daniocode_hub_280355_dopes_all.txt
python fig3c_private_TE_per_fish.py   --root vcf/ --reads reads.txt
python te_open_chromatin_enrichment.py --private-root vcf/private_uniq --dopes daniocode_hub_280355_dopes_all.txt
Rscript active_TE_class_landscape_private.R
Rscript TE_age_transcribed_vs_mobilized.R
```

Common options: `--window` (bp to call two insertions the same, default 100), `--fl` (RetroSeq FL filter, default 8), `--collapse` (merge nearby calls within a fish), `--mode temperature` (map: unique to temperature rather than to one fish).

---

## Statistics

- **Insertion counts (Fig 3C).** Negative binomial GLM per tissue, `private ~ temperature + offset(log(total reads))`, likelihood-ratio test (n = 6 vs 6); rate ratio 34 vs 28 °C with 95 % CI. Benjamini–Hochberg correction across TE classes and tissues. Sequencing depth is included because private-insertion counts track depth closely (Spearman ρ ≈ 0.8).
- **Open chromatin.** Exact binomial test and permutation test (10,000 random placements per chromosome) against the fraction of chr 1–25 covered by open chromatin; heat vs control by per-fish Wilcoxon rank-sum test and Fisher's exact test; BH-adjusted.
- **Class composition.** χ² tests of private insertions vs genome-wide TE content, and 34 °C vs 28 °C.

## Summary of results

| Tissue | Private insertions 28 °C | 34 °C | Depth-adjusted rate ratio (95 % CI) | p |
| --- | --- | --- | --- | --- |
| Ovaries | 1,252 | 1,271 | 1.05 (0.64–1.73) | 0.85 |
| Testes | 1,046 | 1,638 | 1.34 (0.87–2.07) | 0.20 |

Private insertions were enriched near open chromatin in all groups (≈1.5–2.3-fold within 1 kb), driven by DNA transposons; heat did not change insertion location or TE class composition in either tissue.

---

## Requirements

- RetroSeq v1.5, samtools, bedtools (HPC steps)
- Python ≥ 3.9: `numpy`, `pandas`, `scipy`, `matplotlib`
- R ≥ 4.1: `dplyr`, `tidyr`, `ggplot2`, `reshape2`
- Reference: GRCz11; TE annotation BED from RepeatMasker; open chromatin: DANIO-CODE DOPEs track (`daniocode_hub_280355_dopes_all.txt`)

## Contact

Dr Alice M. Godden – University of East Anglia
