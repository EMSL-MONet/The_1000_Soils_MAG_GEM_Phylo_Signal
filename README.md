# Phylogenetic signal in predicted metabolite exchange of soil MAGs

R code for the phylogenetic analyses of flux balance analysis (FBA) predictions
from genome-scale metabolic models (GEMs) of soil metagenome-assembled genomes
(MAGs) recovered from 94 soil samples across the continental United States.

The pipeline starts from compiled FBA exchange fluxes and IQ-TREE phylogenies and
produces the Mantel test results, Figure 1, Figure 2, Supplementary Figure 1,
Supplementary Table 2, and the GC-content comparison reported in the manuscript.

## Folder structure

```
The_1000_Soils_MAG_GEM_Phylo_Signal/
├── run_all.R                 # runs every script in order
├── R/
│   ├── config.R              # input/output paths and parameters
│   └── utils.R               # shared helper functions
├── scripts/
│   ├── 01_filter_fba_tables.R
│   ├── 02_mantel_phylo_signal.R
│   ├── 03_build_tip_metadata.R
│   ├── 04_figure1_tree_with_rings.R
│   ├── 05_figure2_mantel_scatter.R
│   └── 06_gc_content_test.R
├── data/                     # inputs (see below)
└── results/                  # created when the scripts run (not tracked in git)
```

## Input data (`data/`)

| File | Description |
|---|---|
| `fba_uptake_bacteria.tsv` | FBA exchange fluxes for the 498 bacterial MAGs (254 compounds). One row per compound, one column per MAG headed `<phylum>\|\|<MAG name>`. Positive = uptake, negative = secretion, NA = compound absent from that MAG's model. |
| `fba_uptake_archaea.tsv` | Same, for the 53 archaeal MAGs (141 compounds). |
| `bac120.contree` | IQ-TREE ultrafast-bootstrap consensus tree of the bacterial MAGs and 88 GTDB reference genomes (bac120 alignment). |
| `ar53.contree` | Same, for the archaeal MAGs and 18 GTDB reference genomes (ar53 alignment). |
| `checkm2_551_MAGs.tsv` | CheckM2 quality report (completeness, contamination, genome size, GC content) for the 551 dereplicated MAGs. |

The FBA tables were compiled from per-MAG KBase FBA outputs (simulations run under
default complete media) by joining each MAG's `uptake` column on compound ID and
labeling MAGs with their GTDB-Tk phylum. One bacterial MAG (PPRH-TOP_bin.5.o) has
no FBA output: a draft model was reconstructed with the Gram-negative template, but
gapfilling did not yield a solution supporting biomass production and FBA predicted
no growth. It appears as an all-NA column and is removed in step 1.

The per-MAG KBase model and FBA files (550 MAGs) are deposited separately on
Zenodo [DOI]; see the README there for details.

Upstream steps (MAG recovery with metaGEM, dereplication with galah, quality
assessment with CheckM2, classification with GTDB-Tk against GTDB R10-RS226,
phylogenetic inference with IQ-TREE v3.1.3, and model reconstruction and FBA in
KBase) are described in the Supplementary Methods and are not part of this code.

## Running

From the repository root:

```r
setwd("path/to/The_1000_Soils_MAG_GEM_Phylo_Signal")
source("run_all.R")
```

or `Rscript run_all.R` from a terminal. Total runtime is about 30 minutes, most of
it the bacterial Mantel tests (step 2). Each script can also be run on its own,
in order, from the repository root.

## Steps and outputs (`results/`)

| Script | What it does | Output | Manuscript |
|---|---|---|---|
| `01_filter_fba_tables.R` | Removes MAGs without FBA output; drops compounds missing in ≥ one-third of MAGs; sets remaining NAs to zero. Bacteria: 497 MAGs, 50 of 254 compounds. Archaea: 53 MAGs, 67 of 141 compounds. | `fba_filtered_bacteria.tsv`, `fba_filtered_archaea.tsv` | |
| `02_mantel_phylo_signal.R` | Per compound, Mantel test (Spearman, 999 permutations, seed 1) between pairwise patristic distance and pairwise difference in normalized flux; Benjamini–Hochberg correction within each domain. | `mantel_bacteria.tsv`, `mantel_archaea.tsv`, `Supplementary_Table_2_bacteria_mantel_r_Padj.tsv` | Supplementary Table 2 |
| `03_build_tip_metadata.R` | Phylum and MAG/reference status for each tip of the bacterial tree. | `tip_metadata_bacteria.tsv` | |
| `04_figure1_tree_with_rings.R` | Circular bacterial tree with GC content, genome size, and flux rings for compounds with P_adj < 0.05. | `figures/Figure_1_bacteria_tree_with_rings.pdf` | Figure 1 |
| `05_figure2_mantel_scatter.R` | Patristic distance vs. flux difference for every MAG pair, per compound, annotated with Mantel r and P_adj. | `figures/Figure_2_mantel_scatter_significant.pdf`, `figures/Supplementary_Figure_1_mantel_scatter_other.pdf` | Figure 2, Supplementary Figure 1 |
| `06_gc_content_test.R` | Wilcoxon rank-sum (Mann–Whitney U) test of GC content between MAGs combining BCAA uptake with pyruvate and ethanol secretion and all other bacterial MAGs. | `gc_content_exchange_pattern_test.tsv` | Main text (P = 0.33) |

`run_all.R` also writes `results/sessionInfo.txt` with the R and package versions used.

Figures are written as raw ggplot/ggtree output. The published versions were
finished in Adobe Illustrator (legend layout, ring numbering, label formatting);
data, colors, and ordering are unchanged.

## Reproducibility notes

- **Permutation p-values depend on input order.** The random seed is set once
  before looping over compounds, so the p-value for a compound depends on its
  position in the filtered table and on the order of MAG columns. Compounds are
  sorted by mean flux, with ties kept in compound-ID order, and MAG columns keep
  the order of the input table; reordering rows or columns will change
  permutation p-values slightly (Mantel r values are unaffected).
- **Duplicate tests.** Eight trace ions (Mg, Co²⁺, K⁺, Mn²⁺, Cl⁻, Ca²⁺, Cu²⁺, Zn²⁺)
  have identical flux profiles, and enterobactin and Fe-enterochelin have identical
  distance matrices. All are kept as separate tests, which makes the BH correction
  conservative.
- **Untested compounds.** Compounds with zero variance across MAGs (six for
  Bacteria: Cd²⁺, Hg²⁺, Na⁺, Ni²⁺, Pb, glycogen; nine for Archaea) have no Mantel
  r and are excluded from the BH correction.

## Requirements

R (≥ 4.5) with the following packages:

- CRAN: `ape`, `vegan`, `dplyr`, `readr`, `tidyr`, `tibble`, `ggplot2`, `ggnewscale`, `scales`
- Bioconductor: `treeio`, `ggtree`, `ggtreeExtra`

```r
install.packages(c("ape", "vegan", "dplyr", "readr", "tidyr", "tibble",
                   "ggplot2", "ggnewscale", "scales"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("treeio", "ggtree", "ggtreeExtra"))
```

Versions used for the manuscript: R 4.5.1, ape 5.8.1, vegan 2.7.5, treeio 1.32.0,
ggtree 3.16.3, ggtreeExtra 1.18.1, ggplot2 4.0.3, ggnewscale 0.5.2, dplyr 1.2.1,
readr 2.2.0.

## Citation

If you use this code, please cite: [manuscript citation].

## License

[License to be added.]
