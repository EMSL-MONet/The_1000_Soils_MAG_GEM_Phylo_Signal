# Phylogenetic signal in predicted metabolite exchange of soil MAGs

R code and input data for the phylogenetic analyses of flux balance analysis
(FBA) predictions from genome-scale metabolic models (GEMs) of soil
metagenome-assembled genomes (MAGs) recovered from 94 soil samples across the
continental United States.

Everything runs from one R Markdown notebook, `analysis.Rmd`, which

1. **reproduces the published results**: the Mantel tests of phylogenetic
   signal, Figure 1, Figure 2, Supplementary Figure 1, Supplementary Table 2,
   and the statistics quoted in the main text, and then checks each one
   against the value reported in the manuscript; and
2. **lets you explore the data**: what a given MAG takes up and secretes, and
   how a metabolite's exchange varies across phyla.

## Quick start

1. Install the packages listed under [Requirements](#requirements).
2. Clone or download this repository.
3. Open `analysis.Rmd` in RStudio and click **Knit**, or from the repository
   root run:

   ```r
   rmarkdown::render("analysis.Rmd")
   ```

The first run takes about 30 minutes, almost all of it the bacterial Mantel
tests. Their results are cached (in `analysis_cache/`), so later runs take a few
minutes. The rendered report is `analysis.html`; figures and tables are written
to `results/`.

To explore a different MAG or metabolite, change `example_mag` and
`example_compound` at the top of `analysis.Rmd`, or pass them when rendering:

```r
rmarkdown::render("analysis.Rmd",
                  params = list(example_mag = "KONA-TOP_bin.2.s",
                                example_compound = "Ethanol"))
```

## Repository contents

```
The_1000_Soils_MAG_GEM_Phylo_Signal/
├── analysis.Rmd              # the notebook: reproduce + explore
├── R/
│   ├── config.R              # input/output paths and analysis parameters
│   ├── utils.R               # shared helpers (tree distances, flux normalization)
│   ├── analysis.R            # functions that reproduce the published analyses
│   └── explore.R             # functions for exploring the flux data
├── data/                     # inputs (see below)
└── results/                  # created when the notebook runs (not tracked in git)
```

## Input data (`data/`)

| File | Description |
|---|---|
| `fba_uptake_bacteria.tsv` | FBA exchange fluxes for the 498 bacterial MAGs (254 compounds). One row per compound, one column per MAG headed `<phylum>\|\|<MAG name>`. Positive = uptake, negative = secretion, NA = compound absent from that MAG's model. |
| `fba_uptake_archaea.tsv` | Same, for the 53 archaeal MAGs (141 compounds). |
| `bac120.contree` | IQ-TREE ultrafast-bootstrap consensus tree of the bacterial MAGs and 88 GTDB reference genomes (bac120 alignment). |
| `ar53.contree` | Same, for the archaeal MAGs and 18 GTDB reference genomes (ar53 alignment). |
| `checkm2_551_MAGs.tsv` | CheckM2 quality report (completeness, contamination, genome size, GC content) for the 551 dereplicated MAGs. |
| `published_supplementary_table_2.tsv` | Supplementary Table 2 as published, used by the notebook to check its reproduction. |

The FBA tables were compiled from the per-MAG KBase FBA outputs (simulations
under default complete media). One bacterial MAG (PPRH-TOP_bin.5.o) has no FBA
output: a draft model was reconstructed with the Gram-negative template, but
gapfilling did not yield a solution supporting biomass production and FBA
predicted no growth. It appears as an all-NA column and is removed at the start
of the analysis.

The per-MAG model and FBA files for 550 MAGs, and the genome sequences of all
551 MAGs, are on Zenodo: https://doi.org/10.5281/zenodo.23480528

Upstream steps (MAG recovery with metaGEM, dereplication with galah, quality
assessment with CheckM2, classification with GTDB-Tk against GTDB R10-RS226,
phylogenetic inference with IQ-TREE v3.1.3, and model reconstruction and FBA in
KBase) are described in the manuscript's Supplementary Methods and are not part
of this code.

## Outputs (`results/`)

| File | Manuscript |
|---|---|
| `mantel_bacteria.tsv`, `mantel_archaea.tsv` | Mantel r, raw and BH-adjusted p-values for every compound |
| `Supplementary_Table_2_bacteria_mantel_r_Padj.tsv` | Supplementary Table 2 |
| `figures/Figure_1_bacteria_tree_with_rings.pdf` | Figure 1 |
| `figures/Figure_2_mantel_scatter_significant.pdf` | Figure 2 |
| `figures/Supplementary_Figure_1_mantel_scatter_other.pdf` | Supplementary Figure 1 |
| `gc_content_exchange_pattern_test.tsv` | GC-content comparison (main text, P = 0.33) |
| `fba_filtered_bacteria.tsv`, `fba_filtered_archaea.tsv`, `tip_metadata_bacteria.tsv` | Intermediate tables |

Figures are raw ggplot/ggtree output. The published versions were finished in
Adobe Illustrator (legend layout, ring numbering, label formatting); data,
colors, and ordering are unchanged.

## Reproducibility notes

- **Permutation p-values depend on input order.** The random seed is set once
  before looping over compounds, so a compound's p-value depends on its position
  in the filtered table and on the order of MAG columns. Compounds are sorted by
  mean flux (ties kept in compound-ID order) and MAG columns keep the order of
  the input table; reordering either changes permutation p-values slightly
  (Mantel r values are unaffected).
- **Duplicate tests.** Eight trace ions (Mg, Co²⁺, K⁺, Mn²⁺, Cl⁻, Ca²⁺, Cu²⁺,
  Zn²⁺) have identical flux profiles, and enterobactin and Fe-enterochelin have
  identical distance matrices. All are kept as separate tests, which makes the
  BH correction conservative.
- **Untested compounds.** Compounds with zero variance across MAGs (six for
  Bacteria: Cd²⁺, Hg²⁺, Na⁺, Ni²⁺, Pb, glycogen; nine for Archaea) have no Mantel
  r and are excluded from the BH correction.

## Requirements

R (≥ 4.5) with:

- CRAN: `rmarkdown`, `knitr`, `ape`, `vegan`, `dplyr`, `readr`, `tidyr`,
  `tibble`, `ggplot2`, `ggnewscale`, `scales`
- Bioconductor: `treeio`, `ggtree`, `ggtreeExtra`

```r
install.packages(c("rmarkdown", "knitr", "ape", "vegan", "dplyr", "readr",
                   "tidyr", "tibble", "ggplot2", "ggnewscale", "scales"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("treeio", "ggtree", "ggtreeExtra"))
```

Versions used for the manuscript: R 4.5.1, ape 5.8.1, vegan 2.7.5, treeio 1.32.0,
ggtree 3.16.3, ggtreeExtra 1.18.1, ggplot2 4.0.3, ggnewscale 0.5.2, dplyr 1.2.1,
readr 2.2.0, rmarkdown 2.31, knitr 1.51. The rendered notebook ends with the
full `sessionInfo()`.

## Citation

If you use this code, please cite: [manuscript citation].

## License

[License to be added.]
