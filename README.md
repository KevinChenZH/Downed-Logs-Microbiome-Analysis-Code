# Downed-Logs-Microbiome-Analysis-Code

**Project Title**  
Internal–external heterogeneity drives bacterial-fungal assembly in downed logs through environmental co-filtering rather than metabolic coupling

**Authors**  
Zhihao Chen, Yanbing Lin  
Northwest A&F University, College of Life Sciences

## System Requirements

- R version 4.5.3
- RStudio 2026.01.2

## Key R Packages

`phyloseq`, `vegan`, `ALDEx2` (v1.38.0), `plspm`, `SpiecEasi`, `igraph`, `ggplot2`, `dplyr`, `agricolae`, `rcompanion`, `linkET` (v0.0.3; Huang 2021), `ggcor`

## File Inventory & Execution Order

| No. | File Name | Purpose | Input | Output |
|:---|:---|:---|:---|:---|
| 01 | `physicochemical_profiles.R` | Wood chemistry summary and partial eta-squared (η²) effect-size visualization | `physicochemistry.csv` | Table S1, Fig. S1 |
| 02 | `alpha_tukey_posthoc.R` | Tukey HSD and compact letter display (CLD) for alpha diversity | `alpha_diversity.csv` | Table S3 |
| 03 | `alpha_two_way_anova.R` | Two-way ANOVA (Type II/III with robust standard errors when Levene's test fails) for alpha diversity | `alpha_diversity.csv` | Table S5 |
| 04 | `alpha_kw_test.R` | Kruskal–Wallis non-parametric validation of alpha diversity | `alpha_diversity.csv` | Table S4 |
| 05 | `alpha_diversity_violin.R` | Violin plots of Shannon index and observed ASVs | `alpha_diversity.csv` | Fig. 1b |
| 06 | `beta_diversity_pcoa_permanova.R` | PCoA and global PERMANOVA (999 permutations, weighted UniFrac) | `ASV_table.biom`, `tree.nwk`, `metadata.tsv` | Fig. 1c |
| 07 | `beta_pairwise_permanova.R` | Pairwise PERMANOVA (simple effects) | `ASV_table.biom`, `metadata.tsv` | Table S6 |
| 08 | `Stacked_Bar_Community_Composition.R` | Stacked bar plots of bacterial genus-level and fungal phylum-level composition | `taxa_summary.csv` | Fig. S2 |
| 09 | `ALDEx2_Scale_Model_Analysis.R` | Differential abundance using ALDEx2 Scale Model (γ = 0.5, 1,000 Monte Carlo samples) | `ASV_table.csv`, `metadata.tsv` | Table S7, Fig. S3, Fig. 3 |
| 10 | `Cross_Kingdom_Cooccurrence_Networks.R` | Cross-kingdom co-occurrence networks via SpiecEasi (mb, StARS λ = 0.1, 100 subsamples) | `bacteria_ASV.csv`, `fungi_ASV.csv` | Table S10a, S10b, Fig. 4 |
| 11 | `Mantel_Test_Functional_Guilds.R` | Mantel test (999 permutations) between guild dissimilarity and wood-chemistry distance matrices; visualization via `linkET` and `ggcor` | `guild_abundance.csv`, `physicochemistry.csv` | Table S8, Fig. 5a,b |
| 12 | `PLS_SEM_Models_M1_M3.R` | Partial least squares path modeling (1,000 bootstrap resamples) for bacterial B1 guild (M1) and fungal F3 guild (M3) | `sem_dataset.csv` | Table S9a, S9b, S9c, Fig. 5c–f |

## Execution Notes

- Scripts 01–04 are independent and can be run in any order.
- Scripts 05–08 require processed outputs from the QIIME2/DADA2 pipeline (ASV tables, taxonomy assignments, and phylogenetic tree).
- Scripts 09–12 should be executed sequentially because downstream analyses depend on upstream guild assignments and differential-abundance results.
- All scripts use **relative paths** (`./data/`, `./output/`). Update the root directory variable at the top of each script to match your local environment before running.
- For reproducibility, package versions are pinned as listed above. `linkET` (v0.0.3) is specifically required for the Mantel-test network heatmap visualization in Script 11.

## Data Availability

Raw 16S rRNA (V5–V7) and ITS1 amplicon sequencing data are deposited in NCBI SRA under BioProject PRJNAxxxxxx (to be updated upon acceptance). Processed ASV tables, metadata, and physicochemical data are available upon reasonable request.

## References

Huang, H. (2021). *linkET: Everything is Linkable*. R package version 0.0.3.
