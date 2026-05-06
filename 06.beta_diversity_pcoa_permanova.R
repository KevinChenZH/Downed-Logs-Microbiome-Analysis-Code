# beta_diversity_pcoa_permanova.R
# Required: vegan, ape, ggplot2, patchwork, dplyr, readr
# Input: metadata.csv (SampleID, Stage, Position), 16S/ITS distance-matrix.tsv (sample x sample)
# Output: Figure2 PDF/PNG (600dpi), Table2_PERMANOVA.csv

library(vegan)
library(ape)
library(ggplot2)
library(patchwork)
library(dplyr)
library(readr)

STAGE_COLORS <- c(I = "#009E73", III = "#0072B2", V = "#D55E00")
POS_SHAPES <- c(Internal = 16, External = 17)

load_dm <- function(path) {
  df <- read_tsv(path, col_names = TRUE, show_col_types = FALSE)
  mat <- as.matrix(df[, -1])
  rownames(mat) <- df[[1]]
  class(mat) <- "numeric"
  as.dist(mat)
}

run_perm <- function(dm, meta) {
  meta <- meta[match(labels(dm), meta$SampleID), ]
  f <- adonis2(dm ~ Stage * Position, data = meta, permutations = 999, by = "margin")
  data.frame(
    Factor = rownames(f)[1:3],
    R2 = round(f$R2[1:3], 3),
    F = round(f$F[1:3], 2),
    p = ifelse(f$`Pr(>F)`[1:3] < 0.001, "<0.001", format(round(f$`Pr(>F)`[1:3], 4), nsmall = 4))
  )
}

plot_pcoa <- function(dm, meta, title, perm_df, panel_label) {
  meta <- meta[match(labels(dm), meta$SampleID), ]
  ord <- pcoa(dm)
  eig <- ord$values$Relative_eig[1:2] * 100
  df <- data.frame(PC1 = ord$vectors[, 1], PC2 = ord$vectors[, 2], SampleID = labels(dm))
  df <- left_join(df, meta, by = "SampleID")
  
  st <- paste0(
    "Stage: R2=", perm_df$R2[perm_df$Factor == "Stage"], ", p=", perm_df$p[perm_df$Factor == "Stage"], "\n",
    "Position: R2=", perm_df$R2[perm_df$Factor == "Position"], ", p=", perm_df$p[perm_df$Factor == "Position"], "\n",
    "Interaction: R2=", perm_df$R2[perm_df$Factor == "Stage:Position"], ", p=", perm_df$p[perm_df$Factor == "Stage:Position"]
  )
  
  xr <- diff(range(df$PC1, na.rm = TRUE))
  yr <- diff(range(df$PC2, na.rm = TRUE))
  
  ggplot(df, aes(x = PC1, y = PC2, color = Stage, shape = Position)) +
    stat_ellipse(aes(group = Stage, fill = Stage), geom = "polygon", level = 0.95,
                 type = "norm", alpha = 0.25, color = NA) +
    geom_point(size = 3, alpha = 0.9, stroke = 0.8) +
    scale_color_manual(values = STAGE_COLORS) +
    scale_fill_manual(values = STAGE_COLORS, guide = "none") +
    scale_shape_manual(values = POS_SHAPES) +
    labs(title = paste(panel_label, title),
         x = sprintf("PC1 (%.1f%%)", eig[1]),
         y = sprintf("PC2 (%.1f%%)", eig[2])) +
    annotate("label", x = min(df$PC1) + xr * 0.02, y = max(df$PC2) - yr * 0.02,
             label = st, hjust = 0, vjust = 1, size = 3, label.size = 0.4,
             fill = "white", color = "black", alpha = 0.9, parse = FALSE) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major = element_line(color = "grey90", linetype = "dashed"),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      panel.border = element_rect(linewidth = 1),
      axis.line = element_line(linewidth = 0.8)
    )
}

main <- function() {
  meta <- read_csv("path/to/metadata.csv", show_col_types = FALSE)
  dm_b <- load_dm("path/to/16S_distance_matrix.tsv")
  dm_f <- load_dm("path/to/ITS_distance_matrix.tsv")
  
  perm_b <- run_perm(dm_b, meta)
  perm_f <- run_perm(dm_f, meta)
  
  p_b <- plot_pcoa(dm_b, meta, "Bacteria", perm_b, "(A)")
  p_f <- plot_pcoa(dm_f, meta, "Fungi", perm_f, "(B)")
  
  combined <- (p_b + p_f) + plot_layout(ncol = 2, guides = "collect") &
    theme(legend.position = "bottom", legend.box = "horizontal", legend.title = element_blank())
  
  ggsave("path/to/output/Figure2_Beta_Diversity.pdf", combined,
         width = 10, height = 4.5, dpi = 600)
  ggsave("path/to/output/Figure2_Beta_Diversity.png", combined,
         width = 10, height = 4.5, dpi = 600, bg = "white")
  
  bind_rows(perm_b %>% mutate(Group = "Bacteria"), perm_f %>% mutate(Group = "Fungi")) %>%
    select(Group, everything()) %>%
    write_csv("path/to/output/Table2_PERMANOVA.csv")
}

main()