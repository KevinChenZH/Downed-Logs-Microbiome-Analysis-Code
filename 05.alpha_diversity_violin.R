# alpha_diversity_violin.R
# Required: readxl, dplyr, tidyr, ggplot2, patchwork
# Input: totalgroup.xlsx with 16SSample_Name, ITSSample_Name, and metric columns
# Output: 600dpi PDF and PNG

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

COLORS <- c(Internal = "#5DADE2", External = "#F39C12")
POINT_COLORS <- c(Internal = "#AED6F1", External = "#FAD7A0")

load_data <- function(file_path) {
  df <- read_excel(file_path)
  
  df_16s <- data.frame(
    Sample_ID = df[["16SSample_Name"]],
    chao1 = as.numeric(df[["chao1"]]),
    observed_features = as.numeric(df[["observed_features"]]),
    pielou_e = as.numeric(df[["pielou_e"]]),
    shannon = as.numeric(df[["shannon"]]),
    Group = "Bacteria",
    stringsAsFactors = FALSE
  )
  
  its_map <- c("chao1", "observed_features", "pielou_e", "shannon")
  its_df <- data.frame(Sample_ID = df[["ITSSample_Name"]], stringsAsFactors = FALSE)
  for (m in its_map) {
    cand <- c(paste0(m, ".1"), paste0(m, "_1"), m)
    col <- intersect(cand, names(df))[1]
    if (is.na(col)) col <- names(df)[which(names(df) == m)[1] + 4]
    its_df[[m]] <- as.numeric(df[[col]])
  }
  its_df$Group <- "Fungi"
  
  bind_rows(df_16s, its_df) %>%
    filter(!is.na(Sample_ID)) %>%
    mutate(
      Position = ifelse(substr(Sample_ID, 1, 1) == "n", "Internal", "External"),
      Stage = case_when(
        substr(Sample_ID, 2, 2) == "1" ~ "I",
        substr(Sample_ID, 2, 2) == "3" ~ "III",
        substr(Sample_ID, 2, 2) == "5" ~ "V"
      ),
      Stage = factor(Stage, levels = c("I", "III", "V")),
      Position = factor(Position, levels = c("Internal", "External")),
      Group = factor(Group, levels = c("Bacteria", "Fungi"))
    )
}

make_panel <- function(df, grp, metric, title, ylab, show_legend = FALSE) {
  sub <- df %>% filter(Group == grp)
  
  p <- ggplot(sub, aes(x = Stage, y = .data[[metric]], fill = Position)) +
    geom_violin(position = position_dodge(width = 0.7), alpha = 0.75, linewidth = 0.8, width = 0.7, trim = TRUE) +
    geom_boxplot(position = position_dodge(width = 0.7), width = 0.15, outlier.shape = NA, alpha = 0.5, linewidth = 0.5) +
    geom_point(aes(color = Position), position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7),
               size = 2, alpha = 0.7, show.legend = FALSE) +
    scale_fill_manual(values = COLORS, name = "") +
    scale_color_manual(values = POINT_COLORS) +
    labs(title = paste(grp, "-", title), x = "Decay Level", y = ylab) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major.y = element_line(color = "grey90", linetype = "dashed"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
      panel.border = element_rect(linewidth = 1),
      axis.line = element_line(linewidth = 0.8)
    )
  
  if (show_legend) {
    p + theme(
      legend.position = "inside",
      legend.position.inside = c(0.02, 0.98),
      legend.justification = c("left", "top"),
      legend.direction = "horizontal",
      legend.background = element_rect(fill = "white", color = NA),
      legend.key.size = unit(0.4, "cm"),
      legend.text = element_text(size = 10),
      legend.title = element_blank()
    )
  } else {
    p + theme(legend.position = "none")
  }
}

main <- function() {
  df <- load_data("path/to/totalgroup.xlsx")
  
  p1 <- make_panel(df, "Bacteria", "shannon", "Shannon Diversity", "Shannon Diversity", show_legend = TRUE)
  p2 <- make_panel(df, "Bacteria", "observed_features", "Observed ASVs", "Observed ASVs")
  p3 <- make_panel(df, "Fungi", "shannon", "Shannon Diversity", "Shannon Diversity")
  p4 <- make_panel(df, "Fungi", "observed_features", "Observed ASVs", "Observed ASVs")
  
  combined <- (p1 + p2) / (p3 + p4) +
    plot_annotation(
      title = "Alpha Diversity Patterns in Decaying Wood",
      theme = theme(
        plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        legend.position = "none"
      )
    )
  
  ggsave("path/to/output/Alpha_Diversity_CNS_Integrated.pdf", combined,
         width = 10, height = 8, dpi = 600)
  ggsave("path/to/output/Alpha_Diversity_CNS_Integrated.png", combined,
         width = 10, height = 8, dpi = 600, bg = "white")
}

main()