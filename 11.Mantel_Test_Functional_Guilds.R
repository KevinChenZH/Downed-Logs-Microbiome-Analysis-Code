# Mantel test between functional guilds and wood physicochemical properties
# Required input files:
#   1. Wood chemistry CSV with columns matching: sample index, DM, Ash, TN, TP, TK, Cellulose, Lignin, Hemicellulose
#   2. Eight subcommunity CSV files: subcomm_B1.csv through subcomm_B4.csv, subcomm_F1.csv through subcomm_F4.csv
#      Each with ASVs as rows, samples as columns, optional Taxonomy column to drop

suppressPackageStartupMessages({
  library(linkET)
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(writexl)
})

# --- Configuration ---
chem_file   <- "WOOD_CHEMISTRY_CSV_PATH"
subcomm_dir <- "SUBCOMMUNITY_CSV_DIRECTORY"
out_dir     <- "OUTPUT_DIRECTORY"
dir.create(out_dir, showWarnings = FALSE)

# --- Read and standardize chemistry data ---
chem_raw <- read.csv(chem_file, fileEncoding = "GB18030",
                     stringsAsFactors = FALSE, check.names = FALSE)

id_map <- c("n1.1","n1.2","n1.3","n1.4","n3.1","n3.2","n3.3","n3.4",
            "n5.1","n5.2","n5.3","n5.4","w1.1","w1.2","w1.3","w1.4",
            "w3.1","w3.2","w3.3","w3.4","w5.1","w5.2","w5.3","w5.4")

chem_df <- data.frame(
  SampleID = id_map[as.integer(chem_raw[[2]])],
  DM = as.numeric(chem_raw[[4]]),
  Ash = as.numeric(chem_raw[[5]]),
  TN = as.numeric(chem_raw[[6]]),
  TP = as.numeric(chem_raw[[8]]),
  TK = as.numeric(chem_raw[[9]]),
  Cellulose = as.numeric(chem_raw[[10]]),
  Lignin = as.numeric(chem_raw[[11]]),
  Hemicellulose = as.numeric(chem_raw[[12]])
)

chem_scaled <- chem_df
chem_scaled[, 2:9] <- scale(chem_scaled[, 2:9])
rownames(chem_scaled) <- chem_scaled$SampleID

chem_matrix <- as.matrix(chem_scaled[, -1])
chem_cor <- correlate(chem_matrix, method = "pearson")

# --- Mantel test runner ---
run_mantel <- function(sc_name) {
  df <- read.csv(file.path(subcomm_dir, paste0("subcomm_", sc_name, ".csv")),
                 row.names = 1, stringsAsFactors = FALSE)
  
  if("Taxonomy" %in% colnames(df)) df <- df[, !colnames(df) %in% "Taxonomy"]
  df <- df[, id_map]
  
  comm_mat <- t(as.matrix(df))
  valid_samples <- rownames(comm_mat)[rowSums(comm_mat) > 0]
  
  if(length(valid_samples) < 4) return(NULL)
  
  comm_clean <- comm_mat[valid_samples, colSums(comm_mat[valid_samples, ]) > 0, drop = FALSE]
  env_sub <- chem_scaled[valid_samples, c("DM","Ash","TN","TP","TK","Cellulose","Lignin","Hemicellulose")]
  
  comm_dist <- vegdist(comm_clean, method = "bray")
  
  results <- data.frame()
  for(factor in colnames(env_sub)) {
    env_dist <- dist(env_sub[[factor]])
    mt <- mantel(comm_dist, env_dist, method = "pearson", permutations = 999)
    
    results <- rbind(results, data.frame(
      spec = sc_name,
      env = factor,
      r = mt$statistic,
      p = mt$signif,
      n_sample = length(valid_samples)
    ))
  }
  return(results)
}

bacteria_groups <- c("B1","B2","B3","B4")
fungi_groups <- c("F1","F2","F3","F4")

bacteria_results <- do.call(rbind, lapply(bacteria_groups, run_mantel))
fungi_results <- do.call(rbind, lapply(fungi_groups, run_mantel))

# --- Format for plotting ---
format_mantel <- function(df) {
  df %>%
    mutate(
      rd = cut(abs(r), breaks = c(-Inf, 0.3, 0.6, Inf),
               labels = c("< 0.3", "0.3 - 0.6", ">= 0.6")),
      pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
               labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")),
      line_color = case_when(
        p < 0.01 ~ "#08519C",
        p < 0.05 ~ "#6BAED6",
        TRUE ~ "#C6DBEF"
      ),
      line_size = case_when(
        abs(r) >= 0.6 ~ 1.8,
        abs(r) >= 0.3 ~ 1.0,
        TRUE ~ 0.4
      )
    )
}

bacteria_plot <- format_mantel(bacteria_results) %>%
  mutate(spec = factor(spec, levels = bacteria_groups))
fungi_plot <- format_mantel(fungi_results) %>%
  mutate(spec = factor(spec, levels = fungi_groups))

# --- Visualization ---
create_mantel_plot <- function(mantel_data, group_label, title_text) {
  
  color_pal <- c("< 0.01" = "#08519C", "0.01 - 0.05" = "#6BAED6", ">= 0.05" = "#C6DBEF")
  size_vals <- c("< 0.3" = 0.4, "0.3 - 0.6" = 1.0, ">= 0.6" = 1.8)
  
  p <- qcorrplot(chem_cor, type = "lower", diag = FALSE) +
    geom_square(aes(fill = r), colour = "white", size = 0.3) +
    scale_fill_gradient2(name = "Pearson's r",
                         low = "#74A9CF",
                         mid = "#F7F7F7",
                         high = "#C94641",
                         midpoint = 0, limits = c(-1, 1), na.value = "white") +
    geom_couple(aes(colour = pd, size = rd),
                data = mantel_data,
                curvature = 0.08,
                lineend = "round") +
    scale_colour_manual(name = "Mantel's p", values = color_pal) +
    scale_size_manual(name = "Mantel's r", values = size_vals) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "grey70"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
      axis.text.y = element_text(size = 11, face = "bold"),
      legend.position = "right",
      legend.box = "vertical",
      legend.margin = margin(0, 0, 0, 10),
      legend.key = element_rect(fill = NA),
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5, color = "#1a1a1a")
    ) +
    labs(title = title_text)
  
  return(p)
}

p_bacteria <- create_mantel_plot(bacteria_plot, "Bacterial Guilds (B1-B4)", "A  Bacterial Functional Groups")
p_fungi <- create_mantel_plot(fungi_plot, "Fungal Guilds (F1-F4)", "B  Fungal Functional Groups")

ggsave(file.path(out_dir, "FigA_Bacterial_Mantel.pdf"), p_bacteria,
       width = 10, height = 8, device = cairo_pdf)
ggsave(file.path(out_dir, "FigB_Fungal_Mantel.pdf"), p_fungi,
       width = 10, height = 8, device = cairo_pdf)

ggsave(file.path(out_dir, "FigA_Bacterial_Mantel.png"), p_bacteria,
       width = 10, height = 8, dpi = 600)
ggsave(file.path(out_dir, "FigB_Fungal_Mantel.png"), p_fungi,
       width = 10, height = 8, dpi = 600)

# --- Export results ---
write_xlsx(list(
  "Bacteria_B1-B4" = bacteria_plot %>%
    select(spec, env, r, p, n_sample, rd, pd) %>%
    arrange(desc(abs(r))),
  "Fungi_F1-F4" = fungi_plot %>%
    select(spec, env, r, p, n_sample, rd, pd) %>%
    arrange(desc(abs(r))),
  "Summary_Statistics" = data.frame(
    Group = c(rep("Bacteria", 4), rep("Fungi", 4)),
    Guild = c(bacteria_groups, fungi_groups),
    Significant_p0.01 = c(
      sum(bacteria_plot$p < 0.01 & bacteria_plot$spec == "B1"),
      sum(bacteria_plot$p < 0.01 & bacteria_plot$spec == "B2"),
      sum(bacteria_plot$p < 0.01 & bacteria_plot$spec == "B3"),
      sum(bacteria_plot$p < 0.01 & bacteria_plot$spec == "B4"),
      sum(fungi_plot$p < 0.01 & fungi_plot$spec == "F1"),
      sum(fungi_plot$p < 0.01 & fungi_plot$spec == "F2"),
      sum(fungi_plot$p < 0.01 & fungi_plot$spec == "F3"),
      sum(fungi_plot$p < 0.01 & fungi_plot$spec == "F4")
    ),
    Significant_p0.05 = c(
      sum(bacteria_plot$p < 0.05 & bacteria_plot$spec == "B1"),
      sum(bacteria_plot$p < 0.05 & bacteria_plot$spec == "B2"),
      sum(bacteria_plot$p < 0.05 & bacteria_plot$spec == "B3"),
      sum(bacteria_plot$p < 0.05 & bacteria_plot$spec == "B4"),
      sum(fungi_plot$p < 0.05 & fungi_plot$spec == "F1"),
      sum(fungi_plot$p < 0.05 & fungi_plot$spec == "F2"),
      sum(fungi_plot$p < 0.05 & fungi_plot$spec == "F3"),
      sum(fungi_plot$p < 0.05 & fungi_plot$spec == "F4")
    ),
    Max_r = c(
      max(abs(bacteria_plot$r[bacteria_plot$spec == "B1"])),
      max(abs(bacteria_plot$r[bacteria_plot$spec == "B2"])),
      max(abs(bacteria_plot$r[bacteria_plot$spec == "B3"])),
      max(abs(bacteria_plot$r[bacteria_plot$spec == "B4"])),
      max(abs(fungi_plot$r[fungi_plot$spec == "F1"])),
      max(abs(fungi_plot$r[fungi_plot$spec == "F2"])),
      max(abs(fungi_plot$r[fungi_plot$spec == "F3"])),
      max(abs(fungi_plot$r[fungi_plot$spec == "F4"]))
    )
  )
), path = file.path(out_dir, "L1_Mantel_Results_Writing.xlsx"))