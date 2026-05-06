# ALDEx2 Scale Model Analysis for Downed Log Microbiome
# Gamma = 0.5 | MC samples = 1000 | Genus-level aggregation
# Required input files:
#   1. Metadata CSV with columns: Stage (I/III/V), Position (Internal/External)
#   2. 16S feature table in BIOM format (.biom)
#   3. 16S taxonomy TSV with Silva-style taxonomic strings
#   4. ITS feature table in BIOM format (.biom)
#   5. ITS taxonomy TSV with UNITE-style taxonomic strings

suppressPackageStartupMessages({
  library(ALDEx2)
  library(biomformat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggrepel)
  library(cowplot)
  library(stringr)
  library(forcats)
})

# --- Configuration ---
output_dir <- "OUTPUT_DIRECTORY_PATH"
input_metadata <- "METADATA_CSV_PATH"
input_16s_biom <- "16S_FEATURE_TABLE_BIOM_PATH"
input_16s_tax <- "16S_TAXONOMY_TSV_PATH"
input_its_biom <- "ITS_FEATURE_TABLE_BIOM_PATH"
input_its_tax <- "ITS_TAXONOMY_TSV_PATH"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "Paired_Tests"), showWarnings = FALSE)
dir.create(file.path(output_dir, "Main_Figure"), showWarnings = FALSE)
setwd(output_dir)

# --- Utility functions ---
clean_sample_ids <- function(ids) {
  ids <- str_trim(as.character(ids))
  ids <- str_replace_all(ids, "\\s+", "_")
  ids <- str_replace_all(ids, "-", "_")
  ids <- make.names(ids, unique = TRUE)
  return(ids)
}

read_taxonomy_genus <- function(tax_path) {
  lines <- readLines(tax_path)
  data_lines <- lines[!grepl("^#", lines)]
  if(length(data_lines) < 2) stop("Taxonomy file has insufficient data")
  
  header <- strsplit(data_lines[1], "\t")[[1]]
  data <- do.call(rbind, strsplit(data_lines[-1], "\t"))
  tax_df <- as.data.frame(data, stringsAsFactors = FALSE)
  colnames(tax_df) <- header
  
  id_col <- grep("Feature ID|feature-id|#OTU ID|ASV_ID", colnames(tax_df), ignore.case = TRUE)[1]
  tax_col <- grep("Taxon|taxonomy", colnames(tax_df), ignore.case = TRUE)[1]
  colnames(tax_df)[id_col] <- "ASV_ID"
  colnames(tax_df)[tax_col] <- "Taxon"
  
  tax_df <- tax_df %>%
    mutate(
      Taxon_clean = str_remove_all(Taxon, "^[dpcofgs]__"),
      Taxon_clean = str_replace_all(Taxon_clean, ";[dpcofgs]__", ";"),
      Taxon_Split = str_split(Taxon_clean, ";"),
      Kingdom = map_chr(Taxon_Split, ~ ifelse(length(.x) >= 1, .x[1], NA_character_)),
      Phylum = map_chr(Taxon_Split, ~ ifelse(length(.x) >= 2, .x[2], NA_character_)),
      Class = map_chr(Taxon_Split, ~ ifelse(length(.x) >= 3, .x[3], NA_character_)),
      Order = map_chr(Taxon_Split, ~ ifelse(length(.x) >= 4, .x[4], NA_character_)),
      Family = map_chr(Taxon_Split, ~ ifelse(length(.x) >= 5, .x[5], NA_character_)),
      Genus = map_chr(Taxon_Split, ~ ifelse(length(.x) >= 6, .x[6], NA_character_)),
      Genus = ifelse(is.na(Genus) | Genus == "" | Genus == "NA",
                     paste0("Unclassified_", coalesce(Family, Order, "Higher")),
                     Genus),
      Phylum = ifelse(is.na(Phylum), "Unknown_Phylum", Phylum),
      Genus_Label = paste0(Phylum, "|", Genus, " (", ASV_ID, ")")
    ) %>%
    select(-Taxon_clean, -Taxon_Split)
  
  cat("Parsed", nrow(tax_df), "ASVs into", length(unique(tax_df$Genus)), "genera\n")
  return(tax_df)
}

aggregate_to_genus <- function(otu_table, tax_df) {
  common_asvs <- intersect(rownames(otu_table), tax_df$ASV_ID)
  if(length(common_asvs) < 10) stop("Low ASV matching rate:", length(common_asvs))
  
  otu_sub <- otu_table[common_asvs, , drop = FALSE]
  tax_sub <- tax_df %>% filter(ASV_ID %in% common_asvs)
  
  genus_mat <- rowsum(otu_sub, tax_sub$Genus_Label)
  colnames(genus_mat) <- clean_sample_ids(colnames(genus_mat))
  
  sample_sums <- colSums(genus_mat)
  if(any(sample_sums == 0)) {
    zero_samples <- names(sample_sums)[sample_sums == 0]
    cat("Warning: removed zero-count samples:", paste(zero_samples, collapse = ", "), "\n")
    genus_mat <- genus_mat[, sample_sums > 0, drop = FALSE]
  }
  
  prop <- genus_mat / colSums(genus_mat)
  keep <- rowSums(prop > 0.0005) >= 4
  genus_mat_filtered <- genus_mat[keep, , drop = FALSE]
  
  cat("Aggregated", nrow(genus_mat), "genera -> filtered to", nrow(genus_mat_filtered), "genera\n")
  return(genus_mat_filtered)
}

# --- ALDEx2 core analysis ---
run_aldex_comparison <- function(genus_table, meta, group_col, comp_name, prefix) {
  meta_row_names <- clean_sample_ids(rownames(meta))
  rownames(meta) <- meta_row_names
  
  common_samples <- intersect(colnames(genus_table), rownames(meta))
  if(length(common_samples) < 6) {
    cat("Skipped", comp_name, ": insufficient samples (n < 6)\n")
    return(NULL)
  }
  
  common_samples <- sort(common_samples)
  genus_table <- genus_table[, common_samples, drop = FALSE]
  meta <- meta[common_samples, , drop = FALSE]
  conds <- meta[[group_col]]
  
  if(any(is.na(conds))) {
    valid_idx <- !is.na(conds)
    genus_table <- genus_table[, valid_idx, drop = FALSE]
    meta <- meta[valid_idx, , drop = FALSE]
    conds <- conds[valid_idx]
  }
  
  conds <- droplevels(as.factor(conds))
  if(length(unique(conds)) != 2 || min(table(conds)) < 2) {
    cat("Skipped", comp_name, ": invalid grouping\n")
    return(NULL)
  }
  
  cat("Running", comp_name, ":", paste(levels(conds), collapse = " vs "), "\n")
  
  reads_matrix <- as.matrix(genus_table)
  conds_char <- as.character(conds)
  
  tryCatch({
    x <- aldex(reads_matrix, conds_char,
               mc.samples = 1000,
               test = "t",
               effect = TRUE,
               verbose = FALSE,
               gamma = 0.5)
    
    x$Genus_Label <- rownames(x)
    x$Comparison <- comp_name
    x$Group <- prefix
    
    parsed <- strsplit(x$Genus_Label, "[|()]")
    x$Phylum <- sapply(parsed, function(v) v[1])
    x$Genus_Name <- sapply(parsed, function(v) trimws(v[2]))
    x$ASV_ID <- sapply(parsed, function(v) ifelse(length(v) >= 4, v[3], NA))
    
    sig_count <- sum(x$we.eBH < 0.05, na.rm = TRUE)
    cat("Completed", comp_name, ":", sig_count, "significant genera\n")
    return(x)
  }, error = function(e) {
    cat("Error in", comp_name, ":", e$message, "\n")
    return(NULL)
  })
}

# --- Visualization functions ---
save_plot <- function(plot, filename, width, height) {
  if(is.null(plot)) return()
  tryCatch({
    if(capabilities("cairo")[["cairo"]]) {
      ggsave(filename, plot, width = width, height = height, dpi = 600, device = cairo_pdf)
    } else {
      ggsave(filename, plot, width = width, height = height, dpi = 600, device = pdf)
    }
    cat("Saved:", basename(filename), "\n")
  }, error = function(e) {
    cat("Failed to save plot:", e$message, "\n")
  })
}

create_cns_volcano <- function(aldex_result, title, subtitle) {
  if(is.null(aldex_result) || nrow(aldex_result) == 0) return(NULL)
  
  data <- as.data.frame(aldex_result) %>%
    mutate(
      significance = case_when(
        we.eBH < 0.01 & abs(effect) > 1.5 ~ "Highly Sig",
        we.eBH < 0.05 & abs(effect) > 1.0 ~ "Significant",
        TRUE ~ "Non-sig"
      ),
      logFDR = -log10(we.eBH + 1e-10)
    )
  
  top_labels <- data %>%
    filter(we.eBH < 0.05, abs(effect) > 0.8) %>%
    arrange(we.eBH) %>%
    head(10)
  
  colors <- c("Highly Sig" = "#D62728", "Significant" = "#FF7F0E", "Non-sig" = "#7F7F7F")
  
  p <- ggplot(data, aes(x = effect, y = logFDR, color = significance)) +
    geom_point(aes(size = abs(diff.btw)), alpha = 0.7) +
    scale_color_manual(values = colors, name = "Significance") +
    scale_size_continuous(range = c(2, 6), name = "|Diff. btw|") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40", linewidth = 0.5) +
    geom_vline(xintercept = c(-1, 0, 1), linetype = c("dashed", "solid", "dashed"),
               color = c("gray60", "black", "gray60"), linewidth = c(0.5, 0.8, 0.5)) +
    geom_label_repel(data = top_labels, aes(label = Genus_Name),
                     size = 3, box.padding = 0.4, point.padding = 0.3,
                     max.overlaps = 15, fontface = "italic", fill = "white", alpha = 0.9) +
    labs(title = title, subtitle = subtitle,
         x = "Effect Size (log2 FC)", y = expression(-log[10]~(FDR))) +
    theme_bw(base_family = "Arial", base_size = 12) +
    theme(
      panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "gray30"),
      axis.text = element_text(color = "black"),
      legend.position = "bottom",
      legend.box = "horizontal",
      aspect.ratio = 0.85
    )
  return(p)
}

create_effect_plot <- function(aldex_result, title, top_n = 15) {
  if(is.null(aldex_result) || nrow(aldex_result) == 0) return(NULL)
  
  top_data <- as.data.frame(aldex_result) %>%
    filter(we.eBH < 0.1) %>%
    arrange(we.eBH) %>%
    head(top_n) %>%
    mutate(
      Genus_Name = fct_reorder(Genus_Name, effect),
      Color = ifelse(we.eBH < 0.05, "Significant", "Trend"),
      CI_lower = effect - 1.96 * diff.win,
      CI_upper = effect + 1.96 * diff.win
    )
  
  if(nrow(top_data) == 0) return(NULL)
  
  colors_fill <- c("Significant" = "#2C3E50", "Trend" = "#95A5A6")
  
  p <- ggplot(top_data, aes(x = effect, y = Genus_Name, fill = Color)) +
    geom_col(width = 0.7, alpha = 0.9) +
    geom_errorbar(aes(xmin = CI_lower, xmax = CI_upper),
                  orientation = "y", width = 0.2, linewidth = 0.5, color = "gray30") +
    scale_fill_manual(values = colors_fill, name = "Significance") +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.8) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray50", linewidth = 0.4) +
    labs(title = title,
         subtitle = paste0("Top ", min(top_n, nrow(top_data)), " genera (FDR<0.1)"),
         x = "Effect Size (log2 FC, 95% CI)", y = NULL) +
    theme_bw(base_family = "Arial", base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.y = element_text(face = "italic", size = 9),
      plot.title = element_text(face = "bold", size = 12),
      legend.position = "bottom",
      legend.box = "horizontal"
    )
  return(p)
}

create_ma_plot <- function(aldex_result, title) {
  if(is.null(aldex_result) || nrow(aldex_result) == 0) return(NULL)
  
  data <- as.data.frame(aldex_result)
  if("rab.all" %in% colnames(data)) {
    data <- data %>% mutate(mean_abund = rab.all)
  } else if("rab.win" %in% colnames(data) && "rab.btw" %in% colnames(data)) {
    data <- data %>% mutate(mean_abund = (rab.win + rab.btw) / 2)
  } else {
    data <- data %>% mutate(mean_abund = abs(diff.btw))
  }
  
  data <- data %>% mutate(
    significant = ifelse(we.eBH < 0.05, "FDR<0.05", "Not Sig")
  )
  
  p <- ggplot(data, aes(x = mean_abund, y = effect)) +
    geom_point(aes(color = significant), size = 2, alpha = 0.7) +
    scale_color_manual(values = c("FDR<0.05" = "#E74C3C", "Not Sig" = "#3498DB"),
                       name = "Significance") +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed", color = "gray50", linewidth = 0.4) +
    geom_label_repel(data = filter(data, we.eBH < 0.05, abs(effect) > 1) %>% head(8),
                     aes(label = Genus_Name),
                     size = 3, box.padding = 0.3, force = 10, fontface = "italic") +
    labs(title = title,
         subtitle = "MA Plot: Effect vs Mean Abundance",
         x = "Mean Abundance (log2)", y = "Effect Size (log2 FC)") +
    theme_bw(base_family = "Arial", base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
      panel.grid.minor = element_blank()
    )
  return(p)
}

# --- F/B ratio analysis ---
analyze_fb_ratio <- function(meta, otu_16s, otu_its, comparisons) {
  cat("\n========== F/B Ratio Analysis ==========\n")
  
  bact_total <- colSums(otu_16s)
  fung_total <- colSums(otu_its)
  
  common_samples <- intersect(names(bact_total), names(fung_total))
  common_samples <- intersect(common_samples, rownames(meta))
  
  fb_data <- data.frame(
    Sample = common_samples,
    Bacteria_Total = bact_total[common_samples],
    Fungi_Total = fung_total[common_samples],
    F_B_Ratio = fung_total[common_samples] / bact_total[common_samples],
    Log_F_B_Ratio = log2((fung_total[common_samples] + 1) / (bact_total[common_samples] + 1)),
    stringsAsFactors = FALSE
  )
  
  fb_meta <- meta[common_samples, ]
  fb_data$Stage <- fb_meta$Stage
  fb_data$Position <- fb_meta$Position
  
  write.csv(fb_data, file.path(output_dir, "F_B_Ratio_Sample_Level.csv"), row.names = FALSE)
  
  fb_results <- list()
  
  for(comp in comparisons) {
    if(comp$type == "stage") {
      idx <- which(fb_data$Position == comp$pos & fb_data$Stage %in% comp$stages)
      sub_data <- fb_data[idx, ]
      if(nrow(sub_data) < 4) next
      group_vec <- factor(sub_data$Stage, levels = comp$stages)
    } else {
      idx <- which(fb_data$Stage == comp$stage)
      sub_data <- fb_data[idx, ]
      if(nrow(sub_data) < 4) next
      group_vec <- factor(sub_data$Position, levels = c("Internal", "External"))
    }
    
    if(length(unique(group_vec)) != 2 || min(table(group_vec)) < 2) next
    
    t_res <- t.test(Log_F_B_Ratio ~ group_vec, data = sub_data)
    w_res <- wilcox.test(Log_F_B_Ratio ~ group_vec, data = sub_data)
    
    g1 <- sub_data$Log_F_B_Ratio[group_vec == levels(group_vec)[1]]
    g2 <- sub_data$Log_F_B_Ratio[group_vec == levels(group_vec)[2]]
    
    paired_p <- NA
    if(length(g1) == length(g2) && length(g1) >= 3) {
      paired_w <- wilcox.test(g1, g2, paired = TRUE)
      paired_p <- paired_w$p.value
    }
    
    fb_results[[comp$name]] <- data.frame(
      Comparison = comp$name,
      Group1 = levels(group_vec)[1],
      Group2 = levels(group_vec)[2],
      n1 = sum(group_vec == levels(group_vec)[1]),
      n2 = sum(group_vec == levels(group_vec)[2]),
      F_B_Mean_G1 = mean(2^g1),
      F_B_Mean_G2 = mean(2^g2),
      Log2FC_F_B = mean(g2) - mean(g1),
      T_Test_p = t_res$p.value,
      Wilcox_p = w_res$p.value,
      Paired_Wilcox_p = paired_p,
      stringsAsFactors = FALSE
    )
  }
  
  if(length(fb_results) > 0) {
    fb_summary <- do.call(rbind, fb_results)
    fb_summary$FDR_T <- p.adjust(fb_summary$T_Test_p, method = "BH")
    fb_summary$FDR_W <- p.adjust(fb_summary$Wilcox_p, method = "BH")
    fb_summary$FDR_Paired <- ifelse(!is.na(fb_summary$Paired_Wilcox_p),
                                    p.adjust(fb_summary$Paired_Wilcox_p, method = "BH"), NA)
    write.csv(fb_summary, file.path(output_dir, "F_B_Ratio_Results.csv"), row.names = FALSE)
    cat("F/B analysis completed\n")
    print(fb_summary[, c("Comparison", "Log2FC_F_B", "FDR_T", "FDR_Paired")])
  }
  return(fb_results)
}

# --- Paired test function ---
run_paired_analysis <- function(genus_table, meta, group_col, comp_name, prefix, output_subdir) {
  cat("\nPaired test:", comp_name, "\n")
  
  common_samples <- intersect(colnames(genus_table), rownames(meta))
  if(length(common_samples) < 6) return(NULL)
  
  meta <- meta[common_samples, , drop = FALSE]
  genus_table <- genus_table[, common_samples, drop = FALSE]
  
  conds <- meta[[group_col]]
  if(any(is.na(conds))) {
    valid_idx <- !is.na(conds)
    genus_table <- genus_table[, valid_idx, drop = FALSE]
    meta <- meta[valid_idx, , drop = FALSE]
    conds <- conds[valid_idx]
  }
  
  conds <- droplevels(as.factor(conds))
  if(length(unique(conds)) != 2) return(NULL)
  
  g1_samples <- rownames(meta)[conds == levels(conds)[1]]
  g2_samples <- rownames(meta)[conds == levels(conds)[2]]
  
  if(length(g1_samples) != length(g2_samples)) {
    cat("Skipped paired test: unequal sample sizes\n")
    return(NULL)
  }
  
  paired_results <- data.frame(Genus = rownames(genus_table), stringsAsFactors = FALSE)
  paired_pvalues <- c()
  
  for(i in 1:nrow(genus_table)) {
    vals1 <- as.numeric(genus_table[i, g1_samples])
    vals2 <- as.numeric(genus_table[i, g2_samples])
    
    if(sum(vals1 != vals2) > 0) {
      wt <- wilcox.test(vals1, vals2, paired = TRUE)
      paired_pvalues[i] <- wt$p.value
    } else {
      paired_pvalues[i] <- 1
    }
  }
  
  paired_results$Paired_Wilcox_p <- paired_pvalues
  paired_results$FDR <- p.adjust(paired_pvalues, method = "BH")
  paired_results$Comparison <- comp_name
  
  write.csv(paired_results, file.path(output_dir, output_subdir,
                                      paste0(comp_name, "_Paired_Results.csv")), row.names = FALSE)
  cat("Paired test results saved to", output_subdir, "\n")
  return(paired_results)
}

# --- Main workflow ---
main_analysis <- function() {
  # Read metadata
  meta <- read.csv(input_metadata, row.names = 1, stringsAsFactors = TRUE, check.names = FALSE)
  meta$Stage <- factor(meta$Stage, levels = c("I", "III", "V"))
  meta$Position <- factor(meta$Position, levels = c("Internal", "External"))
  rownames(meta) <- clean_sample_ids(rownames(meta))
  
  cat("Metadata samples:", nrow(meta), "\n")
  cat("Experimental design:\n")
  print(table(meta$Stage, meta$Position))
  
  # Read 16S data
  cat("\nReading 16S data\n")
  biom_16s <- read_biom(input_16s_biom)
  otu_16s_raw <- as.matrix(biom_data(biom_16s))
  if(ncol(otu_16s_raw) > nrow(otu_16s_raw)) otu_16s_raw <- t(otu_16s_raw)
  colnames(otu_16s_raw) <- clean_sample_ids(colnames(otu_16s_raw))
  
  tax_16s <- read_taxonomy_genus(input_16s_tax)
  genus_16s <- aggregate_to_genus(otu_16s_raw, tax_16s)
  
  # Read ITS data
  cat("\nReading ITS data\n")
  biom_its <- read_biom(input_its_biom)
  otu_its_raw <- as.matrix(biom_data(biom_its))
  if(ncol(otu_its_raw) > nrow(otu_its_raw)) otu_its_raw <- t(otu_its_raw)
  colnames(otu_its_raw) <- clean_sample_ids(colnames(otu_its_raw))
  
  tax_its <- read_taxonomy_genus(input_its_tax)
  genus_its <- aggregate_to_genus(otu_its_raw, tax_its)
  
  # Define 16 comparisons: 12 stage comparisons + 4 position comparisons
  comparisons <- list(
    # Bacteria Internal (3)
    list(name = "Bac_Int_IvsIII", domain = "Bacteria", type = "stage", pos = "Internal", stages = c("I", "III"), main = FALSE),
    list(name = "Bac_Int_IvsV", domain = "Bacteria", type = "stage", pos = "Internal", stages = c("I", "V"), main = FALSE),
    list(name = "Bac_Int_IIIvsV", domain = "Bacteria", type = "stage", pos = "Internal", stages = c("III", "V"), main = TRUE),
    
    # Bacteria External (3)
    list(name = "Bac_Ext_IvsIII", domain = "Bacteria", type = "stage", pos = "External", stages = c("I", "III"), main = FALSE),
    list(name = "Bac_Ext_IvsV", domain = "Bacteria", type = "stage", pos = "External", stages = c("I", "V"), main = FALSE),
    list(name = "Bac_Ext_IIIvsV", domain = "Bacteria", type = "stage", pos = "External", stages = c("III", "V"), main = FALSE),
    
    # Fungi Internal (3)
    list(name = "Fung_Int_IvsIII", domain = "Fungi", type = "stage", pos = "Internal", stages = c("I", "III"), main = FALSE),
    list(name = "Fung_Int_IvsV", domain = "Fungi", type = "stage", pos = "Internal", stages = c("I", "V"), main = FALSE),
    list(name = "Fung_Int_IIIvsV", domain = "Fungi", type = "stage", pos = "Internal", stages = c("III", "V"), main = TRUE),
    
    # Fungi External (3)
    list(name = "Fung_Ext_IvsIII", domain = "Fungi", type = "stage", pos = "External", stages = c("I", "III"), main = FALSE),
    list(name = "Fung_Ext_IvsV", domain = "Fungi", type = "stage", pos = "External", stages = c("I", "V"), main = FALSE),
    list(name = "Fung_Ext_IIIvsV", domain = "Fungi", type = "stage", pos = "External", stages = c("III", "V"), main = FALSE),
    
    # Position comparisons (4)
    list(name = "Bac_PosIII", domain = "Bacteria", type = "position", stage = "III", main = TRUE),
    list(name = "Bac_PosV", domain = "Bacteria", type = "position", stage = "V", main = FALSE),
    list(name = "Fung_PosIII", domain = "Fungi", type = "position", stage = "III", main = FALSE),
    list(name = "Fung_PosV", domain = "Fungi", type = "position", stage = "V", main = TRUE)
  )
  
  cat("\nStarting ALDEx2 analysis (", length(comparisons), " comparisons)\n")
  
  # F/B ratio analysis
  fb_res <- analyze_fb_ratio(meta, otu_16s_raw, otu_its_raw, comparisons)
  
  all_results <- list()
  main_plots_vol <- list()
  main_plots_eff <- list()
  
  # Execute all comparisons
  for(i in 1:length(comparisons)) {
    comp <- comparisons[[i]]
    cat("\n[", i, "/", length(comparisons), "] ", comp$name, "\n")
    
    genus_table <- if(comp$domain == "Bacteria") genus_16s else genus_its
    
    if(comp$type == "stage") {
      idx <- which(meta$Position == comp$pos & meta$Stage %in% comp$stages)
      sub_meta <- meta[idx, , drop = FALSE]
      sub_meta$Stage <- droplevels(factor(sub_meta$Stage, levels = comp$stages))
      group_col <- "Stage"
      comp_label <- paste0(comp$domain, " ", comp$pos, ": ", paste(comp$stages, collapse = "->"))
    } else {
      idx <- which(meta$Stage == comp$stage)
      sub_meta <- meta[idx, , drop = FALSE]
      group_col <- "Position"
      comp_label <- paste0(comp$domain, " Stage", comp$stage, ": Int vs Ext")
    }
    
    res <- run_aldex_comparison(genus_table, sub_meta, group_col, comp_label, comp$domain)
    
    if(!is.null(res)) {
      all_results[[comp$name]] <- res
      
      write.csv(res, file.path(output_dir, paste0(comp$name, "_Results.csv")), row.names = FALSE)
      
      p_vol <- create_cns_volcano(res, comp$name, paste0("gamma=0.5 | ", comp$domain))
      p_eff <- create_effect_plot(res, paste0(comp$name, " Effect"))
      p_ma <- create_ma_plot(res, paste0(comp$name, " MA"))
      
      save_plot(p_vol, file.path(output_dir, paste0(comp$name, "_Volcano.pdf")), 9, 8)
      save_plot(p_eff, file.path(output_dir, paste0(comp$name, "_Effect.pdf")), 9, 7)
      save_plot(p_ma, file.path(output_dir, paste0(comp$name, "_MA.pdf")), 8, 6)
      
      if(!is.null(comp$main) && comp$main) {
        main_plots_vol[[comp$name]] <- p_vol
        main_plots_eff[[comp$name]] <- p_eff
        
        run_paired_analysis(genus_table, sub_meta, group_col, comp$name, comp$domain, "Paired_Tests")
      }
    }
  }
  
  # Generate main figure (8-in-1)
  cat("\nGenerating main figure (8-in-1)\n")
  
  plot_order <- c("Bac_PosIII", "Bac_Int_IIIvsV", "Fung_PosV", "Fung_Int_IIIvsV")
  
  vol_list <- list()
  eff_list <- list()
  
  for(i in 1:4) {
    key <- plot_order[i]
    if(!is.null(main_plots_vol[[key]])) {
      vol_list[[i]] <- main_plots_vol[[key]] +
        ggtitle(paste0(letters[(i - 1) * 2 + 1], ") ", key))
      eff_list[[i]] <- main_plots_eff[[key]] +
        ggtitle(paste0(letters[(i - 1) * 2 + 2], ") ", key, " Top Genera"))
    } else {
      cat("Warning: missing main plot", key, "\n")
    }
  }
  
  if(length(vol_list) == 4 && length(eff_list) == 4) {
    left_col <- plot_grid(plotlist = vol_list, ncol = 1, rel_heights = c(1, 1, 1, 1))
    right_col <- plot_grid(plotlist = eff_list, ncol = 1, rel_heights = c(1, 1, 1, 1))
    
    main_fig <- plot_grid(left_col, right_col, ncol = 2, rel_widths = c(1, 1))
    main_fig <- ggdraw(main_fig) +
      draw_label("ALDEx2 Scale Model Analysis of Decayed Wood Microbiome (gamma=0.5)",
                 x = 0.5, y = 0.98, size = 14, fontface = "bold", fontfamily = "Arial")
    
    save_plot(main_fig, file.path(output_dir, "Main_Figure", "Main_Figure_8in1.pdf"), 16, 20)
    cat("Main figure saved\n")
  } else {
    cat("Main figure generation failed: insufficient plots\n")
  }
  
  # Generate SCI summary table
  cat("\nGenerating SCI summary table\n")
  if(length(all_results) > 0) {
    summary_df <- data.frame(
      Comparison = names(all_results),
      Domain = sapply(all_results, function(x) unique(x$Group)),
      Total_Genera = sapply(all_results, nrow),
      Sig_FDR05 = sapply(all_results, function(x) sum(x$we.eBH < 0.05)),
      Sig_FDR01 = sapply(all_results, function(x) sum(x$we.eBH < 0.01)),
      Max_Effect = sapply(all_results, function(x) {
        sig <- x[x$we.eBH < 0.05, ]
        ifelse(nrow(sig) > 0, max(abs(sig$effect)), 0)
      }),
      stringsAsFactors = FALSE
    )
    
    if(length(fb_res) > 0) {
      fb_df <- do.call(rbind, fb_res)
      rownames(fb_df) <- fb_df$Comparison
      summary_df$F_B_Log2FC <- fb_df[summary_df$Comparison, "Log2FC_F_B"]
      summary_df$F_B_FDR <- fb_df[summary_df$Comparison, "FDR_T"]
    }
    
    write.csv(summary_df, file.path(output_dir, "ALDEx2_Summary_SCI.csv"), row.names = FALSE)
    cat("Summary table saved\n")
    print(summary_df)
  }
  
  cat("\nAnalysis completed\n")
  cat("Output directory:", output_dir, "\n")
  cat("- Individual results: root directory CSV + PDF\n")
  cat("- Paired tests: Paired_Tests/ subdirectory\n")
  cat("- F/B ratio: F_B_Ratio_Results.csv\n")
  cat("- Main figure: Main_Figure/Main_Figure_8in1.pdf\n")
}

main_analysis()