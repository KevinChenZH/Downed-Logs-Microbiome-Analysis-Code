# Stacked Bar Plot for Microbial Community Composition
# CNS-level output: 600 dpi, editable Arial fonts, italic bacterial genus names
# Required input files:
#   1. Metadata CSV with columns: SampleID, Stage (I/III/V), Position (Internal/External)
#   2. 16S feature table (TSV/CSV/XLS) with QIIME-style taxonomy strings
#   3. ITS feature table (TSV/CSV/XLS) with QIIME-style taxonomy strings

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(cowplot)
})

has_readxl <- requireNamespace("readxl", quietly = TRUE)

# --- Configuration ---
BASE_PATH <- "BASE_DIRECTORY_PATH"
BACTERIA_FILE <- file.path(BASE_PATH, "16SfeatureTable.group.total.relative.xls")
FUNGI_FILE    <- file.path(BASE_PATH, "ITSfeatureTable.group.total.relative.xls")
METADATA_FILE <- file.path(BASE_PATH, "Metadata.csv")
OUTPUT_PDF    <- file.path(BASE_PATH, "Fig3_Community_Composition_CNS_Editable.pdf")
OUTPUT_TIFF   <- file.path(BASE_PATH, "Fig3_Community_Composition_CNS_Editable.tiff")

BACTERIA_PALETTE <- c(
  "Burkholderia-Caballeronia-Paraburkholderia" = "#1f497d",
  "Pseudomonas"      = "#4a7c9b",
  "Sphingomonas"     = "#76a5c4",
  "Luteibacter"      = "#17becf",
  "Nocardioides"     = "#c55a11",
  "Mycobacterium"    = "#e08a4c",
  "Streptomyces"     = "#ffbb78",
  "Bacillus"         = "#2ca02c",
  "unclassified_Bacillales" = "#98df8a",
  "Granulicella"     = "#d62728",
  "Acidibacter"      = "#9467bd",
  "unclassified_Bacteroidales" = "#8c564b",
  "Gryllotalpicola"  = "#bcbd22",
  "Jatrophihabitans" = "#e377c2",
  "unclassified_Micropepsaceae" = "#f7b6d3",
  "Others"           = "#c7c7c7"
)

FUNGI_PALETTE <- c(
  "Ascomycota"            = "#ff7f0e",
  "Basidiomycota"         = "#1f77b4",
  "Mortierellomycota"     = "#8c564b",
  "Chytridiomycota"       = "#2ca02c",
  "Mucoromycota"          = "#d62728",
  "Rozellomycota"         = "#9467bd",
  "Olpidiomycota"         = "#7f7f7f",
  "Glomeromycota"         = "#bcbd22",
  "Fungi_phy_Incertae_sedis" = "#e377c2",
  "Others"                = "#c7c7c7"
)

# --- Utility functions ---

parse_taxonomy <- function(tax_str, level = "genus") {
  if (is.na(tax_str) || tax_str == "") return("Unclassified")
  parts <- strsplit(as.character(tax_str), ";")[[1]]
  prefix <- paste0(substr(level, 1, 1), "__")
  for (part in parts) {
    part <- trimws(part)
    if (grepl(paste0("^", prefix), part)) {
      name <- gsub(prefix, "", part)
      name <- trimws(name)
      if (name %in% c("", "NA", "uncultured", "Unknown")) {
        higher <- "Unknown"
        if (level == "genus") {
          for (p in parts) {
            if (grepl("^f__", p)) { higher <- gsub("^f__", "", p); break }
          }
        }
        return(ifelse(higher != "Unknown", paste0("unclassified_", higher), "Unclassified"))
      }
      return(name)
    }
  }
  return("Unclassified")
}

get_phylum <- function(tax_str) {
  if (is.na(tax_str)) return("Unknown")
  for (part in strsplit(as.character(tax_str), ";")[[1]]) {
    if (grepl("^p__", part)) {
      name <- trimws(gsub("^p__", "", part))
      return(ifelse(name == "", "Unknown", name))
    }
  }
  return("Unknown")
}

read_feature_table <- function(file_path) {
  df <- tryCatch({
    read.delim(file_path, check.names = FALSE, stringsAsFactors = FALSE)
  }, error = function(e) {
    tryCatch({
      read.csv(file_path, check.names = FALSE, stringsAsFactors = FALSE)
    }, error = function(e2) {
      if (has_readxl) {
        as.data.frame(readxl::read_excel(file_path))
      } else {
        stop("Unable to read file (install 'readxl' for Excel support): ", e2$message)
      }
    })
  })
  return(df)
}

load_and_aggregate <- function(file_path, level = "genus") {
  df <- read_feature_table(file_path)
  
  id_col <- grep("OTU|Feature|ASV|#OTU", colnames(df), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(id_col)) id_col <- colnames(df)[1]
  
  tax_col <- grep("Taxon", colnames(df), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(tax_col)) tax_col <- colnames(df)[ncol(df)]
  
  sample_cols <- setdiff(colnames(df), c(id_col, tax_col))
  
  if (level == "genus") {
    df$Taxon  <- sapply(df[[tax_col]], parse_taxonomy, level = "genus")
    df$Phylum <- sapply(df[[tax_col]], get_phylum)
  } else {
    df$Taxon  <- sapply(df[[tax_col]], parse_taxonomy, level = "phylum")
    df$Phylum <- df$Taxon
  }
  
  df <- df %>% filter(Taxon != "Unclassified")
  
  agg <- df %>%
    group_by(Taxon) %>%
    summarise(across(all_of(sample_cols), sum), .groups = "drop") %>%
    as.data.frame()
  rownames(agg) <- agg$Taxon
  agg$Taxon <- NULL
  
  tax_map <- df %>% select(Taxon, Phylum) %>% distinct() %>% as.data.frame()
  rownames(tax_map) <- tax_map$Taxon
  
  # Remove zero-count samples
  cs <- colSums(agg)
  zero_cols <- names(cs)[cs == 0]
  if (length(zero_cols) > 0) {
    agg <- agg[, setdiff(colnames(agg), zero_cols), drop = FALSE]
    cs <- cs[setdiff(names(cs), zero_cols)]
  }
  
  agg <- sweep(agg, 2, cs, "/") * 100
  return(list(agg_df = agg, tax_map = tax_map))
}

calculate_group_means <- function(agg_df, metadata, stage_order = c("I", "III", "V")) {
  positions <- c("Internal", "External")
  res <- list()
  for (stage in stage_order) {
    for (pos in positions) {
      samples <- metadata %>% 
        filter(Stage == stage, Position == pos) %>% 
        pull(SampleID) %>% 
        as.character()
      
      avail <- intersect(samples, colnames(agg_df))
      if (length(avail) == 0) {
        avail <- c()
        for (s in samples) {
          base <- strsplit(s, "\\.")[[1]][1]
          avail <- c(avail, colnames(agg_df)[grepl(paste0("^", base), colnames(agg_df))])
        }
        avail <- unique(avail)
      }
      
      if (length(avail) > 0) {
        res[[paste0(stage, "_", pos)]] <- rowMeans(agg_df[, avail, drop = FALSE])
      } else {
        res[[paste0(stage, "_", pos)]] <- setNames(rep(0, nrow(agg_df)), rownames(agg_df))
      }
    }
  }
  mean_df <- as.data.frame(res)
  mean_df <- sweep(mean_df, 2, colSums(mean_df), "/") * 100
  return(mean_df)
}

select_top_taxa <- function(mean_df, top_n = 15) {
  om <- rowMeans(mean_df)
  top_taxa <- names(sort(om, decreasing = TRUE))[seq_len(top_n)]
  others <- 100 - colSums(mean_df[top_taxa, , drop = FALSE])
  final <- rbind(mean_df[top_taxa, , drop = FALSE], Others = others)
  return(list(final_df = final, top_taxa = top_taxa))
}

prepare_plot_data <- function(final_df, palette) {
  pd <- final_df %>%
    mutate(Taxon = rownames(final_df)) %>%
    pivot_longer(cols = -Taxon, names_to = "Group", values_to = "Abundance") %>%
    mutate(
      Stage    = factor(gsub("_.*", "", Group), levels = c("I", "III", "V")),
      Position = factor(gsub(".*_", "", Group), levels = c("Internal", "External")),
      Taxon    = factor(Taxon, levels = rev(rownames(final_df)))
    )
  
  all_taxa <- unique(as.character(pd$Taxon))
  cmap <- palette
  missing <- setdiff(all_taxa, names(cmap))
  backup <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd")
  for (i in seq_along(missing)) {
    cmap[missing[i]] <- backup[(i - 1) %% length(backup) + 1]
  }
  cmap["Others"] <- "#c7c7c7"
  
  return(list(plot_data = pd, color_map = cmap))
}

make_legend_labels <- function(taxa, group_type = "bacteria") {
  as.expression(lapply(taxa, function(t) {
    if (group_type == "bacteria" && t != "Others" && !grepl("^unclassified", t)) {
      bquote(italic(.(t)))
    } else {
      t
    }
  }))
}

# --- Main workflow ---

main <- function() {
  cat("CNS-level Stacked Bar Plot Generation\n")
  
  metadata <- read.csv(METADATA_FILE, stringsAsFactors = FALSE)
  metadata$Stage    <- factor(metadata$Stage,    levels = c("I", "III", "V"))
  metadata$Position <- factor(metadata$Position, levels = c("Internal", "External"))
  cat("Metadata samples:", nrow(metadata), "\n")
  
  cat("Processing bacteria (genus level)...\n")
  bac_res  <- load_and_aggregate(BACTERIA_FILE, level = "genus")
  bac_mean <- calculate_group_means(bac_res$agg_df, metadata)
  bac_sel  <- select_top_taxa(bac_mean, top_n = 15)
  
  cat("Processing fungi (phylum level)...\n")
  fun_res  <- load_and_aggregate(FUNGI_FILE, level = "phylum")
  fun_mean <- calculate_group_means(fun_res$agg_df, metadata)
  fun_sel  <- select_top_taxa(fun_mean, top_n = 10)
  
  cat("Generating plots...\n")
  
  bac_plot <- prepare_plot_data(bac_sel$final_df, BACTERIA_PALETTE)
  fun_plot <- prepare_plot_data(fun_sel$final_df, FUNGI_PALETTE)
  
  # Bacteria plot
  p_bac <- ggplot(bac_plot$plot_data, aes(x = Stage, y = Abundance, fill = Taxon)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7, colour = "white", linewidth = 0.3) +
    facet_wrap(~Position, ncol = 2, strip.position = "bottom") +
    scale_fill_manual(
      values = bac_plot$color_map,
      name = "Genus",
      labels = make_legend_labels(levels(bac_plot$plot_data$Taxon), "bacteria")
    ) +
    scale_x_discrete(labels = c("I" = "Stage I", "III" = "Stage III", "V" = "Stage V")) +
    labs(title = "Bacterial Community Composition", y = "Relative Abundance (%)") +
    coord_cartesian(ylim = c(0, 100)) +
    theme_bw(base_family = "Arial", base_size = 10) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "#333333", linewidth = 0.8),
      axis.ticks = element_line(colour = "#333333", linewidth = 0.8),
      axis.text = element_text(colour = "#333333", size = 10),
      axis.title.y = element_text(colour = "#333333", size = 11, face = "bold"),
      plot.title = element_text(colour = "#333333", size = 12, face = "bold", hjust = 0.5),
      strip.background = element_blank(),
      strip.text = element_text(colour = "#333333", size = 10, face = "bold"),
      strip.placement = "outside",
      legend.position = "none",
      plot.margin = margin(10, 5, 10, 5)
    )
  
  # Fungi plot
  p_fun <- ggplot(fun_plot$plot_data, aes(x = Stage, y = Abundance, fill = Taxon)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7, colour = "white", linewidth = 0.3) +
    facet_wrap(~Position, ncol = 2, strip.position = "bottom") +
    scale_fill_manual(
      values = fun_plot$color_map,
      name = "Phylum",
      labels = make_legend_labels(levels(fun_plot$plot_data$Taxon), "fungi")
    ) +
    scale_x_discrete(labels = c("I" = "Stage I", "III" = "Stage III", "V" = "Stage V")) +
    labs(title = "Fungal Community Composition (Phylum Level)", y = "Relative Abundance (%)") +
    coord_cartesian(ylim = c(0, 100)) +
    theme_bw(base_family = "Arial", base_size = 10) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "#333333", linewidth = 0.8),
      axis.ticks = element_line(colour = "#333333", linewidth = 0.8),
      axis.text = element_text(colour = "#333333", size = 10),
      axis.title.y = element_text(colour = "#333333", size = 11, face = "bold"),
      plot.title = element_text(colour = "#333333", size = 12, face = "bold", hjust = 0.5),
      strip.background = element_blank(),
      strip.text = element_text(colour = "#333333", size = 10, face = "bold"),
      strip.placement = "outside",
      legend.position = "none",
      plot.margin = margin(10, 5, 10, 5)
    )
  
  # Extract legends
  leg_bac <- get_legend(
    ggplot(bac_plot$plot_data, aes(x = Stage, y = Abundance, fill = Taxon)) +
      geom_bar(stat = "identity", position = "stack", width = 0.7) +
      scale_fill_manual(
        values = bac_plot$color_map, name = "Genus",
        labels = make_legend_labels(levels(bac_plot$plot_data$Taxon), "bacteria")
      ) +
      theme_void(base_family = "Arial", base_size = 9) +
      theme(
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.4, "cm"),
        legend.spacing.y = unit(0.15, "cm")
      )
  )
  
  leg_fun <- get_legend(
    ggplot(fun_plot$plot_data, aes(x = Stage, y = Abundance, fill = Taxon)) +
      geom_bar(stat = "identity", position = "stack", width = 0.7) +
      scale_fill_manual(
        values = fun_plot$color_map, name = "Phylum",
        labels = make_legend_labels(levels(fun_plot$plot_data$Taxon), "fungi")
      ) +
      theme_void(base_family = "Arial", base_size = 9) +
      theme(
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.4, "cm"),
        legend.spacing.y = unit(0.15, "cm")
      )
  )
  
  # Assemble panels
  upper <- plot_grid(p_bac, leg_bac, ncol = 2, rel_widths = c(4, 1), align = "v", axis = "l")
  lower <- plot_grid(p_fun, leg_fun, ncol = 2, rel_widths = c(4, 1), align = "v", axis = "l")
  combined <- plot_grid(upper, lower, ncol = 1, rel_heights = c(1, 1))
  combined <- ggdraw(combined) +
    draw_label("Microbial Community Succession During Wood Decay",
               x = 0.5, y = 0.98, size = 14, fontface = "bold", fontfamily = "Arial")
  
  # Save
  if (capabilities("cairo")) {
    ggsave(OUTPUT_PDF, combined, width = 14, height = 10, dpi = 600, device = cairo_pdf)
  } else {
    ggsave(OUTPUT_PDF, combined, width = 14, height = 10, dpi = 600, device = "pdf")
  }
  ggsave(OUTPUT_TIFF, combined, width = 14, height = 10, dpi = 600, 
         device = "tiff", compression = "lzw")
  
  cat("Saved PDF:", OUTPUT_PDF, "\n")
  cat("Saved TIFF:", OUTPUT_TIFF, "\n")
  cat("Completed. Arial fonts editable; bacterial genus names italicized.\n")
}

main()