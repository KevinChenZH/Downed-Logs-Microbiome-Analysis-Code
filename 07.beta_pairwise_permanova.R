library(vegan)
library(dplyr)
library(readr)

load_dm <- function(path) {
  df <- read_tsv(path, col_names = TRUE, show_col_types = FALSE)
  mat <- as.matrix(df[, -1])
  rownames(mat) <- df[[1]]
  class(mat) <- "numeric"
  as.dist(mat)
}

pairwise_perm <- function(dm, meta, factor, label) {
  ids <- labels(dm)
  meta <- meta[match(ids, meta$SampleID), ]
  lvls <- sort(unique(meta[[factor]]))
  out <- list()
  for (i in seq_len(length(lvls) - 1)) {
    for (j in (i + 1):length(lvls)) {
      l1 <- lvls[i]; l2 <- lvls[j]
      mask <- meta[[factor]] %in% c(l1, l2)
      sub_ids <- ids[mask]
      if (length(sub_ids) < 3) next
      sub_dm <- as.dist(as.matrix(dm)[sub_ids, sub_ids])
      sub_meta <- meta[match(sub_ids, meta$SampleID), ]
      df_tmp <- data.frame(group = sub_meta[[factor]], stringsAsFactors = FALSE)
      res <- adonis2(sub_dm ~ group, data = df_tmp, permutations = 999, by = "margin")
      r2 <- round(res$R2[1], 3)
      p <- res$`Pr(>F)`[1]
      p_str <- ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
      out[[length(out) + 1]] <- data.frame(
        Group = label, Comparison = paste(l1, "vs", l2), Factor = factor,
        R2 = r2, p_value = p_str, stringsAsFactors = FALSE
      )
    }
  }
  bind_rows(out)
}

simple_effect <- function(dm, meta, factor, cond_factor, cond_val, label) {
  ids <- labels(dm)
  meta <- meta[match(ids, meta$SampleID), ]
  mask <- meta[[cond_factor]] == cond_val
  if (sum(mask) == 0) return(NULL)
  sub_ids <- ids[mask]
  sub_dm <- as.dist(as.matrix(dm)[sub_ids, sub_ids])
  sub_meta <- meta[match(sub_ids, meta$SampleID), ]
  lvls <- sort(unique(sub_meta[[factor]]))
  if (length(lvls) < 2) return(NULL)
  out <- list()
  for (i in seq_len(length(lvls) - 1)) {
    for (j in (i + 1):length(lvls)) {
      l1 <- lvls[i]; l2 <- lvls[j]
      pair_mask <- sub_meta[[factor]] %in% c(l1, l2)
      pair_ids <- sub_ids[pair_mask]
      if (length(pair_ids) < 3) next
      pair_dm <- as.dist(as.matrix(sub_dm)[pair_ids, pair_ids])
      pair_meta <- sub_meta[match(pair_ids, sub_meta$SampleID), ]
      df_tmp <- data.frame(group = pair_meta[[factor]], stringsAsFactors = FALSE)
      res <- adonis2(pair_dm ~ group, data = df_tmp, permutations = 999, by = "margin")
      r2 <- round(res$R2[1], 3)
      p <- res$`Pr(>F)`[1]
      p_str <- ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
      out[[length(out) + 1]] <- data.frame(
        Group = label,
        Comparison = paste0(l1, " vs ", l2, " (in ", cond_val, ")"),
        Factor = paste0(factor, "|", cond_factor),
        R2 = r2, p_value = p_str, stringsAsFactors = FALSE
      )
    }
  }
  bind_rows(out)
}

main <- function() {
  meta <- read_csv("path/to/metadata.csv", show_col_types = FALSE)
  dm_b <- load_dm("path/to/16S_distance_matrix.tsv")
  dm_f <- load_dm("path/to/ITS_distance_matrix.tsv")
  
  res <- list(
    pairwise_perm(dm_b, meta, "Stage", "Bacteria"),
    pairwise_perm(dm_b, meta, "Position", "Bacteria"),
    simple_effect(dm_b, meta, "Stage", "Position", "Internal", "Bacteria"),
    simple_effect(dm_b, meta, "Stage", "Position", "External", "Bacteria"),
    simple_effect(dm_b, meta, "Position", "Stage", "I", "Bacteria"),
    simple_effect(dm_b, meta, "Position", "Stage", "III", "Bacteria"),
    simple_effect(dm_b, meta, "Position", "Stage", "V", "Bacteria"),
    pairwise_perm(dm_f, meta, "Stage", "Fungi"),
    pairwise_perm(dm_f, meta, "Position", "Fungi"),
    simple_effect(dm_f, meta, "Stage", "Position", "Internal", "Fungi"),
    simple_effect(dm_f, meta, "Stage", "Position", "External", "Fungi"),
    simple_effect(dm_f, meta, "Position", "Stage", "I", "Fungi"),
    simple_effect(dm_f, meta, "Position", "Stage", "III", "Fungi"),
    simple_effect(dm_f, meta, "Position", "Stage", "V", "Fungi")
  )
  
  bind_rows(res) %>% write_csv("path/to/output/Pairwise_PERMANOVA.csv")
}

main()