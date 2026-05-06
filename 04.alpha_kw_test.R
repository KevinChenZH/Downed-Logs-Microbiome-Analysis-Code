# alpha_kw_test.R
# Required: readxl, dplyr, tidyr, openxlsx
# Input: Excel with sample IDs in first column (e.g., n1.1, w5.4)
# Output: 3-sheet Excel (Mean±SE pivot, raw summary, K-W test results)

library(readxl)
library(dplyr)
library(tidyr)
library(openxlsx)

METRICS <- c("chao1", "observed_features", "shannon", "pielou_e")

parse_name <- function(x) {
  x <- trimws(as.character(x))
  sm <- c("1" = "I", "3" = "III", "5" = "V")
  p <- substr(x, 1, 1)
  s <- sm[substr(x, 2, 2)]
  if (p == "n") return(c(s, "Internal"))
  if (p == "w") return(c(s, "External"))
  c(NA_character_, NA_character_)
}

get_sig <- function(p) {
  if (is.na(p)) return("na")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  "ns"
}

process_diversity <- function(file_path, group_label) {
  df <- read_excel(file_path, col_names = TRUE)
  pr <- t(vapply(df[[1]], parse_name, character(2)))
  df$Decay_Stage <- factor(pr[, 1], levels = c("I", "III", "V"))
  df$Position <- factor(pr[, 2], levels = c("Internal", "External"))
  df <- df %>% filter(!is.na(Decay_Stage), !is.na(Position))
  df <- df %>% mutate(across(all_of(intersect(METRICS, names(df))), as.numeric))
  
  summary_df <- df %>%
    group_by(Decay_Stage, Position) %>%
    summarise(across(all_of(intersect(METRICS, names(.))), list(
      Mean = ~ mean(.x, na.rm = TRUE),
      SE = ~ sd(.x, na.rm = TRUE) / sqrt(sum(!is.na(.x))),
      n = ~ sum(!is.na(.x))
    )), .groups = "drop") %>%
    pivot_longer(cols = -c(Decay_Stage, Position), names_to = c("Index", ".value"), names_sep = "_") %>%
    mutate(Microbial_group = group_label, Mean_SE = sprintf("%.2f±%.2f", Mean, SE))
  
  stats_rows <- list()
  for (col in intersect(METRICS, names(df))) {
    st_list <- list()
    for (st in c("I", "III", "V")) {
      d <- df %>% filter(Decay_Stage == st) %>% pull(!!sym(col)) %>% na.omit()
      if (length(d) > 0) st_list[[st]] <- d
    }
    if (length(st_list) >= 2) {
      kw_s <- kruskal.test(st_list)
      h_s <- kw_s$statistic
      p_s <- kw_s$p.value
    } else {
      h_s <- NA; p_s <- NA
    }
    
    int <- df %>% filter(Position == "Internal") %>% pull(!!sym(col)) %>% na.omit()
    ext <- df %>% filter(Position == "External") %>% pull(!!sym(col)) %>% na.omit()
    if (length(int) > 0 && length(ext) > 0) {
      kw_p <- kruskal.test(list(int, ext))
      h_p <- kw_p$statistic
      p_p <- kw_p$p.value
    } else {
      h_p <- NA; p_p <- NA
    }
    
    stats_rows[[col]] <- data.frame(
      Microbial_group = group_label,
      Index = col,
      Stage_effect_p = ifelse(is.na(p_s), "na", sprintf("%.4f", p_s)),
      Stage_effect_sig = get_sig(p_s),
      Position_effect_p = ifelse(is.na(p_p), "na", sprintf("%.4f", p_p)),
      Position_effect_sig = get_sig(p_p),
      H_value_Stage = ifelse(is.na(h_s), "na", round(h_s, 2)),
      H_value_Pos = ifelse(is.na(h_p), "na", round(h_p, 2)),
      stringsAsFactors = FALSE
    )
  }
  stats_df <- bind_rows(stats_rows)
  
  list(summary = summary_df, stats = stats_df)
}

bacteria <- process_diversity("path/to/16S_alpha_diversity.xlsx", "Bacteria")
fungi <- process_diversity("path/to/ITS_alpha_diversity.xlsx", "Fungi")

summary_all <- bind_rows(bacteria$summary, fungi$summary)
stats_all <- bind_rows(bacteria$stats, fungi$stats)

pivot <- summary_all %>%
  select(Microbial_group, Decay_Stage, Position, Index, Mean_SE) %>%
  pivot_wider(names_from = Index, values_from = Mean_SE) %>%
  select(Microbial_group, Decay_Stage, Position, all_of(intersect(METRICS, names(.))))

wb <- createWorkbook()
addWorksheet(wb, "Table_Mean_SE")
writeData(wb, "Table_Mean_SE", pivot)
addWorksheet(wb, "Raw_Mean_SE_n")
writeData(wb, "Raw_Mean_SE_n", summary_all %>% select(Microbial_group, Decay_Stage, Position, Index, Mean, SE, n))
addWorksheet(wb, "Statistical_Tests")
writeData(wb, "Statistical_Tests", stats_all)
saveWorkbook(wb, "path/to/output/Alpha_Diversity_KW.xlsx", overwrite = TRUE)