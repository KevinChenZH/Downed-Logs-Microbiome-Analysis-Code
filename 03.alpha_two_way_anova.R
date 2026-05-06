# alpha_two_way_anova.R
# Required packages: readxl, dplyr, tidyr, car, openxlsx
# Install via: install.packages(c("readxl", "dplyr", "tidyr", "car", "openxlsx"))
# Input: Excel with first column = sample ID (e.g., n1.1, w3.2), columns: shannon, chao1, observed_features, pielou_e
# Output: Two Excel files (long format + wide format for paper table)

library(readxl)
library(dplyr)
library(tidyr)
library(car)
library(openxlsx)

METRICS <- c("shannon", "chao1", "observed_features", "pielou_e")

parse_name <- function(x) {
  x <- trimws(as.character(x))
  sm <- c("1" = "I", "3" = "III", "5" = "V")
  p <- substr(x, 1, 1)
  s <- sm[substr(x, 2, 2)]
  if (p == "n") return(c(s, "Internal"))
  if (p == "w") return(c(s, "External"))
  c(NA_character_, NA_character_)
}

run_anova <- function(file_path, group_label) {
  df <- read_excel(file_path, col_names = TRUE)
  pr <- t(vapply(df[[1]], parse_name, character(2)))
  df$Stage <- factor(pr[, 1], levels = c("I", "III", "V"))
  df$Position <- factor(pr[, 2], levels = c("Internal", "External"))
  df <- df %>% filter(!is.na(Stage), !is.na(Position))
  df <- df %>% mutate(across(all_of(intersect(METRICS, names(df))), as.numeric))
  
  out <- list()
  for (col in intersect(METRICS, names(df))) {
    sub <- df %>% filter(!is.na(.data[[col]]))
    if (nrow(sub) < 3) next
    
    shp <- sub %>%
      group_by(Stage, Position) %>%
      summarise(p = ifelse(n() >= 3, shapiro.test(.data[[col]])$p.value, NA), .groups = "drop")
    normal_ok <- all(shp$p > 0.05, na.rm = TRUE)
    
    lev <- leveneTest(as.formula(paste(col, "~ Stage * Position")), data = sub)
    equal_var <- lev$`Pr(>F)`[1] > 0.05
    
    fit <- lm(as.formula(paste(col, "~ Stage * Position")), data = sub)
    aov_tbl <- if (equal_var) Anova(fit, type = "II") else Anova(fit, type = "III", white.adjust = "hc3")
    
    ss_res <- aov_tbl["Residuals", "Sum Sq"]
    for (i in seq_len(nrow(aov_tbl))) {
      term <- rownames(aov_tbl)[i]
      if (term == "Residuals") next
      fac <- gsub(":", "×", term)
      f_val <- aov_tbl[i, "F"]
      p_val <- aov_tbl[i, "Pr(>F)"]
      if (is.na(f_val)) next
      sig <- ifelse(p_val < 0.001, "***", ifelse(p_val < 0.01, "**", ifelse(p_val < 0.05, "*", "ns")))
      ss_eff <- aov_tbl[i, "Sum Sq"]
      eta_sq <- ss_eff / (ss_eff + ss_res)
      
      out[[length(out) + 1]] <- data.frame(
        Microbial_group = group_label, Index = col, Factor = fac,
        F_value = round(f_val, 2), p_value = round(p_val, 4),
        Significance = sig, Partial_Eta2 = round(eta_sq, 3),
        Assumptions = ifelse(normal_ok && equal_var, "Met", "Partially Met"),
        Method = ifelse(equal_var, "Two-way ANOVA (Type II)", "Two-way ANOVA (Type III, robust)"),
        stringsAsFactors = FALSE
      )
    }
  }
  bind_rows(out)
}

bacteria <- run_anova("path/to/16S_alpha_diversity.xlsx", "Bacteria")
fungi <- run_anova("path/to/ITS_alpha_diversity.xlsx", "Fungi")
all_res <- bind_rows(bacteria, fungi)

write.xlsx(all_res, "path/to/output/Two_Way_ANOVA_Long.xlsx", rowNames = FALSE)

wide <- all_res %>%
  select(Microbial_group, Index, Factor, F_value, p_value, Significance, Partial_Eta2) %>%
  pivot_wider(
    id_cols = c(Microbial_group, Index),
    names_from = Factor,
    values_from = c(F_value, p_value, Significance, Partial_Eta2),
    names_glue = "{Factor}_{.value}"
  )

write.xlsx(wide, "path/to/output/Two_Way_ANOVA_Wide.xlsx", rowNames = FALSE)