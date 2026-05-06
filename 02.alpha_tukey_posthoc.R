library(readxl)
library(dplyr)
library(openxlsx)
library(agricolae)

METRICS <- c(
  chao1 = "Chao1 index",
  observed_features = "Observed ASVs",
  shannon = "Shannon index",
  pielou_e = "Pielou's evenness"
)

parse_name <- function(x) {
  if (is.na(x)) return(c(NA_character_, NA_character_))
  x <- trimws(as.character(x))
  sm <- c("1" = "I", "3" = "III", "5" = "V")
  p <- substr(x, 1, 1)
  s <- sm[substr(x, 2, 2)]
  if (p == "n") return(c(s, "Internal"))
  if (p == "w") return(c(s, "External"))
  c(NA_character_, NA_character_)
}

run_tukey <- function(file_path, group_label) {
  df <- read_excel(file_path, col_names = TRUE)
  pr <- t(vapply(df[[1]], parse_name, character(2)))
  df$Decay_stage <- factor(pr[, 1], levels = c("I", "III", "V"))
  df$Position <- factor(pr[, 2], levels = c("Internal", "External"))
  df <- df %>% filter(!is.na(Decay_stage), !is.na(Position))
  
  cols <- intersect(names(METRICS), names(df))
  df <- df %>% mutate(across(all_of(cols), as.numeric))
  df$Group <- paste(df$Decay_stage, df$Position, sep = "-")
  
  out <- df %>% distinct(Decay_stage, Position) %>% mutate(Microbial_group = group_label)
  
  for (c in cols) {
    sub <- df %>% filter(!is.na(.data[[c]]))
    if (nrow(sub) < 3 || n_distinct(sub$Group) < 2) {
      out[[METRICS[c]]] <- NA
      next
    }
    
    g <- sub %>%
      group_by(Decay_stage, Position, Group) %>%
      summarise(Mean = mean(.data[[c]]), SE = sd(.data[[c]]) / sqrt(n()), .groups = "drop")
    
    m <- aov(as.formula(paste(c, "~ Group")), data = sub)
    h <- HSD.test(m, "Group", alpha = 0.05, console = FALSE)
    l <- data.frame(Group = rownames(h$groups), Letter = h$groups$groups, stringsAsFactors = FALSE)
    g <- g %>%
      left_join(l, by = "Group") %>%
      mutate(!!METRICS[c] := sprintf("%.2f±%.2f<sup>%s</sup>", Mean, SE, Letter)) %>%
      select(Decay_stage, Position, !!METRICS[c])
    
    out <- out %>% left_join(g, by = c("Decay_stage", "Position"))
  }
  
  out %>%
    select(Microbial_group, Decay_stage, Position, all_of(METRICS[cols])) %>%
    arrange(Microbial_group, Decay_stage, Position)
}

bacteria <- run_tukey("path/to/16S_alpha_diversity.xlsx", "Bacteria")
fungi <- run_tukey("path/to/ITS_alpha_diversity.xlsx", "Fungi")

bind_rows(bacteria, fungi) %>%
  write.xlsx("path/to/output/Alpha_Diversity_Tukey.xlsx", rowNames = FALSE)
#Notes:
#Requires agricolae for Tukey HSD and compact letter display (CLD). Install via install.packages("agricolae") if absent.
#Input Excel files must contain the columns chao1, observed_features, shannon, pielou_e with sample IDs in the first column (e.g., n1.1, w3.2).
#Output strings include HTML <sup> tags for manual superscript formatting in Word.
