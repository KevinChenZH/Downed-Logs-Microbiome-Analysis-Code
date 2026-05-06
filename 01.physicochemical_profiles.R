library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)

# Expected columns: Raw_ID, Mapped_ID, Group, Dry_Matter, Ash, Total_N,
# Crude_Protein, Total_P, Total_K, Cellulose, Lignin, Hemicellulose, Position, Decay_Level

load_data <- function(file_path) {
  df <- read_excel(file_path, sheet = 1)
  df <- df %>% filter(Raw_ID != "Raw_ID")
  
  nums <- c("Dry_Matter", "Ash", "Total_N", "Crude_Protein",
            "Total_P", "Total_K", "Cellulose", "Lignin", "Hemicellulose")
  df <- df %>% mutate(across(all_of(nums), as.numeric))
  
  df %>%
    mutate(
      Position = case_when(trimws(Position) == "N" ~ "Internal",
                           trimws(Position) == "W" ~ "External"),
      Decay_Level = trimws(Decay_Level),
      Decay_Level = ifelse(Decay_Level %in% c("WV", "V"), "V", Decay_Level),
      Decay_Level = factor(Decay_Level, levels = c("I", "III", "V"), ordered = TRUE)
    )
}

plot_one <- function(df, var, ylab, title) {
  smry <- df %>%
    group_by(Position, Decay_Level) %>%
    summarise(
      mean = mean(.data[[var]], na.rm = TRUE),
      sd = sd(.data[[var]], na.rm = TRUE),
      n = sum(!is.na(.data[[var]])),
      .groups = "drop"
    ) %>%
    mutate(sem = sd / sqrt(n))
  
  ggplot(smry, aes(x = Decay_Level, y = mean, group = Position, color = Position)) +
    geom_line(linewidth = 0.8) +
    geom_point(aes(shape = Position), size = 2.5) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.15, linewidth = 0.6) +
    scale_color_manual(values = c(Internal = "#0077BB", External = "#EE7733")) +
    scale_shape_manual(values = c(Internal = 16, External = 17)) +
    labs(title = title, x = "Decay Level", y = ylab) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      axis.title = element_text(face = "bold", size = 10),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linetype = "dotted")
    )
}

main <- function() {
  df <- load_data("path/to/your/physicochemical_data.xlsx")
  
  p1 <- plot_one(df, "Dry_Matter", "Dry Matter (%)", "Dry Matter Content")
  p2 <- plot_one(df, "Total_N", "Total N (g/kg)", "Total Nitrogen")
  p3 <- plot_one(df, "Total_P", "Total P (g/kg)", "Total Phosphorus")
  p4 <- plot_one(df, "Lignin", "Lignin (%)", "Lignin Content")
  p5 <- plot_one(df, "Cellulose", "Cellulose (%)", "Cellulose Content")
  p6 <- plot_one(df, "Hemicellulose", "Hemicellulose (%)", "Hemicellulose Content")
  
  p1 <- p1 + theme(
    legend.position = "inside",
    legend.position.inside = c(0.02, 0.98),
    legend.justification = c("left", "top"),
    legend.title = element_blank()
  )
  
  out <- wrap_plots(p1, p2, p3, p4, p5, p6, ncol = 3)
  ggsave("Fig1_Physicochemical_Profiles.pdf", out, width = 16, height = 10, dpi = 600)
}

main()
#Notes:
#Ensure your input Excel uses the English column names listed in the header comment.
#patchwork is required for panel assembly; install via install.packages("patchwork") if absent.