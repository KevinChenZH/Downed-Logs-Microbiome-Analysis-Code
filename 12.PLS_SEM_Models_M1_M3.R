# PLS-SEM for bacterial (M1: B1 guild) and fungal (M3: F3 guild) functional guilds
# Required input files:
#   1. M1 preprocessed CSV: Stage_III, Position, Cellulose, Hemicellulose, Lignin, B1
#   2. M3 preprocessed CSV: Stage_V, Position, Lignin, TN, TP, TK, F3

library(plspm)
library(ggplot2)
library(dplyr)

BASE_PATH <- "PREPROCESSED_DATA_DIRECTORY"

# --- Utility functions ---

total_effects_mat <- function(P) {
  n <- nrow(P)
  total <- P
  power <- P
  for (i in 2:(n - 1)) {
    power <- power %*% P
    total <- total + power
  }
  return(total)
}

calc_measurement_quality <- function(pls_obj, raw_data) {
  outer <- pls_obj$outer_model
  blocks <- pls_obj$model$blocks
  
  res <- data.frame(
    Latent = character(),
    Indicators = character(),
    Items = integer(),
    Cronbach_Alpha = numeric(),
    Composite_Reliability = numeric(),
    AVE = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (lv in names(blocks)) {
    items <- blocks[[lv]]
    n_items <- length(items)
    loadings <- outer$loading[outer$block == lv]
    
    ave <- mean(loadings^2, na.rm = TRUE)
    sum_load <- sum(loadings, na.rm = TRUE)
    sum_err <- sum(1 - loadings^2, na.rm = TRUE)
    cr <- ifelse(sum_load^2 + sum_err > 0, (sum_load^2) / (sum_load^2 + sum_err), NA)
    
    if (n_items == 1) {
      alpha <- 1.0000
    } else {
      dat <- raw_data[, items, drop = FALSE]
      covmat <- cov(dat, use = "pairwise.complete.obs")
      p <- ncol(covmat)
      col_var <- sum(diag(covmat))
      total_var <- sum(covmat)
      alpha <- ifelse(p > 1, (p / (p - 1)) * (1 - col_var / total_var), 1.0000)
    }
    
    res <- rbind(res, data.frame(
      Latent = lv,
      Indicators = paste(items, collapse = ", "),
      Items = n_items,
      Cronbach_Alpha = round(alpha, 4),
      Composite_Reliability = round(cr, 4),
      AVE = round(ave, 4),
      stringsAsFactors = FALSE
    ))
  }
  return(res)
}

save_full_results <- function(pls_obj, raw_data, model_name, out_dir) {
  
  cat("\n==========", model_name, "==========\n")
  
  # 1. GoF
  gof_val <- pls_obj$gof
  cat("GoF =", round(gof_val, 4), "\n")
  
  # 2. R2
  inner_sum <- pls_obj$inner_summary
  r2_df <- as.data.frame(inner_sum)
  r2_df$Construct <- rownames(r2_df)
  cat("\n--- R2 / Inner Summary ---\n")
  print(r2_df)
  write.csv(r2_df, file.path(out_dir, paste0(model_name, "_inner_summary.csv")), row.names = FALSE)
  
  # 3. Outer loadings
  outer_m <- pls_obj$outer_model
  cat("\n--- Outer Loadings ---\n")
  print(outer_m)
  write.csv(outer_m, file.path(out_dir, paste0(model_name, "_outer_loadings.csv")), row.names = FALSE)
  
  # 4. Measurement quality
  mq <- calc_measurement_quality(pls_obj, raw_data)
  cat("\n--- Measurement Quality ---\n")
  print(mq)
  write.csv(mq, file.path(out_dir, paste0(model_name, "_measurement_quality.csv")), row.names = FALSE)
  
  # 5. Path coefficients
  P <- pls_obj$path_coefs
  cat("\n--- Path Coefficients ---\n")
  print(P)
  write.csv(as.data.frame(P), file.path(out_dir, paste0(model_name, "_path_coefficients.csv")))
  
  # 6. Bootstrap
  boot_paths <- pls_obj$boot$paths
  boot_paths$z_value <- ifelse(boot_paths$Std.Error > 1e-10,
                               abs(boot_paths$Original / boot_paths$Std.Error), NA)
  boot_paths$p_value_approx <- ifelse(is.na(boot_paths$z_value), NA,
                                      2 * (1 - pnorm(boot_paths$z_value)))
  boot_paths$Stars <- ifelse(boot_paths$p_value_approx < 0.001, "***",
                             ifelse(boot_paths$p_value_approx < 0.01, "**",
                                    ifelse(boot_paths$p_value_approx < 0.05, "*", "ns")))
  boot_paths$CI_Significance <- ifelse(boot_paths$perc.025 > 0 | boot_paths$perc.975 < 0, "sig", "ns")
  
  cat("\n--- Bootstrap Paths (95% CI) ---\n")
  print(boot_paths)
  write.csv(boot_paths, file.path(out_dir, paste0(model_name, "_bootstrap_paths.csv")), row.names = FALSE)
  
  # 7. Total effects
  total_mat <- total_effects_mat(P)
  cat("\n--- Total Effects ---\n")
  print(total_mat)
  write.csv(as.data.frame(total_mat), file.path(out_dir, paste0(model_name, "_total_effects.csv")))
  
  # 8. Effect decomposition
  exogenous <- rownames(P)[rowSums(abs(P)) < 1e-10]
  endogenous  <- rownames(P)[rowSums(abs(P)) > 1e-10]
  
  cat("\n--- Effect Decomposition ---\n")
  cat("Exogenous:", paste(exogenous, collapse = ", "), "\n")
  cat("Endogenous:", paste(endogenous, collapse = ", "), "\n\n")
  
  decomp_list <- list()
  for (end in endogenous) {
    for (ex in exogenous) {
      direct <- P[end, ex]
      tot <- total_mat[end, ex]
      indirect <- tot - direct
      decomp_list[[length(decomp_list) + 1]] <- data.frame(
        Endogenous = end,
        Exogenous = ex,
        Direct_Effect = round(direct, 4),
        Indirect_Effect = round(indirect, 4),
        Total_Effect = round(tot, 4),
        stringsAsFactors = FALSE
      )
    }
  }
  decomp_df <- do.call(rbind, decomp_list)
  print(decomp_df)
  write.csv(decomp_df, file.path(out_dir, paste0(model_name, "_effect_decomposition.csv")), row.names = FALSE)
  
  return(list(
    total = total_mat,
    boot = boot_paths,
    gof = gof_val,
    inner = r2_df,
    outer = outer_m,
    measurement = mq,
    decomposition = decomp_df
  ))
}

plot_total_effects_bar <- function(total_mat, target_name, model_label, out_path, gof_val) {
  target_idx <- which(rownames(total_mat) == target_name)
  effs <- total_mat[target_idx, ]
  effs <- effs[abs(effs) > 1e-10]
  if (length(effs) == 0) return(NULL)
  
  df <- data.frame(
    Predictor = names(effs),
    Effect = as.numeric(effs),
    Direction = ifelse(effs > 0, "Positive", "Negative")
  )
  
  p <- ggplot(df, aes(x = reorder(Predictor, -Effect), y = Effect, fill = Direction)) +
    geom_bar(stat = "identity", width = 0.6, color = "black", linewidth = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey30", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.3f", Effect)), 
              vjust = ifelse(df$Effect > 0, -0.6, 1.3),
              size = 4.2, fontface = "bold", color = "grey10") +
    scale_y_continuous(expand = expansion(mult = c(0.22, 0.25))) +
    scale_fill_manual(values = c("Positive" = "#08519C", "Negative" = "#C6DBEF")) +
    labs(title = paste0(model_label, ": Total Effects on ", target_name),
         subtitle = paste("GoF =", round(gof_val, 3), "| PLS-PM with 1000 bootstrap replicates"),
         x = "Predictor Construct", 
         y = "Total Effect (Standardized)") +
    theme_bw(base_size = 14) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 15, margin = margin(b = 8)),
          plot.subtitle = element_text(size = 12, color = "grey30"),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 25, hjust = 1, size = 11),
          plot.margin = margin(20, 40, 20, 20))
  
  ggsave(out_path, plot = p, device = "pdf", width = 9, height = 7, dpi = 600)
  return(p)
}

# --- Model 1: Bacterial B1-Guild (5 latent variables, direct paths included) ---

cat("\n========================================\n")
cat("  Running M1: Bacterial B1-Guild\n")
cat("========================================\n")

dat1 <- read.csv(file.path(BASE_PATH, "sem_model1_bacteria_carbon.csv"), header = TRUE)
dat1 <- na.omit(dat1)
cat("  M1 n =", nrow(dat1), "\n")

latent_m1 <- c("Stage_III", "Position", "Carbon_Polymers", "Lignin", "B1_Guild")

blocks_m1 <- list(
  c("Stage_III"),
  c("Position"),
  c("Cellulose", "Hemicellulose"),
  c("Lignin"),
  c("B1")
)

path_m1 <- matrix(c(
  0, 0, 0, 0, 0,
  0, 0, 0, 0, 0,
  1, 1, 0, 0, 0,
  1, 1, 0, 0, 0,
  1, 1, 1, 1, 0
), nrow = 5, ncol = 5, byrow = TRUE)
dimnames(path_m1) <- list(latent_m1, latent_m1)

pls_m1 <- plspm(dat1, path_matrix = path_m1, blocks = blocks_m1,
                modes = rep("A", 5), boot.val = TRUE, br = 1000)

res_m1 <- save_full_results(pls_m1, dat1, "M1_Bacterial_B1Only", BASE_PATH)

plot_total_effects_bar(res_m1$total, "B1_Guild",
                       "Model 1 (Bacterial B1-Guild)",
                       file.path(BASE_PATH, "Fig5_M1_B1_total_effects.pdf"),
                       gof_val = res_m1$gof)

cat("\n  M1 bar plot saved: Fig5_M1_B1_total_effects.pdf\n")

# --- Model 3: Fungal F3-Guild (5 latent variables, full mediation) ---

cat("\n========================================\n")
cat("  Running M3: Fungal F3-Guild\n")
cat("========================================\n")

dat3 <- read.csv(file.path(BASE_PATH, "sem_model3_fungi_pollution.csv"), header = TRUE)
dat3 <- na.omit(dat3)
cat("  M3 n =", nrow(dat3), "\n")

latent_m3 <- c("Stage_V", "Position", "Lignin", "Nutrient", "F3_Guild")

blocks_m3 <- list(
  c("Stage_V"),
  c("Position"),
  c("Lignin"),
  c("TN", "TP", "TK"),
  c("F3")
)

path_m3 <- matrix(c(
  0, 0, 0, 0, 0,
  0, 0, 0, 0, 0,
  1, 1, 0, 0, 0,
  1, 1, 0, 0, 0,
  0, 0, 1, 1, 0
), nrow = 5, ncol = 5, byrow = TRUE)
dimnames(path_m3) <- list(latent_m3, latent_m3)

pls_m3 <- plspm(dat3, path_matrix = path_m3, blocks = blocks_m3,
                modes = rep("A", 5), boot.val = TRUE, br = 1000)

res_m3 <- save_full_results(pls_m3, dat3, "M3_Fungal_F3", BASE_PATH)

plot_total_effects_bar(res_m3$total, "F3_Guild",
                       "Model 3 (Fungal F3-Guild)",
                       file.path(BASE_PATH, "Fig5_M3_F3_total_effects.pdf"),
                       gof_val = res_m3$gof)

cat("\n  M3 bar plot saved: Fig5_M3_F3_total_effects.pdf\n")

cat("\n========================================\n")
cat("  M1 + M3 execution completed\n")
cat("========================================\n")
cat("\nOutput files saved in:", BASE_PATH, "\n")
cat("M1: inner_summary, outer_loadings, measurement_quality, path_coefficients, bootstrap_paths, total_effects, effect_decomposition, total_effects_bar.pdf\n")
cat("M3: inner_summary, outer_loadings, measurement_quality, path_coefficients, bootstrap_paths, total_effects, effect_decomposition, total_effects_bar.pdf\n")