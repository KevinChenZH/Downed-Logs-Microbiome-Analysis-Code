# Cross-kingdom co-occurrence network analysis for downed log microbiome
# Required input files:
#   1. 16S feature table in BIOM format
#   2. ITS feature table in BIOM format
#   3. Metadata CSV with columns: SampleID, Stage, Position
#   4. 16S taxonomy TSV (QIIME2 format)
#   5. ITS taxonomy TSV (QIIME2 format)
#   6. Physicochemistry summary CSV (optional, for III-Internal null explanation)

suppressPackageStartupMessages({
  library(SpiecEasi)
  library(biomformat)
  library(igraph)
  library(tidyverse)
})

# --- Configuration ---
BIOM_16S_PATH      <- "BIOM_16S_PATH"
BIOM_ITS_PATH      <- "BIOM_ITS_PATH"
METADATA_PATH      <- "METADATA_PATH"
TAXONOMY_16S_PATH  <- "TAXONOMY_16S_PATH"
TAXONOMY_ITS_PATH  <- "TAXONOMY_ITS_PATH"
PHYSICOCHEM_PATH   <- "PHYSICOCHEM_PATH"
OUTPUT_DIR         <- "OUTPUT_DIRECTORY"

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Utility functions ---

read_biom_to_matrix <- function(biom_path) {
  biom_obj <- read_biom(biom_path)
  as(biom_data(biom_obj), "matrix")
}

parse_tax <- function(tax_df) {
  tax_df$Taxon <- as.character(tax_df$Taxon)
  tax_df$Genus <- sapply(strsplit(tax_df$Taxon, ";"), function(x) {
    g <- grep("^g__", x, value = TRUE)
    if (length(g) > 0) {
      genus <- sub("^g__", "", g[1])
      if (genus == "" || genus == "NA" || grepl("Incertae_sedis|unclassified", genus, ignore.case = TRUE)) return(NA)
      return(genus)
    } else return(NA)
  })
  tax_df$Phylum <- sapply(strsplit(tax_df$Taxon, ";"), function(x) {
    p <- grep("^p__", x, value = TRUE)
    if (length(p) > 0) sub("^p__", "", p[1]) else "Unclassified"
  })
  return(tax_df)
}

aggregate_to_genus <- function(otu_mat, tax_df) {
  common_asvs <- intersect(rownames(otu_mat), tax_df[[1]])
  otu_mat <- otu_mat[common_asvs, , drop = FALSE]
  tax_df <- tax_df[match(common_asvs, tax_df[[1]]), ]
  valid_idx <- !is.na(tax_df$Genus)
  otu_mat <- otu_mat[valid_idx, , drop = FALSE]
  tax_df <- tax_df[valid_idx, ]
  genus_list <- split(seq_len(nrow(otu_mat)), tax_df$Genus)
  genus_mat <- do.call(rbind, lapply(names(genus_list), function(g) {
    colSums(otu_mat[genus_list[[g]], , drop = FALSE])
  }))
  rownames(genus_mat) <- names(genus_list)
  genus_tax <- tax_df[!duplicated(tax_df$Genus), c("Genus", "Phylum")]
  rownames(genus_tax) <- genus_tax$Genus
  return(list(otu = genus_mat, tax = genus_tax[rownames(genus_mat), ]))
}

filter_group <- function(otu_mat, samples, prev_ratio = 0.5, min_abund = 0.001) {
  sub_mat <- otu_mat[, samples, drop = FALSE]
  n <- ncol(sub_mat)
  min_prev <- ceiling(n * prev_ratio)
  prevalence <- rowSums(sub_mat > 0)
  rel_abund <- rowSums(sub_mat) / sum(sub_mat)
  keep <- prevalence >= min_prev & rel_abund > min_abund
  return(sub_mat[keep, , drop = FALSE])
}

calculate_zi_pi <- function(g, modules) {
  membership <- modules$membership
  node_stats <- data.frame(
    Node = V(g)$name,
    Label = sub("^[BF]_", "", V(g)$name),
    Kingdom = ifelse(grepl("^B_", V(g)$name), "Bacteria", "Fungi"),
    Module = membership,
    TotalDegree = degree(g),
    stringsAsFactors = FALSE
  )
  
  module_stats <- data.frame(
    Module = 1:length(modules),
    MeanK = sapply(1:length(modules), function(m) mean(degree(g)[membership == m])),
    SDK = sapply(1:length(modules), function(m) sd(degree(g)[membership == m]))
  )
  
  node_stats$Zi <- sapply(1:vcount(g), function(i) {
    m <- membership[i]
    ki <- sum(membership[neighbors(g, i)] == m)
    k_mean <- module_stats$MeanK[m]
    k_sd <- module_stats$SDK[m]
    if (k_sd == 0 || is.na(k_sd)) return(0)
    return((ki - k_mean) / k_sd)
  })
  
  node_stats$Pi <- sapply(1:vcount(g), function(i) {
    m <- membership[i]
    ki <- degree(g)[i]
    if (ki == 0) return(0)
    neighbor_modules <- membership[neighbors(g, i)]
    module_connects <- table(factor(neighbor_modules, levels = unique(membership)))
    return(1 - sum((module_connects / ki)^2))
  })
  
  node_stats$Role <- "Peripheral"
  node_stats$Role[node_stats$Zi > 2.5 & node_stats$Pi <= 0.62] <- "Module Hub"
  node_stats$Role[node_stats$Zi <= 2.5 & node_stats$Pi > 0.62] <- "Connector"
  node_stats$Role[node_stats$Zi > 2.5 & node_stats$Pi > 0.62] <- "Network Hub"
  
  key_genera <- c("Burkholderia", "Caballeronia", "Paraburkholderia", "Pseudomonas",
                  "Granulicella", "Gryllotalpicola", "Luteibacter", "Acidisoma",
                  "Botryobasidium", "Dacryopinax", "Tricholomopsis", "Resinicium",
                  "Lentinellus", "Mortierella", "Cladophialophora", "Mycena",
                  "Lasiosphaeris", "Chaetosphaeria", "Stropharia", "Phallus")
  node_stats$Is_Key <- node_stats$Label %in% key_genera
  
  return(node_stats)
}

build_network <- function(otu_mat, group_name) {
  if (nrow(otu_mat) < 5) return(NULL)
  
  se <- spiec.easi(t(otu_mat), method = "mb", lambda.min.ratio = 0.1, nlambda = 30,
                   sel.criterion = "bstars",
                   pulsar.params = list(rep.num = 100, seed = 2401, ncores = 1, lb.stars = FALSE),
                   verbose = FALSE)
  
  adj <- getRefit(se)
  rownames(adj) <- colnames(adj) <- rownames(otu_mat)
  if (sum(adj) == 0) return(NULL)
  
  g <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
  
  cor_mat <- cor(t(otu_mat), method = "spearman")
  cor_mat[is.na(cor_mat)] <- 0
  edge_list <- as_edgelist(g)
  weights <- apply(edge_list, 1, function(e) cor_mat[e[1], e[2]])
  
  abs_w <- abs(weights)
  w_scaled <- ifelse(max(abs_w) > min(abs_w), (abs_w - min(abs_w)) / (max(abs_w) - min(abs_w)), rep(1, length(weights)))
  
  E(g)$Weight <- w_scaled
  E(g)$Spearman <- weights
  E(g)$Sign <- ifelse(weights > 0, "Positive", "Negative")
  
  modules <- cluster_louvain(g)
  modularity_val <- modularity(modules)
  zi_pi_df <- calculate_zi_pi(g, modules)
  
  edges_df <- data.frame(
    Source = edge_list[, 1],
    Target = edge_list[, 2],
    Sign = E(g)$Sign,
    stringsAsFactors = FALSE
  )
  edges_df$Type <- ifelse(
    substr(edges_df$Source, 1, 1) == substr(edges_df$Target, 1, 1),
    ifelse(substr(edges_df$Source, 1, 1) == "B", "Bac-Bac", "Fun-Fun"),
    "Bac-Fun"
  )
  
  edge_stats <- edges_df %>%
    group_by(Type, Sign) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = Sign, values_from = n, values_fill = 0) %>%
    mutate(Total = Positive + Negative, PosRatio = round(Positive / Total, 3))
  
  return(list(graph = g, zi_pi = zi_pi_df, modules = modules, modularity = modularity_val,
              edge_stats = edge_stats, edges_df = edges_df))
}

# --- Functional guild constants ---
GUILD_MAP <- list(
  B1 = c("Granulicella", "Burkholderia-Caballeronia-Paraburkholderia",
         "Burkholderia", "Caballeronia", "Paraburkholderia",
         "Gryllotalpicola", "Luteibacter"),
  B3 = c("Pseudomonas"),
  F3 = c("Mortierella", "Cladophialophora")
)

# --- Main workflow ---

main <- function() {
  otu_16s <- read_biom_to_matrix(BIOM_16S_PATH)
  otu_its <- read_biom_to_matrix(BIOM_ITS_PATH)
  metadata <- read.csv(METADATA_PATH, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(metadata) <- c("SampleID", "Stage", "Position", "Group")
  
  # Optional physicochemistry for III-Internal null explanation
  tn_iii_int <- 1.84; tp_iii_int <- 0.1
  if (file.exists(PHYSICOCHEM_PATH)) {
    pc <- read.csv(PHYSICOCHEM_PATH, stringsAsFactors = FALSE, skip = 2, header = FALSE)
    colnames(pc) <- c("Variable", "Ext_I", "Ext_III", "Ext_V", "Int_I", "Int_III", "Int_V")
    get_val <- function(var_name, col_name) {
      val <- pc %>% filter(Variable == var_name) %>% pull(!!sym(col_name))
      as.numeric(strsplit(as.character(val), "±")[[1]][1])
    }
    tn_iii_int <- get_val("Total_N", "Int_III")
    tp_iii_int <- get_val("Total_P", "Int_III")
  }
  
  tax_16s <- read.delim(TAXONOMY_16S_PATH, stringsAsFactors = FALSE, check.names = FALSE)
  if (tax_16s[1, 1] == "#q2:types") tax_16s <- tax_16s[-1, ]
  tax_its <- read.delim(TAXONOMY_ITS_PATH, stringsAsFactors = FALSE, check.names = FALSE)
  if (tax_its[1, 1] == "#q2:types") tax_its <- tax_its[-1, ]
  
  tax_16s <- parse_tax(tax_16s)
  tax_its <- parse_tax(tax_its)
  
  agg_16s <- aggregate_to_genus(otu_16s, tax_16s)
  agg_its <- aggregate_to_genus(otu_its, tax_its)
  
  get_samples <- function(stage, position) {
    metadata$SampleID[metadata$Stage == stage & metadata$Position == position]
  }
  
  cross_groups <- list(
    list(name = "Cross_III_int", stage = "III", position = "Internal"),
    list(name = "Cross_III_ext", stage = "III", position = "External"),
    list(name = "Cross_V_int", stage = "V", position = "Internal"),
    list(name = "Cross_V_ext", stage = "V", position = "External")
  )
  
  network_summary <- list()
  connectors_list <- list()
  edge_summary <- list()
  guild_node_list <- list()
  guild_edge_list <- list()
  guild_summary_list <- list()
  
  for (grp in cross_groups) {
    samples <- get_samples(grp$stage, grp$position)
    bac <- filter_group(agg_16s$otu, samples)
    fun <- filter_group(agg_its$otu, samples)
    
    if (nrow(bac) == 0 || nrow(fun) == 0 || (nrow(bac) + nrow(fun)) < 5) {
      next  # III-Internal skipped due to insufficient taxa after filtering
    }
    
    rownames(bac) <- paste0("B_", rownames(bac))
    rownames(fun) <- paste0("F_", rownames(fun))
    bac_tax <- agg_16s$tax[sub("^B_", "", rownames(bac)), ]
    rownames(bac_tax) <- paste0("B_", rownames(bac_tax))
    fun_tax <- agg_its$tax[sub("^F_", "", rownames(fun)), ]
    rownames(fun_tax) <- paste0("F_", rownames(fun_tax))
    bac_tax$Kingdom <- "Bacteria"; fun_tax$Kingdom <- "Fungi"
    
    otu <- rbind(bac, fun)
    res <- build_network(otu, grp$name)
    if (is.null(res)) next
    
    # Network summary
    conn <- res$zi_pi %>% filter(Role == "Connector")
    if (nrow(conn) > 0) {
      conn$Network <- grp$name; conn$Stage <- grp$stage; conn$Position <- grp$position
      connectors_list[[grp$name]] <- conn
    }
    
    es <- res$edge_stats; es$Network <- grp$name
    edge_summary[[grp$name]] <- es
    
    network_summary[[grp$name]] <- data.frame(
      Network = grp$name, Stage = grp$stage, Position = grp$position,
      N_Bacteria = nrow(bac), N_Fungi = nrow(fun),
      Nodes = vcount(res$graph), Edges = ecount(res$graph),
      Modularity = round(res$modularity, 3),
      NegRatio = round(sum(res$edges_df$Sign == "Negative") / ecount(res$graph), 3),
      Connectors = nrow(conn),
      stringsAsFactors = FALSE
    )
    
    # Guild subnetwork analysis
    g <- res$graph; zi_pi_df <- res$zi_pi; edges_df <- res$edges_df
    all_node_names <- V(g)$name
    all_node_labels <- sub("^[BF]_", "", all_node_names)
    all_node_guilds <- rep(NA_character_, length(all_node_names))
    for (guild_name in names(GUILD_MAP)) {
      match_mask <- all_node_labels %in% GUILD_MAP[[guild_name]]
      all_node_guilds[match_mask] <- guild_name
    }
    
    guild_node_mask <- !is.na(all_node_guilds)
    guild_node_names <- all_node_names[guild_node_mask]
    
    if (length(guild_node_names) > 0) {
      otu_rel <- apply(otu, 2, function(x) x / sum(x))
      node_mean_abund <- rowMeans(otu_rel)[guild_node_names]
      
      phylum_vec <- rep(NA_character_, length(guild_node_names))
      names(phylum_vec) <- guild_node_names
      b_mask <- grepl("^B_", guild_node_names)
      phylum_vec[b_mask] <- as.character(bac_tax[guild_node_names[b_mask], "Phylum"])
      phylum_vec[!b_mask] <- as.character(fun_tax[guild_node_names[!b_mask], "Phylum"])
      
      gnode_df <- data.frame(
        Network = grp$name, Stage = grp$stage, Position = grp$position,
        Node = guild_node_names, Label = all_node_labels[guild_node_mask],
        Kingdom = ifelse(grepl("^B_", guild_node_names), "Bacteria", "Fungi"),
        Phylum = phylum_vec, Guild = all_node_guilds[guild_node_mask],
        Mean_Rel_Abundance = as.numeric(node_mean_abund),
        Degree = zi_pi_df$TotalDegree[match(guild_node_names, zi_pi_df$Node)],
        Zi = zi_pi_df$Zi[match(guild_node_names, zi_pi_df$Node)],
        Pi = zi_pi_df$Pi[match(guild_node_names, zi_pi_df$Node)],
        Role = zi_pi_df$Role[match(guild_node_names, zi_pi_df$Node)],
        Is_Key = zi_pi_df$Is_Key[match(guild_node_names, zi_pi_df$Node)],
        stringsAsFactors = FALSE
      )
      guild_node_list[[grp$name]] <- gnode_df
      
      edges_df$Spearman <- E(g)$Spearman
      sub_edges <- edges_df %>%
        filter(Source %in% guild_node_names, Target %in% guild_node_names) %>%
        mutate(
          Network = grp$name, Stage = grp$stage, Position = grp$position,
          Source_Guild = all_node_guilds[match(Source, all_node_names)],
          Target_Guild = all_node_guilds[match(Target, all_node_names)],
          Source_Label = sub("^[BF]_", "", Source),
          Target_Label = sub("^[BF]_", "", Target),
          Guild_Pair = paste0(pmin(Source_Guild, Target_Guild), "-", pmax(Source_Guild, Target_Guild)),
          Cross_Kingdom = substr(Source, 1, 1) != substr(Target, 1, 1)
        ) %>%
        select(Network, Stage, Position, Source, Target, Source_Label, Target_Label,
               Source_Guild, Target_Guild, Guild_Pair, Cross_Kingdom, Sign, Spearman, Edge_Type = Type)
      
      guild_edge_list[[grp$name]] <- sub_edges
      
      if (nrow(sub_edges) > 0) {
        pair_summary <- sub_edges %>%
          group_by(Guild_Pair) %>%
          summarise(
            N_Edges = n(), N_Positive = sum(Sign == "Positive"),
            N_Negative = sum(Sign == "Negative"),
            Pos_Ratio = round(mean(Sign == "Positive"), 3),
            Cross_Kingdom = first(Cross_Kingdom), .groups = "drop"
          ) %>%
          mutate(Network = grp$name, Stage = grp$stage, Position = grp$position) %>%
          select(Network, Stage, Position, Guild_Pair, Cross_Kingdom, everything())
        guild_summary_list[[grp$name]] <- pair_summary
      } else {
        guild_summary_list[[grp$name]] <- data.frame(
          Network = character(0), Stage = character(0), Position = character(0),
          Guild_Pair = character(0), Cross_Kingdom = logical(0), N_Edges = integer(0),
          N_Positive = integer(0), N_Negative = integer(0), Pos_Ratio = numeric(0),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  # Export results
  if (length(network_summary) > 0) {
    write.csv(do.call(rbind, network_summary), file.path(OUTPUT_DIR, "01_Network_Summary.csv"), row.names = FALSE)
  }
  if (length(connectors_list) > 0) {
    cl <- do.call(rbind, connectors_list)
    write.csv(cl, file.path(OUTPUT_DIR, "02_All_Connectors.csv"), row.names = FALSE)
    key_cl <- cl %>% filter(Is_Key == TRUE)
    if (nrow(key_cl) > 0) {
      write.csv(key_cl, file.path(OUTPUT_DIR, "03_Key_Genera_Connectors.csv"), row.names = FALSE)
    }
  }
  if (length(edge_summary) > 0) {
    write.csv(do.call(rbind, edge_summary), file.path(OUTPUT_DIR, "04_Edge_Type_Statistics.csv"), row.names = FALSE)
  }
  
  guild_node_list <- guild_node_list[!sapply(guild_node_list, is.null)]
  guild_edge_list <- guild_edge_list[!sapply(guild_edge_list, is.null)]
  guild_summary_list <- guild_summary_list[!sapply(guild_summary_list, is.null)]
  
  if (length(guild_node_list) > 0) {
    write.csv(do.call(rbind, guild_node_list), file.path(OUTPUT_DIR, "05_Guild_Node_Topology.csv"), row.names = FALSE)
    write.csv(do.call(rbind, guild_edge_list), file.path(OUTPUT_DIR, "06_Guild_Subnetwork_Edges.csv"), row.names = FALSE)
    if (length(guild_summary_list) > 0) {
      write.csv(do.call(rbind, guild_summary_list), file.path(OUTPUT_DIR, "07_Guild_Interaction_Summary.csv"), row.names = FALSE)
    }
    
    # Writing-ready wide-format table
    writing_rows <- list()
    for (nm in names(guild_node_list)) {
      ns <- guild_node_list[[nm]]
      es <- if (nm %in% names(guild_edge_list)) guild_edge_list[[nm]] else NULL
      net_info <- network_summary[[nm]]
      
      row <- data.frame(Network = nm, Stage = ns$Stage[1], Position = ns$Position[1], stringsAsFactors = FALSE)
      for (g in c("B1", "B3", "F3")) row[[paste0(g, "_Nodes")]] <- sum(ns$Guild == g)
      row$Total_Guild_Nodes <- nrow(ns)
      
      if (!is.null(es) && nrow(es) > 0) {
        row$Total_Guild_Edges <- nrow(es)
        row$Guild_Positive_Edges <- sum(es$Sign == "Positive")
        row$Guild_Negative_Edges <- sum(es$Sign == "Negative")
        row$Guild_Pos_Ratio <- round(mean(es$Sign == "Positive"), 3)
        row$Cross_Kingdom_Edges <- sum(es$Cross_Kingdom)
        row$Cross_Kingdom_Positive <- sum(es$Cross_Kingdom & es$Sign == "Positive")
        row$Pct_of_Total_Edges <- round(nrow(es) / net_info$Edges * 100, 1)
        pairs <- c("B1-B1", "B1-B3", "B1-F3", "B3-B3", "B3-F3", "F3-F3")
        for (p in pairs) {
          p_edges <- es %>% filter(Guild_Pair == p)
          row[[paste0(p, "_Edges")]] <- nrow(p_edges)
          row[[paste0(p, "_Pos")]] <- sum(p_edges$Sign == "Positive")
          row[[paste0(p, "_Neg")]] <- sum(p_edges$Sign == "Negative")
        }
      } else {
        row$Total_Guild_Edges <- 0; row$Guild_Positive_Edges <- 0; row$Guild_Negative_Edges <- 0
        row$Guild_Pos_Ratio <- NA; row$Cross_Kingdom_Edges <- 0; row$Cross_Kingdom_Positive <- 0
        row$Pct_of_Total_Edges <- 0
        for (p in c("B1-B1", "B1-B3", "B1-F3", "B3-B3", "B3-F3", "F3-F3")) {
          row[[paste0(p, "_Edges")]] <- 0; row[[paste0(p, "_Pos")]] <- 0; row[[paste0(p, "_Neg")]] <- 0
        }
      }
      writing_rows[[nm]] <- row
    }
    write.csv(do.call(rbind, writing_rows), file.path(OUTPUT_DIR, "08_Guild_Writing_Stats.csv"), row.names = FALSE)
  }
  
  cat("Analysis completed. Output directory:", OUTPUT_DIR, "\n")
  cat("III-Internal network null: TN =", tn_iii_int, ", TP =", tp_iii_int, "(extreme oligotrophy)\n")
}

main()