if(!require("pacman")) {
  install.packages("pacman", repos = "http://cran.us.r-project.org")
}
pacman::p_load("ade4",
               "GUniFrac",
               "phangorn",
               "cluster",
               "fpc",
               "vegan",
               "clusterSim")

# Calculate the species richness in a sample
Species.richness <- function(x)
{
  # Count only the OTUs that are present >0.5 normalized counts (normalization produces real values for counts)
  count = sum(x[x > 0.5] ^ 0)
  return(count)
}

# Calculate the Effective species richness in each individual sample
Eff.Species.richness <- function(x, eff.cutoff = 0.0025)
{
  # Count only the OTUs that are present more than the set proportion
  total = sum(x)
  count = sum(x[x / total > eff.cutoff] ^ 0)
  return(count)
}

# Calculate the Normalized species richness in each individual sample
Norm.Species.richness <- function(x, norm.cutoff = 1000)
{
  # Count only the OTUs that are present >0.5 normalized counts (normalization produces real values for counts)
  # Given a fixed Normalization reads depth
  total = sum(x)
  count = sum(x[norm.cutoff * x / total > 0.5] ^ 0)
  return(count)
}


# Calculate the Shannon diversity index
Shannon.entropy <- function(x)
{
  total = sum(x)
  se = -sum(x[x > 0] / total * log(x[x > 0] / total))
  return(se)
}

# Calculate the effective number of species for Shannon
Shannon.effective <- function(x)
{
  total = sum(x)
  se = round(exp(-sum(x[x > 0] / total * log(x[x > 0] / total))), digits =
               2)
  return(se)
}

# Calculate the Simpson diversity index
Simpson.concentration <- function(x)
{
  total = sum(x)
  si = sum((x[x > 0] / total) ^ 2)
  return(si)
}

# Calculate the effective number of species for Simpson
Simpson.effective <- function(x)
{
  total = sum(x)
  si = round(1 / sum((x[x > 0] / total) ^ 2), digits = 2)
  return(si)
}

Rhea.Alpha_Diversity <- function(ps,
                                 eff.cutoff = 0.0025,
                                 norm.cutoff = 1000) {
  meta_file <- meta_to_df(ps)
  meta_file <- data.frame(meta_file[!apply(is.na(meta_file) |
                                             meta_file == "", 1, all), , drop = FALSE])
  meta_file <- data.frame(meta_file[order(row.names(meta_file)), , drop = FALSE])
  
  # Clean table from empty lines
  otu_table <- ps_to_asvtab(ps)
  otu_table <- otu_table[!apply(is.na(otu_table) |
                                  otu_table == "", 1, all), ]
  
  # Order and transpose OTU-table
  my_otu_table <- otu_table[, order(names(otu_table))]
  my_otu_table <- data.frame(t(my_otu_table))
  
  # Apply diversity functions to table
  otus_div_stats <- data.frame(my_otu_table[, 0])
  otus_div_stats$Richness <- apply(my_otu_table, 1, Species.richness)
  otus_div_stats$Normalized.Richness <- apply(my_otu_table, 1, Norm.Species.richness)
  otus_div_stats$Effective.Richness <- apply(my_otu_table, 1, Eff.Species.richness)
  otus_div_stats$Shannon.Index <- apply(my_otu_table, 1, Shannon.entropy)
  otus_div_stats$Shannon.Effective <- apply(my_otu_table, 1, Shannon.effective)
  otus_div_stats$Simpson.Index <- apply(my_otu_table, 1, Simpson.concentration)
  otus_div_stats$Simpson.Effective <- apply(my_otu_table, 1, Simpson.effective)
  otus_div_stats$Evenness <- otus_div_stats$Shannon.Index / log(otus_div_stats$Richness, 2)
  
  otus_div_stats_meta <- merge(otus_div_stats, meta_file, by = 0) %>%
    column_to_rownames("Row.names")
  return(otus_div_stats_meta)
}


Rhea.Beta_Diversity <- function(ps,
                                group_name,
                                label_samples = 0,
                                label_id = "",
                                kmers_limit = 10) {
  # Read and process the metadata file
  meta_file <- meta_to_df(ps)
  meta_file <- data.frame(meta_file[!apply(is.na(meta_file) |
                                             meta_file == "", 1, all), , drop = FALSE])
  meta_file <- data.frame(meta_file[order(row.names(meta_file)), , drop = FALSE])
  meta_file_pos <- which(colnames(meta_file) == group_name)
  all_groups <- as.factor(meta_file[, meta_file_pos])
  
  # Read and process the OTU table
  otu_file <- ps_to_asvtab(ps)
  otu_file <- otu_file[!apply(is.na(otu_file) |
                                otu_file == "", 1, all), ]
  otu_file <- otu_file[, rownames(meta_file)]
  otu_file <- otu_file[, order(names(otu_file))]
  otu_file <- data.frame(t(otu_file), check.names = FALSE)
  
  # Read and process the phylogenetic tree
  tree_file <- ps@phy_tree
  tree_file$tip.label <- gsub("'", "", tree_file$tip.label)
  rooted_tree <- midpoint(tree_file)
  
  # Calculate UniFrac distance matrix
  unifracs <- GUniFrac(otu_file, rooted_tree, alpha = c(0.0, 0.5, 1.0))$unifracs
  unifract_dist <- unifracs[, , "d_0.5"]
  
  # Generate phylogenetic tree plot
  all_dist_matrix <- as.dist(unifract_dist)
  all_fit <- hclust(all_dist_matrix, method = "ward.D2")
  tree <- as.phylo(all_fit)
  
  plot_color <- rainbow(length(levels(all_groups)))[all_groups]
  
  plot(
    tree,
    type = "phylogram",
    use.edge.length = TRUE,
    tip.color = plot_color,
    label.offset = 0.01
  )
  print.phylo(tree)
  axisPhylo()
  tiplabels(pch = 16, col = plot_color)
  
  all_groups_comp <- all_groups[!is.na(all_groups)]
  unifract_dist_comp <- unifract_dist[!is.na(all_groups), !is.na(all_groups)]
  adonis_result <- adonis2(as.dist(unifract_dist_comp) ~ all_groups_comp)
  permdisp_result <- permutest(betadisper(
    as.dist(unifract_dist_comp),
    as.factor(all_groups_comp),
    type = "median"
  ))
  all_groups_comp <- factor(all_groups_comp, levels(all_groups_comp)[unique(all_groups_comp)])
  
  if (nrow(unifract_dist_comp) > 2) {
    s.class(
      cmdscale(unifract_dist_comp, k = 2),
      col = unique(plot_color),
      cpoint = 2,
      fac = all_groups_comp,
      sub = paste0(
        "MDS plot of Microbial Profiles\nPERMDISP     p=",
        permdisp_result[["tab"]][["Pr(>F)"]][1],
        "\n",
        "PERMANOVA  p=",
        adonis_result[1, 5]
      )
    )
    if (label_samples == 1) {
      lab_samples <- row.names(cmdscale(unifract_dist_comp, k = 2))
      if (label_id != "") {
        lab_samples <- replace(lab_samples, !(lab_samples %in% label_id), "")
      }
      text(
        cmdscale(unifract_dist_comp, k = 2),
        labels = lab_samples,
        cex = 0.7,
        adj = c(-.1, -.8)
      )
    }
    
    meta <- metaMDS(unifract_dist_comp, k = 2)
    s.class(
      meta$points,
      col = unique(plot_color),
      cpoint = 2,
      fac = all_groups_comp,
      sub = paste0(
        "metaNMDS plot of Microbial Profiles\nPERMDISP     p=",
        permdisp_result[["tab"]][["Pr(>F)"]][1],
        "\n",
        "PERMANOVA  p=",
        adonis_result[1, 5]
      )
    )
    if (label_samples == 1) {
      lab_samples <- row.names(meta$points)
      if (label_id != "") {
        lab_samples <- replace(lab_samples, !(lab_samples %in% label_id), "")
      }
      text(
        meta$points,
        labels = lab_samples,
        cex = 0.7,
        adj = c(-.1, -.8)
      )
    }
  }
  
  
  # NMDS for pairwise analysis
  unique_groups <- levels(all_groups_comp)
  if (length(unique_groups) > 2) {
    pVal <- numeric()
    permdisppval <- numeric()
    pairedMatrixList <- list()
    pair_1_list <- character()
    pair_2_list <- character()
    
    combn_groups <- combn(unique_groups, 2)
    for (i in 1:ncol(combn_groups)) {
      pair_1 <- combn_groups[1, i]
      pair_2 <- combn_groups[2, i]
      
      pair_1_list[i] <- pair_1
      pair_2_list[i] <- pair_2
      
      inc_groups <- rownames(subset(meta_file, meta_file[, meta_file_pos] == pair_1 |
                                      meta_file[, meta_file_pos] == pair_2))
      paired_dist <- as.data.frame(unifract_dist_comp)
      row_names <- rownames(paired_dist)
      paired_dist <- cbind(row_names, paired_dist)
      paired_dist <- paired_dist[sapply(paired_dist[, 1], function(x)
        all(x %in% inc_groups)), ]
      paired_dist[, 1] <- NULL
      paired_dist <- rbind(row_names, paired_dist)
      paired_dist <- paired_dist[, sapply(paired_dist[1, ], function(x)
        all(x %in% inc_groups))]
      # Remove first row with unnecessary group information
      paired_dist <- paired_dist[-1, ]
      
      
      
      
      
      # Convert generated distance matrix to data type matrix (needed by multivariate analysis)
      paired_matrix <- as.matrix(paired_dist)
      class(paired_matrix) <- "numeric"
      
      # Save paired matrix in list
      pairedMatrixList[[i]] <- paired_matrix
      
      # Applies multivariate analysis to a pair out of the selected groups
      adonis <- adonis2(paired_matrix ~ all_groups_comp[all_groups_comp == pair_1 |
                                                          all_groups_comp == pair_2])
      
      permdisp <- permutest(betadisper(
        as.dist(paired_matrix),
        as.factor(all_groups_comp[all_groups_comp == pair_1 |
                                    all_groups_comp == pair_2]),
        type = "median"
      ),
      pairwise = T)
      
      # List p-values
      pVal[i] <- adonis[1, 5]
      permdisppval[i] <- permdisp$pairwise[2]
      
    }
    
    # Adjust p-values for multiple testing according to Benjamini-Hochberg method
    pVal_BH <- round(p.adjust(pVal, method = "BH", n = length(pVal)), 4)
    permdisppval_BH <- round(p.adjust(
      permdisppval,
      method = "BH",
      n = length(permdisppval)
    ), 4)
    
    
    
    for (i in 1:length(combn(unique_groups, 2)[1, ])) {
      if (nrow(pairedMatrixList[[i]]) > 2) {
        meta <- metaMDS(pairedMatrixList[[i]], k = 2)
        s.class(
          meta$points,
          col = rainbow(length(levels(
            all_groups_comp
          ))),
          cpoint = 2,
          fac = as.factor(all_groups_comp[all_groups_comp == pair_1_list[i] |
                                            all_groups_comp == pair_2_list[i]]),
          sub = paste(
            "NMDS plot of Microbial Profiles\n ",
            pair_1_list[i],
            " - ",
            pair_2_list[i],
            "\n PERMDISP     p=",
            permdisppval[[i]],
            ",",
            "  p.adj=",
            permdisppval_BH[i],
            "\n",
            " PERMANOVA  p=",
            pVal[i],
            ",",
            " p.adj=",
            pVal_BH[i],
            sep = ""
          )
        )
      }
    }
    
    
  }
  return(unifract_dist_comp)
}

# flatten cor matrix to df 
flattenMatrix <- function(mat) {
  ord <- rownames(mat)
  ordmat <- mat[, ord]
  ut <- upper.tri(ordmat, diag = T)
  data.frame(
    row = rownames(ordmat)[row(ordmat)[ut]],
    column = rownames(ordmat)[col(ordmat)[ut]],
    betadiv  =(ordmat)[ut]
  )
}

# filter square matrix
simul_col_row_filt <- function(df, samples_vector, filt="both"){
  if (filt == "rows"){
    row_filt <- filter_rownames(df, samples_vector)
    return(row_filt)
  }
  else if (filt == "columns"){
    col_filt <- df %>% 
      select(all_of(samples_vector))
    return(col_filt)
  }
  else if (filt == "both"){
    row_filt <- filter_rownames(df, samples_vector)
    col_and_row_filt <- row_filt %>% 
      select(all_of(samples_vector))
    return(col_and_row_filt)
  } 
  else {
    print("Error, select from 'columns', 'rows' or 'both'")
  }
  
}
