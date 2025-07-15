library(tidyverse)
library(reshape2)
library(compositions)
library(ggpubr)
library(coin)
library(docstring)
library(roxygen2)
library(EnhancedVolcano)
library(ggsci)
library(lemon)
library(factoextra)
library(ggfortify)
library(phyloseq)
library(microbiome)
library(cowplot)

taxonomy <- function (ps) {
  return(as.data.frame(tax_table(ps)))
}

meta_to_df <- function(ps) {
  return(as(sample_data(ps), "data.frame"))
}

ps_to_asvtab <- function(ps) {
  return(as.data.frame(ps@otu_table))
}


#########################    IMPORTING DATA #########################
split_taxonomy <- function(asvtab) {
  # select taxa column and split taxonomy column into separate columns using ; delimiter
  taxonomy_phylo <- dplyr::select(asvtab, "taxonomy") %>%
    separate("taxonomy",
             c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
             sep = ";")
  return(taxonomy_phylo)
}

format_taxonomy <- function(ps) {
  ps_tax <- taxonomy(ps) %>%
    mutate_all(na_if, "")
  
  ps_tax[] <- t(na.locf(t(ps_tax))) %>%
    as.data.frame()
  # add unknown to genera which match family
  ps_tax <- ps_tax %>%
    mutate(Genus = case_when(Genus == Family ~ paste0("unknown_", Genus),
                             TRUE ~ Genus)) %>%
    mutate(Family = case_when(Family == Order ~ paste0("unknown_", Family),
                              TRUE ~ Family)) %>%
    mutate(Order = case_when(Order == Class ~ paste0("unknown_", Order),
                             TRUE ~ Order)) %>%
    mutate(Class = case_when(Class == Phylum ~ paste0("unknown_", Class),
                             TRUE ~ Class)) %>%
    mutate(Phylum = case_when(Phylum == Kingdom ~ paste0("unknown_", Phylum),
                              TRUE ~ Phylum)) %>%
    mutate(ASV = rownames(ps_tax)) %>%
    mutate(Highest_classified = paste(rownames(ps_tax), Genus, sep = "; "))
  
  rownames(ps_tax) <- ps_tax[, "Highest_classified"]
  taxa_names(ps) <- ps_tax[, "Highest_classified"]
  tax_table(ps) <- tax_table(as.matrix(ps_tax))
  
  return(ps)
}


read_tab_delim <- function(df) {
  # read all tab delimited files using these params
  df_out <-
    read.table(
      df,
      row.names = 1,
      header = 1,
      sep = "\t",
      comment.char = "",
      check.names = F
    )
  return(df_out)
}

load_phylo <- function(asvtab, taxa, mapping, tree = NULL) {
  # convert to phyloseq and return list
  phylo_asv <- otu_table(asvtab, taxa_are_rows = T)
  
  phylo_tax <- tax_table(as.matrix(taxa))
  
  phylo_map <- sample_data(mapping)
  
  if (exists("tree")) {
    phylo_tree <- read_tree(tree)
    return(merge_phyloseq(phylo_asv, phylo_tax, phylo_tree, phylo_map))
  }
  else {
    return(merge_phyloseq(phylo_asv, phylo_tax, phylo_map))
  }
}
transform <- function(ps, transform = "mss", offset=1) {
  if (length(is(ps)) == 1 && class(ps) == "phyloseq") {
    x <- ps_to_asvtab(ps)
  }
  else {
    print("not a phyloseq object, exiting")
    stop()
  }
  
  if (transform %in% c("mss", "relative", "clr")) {
    if (transform == "mss") {
      ps_t <- t(min(colSums(x)) * t(x) / colSums(x))
      
    } else if (transform == "relative") {
      ps_t <- t(100 * t(x) / colSums(x))
      
    } else if (transform == "clr") {
      ps_t <- t(logratio.transfo(t(x + offset), logratio = "CLR"))
      # fix mixOmics class issue
      class(ps_t) <- "matrix"
    }
    otu_table(ps)@.Data <- ps_t
    
    return(ps)
    
  } else {
    print("Not a valid transform, exiting")
    stop()
  }
  
}

import_as_pseq <- function(asvtab, mapping, tree = NULL) {
  # read files
  asvtab_taxa <- read_tab_delim(asvtab)
  metadata <- read_tab_delim(mapping)
  # convert otu into matrix and drop taxonomy col
  asv_matrix <- as.matrix(subset(asvtab_taxa, select = -taxonomy))
  # split taxonomy
  taxonomy_phylo <- split_taxonomy(asvtab_taxa)
  # make sure taxa table row names match otu
  rownames(taxonomy_phylo) <- rownames(asvtab_taxa)
  # use load phylo function to convert all files to phyloseq objs
  if (exists("tree")) {
    # if user wants to input a phylogenetic tree this code is called
    out <-
      load_phylo(asv_matrix, taxonomy_phylo, metadata, tree = tree)
  }
  else {
    # without a tree
    out <- load_phylo(asv_matrix, taxonomy_phylo, metadata)
  }
  return(out)
}

names_to_vec <- function(mapping, mapping_col, filter_col, filter_by) {
  vec <- filter(mapping, {{ filter_col }} == filter_by) %>% 
    pull({{ mapping_col }})  
  return(vec)
}

round_any <- function(x, accuracy, f=ceiling) {
  f(x/ accuracy) * accuracy
}

filter_cols_from_rows <- function(df_cols, df_rows, df_rows_col, filter_var_col, filter_var, keep_col) {
  vec <- names_to_vec(df_rows, mapping_col = {{ df_rows_col }}, filter_col = {{ filter_var_col }}, filter_var)
  df_filtered <- select(df_cols, vec, one_of(keep_col))
  return(df_filtered)
}

melt_gather <- function(df, col_excl) {
  melt_df <- gather(df, variable, value, -col_excl)
  return(melt_df)
}

calc_gram_stain<- function(df, gather_col, value_col) {
  df_out <- df %>% 
    select({{ gather_col }}, {{ value_col }}) %>% 
    group_by({{ value_col }}) %>% 
    summarise_all(mean)
  return(df_out)
}

normalize <- function(df, method) {
  if (method == "relative") {
    rel_tab <- t(100 * t(df) / colSums(df))
    return(as.data.frame(rel_tab))
  }
  else {
    min_sum <- min(colSums(df))
    norm_tab <- t(min_sum * t(df) / colSums(df))
    return(as.data.frame(norm_tab))
  }
}

filter_prev_abund <- function(asvtab, abund, prev) {
  asvtab_filt <- asvtab[which(rowSums(asvtab > abund) 
                          >= prev * ncol(asvtab)), ]
  return(asvtab_filt)
} 

filter_rownames <- function(df, filt_vector) {
  # wrapper script around filter to handle rownames
  df_filt <- df %>% 
    rownames_to_column(var="id") %>%
    filter(id %in% filt_vector) %>%
    column_to_rownames(var="id")
  return(df_filt)
}


picrust_pca <- function(picrust_matrix, column_1, column_2, meta_variable) {
  # picrust_matrix - picrust data [matrix]
  # column 1 - 1st numeric column [string]
  # last numeric column [string]
  # metadata variable to colour by [string] 
  #' @description This function generates a PCA plot from picrust data
  #' @param picrust_matrix: matrix matrix of picrust data (samples are rows), 
  #' also containing metadata information   
  #' @param column_1 string 1st numeric column
  #' @param column_2 string last numric column
  #' @param meta_variable string metadata variable to color PCA by
  #' @usage picrust_pca(picrust_matrix, column_1, column_2, meta_variable)
  #' @return PCA plot colored by metavariables
  pca <- prcomp(picrust_matrix[c(col1:col2)], scale. = T)
  scree_plot <- fviz_eig(pca)
  pca_plot <- autoplot(pca, data = ko_tumor_filt_meta,
                       colour = meta_variable) + 
    geom_point(aes(color= {{ meta_variable }} ), size=3, alpha=0.75) + 
    scale_color_manual(values = c("#05A8AA", "#E26D5A")) + 
    theme_cowplot()
  pca_plot$layers[[1]] <- NULL
  return(p)
}

fc <- function(seqtab, meta, read=F, log.n0=1e-05, mult.corr="fdr",
               var, case, ctrl, feat_type="taxa") {
  #' @description This function calculates differentially abundant ASVs 
  #' based on fdr adjusted wilcoxon tests and then calculates generalized 
  #' fold change, as introduced in Wirbel et al. 2019 Nat Med. Script adapted from: 
  #' https://github.com/zellerlab/crc_meta/blob/master/src/marker_analysis.R 
  #' @param seqtab: character. File path of ASV/OTU table or object if read=F 
  #' @param meta: character. File path of metadata or object if read=F
  #' @param read: boolean. Whether to read files or not. Defaults to false
  #' @param log.n0 integer. Pseudocount for log(X + 1) transformation. 
  #' Default 1e-05
  #' @param mult.corr character. Multiple testing correction. Passed to p.adjust. 
  #' Default fdr/BH
  #' @param var character. Column variable to test on 
  #' @param case character. Case group 
  #' @param ctrl character. Control/Reference group 
  #' @param feat_type character taxa for ASV/OTU/species abundances or function for 
  #' gene/predicted gene abundances 
  #' @usage gfc(seqtab, meta, var, case, ctrl)
  #' @return List containing two dataframes: one containing pvalues and one
  #' containing generalized fold change
  if (read == T){
    seqtab <- as.matrix(read.table(seqtab, sep='\t', header=TRUE, 
                                   stringsAsFactors = FALSE, 
                                   check.names = FALSE, quote='', row.names = 1))
    meta <- read_tsv(meta)
  }
  
  if (feat_type == "functional"){
    log.n0 = 1e-08
  }
  # infer name of sample column 
  sampleid_col <- colnames(meta)[1]
  
  stopifnot(all(meta[sampleid_col] %>% pull() 
                %in% colnames(seqtab)))
  
  # create pval matrix with row for each feature
  p.val <- matrix(NA, nrow=nrow(seqtab), ncol=1, 
                  dimnames=list(row.names(seqtab), c("pval")))
  # duplicate pval matrix to store gfc 
  fc <- p.val
  colnames(fc) <- c("FC") 
  
  
  cat("Calculating pvalue and FC for each feature...\n")
  pb <- txtProgressBar(max=nrow(seqtab), style=3)
  
  # calculate wilcoxon test and effect size for each feature and study
  for (f in row.names(seqtab)) {
    
    x <- unlist(seqtab[f, meta %>% 
                         filter(.data[[var]]==case) 
                       %>% pull(sampleid_col)])
    
    y <- unlist(seqtab[f, meta %>% 
                         filter(.data[[var]]==ctrl) 
                       %>% pull(sampleid_col)])
    
    # Wilcoxon
    p.val[f, ] <- wilcox.test(x, y, exact=FALSE)$p.value
    
    # FC - doubble check how to calculate FC
    q.p <- log2(x+log.n0)
    q.n <- log2(y+log.n0)
    fc[f,] <- median(q.p - q.n)
    # progressbar
    setTxtProgressBar(pb, (pb$getVal()+1))
  }
  cat('\n')
  
  # multiple hypothesis correction
  p.adj <- data.frame(apply(p.val, MARGIN=2, FUN=p.adjust, method=mult.corr),
                      check.names = FALSE)
  ### Return datafram containing pval and log2FC
  
  df_res <- data.frame(p.adj,fc)
  colnames(df_res) <- c("p.adj", "log2FC")
  
  comparison <- sort(c(case, ctrl))
  
  return_list <- list(res = df_res, groups = comparison)
}



plot_volcano <- function(results_df, pthresh=0.05, FCthresh=0.58, ...) {
  #' @description This function takes the ouput of the GFC function and generates 
  #' a volcano plot
  #' @param results_df: object list ouput of GFC containing dataframe with adjusted
  #' p-values and log2 fold change and vector of the groups being compared  
  #' @details additional arguments are passed to EnhancedVolcano
  #' @usage plot_volcano(results_df, ...)
  #' @return Volcano plot generated by EnhancedVolcano
  res <- results_df$res
  groups <- results_df$groups
  EnhancedVolcano(res, x="log2FC", y="p.adj", 
                  lab=rownames(res),
                  labSize = 4,
                  pCutoff = pthresh, 
                  FCcutoff = FCthresh, 
                  title = paste(groups[1], 
                                "versus", 
                                groups[2], 
                                sep=" "), 
                  subtitle = NULL,...)
}

filter_feat <- function(feat_profile,
                        by = "both",
                        features = NULL,
                        samples = NULL,
                        feat_are_rows = TRUE, 
                        partial_feat_matching=FALSE) {
  if (partial_feat_matching){
    features <- partial_matches(feat_profile, strings=features)
  }
  
  if (by == "features") {
    if (feat_are_rows == T) {
      feat_filt <- filter_rownames(feat_profile, features)
    } else {
      feat_filt <- select(feat_profile, all_of(features))
    }
    
  } else if (by == "samples") {
    if (feat_are_rows == T) {
      feat_filt <- select(feat_profile, all_of(samples))
    } else {
      feat_filt <- filter_rownames(feat_profile, samples)
    }
  } else if (by == "both") {
    feat_filt <-
      filter_feat(feat_profile,
                  by = "features",
                  features = features,
                  feat_are_rows)
    feat_filt <-
      filter_feat(feat_filt,
                  by = "samples",
                  samples = samples,
                  feat_are_rows)
  }
  return(feat_filt)
}

partial_matches <- function(feat_profile, strings) {
  feat_filt <- feat_profile[grep(paste(strings, collapse = "|"),rownames(feat_profile)), ]
  return(rownames(feat_filt))
}

plot_feat <-
  function(feat_string,
           feat_table,
           x,
           y,
           fill,
           comps,
           cols,
           ord=NULL,
           ylab = NULL,
           xlab = NULL) {
    feat <- feat_table %>%
      rownames_to_column("ko") %>%
      select(.data[[ x ]], contains(feat_string))
    p <- ggboxplot(
      feat,
      x = x,
      y = y,
      fill = fill,
      order = ord,
      bxp.errorbar = TRUE,
      bxp.errorbar.width = 0.05,
      ylab = ylab
    ) +
      geom_point(aes(color = .data[[fill]]), position = position_jitterdodge()) +
      stat_compare_means(comparisons = comps, label = "p.format") +
      scale_color_manual(values = cols) +
      scale_fill_manual(values = cols)
    return(p)
  }

get_gene_presence <- function(func_pa, asvtab, func) {
  
  asv_labels <- asvtab %>% 
    rownames()
  
  taxonomy <- asvtab %>% 
    select(taxonomy)
  
  rownames(taxonomy) <- asv_labels
  
  feat_tax <- merge(func_pa, taxonomy, by=0) %>% 
    column_to_rownames("Row.names")
  taxa_func_map <- feat_tax %>% 
    select(all_of(c(func, "taxonomy"))) %>% 
    merge(asvtab %>% select(-taxonomy), by=0) %>% 
    column_to_rownames("Row.names")
  
  return(taxa_func_map)
}

gene_stats <- function(gene_db, gene) {
  encoding_sum <- gene_db %>% filter(.data[[ gene ]] == 1) %>% 
    select(where(is.numeric)) %>% 
    summarise_all(.funs = sum)
  non_encoding_sum <- gene_db %>% filter(.data[[ gene ]] == 0) %>% 
    select(where(is.numeric)) %>% 
    summarise_all(.funs = sum)
  e_non_e_ratio <- encoding_sum / non_encoding_sum 
  
  data_out <- bind_rows(encoding_sum, non_encoding_sum, e_non_e_ratio) %>% 
    add_column(variable =c("encoding", "non-encoding", "ratio")) %>% 
    select(-c(.data[[ gene ]] )) %>% 
    pivot_longer(-variable)
  return(data_out)
}

# create statistical formula
xyform <- function (y_var, x_vars) {
  # y_var: a length-one character vector
  # x_vars: a character vector of object names
  as.formula(sprintf("%s ~ %s", y_var, paste(x_vars, collapse = " + ")))
}


pairwise_wilcox_asv <- function(asv_tab, factor, case, control, 
                                sample_col, meta, p_adjust="BY"){
  
  # initialise pval matrix 
  p.val <- matrix(NA, nrow=nrow(asv_tab), ncol=1, 
                  dimnames=list(row.names(asv_tab), "p.val"))
  
  featmat <- as.matrix(asv_tab)
  
  for (row in rownames(featmat)){
    x <- featmat[row, meta %>% 
                   filter(.data[[ factor ]] ==control) %>% pull(.data[[ sample_col ]])]
    y <- featmat[row, meta %>% 
                   filter(.data[[ factor ]] ==case) %>% pull(.data[[ sample_col ]])]
    
    
    p.val[row, ] <- wilcox.test(x, y, exact=FALSE)$p.value
  }
  
  p.val.adj <- p.val %>% 
    as.data.frame() %>% 
    adjust_pvalue("p.val", method = p_adjust)
  
  return(p.val.adj)
}


add_tax <- function(func_strat, taxonomy) {
  taxa_list <- as.character(rownames(taxonomy))
  names(taxa_list) <- taxonomy$taxonomy
  func_strat$taxonomy <-
    sapply(func_strat$taxon, function(x)
      names(taxa_list)[match(x, taxa_list)])
  return(func_strat)
}

calc_other_site <- function(func_means, taxonomic_level) {
  others <- func_means %>% 
    group_by(site) %>% 
    mutate(other_abund = 100 - sum(Mean_Rel_Abund_Tax)) %>% 
    mutate(other_abund = case_when(other_abund < 0 ~ 0,
                                   TRUE ~ other_abund)) %>% 
    distinct(site, .keep_all = T) %>% 
    pull(other_abund)
  
  n_sites <- length(unique(func_means$site))
  
  sites <- seq(1, n_sites, 1)
  char <- rep("Others", n_sites)
  
  tmp <- data.frame(site = sites,
                    tmp_level = char,
                    Mean_Rel_Abund_Tax = others) 
  names(tmp)[names(tmp) == "tmp_level"] <- taxonomic_level
  
  func_means_out <- func_means %>%
    bind_rows(tmp) %>%
    arrange(site)
  return(func_means_out)
}

calc_other_var <- function(func_means, taxonomic_level, var) {
  others <- func_means %>% 
    group_by(.data[[ var ]]) %>% 
    mutate(other_abund = 100 - sum(Mean_Rel_Abund_Tax)) %>% 
    mutate(other_abund = case_when(other_abund < 0 ~ 0,
                                   TRUE ~ other_abund)) %>% 
    distinct(.data[[ var ]], .keep_all = T) %>% 
    pull(other_abund)
  
  vars <- unique(func_means[[var]])
  
  char <- rep("Others", length(vars))
  
  tmp <- data.frame(tmp_var = vars,
                    tmp_level = char,
                    Mean_Rel_Abund_Tax = others) 
  names(tmp)[names(tmp) == "tmp_level"] <- taxonomic_level
  names(tmp)[names(tmp) == "tmp_var"] <- var
  
  func_means_out <- func_means %>%
    bind_rows(tmp) %>%
    arrange(.data[[ var ]])
  return(func_means_out)
}

taxonomic_contributions_site <- function(func_strat, metadata, samples, 
                                    func, taxonomy, taxonomic_level, n_taxa) {
  
  # filter for samples and gene of interest and add site column for grouping
  func_filt <- func_strat %>% 
    filter(sample %in% samples) %>% 
    filter(`function`==func) %>% 
    mutate(site = as.integer(map_chr(sample, function(s) rev(strsplit(s, "-")[[1]])[1])))
  
  func_contrib_means <- func_filt %>% 
    group_by(site, taxon) %>% 
    summarise(Site_Mean = mean(taxon_function_abun)) %>% 
    group_by(site) %>% 
    mutate(Mean_Rel_Abund = Site_Mean / sum(Site_Mean) * 100)
  
  func_contrib_means_tax <- add_tax(func_contrib_means, taxonomy = taxonomy) %>% 
    group_by(site,taxonomy) %>%
    summarise(Mean_Rel_Abund_Tax = sum(Mean_Rel_Abund)) %>%
    separate(taxonomy,
             into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
             sep = ";") %>%
    mutate(
      across(Kingdom:Genus, ~na_if(., ""))
    )
  
  func_contrib_means_tax[] <- t(na.locf(t(func_contrib_means_tax)))
  func_contrib_means_tax <- func_contrib_means_tax %>%
    # adjust this code here
    mutate(Genus = case_when(Genus == Family ~ paste0("unknown_", Genus),
                             TRUE ~ Genus)) %>%
    mutate(Family = case_when(Family == Order ~ paste0("unknown_", Family),
                              TRUE ~ Family)) %>%
    mutate(Order = case_when(Order == Class ~ paste0("unknown_", Order),
                             TRUE ~ Order)) %>%
    mutate(Class = case_when(Class == Phylum ~ paste0("unknown_", Class),
                             TRUE ~ Class)) %>%
    mutate(Phylum = case_when(Phylum == Kingdom ~ paste0("unknown_", Phylum),
                              TRUE ~ Phylum)) %>%
    mutate(site = as.integer(site)) %>%
    mutate(Mean_Rel_Abund_Tax = as.double(Mean_Rel_Abund_Tax)) %>% 
    # keep taxonomy here
    select(site,Mean_Rel_Abund_Tax, .data[[ taxonomic_level ]]) %>%
    group_by(site) %>%
    slice_max(order_by = Mean_Rel_Abund_Tax, n=n_taxa)
  
  func_contrib_means_out <- calc_other_site(func_contrib_means_tax, 
                                       taxonomic_level = taxonomic_level)
  return(func_contrib_means_out)
}

taxonomic_contributions_var <- function(func_strat_meta, samples, var,
                                         func, taxonomy, taxonomic_level, n_taxa) {
  # filter for samples and gene of interest and add site column for grouping
  func_filt <- func_strat_meta %>% 
    filter(sample %in% samples) %>% 
    filter(`function`==func) 

  func_contrib_means <- func_filt %>% 
    group_by(.data[[ var ]], taxon) %>% 
    summarise(Mean = mean(taxon_function_abun)) %>% 
    group_by(.data[[ var ]]) %>% 
    mutate(Mean_Rel_Abund = Mean / sum(Mean) * 100)
  
  func_contrib_means_tax <- add_tax(func_contrib_means, taxonomy = taxonomy) %>% 
    group_by(.data[[ var ]],taxonomy) %>%
    summarise(Mean_Rel_Abund_Tax = sum(Mean_Rel_Abund)) %>%
    separate(taxonomy,
             into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
             sep = ";") %>%
    mutate(
      across(Kingdom:Genus, ~na_if(., ""))
    )
  
  func_contrib_means_tax[] <- t(na.locf(t(func_contrib_means_tax)))
  func_contrib_means_tax <- func_contrib_means_tax %>%
    # adjust this code here
    mutate(Genus = case_when(Genus == Family ~ paste0("unknown_", Genus),
                             TRUE ~ Genus)) %>%
    mutate(Family = case_when(Family == Order ~ paste0("unknown_", Family),
                              TRUE ~ Family)) %>%
    mutate(Order = case_when(Order == Class ~ paste0("unknown_", Order),
                             TRUE ~ Order)) %>%
    mutate(Class = case_when(Class == Phylum ~ paste0("unknown_", Class),
                             TRUE ~ Class)) %>%
    mutate(Phylum = case_when(Phylum == Kingdom ~ paste0("unknown_", Phylum),
                              TRUE ~ Phylum)) %>%
    mutate(Mean_Rel_Abund_Tax = as.double(Mean_Rel_Abund_Tax)) %>%
    # keep taxonomy here
    select(.data[[ var ]],Mean_Rel_Abund_Tax, .data[[ taxonomic_level ]]) %>%
    group_by(.data[[ var ]]) %>%
    slice_max(order_by = Mean_Rel_Abund_Tax, n=n_taxa)

  func_contrib_means_out <- calc_other_var(func_contrib_means_tax,
                                       taxonomic_level = taxonomic_level,
                                       var = var)
  return(func_contrib_means_out)
}

plot_boxplot <- function(df,
                         variable_col,
                         value_col,
                         comparisons_list,
                         fill_var = variable_col,
                         xlab = variable_col,
                         ylab = value_col,
                         p_title = NULL,
                         multiple_groups = FALSE,
                         cols = NULL,
                         group.order = NULL,
                         paired = FALSE,
                         ...) {
  # extend color palette with transparent value - required due to way we are
  # layering plot
  if (is.null(cols)) {
    cols <- pal_npg()(length(unique(df[, variable_col])))
  }
  cols <- c(cols, "transparent")
  
  if (!is.null(group.order)) {
    df[, variable_col] <-
      factor(df[, variable_col], levels = group.order)
  }
  
  formula <- xyform(value_col, variable_col)
  
  if (multiple_groups == TRUE) {
    if (paired == TRUE) {
      stat_variance <- df %>%
        friedman_test(formula)
      stat_test <- df %>%
        pairwise_wilcox_test(
          formula,
          comparisons = comparisons_list,
          p.adjust.method = "BH",
          paired = TRUE
        ) %>%
        add_significance() %>%
        add_xy_position(x = variable_col) 
    }
    else {
      stat_variance <- df %>%
        rstatix::kruskal_test(formula)
      stat_test <- df %>%
        pairwise_wilcox_test(formula,
                             comparisons = comparisons_list,
                             p.adjust.method = "BH") %>%
        add_significance() %>%
        add_xy_position(x = variable_col) 
    }
  }
  else if (multiple_groups == FALSE) {
    if (paired == TRUE) {
      stat_test <- df %>%
        rstatix::wilcox_test(formula, paired = TRUE) %>%
        add_significance() %>%
        add_xy_position(x = variable_col) 
    }
    else {
      stat_test <- df %>%
        rstatix::wilcox_test(formula) %>%
        add_significance() %>%
        add_xy_position(x = variable_col)
    }
    
  }
  
  # aes string accepts strings as column names, this code plots boxplot and adds error bars
  plot <- ggplot(
    df,
    aes_string(
      x = variable_col,
      y = value_col,
      fill = variable_col,
      color = variable_col
    )
  ) +
    geom_boxplot(
      color = "black",
      alpha = 0.8,
      outlier.shape = 5,
      outlier.size = 1
    ) +
    geom_point(size = 1.5, position = position_jitterdodge()) +
    labs(x = xlab, y = ylab) +
    stat_boxplot(color = "black",
                 geom = "errorbar",
                 width = 0.2)
  # creates new 'finalised plot' and adds statistical significance, labels and adjusts theme and title
  final_plot <- plot +
    theme_classic2() +
    ggtitle(p_title) +
    theme(
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      axis.line = element_line(colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position = "None"
    ) +
    scale_fill_manual(values = cols) +
    scale_color_manual(values = cols) +
    rotate_x_text(angle = 45)
  
  if (dim(stat_test)[1] == 0) {
    plot_out <- final_plot
  }
  else {
    
    
    if (multiple_groups == T) {
      plot_out <- final_plot +
        stat_pvalue_manual(
          stat_test,
          label = "p.adj",
          inherit.aes = FALSE,
          ...
        ) #+
        #coord_capped_cart(left='top', expand = F)
    }
    else {
      plot_out <- final_plot +
        stat_pvalue_manual(
          stat_test,
          label = "p",
          inherit.aes = FALSE,
          ...
        ) #+
        #coord_capped_cart(left='top', expand = F)
      
    }
  }
  
  return(plot_out)
}

# plot scatter plot with correlation if desired
plot_scatter <- function(df,
                         x,
                         y,
                         point_color,
                         line_color,
                         fill_color,
                         xlab,
                         ylab,
                         corr.method = NULL,
                         ...) {
  p <-
    ggplot(data = df, mapping = aes(x = .data[[x]], y = .data[[y]])) +
    geom_point(aes(color = point_color), size = 2.5) +
    geom_smooth(method = "lm",
                color = line_color,
                fill = fill_color) +
    theme_bw() +
    theme(
      legend.position = "None",
      axis.title.x = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.title.y = element_text(size = 14),
      axis.text.y = element_text(size = 12)
    ) +
    xlab(xlab) +
    ylab(ylab)
  
  if (!is.null(corr.method)) {
    p <- p + stat_cor(method = corr.method, ...)
    return(p)
  }
  else {
    return(p)
  }
}


t_df <- function(x) {
  return(as.data.frame(t(x)))
}

run_lefser  <- function(ps, group_col, ...) {
  counts <- unclass(otu_table(ps))
  colData <- as(sample_data(ps), "data.frame")
  ## create a SummarizedExperiment object
  lefse_in <- SummarizedExperiment(assays = list(counts = counts), colData = colData)
  lefse_in_norm <- relativeAb(lefse_in)
  res <- lefser(lefse_in_norm, groupCol = group_col, ...)
  
  
  return(res)
}

filter_ps <- function(ps, prev=0.33, abund=0.25) {
  
  ps_filt <- ps
  
  asvtab <- ps_to_asvtab(ps)
  
  asvtab_filt <- filter_prev_abund(asvtab, abund, prev)
  
  otu_table(ps_filt) <- otu_table(asvtab_filt, taxa_are_rows = T)
  
  return(ps_filt)
}

# should merge selected meta column(s) with ASV tab at first position and transpose
format_lefse <- function(ps, variables) {
  # get metadata
  meta <- meta_to_df(ps)
  # get ASV table
  asvtab <- ps_to_asvtab(ps) %>% 
    t_df()
  # merge metadata with ASV table
  lefse_df <- merge(meta[variables], asvtab, by = 0) %>% 
    column_to_rownames("Row.names")
  # transpose
  lefse_df <- t_df(lefse_df)
  return(lefse_df)
}

### Metabolomics analysis