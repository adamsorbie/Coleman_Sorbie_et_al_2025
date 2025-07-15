library(phyloseq)
library(microbiome)
library(microbiomeutilities)
library(ComplexHeatmap)
library(RColorBrewer)

preprocess_ps <- function(ps, list_taxa, transformation) {
  # this function takes a phyloseq object, prunes taxa, adds a variance 
  # stabilizing transformation and outputs the meta and otu tables as a df
  ps <- prune_taxa(list_taxa, ps)
  ps_f <- format_to_besthit(ps) %>% 
    microbiome::transform(transformation)
  
  return_list <- list("OTU" = abundances(ps_f), "Meta" = meta(ps_f), "Phylo_obj" = ps_f)
  return(return_list)
}

generate_row_ann <- function(ps, level) {
  
  # levels dict 
  levels <- list("Kingdom" = 1, "Phylum" = 2, 
                 "Class" = 3, "Order" = 4, 
                 "Family" = 5, "Genus" = 6)
  # subset tax table by desired taxonomic level
  row_ann <- data.frame(tax_table(ps)[, levels[[level]] ])
  # set rownames as ASV ids 
  rownames(row_ann) <- rownames(abundances(ps))
  
  # get unique taxonomy from row ann df 
  uniq_taxa <- unique(row_ann[,1])
  # get number of unique taxa
  n_taxa <- length(uniq_taxa)
  
  # set color palette based on number of unique taxa
  if (n_taxa <= 2){
    taxa_palette <- brewer.pal(3, "Dark2")[1:n_taxa]
  }
  else if (n_taxa > 2 & n_taxa < 8){
    taxa_palette <- brewer.pal(n_taxa, "Dark2")
  }
  else {
    taxa_palette <- colorRampPalette(brewer.pal(8, "Dark2"))(n_taxa)
  }
  
  # set names 
  names(taxa_palette) <- uniq_taxa
  # generate row annotation
  row_ha <- rowAnnotation(df = row_ann, col=list(Genus = taxa_palette),
                          show_annotation_name=F)
  return(row_ha)
}

generate_col_ann <- function(meta, variable1, name_ann1, multiple_an=F, 
                             variable2=NULL, name_ann2=NULL){
  if (multiple_an == FALSE){
    column_ha <- HeatmapAnnotation(df = meta[c(variable1)], 
                                   col = list(name_ann1 = variable1),
                                   border = T, show_annotation_name = F) 
  } else {
    column_ha <- HeatmapAnnotation(df = meta[c(variable1, variable2)], 
                                   col = list(name_ann1 = variable1,
                                              name_ann2 = variable2),
                                   border = T, show_annotation_name = F)
  }
  
  return(column_ha)
}

plot_heatmap <- function(ps, list_taxa=NULL, 
                         variable, variable.colors, heatmap.colors, 
                         heatmap.type="complexheatmap", split_on=NULL, 
                         transform="Z", ord=F, col_ord=NULL, row_ord=NULL, ...) {
  
  if (!is.null(list_taxa)){
    ps_in <- prune_taxa(list_taxa, ps) 
  } else {
    ps_in <- ps 
  }
 
  ps_z <- transform(ps_in, transform)
  
  otu.mat <- abundances(ps_z)
  
  if (ord == T){
    if (!is.null(col_ord)){
      otu.mat <- otu.mat[,col_ord]
    }
    if (!is.null(row_ord)){
      otu.mat <- otu.mat[row_ord,]
    }
  }
  meta <- meta(ps_z)
  
  subset_meta <- subset(meta, select = c(variable))
  subset_meta_ord <- subset_meta[order(match(rownames(subset_meta), col_ord)), , drop = FALSE]
  
  # Get genus name
  row_ann <- data.frame(tax_table(ps_z)[, 6])
  row_ann_ord <- row_ann[order(match(rownames(row_ann), row_ord)), , drop = FALSE]
  # ASV names as rownames
  rownames(row_ann_ord) <- rownames(otu.mat)
  
  if (heatmap.type == "pheatmap"){
    p <- pheatmap::pheatmap(otu.mat, annotation_col = subset_meta, 
                            annotation_row = row_ann_ord, color=heatmap.colors, 
                            cluster_cols=F, cluster_rpows=F, border_color = "white", 
                            show_colnames = F, cellwidth=12, cellheight=12,
                            show_rownames = F, ...)
    return(p)
  }
  else if (heatmap.type == "complexheatmap"){
    
    
    
    unique_taxa <- unique(row_ann[,1]) 
    n_taxa <- length(unique_taxa)
    #paired or dark color brewer palette 
    # set color palette based on number of unique taxa
    if (n_taxa <= 2){
      taxa_palette <- pals::stepped3(3)[1:n_taxa]
    }
    else {
      taxa_palette <- sample(unique(c(pals::stepped3(), pals::stepped2())), size = n_taxa)
    }
    
    # set names 
    names(taxa_palette) <- unique_taxa
    
    
    
    row_ha <- rowAnnotation(df = row_ann_ord, col=list(Genus = taxa_palette),
                            show_annotation_name=F)
    column_ha <- HeatmapAnnotation(df = subset_meta_ord, 
                                   col = list(variable = variable.colors),
                                   border = T, show_annotation_name = F)

    p <- Heatmap(otu.mat, right_annotation = row_ha, top_annotation = column_ha,
                 col=heatmap.colors, cluster_rows = T,cluster_columns = F,
                 show_row_names = F,
                 show_column_names = F, name = "Z-score abundance",
                 row_title = NULL, row_split = split_on, border = "black",
                 column_gap = unit(1.5, "mm"), ...)

    return(draw(p, heatmap_legend_side="right", annotation_legend_side="right"))
    return(otu.mat)
   
  }
}

plot_multi_heatmap <- function(ps1, ps2, list_taxa1, list_taxa2,
                               heatmap.colors, column_ha1, column_ha2=FALSE, 
                               ord1, ord2, ...) {
  
  otu_meta_ps1 <- preprocess_ps(ps1,list_taxa1, transformation = "Z")
  
  otu_meta_ps2 <- preprocess_ps(ps2,list_taxa2, transformation = "Z")
  
  row_ha1 <- generate_row_ann(otu_meta_ps1$Phylo_obj, "Genus")
  row_ha2 <- generate_row_ann(otu_meta_ps2$Phylo_obj, "Genus")
  
  p1 <- Heatmap(otu_meta_ps1$OTU, right_annotation = row_ha1, top_annotation = column_ha1, 
                col=heatmap.colors, cluster_columns = F, show_row_names = F, 
                show_column_names = F, name = "Z-score abundance",
                row_title = NULL, border = "black", 
                column_gap = unit(1.5, "mm"),column_order = ord1, ...)
  
  p2 <- Heatmap(otu_meta_ps2$OTU, right_annotation = row_ha2, top_annotation = column_ha2, 
                col=heatmap.colors, cluster_columns = F, show_row_names = F, 
                show_column_names = F, name = "Z-score abundance",
                row_title = NULL, border = "black", 
                column_gap = unit(1.5, "mm"),column_order = ord2, ...)
  
  return_list <- list("Heatmap1" = p1, "Heatmap2" = p2)
  return(return_list)
}