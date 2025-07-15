#' Functions to handle imputation, and scaling of metabolome data, while handling
#' feature metadata

metabo_list <- function(df, meta, sampleid_col){
  
  sample_names <- meta %>% pull(.data[[sampleid_col]])
  
  features <- df %>% 
    dplyr::select(all_of(sample_names))
  
  feat_meta <- df %>% 
    dplyr::select(!all_of(sample_names))
  
  return_list <- list("Features" = features, "Feature_Metadata" = feat_meta)
  return(return_list)
}



half_min_imp <- function(metabo_list) {
  
  # transpose
  df_t <- as.data.frame(t(metabo_list$Features))
  
  # colwise min, fill NA with 0.5 * min 
  half_min <- function(x) replace(x, is.na(x), min(x, na.rm = T) * 0.5)
  df_imp <- replace(df_t, TRUE, lapply(df_t, half_min))
  
  # re-tranpose and add back to metabolist
  metabo_list$Features <- as.data.frame(t(df_imp))
  
  return(metabo_list)
}
#######################################

# adjust so you can specify number of samples as well as % 
filter_prevalence <- function(metabo_list, filt_type="prev_trh", prev_trh=0.5, min_samples=0) {
  
  if (filt_type == "prev_trh"){ 
    feat_filt <- metabo_list$Features[which(rowMeans(is.na(metabo_list$Features)) <= prev_trh), ]
  }
  else if (filt_type == "min_samples"){
    
    if (min_samples == 0){
      stop("No value specified for minimum samples")
    } 
    else {
      feat_filt <- metabo_list$Features[which(rowSums(is.na(metabo_list$Features)) <= min_samples), ]
    }
  }
  
  feat_meta_filt <- metabo_list$Feature_Metadata[rownames(feat_filt), ]
  
  return_list <- list("Features" = feat_filt, 
                      "Feature_Metadata" = feat_meta_filt)
  return(return_list)
}