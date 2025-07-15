filter_prev_abund <- function(df, prev, abun) {
  df_filt <- df[rowSums( df >= abun) >= ncol(df) * prev, ]
  return(as.data.frame(df_filt))
}

split_df <- function(df, map, col, var){
  list_samples <- filter(map, map[[col]] == var) %>% 
    pull(var=1)
  filtered_df <- df[colnames(df) %in% list_samples]
  return(filtered_df)
}

normalise <- function(df) {
  normalised <- t(100 * t(df) / colSums(df))
  return(normalised)
}


rem_zero <- function(otu_table) {
  zero_rem <- normalised[apply(normalised, 1, function(X) any(x != 0 | is.na(x))), ]
  return(zero_rem)
}
# define function to filter tables by prevalence and abundance 
filter_otu <- function(df, prev, abund) {
  # this is not the best code
  df_filt <- normalise(df) %>% filter_prev_abund(prev, abund) %>% rem_zero()
  return(df_filt)
}

get_samps <- function(metadata, variable, group) {
  # filter by criteria and return rownames (sample names) as vector
  samps <- metadata %>% filter(.data[[variable]] == group) %>% 
    rownames()
  return(samps)
}

filter_rownames <- function(df, filt_vector) {
  # wrapper script around filter to handle rownames
  df_filt <- df %>% 
    rownames_to_column(var="id") %>%
    filter(id %in% filt_vector) %>%
    column_to_rownames(var="id")
  return(df_filt)
}

# filter based on criteria 
subset_group <- function(df, metadata, variable, group) {
  # get list of samples matching criteria as vector
  samps <- get_samps(metadata, variable, group)
  
  # check if df and convert if not
  if (is.data.frame(df)) {
    # use filter rownames function to filter by rowname
    df_filt <- filter_rownames(df, samps)
  }
  else {
    df <- as.data.frame(df)
    # use filter rownames function to filter by rowname
    df_filt <- filter_rownames(df, samps)
  }
  return(df_filt)
}

colfrom_named_vec <- function(vector, df, col, newcol) {
  df[, newcol] <- lapply(df[, col], function(x) names(vector)[match(x, vector)])
  return(df)
}