# SIAMCAT functions 

library(SIAMCAT)

sc.object <- function(ps, label, case, norm = TRUE) {
  if (norm == FALSE) {
    ps <- transform(ps, transform = "relative")
  }
  label <- create.label(label = label, case = case, sample_data(ps))
  sc.obj <- siamcat(phyloseq = ps, label = label)
  return(sc.obj)
}

sc.preprocess  <-
  function(sc.object,
           abund,
           prev,
           normalize = TRUE,
           norm.method = "log.clr") {
    sc.obj.filt <- sc.object %>%
      filter.features(cutoff = abund, filter.method = 'abundance') %>%
      filter.features(cutoff = prev,
                      filter.method = 'prevalence',
                      feature.type = 'filtered')
    
    if (normalize == TRUE) {
      if (norm.method %in% c("rank.unit",
                             "rank.std",
                             "log.std",
                             "log.unit",
                             "log.clr",
                             "std",
                             "pass")) {
        sc.obj.filt <- normalize.features(sc.obj.filt,
                                          norm.method = "log.clr",
                                          norm.param = list(log.n0 = 1e-6))
      }
    }
    return(sc.obj.filt)
  }

sc.train <- function(sc.obj,
                     folds,
                     resample,
                     stratify = TRUE,
                     model,
                     fs = FALSE) {
  sc.obj <- create.data.split(
    sc.obj,
    num.folds = folds,
    num.resample = resample,
    stratify = stratify
  )
  
  sc.obj <- train.model(sc.obj,
                        method = model,
                        perform.fs = fs)
  return(sc.obj)
}

sc.test  <- function(sc.obj) {
  sc.obj <- make.predictions(sc.obj) %>%
    evaluate.predictions()
  return(sc.obj)
}

sc.subset <- function(sc.obj, keep) {
  warning(
    "samples are filtered after filtering and normalization,
          this may cause issues with certain analyses",
    call. = FALSE
  )
  sc.obj.filt <- sc.obj
  # filt phyloseq object
  sc.obj.filt@phyloseq <- prune_samples(samples = keep, x = sc.obj.filt@phyloseq)
  # filt filtered and normalized feature matrices
  sc.obj.filt@filt_feat$filt.feat <- sc.obj.filt@filt_feat$filt.feat[, keep]
  sc.obj.filt@norm_feat$norm.feat <- sc.obj.filt@filt_feat$filt.feat[, keep]
  # filter label list
  sc.obj.filt@label$label <- sc.obj.filt@label$label[keep]
  return(sc.obj.filt)
}

sc.plot <- function(sc.obj, basefname, plot_feat=F, colours = NULL) {

  if (is.null(colours)) {
    colours <- c("red", "blue")
  }
  model.evaluation.plot(
    sc.obj,
    fn.plot = paste0("figures/", basefname, "_ROC.pdf"),
    colours = colours,
    show.all = T
  )
  
  if (plot_feat == T) {
    model.interpretation.plot(
      sc.obj,
      max.show = 20, 
      fn.plot = paste0("figures/", basefname, "_features.pdf"),
      color.scheme = "RdBu"
    )
  }
}

sc.metrics <- function(sc.obj, metric="AUCROC") {
  
  auc_roc <- sc.obj@eval_data$auroc.all
  auc_pr <- sc.obj@eval_data$auprc.all
  
  if (metric == "AUCROC") {
    return(auc_roc)
  }
  else if (metric == "AUCPR") {
    return(auc_pr)
  }
  else if (metric == "both") {
    return_list <- list("AUCROC" = auc_roc, "AUCPR" = auc_pr)
    return(return_list)
  }
  
}