#source("D:/Users/adam-//Dropbox/Haller_Lab/Data/Metabolomics/scripts/Normalisation.R")

library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(latex2exp)
library(EnhancedVolcano)

get_volcano_data <- function(p.values, log.fc,
                             min.fc=1, min.p=.25) {
  if (is.null(dim(p.values))) rows <- names(p.values)
  else rows <- rownames(p.values)

  sig_name <- sapply(rows,
                     function(r) {
                       pval <- p.values[r] < min.p
                       fc <- abs(log.fc[r]) > min.fc

                       if (pval & fc) {
                         return("p-value and logFC")
                       }
                       if (pval) {
                         return("p-value")
                       }
                       if (fc) {
                         return("logFC")
                       }
                       return("NS")
                     })

  data.frame("logFC"=log.fc[rows],
             "log_pvalue"=-log10(p.values),
             "sig"=sig_name,
             "Feature"=rows) -> volcano_data
  return(volcano_data)
}

# feature.annotation HAS to have a column "Feature" => used for joining
plot_multi_volcano <- function(p.values, log.fc,
                               feature.annotation,
                               colour.name, min.p=.05,
                               min.fc=2) {
  volcano.data <- get_volcano_data(p.values, log.fc,
                                   min.fc, min.p)

  full_join(as_tibble(volcano.data),
            feature.annotation) -> plot.data
    # TODO: use this for showing names for when not facet wrapping
    # add_column(FeatureName=sapply(1:nrow(.),
    #                               function(i) {
    #                                 if (.$sig[i] == "yes") return(.$Name[i])
    #                                 return("")
    #                                })) -> plot.data
  ggplot(data=plot.data,
         aes_string(x="logFC", y="log_pvalue",
                    colour=colour.name)) +
    annotate("rect", xmin=min.fc, ymin=-log10(min.p),
                  xmax=Inf, ymax=Inf,
              fill="#D62728",
              alpha=.3) +
    annotate("rect", xmin=-min.fc, ymin=-log10(min.p),
                  xmax=-Inf, ymax=Inf,
              fill="#D62728",
              alpha=.3) +
      geom_point(data=mutate(plot.data, !!sym(colour.name):=NULL),
                 colour="grey85") +
      geom_point() +
      facet_wrap(as.formula(paste0(colour.name, "~."))) +
      theme_pubr() +
      theme(legend.position="right") +
      xlab(bquote(~log[2]~ "FC")) +
      ylab(bquote(~-log[10]~ "p-value")) -> plot

  return(list(data=plot.data,
              plot=plot))
}

#
plot_volcano_names <- function(p.values, log.fc,
                               feature.names,
                               min.p=.05, min.fc=2,
                               dot.size=3, ...) {
  if (is.null(names(p.values))) stop("p.valuesmust be a named vector!")
  if (is.null(names(log.fc))) stop("log.fc must be a named vector!")
  if (is.null(names(feature.names))) stop("feature.names must be a named vector!")

  volcano.data <- get_volcano_data(p.values, log.fc,
                                   min.fc, min.p)
  volcano.data$Name <- feature.names[rownames(volcano.data)]
  # fa_mask <- grep("^FA |^LP|^P|^DG", volcano.data$Name)

  ggplot(data=volcano.data,
         aes(x=logFC, y=log_pvalue, colour=sig)) +
    # TODO: turn annotate rect into lines
    geom_vline(xintercept=c(-min.fc, min.fc),
               linetype=3) +
    geom_hline(yintercept=-log10(min.p),
               linetype=3) +
    geom_point(size=dot.size) +
    geom_text_repel(data=filter(#volcano.data[fa_mask,],
                                volcano.data,
                                sig == "p-value and logFC"),
                    aes(x=logFC, y=log_pvalue, label=Name),
                    show.legend=FALSE,
                    min.segment.length = unit(0, 'lines'),
                    nudge_y = .2, ...) +
    scale_colour_manual(values=c("NS"="black", "p-value and logFC"="red4",
                                 "p-value"="blue4" , "logFC"="green4")) +
    labs(colour="Significance") +
    theme_pubr() +
    xlab(bquote(~log[2]~ "FC")) +
    ylab(bquote(~-log[10]~ "p-value"))
}

get_padj <- function(data, groups,
                     feature.orient=1,
                     test.method=t.test,
                     p.adjust.method="BH") {
  raw_p <- apply(data, feature.orient,
                 function(x) tryCatch(test.method(x~groups)$p.value, error=function(x) return(1)))
  adj_p <- p.adjust(raw_p, p.adjust.method)
  names(adj_p) <- dimnames(data)[[feature.orient]]

  return(adj_p)
}

pairwise_ttest_metabo <- function(metabo, factor, case, control, 
                                  sample_col, meta, p_adjust="BH"){
  
  # initialise pval matrix 
  p.val <- matrix(NA, nrow=nrow(metabo), ncol=1, 
                  dimnames=list(row.names(metabo), "p.val"))
  
  featmat <- as.matrix(metabo)
  
  for (row in rownames(featmat)){
    x <- featmat[row, meta %>% 
                   filter(.data[[ factor ]] ==control) %>% pull(.data[[ sample_col ]])]
    y <- featmat[row, meta %>% 
                   filter(.data[[ factor ]] ==case) %>% pull(.data[[ sample_col ]])]
    
    
    p.val[row, ] <- t.test(x, y, exact=FALSE)$p.value
  }
  
  p.val.adj <- p.val %>% 
    as.data.frame() %>% 
    rstatix::adjust_pvalue("p.val", method = p_adjust) %>% 
    mutate(signif=case_when(p.val.adj <= 0.05 & p.val.adj > 0.01 ~ "*",
                            p.val.adj <= 0.01 & p.val.adj > 0.001 ~ "**",
                            p.val.adj <= 0.001 & p.val.adj > 0.0001 ~ "***",
                            p.val.adj <= 0.0001 ~ "****",
                            TRUE ~ "ns")) %>% 
    mutate(group1 = case,
           group2 = control)
  
  return(p.val.adj)
}

prettyVolcano <- function(res,
                          lfc_thresh = 1,
                          pval_thresh = 0.05,
                          downreg_col = "dodgerblue",
                          upreg_col = "firebrick1",
                          contrast,
                          gene_names,
                          xlim,
                          ylim,
                          ...) {
  keyvals <- ifelse(res$log2FoldChange < (lfc_thresh* -1) & res$padj < pval_thresh,
                    downreg_col, 
                    ifelse(
                      res$log2FoldChange > lfc_thresh & res$padj < pval_thresh,
                      upreg_col,
                      'black'
                    )
  )
  keyvals[is.na(keyvals)] <- "black"
  names(keyvals)[keyvals == upreg_col] <- contrast[[2]]
  names(keyvals)[keyvals == downreg_col] <- contrast[[1]]
  
  
  
  p <- EnhancedVolcano(
    res,
    lab = res[[gene_names]],
    x = "log2FoldChange",
    y = "padj",
    pCutoff = 0.1,
    FCcutoff = 1,
    xlim = xlim,
    ylim = ylim,
    title = NULL,
    subtitle = NULL,
    caption = NULL,
    colCustom = keyvals,
    colAlpha = 1,
    legendPosition = "None",
    ...
  )
  return(p)
}

# if (sys.nframe() == 0) {
#   setwd("~/PhD/ATF6")
#   load("ProcessedMetabolomeMicrobiome.RData")
#
#   total_samples <- c(colnames(fc.scaled), colnames(cc.scaled))
#
#   line <- sample.info[total_samples, "Mouse_Line"]
#   names(line) <- total_samples
#   age <- sample.info[total_samples,"Age"]
#   names(age) <- total_samples
#   genotype <- sample.info[total_samples,"Genotype"]
#   names(genotype) <- total_samples
#   phenotype <- sample.info[total_samples,"Phenotype"]
#   names(phenotype) <- total_samples
#   responder <- sample.info[total_samples,"Responder_Status"]
#   names(responder) <- total_samples
#
#   cc.av.mask <- line[colnames(cc.scaled)] == "AV"
#   cc.avi.mask <- line[colnames(cc.scaled)] == "AVI"
#   fc.av.mask <- line[colnames(fc.scaled)] == "AV"
#   fc.avi.mask <- line[colnames(fc.scaled)] == "AVI"
#
#   cc_av_5wk <- cc.scaled[,cc.av.mask & age[colnames(cc.scaled)] == "5wk"]
#   cc_av_tumour <- cc.scaled[,cc.av.mask & age[colnames(cc.scaled)] != "5wk"]
#
#   tum_adj_p <- get_padj(cc_av_tumour, phenotype[colnames(cc_av_tumour)],
#                         test.method=wilcox.test, p.adjust.method="BH")
#   tum_fc <- fold.change(data=cc.log.normed[,cc.av.mask & age[colnames(cc.log.normed)] != "5wk"],
#                         groups=phenotype[colnames(cc_av_tumour)])
#
#   ### Testing multi panel volcano
#   add_features <- metabolite.info[,c("Ontology", "Metabolite name")] %>%
#                     add_column(Feature=rownames(metabolite.info))
#   plot_multi_volcano(tum_adj_p, tum_fc$fold.changes,
#                      add_features,
#                      "Ontology") -> pmv
#
#   ### Testing volcano with labels
#   fn <- metabolite.info$`Metabolite name`
#   names(fn) <- rownames(metabolite.info)
#   plot_volcano_names(tum_adj_p, tum_fc$fold.changes,
#                      fn)
# }
