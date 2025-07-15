library(tidyverse)
library(readr)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(factoextra)

# ================================== #
#### Metabolomics Data Processing ####
# ================================== #

### Reading and organising => dataset specific ###
#' Adjusting sample names to naming convention in "Einwagen.xlsx"
adjust.samples.to.weights <- function(sample_name) {
  new_name <- gsub("^[0-9]+_" ,"", sample_name)
  new_name <- gsub("_+[a-z]+-[0-9]+$", "", new_name)
  new_name <- gsub("_[0-9]$", "", new_name)
  return(toupper(new_name))
}

# NOTE: out dated since Adam send a new metadata file with matching names
#' Adjusting sample names to naming convention in "combined_all_susceptible.csv"
adjust.samples.to.meta <- function(sample_name) {
  new_name <- str_replace(sample_name, "FC", "fp")
  parts <- strsplit(new_name, "_")[[1]]
  return(paste(toupper(parts[2]), tolower(parts[1]),
               sep="_"))
}

#'
read.raw.data <- function(file, delim="\t", skip=3,
                          sample_pattern = "_AV[0-9]+_|_av[0-9]+_|_AVI[0-9]+_") {
  full <- read_delim(file, delim=delim,
                     skip=skip)
  sample.cols <- grepl(sample_pattern, colnames(full))
  return(list(info = full[,!sample.cols],
              data = full[,sample.cols]))
}


### Normalisation => general ###
# Normalising to sample weight
weight.normalise <- function(data, weights, factor=1,
                             sample.orient="column") {
  switch(sample.orient,
         "column"={data_ <- as.matrix(data)},
         "row"={data_ <- t(as.matrix(data))},
         stop(paste("sample.orient must be 'column' or 'row' not", sample.orient))
  )

  if (ncol(data_) != length(weights)) stop("Number of weights must be equal to the number of samples")
  if (sum(names(weights) %in% colnames(data_)) != length(weights)) {
    stop("Sample names in the weight vector must match sample names in the data")
  }

  normed <- as.matrix(sapply(names(weights), function(w) data_[,w]/weights[w]*factor))
  colnames(normed) <- names(weights)
  return(normed)
}

# Actual Data Processing
#'
quotient.normalise <- function(data, sample.orient="column",
                               reference=NULL, na.rm=TRUE) {
  switch(sample.orient,
         "column"={raw <- data},
         "row"={raw <- t(data)},
         stop(paste("sample.orient must be 'column' or 'row' not", sample.orient))
  )
  # reference sample
  if (is.null(reference)) {
    ref <- apply(raw, 1, median, na.rm=na.rm)
  } else {
    if (is.null(dim(reference))) {
      assertthat::are_equal(
        length(reference) == nrow(raw),
        "'reference' must have the same number of elements than the samples in data"
      )
      ref <- reference
    } else {
      if (sample.orient == "column") {
        assertthat::are_equal(
          nrow(reference) == nrow(raw),
          "'reference' must have the same number of elements than the samples in data"
        )
        ref <- apply(reference, 1, median, na.rm=na.rm)
      } else {
        assertthat::are_equal(
          ncol(reference) == nrow(raw),
          "'reference' must have the same number of elements than the samples in data"
        )
        ref <- apply(reference, 2, median, na.rm=na.rm)
      }
    }
  }
  # dilution factors
  dilutions <- apply(raw, 2, function(x){median(x/ref, na.rm=na.rm)})
  # applying dilutions
  diluted <- as.data.frame(sapply(1:ncol(raw), function(i){raw[,i]/dilutions[i]}))
  dimnames(diluted) <- dimnames(raw)

  if (sample.orient == "column") {
    return(diluted)
  }
  return(t(diluted))
}

#'
safe.log2 <- function(data, constant=1e-10) {
  return(log2(data + constant))
}

# Genearlised logarithm
glog <- function(data, l=1e-5, log.fun=log2) {
  log.fun(data + sqrt(data^2 + l))
}

#'
scale.z.scores <- function(data, sample.orient="column",
                           na.rm=TRUE) {
  switch(sample.orient,
         "column"={
           return(t(apply(data, 1,
                          function(x){
                            (x-mean(x, na.rm=na.rm))/sd(x, na.rm=na.rm)
                          }))
                  )
         },
         "row"={
           return(apply(data, 2,
                        function(x){
                          (x-mean(x, na.rm=na.rm))/sd(x, na.rm=na.rm)
                        })
                  )
           },
         stop(paste("sample.orient must be 'column' or 'row' not", sample.orient))
  )
}

### Missing Value Imputation ###
#'
impute.missing.values <- function(data, remove.by.missing=.5,
                                  min.fraction=5,
                                  sample.orient="column") {
  switch(sample.orient,
         "column"={imp <- data[rowSums(is.na(data) | data == 0) < remove.by.missing,]},
         "row"={imp <- t(data[colSums(is.na(data) | data == 0) < remove.by.missing,])},
         stop(paste("sample.orient must be 'column' or 'row' not", sample.orient))
  )

  # TODO: add/use more sophisticated methods here?
  # e.g. correlation/regression based, kNN, PPCA/BPCA etc.?
  imp <- apply(data, 1,
               function(x){
                 minimum <- min(x[x > 0], na.rm=TRUE)
                 x[is.na(x) | x == 0] <- minimum/min.fraction
                 return(x)
               }
          )
  return(t(imp))
}

### Feature Filtering ###
# Auxiliaries for filter.features #
# Relative Standard Deviation
rsd <- function(x) abs(sd(x, na.rm=TRUE)/mean(x, na.rm=TRUE))
# Non-Parametric Standard Deviation
nprsd <- function(x) abs(mad(x, na.rm=TRUE)/median(x, na.rm=TRUE))

#'
filter.features <- function(data, method="nprsd", filter.absolute=NULL,
                            filter.quantile=.2, sample.orient="column",
                            return.mask=FALSE) {
  switch(sample.orient,
         "column"={data_ <- data},
         "row"={data_ <- t(data)},
         stop(paste("sample.orient must be 'column' or 'row' not", sample.orient))
  )

  # TODO: other methods? => these methods are all somehow variation-based
  switch(method,
         "rsd"={filters <- apply(data_, 1, rsd)},
         "nprsd"={filters <- apply(data_, 1, nprsd)},
         # Interquantile Range
         "iqr"={filters <- apply(data_, 1, IQR)},
         stop(paste("sample.orient must be one of [rsd, nprsd, iqr] not", sample.orient))
  )

  filter.ranks <- rank(filters)
  if (!is.null(filter.absolute)) {
    # filtering for the n most variable features
    mask <- filter.ranks > filter.absolute
  } else {
    # filtering for the q*N (feature number) most variable features
    mask <- filter.ranks > filter.quantile*nrow(data_)
  }
  if (return.mask) return(mask)
  return(data_[mask,])
}


### Plotting to compare raw to processed data ###
#' @param data : data.frame or matrix with sample in columns and feature names in rownames
#' @param groups : named vector of sample groups, names must be equal to colnames of data
plot.sample.distributions <- function(data, subset=NULL, groups=NULL,
                                      ylab="Concentration") {
  g_art <- FALSE
  # subset: only show random samples
  if (!is.null(subset)) {
    if (is.integer(subset)) {
      mask <- sample(1:ncol(data), size=subset)
    } else {
      mask <- subset
    }
  } else {
    mask <- rep(TRUE, ncol(data))
  }

  if (is.null(groups)) {
    groups <- rep("A", ncol(data))
    names(groups) <- colnames(data)
    g_art <- TRUE
  }

  pivot_longer(data[,mask], colnames(data[,mask]),
               names_to="Sample", values_to="Concentration") %>%
    add_column(Group=groups[.$Sample]) -> plot.data

  plot <- ggplot(plot.data, aes(x=Sample, y=Concentration, fill=Group)) +
            geom_boxplot() +
            theme_pubr() +
            scale_fill_d3("category10") +
            ylab(ylab)

  if (g_art) return(plot + theme(legend.position="none"))
  return(plot)
}

#' @param data : data.frame or matrix with sample in columns and feature names in rownames
plot.feature.distributions <- function(data, subset=NULL, ylab="Concentration") {
  if (is.null(subset)) {
    mask <- rep(TRUE, nrow(data))
  } else {
    if (is.integer(subset)) {
      mask <- sample(1:nrow(data), size=subset)
    } else {
      mask <- subset
    }
  }
  plot.data <- as.data.frame(t(data[mask,]))
  return(pivot_longer(plot.data, colnames(plot.data),
                      names_to="Feature", values_to="Concentration") %>%
            ggplot(aes(x=Feature, y=Concentration)) +
              geom_boxplot() +
              theme_pubr() +
              ylab(ylab))
}

#' @param groups : named vector of sample groups, names must be equal to colnames of data
#' @param shape : list, first element variable name, second element analogous to groups
plot.pca <- function(pca, groups, shape=NULL, pcs=c(1,2), ellipse=TRUE,
                     show.labels=FALSE, size=3, facet_var=NULL,
                     multi_facet=FALSE, add_pval=FALSE) {
  var.expl <- pca$sdev^2/sum(pca$sdev^2)
  
  if (is.null(facet_var)) {
    facet_var <- rep(NA, nrow(pca$x))
    names(facet_var) <- rownames(pca$x)
    facet_ <- FALSE
  } else {
    facet_var <- facet_var[rownames(pca$x)]
    facet_ <- TRUE
  }
  if (is.null(shape)) {
    plot_data <- data.frame(x=pca$x[,pcs[1]], y=pca$x[,pcs[2]],
                            Group=groups[rownames(pca$x)],
                            label=rownames(pca$x),
                            facet_var=facet_var[rownames(pca$x)])
    point_aes <- aes(x=x, y=y, colour=Group, fill=Group)
    p <- ggplot()
  } else {
    plot_data <- data.frame(x=pca$x[,pcs[1]], y=pca$x[,pcs[2]],
                            Group=groups[rownames(pca$x)],
                            Shape=shape[[2]][rownames(pca$x)],
                            label=rownames(pca$x),
                            facet_var=facet_var[rownames(pca$x)])
    point_aes <- aes(x=x, y=y, colour=Group, fill=Group,
                    shape=Shape)
    p <- ggplot() +
          scale_shape_discrete(name=shape[[1]])
  }

  if (ellipse) {
    p <- p + stat_ellipse(data=plot_data, aes(x=x, y=y, fill=Group),
                          inherit.aes=FALSE,
                          geom="polygon", alpha=.2) +
             scale_fill_d3("category10")
  }

  if (show.labels) p <- p + geom_text(data=plot_data, aes(label=label, x=x, y=y))
  else {
    if (!multi_facet) {
        p <- p + geom_point(data=plot_data, size=size,
                            mapping=point_aes)
    }
    else {
      if (is.null(facet_var) | all(is.na(facet_var))) {
        stop("'facet_var' cannot be null when 'multi_facet' is TRUE")
      }
      p <- p +
        geom_point(data=mutate(plot_data, Group=NULL, facet_var=NULL),
                   colour="grey85", size=size,
                   aes(x=x, y=y),
                   inherit.aes=FALSE) +
        geom_point(data=plot_data, size=size,
                   mapping=point_aes) +
        facet_wrap("facet_var")
      facet_ <- FALSE
    }
  }
  p +
    xlab(paste0("PC", pcs[1], " (", round(var.expl[pcs[1]], 3)*100, "%)")) +
    ylab(paste0("PC", pcs[2], " (", round(var.expl[pcs[2]], 3)*100, "%)")) +
    scale_color_d3("category10") +
    theme_pubr() +
    theme(legend.position="right") -> p

  if (add_pval) {
    if (multi_facet) stop("'add_pval' cannot be used with 'multi_facet'")
    pdata <- plot_data[,c(1,2)]
    # TODO: adapt for multiple groups
    groups <- levels(plot_data$Group)
    if (is.null(groups)) groups <- unique(plot_data$Group)
    htest <- Hotelling::hotelling.test(
      x=plot_data$x[plot_data$Group == groups[1]],
      y=plot_data$y[plot_data$Group == groups[2]]
    )
    p +
      annotate(geom='text', label=paste("p-value:", round(htest$pval, 4)),
               x=-Inf, y=Inf,
               hjust=0, vjust=1)
  }
  
  if (facet_) p <- p + facet_wrap("facet_var")

  return(p)
}

#' @param show.type : "coord" or "cor" - determines whether correlations between PCs or coords on PCs are shown
plot.pca.vars <- function(pca, show.type="coord",
                          f.names=NULL, colour=NULL,
                          pcs=c(1,2), scale=1,
                          return.data=FALSE) {
  if (!(show.type %in% c("coord", "cor"))) stop(paste("show.type must be 'coord' or 'cor', not", show.type))

  var.expl <- pca$sdev^2/sum(pca$sdev^2)

  var <- get_pca_var(pca)
  data <- as.data.frame(var[[show.type]])*scale
  data <- data[,pcs]
  colnames(data) <- c("x", "y")

  if (is.null(f.names)) data$Feature <- rownames(data)
  else data$Feature <- f.names

  if (!is.null(colour)) {
    data$colour <- colour[[2]][rownames(data)]
    ggplot(data, aes(x=x, y=y, colour=colour,
                     label=Feature)) +
      geom_text() +
      geom_segment(aes(x=0, y=0, xend=x, yend=y),
                   arrow=arrow(length=unit(0.25, 'cm'))) -> p
  } else{
    ggplot(data, aes(x=x, y=y, label=Feature)) +
      geom_text() +
      geom_segment(aes(x=0, y=0, xend=x, yend=y),
                   arrow=arrow(length=unit(0.25, 'cm'))) -> p
  }

  if (return.data) return(data)

  p <- p +
        xlab(paste0("PC", pcs[1], " (", round(var.expl[pcs[1]], 3)*100, "%)")) +
        ylab(paste0("PC", pcs[2], " (", round(var.expl[pcs[2]], 3)*100, "%)")) +
        theme_pubr()
  return(p)
}


#' @param groups : named vector of sample groups, names must be equal to colnames of data
#' @param shape : list, first element variable name, second element analogous to groups
plot.umap <- function(um, groups, shape=NULL) {
  if (is.null(shape)) {
    p <- ggplot(data.frame(x=um$layout[,1], y=um$layout[,2],
                           Group=groups[rownames(um$layout)]),
                aes(x=x, y=y))
  } else {
    p <- ggplot(data.frame(x=um$layout[,1], y=um$layout[,2],
                           Group=groups[rownames(um$layout)],
                           Shape=shape[[2]][rownames(um$layout)]),
                aes(x=x, y=y, shape=Shape)) +
      scale_shape_discrete(name=shape[[1]])
  }

  p +
    geom_point(aes(color=Group)) +
    stat_ellipse(aes(fill=Group),
                 geom="polygon", alpha=.2) +
    xlab("UMAP1") +
    ylab("UMAP2") +
    scale_color_d3("category10") +
    scale_fill_d3("category10") +
    theme_pubr() +
    theme(legend.position="right")
}

# NOTE: this computes the fold-change on log-data NOT on un-transformed data!
fold.change <- function(data, groups, group.order=NULL) {
  if (is.null(names(groups))) stop("groups must be a named vector")
  if (is.null(group.order)) {
    uni.groups <- sort(unique(groups))
  } else {
    uni.groups <- group.order
  }
  n.groups <- length(uni.groups)

  # NOTE: only geometric mean if log-data
  if (is.null(dim(data))) {
      geom.means <- sapply(uni.groups,
                           function(g) {
                             sub <- data[names(groups[groups == g])]
                             mean(sub, na.rm=TRUE)
                           })

  } else {
      geom.means <- sapply(uni.groups,
                           function(g) {
                             sub <- data[,names(groups[groups == g])]
                             if (is.null(dim(sub))) {
                               warning(
                                 paste0("Only one sample in group ", g, "!")
                               )
                               return(sub)
                             }
                             rowMeans(sub, na.rm=TRUE)
                           })
  }
  if (n.groups == 2) {
    fc <- geom.means[,1] - geom.means[,2]
    comp <- paste(colnames(geom.means), collapse="_")
    return(list(fold.changes=fc,
                comparison=comp))
  } else {
    comps <- sapply(1:(n.groups-1),
                    function(i) {
                      sapply((i+1):n.groups,
                             function(j) return(geom.means[,i] - geom.means[,j])
                      )
                    })

    fcs <- comps[[1]]
    for (i in 2:length(comps)) fcs <- cbind(fcs, comps[[i]])

    comp.names <- sapply(1:(n.groups-1),
                         function(i) {
                           i.col <- colnames(geom.means)[i]
                           paste(i.col, colnames(geom.means)[(i+1):n.groups],
                                 sep="_")
                         })
    colnames(fcs) <- unlist(comp.names)
    rownames(fcs) <- rownames(data)

    return(list(fold.changes=fcs,
                comparison=unlist(comp.names)))
  }
}

log.fold.change <- function(data, groups, log.fun=glog) {
  log.data <- log.fun(data)
  return(fold.change(log.data, groups))
}

# ============================= #
#### OTU Count Normalisation ####
# ============================= #
### Naming and Stuff ###
#' Filtering based on minimum relative abundance in a minimum fraction of samples
filter.otu <- function(comp.data, abundance=.01, frac=.3) {
  rownames(comp.data)[rowSums(comp.data*100 >= abundance) >= frac*ncol(comp.data)]
}

# NOTE: tailored to given data format
read.otu.data <- function(otu.file, meta.file=NULL, filter=TRUE,
                          abundance=.01, frac=.3, verbose=FALSE) {
  # OTU counts
  otu <- read_delim(otu.file, delim="\t")
  colnames(otu)[2:ncol(otu)] <- gsub("^[0-9]+\\.", "", colnames(otu)[2:ncol(otu)])
  otu <- as.data.frame(otu)
  rownames(otu) <- otu$`#OTUId`

  otu.taxa <- data.frame(t(sapply(otu$taxonomy, function(x) strsplit(x, ";")[[1]])))
  colnames(otu.taxa) <- c("Kingdom", "Phylum", "Class",
                          "Order", "Family", "Genus")
  rownames(otu.taxa) <- rownames(otu)

  # extracting counts only
  counts <- otu[,startsWith(colnames(otu), "AV")]
  # turning into compositions
  composition <- data.frame(apply(counts, 2, function(x) x/sum(x)))

  if (filter) {
    otus.kept <- filter.otu(composition, abundance, frac)

    cat(length(otus.kept), "of", nrow(composition),
        "OTUs retained after filtering\n")
    if (verbose) {
      cat("Discarded OTUs:\n", paste(setdiff(rownames(composition), otus.kept), collapse=","))
    }

    counts <- counts[otus.kept,]
    # re-compute compositions => changed by OTU removal
    composition <- data.frame(apply(counts, 2, function(x) x/sum(x)))
    otu.taxa <- otu.taxa[otus.kept,]
  }

  if (is.null(meta.file)) {
    return(list(counts=counts,
                composition=composition,
                taxonomy=otu.taxa))
  } else {
    # metadata
    otu.meta <- as.data.frame(read_delim(meta.file, delim="\t"))
    rownames(otu.meta) <- otu.meta$Sample
    return(list(counts=counts,
                composition=composition,
                taxonomy=otu.taxa,
                meta=otu.meta))
  }
}

subset.samples <- function(data, to.remove) {
  mask <- !(colnames(data) %in% to.remove)
  return(data[,mask])
}
