#source("../scripts/Normalisation.R")
library(ggforce)
library(ggrepel)

# scaling is adapted from the factoextra package
variable.scaling <- function(pca, var.type="coord", pcs=c(1,2)) {
  sams <- pca$x[,pcs]
  # vars <- get_pca_var(pca)[[var.type]][,pcs]
  vars <- pca[[var.type]]
  
  x.fact <- (max(sams[,1])-min(sams[,1])/(max(vars[,1])-min(vars[,1])))
  y.fact <- (max(sams[,2])-min(sams[,2])/(max(vars[,2])-min(vars[,2])))
  
  return(min(x.fact, y.fact))
}

circle.path <- function(r=1, n=100){
  p <- seq(0, 2*pi, length.out=n)
  return(data.frame(x=r*cos(p),
                    y=r*sin(p)))
}

#' @param show.type : "coord" or "cor" - determines whether correlations between PCs or coords on PCs are shown
plot.pca.vars <- function(pca, show.type="coord",
                          f.names=NULL, colour=NULL,
                          pcs=c(1,2), scale=1,
                          return.data=FALSE,
                          add.columns=NULL) {
  if (!(show.type %in% c("coord", "cor", "contrib"))) stop(paste("show.type must be 'coord' or 'cor', not", show.type))
  
  var.expl <- pca$sdev^2/sum(pca$sdev^2)
  
  var <- get_pca_var(pca)
  data <- as.data.frame(var[[show.type]])
  if(show.type == "coord") data <- data*scale
  data <- data[,pcs]
  colnames(data) <- c("x", "y")
  
  if (is.null(f.names)) data$Feature <- rownames(data)
  else data$Feature <- f.names
  
  if (!is.null(add.columns)) {
    data <- cbind(data, add.columns[rownames(data),])
  }
  
  if (return.data) return(data)
  
  if (!is.null(colour)) {
    data$colour <- colour[[2]][rownames(data)]
    ggplot(data, aes(x=x, y=y, colour=colour)) +
      geom_text_repel(label=data$Feature) +
      geom_segment(aes(x=0, y=0, xend=x, yend=y),
                   arrow=arrow(length=unit(0.25, 'cm'))) -> p
  } else{
    ggplot(data, aes(x=x, y=y)) +
      geom_text_repel(label=data$Feature) +
      geom_segment(aes(x=0, y=0, xend=x, yend=y),
                   arrow=arrow(length=unit(0.25, 'cm'))) -> p
  }
  
  p <- p +
    xlab(paste0("PC", pcs[1], " (", round(var.expl[pcs[1]], 3)*100, "%)")) +
    ylab(paste0("PC", pcs[2], " (", round(var.expl[pcs[2]], 3)*100, "%)")) +
    theme_pubr()
  
  return(p)
}

euclid.dist <- function(x) sqrt(sum(x^2, na.rm=TRUE))
plot.feature.contrib <- function(pca, pcs=c(1,2), n=NULL) {
  contrib <- get_pca_var(pca)$contrib
  contrib.dist <- apply(contrib[,pcs], 1, euclid.dist)
  
  if (!is.null(n)) mask <- rank(-contrib.dist) <= n
  else mask <- rep(TRUE, nrow(contrib))
  ggplot(data.frame(Feature=rownames(contrib)[mask],
                    Contrib=contrib.dist[mask])) +
    geom_bar(aes(x=reorder(Feature, -Contrib), y=Contrib), stat="identity") +
    xlab("Feature") +
    ylab("Contribution")
}


### Heatmap Stuff ###
colour.pal.d3 <- function(groups, return.palette=FALSE) {
  uni.g <- unique(groups)
  n <- length(uni.g)
  if (n > 10) {
    if (n > 20) stop("max. 20 groups allowed")
    cs <- sapply(pal_d3("category20")(n), substr,
                 start=1, stop=7)
  } else {
    cs <- sapply(pal_d3("category10")(n), substr,
                 start=1, stop=7)
  }
  
  names(cs) <- uni.g
  
  if (return.palette) return(cs)
  return(cs[groups])
}

plot.heatmap <- function(data, sample.groups,
                         sample.title,
                         feature.groups=NULL,
                         feature.title=NULL,
                         colours=gplots::bluered(100),
                         mar=c(5, 2, 4, 4),
                         col.legend.pos=c(.9, 1.1),
                         row.legend.pos=c(.975, .7),
                         legend.cex=.8,
                         ...) {
  samp.groups <- as.character(sample.groups)
  feat.groups <- as.character(feature.groups)
  
  samp.pal <- colour.pal.d3(samp.groups, return.palette=TRUE)
  
  par(xpd=TRUE, mar=mar)
  
  if (!is.null(feature.groups)) {
    feat.pal <- colour.pal.d3(feat.groups, return.palette=TRUE)
    gplots::heatmap.2(data,
                      col=colours,
                      RowSideColors=feat.pal[feat.groups],
                      ColSideColors=samp.pal[samp.groups],
                      margins=c(5,11),
                      trace="none",
                      ...)
    legend(x=row.legend.pos[1], y=row.legend.pos[2],
           legend=names(samp.pal), fill=samp.pal, bty="n",
           title=sample.title, cex=legend.cex)
    legend(x=col.legend.pos[1], y=col.legend.pos[2],
           legend=names(feat.pal), fill=feat.pal,
           cex=legend.cex, bty="n", title=feature.title)
  } else {
    gplots::heatmap.2(data,
                      col=colours,
                      ColSideColors=samp.pal[samp.groups],
                      margins=c(5,11),
                      trace="none",
                      ...)
    legend(x=col.legend.pos[1], y=col.legend.pos[2],
           legend=names(samp.pal), fill=samp.pal, bty="n",
           title=sample.title)
  }
}


colour.pal.jco <- function(sample.groups, return.palette=TRUE) {
  uni.g <- unique(sample.groups)
  n <- length(uni.g)
  if (n > 10) {
    stop("Maximum of 10 values allowed")
  } else {
    cs <- sapply(pal_jco()(n), substr,
                 start=1, stop=7)
  }
  
  names(cs) <- uni.g
  
  if (return.palette) return(cs)
  return(cs[sample.groups])
}

colour.pal.nejm <- function(sample.groups, return.palette=TRUE) {
  uni.g <- unique(sample.groups)
  n <- length(uni.g)
  if (n > 10) {
    stop("Maximum of 10 values allowed")
  } else {
    cs <- sapply(pal_nejm()(n), substr,
                 start=1, stop=7)
  }
  
  names(cs) <- uni.g
  
  if (return.palette) return(cs)
  return(cs[sample.groups])
}


library(gtools)
# NOTE: gtools is required for heatmap.2x => throws error if not loaded
plot.multi.heatmap <- function(data, sample.groups,
                               feature.groups=NULL,
                               feature.title=NULL,
                               colours=gplots::bluered(100),
                               mar=c(5, 2, 4, 4),
                               col1.legend.pos=c(.9, 1.1),
                               col2.legend.pos=c(1.05, 1.1),
                               row.legend.pos=c(.975, .7),
                               ...) {
  
  if (length(sample.groups) != 2) stop("sample groups must be a list with two entries")
  
  if (any(names(sample.groups[[1]]) != names(sample.groups[[2]]))) {
    sample.groups[[2]] <- sample.groups[[2]][names(sample.groups[[1]])]
  }
  
  rnames <- names(sample.groups[[1]])
  
  sample.groups[[1]] <- as.character(sample.groups[[1]])
  sample.groups[[2]] <- as.character(sample.groups[[2]])
  
  samp.pal.1 <- colour.pal.d3(sample.groups[[1]], return.palette=TRUE)
  samp.pal.2 <- colour.pal.jco(sample.groups[[2]], return.palette=TRUE)
  
  samp.cols <- rbind(samp.pal.1[sample.groups[[1]]],
                     samp.pal.2[sample.groups[[2]]])

  dimnames(samp.cols) <- list(names(sample.groups), rnames)
  
  par(xpd=TRUE, mar=mar)
  
  if (!is.null(feature.groups)) {
    feat.groups <- as.character(feature.groups)
    feat.pal <- colour.pal.d3(feature.groups, return.palette=TRUE)
    heatmap.2x::heatmap.2x(data,
                           col=colours,
                           RowSideColors=as.matrix(feat.pal[feat.groups],
                                                   nrow=1),
                           ColSideColors=samp.cols,
                           margins=c(5,11),
                           trace="none",
                           ...)
    legend(x=col1.legend.pos[1], y=col1.legend.pos[2],
           legend=names(samp.pal.1), fill=samp.pal.1, bty="n",
           title=names(sample.groups)[1])
    legend(x=col2.legend.pos[1], y=col2.legend.pos[2],
           legend=names(samp.pal.2), fill=samp.pal.2, bty="n",
           title=names(sample.groups)[2])
    legend(x=row.legend.pos[1], y=row.legend.pos[2],
           legend=names(feat.pal), fill=feat.pal,
           cex=.8, bty="n", title=feature.title)
  } else {
    heatmap.2x::heatmap.2x(data,
                           col=colours,
                           ColSideColors=samp.cols,
                           margins=c(5,11),
                           trace="none",
                           ...)
    legend(x=col1.legend.pos[1], y=col1.legend.pos[2],
           legend=names(samp.pal.1), fill=samp.pal.1, bty="n",
           title=names(sample.groups)[1])
    legend(x=col2.legend.pos[1], y=col2.legend.pos[2],
           legend=names(samp.pal.2), fill=samp.pal.2, bty="n",
           title=names(sample.groups)[2])
  }
}

plot.3.heatmap <- function(data, sample.groups,
                               feature.groups=NULL,
                               feature.title=NULL,
                               colours=gplots::bluered(100),
                               mar=c(5, 2, 4, 4),
                               margins=c(5, 11),
                               col1.legend.pos=c(.9, 1.1),
                               col2.legend.pos=c(1.05, 1.1),
                               col3.legend.pos=c(1.05, .9),
                               row.legend.pos=c(.975, .7),
                               ...) {
  
  if (length(sample.groups) != 3) stop("sample groups must be a list with two entries")
  
  if (any(names(sample.groups[[1]]) != names(sample.groups[[2]]))) {
    sample.groups[[2]] <- sample.groups[[2]][names(sample.groups[[1]])]
  }
  
  rnames <- names(sample.groups[[1]])
  
  sample.groups[[1]] <- as.character(sample.groups[[1]])
  sample.groups[[2]] <- as.character(sample.groups[[2]])
  sample.groups[[3]] <- as.character(sample.groups[[3]])
  
  samp.pal.1 <- colour.pal.d3(sample.groups[[1]], return.palette=TRUE)
  samp.pal.2 <- colour.pal.jco(sample.groups[[2]], return.palette=TRUE)
  samp.pal.3 <- colour.pal.nejm(sample.groups[[3]], return.palette=TRUE)
  
  samp.cols <- rbind(samp.pal.1[sample.groups[[1]]],
                     samp.pal.2[sample.groups[[2]]],
                     samp.pal.3[sample.groups[[3]]])
  
  dimnames(samp.cols) <- list(names(sample.groups), rnames)
  
  par(xpd=TRUE, mar=mar)
  
  if (!is.null(feature.groups)) {
    feat.groups <- as.character(feature.groups)
    feat.pal <- colour.pal.d3(feature.groups, return.palette=TRUE)
    heatmap.2x::heatmap.2x(data,
                           col=colours,
                           RowSideColors=as.matrix(feat.pal[feat.groups],
                                                   nrow=1),
                           ColSideColors=samp.cols,
                           margins=margins,
                           trace="none",
                           ...)
    legend(x=col1.legend.pos[1], y=col1.legend.pos[2],
           legend=names(samp.pal.1), fill=samp.pal.1, bty="n",
           title=names(sample.groups)[1])
    legend(x=col2.legend.pos[1], y=col2.legend.pos[2],
           legend=names(samp.pal.2), fill=samp.pal.2, bty="n",
           title=names(sample.groups)[2])
    legend(x=col3.legend.pos[1], y=col3.legend.pos[2],
           legend=names(samp.pal.3), fill=samp.pal.3, bty="n",
           title=names(sample.groups)[3])
    legend(x=row.legend.pos[1], y=row.legend.pos[2],
           legend=names(feat.pal), fill=feat.pal,
           cex=.8, bty="n", title=feature.title)
  } else {
    heatmap.2x::heatmap.2x(data,
                           col=colours,
                           ColSideColors=samp.cols,
                           margins=margins,
                           trace="none",
                           ...)
    legend(x=col1.legend.pos[1], y=col1.legend.pos[2],
           legend=names(samp.pal.1), fill=samp.pal.1, bty="n",
           title=names(sample.groups)[1])
    legend(x=col2.legend.pos[1], y=col2.legend.pos[2],
           legend=names(samp.pal.2), fill=samp.pal.2, bty="n",
           title=names(sample.groups)[2])
    legend(x=col3.legend.pos[1], y=col3.legend.pos[2],
           legend=names(samp.pal.3), fill=samp.pal.3, bty="n",
           title=names(sample.groups)[3])
  }
}

rotation.to.coord <- function(pca.obj) {
  pca.obj$coord <- sapply(1:ncol(pca.obj$rotation),
                          function(i) pca.obj$rotation[,i]*pca.obj$sdev[i])
  dimnames(pca.obj$coord) <- dimnames(pca.obj$rotation)
  return(pca.obj)
}


plot_multi_pca <- function(pca, groups, shape=NULL, pcs=c(1,2),
                           facet_var=NULL, size=3, group_name=NULL) {
    var.expl <- pca$sdev^2/sum(pca$sdev^2)
    
    plot_data <- data.frame(x=pca$x[,pcs[1]], y=pca$x[,pcs[2]],
                            Group=groups[rownames(pca$x)],
                            label=rownames(pca$x))
    
    if (!is.null(facet_var)) {
      plot_data[[facet_var[[1]]]] <- facet_var[[2]][rownames(pca$x)]
    }
    
    if (is.null(shape)) {
      p <- ggplot(plot_data, aes(x=x, y=y, label=label))
    } else {
      plot_data[[shape[[1]]]] <- shape[[2]][rownames(pca$x)]
      p <- ggplot(plot_data,
                  aes_string(x="x", y="y", shape=shape[[1]],
                      label="label")) +
        scale_shape_discrete(name=shape[[1]])
    }
    # browser()
    p <- p +
      geom_point(data=mutate(plot_data, !!sym(facet_var[[1]]):=NULL),
                 colour="grey85", size=size) +
      geom_point(aes(colour=Group), size=size)
    
    if (!is.null(facet_var)) p <- p + facet_wrap(facet_var[[1]])
    
    p +
      xlab(paste0("PC", pcs[1], " (", round(var.expl[pcs[1]], 3)*100, "%)")) +
      ylab(paste0("PC", pcs[2], " (", round(var.expl[pcs[2]], 3)*100, "%)")) +
      scale_color_d3("category10") +
      theme_pubr() + 
      theme(legend.position="right") -> p
    
    if (!is.null(group_name)) p <- p + guides(colour=guide_legend(title=group_name))
    
    return(p)
  }


plot_multi_umap <- function(um, groups, shape=NULL, group_name=NULL,
                            facet_var=NULL, size=4) {
  plot_data <- data.frame(x=um$layout[,1], y=um$layout[,2],
                          Group=groups[rownames(um$layout)],
                          label=rownames(um$layout))
  
  if (!is.null(facet_var)) {
    plot_data[[facet_var[[1]]]] <- facet_var[[2]][rownames(um$layout)]
  }
  
  if (is.null(shape)) {
    p <- ggplot(plot_data, aes(x=x, y=y, label=label))
  } else {
    plot_data[[shape[[1]]]] <- shape[[2]][rownames(um$layout)]
    p <- ggplot(plot_data,
                aes_string(x="x", y="y", shape=shape[[1]],
                           label="label")) +
      scale_shape_discrete(name=shape[[1]])
  }
  
  p <- p +
    geom_point(data=mutate(plot_data, !!sym(facet_var[[1]]):=NULL),
               colour="grey85", size=size) +
    geom_point(aes(colour=Group), size=size)
  
  if (!is.null(facet_var)) p <- p + facet_wrap(facet_var[[1]])
  
  p +
    xlab("UMAP1") +
    ylab("UMAP2") +
    scale_color_d3("category10") +
    theme_pubr() + 
    theme(legend.position="right") -> p
  
  if (!is.null(group_name)) p <- p + guides(colour=guide_legend(title=group_name))
  
  return(p)
}
