library(mixOmics)

plot_samples <- function(plsda.obj, group.name="Group",
                         shape=NULL, plot.ellipse=TRUE,
                         size=2, geom="text", ...) {
  pdf(file=NULL)
  sample.plot <- plotIndiv(plsda.obj, ...)
  dev.off()
  
  if (!is.null(shape)) {
    sample.plot$df$shape <- shape[[2]]
    p <- ggplot(data=sample.plot$df,
                aes_string(x="x", y="y",
                           colour="group",
                           shape="shape"),
                size=size) +
            labs(shape=shape[[1]],
                 colour=group.name)
  } else {
    p <- ggplot(data=sample.plot$df,
                aes_string(x="x", y="y",
                           colour="group",
                           label="names"),
                size=size) +
            labs(colour=group.name)
  }
  
  if (plot.ellipse) {
    p <- p + stat_ellipse(data=sample.plot$df,
                          geom="polygon",
                          alpha=.15,
                          aes_string(x="x", y="y",
                                     fill="group"),
                          inherit.aes = FALSE,
                          show.legend=FALSE)
  }
  
  if (geom == "point") p <- p + geom_point(size=size)
  if (geom == "text") p <- p + geom_text(size=size)

  
  p +
    scale_colour_d3("category10") +
    scale_fill_d3("category10") +
    theme_pubr() +
    xlab(sample.plot$graph$labels$x) +
    ylab(sample.plot$graph$labels$y)
}

plot_variable_subset <- function(plsda.obj, add.columns=NULL,
                                 colour=NULL, shape=NULL,
                                 label=NULL, geom="point",
                                 x_threshold=.975, y_threshold=.975,
                                 quantiles=TRUE, relation="OR",
                                 return.full=FALSE, size=2) {
  
  pdf(file=NULL)
  var.plot <- plotVar(plsda.obj)
  dev.off()
  
  if (quantiles) {
    x.cutoff <- quantile(abs(var.plot$x), x_threshold)
    y.cutoff <- quantile(abs(var.plot$y), y_threshold)
  } else {
    x.cutoff <- x_threshold
    y.cutoff <- y_threshold
  }
  
  
  if (!is.null(add.columns)) {
    var.plot <- cbind(var.plot, add.columns[rownames(var.plot),])
  }
  
  if (relation == "AND") {
    select.var.plot <- filter(var.plot,
                              abs(x) >= x.cutoff & abs(y) > y.cutoff)
  } else {
    select.var.plot <- filter(var.plot,
                              abs(x) >= x.cutoff | abs(y) > y.cutoff)
  }
  
  if (dim(select.var.plot)[1] == 0) stop("No features left after filtering. Choose different cutoffs or ")

  ggplot() +
    geom_path(data=circle.path(), aes(x=x, y=y)) +
    geom_path(data=circle.path(r=.5), aes(x=x, y=y)) -> p
  
  if (geom == "point") {
   p <- p + geom_point(data=select.var.plot,
                       aes_string(x="x", y="y", colour=colour,
                                  shape=shape),
                       size=size)
  }
  if (geom == "text") {
    p <- p + geom_text(data=select.var.plot,
                        aes_string(x="x", y="y", colour=colour,
                                   label=label),
                        size=size)
  }
  
  p <- p +
        theme_pubr() +
        theme(legend.position="right")
  
  if (return.full) return(list(plot=p,
                               filtered.data=select.var.plot,
                               full.data=var.plot))
  return(list(plot=p,
              data=select.var.plot))
}

features_from_subset <- function(features, meta) {
  feat_sub <- features[features %in% rownames(meta)]
  return(meta[feat_sub,])
}
