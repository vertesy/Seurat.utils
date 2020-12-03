######################################################################
# Custom.Functions.for.Slingshot.R
######################################################################
# source('~/GitHub/Packages/Seurat.utils/Custom.Functions.for.Slingshot.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Custom.Functions.for.Slingshot.R"))

# ------------------------

# require(ggbeeswarm)
# require(ggthemes)

#' Assign a color to each cell based on some value
#'
#' @param cell_vars Vector indicating the value of a variable associated with cells.
#' @param pal_fun Palette function that returns a vector of hex colors, whose
#' argument is the length of such a vector.
#' @param ... Extra arguments for pal_fun.
#' @return A vector of hex colors with one entry for each cell.

cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}


# ggplot for slinshot by @HectorRDB ------------------------
# https://github.com/kstreet13/slingshot/issues/73#issuecomment-585376827

### Point on curve function ----
points_on_curve <- function(curve, lambda, ...) {
  UseMethod("points_on_curve", curve)
}

points_on_curve.principal_curve <- function(curve, lambda, ...) {
  if (nrow(curve$s) == length(curve$lambda)) { # didn't use approx_points
    S <- apply(curve$s, 2, function(sjj) {
      return(approx(
        x = curve$lambda[curve$ord],
        y = sjj[curve$ord],
        xout = lambda, ties = "ordered"
      )$y)
    })
  } else {
    if (all(curve$ord == seq_along(curve$lambda))) { # used approx_points
      curvelambda <- seq(min(curve$lambda), max(curve$lambda), length.out = nrow(curve$s))
      S <- apply(curve$s, 2, function(sjj) {
        return(approx(
          x = curvelambda,
          y = sjj,
          xout = lambda, ties = "ordered"
        )$y)
      })
    }
  }
  return(S)
}

points_on_curve.SlingshotDataSet <- function(curve, lambda, ...) {
  locs <- lapply(slingCurves(curve), function(crv) {
    points_on_curve(crv, lambda, ...)
  })
  locs <- do.call('rbind', locs)
  colnames(locs) <- paste0("dim", seq_len(ncol(locs)))
  return(as.data.frame(locs))
}

### Extend ggplot function

#' Plot the gene in reduced dimension space
#'
#' @param sds The output from a lineage computation
#' @param col The assignation of each cell to a label. If none is provided,
#' default to the cluster labels provided when creating the \code{\link{SlingshotDataSet}}
#' @param title Title for the plot.
#' @param reduction "UMAP"
#' @param titleFsize title font size, def: 20
#' @param line.colors line color
#' @param lineSize Size of the curve lineages. Default to 1.
#' @param ... Other options passed to \code{\link{geom_point}}
#' @return A \code{\link{ggplot}} object
#' @examples
#' data('slingshotExample')
#' sds <- slingshot(rd, cl)
#' gg_plot(sds)
#'
#' ## Change point size and transparency
#' gg_plot(sds, size = 2, alpha = .5)
#'
#' ## Use grey background
#'
#' gg_plot(sds) + theme_grey()
#'
#' ## Color by gene expression
#' gene_count <- sample(0:10, nrow(reducedDims(sds)), replace = TRUE)
#' gg_plot(sds, col = gene_count)
#'
#' ## Add a marker of pseudotime
#' gg_plot(sds) + geom_point(data = points_on_curve(sds, 10), size = 3)
#' @importFrom slingshot slingPseudotime slingCurves reducedDim slingClusterLabels
#' @import ggplot2
#' @export

gg_plot <- function(sds, col = NULL, title = NULL, lineSize = 1, reduction = "UMAP"
                    , titleFsize = 20
                    , line.colors = gray.colors(n = length(sds@curves), start = 0, end = .6, alpha = 1 )
                    , ...) {
  rd <- reducedDim(sds)

  if (is.null(col)) {
    cl <- slingClusterLabels(sds)
    if ("matrix" %in% is(cl)) {
      cl <- apply(cl, 1, which.max)
      cl <- as.character(cl)
    }
  } else {
    cl <- col
  }

  # Getting the main plot
  df <- data.frame(dim1 = rd[, 1], dim2 = rd[, 2], col = cl)
  p <- ggplot(df, aes(x = dim1, y = dim2, col = col)) +
    geom_point(...) +
    theme_classic() +
    labs(title = title, col = "") +
    theme(plot.title = element_text(size = titleFsize)) # , face = "bold"
  # Adding the curves
  for (i in seq_along(slingCurves(sds))) {
    curve_i <- slingCurves(sds)[[i]]
    curve_i <- curve_i$s[curve_i$ord, ]
    colnames(curve_i) <- c("dim1", "dim2")
    p <- p + geom_path(data = as.data.frame(curve_i), col = line.colors[i], size = 1) +
      labs(x = paste(reduction, 1), y = paste(reduction, 1))
  }
  return(p)
}


# ------------------------

# plotting
#' @title Plot Gene Expression by Pseudotime
#' @name plotFittedGenePseudotime
#' @aliases plotFittedGenePseudotime
#'
#' @description Show the gene expression pattern for an individual gene along
#' lineages inferred by \code{\link{slingshot}}.
#'
#' @param data an object containing \code{\link{slingshot}} output, either a
#'   \code{\link{SlingshotDataSet}} or a \code{\link{SingleCellExperiment}}
#'   object.
#'
#' @export
setGeneric(
  name = "plotFittedGenePseudotime",
  signature = c('data'),
  def = function(data, ...) {
    standardGeneric("plotFittedGenePseudotime")
  }
)

setMethod(
  f = "plotFittedGenePseudotime",
  signature = signature(data = "SlingshotDataSet"),
  definition = function(data, gene, exprs, lcol = 1:4,
                        loess = TRUE, loessCI = TRUE, ...) {
    if (length(gene) > 1 & is.numeric(gene)){
      y <- gene
      genename <- deparse(substitute(gene))
    }
    if (length(gene) == 1){
      y <- exprs[gene, ,drop=FALSE][1,]
      genename <- gene
    }
    pst <- slingPseudotime(data)
    w <- slingCurveWeights(data)
    L <- length(slingLineages(data))

    # par(mfrow = c(L,1))
    i = 0
    for(l in seq_len(L)){
      i = i +1
      # print(l)
      if (l == 1) {
        plot(pst[,l], y, xlab = 'Pseudotime', ylab = 'Expression', cex = 0,
             main=paste(genename, ', Lineage ',l, sep=''), ...)
      }
      if (loess | loessCI){
        l <- loess(y ~ pst[,l], weights = w[,l])
      }
      if (loessCI){
        pl <- predict(l, se=TRUE, )
        polygon(c(l$x[order(l$x)],rev(l$x[order(l$x)])),
                c((pl$fit+qt(0.975,pl$df)*pl$se)[order(l$x)],
                  rev((pl$fit-qt(0.975,pl$df)*pl$se)[order(l$x)])),
                border = NA, col = rgb(0,0,0,.3))
      }
      if (loess){
        lines(l$x[order(l$x)], l$fitted[order(l$x)], lwd=2, col = lcol[i])
      }
    }
    # par(mfrow = c(1,1))
    invisible(NULL)
  }
)
# plotFittedGenePseudotime(data = sds, gene ="SST", expr = EXPR, loessCI=T
#                          , col = colz, pch = 20, panel_first = grid(NULL) )



# ------------------------
# ------------------------


