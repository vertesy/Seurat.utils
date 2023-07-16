# _________________________________________________________________________________________________
# Slingshot.Utils.R ----
# _________________________________________________________________________________________________
# source('~/GitHub/Packages/Seurat.utils/R/Slingshot.Utils.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Custom.Functions.for.Slingshot.R"))

# _________________________________________________________________________________________________

# require(ggbeeswarm)
# require(ggthemes)

# _________________________________________________________________________________________________
# https://github.com/kstreet13/slingshot/issues/73#issuecomment-585376827




# _________________________________________________________________________________________________
#' @title points_on_curve
#'
#' @description Helper function to visualize points on a curve in Slingshot.
#' @param curve The curve on which to visualize points. Default: PARAM_DEFAULT_VALUE.
#' @param lambda The parametric location of points along the curve. Default: PARAM_DEFAULT_VALUE.
#' @return A function call to points_on_curve based on the class of the curve.
#' @details This function uses method dispatch based on the class of the 'curve' argument to visualize points on the curve.
#' @export
points_on_curve <- function(curve, lambda, ...) {
  UseMethod("points_on_curve", curve)
}

# _________________________________________________________________________________________________
#' @title points_on_curve.principal_curve
#' @description Helper function to visualize points on a principal curve in Slingshot.
#' @param curve The principal curve on which to visualize points. Default: PARAM_DEFAULT_VALUE.
#' @param lambda The parametric location of points along the principal curve. Default: PARAM_DEFAULT_VALUE.
#' @return Coordinates of points on the principal curve.
#' @details This function interpolates the coordinates of points along the principal curve in Slingshot.
#' @seealso
#' \code{\link[slingshot]{SlingshotDataSet}}
#' @export
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

# _________________________________________________________________________________________________
#' @title points_on_curve.SlingshotDataSet
#' @description Helper function to visualize points on a SlingshotDataSet in Slingshot.
#' @param curve The SlingshotDataSet on which to visualize points. Default: PARAM_DEFAULT_VALUE.
#' @param lambda The parametric location of points along the curves in the SlingshotDataSet. Default: PARAM_DEFAULT_VALUE.
#' @return A dataframe of coordinates of points on the curves in the SlingshotDataSet.
#' @details This function interpolates the coordinates of points along the curves in a SlingshotDataSet.
#' @seealso
#' \code{\link[slingshot]{SlingshotDataSet}}
#' @export
points_on_curve.SlingshotDataSet <- function(curve, lambda, ...) {
  locs <- lapply(slingCurves(curve), function(crv) {
    points_on_curve(crv, lambda, ...)
  })
  locs <- do.call('rbind', locs)
  colnames(locs) <- paste0("dim", seq_len(ncol(locs)))
  return(as.data.frame(locs))
}


# _________________________________________________________________________________________________
# Colors and visuals ----
# _________________________________________________________________________________________________


#' @title ggplotColours
#' @description Generate ggplot colours for slingshot
#' @param n PARAM_DESCRIPTION, Default: 6
#' @param h height of the plot, Default: c(0, 360) + 15
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}



# _________________________________________________________________________________________________
### Extend ggplot function

#' @title gg_plot
#' @description An adjusted ggplot function for visualizing Slingshot datasets.
#' @param sds A Slingshot dataset. Default: PARAM_DEFAULT_VALUE.
#' @param col The color to be used for points. Default: NULL.
#' @param title The title of the plot. Default: NULL.
#' @param lineSize The size of the lines drawn. Default: 1.
#' @param reduction The dimension reduction method used ('UMAP', 'tSNE', or 'PCA'). Default: 'UMAP'.
#' @param titleFsize The font size for the title. Default: 20.
#' @param line.colors The colors to be used for lines, specified as a vector of colors. Default: gray scale colors.
#' @return A ggplot object.
#' @details This function creates a ggplot object that visualizes a Slingshot dataset, with curves representing the pseudotime trajectories.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  plotFittedGenePseudotime(data = sds, gene ="SST", expr = EXPR, loessCI = T, col = colz, pch = 20,
#'  panel_first = grid(NULL) )
#'  }
#' }
#' @seealso
#' \code{\link[slingshot]{SlingshotDataSet}}, \code{\link[ggplot2]{ggplot}}
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


# _________________________________________________________________________________________________

# plotting
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
      y <- exprs[gene, ,drop = FALSE][1,]
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
             main = paste(genename, ', Lineage ',l, sep=''), ...)
      }
      if (loess | loessCI){
        l <- loess(y ~ pst[,l], weights = w[,l])
      }
      if (loessCI){
        pl <- predict(l, se = TRUE, )
        polygon(c(l$x[order(l$x)],rev(l$x[order(l$x)])),
                c((pl$fit+qt(0.975,pl$df)*pl$se)[order(l$x)],
                  rev((pl$fit-qt(0.975,pl$df)*pl$se)[order(l$x)])),
                border = NA, col = rgb(0,0,0,.3))
      }
      if (loess){
        lines(l$x[order(l$x)], l$fitted[order(l$x)], lwd = 2, col = lcol[i])
      }
    }
    # par(mfrow = c(1,1))
    invisible(NULL)
  }
)




# _________________________________________________________________________________________________
#' @title cell_pal
#' @description Generate palette
#' @param cell_vars PARAM_DESCRIPTION
#' @param pal_fun PARAM_DESCRIPTION
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
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




