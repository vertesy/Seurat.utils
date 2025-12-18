######################################################################
# Cluster.Auto-naming.DE.R
######################################################################
# source('~/GitHub/Packages/Seurat.utils/Functions/Cluster.Auto-naming.DE.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Cluster.Auto-naming.DE.R"))

# ------------------------------------------------------------------------
# require(princurve) # only for AutoNumber.by.PrinCurve


# SmallestNonAboveX ------------------------------------------------------------------------
SmallestNonAboveX <- function(vec, X = 0) { # replace small values with the next smallest value found, which is >X.
  newmin <- min(vec[vec > X])
  vec[vec <= X] <- newmin
  vec
}
# SmallestNonZero(vec = df.markers$"p_val")


# AddCombinedScore2DEGResults ------------------------------------------------------------------------
AddCombinedScore2DEGResults <- function(df=df.markers, p_val_min=1e-25, pval_scaling = 0.001, colP = "p_val"
                                  , colLFC = grepv(pattern = c("avg_logFC|avg_log2FC"), x = colnames(df), perl = T)
                                  # , colLFC = "avg_log2FC"
                                  ) { # Score = -LOG10(p_val) * avg_log2FC
  p_cutoff <- SmallestNonAboveX(vec = df[[colP]], X = p_val_min)
  df$'combined.score' <- round(df[[colLFC]] * -log10( p_cutoff / pval_scaling ) )
  return(df)
}
# df.markers <- AddCombinedScore2DEGResults(df.markers)



# ------------------------------------------------------------------------------------
StoreTop25Markers <- function(obj = combined.obj # Save the top 25 markers based on `avg_log2FC` output table of `FindAllMarkers()` (df_markers) under `@misc$df.markers$res...`. By default, it rounds up insignificant digits up to 3.
                              , df_markers = df.markers, res = 0.5) {
  top25.markers <-
    df_markers %>%
    group_by(cluster) %>%
    top_n(n = 25, wt = avg_2logFC) %>%
    dplyr::select(gene) %>%
    col2named.vec.tbl() %>%
    splitbyitsnames()

  obj@misc$'top25.markers'[[ppp('res',res)]] <- top25.markers
  return(obj)
}
# combined.obj <- StoreTop25Markers(df_markers = df.markers, res = 0.)

# ------------------------------------------------------------------------------------
StoreAllMarkers <- function(obj = combined.obj # Save the output table of `FindAllMarkers()` (df_markers) under `@misc$df.markers$res...`. By default, it rounds up insignificant digits up to 3.
                            , df_markers = df.markers, res = 0.5, digit=c(0,3)[2]) {
  if (digit) df_markers[,1:5] <- signif(df_markers[,1:5], digits = digit)
  obj@misc$'df.markers'[[ppp('res',res)]] <- df_markers
  iprint("DF markers are stored under:", 'obj@misc$df.markers$', ppp('res',res))
  return(obj)
}
# combined.obj <- StoreAllMarkers(df_markers = df.markers, res = 0.5)


# GetTopMarkersDF ------------------------------------------------------------------------------------
GetTopMarkersDF <- function(dfDE = df.markers # Get the vector of N most diff. exp. genes.
                            , n = p$'n.markers', order.by = c("avg_log2FC", "p_val_adj")[1]) {
  'Works on active Idents() -> thus we call cluster'
  TopMarkers <- dfDE %>%
    arrange(desc(!!as.name(order.by))) %>%
    group_by(cluster) %>%
    dplyr::slice(1:n) %>%
    dplyr::select(gene)

  return(TopMarkers)
}
# GetTopMarkers(df = df.markers, n=3 )

# GetTopMarkers ------------------------------------------------------------------------------------
GetTopMarkers <- function(dfDE = df.markers # Get the vector of N most diff. exp. genes.
                            , n = p$'n.markers', order.by = c("combined.score", "avg_log2FC", "p_val_adj")[2]) {
  'Works on active Idents() -> thus we call cluster'
  TopMarkers <- dfDE %>%
    arrange(desc(!!as.name(order.by))) %>%
    group_by(cluster) %>%
    dplyr::slice(1:n) %>%
    dplyr::select(gene) %>%
    col2named.vec.tbl()

  return(TopMarkers)
}
# GetTopMarkers(df = df.markers, n=3 )



# ------------------------------------------------------------------------------------
AutoLabelTop.logFC <- function(obj = combined.obj # Create a new "named identity" column in the metadata of a Seurat object, with `Ident` set to a clustering output matching the `res` parameter of the function. It requires the output table of `FindAllMarkers()`. If you used `StoreAllMarkers()` is stored under `@misc$df.markers$res...`, which location is assumed by default.
                               , res = 0.2, plot.top.genes = T
                               , order_by = c("combined.score", "avg_logFC", "p_val_adj")[1]
                               , df_markers = combined.obj@misc$"df.markers"[[paste0("res.",res)]] ) {
  stopifnot(!is.null("df_markers"))
  stopifnot(order_by %in% colnames(df_markers))

  top.markers <-
    GetTopMarkersDF(df = df_markers, order.by = order_by, n = 1) %>%
    col2named.vec.tbl()

  obj@misc[[ppp("top.markers.res",res)]] <- top.markers

  ids <- unique(Idents(object = obj))
  if(length(ids) != length(top.markers)) {
    warning("Not all clusters returned DE-genes!")
    missing <- setdiff(ids, names(top.markers));  names(missing) <- missing
    iprint("missing:", missing)
    top.markers <- sortbyitsnames(c(top.markers, missing))
  }

  (top.markers.ID <- ppp(names(top.markers), top.markers))
  names(top.markers.ID) <- names(top.markers)
  named.ident <- top.markers.ID[Idents(object = obj)]

  namedIDslot <- ppp('cl.names.top.gene.res', res )
  obj[[namedIDslot]] = named.ident

  if (plot.top.genes) multiFeaturePlot.A4(list.of.genes = top.markers)

  return(obj)
}

# combined.obj <- AutoLabelTop.logFC(); combined.obj$"cl.names.top.gene.res.0.5"



# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------
AutoLabel.KnownMarkers <- function(obj = combined.obj, topN =1, res = 0.5 # Create a new "named identity" column in the metadata of a Seurat object, with `Ident` set to a clustering output matching the `res` parameter of the function. It requires the output table of `FindAllMarkers()`. If you used `StoreAllMarkers()` is stored under `@misc$df.markers$res...`, which location is assumed by default.
                                   , KnownMarkers=c("TOP2A", "EOMES", "SLA", "HOPX", "S100B", "DLX6-AS1", "POU5F1","SALL4","DDIT4", "PDK1", "SATB2", "FEZF2")
                                   , order.by = c("combined.score", "avg_log2FC", "p_val_adj")[1]

                                   , df_markers = combined.obj@misc$"df.markers"[[paste0("res.",res)]] ) {
  stopifnot(!is.null("df_markers"))

  lfcCOL <- grepv(pattern = c("avg_logFC|avg_log2FC"), x = colnames(df_markers), perl = T)
  keep <- unique(c(lfcCOL, 'p_val_adj', 'cluster', order.by, 'gene'  ))


  matching.clusters <-
    df_markers %>%
    dplyr::select(keep) %>%
    arrange(desc(!!as.name(order.by))) %>%
    filter(gene %in%  KnownMarkers) %>%
    group_by(gene) %>%
    dplyr::slice(1:topN) %>%
    arrange(desc(!!as.name(order.by))) %>%
    # top_n(n=1, wt=avg_log2FC) %>% # Select the top cluster for each gene
    arrange(cluster)

  print(matching.clusters)

  unique.matches <-
    matching.clusters %>%
    group_by(cluster) %>% # Select rows with unique values based on column "cluster"
    distinct(cluster,.keep_all = T)  %>%
    dplyr::select(gene)

  print("Best matches:")
  print(unique.matches)

  top.markers.df <- GetTopMarkersDF(dfDE = df_markers, order.by = lfcCOL, n = 1)
  top.markers <- top.markers.df %>% col2named.vec.tbl()

  missing.annotations <-
    top.markers.df %>%
    filter(!cluster %in%  unique.matches$cluster) # filter for clusters that do not have a unique label already

  named.annotations <-
    rbind(unique.matches, missing.annotations) %>%  # merge the 2 df's
    arrange(cluster) %>%
    col2named.vec.tbl() # requires github.com/vertesy/CodeAndRoll

  (top.markers.ID <- ppp(names(named.annotations), named.annotations))
  names(top.markers.ID) <- names(top.markers)
  named.ident <- top.markers.ID[Idents(object = obj)]

  namedIDslot <- ppp('cl.names.KnownMarkers', res )
  obj[[namedIDslot]] = named.ident
  return(obj)
}
# combined.obj <- AutoLabel.KnownMarkers(); # combined.obj$cl.names.KnownMarkers.0.5
# DimPlot.ClusterNames(ident = "cl.names.KnownMarkers.0.5")
# qUMAP("XACT")


# ------------------------------------------------------------------------------------
DimPlot.ClusterNames <- function(obj = combined.obj # Plot UMAP with Cluster names.
                                 , ident = "cl.names.top.gene.res.0.5", reduct = "umap", title = ident, ...) {
  DimPlot(object = obj, reduction = reduct, group.by = ident, label = T, repel = T, ...) + NoLegend() + ggtitle(title)
}
# DimPlot.ClusterNames()


# ------------------------------------------------------------------------------------
AutoNumber.by.UMAP <- function(obj = combined.obj # Relabel cluster numbers along a UMAP (or tSNE) axis
                               , dimension=1, swap= F, reduction="umap", res = "integrated_snn_res.0.5" ) {

  dim_name <- kppu(toupper(reduction),dimension)
  coord.umap <- as.named.vector(FetchData(object = obj, vars = dim_name))
  ls.perCl <- split(coord.umap, f = obj[[res]])
  MedianClusterCoordinate <- unlapply(ls.perCl, median)
  OldLabel <- names(sort(MedianClusterCoordinate, decreasing = swap))
  NewLabel <- as.character(0:(length(MedianClusterCoordinate) - 1))
  NewMeta <- translate(vec = obj[[res]], old = OldLabel, new = NewLabel)
  NewMetaCol <- kpp(res,"ordered")
  iprint("NewMetaCol:",NewMetaCol)
  obj[[NewMetaCol]] <- NewMeta
  return(obj)
}

# combined.obj <- AutoNumber.by.UMAP(obj = combined.obj, dimension=1, reduction="umap", res = "integrated_snn_res.0.5" )
# DimPlot.ClusterNames(ident = "integrated_snn_res.0.5.ordered")

# ------------------------------------------------------------------------------------
AutoNumber.by.PrinCurve <- function(obj = combined.obj # Relabel cluster numbers along the principal curve of 2 UMAP (or tSNE) dimensions.
                                    , dimension=1:2, plotit=T, swap= -1
                                    , reduction="umap", res = "integrated_snn_res.0.5" ) {
  # require(princurve)
  dim_name <- ppu(toupper(reduction),dimension)
  coord.umap <- FetchData(object = obj, vars = dim_name)
  fit <- princurve::principal_curve(x = as.matrix(coord.umap))
  if (plotit) {
    plot(fit, xlim = range(coord.umap[, 1]), ylim = range(coord.umap[, 2])
         , main = "principal_curve")
    # points(fit)
    points(coord.umap, pch = 18, cex = .25)
    whiskers(coord.umap, fit$s, lwd = .1)
    wplot_save_this(plotname = "principal_curve")
  }

  ls.perCl <- split(swap * fit$lambda, f = obj[[res]])
  MedianClusterCoordinate <- unlapply(ls.perCl, median)
  OldLabel <- names(sort(MedianClusterCoordinate))
  NewLabel <- as.character(0:(length(MedianClusterCoordinate) - 1))
  NewMeta <- translate(vec = obj[[res]], old = OldLabel, new = NewLabel)
  NewMetaCol <- kpp(res,"prin.curve")
  iprint("NewMetaCol:",NewMetaCol)
  obj[[NewMetaCol]] <- NewMeta
  return(obj)
}
# DimPlot.ClusterNames(ident = "integrated_snn_res.0.5")
# combined.obj <- AutoNumber.by.PrinCurve(obj = combined.obj, dimension=1:2, reduction="umap", plotit=T, swap= -1, res = "integrated_snn_res.0.5" )
# DimPlot.ClusterNames(ident = "integrated_snn_res.0.5.prin.curve")

######################################################################
# Custom.Functions.for.Slingshot.R
######################################################################
# source('~/GitHub/Packages/Seurat.utils/Functions/Custom.Functions.for.Slingshot.R')
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


######################################################################
# Jaccard.toolkit.R
######################################################################
# try(source('~/GitHub/Packages/Seurat.utils/Functions/Jaccard.toolkit.R'))
# try(source('https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Jaccard.toolkit.R'))


# Functions ------------------------
# try(source('~/GitHub/Packages/CodeAndRoll/CodeAndRoll.R'),silent= F)
# try(require('MarkdownReportsDev'),  silent = T)
# try(require('tidyverse'),  silent = T)
# source('~/Github/TheCorvinas/R/DatabaseLinke.r')


# --------------------------------------------------------------------------------------------
# Fast direct calculation from a list --------------------------------------------------------
# --------------------------------------------------------------------------------------------


# jJaccardIndexVec ----------------------------------------
jJaccardIndexVec <- function(A = 1:3, B = 2:4) length(intersect(A,B)) / length(union(A,B))

# jPairwiseJaccardIndexList ----------------------------------------

jPairwiseJaccardIndexList <- function(lsG = ls_genes) { # Create a pairwise jaccard similarity matrix across all combinations of columns in binary.presence.matrix. Modified from: https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/
  if (l(names(lsG)) < l(lsG)) {
    iprint("Gene lists were not (all) named, now renamed as:")
    names(lsG) <- ppp("dataset", 1:l(lsG))
    print(names(lsG))
  }
  m = matrix.fromNames(rowname_vec = names(lsG), colname_vec = names(lsG))
  n.sets <- length(lsG)
  for (r in 1:n.sets) {
    # print(percentage_formatter(r/n.sets))
    for (c in 1:n.sets) {
      if (c == r) {
        m[r,c] = 1
      } else {
        m[r,c] =signif(jJaccardIndexVec(lsG[[r]], lsG[[c]]), digits = 2)
      }
    }
  }
  return(m)
}
# jPairwiseJaccardIndexList(lsG = ls_genes)






# --------------------------------------------------------------------------------------------
# Much slower Indirect calculation via PresenceMatrix ----------------------------------------
# --------------------------------------------------------------------------------------------


# jPresenceMatrix ----------------------------------------
jPresenceMatrix <- function(string_list = lst(a=1:3, b=2:5,c=4:9, d=-1:4) ) { # Make a binary presence matrix from a list. Source: https://stackoverflow.com/questions/56155707/r-how-to-create-a-binary-relation-matrix-from-a-list-of-strings
  df.presence <- string_list %>%
    enframe %>%
    unnest(cols = "value") %>%
    count(name, value) %>%
    spread(value, n, fill = 0)
  df.presence2 <- FirstCol2RowNames(df.presence)
  return(t(df.presence2))
}
# df.presence <- jPresenceMatrix(string_list = lst(a=1:3, b=2:5,c=4:9, d=-1:4))

# jJaccardIndexBinary ----------------------------------------
jJaccardIndexBinary = function (x, y) { # Calculate Jaccard Index. Modified from: https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/
  elements.found <- sort(unique(union(x, y)))
  stopifnot(l(elements.found) == 2) # check if you only have [0,1]
  stopifnot(as.numeric(elements.found) == 0:1) # check if you only have [0,1]

  M.11 = sum(x == 1 & y == 1)
  M.10 = sum(x == 1 & y == 0)
  M.01 = sum(x == 0 & y == 1)
  return (M.11 / (M.11 + M.10 + M.01))
}
# JaccardSimilarity <- jJaccardIndexBinary(  x=sample(x = 0:1, size = 100, replace = T)
#               , y=sample(x = 0:1, size = 100, replace = T))



# jPairwiseJaccardIndex ----------------------------------------
jPairwiseJaccardIndex <- function(binary.presence.matrix = df.presence) { # Create a pairwise jaccard similarity matrix across all combinations of columns in binary.presence.matrix. Modified from: https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/
  m = matrix.fromNames(rowname_vec = colnames(binary.presence.matrix), colname_vec = colnames(binary.presence.matrix) )
  n.sets <- ncol(binary.presence.matrix)
  for (r in 1:n.sets) {
    print(percentage_formatter(r/n.sets))
    for (c in 1:n.sets) {
      if (c == r) {
        m[r,c] = 1
      } else {
        m[r,c] = signif(jJaccardIndexBinary(binary.presence.matrix[,r], binary.presence.matrix[,c]), digits = 2)
      }
    }
  }
  return(m)
}
# PairwiseJaccardIndices <- jPairwiseJaccardIndex(binary.presence.matrix = df.presence)

######################################################################
# metadata.manipulation.R
######################################################################
# source('~/GitHub/Packages/Seurat.utils/Functions/Seurat.object.manipulations.etc.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Metadata.manipulation.R"))
# Source: self + web

# - getMedianMetric.lsObj
# - add.meta.tags
# - addMetaFraction
# - GetClusteringRuns
# - GetNamedClusteringRuns
# - GetOrderedClusteringRuns
# - GetNumberOfClusters
# - getMetadataColumn <- mmeta
# - getCellIDs.from.meta
# - seu.add.meta.from.vector
# - seu.map.and.add.new.ident.to.meta
# - calc.cluster.averages
# - seu.add.meta.from.table
# - sampleNpc
# - calc.q90.Expression.and.set.all.genes
# - PlotTopGenes
# - fix.orig.ident
# - set.all.genes
# - recallAllGenes
# - recall.meta.tags.n.datasets
# - recall.parameters
# - recall.genes.l
# - save.parameters
# - plot.expression.rank.q90
# - FlipReductionCoordinates
# - SeuratColorVector
# - getClusterColors
# - get.clustercomposition

# getMedianMetric.lsObj ------------------------------------------------------------------------------------------------
getMedianMetric.lsObj <- function(ls.Obj = ls.Seurat, n.datasets = length(ls.Seurat), mColname = "percent.mito") {
  medMetric <- vec.fromNames(names(ls.Seurat))
  for(i in 1:n.datasets ) {
    medMetric[i] <- median(ls.Seurat[[i]]@meta.data[,mColname])
  }
  return(medMetric)
}
# ls.Seurat <- getMedianMetric.lsObj(ls.Obj = ls.Seurat, n.datasets = length(ls.Seurat), mColname = "percent.mito")


# add.meta.tags ------------------------------------------------------------------------------------------------
add.meta.tags <- function(list.of.tags = tags, obj = ls.Seurat[[1]], n = 1) {  # N is the for which dataset
  stopifnot( length(names(tags)) == length(tags) )
  nCells = nrow(obj@meta.data)
  for (i in 1:length(list.of.tags)) {
    tagX <- list.of.tags[[i]]
    new.meta.tag.per.cell <- rep(tagX[n], nCells)
    obj <- AddMetaData(object = obj, metadata = new.meta.tag.per.cell, col.name = names(tags)[i])
  }
  return(obj)
}
# ls.Seurat[[1]] <- add.meta.tags(list.of.tags = tags, obj = ls.Seurat[[1]], n = 1)


# addMetaFraction ------------------------------------------------------------------------------------------------
addMetaFraction <- function(col.name = "percent.mito", gene.symbol.pattern = c("^MT\\.|^MT-", F)[1]
                              , gene.set = F, obj = ls.Seurat[[1]], verbose = T) {
  stopif(condition = isFALSE(gene.set) && isFALSE(gene.symbol.pattern), "Either gene.set OR gene.symbol.pattern has to be defined (!= FALSE).")
  if (!isFALSE(gene.set) && !isFALSE(gene.symbol.pattern) && verbose) print("Both gene.set AND gene.symbol.pattern are defined. Only using gene.set.")

  geneset <- check.genes(list.of.genes = gene.set, obj = obj)
  total_expr <- Matrix::colSums(GetAssayData(object = obj))
  genes.matching <- if (!isFALSE(gene.set)) intersect(gene.set, rownames(obj)) else grepv(pattern = gene.symbol.pattern, x = rownames(obj))

  genes.expr = GetAssayData(object = obj)[genes.matching, ]
  target_expr <- if (l(genes.matching) >1) Matrix::colSums(genes.expr) else genes.expr
  obj <- AddMetaData(object = obj, metadata = target_expr / total_expr, col.name = col.name)
  colnames(obj@meta.data)
  return(obj)
}


# ls.Seurat[[1]] <- addMetaFraction(col.name = "percent.mito", gene.symbol.pattern = "^MT\\.|^MT-")
# ls.Seurat[[1]] <- addMetaFraction(col.name = "percent.ribo", gene.symbol.pattern = "^RPL|^RPS")
# ls.Seurat[[1]] <- addMetaFraction(col.name = "percent.AC.GenBank", gene.symbol.pattern = "^AC[0-9]{6}\\.")
# ls.Seurat[[1]] <- addMetaFraction(col.name = "percent.AL.EMBL", gene.symbol.pattern = "^AL[0-9]{6}\\.")
# ls.Seurat[[1]] <- addMetaFraction(col.name = "percent.LINC", gene.symbol.pattern = "^LINC0")
# ls.Seurat[[1]] <- addMetaFraction(col.name = "percent.MALAT1", gene.symbol.pattern = "^MALAT1")
# colnames(ls.Seurat[[1]]@meta.data)
# HGA_MarkerGenes <- c("ENO1", "IGFBP2", "WSB1", "DDIT4", "PGK1", "BNIP3", "FAM162A", "TPI1", "VEGFA", "PDK1", "PGAM1", "IER2", "FOS", "BTG1", "EPB41L4A-AS1","NPAS4", "HK2", "BNIP3L", "JUN", "ENO2", "GAPDH", "ANKRD37", "ALDOA", "GADD45G", "TXNIP")
# sobj <- addMetaFraction(col.name = "percent.HGA", gene.set = HGA_MarkerGenes, obj =  sobj)


# ------------------------------------------------------------------------------------
GetClusteringRuns <- function(obj = combined.obj, res = F, pat = "*snn_res.*[0-9]$") { # Get Clustering Runs: metadata column names
  if (res) pat = gsub(x = pat, pattern = '\\[.*\\]', replacement = res)
  clustering.results <- grepv(x = colnames(obj@meta.data), pattern = pat)
  if ( identical(clustering.results, character(0)) ) warning("No matching column found!")
  return(clustering.results)
}
# GetClusteringRuns()


# ------------------------------------------------------------------------------------
GetNamedClusteringRuns <- function(obj = combined.obj  # Get Clustering Runs: metadata column names
                                   , res = c(F, 0.5)[1], topgene = F, pat = "^cl.names.Known.*[0,1]\\.[0-9]$") {
  if (res) pat = gsub(x = pat, pattern = '\\[.*\\]', replacement = res)
  if (topgene) pat = gsub(x = pat, pattern = 'Known', replacement = 'top')
  clustering.results <- grepv(x = colnames(obj@meta.data), pattern = pat)
  if ( identical(clustering.results, character(0)) ) {
    print("Warning: NO matching column found! Trying GetClusteringRuns(..., pat = '*_res.*[0,1]\\.[0-9]$)")
    clustering.results <- GetClusteringRuns(obj = obj, res = F, pat = "*_res.*[0,1]\\.[0-9]$")
  }
  return(clustering.results)
}
# GetNamedClusteringRuns()


# ------------------------------------------------------------------------------------
GetOrderedClusteringRuns <- function(obj = combined.obj, res = F, pat = "*snn_res.*[0,1]\\.[0-9]\\.ordered$") { # Get Clustering Runs: metadata column names
  if (res) pat = gsub(x = pat, pattern = '\\[.*\\]', replacement = res)
  clustering.results <- grepv(x = colnames(obj@meta.data), pattern = pat)
  if ( identical(clustering.results, character(0)) ) warning("No matching column found!")
  return(clustering.results)
}
# GetOrderedClusteringRuns(); GetOrderedClusteringRuns(res = 0.5)


# ------------------------------------------------------------------------------------
GetNumberOfClusters <- function(obj = combined.obj) { # Get Number Of Clusters
  clustering.results <- GetClusteringRuns(obj)
  print("## Number of clusters: ---------")
  for (cc in clustering.results) {
    NrCl <- length(unique(obj@meta.data[[cc]]))
    iprint( cc, "   ", NrCl)
  }
}
# GetNumberOfClusters()

# get Cells from metadata  ------------------------------------------------
getMetadataColumn <- mmeta <- function(ColName.metadata = 'batch', obj = combined.obj, as_numeric =F) { # Get a metadata column from a Seurat object as a named vector
  stopifnot(ColName.metadata %in% colnames(obj@meta.data))

  x = as.named.vector(obj@meta.data[ ,ColName.metadata, drop=F])
  if (as_numeric) {
    as.numeric.wNames(x)+1
  } else {x}
}

# GetCellIDs from metadata ---------------
getCellIDs.from.meta <- function(ident = 'res.0.6', values = NA, obj=combined.obj, inverse = F ) { # Get cellIDs from a metadata column, matching a list of values (using %in%).
  mdat <- obj@meta.data[ , ident]
  cells <- if (inverse) {mdat %!in% values} else {mdat %in% values}
  idx.matching.cells = which(cells)
  iprint(length(idx.matching.cells), 'cells found.')
  return(rownames(obj@meta.data)[idx.matching.cells])
}
# getCellIDs.from.meta()

# seu.add.meta.from.vector ------------------------------------------------------------------------
seu.add.meta.from.vector <- function(obj = combined.obj, metaD.colname = metaD.colname.labeled, Label.per.cell=Cl.Label.per.cell ) { # Add a new metadata column to a Seurat  object
  obj@meta.data[, metaD.colname ] = Label.per.cell
  iprint(metaD.colname, "contains the named identities. Use Idents(combined.obj) = '...'. The names are:", unique(Label.per.cell))
  return(obj)
}
# combined.obj <- add.Cl.Label.2.Metadata(obj = combined.obj, metaD.colname = metaD.colname.labeled, Label.per.cell=Cl.Label.per.cell )
# formerly add.Cl.Label.2.Metadata


# seu.map.and.add.new.ident.to.meta ------------------------------------------------------------------------

seu.map.and.add.new.ident.to.meta <- function(obj = combined.obj, ident.table = clusterIDs.GO.process
                                              , metaD.colname = substitute(ident.table) ) { # Add a new metadata column to a Seurat  object
  ident.vec <- as.named.vector(ident.table)

  # identities should match ----------------
  ident.X <- names(ident.vec)
  ident.Y <- as.character(ident.vec)
  ident.Seu <- gtools::mixedsort(levels(Idents(obj)))
  iprint("ident.Seu: ", ident.Seu)

  OnlyInIdentVec      <- setdiff(ident.X, ident.Seu)
  OnlyInSeuratIdents  <- setdiff(ident.Seu, ident.X)

  msg.IdentVec <- kollapse("Rownames of 'ident.table' have entries not found in 'Idents(obj)':"
                           , OnlyInIdentVec, " not found in ", ident.Seu, collapseby = " ")

  msg.Seu <- kollapse("Rownames of 'Idents(obj)' have entries not found in 'ident.table':"
                      , OnlyInSeuratIdents, " not found in ", ident.X, collapseby = " ")

  stopif (l(OnlyInIdentVec), message = msg.IdentVec)
  stopif (l(OnlyInSeuratIdents), message = msg.Seu)

  # identity mapping ----------------
  new.ident <- translate(vec = as.character(Idents(obj)), old = ident.X, new = ident.Y)
  obj@meta.data[[metaD.colname]] = new.ident
  iprint(metaD.colname, "contains the named identities. Use Idents(combined.obj) = '...'. The names are:"); cat(paste0("\t", ident.Y, "\n"))
}
# combined.obj <- seu.map.and.add.new.ident.to.meta(obj = combined.obj, ident.table = clusterIDs.GO.process)




# calc.cluster.averages ------------------------------------------------
calc.cluster.averages <- function(col_name = "Score.GO.0006096"
                                  , plot.UMAP.too = TRUE
                                  , return.plot = F
                                  , obj =  combined.obj
                                  , split_by = GetNamedClusteringRuns()[1]
                                  , scale.zscore = FALSE
                                  , simplify=T, plotit = T
                                  , histogram = FALSE, nbins = 50
                                  , suffix = NULL
                                  , stat = c("mean", "median")[2]
                                  , quantile.thr = 0.9
                                  , absolute.thr = FALSE
                                  , filter = c(FALSE, 'above', 'below')[1]
                                  , ylab.text = paste("Cluster", stat, "score")
                                  , title = paste("Cluster", stat, col_name)
                                  , subtitle = NULL
                                  , width = 8, height =6
                                  , ...
                                  # , ylb = paste(ylab.text, col_name)
                                  # , xlb = paste("Clusters >",percentage_formatter(quantile.thr),"quantile are highlighted. |", split_by)
                                  , xlb = if (absolute.thr) paste("Threshold at", absolute.thr) else paste(
                                    "Black lines: " , kppd(percentage_formatter(c(1-quantile.thr, quantile.thr))) ,"quantiles |"
                                    , "Cl. >",percentage_formatter(quantile.thr),"are highlighted. |", split_by
                                  )

                                  , fname = ppp(col_name,split_by,"cluster.average.barplot.pdf", ...)
) { # calc.cluster.averages of a m
  iprint(substitute(obj), "split by", split_by)
  if(absolute.thr) iprint('In case of the absolute threshold, only the returned values are correct, the plot annotations are not!')

  if (plot.UMAP.too) qUMAP(obj = obj, feature = col_name)

  df.summary <-
    obj@meta.data %>%
    select_at(c(col_name, split_by)) %>%
    group_by_at(split_by) %>%
    summarize('nr.cells' = n()
              , 'mean' = mean(!!sym(col_name), na.rm = TRUE)
              , 'SEM' = sem(!!sym(col_name), na.rm = TRUE)
              , 'median' = median(!!sym(col_name), na.rm = TRUE)
              , 'SE.median' = 1.2533 * sem(!!sym(col_name), na.rm = TRUE)
    )

  if (simplify) {
    av.score <- df.summary[[stat]]
    names(av.score) <- ppp("cl",df.summary[[1]])
    av.score <- sortbyitsnames(av.score)
    if (scale.zscore) av.score <- (scale(av.score)[,1])

    cutoff <- if(absolute.thr) absolute.thr else quantile(av.score, quantile.thr)
    cutoff.low <- if(absolute.thr) NULL else  quantile(av.score, (1-quantile.thr) )

    if (plotit) {
      if (histogram) {
        p <- qhistogram(vec = as.numeric(av.score), save = F
                        , vline = cutoff
                        , plotname = ppp(title, quantile.thr)
                        , bins = nbins
                        , subtitle = paste(subtitle, "| median in blue/dashed")
                        , ylab = ylab.text
                        , xlab = xlb # Abused
                        , xlab.angle = 45
                        # , ylim=c(-1,1)
                        , ...
                        # , ext = "png", w = 7, h = 5
        ) + geom_vline(xintercept = cutoff.low, lty=2)
        print(p)
        title_ <- ppp(title, suffix, flag.nameiftrue(scale.zscore))
        qqSave(ggobj = p, title = title_, ext = "png", w = width, h = height)
      } else {
        p <- qbarplot(vec = av.score, save = F
                      , hline = cutoff
                      , title = title
                      , suffix = quantile.thr
                      , subtitle = subtitle
                      , ylab = ylab.text
                      , xlab = xlb # Abused
                      , xlab.angle = 45
                      # , ylim=c(-1,1)
                      , ...
                      # , ext = "png", w = 7, h = 5
        ) + geom_hline(yintercept = cutoff.low , lty = 2)

        print(p)
        title_ <- ppp(title, suffix, flag.nameiftrue(scale.zscore))
        qqSave(ggobj = p, title = title_, fname = ppp(title_, split_by, "png"),  w = width, h = height)
      }
    }
    print(quantile.thr)

    if (filter == 'below') {
      return(filter_LP(av.score, threshold = cutoff, plot.hist = F))
    } else if (filter == 'above') {
      return(filter_HP(av.score, threshold = cutoff, plot.hist = F))
    } else {
      return(av.score)
    }
  } else if (return.plot) { # if /not/ simplify
    return(p)
  } else {
    return(df.summary)
  }
}

# Add to obj@metadata from an external table ------------------------------------------------------------------------
seu.add.meta.from.table <- function(obj = seu.ORC, meta = MetaData.ORC, suffix = ".fromMeta") { # Add multiple new metadata columns to a Seurat object from a table.
  NotFound  = setdiff(colnames(obj), rownames(meta))
  Found     = intersect(colnames(obj), rownames(meta))
  if (length(NotFound)) iprint(length(NotFound), 'cells were not found in meta, e.g.: ', trail(NotFound, N = 10))

  mCols.new = colnames(meta)
  mCols.old = colnames(obj@meta.data)
  overlap = intersect(mCols.new, mCols.old)
  if (length(overlap)) {
    iprint(length(overlap), 'metadata columns already exist in the seurat object: ', overlap, '. These are tagged as: *', suffix)
    colnames(meta)[overlap] = paste0(overlap, suffix)
  }
  mCols.add = colnames(meta)
  obj@meta.data[Found, mCols.add] = meta[ Found,]

  return(obj)
} # x=seu.add.meta.from.table()


# sampleNpc ------------------------------------------------------------------------
sampleNpc <- function(metaDF = MetaData[which(Pass),], pc=0.1) { # Sample N % of a dataframe (obj@metadata), and return the cell IDs.
  cellIDs = rownames(metaDF)
  nr_cells = floor(length(cellIDs) * pc)
  cellIDs.keep = sample(cellIDs, size = nr_cells, replace = FALSE)
  return(cellIDs.keep)
}


# calc.q90.Expression.and.set.all.genes ------------------------------------------------------------------------
calc.q90.Expression.and.set.all.genes <- function(obj = combined.obj # Calculate the gene expression of the e.g.: 90th quantile (expression in the top 10% cells).
                              , quantileX=0.9, max.cells =  100000
                              , slot = "data", assay = c("RNA", "integrated")[1]
                              , set.all.genes = TRUE, show = TRUE) {
  tic()
  x = GetAssayData(object = obj, assay = assay, slot = slot) #, assay = 'RNA'
  if (ncol(x) > max.cells) {
    dsampled = sample(x = 1:ncol(x), size = max.cells)
    x = x[ , dsampled]
  }
  qname = p0("q", quantileX * 100)
  slot_name = kpp("expr", qname)

  # expr.q90 = iround(apply(x, 1, quantile, probs = quantileX) )
  expr.q90.df = sparseMatrixStats::rowQuantiles(x, probs = quantileX)
  expr.q90 = iround(as.named.vector(expr.q90.df))
  toc();

  log2.gene.expr.of.the.90th.quantile <- as.numeric(log2(expr.q90 + 1)) # strip names
  try(
    qhistogram(log2.gene.expr.of.the.90th.quantile, ext = "pdf", breaks = 30
               , plotname = kpp("log2.gene.expr.of.the ", qname," quantile")
               , subtitle = kollapse(pc_TRUE(expr.q90 > 0, NumberAndPC = T), " genes have ", qname ," expr. > 0.")
               , xlab = p0("log2(expr.",qname,"+1) [UMI]"), ylab = "Genes"
               , plot = show, save = TRUE, vline  = .2)
    , silent = TRUE)

  {
    all.genes = percent_rank(expr.q90); names(all.genes) = names(expr.q90); all.genes <- sort.decreasing(all.genes)
    if (set.all.genes) obj@misc$'all.genes' = all.genes = as.list(all.genes)
    assign('all.genes', all.genes, envir = as.environment(1))
  }

  obj@misc[[slot_name]] <-  expr.q90

  iprint('Quantile', quantileX ,'is now stored under obj@misc$all.genes and $', slot_name, ' Please execute all.genes <- obj@misc$all.genes.')
  return(obj)
}

# combined.obj <- calc.q90.Expression.and.set.all.genes(obj = combined.obj, quantileX=0.9, max.cells =  25000)
# head(sort(as.numeric.wNames(obj@misc$expr.q90), decreasing = T))
# combined.obj <- calc.q90.Expression.and.set.all.genes(obj = combined.obj, quantileX=0.95, max.cells =  25000, set.all.genes = FALSE)

# PlotTopGenes ------------------------------------------------------------------------
PlotTopGenes <- function(obj = combined.obj, n=32 ){ # Plot the highest expressed genes on umaps, in a subfolder. Requires calling calc.q90.Expression.and.set.all.genes before.
  Highest.Expressed.Genes = names(head(sort(obj@misc$expr.q90, decreasing = T), n = n))
  multiFeaturePlot.A4(list.of.genes = Highest.Expressed.Genes, foldername = "Highest.Expressed.Genes" )
}
# PlotTopGenes()


# fix.orig.ident ------------------------------------------------------------------------
fix.orig.ident <- function(obj = merged.obj) {
  fixed <- sub(obj$'orig.ident', pattern = 'filtered_feature_bc_matrix.', replacement = '')
  return(fixed)
}
# merged.obj$orig.ident <- fix.orig.ident(obj = merged.obj); table(merged.obj$orig.ident)

# set.all.genes ------------------------------------------------------------------------
set.all.genes <- function(obj = combined.obj) iprint("Use calc.q90.Expression.and.set.all.genes()")
# set.all.genes(); all.genes

# recallAllGenes ------------------------------------------------------------------------
recallAllGenes <- function(obj = combined.obj) { # all.genes set by calc.q90.Expression.and.set.all.genes()
  if (!exists('all.genes')) {
    all.genes <- obj@misc$all.genes
    print(head(unlist(all.genes)))
    ww.assign_to_global(name = "all.genes", value = all.genes)
  } else {print("variable 'all.genes' exists in the global namespace")}
}
# recallAllGenes(); all.genes

# recall.meta.tags.n.datasets ------------------------------------------------------------------------
recall.meta.tags.n.datasets <- function(obj = combined.obj) {
  if (!exists('n.datasets')) {
    n.datasets <- obj@misc$n.datasets
    print(head(unlist(n.datasets)))
    ww.assign_to_global(name = "n.datasets", value = n.datasets)
  } else {print("variable 'n.datasets' exists in the global namespace")}

  if (!exists('meta.tags')) {
    meta.tags <- obj@misc$meta.tags
    print(head(unlist(meta.tags)))
    ww.assign_to_global(name = "meta.tags", value = meta.tags)
  } else {print("variable 'meta.tags' exists in the global namespace")}

}
# recall.n.datasets(); n.datasets

# recall.parameters ------------------------------------------------------------------------
recall.parameters <- function(obj = combined.obj, overwrite = FALSE) {
  if (is.null(obj@misc$'p')) {
    print("obj does not have: obj@misc$p")
  } else {
    p <- obj@misc$'p'
    print(head(p))
    if (exists('p')) iprint("variable 'p' exists in the global namespace:");

    if (!exists('p') | (exists('p') & overwrite == TRUE) ) {
      ww.assign_to_global(name = "p", value = p); print("Overwritten.")
    } else {
      print("Not overwritten.")
    }
  } # else if obj@misc$'p'
}
# recall.parameters(); p


# recall.genes.ls ------------------------------------------------------------------------
recall.genes.ls<- function(obj = combined.obj) { # genes.ls
  if (!exists('genes.ls')) {
    genes.ls <- obj@misc$genes.ls
    print(head(unlist(genes.ls)))
    ww.assign_to_global(name = "genes.ls", value = genes.ls)
  } else {print("variable 'genes.ls' exists in the global namespace")}
}
# recall.genes.ls(); genes.ls


# save.parameters ------------------------------------------------------------------------
save.parameters <- function(obj = combined.obj, params = p) {
  if (!is.null(obj@misc$'p')) print("Overwriting already existing obj@misc$p. Old version:") ; print(head(unlist(obj@misc$'p')))
  obj@misc$p <- params
}
# save.parameters(obj = combined.obj, params = p);


# plot.expression.rank.q90 ------------------------------------------------------------------------
plot.expression.rank.q90 <- function(obj = combined.obj, gene="ACTB", filterZero=T) {
  expr.GOI <- obj@misc$expr.q90[gene]
  expr.all <- unlist(obj@misc$expr.q90)
  gene.found <- gene %in% names(expr.all)
  stopifnot(gene.found)

  if (expr.GOI == 0) iprint(gene, "is not expressed. q90-av.exp:", expr.GOI) else
    if (expr.GOI < 0.05) iprint(gene, "is lowly expressed. q90-av.exp:", expr.GOI)
    if (filterZero) {
      iprint("Zero 'q90 expression' genes (", pc_TRUE(expr.all == 0), ") are removed.")
      expr.all <- expr.all[expr.all > 0]
    }
  counts <- sum(obj@assays$RNA@counts[gene,])
  if (expr.GOI == 0) {
    quantile.GOI <- 0
    title <- paste(gene, "is too lowly expressed: q90-av.exp is zero. \n There are", counts,"counts." )
  } else {
    pos.GOI <- which(names(expr.all) == gene)
    quantile.GOI <- ecdf(expr.all)(expr.all)[pos.GOI]
    title <- paste(gene, "is in the", percentage_formatter(quantile.GOI), "quantile of 'q90-av' expression. \n There are", counts,"counts" )
  }
  suppressWarnings(
    whist(expr.all, vline = expr.GOI, breaks = 100, main = title, plotname =   make.names(title)
          , ylab = "Genes"
          , xlab = "Av. mRNA in the 10% top expressing cells (q90 av.exp.)")
  )
}
# plot.expression.rank.q90(gene = "SATB2")



# FlipReductionCoordinates ------------------------------------------------------------------------
FlipReductionCoordinates <- function(obj = combined.obj, dim=2, reduction="umap"
                                     , flip=c('x', 'y', 'xy', NULL)[1], FlipReductionBackupToo = TRUE) { # Set active UMAP to `obj@reductions$umap` from `obj@misc$reductions.backup`.
  coordinates <- Embeddings(obj, reduction = reduction)
  stopifnot(ncol(coordinates) == dim )

  if (flip %in% c('x', 'xy')) coordinates[,1] = coordinates[,1] * -1
  if (flip %in% c('y', 'xy')) coordinates[,2] = coordinates[,2] * -1
  obj@reductions[[reduction]]@cell.embeddings <- coordinates

  if (FlipReductionBackupToo) {
    bac.slot <- p0(reduction,dim,"d")
    if (length(obj@misc$reductions.backup[[bac.slot]])) {
      obj@misc$reductions.backup[[bac.slot]]@cell.embeddings <- coordinates
      iprint(dim, "dimensional",reduction,"backup flipped too.")
    }
  }
  return(obj)
}
# clUMAP(); combined.obj <- FlipReductionCoordinates(combined.obj); clUMAP()



# SeuratColorVector ------------------------------------------------------------------------
SeuratColorVector <- function(ident = NULL, obj = combined.obj, plot.colors = F) {
  if (!is.null(ident)) {
    print(ident)
    ident.vec <- obj[[ident]][,1]
  } else {
    ident.vec <- obj@active.ident
  }
  ident.vec <- as.factor(ident.vec)
  print(table(ident.vec))
  colorlevels <- scales::hue_pal()(length(levels(ident.vec)))
  if (plot.colors) color_check(colorlevels)
  translate(vec = as.character(ident.vec)
            , old = levels(ident.vec)
            , new = colorlevels)
}
# SeuratColorVector()
# SeuratColorVector(ident = GetNamedClusteringRuns()[2], plot.colors = T)

# getClusterColors ------------------------------------------------------------------------
getClusterColors <- function(obj = combined.obj
                             , ident = GetClusteringRuns()[1]
                             , show = T) {
  (identities <- levels(as.factor(obj[[ident]][,1])))
  color_palette <- scales::hue_pal()(length(identities))
  # color_check(color_palette)
  # names(color_palette) <- sort(as.factor(identities))
  names(color_palette) <- ((identities))
  identvec <- obj[[ident]][,1]
  colz <- color_palette[identvec]
  names(colz) <- identvec
  if (show) color_check(unique(colz)) # MarkdownReports
  colz
}
# getClusterColors(obj = combined.obj, ident = GetClusteringRuns()[2] )


#  get.clustercomposition ------------------------------------------------------------------------
get.clustercomposition <- function(obj = combined.obj, x = 'integrated_snn_res.0.3', y = 'project', color = y, ...) {
  setwd(OutDir)
  clUMAP(obj = obj, ident = x, save.plot = T, suffix = "as.in.barplot")
  categ.per.cluster <- ggbarplot(obj@meta.data
                                 , x = x
                                 , y = y
                                 , color = y
                                 , ...
  )
  qqSave(categ.per.cluster)
}

# get.clustercomposition()
# get.clustercomposition(, ylim=c(0,10000))


#  ------------------------------------------------------------------------
#  ------------------------------------------------------------------------

######################################################################
# Monocle.Utils.R
######################################################################
# source('~/GitHub/Packages/Seurat.utils/Functions/Monocle.Utils.R')
# rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)

# Functions ------------------------
# source('https://raw.githubusercontent.com/vertesy/Seurat.Pipeline/main/elements/Load.packages.local.R')
# try(source("~/GitHub/Packages/Seurat.multicore/00.Load.Seurat3.Multicore.LOCAL.R"));


# Metadata ------------------------

# Parameters ------------------------


# ------------------------
# ------------------------
# ------------------------
# ------------------------

# mplotGene ------------------------
mplotGene <- function(gene = "PGK1", reduc = "UMAP", obj = cds_from_seurat) {
  pl1 <- plot_cells(cds = obj,
                    # color_cells_by = 'clusters_low',
                    genes = gene,
                    reduction_method =reduc,
                    # color_cells_by = 'clusters_superlow',
                    label_groups_by_cluster = F,
                    label_leaves = F,
                    label_branch_points = F, label_cell_groups = F
                    , group_label_size = 10
                    , cell_size = 1, alpha = .5)
  qqSave(pl1, fname = ppp(reduc, gene, ".png"), w = 14, h=7)
}
# mplotGene()
# mplotGene("GAPDH")


# mplotManyGenes ------------------------
mplotManyGenes <- function(ls.genes = c(
  `S-phase` = "TOP2A", `G2M-phase` = "HIST1H4C"
  , `oRG` = "ID4", `oRG` = "HOPX"
  , `Intermediate progenitor` = "EOMES",  `Intermediate progenitor1` = "TAC3"
  , Astroglia = "GFAP", Astrocyte = "S100B"
  , `Immature neurons` = "SLA", Interneurons = "DLX6-AS1"
  , `Hypoxia/Stress` = "DDIT4", Glycolytic = "PDK1"
  , `Low-Quality` = "POLR2A", `PGC` = "DCN"
  , `dl-EN` = "KAZN", `ul-EN` = "SATB2")
) {
  create_set_SubDir("gene.expression")
  for (g in ls.genes) {
    mplotGene(g)
  }
  create_set_OutDir(ParentDir)
}
# mplotManyGenes()


m3DplotGene <- function(gene = "PGK1", reduc = "UMAP", obj = cds.10pc, ttl.suffix = "expression", suffix = "", cex = 5) {
  gene <- intersect(rownames(obj), gene)
  if (l(gene)) {
    pl1 <- plot_cells_3d(cds = obj
                         # color_cells_by = 'clusters_low',
                         , genes = gene
                         , reduction_method = reduc
                         # color_cells_by = 'clusters_superlow',
                         # label_groups_by_cluster = F,
                         # label_leaves = F,
                         # label_branch_points = F, label_cell_groups = F
                         # , group_label_size = 10
                         , cell_size = cex
                         , alpha = .5) %>% layout(title=p0(gene, ttl.suffix))
    SavePlotlyAsHtml(pl1, category. = gene, suffix. = suffix)
  } else { iprint("gene not found") }
}
# m3DplotGene("DDIT4")


# m3DplotKeyGenes ------------------------
m3DplotKeyGenes <- function(obj = cds.10pc, cex = iround(log10(idim(obj)[2]))
                            , ls.genes =  c(
                              `S-phase` = "TOP2A", `G2M-phase` = "HIST1H4C"
                              , `oRG` = "ID4", `oRG` = "HOPX"
                              , `Intermediate progenitor` = "EOMES",  `Intermediate progenitor1` = "TAC3"
                              , Astroglia = "GFAP", Astrocyte = "S100B"
                              , `Immature neurons` = "SLA", Interneurons = "DLX6-AS1"
                              , `Hypoxia/Stress` = "DDIT4", Glycolytic = "PDK1"
                              , `Low-Quality` = "POLR2A", `PGC` = "DCN"
                              , `dl-EN` = "KAZN", `ul-EN` = "SATB2")
                            , reduc = "UMAP", suffix = "") {

  ls.genes <- intersect(rownames(obj), ls.genes)
  iprint(l(ls.genes), "genes found.")
  create_set_SubDir(ppp("3D.gex.plots", substitute(obj)))
  for (g in ls.genes) {
    m3DplotGene(gene = g, obj = obj, cex = cex)
  }
  create_set_OutDir(ParentDir)
}
# m3DplotKeyGenes()

# subsetMonocleObject ------------------------
subsetMonocleObject <-function(obj = cds_from_seurat, fraction_ = 0.1, nCells = F, seed_ = 1989 ) { # Subset a compressed Seurat Obj and save it in wd.
  set.seed(seed_)
  if (isFALSE(nCells)) {
    all.cells <- colnames(obj)
    n.keep <- floor(l(all.cells) * fraction_)
    cellIDs.keep <- sample(all.cells, size = n.keep, replace = FALSE)
    iprint(length(cellIDs.keep), "or",percentage_formatter(fraction_),"of the cells are kept. Seed:", head(cellIDs.keep), seed_)
    # cellIDs.keep
  } else if (nCells > 1) {
    nKeep = min(ncol(obj), nCells)
    # print(nKeep)
    cellIDs.keep = sample(colnames(obj), size = nKeep, replace = F)
    if (nKeep < nCells) iprint("Only",nCells,"cells were found in the object, so downsampling is not possible.")
  }
  obj <- obj[,cellIDs.keep] # downsample
  return(obj)
}
# cds.10pc <- subsetMonocleObject(cds_from_seurat, fraction_ = .1)
# cds.25pc <- subsetMonocleObject(cds_from_seurat, fraction_ = .25)



# m3.get.umap ------------------------
m3.get.umap <- function(obj = cds_from_seurat, slot = 'UMAP', dim = (2:3)[2]) {
  iprint(dim, 'dimensional', slot)
  df.reduc <- obj@int_colData@listData$reducedDims[[slot]]
  stopif(condition = is_null(df.reduc), message = paste('slot not found', slot))
  stopifnot(ncol(df.reduc) == dim)
  return(df.reduc)
}
# m3.get.umap(obj = cds_from_seurat, slot = 'UMAP', dim = (2:3)[2])

# ------------------------
m3.backup.umap <- function(obj = cds_from_seurat, slot = 'UMAP', dim = (2:3)[2], ...) {
  new.slot <- p0(slot,'.',dim, 'D')
  obj@int_colData@listData$reducedDims[[new.slot]] <- m3.get.umap(obj = obj, slot = slot, dim = dim, ... )
  iprint('obj@int_colData@listData$reducedDims$', new.slot)
  return(obj)
}
# cds_from_seurat <- m3.backup.umap(obj = cds_from_seurat, slot = 'UMAP', dim = (2:3)[2])


# ------------------------
m3.recall.umap <- function(obj = cds_from_seurat, slot = 'UMAP', dim = (2:3)[2], ...) {
  backup.slot <- p0(slot,'.',dim, 'D')
  old.dim <- ncol(obj@int_colData@listData$reducedDims[[slot]])
  obj@int_colData@listData$reducedDims[[slot]] <- obj@int_colData@listData$reducedDims[[backup.slot]]
  iprint(old.dim, 'dimensional', slot, 'replaced by', dim, slot, 'reduction.')
  return(obj)
}
# cds_from_seurat <- m3.recall.umap(obj = cds_from_seurat, slot = 'UMAP', dim = (2:3)[2])
# cds_from_seurat.bac <- cds_from_seurat


# ------------------------
m3.export.umap.2.Seurat <- function(mobj = cds_from_seurat, sobj = combined.obj, def.dim = 2) {

  reduc_bac <- sobj@misc$reductions.backup
  reduc_bac$'umap2d.Seurat' <- reduc_bac$'umap2d'
  reduc_bac$'umap2d'@cell.embeddings <- m3.get.umap(obj = cds_from_seurat, slot = 'UMAP.2D', dim = 2)
  colnames(reduc_bac$'umap2d'@cell.embeddings) <- paste('UMAP', 1:2, sep = '_')

  reduc_bac$'umap3d.Seurat' <- reduc_bac$'umap3d'
  reduc_bac$'umap3d'@cell.embeddings <- m3.get.umap(obj = cds_from_seurat, slot = 'UMAP.3D', dim = 3)
  colnames(reduc_bac$'umap3d'@cell.embeddings) <- paste('UMAP', 1:3, sep = '_')

  sobj@misc$reductions.backup <- reduc_bac
  sobj@reductions$umap <-  reduc_bac$'umap2d'  # if (def.dim == 2) { reduc_bac$'umap2d' } else if (def.dim == 3) { reduc_bac$'umap3d' }

  return(sobj)
}
# combined.obj.bac <- combined.obj
# combined.obj <- combined.obj.bac
# combined.obj <- m3.export.umap.2.Seurat(mobj = cds_from_seurat, sobj = combined.obj); qUMAP("DDIT4")

# ------------------------






######################################################################
# MULTI-seq.functions.R
######################################################################
# source('~/GitHub/Packages/Seurat.utils/Functions/MULTI-seq.functions.R')
# try(source('https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/MULTI-seq.functions.R'))

# Requirements ------------------------
# try(require(MarkdownReportsDev),  silent = T)
# try(require(pheatmap),  silent = T)
# May also require


# BarTableSweepList -----------------------------------------------------------------------
BarTableSweepList <- function(min=0.01, max=0.99, step=0.02, bar_table =bar.table) {
  bar.table_sweep.list <- list()
  n <- 0
  Quantiles = seq(from = min, to = max, by=step)
  for (n in 1:length(Quantiles)) { # print(q)
    bar.table_sweep.list[[n]] <- classifyCells(bar_table, q=Quantiles[n])
    names(bar.table_sweep.list)[n] <- paste("q=",Quantiles[n], sep="")
  }
  return(bar.table_sweep.list)
}
# bar.table_sweep.list <- BarTableSweepList(bar_table = bar.table.solo) # Quantile Sweep List




# mSeq.map.all96.BCs ------------------------------------------------------------------------
mSeq.map.all96.BCs <- function(readTable = readTable, CellIDs = CellIDs
                               , path2allBCs = '~/Google_Drive/Science/IMBA/MULTI.seq/from.US/All.MULTI-seq_barcodes.Mar2019.tsv'
) {
  (bar.ref <- read_tsv(path2allBCs)[[1]]) # Vector of reference all MULTI-seq sample barcode sequences.
  MULTIseq.align(readTable = readTable, cellIDs = CellIDs, ref = bar.ref)
}
# bar.table <-mSeq.map.all96.BCs(readTable = readTable, CellIDs = CellIDs)


# aux.plotAllMseqBCs ------------------------------------------------------------------------

aux.plotAllMseqBCs <- function(bar.table = bar.table[,1:96], barcodes.used = BCs.used
                               , plotname = "Barcode seq depth") {
  stopifnot(is.numeric(BCs.used))
  BC.depth <- colSums(bar.table)[1:96]
  if (min(BC.depth) < 1) { BC.depth <- BC.depth+1 }
  log.depth <- log10(BC.depth); range(log.depth)
  ccc <- colnames(bar.table) %in% BCs.used

  wbarplot(log.depth, col = ccc, lwd=1, lcol=1, lty=2, plotname = plotname
           , vline = (range(BCs.used)+c(-1,1))
           , hline = quantile(log.depth[setdiff(1:69, BCs.used)], .95)
           , ylab = "log10(total reads / BC)", main =  plotname)
  wlegend.label(
    "    Horiz. line at 95%
    of unused BC's.
    Vertical line: range
    of used BCs.",poz = 2, cex=1)
}
# aux.plotAllMseqBCs(bar.table = bar.table[,1:96], barcodes.used = BCs.used, plotname = "Barcode seq depth")

# ------------------------------------------------------------------------
# bar.table.log <- t(log10(bar.table[,BCs.used]+1))
# bar.table.log <- clip.outliers(bar.table.log)

# p$'dist' <- c("correlation", "manhattan") # 'canberra', 'binary', 'minkowski', "euclidean")
# for (i in 1:l(p$'dist')) {
#   distX <- p$'dist'[i]; print(distX)
#   pheatmap(bar.table.log, clustering_distance_cols = distX# , cluster_rows = FALSE
#            , show_colnames = F, cutree_cols = nrow(bar.table.log), treeheight_col = 0
#            , filename = kpp("Barcode.log10.RC", distX,"heatmap", idate(Format = "%Y.%m.%d_%H.%M.%S"), "png"))
# }
# oo()


# subset.ReadTable <- function(variables) {
#
# }
# ------------------------------------------------------------------------
######################################################################
# plotting.dim.reduction.2D.R
######################################################################
# source('~/GitHub/Packages/Seurat.utils/Functions/plotting.dim.reduction.2D.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Plotting.dim.reduction.2D.R"))
# Source: self + web

# Requirements ------------------------
# library(plotly)
# try(source("~/GitHub/Packages/ggExpressDev/ggExpress.functions.R"), silent = T)
# try(source("https://raw.githubusercontent.com/vertesy/ggExpressDev/main/ggExpress.functions.R"), silent = T)

# May also require
# try (source('/GitHub/Packages/CodeAndRoll/CodeAndRoll.R'),silent= F) # generic utilities functions
# require('MarkdownReportsDev') # require("devtools") # plotting related utilities functions # devtools::install_github(repo = "vertesy/MarkdownReportsDev")


# Quick gene expression umap ------------------------------------------------------------------------
qUMAP <- function( feature= 'TOP2A', obj =  combined.obj  # The quickest way to draw a gene expression UMAP
                   , title = feature, sub =NULL
                   , reduct ="umap", splitby = NULL
                   , suffix = sub
                   , save.plot=T, PNG = T
                   , h=7, w=NULL, nr.cols = NULL
                   , assay = c("RNA","integrated")[1]
                   , HGNC.lookup= TRUE, make.uppercase = TRUE
                   , qlow = "q10", qhigh = "q90", ...) {

  if ( !(feature %in% colnames(obj@meta.data))) {
    feature <- check.genes(feature, verbose = F, HGNC.lookup = HGNC.lookup, makeuppercase = make.uppercase)
  }

  DefaultAssay(obj) <- assay
  ggplot.obj <- FeaturePlot(obj, features = feature
                            , reduction = reduct
                            , min.cutoff = qlow, max.cutoff = qhigh
                            # , plotname = ppp(toupper(reduct), feature)
                            , ncol = nr.cols
                            , split.by = splitby
                            , ...) + ggtitle(label = title, subtitle = sub)
  if (save.plot) {
    plotname <- ppp(toupper(reduct), feature)
    fname = ww.FnP_parser(ppp(plotname, assay, suffix), if (PNG) "png" else "pdf")
    try(save_plot(filename = fname, plot = ggplot.obj, base_height=h, base_width = w)) #, ncol=1, nrow=1
  }
  return(ggplot.obj)
}
# qUMAP('nFeature_RNA')
# qUMAP('VGLUT') # old name


# Quick clustering result or categorical umap  ------------------------------------------------------------------------
clUMAP <- function(ident = "integrated_snn_res.0.5", obj =  combined.obj   # The quickest way to draw a clustering result  UMAP
                   , reduct ="umap", splitby = NULL
                   , title = ident, sub =NULL, suffix = sub
                   , label.cex = 7
                   , h=7, w=NULL, nr.cols = NULL
                   , plotname = ppp(toupper(reduct), ident)
                   , cols = getDiscretePalette(ident.used = ident, show.colors = F)
                   , highlight.clusters = NULL, cells.highlight = NULL
                   , label = T, repel = T, legend = !label, MaxCategThrHP = 200
                   , save.plot=T, PNG = T
                   , save.object = F, ...) {
  IdentFound <- (ident %in%  colnames(obj@meta.data))

  if (!IdentFound) {
    ident <- GetClusteringRuns(obj = obj, pat = "_res.*[0,1]\\.[0-9]$")[1]
    iprint("Identity not found. Plotting", ident)
  }

  if ( !missing(highlight.clusters)) {
    x <- obj[[ident]]
    idx.ok <- x[,1] %in% highlight.clusters
    highlight.these <- rownames(x)[idx.ok]
  } else { highlight.these <- NULL}
  if ( !missing(cells.highlight)) {highlight.these <- cells.highlight} # overwrite, if directly defined



  NtCategs <- length(unique(obj[[ident]][,1]))
  if( NtCategs > MaxCategThrHP ) {
    iprint("Too many categories (",NtCategs,") in ", ident, "- use qUMAP for continous variables.")
  } else {
    if( length(unique(obj[[ident]])) < MaxCategThrHP )
      ggplot.obj <-
        DimPlot(object = obj, group.by = ident
                , cols = cols
                , reduction = reduct, split.by = splitby
                , ncol = nr.cols, cells.highlight = highlight.these
                , label = label, repel = repel, label.size = label.cex, ...) +
        ggtitle(label = title, subtitle = sub) +
        if (!legend) NoLegend() else NULL

    if (save.plot) {
      fname = ww.FnP_parser(ppp(plotname, suffix, kpp(highlight.clusters)), if (PNG) "png" else "pdf")
      try(save_plot(filename = fname, plot = ggplot.obj, base_height=h, base_width = w)) #, ncol=1, nrow=1
    }
    if(save.object) saveRDS(object = ggplot.obj, file = ppp(fname, 'ggobj.RDS'))
    return(ggplot.obj)
  } # if not too many categories
}
# clUMAP('cl.names.KnownMarkers.0.5' )
# clUMAP('cl.names.KnownMarkers.0.5', cols = NULL)




# ------------------------------------------------------------------------
gg_color_hue <- function(n) { # reproduce the ggplot2 default color palette
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
# https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette

# ------------------------------------------------------------------------
save2umaps.A4 <- function(plot_list, pname = F, suffix = NULL, scale = 1
                          , nrow = 2, ncol = 1
                          , h = hA4 * scale, w = wA4 * scale, ...) { # Save 2 umaps on an A4 page.
  if (pname ==F) pname = sppp(substitute(plot_list), suffix)
  p1 = plot_grid(plotlist = plot_list, nrow = nrow, ncol = ncol, labels = LETTERS[1:length(plot_list)], ...  )
  save_plot(plot = p1, filename = extPNG(pname), base_height = h, base_width = w)
}

# ------------------------------------------------------------------------
save4umaps.A4 <- function(plot_list, pname = F, suffix = NULL, scale = 1
                          , nrow = 2, ncol = 2
                          , h = wA4 * scale, w = hA4 * scale, ...) { # Save 4 umaps on an A4 page.
  if (pname==F) pname =  sppp(substitute(plot_list), suffix)
  p1 = plot_grid(plotlist = plot_list, nrow = nrow, ncol = ncol, labels = LETTERS[1:length(plot_list)], ...  )
  save_plot(plot = p1, filename = extPNG(pname), base_height = h, base_width = w)
}

# # ------------------------------------------------------------------------
# save4umaps.A4 <- function(plot_list, pname = F, suffix = NULL, scale = 1
#                           , nrow = 2, ncol = 2
#                           , h = wA4 * scale, w = hA4 * scale, ...) { # Save 4 umaps on an A4 page.
#   ww.saveXumaps(plot_list = plot_list, pname = pname, suffix = suffix, scale = scale
#                 , nrow = nrow, ncol = ncol
#                 , h = h, w =w, ...)
# }
#
# # ------------------------------------------------------------------------
# ww.saveXumaps <- function(plot_list, pname = F, suffix = NULL, scale = 1
#                           , nrow = 2, ncol = 2
#                           , h = wA4 * scale, w = hA4 * scale, ...) { # Save 4 umaps on an A4 page.
#   if (pname==F) pname =  sppp(substitute(plot_list), suffix)
#   p1 = plot_grid(plotlist = plot_list, nrow = nrow, ncol = ncol, labels = LETTERS[1:length(plot_list)], ...  )
#   save_plot(plot = p1, filename = extPNG(pname), base_height = h, base_width = w)
# }
#


# ------------------------------------------------------------------------
umapNamedClusters <- function(obj = combined.obj, metaD.colname = metaD.colname.labeled, ext = "png", ...) { # Plot and save umap based on a metadata column.
  fname = ppp("Named.clusters", metaD.colname, ext)
  p.named =
    DimPlot(obj, reduction = "umap", group.by = metaD.colname, label = T, ...) +
    NoLegend() +
    ggtitle(metaD.colname)
  save_plot(p.named, filename = fname); p.named
}
# umapNamedClusters(obj = combined.obj, metaD.colname = metaD.colname.labeled)


# ------------------------------------------------------------------------

# qqsave <- function(ggplot.obj# Quickly save a ggplot object, and optionally display it.
#                    , h=7, PNG =F, plotname = substitute(ggplot.obj), title=NULL, plotit = F) {
#   fname = ww.FnP_parser(plotname, if (PNG) "png" else "pdf")
#   save_plot(filename =fname, plot = ggplot.obj, base_height=h) #, ncol=1, nrow=1
#   if (plotit) ggplot.obj
# }

# ------------------------------------------------------------------------
qqSaveGridA4 <- function(plotlist= pl # Save 2 or 4 ggplot objects using plot_grid() on an A4 page
                         , plots = 1:2, NrPlots = length(plots), height = hA4, width = wA4
                         , fname = "Fractions.Organoid-to-organoid variation.png", ...) {
  stopifnot(NrPlots %in% c(2,4))
  iprint(NrPlots,"plots found,", plots,"are saved.")
  pg.cf = plot_grid(plotlist = plotlist[plots], nrow = 2, ncol = NrPlots/2, labels = LETTERS[1:NrPlots], ...  )
  if (NrPlots == 4) list2env(list(height = width, width = height), envir=as.environment(environment()))
  save_plot(filename = fname,
            plot = pg.cf, base_height = height, base_width = width)
  ww.FnP_parser(fname)
}
# qqSaveGridA4(plotlist= pl, plots = 1:2, fname = "Fractions.per.Cl.png")
# qqSaveGridA4(plotlist= pl, plots = 1:4, fname = "Fractions.per.Cl.4.png")

# ------------------------------------------------------------------------

# umapHiLightSel highlight a set of cells based on clusterIDs provided---------------
umapHiLightSel <- function(obj = combined.obj, # Highlight a set of cells based on clusterIDs provided.
                           COI =  c("0", "2", "4", "5",  "11"), res.cl = 'integrated_snn_res.0.3') {
  cellsSel = getCellIDs.from.meta(obj, values = COI, ident = res.cl)
  DimPlot(obj, reduction = "umap", group.by = res.cl,
          label = T, cells.highlight = cellsSel)
  ggsave(filename = extPNG(kollapse("cells",COI, collapseby = '.')))
}




# Save multiple FeaturePlot from a list of genes on A4 jpeg ------------------------
multiFeaturePlot.A4 <- function(list.of.genes # Save multiple FeaturePlots, as jpeg, on A4 for each gene, which are stored as a list of gene names.
                                , obj = combined.obj
                                , foldername = substitute(list.of.genes), plot.reduction='umap'
                                , intersectionAssay = c('RNA', 'integrated')[1]
                                , layout = c('tall', 'wide', FALSE )[2]
                                , colors=c("grey", "red"), nr.Col=2, nr.Row =4, cex = round(0.1/(nr.Col*nr.Row), digits = 2)
                                , gene.min.exp = 'q01', gene.max.exp = 'q99', subdir =T
                                , prefix = NULL , suffix = NULL
                                , saveGeneList = FALSE
                                , w = wA4, h = hA4, scaling = 1
                                , format = c('jpg', 'pdf', 'png')[1]
                                , ...
                                # , jpeg.res = 225, jpeg.q = 90
) {
  tictoc::tic()
  ParentDir = OutDir
  if (is.null(foldername)) foldername = "genes"
  if (subdir) create_set_SubDir( paste0(foldername,'-', plot.reduction),'/')
  list.of.genes.found = check.genes(list.of.genes = list.of.genes, obj = obj, assay.slot = intersectionAssay, makeuppercase = F)
  DefaultAssay(obj) <- intersectionAssay

  if (layout == 'tall') { w = wA4 * scaling; h = hA4 * scaling; nr.Col = 2; nr.Row = 4; print('layout active, nr.Col ignored.') }
  if (layout == 'wide') { w = hA4 * scaling; h = wA4 * scaling; nr.Col = 2; nr.Row = 2; print('layout active, nr.Col ignored.') }

  lsG = CodeAndRoll2::split_vec_to_list_by_N(1:length(list.of.genes.found), by = nr.Row * nr.Col)
  for (i in 1:length(lsG)) {
    genes = list.of.genes.found[lsG[[i]]]
    iprint(i,genes )
    plotname = kpp(c(prefix, plot.reduction,i, genes, suffix, format ))

    plot.list = FeaturePlot(object = obj, features = genes, reduction = plot.reduction, combine = F
                            , ncol = nr.Col, cols = colors
                            , min.cutoff = gene.min.exp, max.cutoff = gene.max.exp
                            , pt.size = cex, ...)

    for (i in 1:length(plot.list)) {
      plot.list[[i]] <- plot.list[[i]] + NoLegend() + NoAxes()
    }

    ggsave(filename = plotname, width = w, height = h,
           plot = cowplot::plot_grid(plotlist = plot.list, ncol = nr.Col, nrow = nr.Row)
    )
  }

  if (subdir) create_set_OutDir(... = ParentDir)
  if (saveGeneList) {
    if (is.null(obj@misc$gene.lists)) obj@misc$gene.lists <- list()
    obj@misc$gene.lists[[substitute(list.of.genes)]] <- list.of.genes.found
    print("Genes saved under: obj@misc$gene.lists")
    return(obj)
  }
  tictoc::toc()
};


# Save multiple FeatureHeatmaps from a list of genes on A4 jpeg -----------------------
# code for quantile: https://github.com/satijalab/seurat/blob/master/R/plotting_internal.R

multiSeuratHeatmap.A4 <- function(obj = combined.obj # Save multiple FeatureHeatmaps from a list of genes on A4 jpeg
                                   , list.of.genes, gene.per.page=5
                                   , group.cells.by= "batch", plot.reduction='umap'
                                   , cex = iround(3/gene.per.page), sep_scale = F
                                   , gene.min.exp = 'q5', gene.max.exp = 'q95'
                                   , jpeg.res = 225, jpeg.q = 90, ...) {

  tictoc::tic()
  list.of.genes = check.genes(list.of.genes, obj = obj)

  lsG = CodeAndRoll2::split_vec_to_list_by_N(1:length(list.of.genes), by=gene.per.page)
  for (i in 1:length(lsG)) { print(i )
    genes = list.of.genes[lsG[[i]]]
    plotname = kpp(c("FeatureHeatmap",plot.reduction,i, genes, 'jpg' ))
    print(plotname)
    jjpegA4(plotname, r = jpeg.res, q = jpeg.q)
    try(
      FeatureHeatmap(obj, features.plot =genes , group.by = group.cells.by
                     , reduction.use = plot.reduction, do.return = F
                     , sep.scale = sep_scale, min.exp = gene.min.exp, max.exp = gene.max.exp
                     , pt.size = cex, key.position = "top", ...)
      , silent = F
    )
    try.dev.off()
  }
  tictoc::toc()
}


# plot.UMAP.tSNE.sidebyside ---------------------------------------------------------------------

plot.UMAP.tSNE.sidebyside <- function(obj = combined.obj, grouping = 'res.0.6',  # Plot a UMAP and tSNE sidebyside
                                      no_legend = F,
                                      do_return = TRUE,
                                      do_label = T,
                                      label_size = 10,
                                      vector_friendly = TRUE,
                                      cells_use = NULL,
                                      no_axes = T,
                                      pt_size = 0.5,
                                      name.suffix = NULL,
                                      width = hA4, heigth = 1.75*wA4, filetype = "pdf", ...) {

  p1 <- DimPlot(object = obj, reduction.use = "tsne", no.axes = no_axes, cells.use = cells_use
                , no.legend = no_legend, do.return = do_return, do.label = do_label, label.size = label_size
                , group.by = grouping, vector.friendly = vector_friendly, pt.size = pt_size, ...) +
    ggtitle("tSNE") + theme(plot.title = element_text(hjust = 0.5))

  p2 <- DimPlot(object = obj, reduction.use = "umap", no.axes = no_axes, cells.use = cells_use
                , no.legend = T, do.return = do_return, do.label = do_label, label.size = label_size
                , group.by = grouping, vector.friendly = vector_friendly, pt.size = pt_size, ...) +
    ggtitle("UMAP") + theme(plot.title = element_text(hjust = 0.5))

  plots = plot_grid(p1, p2, labels=c("A", "B"), ncol = 2)
  plotname=kpp( 'UMAP.tSNE', grouping, name.suffix, filetype)

  cowplot::save_plot(filename = plotname, plot = plots
                     , ncol = 2 # we're saving a grid plot of 2 columns
                     , nrow = 1 # and 2 rows
                     , base_width = width
                     , base_height = heigth
                     # each individual subplot should have an aspect ratio of 1.3
                     # , base_aspect_ratio = 1.5
  )
}

# PlotTopGenesPerCluster --------------------------------------------------------------------------------
PlotTopGenesPerCluster <- function(obj = combined.obj, cl_res = res, nrGenes = p$'n.markers'
                                   , order_by = c("combined.score","avg_log2FC", "p_val_adj")[1]
                                   , df_markers = combined.obj@misc$"df.markers"[[paste0("res.",cl_res)]]) {
  topX.markers <- GetTopMarkers(df = df_markers,  n= nrGenes
                                , order.by = order_by )
  ls.topMarkers <-  splitbyitsnames(topX.markers)
  for (i in 1:l(ls.topMarkers)) {
    multiFeaturePlot.A4(list.of.genes = ls.topMarkers[[i]], obj = obj, subdir = F
                        , prefix = ppp("DEG.markers.res",cl_res,"cluster",names(ls.topMarkers)[i]))
  }

}
# PlotTopGenesPerCluster(obj = combined.obj, cl_res = 0.5, nrGenes = p$'n.markers')



# qFeatureScatter --------------------------------------------------------------------------------

qFeatureScatter <- function(feature1 = "M.RabV.N2c", feature2 = "P.RabV.N2c", obj = combined.obj
                            , ext ="png", plot = TRUE, ...) {
  plotname <- kpp(feature1,"VS", feature2)
  p <- FeatureScatter(object = obj, feature1 = feature1, feature2 = feature2, ...)
  fname = kpp("FeatureScatter", plotname)
  qqSave(ggobj = p, title = plotname, ext = ext)
  if (plot) p
}

# qMarkerCheck.BrainOrg --------------------------------------------------------------------------------
qMarkerCheck.BrainOrg <- function(obj = combined.obj, custom.genes = F) {
  Signature.Genes.Top16 <- if (custom.genes) custom.genes else
  {
    Signature.Genes.Top16  <- c(
      `S-phase` = "TOP2A", `G2M-phase` = "HIST1H4C"
      , `oRG` = "ID4", `oRG` = "HOPX" # oRG outer radial glia
      , `Intermediate progenitor` = "EOMES",  `Intermediate progenitor1` = "TAC3"
      , Astroglia = "GFAP", Astrocyte = "S100B"
      , `Immature neurons` = "SLA", Interneurons = "DLX6-AS1"
      , `Hypoxia/Stress` = "DDIT4", Glycolytic = "PDK1"
      , `Low-Quality` = "POLR2A", `Choroid.Plexus` = "DCN"
      , `dl-EN` = "KAZN", `ul-EN` = "SATB2" # dl-EN = deep layer excitatory neuron
      # , `Choroid.Plexus` = "OTX2", `Choroid.Plexus` = "BMP4"
    )
  }
  print(as_tibble_from_named_vec(Signature.Genes.Top16))
  multiFeaturePlot.A4(obj = obj, list.of.genes = Signature.Genes.Top16, layout = "tall")
}
# qMarkerCheck.BrainOrg(combined.obj)



# getDiscretePalette --------------------------------------------------------------------------------
getDiscretePalette <- function(ident.used = GetClusteringRuns()[1]
                               , obj = combined.obj
                               , palette.used = c("alphabet", "alphabet2", "glasbey", "polychrome", "stepped")[1]
                               , show.colors = F) {
  n.clusters <-  nrow(unique(obj[[ident.used]]))
  colz <- DiscretePalette(n = n.clusters, palette = palette.used)
  if (show.colors) MarkdownHelpers::color_check(colz)
  return(colz)
}
# getDiscretePalette()

######################################################################
# plotting.dim.reduction.3D.R
######################################################################
# source('~/GitHub/Packages/Seurat.utils/Functions/plotting.dim.reduction.3D.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Plotting.dim.reduction.3D.R"))
# Source: self + https://github.com/Dragonmasterx87/Interactive-3D-Plotting-in-Seurat-3.0.0

# Requirements ------------------------
# try(library(plotly), silent = T)
# try(library(MarkdownReportsDev), silent = T)
# try(library(htmlwidgets), silent = T)

# May also require
# try (source('~/GitHub/Packages/CodeAndRoll/CodeAndRoll.R'),silent= T) # generic utilities functions
# require('MarkdownReportsDev') # require("devtools") # plotting related utilities functions # devtools::install_github(repo = "vertesy/MarkdownReportsDev")


# ------------------------------------------------------------------------
ww.check.if.3D.reduction.exist <- function(obj = obj) { # ww.check.if.3D.reduction.exist in backup slot
  if( !("UMAP_3" %in% colnames(obj@reductions$'umap'))) {
    stopif( is.null(combined.obj@misc$reductions.backup$'umap3d')
             , "No 3D umap found in backup slot, @misc$reductions.backup. Run SetupReductionsNtoKdimensions() first.")
    RecallReduction(obj = obj, dim = 3, reduction = "umap")
  } else { # Reduction found in normal UMAP slot
    obj
  }
}

# ww.check.quantile.cutoff ------------------------------------------------------------------------
ww.check.quantile.cutoff.and.clip.outliers <- function(expr.vec = plotting.data[,gene], quantileCutoffX = quantileCutoff, min.cells.expressing = 10) {
  expr.vec.clipped <- clip.outliers(expr.vec, probs = c(1 - quantileCutoffX, quantileCutoffX))
  if( sum(expr.vec.clipped > 0) > min.cells.expressing ){
    expr.vec <- expr.vec.clipped
  } else {
    iprint("WARNING: quantile.cutoff too stringent, would leave <", min.cells.expressing, "cells. It is NOT applied.")
  }
  return(expr.vec)
}

# plot3D.umap.gene ------------------------------------------------------------------------
plot3D.umap.gene <- function(gene="TOP2A", obj=combined.obj # Plot a 3D umap with gene expression. Uses plotly. Based on github.com/Dragonmasterx87.
                             , quantileCutoff = .99, def.assay = c("integrated", "RNA")[2]
                             , suffix = NULL, AutoAnnotBy = GetNamedClusteringRuns(obj)[1]
                             , alpha = .5, dotsize=1.25 ){
  # stopifnot(AutoAnnotBy %in% colnames(obj@meta.data) | AutoAnnotBy = FALSE)

  obj <- ww.check.if.3D.reduction.exist(obj = obj)
  stopifnot((gene %in% rownames(obj) | gene %in% colnames(obj@meta.data)))
  DefaultAssay(object = obj) <- def.assay; iprint(DefaultAssay(object = obj), "assay")

  plotting.data <- FetchData(object = obj, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "Expression" = gene), slot = 'data')

  plotting.data$'Expression' <- ww.check.quantile.cutoff.and.clip.outliers(expr.vec = plotting.data[,gene], quantileCutoffX = quantileCutoff, min.cells.expressing = 10)
  clip.outliers(plotting.data[,gene], probs = c(1 - quantileCutoff, quantileCutoff))
  plotting.data$'label' <- paste(rownames(plotting.data), " - ", plotting.data[,gene], sep = "")

  ls.ann.auto <- if (AutoAnnotBy != FALSE) {
    Annotate4Plotly3D(obj. = obj, plotting.data. = plotting.data, AnnotCateg = AutoAnnotBy)
  } else { NULL }

  plt <- plot_ly(data = plotting.data
                 , x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3
                 , type = "scatter3d"
                 , mode = "markers"
                 , marker = list(size = dotsize)
                 , text =  ~label
                 , color = ~Expression
                 , opacity = alpha
                 # , colors = c('darkgrey', 'red')
                 , colorscale='Viridis'
                 #, hoverinfo="text"
  ) %>% layout(title = gene, scene = list(annotations = ls.ann.auto))
  SavePlotlyAsHtml(plt, category. = gene, suffix. = suffix)
  return(plt)
}
# plot3D.umap.gene(obj = combined.obj, gene = "DDIT4", quantileCutoff = .95)
# plot3D.umap.gene(obj = combined.obj, gene = "percent.mito", quantileCutoff = .95) # for continous meta variables
# plot3D.umap.gene(obj = combined.obj, gene = "nFeature_RNA", quantileCutoff = .95) # for continous meta variables


# plot3D.umap ------------------------------------------------------------------------
plot3D.umap <- function(category="v.project", obj=combined.obj # Plot a 3D umap based on one of the metadata columns. Uses plotly. Based on github.com/Dragonmasterx87.
                        , suffix = NULL, AutoAnnotBy = GetNamedClusteringRuns(obj)[1]
                        , dotsize = 1.25) {

  stopifnot(category %in% colnames(obj@meta.data))
  obj <- ww.check.if.3D.reduction.exist(obj = obj)

  plotting.data <- FetchData(object = obj, vars = c("UMAP_1", "UMAP_2", "UMAP_3", category))
  colnames(plotting.data)[4] = "category"
  plotting.data$label <- paste(rownames(plotting.data))   # Make a column of row name identities (these will be your cell/barcode names)

  ls.ann.auto <- if (AutoAnnotBy != FALSE) {
    Annotate4Plotly3D(obj. = obj, plotting.data. = plotting.data, AnnotCateg = AutoAnnotBy)
  } else { NULL }

  plt <- plot_ly(data = plotting.data
                 , x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3
                 , type = "scatter3d"
                 , mode = "markers"
                 , marker = list(size = dotsize)
                 , text = ~label
                 , color = ~category
                 , colors = gg_color_hue(length(unique(plotting.data$'category')))
                 # , hoverinfo="text"
  ) %>% layout(title = category, scene = list(annotations = ls.ann.auto))
  SavePlotlyAsHtml(plt, category. = category, suffix. = suffix)
  return(plt)
}
# plot3D.umap(combined.obj, category = "Phase")

# ------------------------------------------------------------------------
SavePlotlyAsHtml <- function(plotly_obj, category.=category, suffix. = NULL) { # Save Plotly 3D scatterplot as an html file.
  OutputDir <- if (exists("OutDir")) OutDir else getwd()
  name.trunk <- kpp("umap.3D", category., suffix., idate(), "html")
  fname <- kpps(OutputDir, name.trunk)
  iprint("Plot saved as:", fname)
  htmlwidgets::saveWidget(plotly_obj, file = fname, selfcontained = TRUE, title = category.)
}


# ------------------------------------------------------------------------
BackupReduction <- function(obj = combined.obj, dim=2, reduction="umap") { # Backup UMAP to `obj@misc$reductions.backup` from `obj@reductions$umap`.
  if (is.null(obj@misc$"reductions.backup")) obj@misc$"reductions.backup" <- list()
  dslot = paste0(reduction,dim,"d")
  obj@misc$reductions.backup[[dslot]] <- obj@reductions[[reduction]]
  return(obj)
}
# Example
# obj <- BackupReduction(obj = obj, dim=2, reduction=umap"")

# ------------------------------------------------------------------------
SetupReductionsNtoKdimensions <- function(obj = combined.obj, nPCs = p$'n.PC', dimensions=3:2, reduction="umap", ...) { # Calculate N-to-K dimensional umaps (default = 2:3); and back them up UMAP to `obj@misc$reductions.backup` from @reductions$umap
  red <- reduction
  for (d in dimensions) {
    iprint(d, "dimensional", red, "is calculated")
    obj <- if (reduction == "umap") {
      RunUMAP(obj, dims = 1:nPCs, n.components = d, ...)
    } else if (reduction == "tsne") {
      RunTSNE(obj, dims = 1:nPCs, n.components = d, ...)
    } else if (reduction == "pca") {
      RunPCA(obj, dims = 1:nPCs, n.components = d, ...)
    }
    obj <- BackupReduction(obj = obj, dim = d, reduction = red)
  }
  return(obj)
}
# Example
# combined.obj <- SetupReductionsNtoKdimensions(obj = combined.obj, nPCs = p$'n.PC', dimensions=2:3, reduction="umap")
# qUMAP()

# ------------------------------------------------------------------------
RecallReduction <- function(obj = combined.obj, dim=2, reduction="umap") { # Set active UMAP to `obj@reductions$umap` from `obj@misc$reductions.backup`.
  dslot = paste0(reduction,dim,"d")
  reduction.backup <- obj@misc$reductions.backup[[dslot]]
  msg <-  paste(dim, "dimensional", reduction, "from obj@misc$reductions.backup" )
  stopif(is.null(reduction.backup), message = p0(msg," is NOT FOUND")); iprint(msg, "is set active. " )
  stopifnot(dim == ncol(reduction.backup))
  obj@reductions[[reduction]] <- reduction.backup
  return(obj)
}
# Example
# combined.obj <- RecallReduction(obj = combined.obj, dim=2, reduction="umap")
# qUMAP()
# combined.obj <- RecallReduction(obj = combined.obj, dim=3, reduction="umap")
# qUMAP()


# ------------------------------------------------------------------------
Annotate4Plotly3D <- function(obj. = combined.obj # Create annotation labels for 3D plots. Source https://plot.ly/r/text-and-annotations/#3d-annotations
                              , plotting.data. = plotting.data
                              , AnnotCateg = AutoAnnotBy) {
  stopifnot(AnnotCateg %in% colnames(obj.@meta.data))

  plotting.data.$'annot' <- FetchData(object = obj., vars = c(AnnotCateg))[,1]
  auto_annot <-
    plotting.data. %>%
    group_by(annot) %>%
    summarise(showarrow = F
              , xanchor = "left"
              , xshift = 10
              , opacity = 0.7
              , "x" = mean(UMAP_1)
              , "y" = mean(UMAP_2)
              , "z" = mean(UMAP_3)
    )
  names(auto_annot)[1] = "text"
  ls.ann.auto = apply(auto_annot, 1, as.list)
  return(ls.ann.auto)
}

# ------------------------------------------------------------------------
Plot3D.ListOfGenes <- function(obj = combined.obj # Plot and save list of 3D UMAP ot tSNE plots using plotly.
                               , annotate.by = "integrated_snn_res.0.7", opacity = 0.5, cex = 1.25, default.assay = c("integrated", "RNA")[2]
                               , ListOfGenes = c("BCL11B" , "FEZF2", "EOMES", "DLX6-AS1", "HOPX", "DDIT4")
                               , SubFolderName=ppp("plot3D", substitute(ListOfGenes))) {


  try(create_set_SubDir(SubFolderName))
  obj. <- obj; rm("obj")
  stopifnot(annotate.by %in% c(colnames(obj.@meta.data), FALSE))

  DefaultAssay(object = obj.) <- default.assay
  MissingGenes <- setdiff(ListOfGenes, rownames(obj.))
  if ( length(MissingGenes)) iprint("These genes are not found, and omitted:", MissingGenes, ". Try to change default assay.")
  ListOfGenes <- intersect(ListOfGenes, rownames(obj.))

  for (i in 1:length(ListOfGenes)) {
    g <- ListOfGenes[i]; print(g)
    plot3D.umap.gene(obj = obj., gene = g, AutoAnnotBy = annotate.by, alpha = opacity, def.assay = default.assay, dotsize = cex)
  }
  try(oo())
  try(create_set_Original_OutDir(NewOutDir = ParentDir))
}
# CellTypeMarkers <- c(  "PGK1", "CTIP2" = "BCL11B" , "FEZF2", "EOMES", "DLX6-AS1", "HOPX", "DDIT4","TOP2A", "PTGDS", "EDNRB", "EGFR", "SCGN", "NR2F2", "EMX2", "GAD2", "DLX2", "SATB2")
# Plot3D.ListOfGenes(obj = combined.obj, ListOfGenes = CellTypeMarkers)


# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
Plot3D.ListOfCategories <- function(obj = combined.obj # Plot and save list of 3D UMAP ot tSNE plots using plotly.
                                    , annotate.by = "integrated_snn_res.0.7", cex = 1.25, default.assay = c("integrated", "RNA")[2]
                                    , ListOfCategories=c("v.project","experiment", "Phase", "integrated_snn_res.0.7")
                                    , SubFolderName=ppp("plot3D", substitute(ListOfCategories))) {

  try(create_set_SubDir(SubFolderName))
  obj. <- obj; rm("obj")
  stopifnot(annotate.by %in% colnames(obj.@meta.data))
  DefaultAssay(object = obj.) <- default.assay

  MissingCateg <- setdiff(ListOfCategories, colnames(obj.@meta.data))
  if ( length(MissingCateg)) iprint("These metadata categories are not found, and omitted:", MissingCateg, ". See colnames(obj@meta.data).")
  ListOfCategories <- intersect(ListOfCategories, colnames(obj.@meta.data))

  for (i in 1:length(ListOfCategories)) {
    categ <- ListOfCategories[i]; print(categ)
    plot3D.umap(obj = obj., category = categ, AutoAnnotBy = annotate.by, dotsize = cex)
  }
  try(oo())
  try(create_set_Original_OutDir(NewOutDir = ParentDir))
}
# categ3Dplots <- c("v.project","experiment", "Phase", "integrated_snn_res.0.7", "Area", "Individual", "Type")
# Plot3D.ListOfCategories(obj = combined.obj, ListOfCategories = categ3Dplots)


# ------------------------------------------------------------------------
# ------------------------------------------------------------------------


######################################################################
# plotting.filtering.R
######################################################################
# source('~/GitHub/Packages/Seurat.utils/Functions/plotting.filtering.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Plotting.filtering.R"))


# PlotFilters ------------------------------------------------------------------------------------
PlotFilters <- function(ls.obj = ls.Seurat # Plot filtering threshold and distributions, using four panels to highlight the relation between Gene- and UMI-count, ribosomal- and mitochondrial-content.
                        , parentdir= OutDirOrig
                        , suffices = names(ls.obj)
                        , filetype='.png'
                        , below.mito = p$"thr.lp.mito"
                        , above.mito = p$"thr.hp.mito"
                        , below.ribo = p$"thr.lp.ribo"
                        , above.ribo = p$"thr.hp.ribo"
                        , below.nFeature_RNA = p$"thr.lp.nFeature_RNA"
                        , above.nFeature_RNA = p$"thr.hp.nFeature_RNA"
                        , subdir= kpp("Filtering.plots"
                                      , "mito", p$"thr.hp.mito", p$"thr.lp.mito"
                                      , "ribo", p$"thr.hp.ribo", p$"thr.lp.ribo"
                                      , "nFeature", p$"thr.hp.nFeature_RNA", p$"thr.lp.nFeature_RNA", "/")
                        , transparency = 0.25
                        , cex = 0.75
                        , theme.used = theme_bw(base_size = 18)
                        , LabelDistFromTop = 200 # for barplot_label
) {

  llprint(
    "We filtered for high quality cells based on the number of genes detected [", above.nFeature_RNA, ";" ,below.nFeature_RNA
    , "] and the fraction of mitochondrial [", percentage_formatter(above.mito), ";" ,percentage_formatter(below.mito)
    , "] and ribosomal [",percentage_formatter(above.ribo), ";" ,percentage_formatter(below.ribo), "] reads."
  )


  theme_set(theme.used)
  create_set_OutDir(parentdir, subdir)
  # require(ggplot2)
  if (suffices == l(ls.obj)) print("ls.Obj elements have no names (required).")

  for (i in 1:l(ls.obj)) {
    print(suffices[i])
    mm =  ls.obj[[i]]@meta.data

    AllMetaColumnsPresent <- all(c('nFeature_RNA', 'percent.mito', 'percent.ribo') %in% colnames(mm))
    if (!AllMetaColumnsPresent) {
      print(c('nFeature_RNA', 'percent.mito', 'percent.ribo'))
      print(c('nFeature_RNA', 'percent.mito', 'percent.ribo') %in% colnames(mm))
      print("Try to run:")
      print('objX <- addMetaFraction(obj = objX, col.name = "percent.mito", gene.symbol.pattern =  "^MT\\.|^MT-")')
      print('objX <- addMetaFraction(obj = objX, col.name = "percent.ribo", gene.symbol.pattern =  "^RPL|^RPS")')
      stop()
    }



    filt.nFeature_RNA = (mm$'nFeature_RNA' < below.nFeature_RNA & mm$'nFeature_RNA' > above.nFeature_RNA)
    filt.below.mito = (mm$'percent.mito' < below.mito & mm$'percent.mito' > above.mito)

    # filt.below.mito = (mm$'percent.mito' < below.mito)
    filt.below.ribo = (mm$'percent.ribo' < below.ribo & mm$'percent.ribo' > above.ribo)

    mm =  cbind(mm, filt.nFeature_RNA, filt.below.mito, filt.below.ribo)

    mm$colour.thr.nFeature <- cut(mm$'nFeature_RNA',
                                  breaks = c(-Inf, above.nFeature_RNA, below.nFeature_RNA, Inf),
                                  labels = c(p0("LQ (<", above.nFeature_RNA,")"),
                                             p0("HQ (", above.nFeature_RNA,"< X <", below.nFeature_RNA,")"),
                                             p0("Dbl/Outlier (>", below.nFeature_RNA,")")
                                  )
    )

    A = ggplot(data = mm, aes(x = nFeature_RNA, fill = colour.thr.nFeature)) +
      geom_histogram(binwidth = 100) +
      ggtitle(paste("Cells between", above.nFeature_RNA,"and",below.nFeature_RNA, " UMIs are selected (", pc_TRUE(filt.nFeature_RNA), ")")) +
      geom_vline(xintercept = below.nFeature_RNA) +
      geom_vline(xintercept = above.nFeature_RNA);
    # A

    B = ggplot2::ggplot(mm, aes(x = nFeature_RNA, y = percent.mito)) +
      ggplot2::ggtitle(paste("Cells below", percentage_formatter(below.mito),
                             "mito reads are selected (with A:", pc_TRUE(filt.nFeature_RNA & filt.below.mito), ")")) +
      ggplot2::geom_point(alpha = transparency, size = cex,  show.legend = FALSE,
                          aes(color = filt.nFeature_RNA & filt.below.mito)  ) +
      scale_x_log10() + # scale_y_log10() +
      # annotation_logticks() +
      geom_hline(yintercept = below.mito) +
      geom_hline(yintercept = above.mito) +
      geom_vline(xintercept = below.nFeature_RNA) +
      geom_vline(xintercept = above.nFeature_RNA);
    # B


    C = ggplot(mm, aes(x = nFeature_RNA, y = percent.ribo)) +
      ggtitle(paste("Cells below", percentage_formatter(below.ribo),
                    "ribo reads are selected (with A:"
                    , pc_TRUE(filt.nFeature_RNA & filt.below.ribo), ")")) +
      geom_point(alpha = transparency, size = cex,   show.legend = FALSE,
                 aes(color = filt.nFeature_RNA & filt.below.ribo)  ) +
      scale_x_log10() + # scale_y_log10() +
      annotation_logticks() +
      geom_hline(yintercept = below.ribo) +
      geom_hline(yintercept = above.ribo) +
      geom_vline(xintercept = below.nFeature_RNA) +
      geom_vline(xintercept = above.nFeature_RNA);
    # C


    D = ggplot(mm, aes(x = percent.ribo, y = percent.mito)) +
      ggtitle(paste("Cells w/o extremes selected (with A,B,C:"
                    , pc_TRUE(filt.nFeature_RNA & filt.below.mito & filt.below.ribo), ")")) +

      geom_point(alpha = transparency, size =  cex,  show.legend = FALSE,
                 aes(color = filt.nFeature_RNA & filt.below.mito & filt.below.ribo)  ) +
      scale_x_log10() + scale_y_log10() +
      annotation_logticks() +
      geom_hline(yintercept = below.mito) +
      geom_hline(yintercept = above.mito) +
      geom_vline(xintercept = below.ribo) +
      geom_vline(xintercept = above.ribo);
    # D


    plot_list = list(A,B,C,D)
    px = plot_grid(plotlist = plot_list, nrow = 2, ncol = 2, labels = LETTERS[1:4])
    fname = ppp("Filtering.thresholds", suffices[i], filetype)
    save_plot(filename = fname, plot = px, base_height = 12, ncol = 1, nrow = 1) #Figure 2

  } # for

  # End ------------------------------------------------------------------------
  create_set_Original_OutDir()
}
# PlotFilters(ls.Seurat)


######################################################################
# plotting.statistics.and.QC.R
######################################################################
# source('~/GitHub/Packages/Seurat.utils/Functions/Plotting.statistics.and.QC.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Plotting.statistics.and.QC.R"))

# Source: self + web

# Requirements ------------------------
# require(Seurat)
# require(ggplot2)
# tools for tools::toTitleCase

# May also require
# try (source('/GitHub/Packages/CodeAndRoll/CodeAndRoll.R'),silent= F) # generic utilities functions
# require('MarkdownReportsDev') # require("devtools") # plotting related utilities functions # devtools::install_github(repo = "vertesy/MarkdownReportsDev")

# PCA percent of variation associated with each PC ------------------------------------------------------------
scCalcPCAVarExplained <- function(obj =  combined.obj) { # Determine percent of variation associated with each PC.
  pct <- obj@reductions$pca@stdev / sum(obj@reductions$pca@stdev) * 100
  names(pct) =1:length(obj@reductions$pca@stdev)
  return(pct)
}

# plot percent of variation associated with each PC ------------------------------------------------------------
scPlotPCAvarExplained <- function(obj =  combined.obj, use.MDrep = F) { # Plot the percent of variation associated with each PC.
  pct <- scCalcPCAVarExplained(obj)
  if (use.MDrep) {
    wbarplot(pct , xlab = "Principal Components", ylab = "% of variation explained")
    barplot_label(round(pct, digits = 2), barplotted_variable = pct, cex = .5 )
  } else {
    qbarplot(vec = pct, xlab = "Principal Components", ylab =  "% of variation explained", w = 10, h = 5, hline = 1 )
  }
}


# BarplotCellsPerObject ------------------------------------------------------------

BarplotCellsPerObject <- function(ls.Seu = ls.Seurat, # Take a List of Seurat objects and draw a barplot for the number of cells per object.
                                  plotname="Nr.Cells.After.Filtering", names=F ) {
  cellCounts = unlapply(ls.Seu, ncol)
  names(cellCounts) = if (length(names) == length(ls.Seurat)) names else names(ls.Seurat)
  wbarplot(cellCounts, plotname = plotname,tilted_text = T, ylab="Cells")
  barplot_label(cellCounts, TopOffset = 500, w = 4)
}

# CellFractionsBarplot2 ------------------------------------------------------------
CellFractionsBarplot2 <- function(obj = combined.obj
                                  , group.by = "integrated_snn_res.0.5.ordered", fill.by = "age", downsample = T
                                  , plotname = paste(tools::toTitleCase(fill.by), "proportions"), hlines = c(.25, .5, .75), seedNr = 1989) {
  set.seed(seedNr)
  pname.suffix <- capt.suffix <- NULL
  if (downsample) {
    downsample <- min (table(obj@meta.data[[fill.by]]))
    pname.suffix <- "(downsampled)"
    capt.suffix <- paste("Downsampled to", downsample, "cells in the smallest", fill.by, "group.")
  }
  caption_ <- paste("Numbers denote # cells.", capt.suffix)
  pname_ <- paste(plotname, pname.suffix)

  obj@meta.data %>%
    group_by( (!!as.name(fill.by)) ) %>%
    { if (downsample) sample_n(., downsample) else . } %>%
    group_by( (!!as.name(group.by)) ) %>%

    ggplot( aes(fill = (!!(as.name(fill.by))),  x = (!!(as.name(group.by)))) ) +
    geom_hline( yintercept = hlines, lwd=1.5)  +
    geom_bar( position = "fill" ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_text(aes(label = ..count..), stat='count',position = position_fill(vjust = 0.5)) +
    labs(title = pname_,  x = "Clusters", y = "Fraction", caption = caption_) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
}
# CellFractionsBarplot2(obj = combined.obj, group.by = "integrated_snn_res.0.1", fill.by = "Phase", downsample = T)
# CellFractionsBarplot2(obj = combined.obj, group.by = "integrated_snn_res.0.1", fill.by = "Phase", downsample = F)



#  ------------------------------------------------
barplot.cells.per.cluster <- function(obj = combined.obj, ident =  "cl.names.KnownMarkers.0.5", sort = F) {
  cell.per.cluster <- (table(obj[[ident]][,1]))
  if (sort) cell.per.cluster <- sort(cell.per.cluster)
  qbarplot(cell.per.cluster, subtitle = ident, suffix = ident
           , col = rainbow(l(cell.per.cluster))
           , xlab.angle = 45
           # , col = getClusterColors(ident = ident, show = T)
           , palette_use = NULL, )
}
# barplot.cells.per.cluster()
# barplot.cells.per.cluster(sort=T)



# BulkGEScatterPlot ------------------------------------------------------------------------
BulkGEScatterPlot <- function(obj = combined.obj # Plot bulk scatterplots to identify differential expressed genes across conditions
                              , clusters = "cl.names.KnownMarkers.0.2", TwoCategIdent = 'age', genes.from.bulk.DE = rownames(df.markers.per.AGE)) {

  (SplitIdents <- unique(obj[[TwoCategIdent]][,1]))
  stopifnot(length(SplitIdents) == 2)

  Idents(obj) <- clusters
  IdentsUsed <- gtools::mixedsort(as.character(unique(Idents(obj))))
  NrPlots <- length(IdentsUsed)
  p.clAv <- p.clAv.AutoLabel <- genes.to.label <- list.fromNames(IdentsUsed)

  # i = 1
  for (i in 1:NrPlots) {
    print(IdentsUsed[i])
    ClX <- subset(obj, idents = IdentsUsed[i])
    Idents(ClX) <- TwoCategIdent
    avg.ClX.cells <- log2(AverageExpression(ClX, verbose = FALSE)$RNA + 1)
    avg.ClX.cells$gene <- rownames(avg.ClX.cells)

    # plot ----
    p.clAv[[i]] <- p.clAv.AutoLabel[[i]] <-
      ggplot(avg.ClX.cells, aes(x = !!as.name(SplitIdents[1]), y = !!as.name(SplitIdents[2]) )) +
      geom_point(data = avg.ClX.cells, color = rgb(0, .5, 0, 0.25), size = 1) +
      FontSize(x.title = 8, x.text = 8, y.title = 8, y.text = 8)+
      geom_abline(slope = 1, intercept = 0, color = 'grey') +
      ggtitle(paste("Cluster", IdentsUsed[i] )) +
      # ggtitle(paste0("Cluster ", i) ) +
      scale_x_log10() + scale_y_log10() + annotation_logticks()
    # p.clAv[[i]]

    "Auto identify divergent genes"
    dist.from.axis = eucl.dist.pairwise(avg.ClX.cells[,1:2])
    genes.to.label[[i]] = names(head(sort(dist.from.axis, decreasing = T),n = 20))
    p.clAv.AutoLabel[[i]] <- LabelPoints(plot = p.clAv[[i]], points = genes.to.label[[i]], xnudge = 0, ynudge = 0, repel = TRUE, size=2);
    p.clAv.AutoLabel[[i]]

    "Pre-identified genes"
    p.clAv[[i]] <- LabelPoints(plot = p.clAv[[i]], points = genes.from.bulk.DE, repel = TRUE, size=2);
  }

  PlotIter <- CodeAndRoll2::split_vec_to_list_by_N(1:NrPlots, by = 4)
  for (i in 1:length(PlotIter)) {
    plotLS = p.clAv.AutoLabel[PlotIter[[i]]]
    qqSaveGridA4(plotlist = plotLS, plots = 1:4, fname = ppp("BulkGEScatterPlot.AutoGenes",kpp(PlotIter[[i]]), "png"))

    plotLS = p.clAv[PlotIter[[i]]]
    qqSaveGridA4(plotlist = plotLS, plots = 1:4, fname= ppp("BulkGEScatterPlot.BulkGenes",kpp(PlotIter[[i]]), "png"))
  }
}
# BulkGEScatterPlot(obj = combined.obj, clusters = "cl.names.KnownMarkers.0.2", TwoCategIdent = 'age', genes.from.bulk.DE = rownames(df.markers.per.AGE))


# qqSaveGridA4 ------------------------------------------------------------------------------------
qqSaveGridA4 <- function(plotlist= pl # Save 2 or 4 ggplot objects using plot_grid() on an A4 page
                         , plots = 1:2, NrPlots = length(plots), height = hA4, width = wA4
                         , fname = "Fractions.Organoid-to-organoid variation.png") {
  stopifnot(NrPlots %in% c(2,4))
  iprint(NrPlots,"plots found,", plots,"are saved.")
  pg.cf = plot_grid(plotlist = plotlist[plots], nrow = 2, ncol = NrPlots/2, labels = LETTERS[1:NrPlots]  )
  if (NrPlots == 4) list2env(list(height = width, width = height), envir=as.environment(environment()))
  save_plot(filename = fname,
            plot = pg.cf, base_height = height, base_width = width)
  ww.FnP_parser(fname)
}
# qqSaveGridA4(plotlist= pl, plots = 1:2, fname = "Fractions.per.Cl.png")
# qqSaveGridA4(plotlist= pl, plots = 1:4, fname = "Fractions.per.Cl.4.png")



# ------------------------
#' sparse.cor
#' Correlation calculation for sparse matrices. From https://stackoverflow.com/questions/5888287/running-cor-or-any-variant-over-a-sparse-matrix-in-r
#' @param smat sparse matrix
#'
#' @return
#' @export
#'
#' @examples

sparse.cor <- function(smat){
  n <- nrow(smat)
  cMeans <- colMeans(smat)
  covmat <- (as.matrix(crossprod(smat)) - n * tcrossprod(cMeans))/(n - 1)
  sdvec <- sqrt(diag(covmat))
  cormat <- covmat / tcrossprod(sdvec)
  list(cov = covmat, cor = cormat)
}


# Calc.Cor.Seurat ------------------------------------------------------------------------

Calc.Cor.Seurat <- function(assay = "RNA", slot = "data"
                            , digits = 2, obj = combined.obj, ...) {
  expr.mat <- GetAssayData(slot = slot, assay = assay, object = obj)
}


# plot.Metadata.Cor.Heatmap ------------------------------------------------------------------------
plot.Metadata.Cor.Heatmap <- function(
  columns = c( "nCount_RNA", "nFeature_RNA", "percent.mito", "percent.ribo")
  , cormethod = c('pearson', 'spearman')[1]
  , main =paste( "Metadata", cormethod,"correlations")
  , obj = combined.obj
  , w = 10, h = w
  , ...){
  # library(ggcorrplot)


  meta.data <- obj@meta.data
  columns.found <- intersect(colnames(meta.data), columns)

  corX <- cor(meta.data[ , columns.found], method = cormethod)
  pl <- ggcorrplot::ggcorrplot(corX, hc.order = TRUE, title = main
                               , type = "full", lab = T)
  qqSave(pl, fname = ppp(make.names(main),'pdf'), w = w, h = h)
  pl
}




# plot.Metadata.median.fraction.barplot ------------------------------------------------------------------------
plot.Metadata.median.fraction.barplot <- function(
  columns = c(  "percent.mito", "percent.ribo")
  , suffix =  NULL
  , group.by = GetClusteringRuns(obj = obj)[2]
  , method = c('median', 'mean' )[1]
  , min.thr = 2.5 # At least this many percent in at least 1 cluster
  , return.matrix = F
  , main = paste( method, "read fractions per transcript class and cluster", suffix)
  , ylab = "Fraction of transcriptome (%)"
  , percentify = T
  , subt = NULL
  , position = position_stack()
  , w = 10, h = 6
  , obj = combined.obj
  , ...){

  meta.data <- obj@meta.data
  stopifnot(group.by %in% colnames(meta.data))
  columns.found <- intersect(colnames(meta.data), c(group.by, columns) )

  (mat.cluster.medians1 <- meta.data[ , columns.found] %>%
      group_by_at(group.by) %>%
      dplyr::summarize_all(median)
  )
  if (min.thr>0) {
    pass.cols <- colMax(mat.cluster.medians1[,-1]) > (min.thr/100)
    cols.OK <- which_names(pass.cols)
    cols.FAIL <- which_names(!pass.cols)
    subt = paste(length(cols.FAIL), "classed do not reach", min.thr, "% :", kpps(cols.FAIL))
    iprint(subt)
    mat.cluster.medians1 <- mat.cluster.medians1[ , c( group.by, cols.OK) ]
  }


  mat.cluster.medians <- mat.cluster.medians1 %>%
    reshape2::melt(id.vars = c(group.by), value.name = "Fraction")


  if (percentify)  mat.cluster.medians$'Fraction' = 100*mat.cluster.medians$'Fraction'

  pl <- ggbarplot(mat.cluster.medians, x = group.by, y = 'Fraction', fill = 'variable'
                  , position = position
                  , title = main, subtitle = subt ,ylab = ylab)
  qqSave(pl, fname = ppp(make.names(main),'pdf'), w = w, h = h)
  pl
  if (return.matrix) mat.cluster.medians1 else pl
}

# plot.Metadata.median.fraction.barplot()



# plot.Gene.Cor.Heatmap ------------------------------------------------------------------------
plot.Gene.Cor.Heatmap <- function(genes = WU.2017.139.IEGsf
                                  , assay.use = "RNA", slot.use = c("data", "scale.data", "data.imputed")[1], quantileX = 0.95
                                  , min.g.cor =  0.3, calc.COR = FALSE
                                  , cutRows = NULL, cutCols = cutRows
                                  , obj = combined.obj, ...){
  expr.mat <- GetAssayData(slot = slot.use, assay = assay.use, object = obj)
  if (slot.use == c("data.imputed")) {
    "WIP"
  }
  expr.mat <- GetAssayData(slot = slot.use, assay = assay.use, object = obj)

  qname = p0("expr.q", quantileX * 100)
  slotname_cor.mat <- kpp('cor', slot.use, assay.use, qname)
  cor.mat <- obj@misc[[slotname_cor.mat]]

  if (is_null(cor.mat)) {
    iprint(slotname_cor.mat, " not found in @misc.")
    iprint("Correlation slots present in @misc:",grepv(names(combined.obj@misc), pattern = "^cor"))

    # Calculate ------------------------------------
    if (calc.COR) {
      print("Calculating correlation now.")
      genes.found <- check.genes(genes)
      iprint(l(genes.found), "genes are found in the object.")
      if (l(genes.found) > 200) iprint("Too many genes found in data, cor will be slow: ", l(genes.found))
      ls.cor <- sparse.cor(t(expr.mat[genes.found,]))
      cor.mat <- ls.cor$cor
    } else { stop() }
  } else {
    print("Correlation is pre-calculated")
    genes.found <- intersect(genes, rownames(cor.mat))
    iprint(l(genes.found), "genes are found in the correlation matrix.")
    cor.mat <- cor.mat[genes.found, genes.found]
  }


  # Filter ------------------------------------
  diag(cor.mat) <- NaN
  corgene.names <- union(
    which_names(rowMax(cor.mat) >= min.g.cor),
    which_names(rowMin(cor.mat) <= -min.g.cor)
  )
  iprint(l(corgene.names), "genes are more (anti-)correlated than +/-:", min.g.cor)

  pname = p0("Pearson correlations of ", substitute(genes),"\n min.cor:", min.g.cor, " | ",  assay.use ,'.', slot.use )
  o.heatmap <- pheatmap(cor.mat[corgene.names,corgene.names],main = pname, cutree_rows = cutRows, cutree_cols = cutCols, ...)
  wplot_save_pheatmap(o.heatmap, filename = make.names(pname))

  # return values
  maxCorrz <- rowMax(cor.mat)[corgene.names]; names(maxCorrz) <- corgene.names
  dput(maxCorrz)
}

# Calc.Cor.Seurat ------------------------------------------------
Calc.Cor.Seurat <- function(assay.use = "RNA", slot.use = "data"
                            , quantileX = 0.95, max.cells =  40000, seed = p$"seed"
                            , digits = 2, obj = combined.obj) {
  expr.mat <- GetAssayData(slot = slot.use, assay = assay.use, object = obj)
  if (ncol(expr.mat) > max.cells) {
    set.seed(seed = seed)
    cells.use <- sample(x = colnames(expr.mat), size = max.cells)
  }

  qname = p0("q", quantileX * 100)
  quantile_name = kpp("expr", qname)

  if (is.null(obj@misc[[quantile_name]])) iprint("Call: combined.obj <- calc.q90.Expression.and.set.all.genes(combined.obj, quantileX =",quantileX," first )")
  genes.HE = which_names(obj@misc[[quantile_name]] > 0)
  iprint("Pearson correlation is calculated for", l(genes.HE), "HE genes with expr.",qname,": > 0.")
  tic(); ls.cor <- sparse.cor(smat = t(expr.mat[genes.HE, cells.use])); toc()
  ls.cor <- lapply(ls.cor, round, digits = 2)

  slot__name <- kpp(slot.use, assay.use, quantile_name)
  obj@misc[[kpp('cor', slot__name)]] <- ls.cor$'cor'
  obj@misc[[kpp('cov', slot__name)]] <- ls.cor$'cov'
  iprint("Stored under obj@misc$", kpp('cor', slot.use, assay.use), "or cov... ." )
  return(obj)
}
# combined.obj <- calc.q90.Expression.and.set.all.genes(combined.obj, quantileX = 0.99, max.cells =  400000, set.all.genes = F)
# combined.obj <- Calc.Cor.Seurat(assay.use = "RNA", slot.use = "data", digits = 2, obj = combined.obj, quantile = 0.99, max.cells = 40000)

# plot.clust.size.distr ------------------------------------------------
plot.clust.size.distr <- function(obj = combined.obj, ident = GetClusteringRuns()[2]
                                  , plot = T, thr.hist = 30, ...) {
  clust.size.distr <- table(obj@meta.data[,ident])
  print(clust.size.distr)
  resX <- gsub(pattern = ".*res\\.", replacement = '',x = ident)
  ptitle <- ppp('clust.size.distr', ident)
  psubtitle <- paste("Nr.clusters:", l(clust.size.distr)
                     , "| median:", median(clust.size.distr)
                     , "| CV:", percentage_formatter(cv(clust.size.distr))
  )
  xlb = "Cluster size (cells)"
  xlim = c(0, max(clust.size.distr))

  if (plot) {
    if (l(clust.size.distr) < thr.hist) {
      qbarplot(clust.size.distr, plotname = ptitle, subtitle = psubtitle, xlab = xlb, ...)
    } else {
      qhistogram(vec = clust.size.distr, plotname = ptitle, subtitle = psubtitle, xlab = xlb, xlim = xlim, ...)
    }
  } else {    "return vector"
    clust.size.distr
  }

}
# plot.clust.size.distr()



#  ------------------------------------------------
geneExpressionLevelPlots <- function(gene = 'TOP2A', obj = ls.Seurat[[1]], slot = c('counts', 'data')[2] ) {
  slot = 'data'
  print(gene)
  if (gene %in% rownames(obj)) {
    GEX.Counts <- GetAssayData(object = obj, assay = 'RNA', slot = slot)

    GEX.Counts.total <- rowSums(GEX.Counts)
    genes.expression <- GEX.Counts.total[gene]
    mean.expr <- iround(mean(GEX.Counts[gene,]))

    suffx = if (slot == 'counts') 'raw' else 'normalised, logtransformed'
    (pname = paste(gene, 'and the', suffx,'transcript count distribution'))

    qhistogram(GEX.Counts.total, vline = genes.expression, logX = T, w = 6, h = 4
               , subtitle = paste('It belong to the top', pc_TRUE(GEX.Counts.total > genes.expression), 'of genes (black line). Mean expr:', mean.expr)
               , plotname = pname, xlab = 'Total Transcripts in Dataset', ylab = 'Number of Genes')
  } else { print("     !!! Gene not found in object!")}
}

#  ------------------------------------------------
PrctCellExpringGene <- function(genes, group.by = "all", object = combined.obj){ # From Github/Ryan-Zhu https://github.com/satijalab/seurat/issues/371
  if(group.by == "all"){
    prct = unlist(lapply(genes, ww.calc_helper, object = object))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }

  else{
    list = SplitObject(object, group.by)
    factors = names(list)
    results = lapply(list, PrctCellExpringGene, genes = genes)
    for (i in 1:length(factors)) {
      results[[i]]$Feature = factors[i]
    }
    combined = do.call("rbind", results)
    return(combined)
  }
}


#  ------------------------------------------------
ww.calc_helper <- function(object, genes){ # From Github/Ryan-Zhu https://github.com/satijalab/seurat/issues/371
  counts = object[['RNA']]@counts
  ncells = ncol(counts)
  if (genes %in% row.names(counts)) {
    sum(counts[genes, ] > 0) / ncells
  } else{
    return(NA)
  }
}

#  ------------------------------------------------

# scBarplotFractionAboveThr <- function(thrX = 0.01, columns.used = c('cl.names.top.gene.res.0.3', 'percent.ribo')
#                                       , obj = combined.obj, ) { # Calculat the fraction of cells per cluster above a certain threhold
#   meta = obj@meta.data
#   (fr_ribo_low_cells <- meta %>%
#       dplyr::select(columns.used)  %>%
#       dplyr::group_by_(columns.used[1])  %>%
#       summarize(n_cells = n(),
#                 n_ribo_low_cells = sum(!!as.name(columns.used[2]) < thrX),
#                 fr_ribo_low_cells = n_ribo_low_cells / n_cells) %>%
#       FirstCol2RowNames())
#   print(fr_ribo_low_cells)
#
#   (v.fr_ribo_low_cells <- 100* as.named.vector(fr_ribo_low_cells[3]))
#   qbarplot(v.fr_ribo_low_cells, xlab.angle = 45, xlab = 'Clusters', ylab = '% Cells')
# }




scBarplotFractionAboveThr <- function(thrX = 0., value.col = 'TVA', id.col =  'cl.names.top.gene.res.0.3'
                                      , obj = combined.obj, return.df = F) { # Calculat the fraction of cells per cluster above a certain threhold
  meta = obj@meta.data
  (df_cells_above <- meta %>%
      dplyr::select(c(id.col, value.col))  %>%
      dplyr::group_by_(id.col)  %>%
      summarize(n_cells = n(),
                n_cells_above = sum(!!as.name(value.col) > thrX),
                fr_n_cells_above = n_cells_above / n_cells) %>%
      FirstCol2RowNames())

  (v.fr_n_cells_above <- 100* as.named.vector(df_cells_above[3]))
  ggobj <- qbarplot(v.fr_n_cells_above, xlab = 'Clusters', ylab = '% Cells'
                    , plotname = paste('Cells with', value.col, '>', thrX)
                    , subtitle = id.col, xlab.angle = 45)
  if (return.df) return(df_cells_above) else ggobj
}

# combined.obj$'TVA' = combined.obj@assays$RNA@data['AAV.TVA.p2a.G.p2a.GFP.WPRE.bGH',]
# scBarplotFractionAboveThr(id.col =  'cl.names.top.gene.res.0.3', value.col = 'TVA', thrX = 0)


#  ------------------------------------------------
scBarplotFractionBelowThr <- function(thrX = 0.01, value.col = 'percent.ribo', id.col =  'cl.names.top.gene.res.0.3'
                                      , obj = combined.obj, return.df = F) { # Calculat the fraction of cells per cluster below a certain threhold
  meta = obj@meta.data
  (df_cells_below <- meta %>%
      dplyr::select(c(id.col, value.col))  %>%
      dplyr::group_by_(id.col)  %>%
      summarize(n_cells = n(),
                n_cells_below = sum(!!as.name(value.col) < thrX),
                fr_n_cells_below = n_cells_below / n_cells) %>%
      FirstCol2RowNames())

  (v.fr_n_cells_below <- 100* as.named.vector(df_cells_below[3]))
  ggobj <- qbarplot(v.fr_n_cells_below, xlab = 'Clusters', ylab = '% Cells'
                    , plotname = make.names(paste('Cells with', value.col, '<', thrX))
                    , subtitle = id.col, xlab.angle = 45)
  if (return.df) return(df_cells_below) else ggobj
}

# scBarplotFractionBelowThr(id.col =  'cl.names.top.gene.res.0.3', value.col = 'percent.ribo', thrX = 0.01, return.df = T)



#  ------------------------------------------------

#  ------------------------------------------------

######################################################################
# Read.Write.Save.Load.functions.R
######################################################################
# source('~/GitHub/Packages/Seurat.utils/Functions/Read.Write.Save.Load.functions.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Read.Write.Save.Load.functions.R"))

"Multicore read / write (I/O) functions are https://github.com/vertesy/Seurat.multicore"
"Single core read / write (I/O) functions are in https://github.com/vertesy/Seurat.utils/"


# Convert10Xfolders ------------------------------------------------------------------------
Convert10Xfolders <- function(InputDir # Take a parent directory with a number of subfolders, each containing the standard output of 10X Cell Ranger. (1.) It loads the filtered data matrices; (2.) converts them to Seurat objects, and (3.) saves them as *.RDS files.
                              , regex = F, folderPattern = c("filtered_feature", "SoupX_decont")[1]
                              , min.cells = 5, min.features = 200
                              , updateHGNC = T, ShowStats = T) {

  # finOrig <- list.dirs(InputDir, recursive = subdirs)
  finOrig <- list.dirs.depth.n(InputDir, depth = 2)
  fin <- grepv(x = finOrig, pattern = folderPattern, perl = regex)

  iprint(length(fin), "samples found.")
  if (l(fin)) {
    for (i in 1:length(fin)) { print(i)
      pathIN = fin[i]; print(pathIN)
      # fnameIN = basename(dirname(xx))
      fnameIN = strsplit(basename(dirname(fin[i])),split = "_")[[1]][1]
      print(fnameIN)
      fnameOUT = ppp(paste0(InputDir, '/', fnameIN), 'min.cells', min.cells, 'min.features', min.features,"Rds")
      print(fnameOUT)
      count_matrix <- Read10X(pathIN)

      if ( !is.list(count_matrix) | length(count_matrix) == 1) {
        seu <- CreateSeuratObject(counts = count_matrix, project = fnameIN,
                                  min.cells = min.cells, min.features = min.features)
      } else if (is.list(count_matrix) & length(count_matrix) == 2)  {
        seu <- CreateSeuratObject(counts = count_matrix[[1]], project = fnameIN,
                                  min.cells = min.cells, min.features = min.features)

        # LSB, Lipid Sample barcode (Multi-seq) --------------------
        LSB <- CreateSeuratObject(counts = count_matrix[[2]], project = fnameIN)
        LSBnameOUT = ppp(paste0(InputDir, '/LSB.', fnameIN),"Rds")
        saveRDS(LSB, file = LSBnameOUT)
      } else {
        print('More than 2 elements in the list of matrices')
      }
      # update----
      if (updateHGNC) seu <- UpdateGenesSeurat(seu, EnforceUnique = T, ShowStats = T)
      saveRDS(seu, file = fnameOUT)
    }
  } else { iprint("No subfolders found with pattern", folderPattern, "in dirs like: ", finOrig[1:3]) }
}
# Convert10Xfolders(InputDir)


# Convert10Xfolders.old ------------------------------------------------------------------------
Convert10Xfolders.old <- function(InputDir # Take a parent directory with a number of subfolders, each containing the standard output of 10X Cell Ranger. (1.) It loads the filtered data matrices; (2.) converts them to Seurat objects, and (3.) saves them as *.RDS files.
                                  , folderPattern = c("filtered", "SoupX_decont")[1]
                                  , min.cells=10, min.features=200, updateHGNC=T, ShowStats=T) {
  fin <- list.dirs(InputDir, recursive = F)
  fin <- grepv(x = fin, pattern = folderPattern, perl = F)

  for (i in 1:length(fin)) {
    pathIN = fin[i]; print(pathIN)
    fnameIN = basename(fin[i])
    fnameOUT = ppp(paste0(InputDir, '/', fnameIN), 'min.cells', min.cells, 'min.features', min.features,"Rds")
    count_matrix <- Read10X(pathIN)

    if ( !is.list(count_matrix) | length(count_matrix) == 1) {
      seu <- CreateSeuratObject(counts = count_matrix, project = fnameIN,
                                min.cells = min.cells, min.features = min.features)
    } else if (is.list(count_matrix) & length(count_matrix) == 2)  {
      seu <- CreateSeuratObject(counts = count_matrix[[1]], project = fnameIN,
                                min.cells = min.cells, min.features = min.features)

      # LSB, Lipid Sample barcode (Multi-seq) --------------------
      LSB <- CreateSeuratObject(counts = count_matrix[[2]], project = fnameIN)
      LSBnameOUT = ppp(paste0(InputDir, '/LSB.', fnameIN),"Rds")
      saveRDS(LSB, file = LSBnameOUT)
    } else {
      print('More than 2 elements in the list of matrices')
    }
    # update----
    if (updateHGNC) seu <- UpdateGenesSeurat(seu, EnforceUnique = T, ShowStats = T)
    saveRDS(seu, file = fnameOUT)
  }
}
# Convert10Xfolders(InputDir = InputDir)

# ConvertDropSeqfolders ------------------------------------------------------------------------
ConvertDropSeqfolders <- function(InputDir # Take a parent directory with a number of subfolders, each containing the standard output of 10X Cell Ranger. (1.) It loads the filtered data matrices; (2.) converts them to Seurat objects, and (3.) saves them as *.RDS files.
                                  , folderPattern = "SRR*", filePattern = "expression.tsv.gz"
                                  , useVroom = T, col_types.vroom = list("GENE" = "c", .default = "d")
                                  , min.cells=10, min.features=200, updateHGNC=T, ShowStats=T, minDimension = 10, overwrite = FALSE) {
  InputDir <- FixPath(InputDir)
  fin <- list.dirs(InputDir, recursive = F)
  fin <- grepv(x = fin, pattern = folderPattern, perl = F)

  for (i in 1:length(fin)) { print(i)
    pathIN <- FixPath(fin[i]); print(pathIN)
    fnameIN <- basename(fin[i])
    subdir <- p0(InputDir, fnameIN)
    fnameOUT <- ppp(subdir, 'min.cells', min.cells, 'min.features', min.features,"Rds"); print(fnameOUT)
    if (!overwrite) {
      OutFile <- list.files(InputDir, pattern = basename(fnameOUT), recursive = T)
      if (length(OutFile) > 0) {
        if (grepl(pattern = ".Rds$", OutFile, perl = T)) {
          iprint("      RDS OBJECT ALREADY EXISTS.");
          next
        }
      } # if length
    }
    CountTable <- list.files(subdir, pattern = filePattern,recursive = F)
    stopifnot(length(CountTable) == 1 )
    count_matrix <- if (useVroom) {
      vroom::vroom(file = kpps(subdir, CountTable), col_types = col_types.vroom)
    } else {
      readr::read_tsv(file = kpps(subdir, CountTable))
    }

    if (nrow(count_matrix) < minDimension | ncol(count_matrix) < minDimension ) {
      iprint(""); iprint("      EXPRESSION MATRIX TOO SMALL.", nrow(count_matrix), "x", ncol(count_matrix),". Not processed.");
    } else {
      count_matrix <- FirstCol2RowNames(count_matrix)[,-1] # remove 1st "Cell column" # https://github.com/vertesy/SEO/issues/63
      seu <- CreateSeuratObject(counts = count_matrix, project = fnameIN,
                                min.cells = min.cells, min.features = min.features)
      if (ncol(seu) < 1000) print("Only", ncol(seu), "cells survived filtering in the Seurat obj!")
      if (nrow(seu) < 1000) print("Only", nrow(seu), "genes found in the Seurat obj!")

      # update HGNC ----
      Sys.setenv('R_MAX_VSIZE' = 32000000000)
      if (updateHGNC) seu <- UpdateGenesSeurat(seu, EnforceUnique = T, ShowStats = T)
      saveRDS(seu, file = fnameOUT)
    }
  }
}
# ConvertDropSeqfolders(InputDir)

# LoadAllSeurats ------------------------------------------------------------------------
LoadAllSeurats <- function(InputDir # Load all Seurat objects found in a directory. Also works with symbolic links (but not with aliases).
                           , file.pattern = "^filtered.+Rds$"
                           , string.remove1 = c(F, "filtered_feature_bc_matrix.", "raw_feature_bc_matrix." )[2]
                           , string.remove2 = c(F, ".min.cells.10.min.features.200.Rds")[2]) {
  tic()
  InputDir <- AddTrailingSlash(InputDir) # add '/' if necessary

  fin.orig <- list.files(InputDir, include.dirs = F, pattern = file.pattern)
  print(fin.orig)
  fin <- if (!isFALSE(string.remove1)) sapply(fin.orig, gsub, pattern = string.remove1, replacement = "") else fin.orig
  fin <- if (!isFALSE(string.remove2)) sapply(fin, gsub, pattern = string.remove2, replacement = "") else fin

  ls.Seu <- list.fromNames(fin)
  for (i in 1:length(fin)) {print(fin[i]); ls.Seu[[i]] <- readRDS(paste0(InputDir, fin.orig[i]))}
  print(toc())
  return(ls.Seu)
}
# ls.Seurat <- LoadAllSeurats(InputDir)


# ------------------------------------------------------------------------------------------------
read10x <- function(dir) { # read10x from gzipped matrix.mtx, features.tsv and barcodes.tsv
  tictoc::tic()
  names <- c("barcodes.tsv", "features.tsv", "matrix.mtx")
  for (i in 1:length(names)) {
    R.utils::gunzip(paste0(dir, "/", names[i], ".gz"))
  }
  file.copy(paste0(dir, "/features.tsv"), paste0(dir, "/genes.tsv"))
  mat <- Seurat::Read10X(dir)
  file.remove(paste0(dir, "/genes.tsv"))
  for (i in 1:length(names)) {
    R.utils::gzip(paste0(dir, "/", names[i]))
  }
  tictoc::toc()
  mat
}

#### Functions in Saving.and.loading.R


# saveRDS.compress.in.BG ------------------------------------------------------------------------
saveRDS.compress.in.BG <- function(obj, compr = FALSE, fname) {
  try(tictoc::tic(), silent = T)
  saveRDS(object = obj, compress = compr, file = fname)
  try(tictoc::toc(), silent = T)
  print(paste("Saved, being compressed", fname))
  system(paste("gzip", fname),  wait = FALSE) # execute in the background
  try(say(), silent = T)
}

# Save an object -----------------------------------------------
isave.RDS <- function(object, prefix =NULL, suffix=NULL, inOutDir = F
                      , alternative_path_rdata = paste0("~/Dropbox/Abel.IMBA/AnalysisD/_RDS.files/", basename(OutDir))
                      , showMemObject=T, saveParams =T){ # Faster saving of workspace, and compression outside R, when it can run in the background. Seemingly quite CPU hungry and not very efficient compression.
  path_rdata = if (inOutDir) OutDir else alternative_path_rdata
  dir.create(path_rdata)

  if (showMemObject) { memory.biggest.objects() }
  if ( "seurat" %in% is(object) & saveParams) {
    try(object@misc$p <- p, silent = T)
    try(object@misc$all.genes  <- all.genes, silent = T)
  }
  fnameBase = kppu(prefix, substitute(object), suffix, idate(Format = "%Y.%m.%d_%H.%M"))
  fnameBase = trimws(fnameBase, whitespace = '_')
  saveRDS.compress.in.BG(obj = object, fname = paste0(path_rdata, "/",fnameBase , ".Rds") )
}



# downsampleSeuObj.and.Save ------------------------------------------------------------------------
downsampleSeuObj.and.Save <- function(obj=ORC, fraction = 0.25, seed = 1989, dir = OutDir, min.features = p$'min.features', suffix = '') { # Subset a compressed Seurat Obj and save it in wd.
  obj_Xpc <- downsampleSeuObj(obj = obj, fraction_ =  fraction, seed_ = seed)
  nr.cells.kept <- ncol(obj_Xpc)
  saveRDS.compress.in.BG(obj = obj_Xpc, fname = ppp(paste0(dir, substitute(obj)),suffix, nr.cells.kept, 'cells.with.min.features', min.features,"Rds" ) )
}


# Save workspace -----------------------------------------------
# requires MarkdownReportsDev (github) and defining OutDir
# requires github/vertesy/CodeAndRoll.r

isave.image <- function(..., path_rdata = paste0("~/Dropbox/Abel.IMBA/AnalysisD/_Rdata.files/", basename(OutDir))
                        , showMemObject=T, options=c("--force", NULL)[1]
){ # Faster saving of workspace, and compression outside R, when it can run in the background. Seemingly quite CPU hungry and not very efficient compression.

  dir.create(path_rdata)

  if (showMemObject) { memory.biggest.objects() }
  fname = MarkdownReportsDev::kollapse(path_rdata, "/",idate(),...,".Rdata")
  print(fname)
  if (nchar(fname) > 2000) stop()
  save.image(file = fname, compress = F)
  MarkdownReportsDev::iprint("Saved, being compressed", fname)
  system(paste("gzip", options, fname),  wait = FALSE) # execute in the background
}


# Save workspace -----------------------------------------------
# requires MarkdownReportsDev (github) and defining OutDir
# requires github/vertesy/CodeAndRoll.r

qsave.image <- function(..., showMemObject=T, options=c("--force", NULL)[1]){ # Faster saving of workspace, and compression outside R, when it can run in the background. Seemingly quite CPU hungry and not very efficient compression.
  fname = MarkdownReportsDev::kollapse(getwd(), "/",basename(OutDir),idate(),...,".Rdata")
  print(fname)
  if (nchar(fname) > 2000) stop()
  tic()
  save.image(file = fname, compress = F)
  MarkdownReportsDev::iprint("Saved, being compressed", fname)
  system(paste("gzip", options, fname),  wait = FALSE) # execute in the background
  cat(toc)
}


# downsampleSeuObj -----------------------------------------------------------------------
downsampleSeuObj <- function(obj=ls.Seurat[[i]], fraction_ = 0.25, nCells = F, seed_ = 1989 ) { # Subset a compressed Seurat Obj and save it in wd.
  set.seed(seed_)
  if (isFALSE(nCells)) {
    cellIDs.keep = sampleNpc(metaDF = obj@meta.data, pc = fraction_)
    iprint(length(cellIDs.keep), "or",percentage_formatter(fraction_),"of the cells are kept. Seed:", seed_)
  } else if (nCells > 1) {
    nKeep = min(ncol(obj), nCells)
    # print(nKeep)
    cellIDs.keep = sample(colnames(obj), size = nKeep, replace = F)
    if (nKeep < nCells) iprint("Only",nCells,"cells were found in the object, so downsampling is not possible.")
  }
  obj <- subset(x = obj, cells = cellIDs.keep) # downsample
  return(obj)
}

# downsampleSeuObj.and.Save ------------------------------------------------------------------------
downsampleSeuObj.and.Save <- function(obj=ORC, fraction = 0.25, seed = 1989, dir = OutDir, suffix = '') { # Subset a compressed Seurat Obj and save it in wd.
  obj_Xpc <- downsampleSeuObj(obj = obj, fraction_ =  fraction, seed_ = seed)
  nr.cells.kept <- ncol(obj_Xpc)
  saveRDS.compress.in.BG(obj = obj_Xpc, fname = ppp(paste0(dir, substitute(obj)),suffix, nr.cells.kept, 'cells.with.min.features', p$min.features,"Rds" ) )
}


# Downsample.Seurat.Objects ------------------------------------------------------------------------
Downsample.Seurat.Objects <- function(ls.obj = ls.Seurat, NrCells = p$"dSample.Organoids") {
  names.ls = names(ls.obj)
  n.datasets = length(ls.obj)
  iprint(NrCells, "cells")
  tic()
  if (getDoParRegistered() ) {
    ls.obj.downsampled <- foreach(i = 1:n.datasets ) %dopar% {
      iprint(names(ls.obj)[i], percentage_formatter(i/n.datasets, digitz = 2))
      downsampleSeuObj(obj = ls.obj[[i]], nCells = NrCells)
    }; names(ls.obj.downsampled)  <- names.ls
  } else {
    ls.obj.downsampled <- list.fromNames(names.ls)
    for (i in 1:n.datasets ) {
      iprint(names(ls.obj)[i], percentage_formatter(i/n.datasets, digitz = 2))
      ls.obj.downsampled[[i]] <- downsampleSeuObj(obj = ls.obj[[i]], nCells = NrCells)
    };
  } # else
  toc();

  print(head(unlapply(ls.obj, ncol)))
  print(head(unlapply(ls.obj.downsampled, ncol)))

  isave.RDS(object = ls.obj.downsampled, suffix = ppp(NrCells, "cells"), inOutDir = T)

}
# Downsample.Seurat.Objects(NrCells = 2000)
# Downsample.Seurat.Objects(NrCells = 1000)
# Downsample.Seurat.Objects(NrCells = 500)
# Downsample.Seurat.Objects(NrCells = 200)


######################################################################
# Seurat.object.manipulations.etc.R
######################################################################
# source('~/GitHub/Packages/Seurat.utils/Functions/Seurat.object.manipulations.etc.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Seurat.object.manipulations.etc.R"))

# ------------------------------------------------------------------------
clip10Xcellname <- function(cellnames) str_split_fixed(cellnames, "_", n = 2)[,1] # Clip all suffices after underscore (10X adds it per chip-lane, Seurat adds in during integration).

# ------------------------------------------------------------------------
make10Xcellname <- function(cellnames, suffix="_1") paste0(cellnames, suffix) # Add a suffix to cell names, so that it mimics the lane-suffix, e.g.: "_1".


# ------------------------------------------------------------------------
seu.Make.Cl.Label.per.cell <- function(TopGenes, clID.per.cell) { # Take a named vector (of e.g. values ="gene names", names = clusterID), and a vector of cell-IDs and make a vector of "GeneName.ClusterID".
  Cl.names_class= TopGenes[ clID.per.cell ]
  Cl.names_wNr = paste0(Cl.names_class,' (',names(Cl.names_class),')')
  return(Cl.names_wNr)
}
# seu.Make.Cl.Label.per.cell(TopGenes = TopGenes.Classic,
#                            clID.per.cell = getMetadataColumn(ColName.metadata = metaD.CL.colname)  )

# FeaturePlot with different defaults ------------------------------------------------------------------
GetMostVarGenes <- function(obj=org, nGenes = p$nVarGenes) { # Get the most variable rGenes
  head(rownames(slot(object = obj, name = "hvg.info")), n = nGenes)
}

# gene.name.check for read .mtx /write .rds script ---------------------------------------
gene.name.check <- function(Seu.obj = ls.Seurat[[1]] ) { # Check gene names in a seurat object, for naming conventions (e.g.: mitochondrial reads have - or .). Use for reading .mtx & writing .rds files.
  rn = rownames(GetAssayData(object = Seu.obj, slot = "counts"))
  llprint("### Gene name pattern")

  llogit('`rn = rownames(GetAssayData(object = ls.Seurat[[1]], slot = "counts"))`')
  llogit('`head(grepv(rn, pattern = "-"), 10)`')
  print('pattern = -')
  llprint(head(grepv(rn, pattern = "-"), 10))

  llogit('`head(grepv(rn, pattern = "_"), 10)`')
  print('pattern = _')
  llprint(head(grepv(rn, pattern = "_"), 10))

  llogit('`head(grepv(rn, pattern = "\\."), 10)`')
  print('pattern = \\.')
  llprint(head(grepv(rn, pattern = "\\."), 10))

  llogit('`head(grepv(rn, pattern = "\\.AS[1-9]"), 10)`')
  print('pattern = \\.AS[1-9]')
  llprint(head(grepv(rn, pattern = "\\.AS[1-9]"), 10))
}


# check.genes ---------------------------------------
check.genes <- function(list.of.genes = ClassicMarkers, makeuppercase = FALSE, verbose  = TRUE, HGNC.lookup = FALSE
                        , obj = combined.obj
                        , assay.slot=c('RNA', 'integrated')[1]
                        , dataslot = c("counts", "data")[2]) { # Check if genes exist in your dataset.
  if (makeuppercase) list.of.genes <- toupper(list.of.genes)
  all_genes = rownames(GetAssayData(object = obj, assay = assay.slot, slot = dataslot)); length(all_genes)
  missingGenes = setdiff(list.of.genes, all_genes)
  if (length(missingGenes) > 0) {
    if (verbose) { iprint(l(missingGenes), "or", percentage_formatter(l(missingGenes) / l(list.of.genes)), "genes not found in the data, e.g:", head(missingGenes, n = 10))  }
    if (HGNC.lookup) {
      if (exists('qHGNC', mode='function')) { try(qHGNC(missingGenes)) } else { print("load qHGNC() function, see database.linker")}
    }
  }
  intersect(list.of.genes, all_genes)
}
# check.genes("top2a", makeuppercase = TRUE)
# check.genes("VGLUT2", verbose = F, HGNC.lookup = T)


# replace zero indexed clusternames ------------------------------------------------
fixZeroIndexing.seurat <- function(ColName.metadata = 'res.0.6', obj=org) { # Fix zero indexing in seurat clustering, to 1-based indexing
  obj@meta.data[ ,ColName.metadata] =  as.numeric(obj@meta.data[ ,ColName.metadata])+1
  print(obj@meta.data[ ,ColName.metadata])
  return(obj)
}


# CalculateFractionInTrome ------------------------------------------------
CalculateFractionInTrome <- function(geneset = c("MALAT1") # Calculate the fraction of a set of genes within the full Transcriptome of each cell.
                                     , obj = combined.obj
                                     , dataslot = c("counts", "data")[2]
) {
  print("    >>>> Use addMetaFraction() <<<<")
  geneset <- check.genes(geneset)
  stopifnot(length(geneset)>0)

  mat <- as.matrix(slot(obj@assays$RNA, name=dataslot))
  mat.sub <- mat[geneset,,drop=F]
  RC.per.cell.geneset <- colSums(mat.sub)

  RC.per.cell <- colSums(mat)
  gene.fraction.per.cell <- 100*RC.per.cell.geneset / RC.per.cell
  return(gene.fraction.per.cell)
}

# ------------------------------------------------------------------------
AddNewAnnotation <- function(obj = obj # Create a new metadata column based on an exisiting metadata column and a list of mappings (name <- IDs).
                             , source = "RNA_snn_res.0.5", named.list.of.identities = ls.Subset.ClusterLists) {
  NewID <- as.named.vector(obj[[source]])

  for (i in 1:length(named.list.of.identities)) {
    lx <- as.character(named.list.of.identities[[i]])
    name.lx <- names(named.list.of.identities)[i]
    NewID <- translate(vec = NewID, old = lx, new = name.lx)
  }
  print(table(NewID))
  return(NewID)
}
# ls.Subset.ClusterLists = list( "hESC.h9" = c("4", "10", "14"), "hESC.176" = c("0", "1", "2")); AddNewAnnotation()

# whitelist.subset.ls.Seurat ------------------------------------------------------------------------
whitelist.subset.ls.Seurat <- function(ls.obj = ls.Seurat
                                       , metadir = p$'cellWhiteList' #  '~/Dropbox/Abel.IMBA/MetadataD/POL.meta/cell.lists/'
                                       , whitelist.file = "NonStressedCellIDs.2020.10.21_18h.tsv"
) {
  cells.before <- unlapply(ls.obj, ncol)
  # Find file
  df.cell.whitelist <- read.simple.tsv(metadir, whitelist.file)
  dsets <- table(df.cell.whitelist[,1])

  ls.orig.idents <- lapply(lapply(ls.Seurat, getMetadataColumn, ColName.metadata = "orig.ident"), unique)
  stopif(any(unlapply(ls.orig.idents, l) == l(ls.Seurat)), message = "Some ls.Seurat objects have 1+ orig identity.")

  dsets.in.lsSeu <- unlist(ls.orig.idents)
  isMathced <- all(dsets.in.lsSeu == names(dsets)) # Stop if either ls.Seurat OR the metadata has identities not found in the other, in the same order.
  stopif(!isMathced, message = paste("either ls.Seurat OR the metadata has identities not found in the other, or they are not in same order."
                                     , kpps(dsets.in.lsSeu),"vs.", kpps(names(dsets) ) )
  )

  # identX <- ls.orig.idents[[1]]
  for (i in 1:l(ls.orig.idents)) {
    identX <- ls.orig.idents[[i]]; print(identX)

    # Extract and process cellIDs ----
    idx.match <- which(df.cell.whitelist[,1] == identX)
    cell.whitelist <- rownames(df.cell.whitelist)[idx.match]
    cell.whitelist <- substr(x = cell.whitelist
                             , start = 1 ,stop = nchar(cell.whitelist)-2)

    # Extract and process cellIDs ----
    ls.obj[[i]] <- subset(x = ls.obj[[i]], cells = cell.whitelist)
  }
  cells.after <- unlapply(ls.obj, ncol)
  iprint("cells.before",cells.before,"cells.after",cells.after)
  return(ls.obj)
}

# FindCorrelatedGenes ------------------------------------------------------------------------
FindCorrelatedGenes <- function(gene ="N.RabV.N2c", obj = combined.obj, assay = "RNA", slot = "data"
                                , HEonly =F , minExpr = 1, minCells = 1000
                                , trailingNgenes = 1000) {
  tic()
  AssayData <- GetAssayData(object = obj, assay = assay, slot = slot )
  matrix_mod <- iround(as.matrix(AssayData))
  if (HEonly) {
    idx.pass <- (matrixStats::rowSums2(matrix_mod > minExpr) > minCells)
    pc_TRUE(idx.pass)
    matrix_mod <- matrix_mod[ which(idx.pass), ]
  }
  geneExpr <- as.numeric(matrix_mod[gene, ])
  correlations <- apply(matrix_mod, 1, cor, geneExpr)
  topGenes <- trail(sort(correlations, decreasing = T), N = trailingNgenes)
  toc()
  wbarplot(head(topGenes, n =25))
  topGenes
}
# FindCorrelatedGenes(gene ="N.RabV.N2c", obj = combined.obj)
# write_clip(names(head(topGenes[-(1:6)], n=50)))


# ------------------------------------------------------------------------
Calc.Cor.Seurat <- function(assay.use = "RNA", slot.use = "data", geneset = FALSE
                            , quantileX = 0.95, max.cells =  10000, seed = p$"seed"
                            , digits = 2, obj = combined.obj) {
  expr.mat <- GetAssayData(slot = slot.use, assay = assay.use, object = obj)
  if (ncol(expr.mat) > max.cells) {
    set.seed(seed = seed)
    cells.use <- sample(x = colnames(expr.mat), size = max.cells)
  } else {
    cells.use <- colnames(obj)
  }

  qname = p0("q", quantileX * 100)
  quantile_name = kpp("expr", qname)
  if (is.null(obj@misc[[quantile_name]])) { iprint("Quantile data missing! Call: combined.obj <- calc.q90.Expression.and.set.all.genes(combined.obj, quantileX =",quantileX,") first!"); stop()}

  genes.HE  <- if (isFALSE(geneset)) {  which_names(obj@misc[[quantile_name]] > 0) } else {
    check.genes(geneset)  }
  iprint("Pearson correlation is calculated for", l(genes.HE), "HE genes with expr."
         , qname,": > 0 on a sample of", max.cells, " cells.")
  tic(); ls.cor <- sparse.cor(smat = t(expr.mat[genes.HE, cells.use])); toc()

  ls.cor <- lapply(ls.cor, round, digits = 2)

  slot__name <- kpp(slot.use, assay.use, quantile_name)
  obj@misc[[kpp('cor', slot__name)]] <- ls.cor$'cor'
  obj@misc[[kpp('cov', slot__name)]] <- ls.cor$'cov'
  iprint("Stored under obj@misc$", kpp('cor', slot.use, assay.use), "or cov... ." )
  return(obj)
}



# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------

######################################################################
# Seurat.update.gene.symbols.HGNC.R
######################################################################
# source('~/GitHub/Packages/Seurat.utils/Functions/Seurat.update.gene.symbols.HGNC.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Seurat.update.gene.symbols.HGNC.R"))
# require(HGNChelper)

# updateHGNC ------------------------------------------------------------------------------------
UpdateGenesSeurat <- function(obj = ls.Seurat[[i]], species_="human", EnforceUnique = T, ShowStats=F ) { # Update genes symbols that are stored in a Seurat object. It returns a data frame. The last column are the updated gene names.
  HGNC.updated <- HGNChelper::checkGeneSymbols(rownames(obj), unmapped.as.na = FALSE, map = NULL, species = species_)
  if (EnforceUnique) HGNC.updated <- HGNC.EnforceUnique(HGNC.updated)
  if (ShowStats) print(GetUpdateStats(HGNC.updated))
  obj <- RenameGenesSeurat(obj, newnames = HGNC.updated$Suggested.Symbol)
  return(obj)
}
# UpdateGenesSeurat()

# HELPER RenameGenesSeurat  ------------------------------------------------------------------------------------
RenameGenesSeurat <- function(obj = ls.Seurat[[i]], newnames = HGNC.updated[[i]]$Suggested.Symbol) { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  RNA <- obj@assays$RNA

  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
    if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- newnames
    # if (length(obj@meta.data)) rownames(obj@meta.data)          <- newnames
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}
# RenameGenesSeurat(obj = SeuratObj, newnames = HGNC.updated.genes$Suggested.Symbol)


# RemoveGenesSeurat ------------------------------------------------------------------------------------
RemoveGenesSeurat <- function(obj = ls.Seurat[[i]], symbols2remove = c("TOP2A")) { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes metadata; obj@assays$RNA@counts, @data and @scale.data.
  print("Run this as the first thing after creating the Seurat object. It only removes genes from: metadata; obj@assays$RNA@counts, @data and @scale.data.")
  RNA <- obj@assays$RNA

  if (length(RNA@counts)) {
    NotFound <- setdiff(symbols2remove, RNA@counts@Dimnames[[1]])
    if (length(NotFound) == 0)  {
      RNA@counts@Dimnames[[1]] <- symbols2remove
      print("Genes removed from RNA@counts")
    } else {print("Not All Genes Found in RNA@counts. Missing:"); print(NotFound)}
  }
  if (length(RNA@data)) {
    if (length(setdiff(symbols2remove, RNA@data@Dimnames[[1]]) ) == 0)  {
      RNA@data@Dimnames[[1]] <- symbols2remove
      print("Genes removed from RNA@data.")
    } else {print("Not All Genes Found in RNA@data")}
  }
  if (length(RNA@scale.data)) {
    if (length(setdiff(symbols2remove, RNA@scale.data@Dimnames[[1]]) ) == 0)  {
      RNA@scale.data@Dimnames[[1]] <- symbols2remove
      print("Genes removed from RNA@scale.data.")
    } else {print("Not All Genes Found in RNA@scale.data")}
  }
  if (length(obj@meta.data)) {
    if (length(setdiff(symbols2remove, rownames(obj@meta.data)) ) == 0)  {
      rownames(obj@meta.data) <- symbols2remove
      print("Genes removed from @meta.data.")
    } else {print("Not All Genes Found in @metadata")}
  }
  obj@assays$RNA <- RNA
  return(obj)
}
# RemoveGenesSeurat(obj = SeuratObj, symbols2remove = "TOP2A")


# HELPER Enforce Unique names ------------------------------------------------------------------------------------
HGNC.EnforceUnique <- function(updatedSymbols) { # Enforce Unique names after HGNC symbol update. updatedSymbols is the output of HGNChelper::checkGeneSymbols.
  NGL <- updatedSymbols[,3]
  if (any.duplicated(NGL)) {
    updatedSymbols[,3] <- make.unique(NGL); "Unique names are enforced by suffixing .1, .2, etc."
  }
  return(updatedSymbols)
}
# x <- HGNC.EnforceUnique(updatedSymbols = SymUpd)
# While "make.unique" is not the ideal solution, because it generates mismatched, in my integration example it does reduce the mismatching genes from ~800 to 4


# update stats HGNC  ------------------------------------------------------------------------------------
GetUpdateStats <- function(genes = HGNC.updated[[i]]) { # Plot the Symbol-update statistics. Works on the data frame returned by `UpdateGenesSeurat()`.
  (MarkedAsUpdated <- genes[genes$Approved == FALSE, ])
  (AcutallyUpdated <- sum(MarkedAsUpdated[,1] != MarkedAsUpdated[,3]))
  (UpdateStats = c("Updated (%)"=percentage_formatter(AcutallyUpdated / nrow(genes)), "Updated Genes"=floor(AcutallyUpdated), "Total Genes"=floor(nrow(genes))))
  return(UpdateStats)
}
# GetUpdateStats(genes = HGNC.updated.genes)

# update stats HGNC plot ------------------------------------------------------------------------------------
PlotUpdateStats <- function(mat = UpdateStatMat, column.names = c("Updated (%)",  "Updated (Nr.)")) { # Scatter plot of update stats.
  stopifnot(column.names %in% colnames(UpdateStatMat))
  HGNC.UpdateStatistics <- mat[, column.names]
  HGNC.UpdateStatistics[, "Updated (%)"] <- 100*HGNC.UpdateStatistics[, "Updated (%)"]
  colnames(HGNC.UpdateStatistics) <-  c("Gene Symbols updated (% of Total Genes)",  "Number of Gene Symbols updated")
  lll <- wcolorize(vector = rownames(HGNC.UpdateStatistics))
  wplot(HGNC.UpdateStatistics, col = lll
        , xlim = c(0,max(HGNC.UpdateStatistics[,1]))
        , ylim = c(0,max(HGNC.UpdateStatistics[,2])) )
  wlegend(NamedColorVec = lll, poz = 1)
}
# PlotUpdateStats(mat = result.of.GetUpdateStats)

######################################################################
# Soup.Analysis.of.ambient.RNA.R
######################################################################
# source('~/GitHub/Packages/Seurat.utils/Functions/Soup.Analysis.of.ambient.RNA.R')
# try (source('https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Soup.Analysis.of.ambient.RNA.R'))
# Source: self + web

# Requirements ------------------------
# require(tibble)

# plotTheSoup ------------------------------------------------------------------------
plotTheSoup <- function(CellRangerOutputDir = "~/Data/114593/114593"
                        , SeqRun = gsub('*([0-9]+).*','\\1', x = basename(CellRangerOutputDir))) { # Plot the ambient RNA content of droplets without a cell (background droplets).

  ls.Alpha=1
  # Setup ------------------------
  # require(Matrix); require(ggrepel)

  dirz <- list.dirs(CellRangerOutputDir, full.names = F, recursive = F)
  path.raw <- file.path(CellRangerOutputDir, grep(x = dirz, pattern = "^raw_*", value = T))
  path.filt <- file.path(CellRangerOutputDir, grep(x = dirz, pattern = "^filt_*", value = T))
  CR.matrices <- list.fromNames(c("raw", "filt"))

  # Adapter for Markdownreports background variable "OutDir" ----------------------------------------------------------------
  if (exists('OutDir')) OutDirBac <- OutDir
  OutDir <- file.path(CellRangerOutputDir,p0(kpp("SoupStatistics", SeqRun)))
  try(dir.create(OutDir))
  ww.assign_to_global("OutDir", OutDir, 1)

  # Read In ------------------------
  print("Reading raw CellRanger output matrices")
  CR.matrices$'raw' <- Read10X(path.raw)
  if (length(CR.matrices$'raw') == 2 ) { CR.matrices$'raw' <- CR.matrices$'raw'[[1]] } # Maybe AB table is present too at slot 2!

  print("Reading filtered CellRanger output matrices")
  CR.matrices$'filt' <- Read10X(path.filt)
  if (length(CR.matrices$'filt') == 2 ) { CR.matrices$'filt' <- CR.matrices$'filt'[[1]] } # Maybe AB table is present too at slot 2!

  # Profiling the soup ------------------------
  print("Profiling the soup")
  GEMs.all <- CR.matrices$'raw'@Dimnames[[2]]
  GEMs.cells <- CR.matrices$'filt'@Dimnames[[2]]
  iprint("There are", l(GEMs.all), "GEMs sequenced, and",l(GEMs.cells), "are cells among those." )

  GEMs.soup <- setdiff(GEMs.all, GEMs.cells)
  CR.matrices$'soup' <- CR.matrices$'raw'[,GEMs.soup]
  CR.matrices$'soup.total.RC' <- Matrix::rowSums(CR.matrices$'soup')
  CR.matrices$'soup.total.sum' <- sum(CR.matrices$'soup')
  CR.matrices$'cells.total.sum' <- sum(CR.matrices$'filt')

  CR.matrices$'soup.rel.RC'  <- CR.matrices$'soup.total.RC' / CR.matrices$'soup.total.sum'

  # Diff Exp ----------------------------------------------------------------
  Soup.VS.Cells.Av.Exp <- cbind(
    'Soup' = Matrix::rowSums(CR.matrices$'soup'),
    'Cells' = Matrix::rowSums(CR.matrices$'filt')
  )
  colnames(Soup.VS.Cells.Av.Exp)
  idx.HE <- rowSums(Soup.VS.Cells.Av.Exp)>10; pc_TRUE(idx.HE)
  Soup.VS.Cells.Av.Exp <- Soup.VS.Cells.Av.Exp[idx.HE,]; idim(Soup.VS.Cells.Av.Exp)
  Soup.VS.Cells.Av.Exp.log10 <- log10(Soup.VS.Cells.Av.Exp+1)
  # wplot(Soup.VS.Cells.Av.Exp.log10, col = rgb(0,0,0,.25), PNG = T
  #       , xlab = "Total Expression in Soup [log10(mRNA+1)]"
  #       , ylab = "Total Expression in Cells [log10(mRNA+1)]"
  # )

  # ggplot prepare ----------------------------------------------------------------
  Soup.VS.Cells.Av.Exp.gg <- tibble::rownames_to_column(as.data.frame(Soup.VS.Cells.Av.Exp.log10), "gene")
  (Soup.VS.Cells.Av.Exp.gg <- as_tibble(Soup.VS.Cells.Av.Exp.gg))
  soup.rate <- Soup.VS.Cells.Av.Exp.gg$Soup / (Soup.VS.Cells.Av.Exp.gg$Cells + Soup.VS.Cells.Av.Exp.gg$Soup)
  cell.rate <- Soup.VS.Cells.Av.Exp.gg$Cells / (Soup.VS.Cells.Av.Exp.gg$Cells + Soup.VS.Cells.Av.Exp.gg$Soup)

  axl.pfx <- "Total Expression in"
  axl.sfx <- "[log10(mRNA+1)]"



  HGNC <- Soup.VS.Cells.Av.Exp.gg$gene
  Class <- rep("Other", times=nrow(Soup.VS.Cells.Av.Exp.gg))
  Class[grep('^RPL|^RPS', HGNC)]  <- "RP"
  Class[grep('^MT-', HGNC)] <- "MT"
  Class[grep('^LINC', HGNC)]  <- "LINC"
  Class[grep('^AC', HGNC)]  <- "AC"
  Class[grep('^AL', HGNC)]  <- "AL"
  wpie(table(Class))
  Soup.VS.Cells.Av.Exp.gg$Class <- Class

  fname <- kpp("Soup.VS.Cells.Av.Exp.GeneClasses",SeqRun,"pdf")
  pgg <-
    ggplot(Soup.VS.Cells.Av.Exp.gg %>% arrange(-nchar(Class) )
           , aes(x= Soup, y= Cells, label=gene, col= Class))  +
    geom_abline(slope=1, col='darkgrey') + geom_point()+
    scale_alpha_manual(guide='none', values = ls.Alpha) +
    xlab(paste(axl.pfx, "Soup", axl.sfx)) + ylab(paste(axl.pfx, "Cells", axl.sfx)) +
    ggtitle("Soup VS. Cells | gene classes")

  ggsave(pgg, filename = file.path(OutDir, fname))

  # ggplot ----------------------------------------------------------------
  quantiles <- c(0.025, 0.01, 0.0025)

  i=1
  for (i in 1:l(quantiles)) {
    pr <- quantiles[i]; print(pr)
    HP.thr <- 200*pr/quantiles[2]
    idx.HE2 <- rowSums(Soup.VS.Cells.Av.Exp) > HP.thr
    pc_TRUE(idx.HE2)

    fname <- kpp("Soup.VS.Cells.Av.Exp.quantile",pr,SeqRun,"pdf")

    Outlier <- idx.HE2 &
      (cell.rate < quantile(cell.rate, probs = pr) |
         soup.rate < quantile(soup.rate, probs = pr))

    pc_TRUE(Outlier); sum(Outlier)
    HP.thr.mod <- HP.thr
    while (sum(Outlier) > 40) {
      HP.thr.mod <- HP.thr.mod *2
      Outlier <- Outlier &  rowSums(Soup.VS.Cells.Av.Exp) > HP.thr.mod
    }
    sum(Outlier)


    pgg <-
      ggplot(Soup.VS.Cells.Av.Exp.gg, aes(x= Soup, y= Cells, label=gene,
                                          col= Outlier))  +
      geom_point() + theme(legend.position = "none") +
      xlab(paste(axl.pfx, "Soup", axl.sfx)) + ylab(paste(axl.pfx, "Cells", axl.sfx)) +
      ggtitle("Soup VS. Cells", subtitle = pr) +
      ggrepel::geom_text_repel(aes(label= ifelse(Outlier
                                                 , as.character(gene),'')))
    ggsave(pgg, filename = file.path(OutDir, fname))
  }


  # Per Gene ----------------------------------------------------------------
  PC.mRNA.in.Soup <- sum(CR.matrices$'soup')/sum(CR.matrices$'raw')
  PC.mRNA.in.Cells <- 100*sum(CR.matrices$'filt')/sum(CR.matrices$'raw')
  wbarplot(variable = PC.mRNA.in.Cells, col ="seagreen", plotname = kppd("PC.mRNA.in.Cells", SeqRun)
           , ylim = c(0,100), ylab = "% mRNA in cells"
           , sub = "% mRNA is more meaningful than % reads reported by CR")
  barplot_label(barplotted_variable = PC.mRNA.in.Cells
                , labels = percentage_formatter(PC.mRNA.in.Cells/100, digitz = 2)
                , TopOffset = 10)


  # Plot top gene's expression ----------------------------------------------------------------
  Soup.GEMs.top.Genes = 100*head(sort(CR.matrices$'soup.rel.RC', decreasing = T), n = 20)

  wbarplot(Soup.GEMs.top.Genes, plotname = kppd("Soup.GEMs.top.Genes", SeqRun)
           , ylab="% mRNA in the Soup"
           , sub = paste("Within the", SeqRun, "dataset")
           , tilted_text = T
           , ylim = c(0, max(Soup.GEMs.top.Genes)*1.5))
  barplot_label(barplotted_variable = Soup.GEMs.top.Genes
                , labels = percentage_formatter(Soup.GEMs.top.Genes/100, digitz = 2)
                , TopOffset = -.5, srt = 90, cex=.75)

  # Plot summarize expression ----------------------------------------------------------------
  soupProfile <- CR.matrices$'soup.total.RC'
  {
    soup.RP.sum   <- sum(soupProfile[grep('^RPL|^RPS', names(soupProfile))])
    soup.RPL.sum   <- sum(soupProfile[grep('^RPL', names(soupProfile))])
    soup.RPS.sum   <- sum(soupProfile[grep('^RPS', names(soupProfile))])
    soup.mito.sum <- sum(soupProfile[grep('^MT-', names(soupProfile))])
    soup.LINC.sum <- sum(soupProfile[grep('^LINC', names(soupProfile))])
    soup.AC.sum <- sum(soupProfile[grep('^AC', names(soupProfile))])
    soup.AL.sum <- sum(soupProfile[grep('^AL', names(soupProfile))])
    genes.non.Above <- soupProfile[grepv('^RPL|^RPS|^MT-|^LINC|^AC|^AL', names(soupProfile), invert = T)]
  }
  head(sort(genes.non.Above), n=50)


  soupProfile.summarized <- c(
    'Mitochondial' = soup.mito.sum,
    'Ribosomal' = soup.RP.sum,
    'Ribosomal.L' = soup.RPL.sum,
    'Ribosomal.S' = soup.RPS.sum,
    'GenBank (AC)' = soup.AC.sum,
    'EMBL (AL)' = soup.AL.sum,
    'LINC' = soup.LINC.sum,
    sort(genes.non.Above, decreasing = T)
  )
  NrColumns2Show  = min(10, nrow(soupProfile.summarized))
  ccc <- c("#FF4E00","#778B04","#8ea604","#8ea604","#F5BB00","#F5BB00","#EC9F05",rep(x = "#BF3100", times=NrColumns2Show-6)) # ,"#"


  Soup.GEMs.top.Genes.summarized = 100 * soupProfile.summarized[1:NrColumns2Show] / CR.matrices$'soup.total.sum'
  maxx <- max(Soup.GEMs.top.Genes.summarized)
  wbarplot(Soup.GEMs.top.Genes.summarized, plotname = kppd("Soup.GEMs.top.Genes.summarized", SeqRun)
           , ylab="% mRNA in the Soup", ylim = c(0, maxx+3)
           , sub = paste("Within the", SeqRun, "dataset")
           , tilted_text = T, col = ccc)
  barplot_label(barplotted_variable = Soup.GEMs.top.Genes.summarized
                , srt = 45, labels = percentage_formatter(Soup.GEMs.top.Genes.summarized/100, digitz = 2)
                , TopOffset = -1.5)

  # Absolute.fraction ---------------------------
  Absolute.fraction.soupProfile.summarized <- Soup.GEMs.top.Genes.summarized * PC.mRNA.in.Soup

  maxx <- max(Absolute.fraction.soupProfile.summarized)
  wbarplot(Absolute.fraction.soupProfile.summarized, plotname = kppd("Absolute.fraction.soupProfile.summarized", SeqRun)
           , ylab="% of mRNA in cells", ylim = c(0, maxx*1.33)
           , sub = paste(percentage_formatter(PC.mRNA.in.Soup), "of mRNA counts are in the Soup, in the dataset ", SeqRun)
           , tilted_text = T, col = ccc)
  barplot_label(barplotted_variable = Absolute.fraction.soupProfile.summarized
                , srt = 45, labels = percentage_formatter(Absolute.fraction.soupProfile.summarized/100, digitz = 2)
                # formatC(Absolute.fraction.soupProfile.summarized, format="f", big.mark = " ", digits=0)
                , TopOffset = -maxx*0.15)

  # -----
  Soup.GEMs.top.Genes.non.summarized <- 100* sort(genes.non.Above, decreasing = T)[1:20]/ CR.matrices$'soup.total.sum'
  maxx <- max(Soup.GEMs.top.Genes.non.summarized)
  wbarplot(Soup.GEMs.top.Genes.non.summarized, plotname = kppd("Soup.GEMs.top.Genes.non.summarized", SeqRun)
           , ylab="% mRNA in the Soup"
           , sub = paste("Within the", SeqRun, "dataset")
           , tilted_text = T, col = "#BF3100"
           , ylim = c(0, maxx*1.5))
  barplot_label(barplotted_variable = Soup.GEMs.top.Genes.non.summarized
                , labels = percentage_formatter(Soup.GEMs.top.Genes.non.summarized/100, digitz = 2)
                # , labels = p0(round(1e6 * Soup.GEMs.top.Genes.non.summarized), " ppm")
                , TopOffset = -maxx*0.2, srt = 90, cex=.75)

  # Diff Exp ----------------------------------------------------------------


  # Adapter for Markdownreports background variable "OutDir" ----------------------------------------------------------------
  if (exists('OutDirBac'))  ww.assign_to_global("OutDir", OutDirBac, 1)


} # plotTheSoup
# plotTheSoup(CellRangerOutputDir = "~/Data/114593/114593" , SeqRun = gsub('*([0-9]+).*','\\1', x = basename(CellRangerOutputDir)))



load10Xv3 <- function(dataDir, cellIDs = NULL, channelName = NULL, readArgs = list(),
                      includeFeatures = c("Gene Expression"), verbose = TRUE,
                      ...)
{

  # include
  dirz <- list.dirs(dataDir, full.names = F, recursive = F)
  path.raw <- file.path(dataDir, grep(x = dirz, pattern = "^raw_*", value = T))
  path.filt <- file.path(dataDir, grep(x = dirz, pattern = "^filt_*", value = T))
  CR.matrices <- list.fromNames(c("raw", "filt"))


  (isV3 = any(grepl(x = dirz, pattern = "^raw_feature_bc*")))
  tgt = path.raw

  if (!isV3)
    tgt = file.path(tgt, list.files(tgt))
  if (verbose)
    message(sprintf("Loading raw count data"))
  dat = do.call(Read10X, c(list(data.dir = tgt), readArgs))
  if (verbose)
    message(sprintf("Loading cell-only count data"))
  if (!is.null(cellIDs)) {
    if (all(grepl("\\-1$", cellIDs)))
      cellIDs = gsub("\\-1$", "", cellIDs)
    if (!all(cellIDs %in% colnames(dat)))
      stop("Not all supplied cellIDs found in raw data.")
    datCells = dat[, match(cellIDs, colnames(dat))]
  }
  else {
    tgt = path.filt
    if (!isV3)
      tgt = file.path(tgt, list.files(tgt))
    datCells = do.call(Read10X, c(list(data.dir = tgt),
                                  readArgs))
    if (is.list(dat)) {
      dat = do.call(rbind, dat[includeFeatures])
      datCells = do.call(rbind, datCells[includeFeatures])
    }
  }
  if (verbose)
    message(sprintf("Loading extra analysis data where available"))
  mDat = NULL
  tgt = file.path(dataDir, "analysis", "clustering", "graphclust",
                  "clusters.csv")
  if (file.exists(tgt)) {
    clusters = read.csv(tgt)
    mDat = data.frame(clusters = clusters$Cluster, row.names = clusters$Barcode)
  }
  tgt = file.path(dataDir, "analysis", "clustering", "kmeans_10_clusters",
                  "clusters.csv")
  if (file.exists(tgt)) {
    clusters = read.csv(tgt)
    mDat$clustersFine = clusters$Cluster
  }
  tgt = file.path(dataDir, "analysis", "tsne", "2_components",
                  "projection.csv")
  if (file.exists(tgt)) {
    tsne = read.csv(tgt)
    if (is.null(mDat)) {
      mDat = data.frame(tSNE1 = tsne$TSNE.1, tSNE2 = tsne$TSNE.2,
                        row.names = tsne$Barcode)
    }
    else {
      mDat$tSNE1 = tsne$TSNE.1[match(rownames(mDat), tsne$Barcode)]
      mDat$tSNE2 = tsne$TSNE.2[match(rownames(mDat), tsne$Barcode)]
    }
    DR = c("tSNE1", "tSNE2")
  }
  else {
    DR = NULL
  }
  if (!is.null(mDat) && any(rownames(mDat) != colnames(datCells))) {
    rownames(mDat) = gsub("-1$", "", rownames(mDat))
    if (any(rownames(mDat) != colnames(datCells)))
      stop("Error matching meta-data to cell names.")
  }
  if (is.null(channelName))
    channelName = ifelse(is.null(names(dataDir)), dataDir,
                         names(dataDir))
  channel = SoupX::SoupChannel(tod = dat, toc = datCells, metaData = mDat,
                               channelName = channelName, dataDir = dataDir, dataType = "10X",
                               isV3 = isV3, DR = DR, ...)
  return(channel)
}




# dataDir="/Volumes/copy.your.own.data.here/A.Vertesy/SEO/HNV73DRXX_R10015/HNV73DRXX_R10015/aligned_rna/124851_rnacount"

