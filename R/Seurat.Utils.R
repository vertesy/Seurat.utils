
# source("/Users/abel.vertesy/GitHub/Packages/Seurat.utils/R/Seurat.Utils.R")
# _________________________________________________________________________________________________
# Cluster.Auto-naming.DE.R
# _________________________________________________________________________________________________
# source('~/GitHub/Packages/Seurat.utils/Functions/Cluster.Auto-naming.DE.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Cluster.Auto-naming.DE.R"))

# _________________________________________________________________________________________________
# require(princurve) # only for AutoNumber.by.PrinCurve



# _________________________________________________________________________________________________
#' @title SmallestNonAboveX
#' @description replace small values with the next smallest value found, which is >X. #
#' @param vec Numeric input vector
#' @param X Threshold, Default: 0
#' @examples
#' \dontrun{
#' if(interactive()){
#'  SmallestNonZero(vec = df.markers$"p_val")
#'  }
#' }
#' @export
SmallestNonAboveX <- function(vec, X = 0) { # replace small values with the next smallest value found, which is >X.
  newmin <- min(vec[vec > X])
  vec[vec <= X] <- newmin
  vec
}



# _________________________________________________________________________________________________
#' @title Add.DE.combined.score
#' @description Add combined score to DE results. (LFC * -log10( p_cutoff / pval_scaling ) )
#' @param df Data frame, result of DGEA analysis (FindAllMarkers), Default: df.markers
#' @param p_val_min PARAM_DESCRIPTION, Default: 1e-25
#' @param pval_scaling PARAM_DESCRIPTION, Default: 0.001
#' @param colP PARAM_DESCRIPTION, Default: 'p_val'
#' @param colLFC PARAM_DESCRIPTION, Default: CodeAndRoll2::grepv(pattern = c("avg_logFC|avg_log2FC"), x = colnames(df),
#'    perl = T)
#' @examples
#' \dontrun{
#' if(interactive()){
#'  df.markers <- Add.DE.combined.score(df.markers)
#'  }
#' }
#' @export
Add.DE.combined.score <- function(df = df.markers, p_val_min = 1e-25, pval_scaling = 0.001, colP = "p_val"
                                  , colLFC = CodeAndRoll2::grepv(pattern = c("avg_logFC|avg_log2FC"), x = colnames(df), perl = T)
                                  # , colLFC = "avg_log2FC"
) { # Score = -LOG10(p_val) * avg_log2FC
  p_cutoff <- SmallestNonAboveX(vec = df[[colP]], X = p_val_min)
  df$'combined.score' <- round(df[[colLFC]] * -log10( p_cutoff / pval_scaling ) )
  return(df)
}




# _________________________________________________________________________________________________
#' @title StoreTop25Markers
#' @description Save the top 25 makers based on `avg_log2FC` output table of `FindAllMarkers()` (df_markers) under `@misc$df.markers$res...`. By default, it rounds up insignificant digits up to 3. #
#' @param obj Seurat object, Default: combined.obj
#' @param df_markers Data frame, result of DGEA analysis (FindAllMarkers), Default: df.markers
#' @param res Clustering resoluton to use, Default: 0.5
#' @examples
#' \dontrun{
#' if(interactive()){
#'  combined.obj <- StoreTop25Markers(df_markers = df.markers, res = 0.)
#'  }
#' }
#' @seealso
#'  \code{\link[dplyr]{select}}
#' @export
#' @importFrom dplyr select
StoreTop25Markers <- function(obj = combined.obj # Save the top 25 makers based on `avg_log2FC` output table of `FindAllMarkers()` (df_markers) under `@misc$df.markers$res...`. By default, it rounds up insignificant digits up to 3.
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


# _________________________________________________________________________________________________
#' @title StoreAllMarkers
#' @description Save the output table of `FindAllMarkers()` (df_markers) under `@misc$df.markers$res...`. By default, it rounds up insignificant digits up to 3. #
#' @param obj Seurat object, Default: combined.obj
#' @param df_markers Data frame, result of DGEA analysis (FindAllMarkers), Default: df.markers
#' @param res Clustering resoluton to use, Default: 0.5
#' @param digit Number of digits to keep, Default: c(0, 3)[2]
#' @examples
#' \dontrun{
#' if(interactive()){
#'  combined.obj <- StoreAllMarkers(df_markers = df.markers, res = 0.5)
#'  }
#' }
#' @export
StoreAllMarkers <- function(obj = combined.obj # Save the output table of `FindAllMarkers()` (df_markers) under `@misc$df.markers$res...`. By default, it rounds up insignificant digits up to 3.
                            , df_markers = df.markers, res = 0.5, digit = c(0,3)[2]) {
  if (digit) df_markers[,1:5] <- signif(df_markers[,1:5], digits = digit)
  obj@misc$'df.markers'[[ppp('res',res)]] <- df_markers
  iprint("DF markers are stored under:", 'obj@misc$df.markers$', ppp('res',res))
  return(obj)
}


# _________________________________________________________________________________________________
#' @title GetTopMarkersDF
#' @description Get the vector of N most diff. exp. genes. #
#' @param dfDE Data frame, result of DGEA analysis (FindAllMarkers), Default: df.markers
#' @param n Number of markers to return, Default: p$n.markers
#' @param order.by Sort output tibble by which column, Default: c("avg_log2FC", "p_val_adj")[1]
#' @examples
#' \dontrun{
#' if(interactive()){
#'  GetTopMarkers(df = df.markers, n = 3 )
#'  }
#' }
#' @seealso
#'  \code{\link[dplyr]{slice}}, \code{\link[dplyr]{select}}
#' @export
#' @importFrom dplyr slice select
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


# _________________________________________________________________________________________________
#' @title GetTopMarkers
#' @description Get the vector of N most diff. exp. genes. #
#' @param dfDE Data frame, result of DGEA analysis (FindAllMarkers), Default: df.markers
#' @param n Number of markers to return, Default: p$n.markers
#' @param order.by Sort output tibble by which column, Default: c("combined.score", "avg_log2FC", "p_val_adj")[2]
#' @examples
#' \dontrun{
#' if(interactive()){
#'  GetTopMarkers(df = df.markers, n = 3 )
#'  }
#' }
#' @seealso
#'  \code{\link[dplyr]{slice}}, \code{\link[dplyr]{select}}
#' @export
#' @importFrom dplyr slice select
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




# _________________________________________________________________________________________________
#' @title AutoLabelTop.logFC
#' @description Create a new "named identity" column in the metadata of a Seurat object, with `Ident` set to a clustering output matching the `res` parameter of the function. It requires the output table of `FindAllMarkers()`. If you used `StoreAllMarkers()` is stored under `@misc$df.markers$res...`, which location is assumed by default. #
#' @param obj Seurat object, Default: combined.obj
#' @param res Clustering resoluton to use, Default: 0.2
#' @param plot.top.genes PARAM_DESCRIPTION, Default: T
#' @param order.by Sort output tibble by which column, Default: c("combined.score", "avg_logFC", "p_val_adj")[1]
#' @param df_markers Data frame, result of DGEA analysis (FindAllMarkers), Default: combined.obj@misc$df.markers[[paste0("res.", res)]]
#' @examples
#' \dontrun{
#' if(interactive()){
#'  combined.obj <- AutoLabelTop.logFC(); combined.obj$"cl.names.top.gene.res.0.5"
#'  }
#' }
#' @export
AutoLabelTop.logFC <- function(obj = combined.obj # Create a new "named identity" column in the metadata of a Seurat object, with `Ident` set to a clustering output matching the `res` parameter of the function. It requires the output table of `FindAllMarkers()`. If you used `StoreAllMarkers()` is stored under `@misc$df.markers$res...`, which location is assumed by default.
                               , ident = GetClusteringRuns()[1]
                               , res = 0.2, plot.top.genes = T
                               , suffix = res
                               , order.by = c("combined.score", "avg_log2FC", "p_val_adj")[2]
                               , df_markers = obj@misc$"df.markers"[[paste0("res.",res)]] ) {
  stopifnot(!is.null("df_markers"))
  stopifnot(order.by %in% colnames(df_markers))

  top.markers <-
    GetTopMarkersDF(df = df_markers, order.by = order.by, n = 1) %>%
    col2named.vec.tbl()

  obj@misc[[ppp("top.markers.res",res)]] <- top.markers


  ids <- deframe(unique(obj[[ident]]))
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

  if (plot.top.genes) multiFeaturePlot.A4(list.of.genes = top.markers, suffix = suffix)

  return(obj)
}






# _________________________________________________________________________________________________
#' @title AutoLabel.KnownMarkers
#' @description Create a new "named identity" column in the metadata of a Seurat object, with `Ident` set to a clustering output matching the `res` parameter of the function. It requires the output table of `FindAllMarkers()`. If you used `StoreAllMarkers()` is stored under `@misc$df.markers$res...`, which location is assumed by default. #
#' @param obj Seurat object, Default: combined.obj
#' @param topN Take the top N genes, Default: 1
#' @param res Clustering resolution to use, Default: 0.5
#' @param KnownMarkers PARAM_DESCRIPTION, Default: c("TOP2A", "EOMES", "SLA", "HOPX", "S100B", "DLX6-AS1", "POU5F1",
#'    "SALL4", "DDIT4", "PDK1", "SATB2", "FEZF2")
#' @param order.by Sort output tibble by which column, Default: c("combined.score", "avg_log2FC", "p_val_adj")[1]
#' @param df_markers PARAM_DESCRIPTION, Default: combined.obj@misc$df.markers[[paste0("res.", res)]]
#' @examples
#' \dontrun{
#' if(interactive()){
#'  combined.obj <- AutoLabel.KnownMarkers(); # combined.obj$cl.names.KnownMarkers.0.5
#'  DimPlot.ClusterNames(ident = "cl.names.KnownMarkers.0.5")
#'  }
#' }
#' @seealso
#'  \code{\link[dplyr]{select}}, \code{\link[dplyr]{slice}}
#' @export
#' @importFrom dplyr select slice
AutoLabel.KnownMarkers <- function(obj = combined.obj, topN =1, res = 0.5 # Create a new "named identity" column in the metadata of a Seurat object, with `Ident` set to a clustering output matching the `res` parameter of the function. It requires the output table of `FindAllMarkers()`. If you used `StoreAllMarkers()` is stored under `@misc$df.markers$res...`, which location is assumed by default.
                                   , KnownMarkers = c("TOP2A", "EOMES", "SLA", "HOPX", "S100B", "DLX6-AS1", "POU5F1","SALL4","DDIT4", "PDK1", "SATB2", "FEZF2")
                                   , order.by = c("combined.score", "avg_log2FC", "p_val_adj")[1]

                                   , df_markers = obj@misc$"df.markers"[[paste0("res.",res)]] ) {
  stopifnot(!is.null("df_markers"))

  lfcCOL <- CodeAndRoll2::grepv(pattern = c("avg_logFC|avg_log2FC"), x = colnames(df_markers), perl = T)
  keep <- unique(c(lfcCOL, 'p_val_adj', 'cluster', order.by, 'gene'  ))


  matching.clusters <-
    df_markers %>%
    dplyr::select(keep) %>%
    arrange(desc(!!as.name(order.by))) %>%
    filter(gene %in%  KnownMarkers) %>%
    group_by(gene) %>%
    dplyr::slice(1:topN) %>%
    arrange(desc(!!as.name(order.by))) %>%
    # top_n(n = 1, wt = avg_log2FC) %>% # Select the top cluster for each gene
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





# _________________________________________________________________________________________________
#' @title DimPlot.ClusterNames
#' @description Plot UMAP with Cluster names. #
#' @param obj Seurat object, Default: combined.obj
#' @param ident identity used, Default: 'cl.names.top.gene.res.0.5'
#' @param reduction UMAP, tSNE, or PCA (Dim. reduction to use), Default: 'umap'
#' @param title Title of the plot, Default: ident
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @examples
#' \dontrun{
#' if(interactive()){
#'  DimPlot.ClusterNames()
#'  }
#' }
#' @export
DimPlot.ClusterNames <- function(obj = combined.obj # Plot UMAP with Cluster names.
                                 , ident = "cl.names.top.gene.res.0.5", reduction = "umap", title = ident, ...) {
  Seurat::DimPlot(object = obj, reduction = reduction, group.by = ident, label = T, repel = T, ...) + NoLegend() + ggtitle(title)
}



# _________________________________________________________________________________________________
#' @title AutoNumber.by.UMAP
#' @description Relabel cluster numbers along a UMAP (or tSNE) axis #
#' @param obj Seurat object, Default: combined.obj
#' @param dim Which dimension? Default: 1
#' @param swap Swap direction? Default: F
#' @param reduction UMAP, tSNE, or PCA (Dim. reduction to use), Default: 'umap'
#' @param res Clustering resoluton to use, Default: 'integrated_snn_res.0.5'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  combined.obj <- AutoNumber.by.UMAP(obj = combined.obj, dim = 1, reduction="umap", res = "integrated_snn_res.0.5" );
#'  DimPlot.ClusterNames(ident = "integrated_snn_res.0.5.ordered")
#'  }
#' }
#' @export
AutoNumber.by.UMAP <- function(obj = combined.obj # Relabel cluster numbers along a UMAP (or tSNE) axis
                               , dim = 1, swap= F, reduction="umap", res = "RNA_snn_res.0.5" ) {

  dim_name <- kppu(toupper(reduction),dim)
  coord.umap <- as.named.vector(FetchData(object = obj, vars = dim_name))
  identX <- as.character(obj@meta.data[[res]])

  ls.perCl <- split(coord.umap, f = identX)
  MedianClusterCoordinate <- unlapply(ls.perCl, median)

  OldLabel <- names(sort(MedianClusterCoordinate, decreasing = swap))
  NewLabel <- as.character(0:(length(MedianClusterCoordinate) - 1))
  NewMeta <- translate(vec = identX, oldvalues = OldLabel, newvalues = NewLabel)
  NewMetaCol <- kpp(res,"ordered")
  iprint("NewMetaCol:",NewMetaCol)

  obj[[NewMetaCol]] <- NewMeta
  return(obj)
}



# _________________________________________________________________________________________________
#' @title AutoNumber.by.PrinCurve
#' @description Relabel cluster numbers along the principal curve of 2 UMAP (or tSNE) dimensions. #
#' @param obj Seurat object, Default: combined.obj
#' @param dim Dimensions to use, Default: 1:2
#' @param plotit Plot results (& show it), Default: T
#' @param swap Swap Lambda paramter (multiplied with this) , Default: -1
#' @param reduction UMAP, tSNE, or PCA (Dim. reduction to use), Default: 'umap'
#' @param res Clustering resoluton to use, Default: 'integrated_snn_res.0.5'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  DimPlot.ClusterNames(ident = "integrated_snn_res.0.5")
#'  combined.obj <- AutoNumber.by.PrinCurve(obj = combined.obj, dim = 1:2, reduction="umap", plotit = T,
#'  swap= -1, res = "integrated_snn_res.0.5" )
#'  DimPlot.ClusterNames(ident = "integrated_snn_res.0.5.prin.curve")
#'  }
#' }
#' @seealso
#'  \code{\link[princurve]{principal_curve}}
#' @export
#' @importFrom princurve principal_curve
AutoNumber.by.PrinCurve <- function(obj = combined.obj # Relabel cluster numbers along the principal curve of 2 UMAP (or tSNE) dimensions.
                                    , dim = 1:2, plotit = T, swap= -1
                                    , reduction="umap", res = "integrated_snn_res.0.5" ) {
  # require(princurve)
  dim_name <- ppu(toupper(reduction),dim)
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
  NewMeta <- translate(vec = obj[[res]], oldvalues = OldLabel, newvalues = NewLabel)
  NewMetaCol <- kpp(res,"prin.curve")
  iprint("NewMetaCol:",NewMetaCol)
  obj[[NewMetaCol]] <- NewMeta
  return(obj)
}



# _________________________________________________________________________________________________
# _________________________________________________________________________________________________


# _________________________________________________________________________________________________
# Jaccard.toolkit.R
# _________________________________________________________________________________________________
# try(source('~/GitHub/Packages/Seurat.utils/Functions/Jaccard.toolkit.R'))
# try(source('https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Jaccard.toolkit.R'))


#  __________________________________________
# Fast direct calculation from a list


# _________________________________________________________________________________________________
#' @title jJaccardIndexVec
#' @description Calculate jaccard similarity for 2 vecotrs. Helper to jPairwiseJaccardIndexList.
#' @param A Set A, Default: 1:3
#' @param B Set B, Default: 2:4
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
jJaccardIndexVec <- function(A = 1:3, B = 2:4) length(intersect(A,B)) / length(union(A,B))

# _________________________________________________________________________________________________

#' @title jPairwiseJaccardIndexList
#' @description Create a pairwise jaccard similarity matrix across all combinations of columns in binary.presence.matrix. Modified from: https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/ #
#' @param lsG List of genes, Default: ls_genes
#' @examples
#' \dontrun{
#' if(interactive()){
#'  jPairwiseJaccardIndexList(lsG = ls_genes)
#'  }
#' }
#' @export
#' @importFrom Stringendo percentage_formatter
jPairwiseJaccardIndexList <- function(lsG = ls_genes) { # Create a pairwise jaccard similarity matrix across all combinations of columns in binary.presence.matrix. Modified from: https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/
  if (length(names(lsG)) < length(lsG)) {
    iprint("Gene lists were not (all) named, now renamed as:")
    names(lsG) <- ppp("dataset", 1:length(lsG))
    print(names(lsG))
  }
  m = matrix.fromNames(rowname_vec = names(lsG), colname_vec = names(lsG))
  n.sets <- length(lsG)
  for (r in 1:n.sets) {
    # print(Stringendo::percentage_formatter(r/n.sets))
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







# Much slower Indirect calculation via PresenceMatrix
# _________________________________________________________________________________________________


# _________________________________________________________________________________________________
#' @title jPresenceMatrix
#' @description Make a binary presence matrix from a list. Source: https://stackoverflow.com/questions/56155707/r-how-to-create-a-binary-relation-matrix-from-a-list-of-strings #
#' @param string_list PARAM_DESCRIPTION, Default: lst(a = 1:3, b = 2:5, c = 4:9, d = -1:4)
#' @examples
#' \dontrun{
#' if(interactive()){
#'  df.presence <- jPresenceMatrix(string_list = lst(a = 1:3, b = 2:5,c = 4:9, d=-1:4))
#'  }
#' }
#' @export
jPresenceMatrix <- function(string_list = lst(a = 1:3, b = 2:5,c = 4:9, d=-1:4) ) { # Make a binary presence matrix from a list. Source: https://stackoverflow.com/questions/56155707/r-how-to-create-a-binary-relation-matrix-from-a-list-of-strings
  df.presence <- string_list %>%
    enframe %>%
    unnest(cols = "value") %>%
    count(name, value) %>%
    spread(value, n, fill = 0)
  df.presence2 <- FirstCol2RowNames(df.presence)
  return(t(df.presence2))
}


# _________________________________________________________________________________________________
#' @title jJaccardIndexBinary
#' @description Calculate Jaccard Index. Modified from: https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/ #
#' @param x Set X
#' @param y Set Y
#' @examples
#' \dontrun{
#' if(interactive()){
#'  JaccardSimilarity <- jJaccardIndexBinary(  x = sample(x = 0:1, size = 100, replace = T),
#'  y = sample(x = 0:1, size = 100, replace = T))
#'  }
#' }
#' @export
jJaccardIndexBinary <- function(x, y) { # Calculate Jaccard Index. Modified from: https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/
  elements.found <- sort(unique(union(x, y)))
  stopifnot(length(elements.found) == 2) # check if you only have [0,1]
  stopifnot(as.numeric(elements.found) == 0:1) # check if you only have [0,1]

  M.11 = sum(x == 1 & y == 1)
  M.10 = sum(x == 1 & y == 0)
  M.01 = sum(x == 0 & y == 1)
  return (M.11 / (M.11 + M.10 + M.01))
}




# _________________________________________________________________________________________________
#' @title jPairwiseJaccardIndex
#' @description Create a pairwise jaccard similarity matrix across all combinations of columns in binary.presence.matrix. Modified from: https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/ #
#' @param binary.presence.matrix PARAM_DESCRIPTION, Default: df.presence
#' @examples
#' \dontrun{
#' if(interactive()){
#'  PairwiseJaccardIndices <- jPairwiseJaccardIndex(binary.presence.matrix = df.presence)
#'  }
#' }
#' @export
#' @importFrom Stringendo percentage_formatter
jPairwiseJaccardIndex <- function(binary.presence.matrix = df.presence) { # Create a pairwise jaccard similarity matrix across all combinations of columns in binary.presence.matrix. Modified from: https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/
  m = matrix.fromNames(rowname_vec = colnames(binary.presence.matrix), colname_vec = colnames(binary.presence.matrix) )
  n.sets <- ncol(binary.presence.matrix)
  for (r in 1:n.sets) {
    print(Stringendo::percentage_formatter(r/n.sets))
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

# _________________________________________________________________________________________________
# metadata.manipulation.R
# ____________________________________________________________________ ----
# source('~/GitHub/Packages/Seurat.utils/Functions/Seurat.object.manipulations.etc.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Metadata.manipulation.R"))
# Source: self + web

# - getMedianMetric
# - add.meta.tags
# - add.meta.fraction
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
# - calc.q99.Expression.and.set.all.genes
# - PlotTopGenes
# - fix.orig.ident
# - set.all.genes
# - recall.all.genes
# - recall.meta.tags.n.datasets
# - recall.parameters
# - recall.genes.l
# - save.parameters
# - plot.expression.rank.q90
# - FlipReductionCoordinates
# - SeuratColorVector
# - getClusterColors
# - get.clustercomposition


# _________________________________________________________________________________________________
#' @title getMedianMetric
#' @description Get the median values of different columns in meta.data, can iterate over a list of Seurat objects.
#' @param ls.obj List of Seurat objects, Default: ls.Seurat
#' @param n.datasets lenght of list (n objects), Default: length(ls.Seurat)
#' @param mColname PARAM_DESCRIPTION, Default: 'percent.mito'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  ls.Seurat <- getMedianMetric(ls.obj = ls.Seurat, n.datasets = length(ls.Seurat),
#'  mColname = "percent.mito")
#'  }
#' }
#' @export
getMedianMetric <- function(ls.obj = ls.Seurat, n.datasets = length(ls.Seurat), mColname = "percent.mito") {
  medMetric <- vec.fromNames(names(ls.obj))
  for(i in 1:n.datasets ) {
    medMetric[i] <- median(ls.obj[[i]]@meta.data[,mColname])
  }
  return(medMetric)
}




# _________________________________________________________________________________________________
#' @title add.meta.tags
#' @description N is the for which dataset #
#' @param list.of.tags PARAM_DESCRIPTION, Default: tags
#' @param obj Seurat object, Default: ls.Seurat[[1]]
#' @param n PARAM_DESCRIPTION, Default: 1
#' @examples
#' \dontrun{
#' if(interactive()){
#'  ls.Seurat[[1]] <- add.meta.tags(list.of.tags = tags, obj = ls.Seurat[[1]], n = 1)
#'  }
#' }
#' @export
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



# _________________________________________________________________________________________________
#' @title add.meta.fraction
#' @description Add a new metadata column, with the fraction of gene set in the transcripome (percentage).
#' @param col.name PARAM_DESCRIPTION, Default: 'percent.mito'
#' @param gene.symbol.pattern PARAM_DESCRIPTION, Default: c("^MT\\.|^MT-", F)[1]
#' @param gene.set PARAM_DESCRIPTION, Default: F
#' @param obj Seurat object, Default: ls.Seurat[[1]]
#' @param verbose PARAM_DESCRIPTION, Default: T
#' @examples
#' \dontrun{
#' if(interactive()){
#'  ls.Seurat[[1]] <- add.meta.fraction(col.name = "percent.mito", gene.symbol.pattern = "^MT\\.|^MT-")
#'  ls.Seurat[[1]] <- add.meta.fraction(col.name = "percent.ribo", gene.symbol.pattern = "^RPL|^RPS")
#'  ls.Seurat[[1]] <- add.meta.fraction(col.name = "percent.AC.GenBank", gene.symbol.pattern = "^AC[0-9]{6}\\.")
#'  ls.Seurat[[1]] <- add.meta.fraction(col.name = "percent.AL.EMBL", gene.symbol.pattern = "^AL[0-9]{6}\\.")
#'  ls.Seurat[[1]] <- add.meta.fraction(col.name = "percent.LINC", gene.symbol.pattern = "^LINC0")
#'  ls.Seurat[[1]] <- add.meta.fraction(col.name = "percent.MALAT1", gene.symbol.pattern = "^MALAT1")
#'  colnames(ls.Seurat[[1]]@meta.data)
#'  HGA_MarkerGenes <- c("ENO1", "IGFBP2", "WSB1", "DDIT4", "PGK1", "BNIP3", "FAM162A", "TPI1",
#'  "VEGFA", "PDK1", "PGAM1", "IER2", "FOS", "BTG1", "EPB41L4A-AS1","NPAS4", "HK2", "BNIP3L",
#'  "JUN", "ENO2", "GAPDH", "ANKRD37", "ALDOA", "GADD45G", "TXNIP")
#'  sobj <- add.meta.fraction(col.name = "percent.HGA", gene.set = HGA_MarkerGenes, obj =  sobj)
#'  }
#' }
#' @seealso
#'  \code{\link[Matrix]{colSums}}
#' @export
#' @importFrom Matrix colSums
add.meta.fraction <- function(col.name = "percent.mito", gene.symbol.pattern = c("^MT\\.|^MT-", F)[1]
                              , gene.set = F, obj = ls.Seurat[[1]], verbose = T) {
  stopif2(condition = isFALSE(gene.set) && isFALSE(gene.symbol.pattern), "Either gene.set OR gene.symbol.pattern has to be defined (!= FALSE).")
  if (!isFALSE(gene.set) && !isFALSE(gene.symbol.pattern) && verbose) print("Both gene.set AND gene.symbol.pattern are defined. Only using gene.set.")

  if (!isFALSE(gene.set)) geneset <- check.genes(list.of.genes = gene.set, obj = obj)
  total_expr <- Matrix::colSums(GetAssayData(object = obj))
  genes.matching <- if (!isFALSE(gene.set)) intersect(gene.set, rownames(obj)) else CodeAndRoll2::grepv(pattern = gene.symbol.pattern, x = rownames(obj))

  genes.expr = GetAssayData(object = obj)[genes.matching, ]
  target_expr <- if (length(genes.matching) >1) Matrix::colSums(genes.expr) else genes.expr

  iprint(length(genes.matching), "genes found, :", head(genes.matching))

  obj <- AddMetaData(object = obj, metadata = target_expr / total_expr, col.name = col.name)
  colnames(obj@meta.data)
  return(obj)
}




# _________________________________________________________________________________________________
#' @title GetClusteringRuns
#' @description Get Clustering Runs: metadata column names #
#' @param obj Seurat object, Default: combined.obj
#' @param res Clustering resoluton to use, Default: F
#' @param pat Pettern to match, Default: '*snn_res.*[0-9]$'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  GetClusteringRuns()
#'  }
#' }
#' @export
GetClusteringRuns <- function(obj = combined.obj, res = F, pat = "*snn_res.*[0-9]$") { # Get Clustering Runs: metadata column names
  if (res) pat = gsub(x = pat, pattern = '\\[.*\\]', replacement = res)
  clustering.results <- CodeAndRoll2::grepv(x = colnames(obj@meta.data), pattern = pat)
  if ( identical(clustering.results, character(0)) ) warning("No matching column found!")
  return(clustering.results)
}



# _________________________________________________________________________________________________
#' @title GetNamedClusteringRuns
#' @description Get Clustering Runs: metadata column names #
#' @param obj Seurat object, Default: combined.obj
#' @param res Clustering resoluton to use, Default: c(F, 0.5)[1]
#' @param topgene Match clustering named after top expressed gene (see vertesy/Seurat.pipeline/~Diff gene expr.), Default: F
#' @param pat Pettern to match, Default: '^cl.names.Known.*[0,1]\.[0-9]$'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  GetNamedClusteringRuns()
#'  }
#' }
#' @export
GetNamedClusteringRuns <- function(obj = combined.obj  # Get Clustering Runs: metadata column names
                                   , res = c(F, 0.5)[1], topgene = F, pat = "^cl.names.Known.*[0,1]\\.[0-9]$") {
  if (res) pat = gsub(x = pat, pattern = '\\[.*\\]', replacement = res)
  if (topgene) pat = gsub(x = pat, pattern = 'Known', replacement = 'top')
  clustering.results <- CodeAndRoll2::grepv(x = colnames(obj@meta.data), pattern = pat)
  if ( identical(clustering.results, character(0)) ) {
    print("Warning: NO matching column found! Trying GetClusteringRuns(..., pat = '*_res.*[0,1]\\.[0-9]$)")
    clustering.results <- GetClusteringRuns(obj = obj, res = F, pat = "*_res.*[0,1]\\.[0-9]$")
  }
  return(clustering.results)
}



# _________________________________________________________________________________________________
#' @title GetOrderedClusteringRuns
#' @description Get Clustering Runs: metadata column names #
#' @param obj Seurat object, Default: combined.obj
#' @param res Clustering resoluton to use, Default: F
#' @param pat Pettern to match, Default: '*snn_res.*[0,1]\.[0-9]\.ordered$'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  GetOrderedClusteringRuns(); GetOrderedClusteringRuns(res = 0.5)
#'  }
#' }
#' @export
GetOrderedClusteringRuns <- function(obj = combined.obj, res = F, pat = "*snn_res.*[0,1]\\.[0-9]\\.ordered$") { # Get Clustering Runs: metadata column names
  if (res) pat = gsub(x = pat, pattern = '\\[.*\\]', replacement = res)
  clustering.results <- CodeAndRoll2::grepv(x = colnames(obj@meta.data), pattern = pat)
  if ( identical(clustering.results, character(0)) ) warning("No matching column found!")
  return(clustering.results)
}



# _________________________________________________________________________________________________
#' @title GetNumberOfClusters
#' @description Get Number Of Clusters #
#' @param obj Seurat object, Default: combined.obj
#' @examples
#' \dontrun{
#' if(interactive()){
#'  GetNumberOfClusters()
#'  }
#' }
#' @export
GetNumberOfClusters <- function(obj = combined.obj) { # Get Number Of Clusters
  clustering.results <- GetClusteringRuns(obj)
  print("## Number of clusters: ---------")
  for (cc in clustering.results) {
    NrCl <- length(unique(obj@meta.data[[cc]]))
    iprint( cc, "   ", NrCl)
  }
}


# _________________________________________________________________________________________________
#' @title getMetadataColumn
#' @description Get a metadata column from a Seurat object as a named vector. get Cells from metadata.
#' @param ColName.metadata PARAM_DESCRIPTION, Default: 'batch'
#' @param obj Seurat object, Default: combined.obj
#' @param as_numeric PARAM_DESCRIPTION, Default: F
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
getMetadataColumn <- mmeta <- function(ColName.metadata = 'batch', obj = combined.obj, as_numeric =F) { # Get a metadata column from a Seurat object as a named vector
  stopifnot(ColName.metadata %in% colnames(obj@meta.data))

  x = as.named.vector(obj@meta.data[ ,ColName.metadata, drop = F])
  if (as_numeric) {
    as.numeric.wNames(x)+1
  } else {x}
}

# _________________________________________________________________________________________________
#' @title getCellIDs.from.meta
#' @description Get cellIDs from a metadata column, matching a list of values (using %in%). #
#' @param ColName.meta PARAM_DESCRIPTION, Default: 'res.0.6'
#' @param values PARAM_DESCRIPTION, Default: NA
#' @param obj Seurat object, Default: combined.obj
#' @param inverse PARAM_DESCRIPTION, Default: F
#' @examples
#' \dontrun{
#' if(interactive()){
#'  getCellIDs.from.meta()
#'  }
#' }
#' @export
getCellIDs.from.meta <- function(ColName.meta = 'res.0.6', values = NA, obj = combined.obj, inverse = F ) { # Get cellIDs from a metadata column, matching a list of values (using %in%).
  mdat <- obj@meta.data[ , ColName.meta]
  cells <- if (inverse) {mdat %!in% values} else {mdat %in% values}
  idx.matching.cells = which(cells)
  iprint(length(idx.matching.cells), 'cells found.')
  return(rownames(obj@meta.data)[idx.matching.cells])
}


# _________________________________________________________________________________________________
#' @title seu.add.meta.from.vector
#' @description Add a new metadata column to a Seurat  object
#' @param obj Seurat object, Default: combined.obj
#' @param metaD.colname PARAM_DESCRIPTION, Default: metaD.colname.labeled
#' @param Label.per.cell PARAM_DESCRIPTION, Default: Cl.Label.per.cell
#' @examples
#' \dontrun{
#' if(interactive()){
#'  combined.obj <- add.Cl.Label.2.Metadata(obj = combined.obj,  # formerly add.Cl.Label.2.Metadata
#'  metaD.colname = metaD.colname.labeled, Label.per.cell = Cl.Label.per.cell )
#'  }
#' }
#' @export
seu.add.meta.from.vector <- function(obj = combined.obj, metaD.colname = metaD.colname.labeled, Label.per.cell = Cl.Label.per.cell ) { # Add a new metadata column to a Seurat  object
  obj@meta.data[, metaD.colname ] = Label.per.cell
  iprint(metaD.colname, "contains the named identitites. Use Idents(combined.obj) = '...'. The names are:", unique(Label.per.cell))
  return(obj)
}



# _________________________________________________________________________________________________

#' @title seu.map.and.add.new.ident.to.meta
#' @description Add a new metadata column to a Seurat  object
#' @param obj Seurat object, Default: combined.obj
#' @param ident.table PARAM_DESCRIPTION, Default: clusterIDs.GO.process
#' @param metaD.colname PARAM_DESCRIPTION, Default: substitute(ident.table)
#' @examples
#' \dontrun{
#' if(interactive()){
#'  combined.obj <- seu.map.and.add.new.ident.to.meta(obj = combined.obj,
#'  dent.table = clusterIDs.GO.process)
#'  }
#' }
#' @export
seu.map.and.add.new.ident.to.meta <- function(obj = combined.obj, ident.table = clusterIDs.GO.process, orig.ident = Idents(obj)
                                              , metaD.colname = substitute(ident.table) ) { # Add a new metadata column to a Seurat  object
  # identities should match
  {
    Idents(obj) <- orig.ident
    ident.vec <- as.named.vector(ident.table)
    ident.X <- names(ident.vec)
    ident.Y <- as.character(ident.vec)
    ident.Seu <- sort.natural(levels(Idents(obj)))
    iprint("ident.Seu: ", ident.Seu)

    OnlyInIdentVec      <- setdiff(ident.X, ident.Seu)
    OnlyInSeuratIdents  <- setdiff(ident.Seu, ident.X)

    msg.IdentVec <- kollapse("Rownames of 'ident.table' have entries not found in 'Idents(obj)':"
                             , OnlyInIdentVec, " not found in ", ident.Seu, collapseby = " ")

    msg.Seu <- kollapse("Rownames of 'Idents(obj)' have entries not found in 'ident.table':"
                        , OnlyInSeuratIdents, " not found in ", ident.X, collapseby = " ")

    stopif (length(OnlyInIdentVec), message = msg.IdentVec)
    stopif (length(OnlyInSeuratIdents), message = msg.Seu)
  }
  # identity mapping
  {
    new.ident <- translate(vec = as.character(Idents(obj)), oldvalues = ident.X, newvalues = ident.Y)
    obj@meta.data[[metaD.colname]] = new.ident
    iprint(metaD.colname, "contains the named identitites. Use Idents(combined.obj) = '...'. The names are:"); cat(paste0("\t", ident.Y, "\n"))
  }
}





# ------------------------------------------------
#' @title calc.cluster.averages
#' @description Calculate the average of a metadata column (numeric) per cluster.
#' @param col_name PARAM_DESCRIPTION, Default: 'Score.GO.0006096'
#' @param plot.UMAP.too PARAM_DESCRIPTION, Default: TRUE
#' @param return.plot PARAM_DESCRIPTION, Default: F
#' @param obj Seurat object, Default: combined.obj
#' @param split_by PARAM_DESCRIPTION, Default: GetNamedClusteringRuns()[1]
#' @param scale.zscore PARAM_DESCRIPTION, Default: FALSE
#' @param simplify PARAM_DESCRIPTION, Default: T
#' @param plotit Plot results (& show it), Default: T
#' @param histogram PARAM_DESCRIPTION, Default: FALSE
#' @param nbins PARAM_DESCRIPTION, Default: 50
#' @param report Summary sentence, Default: F
#' @param suffix A suffix added to the filename, Default: NULL
#' @param stat PARAM_DESCRIPTION, Default: c("mean", "median")[2]
#' @param quantile.thr PARAM_DESCRIPTION, Default: 0.9
#' @param absolute.thr PARAM_DESCRIPTION, Default: FALSE
#' @param filter PARAM_DESCRIPTION, Default: c(FALSE, "above", "below")[1]
#' @param ylab.text PARAM_DESCRIPTION, Default: paste("Cluster", stat, "score")
#' @param title Title of the plot, Default: paste("Cluster", stat, col_name)
#' @param subtitle PARAM_DESCRIPTION, Default: NULL
#' @param width PARAM_DESCRIPTION, Default: 8
#' @param height PARAM_DESCRIPTION, Default: 6
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @param xlb PARAM_DESCRIPTION, Default: if (absolute.thr) paste("Threshold at", absolute.thr) else paste("Black lines: ",
#'    kppd(Stringendo::percentage_formatter(c(1 - quantile.thr, quantile.thr))),
#'    "quantiles |", "Cl. >", Stringendo::percentage_formatter(quantile.thr),
#'    "are highlighted. |", split_by)
#' @param fname File name, Default: ppp(col_name, split_by, "cluster.average.barplot.pdf", ...)
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @importFrom Stringendo percentage_formatter

calc.cluster.averages <- function(col_name = "Score.GO.0006096"
                                  , plot.UMAP.too = TRUE
                                  , return.plot = F
                                  , obj =  combined.obj
                                  , split_by = GetNamedClusteringRuns()[1]
                                  , scale.zscore = FALSE
                                  , simplify = T, plotit = T
                                  , histogram = FALSE, nbins = 50
                                  , suffix = NULL
                                  , stat = c("mean", "median")[2]
                                  , quantile.thr = 0.9
                                  , absolute.thr = FALSE
                                  , filter = c(FALSE, 'above', 'below')[1]
                                  , ylab.text = paste("Cluster", stat, "score")
                                  , title = paste("Cluster", stat, col_name)
                                  , prefix.cl.names= FALSE
                                  , report = TRUE
                                  , subtitle = NULL
                                  , width = 8, height =6
                                  , ...
                                  # , ylb = paste(ylab.text, col_name)
                                  # , xlb = paste("Clusters >",Stringendo::percentage_formatter(quantile.thr),"quantile are highlighted. |", split_by)
                                  , xlb = if (absolute.thr) paste("Threshold at", absolute.thr) else paste(
                                    "Black lines: " , kppd(Stringendo::percentage_formatter(c(1-quantile.thr, quantile.thr))) ,"quantiles |"
                                    , "Cl. >",Stringendo::percentage_formatter(quantile.thr),"are highlighted. |", split_by
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
    names(av.score) <- if ( !isFALSE(prefix.cl.names)) ppp("cl",df.summary[[1]]) else df.summary[[1]]
    av.score <- sortbyitsnames(av.score)
    if (scale.zscore) av.score <- (scale(av.score)[,1])

    cutoff <- if(absolute.thr) absolute.thr else quantile(av.score, quantile.thr)
    cutoff.low <- if(absolute.thr) NULL else  quantile(av.score, (1-quantile.thr) )

    iprint('quantile.thr:', quantile.thr)
    if (plotit) {
      if (histogram) {
        p <- ggExpress::qhistogram(vec = as.numeric(av.score), save = F
                        , vline = cutoff
                        , plotname = ppp(title, quantile.thr)
                        , bins = nbins
                        , subtitle = paste(subtitle, "| median in blue/dashed")
                        , ylab = ylab.text
                        , xlab = xlb # Abused
                        , xlab.angle = 45
                        # , ylim = c(-1,1)
                        , ...
                        # , ext = "png", w = 7, h = 5
        ) + geom_vline(xintercept = cutoff.low, lty = 2)
        print(p)
        title_ <- ppp(title, suffix, flag.nameiftrue(scale.zscore))
        ggExpress::qqSave(ggobj = p, title = title_, ext = "png", w = width, h = height)
      } else {
        p <- ggExpress::qbarplot(vec = av.score, save = F
                      , hline = cutoff
                      , plotname = title
                      , suffix = quantile.thr
                      , subtitle = subtitle
                      , ylab = ylab.text
                      , xlab = xlb # Abused
                      , xlab.angle = 45
                      # , ylim = c(-1,1)
                      , ...
                      # , ext = "png", w = 7, h = 5
        ) + geom_hline(yintercept = cutoff.low , lty = 2)

        print(p)
        title_ <- ppp(title, suffix, flag.nameiftrue(scale.zscore))
        qqSave(ggobj = p, title = title_, fname = ppp(title_, split_by, "png"),  w = width, h = height)
      }
    }

    if (report) print(paste0(col_name, ": ", paste(iround(av.score), collapse = " vs. ")))
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


# _________________________________________________________________________________________________
#' @title seu.add.meta.from.table
#' @description Add multiple new metadata columns to a Seurat object from a table. #
#' @param obj Seurat object, Default: seu.ORC
#' @param meta PARAM_DESCRIPTION, Default: MetaData.ORC
#' @param suffix A suffix added to the filename, Default: '.fromMeta'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  combined.obj <- seu.add.meta.from.table()
#'  }
#' }
#' @export
seu.add.meta.from.table <- function(obj = combined.obj, meta = MetaData.ORC, suffix = ".fromMeta") { # Add multiple new metadata columns to a Seurat object from a table.
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
}


# _________________________________________________________________________________________________
#' @title sampleNpc
#' @description Sample N % of a dataframe (obj@metadata), and return the cell IDs. #
#' @param metaDF PARAM_DESCRIPTION, Default: MetaData[which(Pass), ]
#' @param pc PARAM_DESCRIPTION, Default: 0.1
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
sampleNpc <- function(metaDF = MetaData[which(Pass),], pc = 0.1) { # Sample N % of a dataframe (obj@metadata), and return the cell IDs.
  cellIDs = rownames(metaDF)
  nr_cells = floor(length(cellIDs) * pc)
  cellIDs.keep = sample(cellIDs, size = nr_cells, replace = FALSE)
  return(cellIDs.keep)
}


# _________________________________________________________________________________________________
#' @title calc.q99.Expression.and.set.all.genes
#' @description Calculate the gene expression of the e.g.: 90th quantile (expression in the top 10% cells). #
#' @param obj Seurat object, Default: combined.obj
#' @param quantileX Quantile level, Default: 0.9
#' @param max.cells PARAM_DESCRIPTION, Default: 1e+05
#' @param slot slot in the Seurat object. Default: 'data'
#' @param assay RNA or integrated assay, Default: c("RNA", "integrated")[1]
#' @param set.all.genes Create the "all.genes" variable in the global env?, Default: TRUE
#' @param show PARAM_DESCRIPTION, Default: TRUE
#' @examples
#' \dontrun{
#' if(interactive()){
#'  combined.obj <- calc.q99.Expression.and.set.all.genes(obj = combined.obj, quantileX = 0.9,
#'  max.cells =  25000)
#'  head(sort(as.numeric.wNames(obj@misc$expr.q90), decreasing = T))
#'  combined.obj <- calc.q99.Expression.and.set.all.genes(obj = combined.obj, quantileX = 0.95,
#'  max.cells =  25000, set.all.genes = FALSE)
#'  }
#' }
#' @seealso
#'  \code{\link[sparseMatrixStats]{character(0)}}
#' @export
#' @importFrom tictoc tic toc
#' @importFrom sparseMatrixStats rowQuantiles
calc.q99.Expression.and.set.all.genes <- function(obj = combined.obj # Calculate the gene expression of the e.g.: 90th quantile (expression in the top 10% cells).
                                                  , quantileX = 0.99, max.cells =  1e5
                                                  , slot = "data", assay = c("RNA", "integrated")[1]
                                                  , set.all.genes = TRUE, show = TRUE) {
  tictoc::tic()
  x = GetAssayData(object = obj, assay = assay, slot = slot) #, assay = 'RNA'
  if (ncol(x) > max.cells) {
    dsampled = sample(x = 1:ncol(x), size = max.cells)
    x = x[ , dsampled]
  }
  qname = paste0("q", quantileX * 100)
  slot_name = kpp("expr", qname)

  # expr.q99 = iround(apply(x, 1, quantile, probs = quantileX) )
  print("Calculating Gene Quantiles")
  expr.q99.df = sparseMatrixStats::rowQuantiles(x, probs = quantileX)
  expr.q99 = iround(expr.q99.df)
  # expr.q99 = iround(as.named.vector(expr.q99.df))
  toc();

  log2.gene.expr.of.the.90th.quantile <- as.numeric(log2(expr.q99 + 1)) # strip names
  n.cells <- floor(ncol(obj) * (1-quantileX) )
  qnameP <- p0(100*quantileX,'th quantile')
  try(
    ggExpress::qhistogram(log2.gene.expr.of.the.90th.quantile, ext = "pdf", breaks = 30
               , plotname = paste("Gene expression in the", qnameP )
               , subtitle = kollapse(pc_TRUE(expr.q99 > 0, NumberAndPC = T), " genes have ", qname ," expr. > 0.")
               , caption = paste(n.cells, 'cells in', qnameP)
               , xlab = paste0("log2(expr. in the ", qnameP,"quantile+1) [UMI]")
               , ylab = "Nr. of genes"
               , plot = show, save = TRUE
               , vline  = .15
               , filtercol = T
               , palette_use = 'npg'
    )
    , silent = TRUE)

  {
    all.genes = percent_rank(expr.q99);
    names(all.genes) = names(expr.q99);
    all.genes <- sort.decreasing(all.genes)
    if (set.all.genes) obj@misc$'all.genes' = all.genes = as.list(all.genes)
    assign('all.genes', all.genes, envir = as.environment(1))
  }

  obj@misc[[slot_name]] <-  expr.q99

  iprint('Quantile', quantileX ,'is now stored under obj@misc$all.genes and $', slot_name, ' Please execute all.genes <- obj@misc$all.genes.')
  return(obj)
}




# _________________________________________________________________________________________________
#' @title PlotTopGenes
#' @description Plot the highest expressed genes on umaps, in a subfolder. Requires calling calc.q99.Expression.and.set.all.genes before. #
#' @param obj Seurat object, Default: combined.obj
#' @param n PARAM_DESCRIPTION, Default: 32
#' @examples
#' \dontrun{
#' if(interactive()){
#'  PlotTopGenes()
#'  }
#' }
#' @export
PlotTopGenes <- function(obj = combined.obj, n = 32, exp.slot = 'expr.q99'){ # Plot the highest expressed genes on umaps, in a subfolder. Requires calling calc.q99.Expression.and.set.all.genes before.
  Highest.Expressed.Genes = names(head(sort(obj@misc[[exp.slot]], decreasing = T), n = n))
  multiFeaturePlot.A4(list.of.genes = Highest.Expressed.Genes, foldername = "Highest.Expressed.Genes" )
}



# _________________________________________________________________________________________________
#' @title fix.orig.ident
#' @description Remove the string "filtered_feature_bc_matrix." from "orig.ident". Helper function.
#' @param obj Seurat object, Default: merged.obj
#' @examples
#' \dontrun{
#' if(interactive()){
#'  merged.obj$orig.ident <- fix.orig.ident(obj = merged.obj); table(merged.obj$orig.ident)
#'  }
#' }
#' @export
fix.orig.ident <- function(obj = merged.obj) {
  fixed <- sub(obj$'orig.ident', pattern = 'filtered_feature_bc_matrix.', replacement = '')
  return(fixed)
}


# _________________________________________________________________________________________________
#' @title set.all.genes
#' @description It is just a reminder to use calc.q99.Expression.and.set.all.genes to create the all.genes variable
#' @param obj Seurat object, Default: combined.obj
#' @examples
#' \dontrun{
#' if(interactive()){
#'  set.all.genes(); all.genes
#'  }
#' }
#' @export
set.all.genes <- function(obj = combined.obj) iprint("Use calc.q99.Expression.and.set.all.genes()")



# _________________________________________________________________________________________________
#' @title set.mm
#' @description Helps to find metadata columns. It creates a list with the names of of 'obj@meta.data'.
#' @param obj Seurat object, Default: combined.obj
#' @examples
#' \dontrun{
#' if(interactive()){
#'  set.mm(); mm
#'  }
#' }
#' @export
set.mm <- function(obj = combined.obj) {
  mm <- CodeAndRoll2::list.fromNames(colnames(obj@meta.data))
  assign(x = 'mm', value = mm, envir = as.environment(1))
}



# _________________________________________________________________________________________________
#' @title recall.all.genes
#' @description all.genes set by calc.q99.Expression.and.set.all.genes() #
#' @param obj Seurat object, Default: combined.obj
#' @examples
#' \dontrun{
#' if(interactive()){
#'  recall.all.genes(); all.genes
#'  }
#' }
#' @export
recall.all.genes <- function(obj = combined.obj) { # all.genes set by calc.q99.Expression.and.set.all.genes()
  if (!exists('all.genes')) {
    all.genes <- obj@misc$all.genes
    print(head(unlist(all.genes)))
    ww.assign_to_global(name = "all.genes", value = all.genes, verbose = F)
  } else {print("variable 'all.genes' exits in the global namespace")}
}


# _________________________________________________________________________________________________
#' @title recall.meta.tags.n.datasets
#' @description Recall  meta.tags from obj@misc to "meta.tags" in the global environment.
#' @param obj Seurat object, Default: combined.obj
#' @examples
#' \dontrun{
#' if(interactive()){
#'  recall.n.datasets(); n.datasets
#'  }
#' }
#' @export
recall.meta.tags.n.datasets <- function(obj = combined.obj) {
  if (!exists('n.datasets')) {
    n.datasets <- obj@misc$n.datasets
    print(head(unlist(n.datasets)))
    ww.assign_to_global(name = "n.datasets", value = n.datasets)
  } else {print("variable 'n.datasets' exits in the global namespace")}

  if (!exists('meta.tags')) {
    meta.tags <- obj@misc$meta.tags
    print(head(unlist(meta.tags)))
    ww.assign_to_global(name = "meta.tags", value = meta.tags)
  } else {print("variable 'meta.tags' exits in the global namespace")}

}


# _________________________________________________________________________________________________
#' @title recall.parameters
#' @description Recall parameters from obj@misc to "p" in the global environment.
#' @param obj Seurat object, Default: combined.obj
#' @param overwrite PARAM_DESCRIPTION, Default: FALSE
#' @examples
#' \dontrun{
#' if(interactive()){
#'  recall.parameters(); p
#'  }
#' }
#' @export
recall.parameters <- function(obj = combined.obj, overwrite = FALSE) {
  if (is.null(obj@misc$'p')) {
    print("obj does not have: obj@misc$p")
  } else {
    p <- obj@misc$'p'
    print(head(p))
    if (exists('p')) iprint("variable 'p' exits in the global namespace:");

    if (!exists('p') | (exists('p') & overwrite == TRUE) ) {
      ww.assign_to_global(name = "p", value = p); print("Overwritten.")
    } else {
      print("Not overwritten.")
    }
  } # else if obj@misc$'p'
}



# _________________________________________________________________________________________________
#' @title recall.genes.ls
#' @description Recall genes.ls from obj@misc to "genes.ls" in the global environment.
#' @param obj Seurat object, Default: combined.obj
#' @examples
#' \dontrun{
#' if(interactive()){
#'  recall.genes.ls(); genes.ls
#'  }
#' }
#' @export
recall.genes.ls<- function(obj = combined.obj) { # genes.ls
  if (!exists('genes.ls')) {
    genes.ls <- obj@misc$genes.ls
    print(head(unlist(genes.ls)))
    ww.assign_to_global(name = "genes.ls", value = genes.ls)
  } else {print("variable 'genes.ls' exits in the global namespace")}
}



# _________________________________________________________________________________________________
#' @title save.parameters
#' @description Save parameters to obj@misc$p
#' @param obj Seurat object, Default: combined.obj
#' @param params PARAM_DESCRIPTION, Default: p
#' @examples
#' \dontrun{
#' if(interactive()){
#'  save.parameters(obj = combined.obj, params = p);
#'  }
#' }
#' @export
save.parameters <- function(obj = combined.obj, params = p) {
  if (!is.null(obj@misc$'p')) print("Overwriting already existing obj@misc$p. Old version:") ; print(head(unlist(obj@misc$'p')))
  obj@misc$p <- params
}



# _________________________________________________________________________________________________
#' @title plot.expression.rank.q90
#' @description Plot gene expression based on the expression at the 90th quantile (so you will not lose genes expressed in few cells).
#' @param obj Seurat object, Default: combined.obj
#' @param gene gene of interest, Default: 'ACTB'
#' @param filterZero PARAM_DESCRIPTION, Default: T
#' @examples
#' \dontrun{
#' if(interactive()){
#'  plot.expression.rank.q90(gene = "SATB2")
#'  }
#' }
#' @export plot.expression.rank.q90
#' @importFrom Stringendo percentage_formatter
plot.expression.rank.q90 <- function(obj = combined.obj, gene="ACTB", filterZero = T) {
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
    title <- paste(gene, "is in the", Stringendo::percentage_formatter(quantile.GOI), "quantile of 'q90-av' expression. \n There are", counts,"counts" )
  }
  suppressWarnings(
    whist(expr.all, vline = expr.GOI, breaks = 100, main = title, plotname =   make.names(title)
          , ylab = "Genes"
          , xlab = "Av. mRNA in the 10% top expressing cells (q90 av.exp.)")
  )
}




# _________________________________________________________________________________________________
#' @title FlipReductionCoordinates
#' @description Flip reduction coordinates (like UMAP upside down).
#' @param obj Seurat object, Default: combined.obj
#' @param dim Numer of dimensions used, Default: 2
#' @param reduction UMAP, tSNE, or PCA (Dim. reduction to use), Default: 'umap'
#' @param flip PARAM_DESCRIPTION, Default: c("x", "y", "xy", NULL)[1]
#' @param FlipReductionBackupToo PARAM_DESCRIPTION, Default: TRUE
#' @examples
#' \dontrun{
#' if(interactive()){
#'  clUMAP(); combined.obj <- FlipReductionCoordinates(combined.obj); clUMAP()
#'  }
#' }
#' @export
FlipReductionCoordinates <- function(obj = combined.obj, dim = 2, reduction="umap"
                                     , flip = c('x', 'y', 'xy', NULL)[1], FlipReductionBackupToo = TRUE) { # Set active UMAP to `obj@reductions$umap` from `obj@misc$reductions.backup`.
  coordinates <- Embeddings(obj, reduction = reduction)
  stopifnot(ncol(coordinates) == dim )

  if (flip %in% c('x', 'xy')) coordinates[,1] = coordinates[,1] * -1
  if (flip %in% c('y', 'xy')) coordinates[,2] = coordinates[,2] * -1
  obj@reductions[[reduction]]@cell.embeddings <- coordinates

  if (FlipReductionBackupToo) {
    bac.slot <- paste0(reduction,dim,"d")
    if (length(obj@misc$reductions.backup[[bac.slot]])) {
      obj@misc$reductions.backup[[bac.slot]]@cell.embeddings <- coordinates
      iprint(dim, "dimensional",reduction,"backup flipped too.")
    }
  }
  return(obj)
}




# _________________________________________________________________________________________________
#' @title SeuratColorVector
#' @description Recall a Seurat color vector.
#' @param ident identity used, Default: NULL
#' @param obj Seurat object, Default: combined.obj
#' @param plot.colors PARAM_DESCRIPTION, Default: F
#' @param simple Return simply the unique colors, in order? Default: F
#' @examples
#' \dontrun{
#' if(interactive()){
#'  SeuratColorVector(); SeuratColorVector(ident = GetNamedClusteringRuns()[2], plot.colors = T)
#'  }
#' }
#' @seealso
#'  \code{\link[scales]{hue_pal}}
#' @export
#' @importFrom scales hue_pal

SeuratColorVector <- function(ident = NULL, obj = combined.obj, plot.colors = F, simple = F) {
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
  if (simple) {
    colorlevels
  } else {
    translate(vec = as.character(ident.vec)
              , oldvalues = levels(ident.vec)
              , newvalues = colorlevels)
  }
}


# _________________________________________________________________________________________________
#' @title getClusterColors
#' @description get Seurat's cluster colors.
#' @param obj Seurat object, Default: combined.obj
#' @param ident identity used, Default: GetClusteringRuns()[1]
#' @param show PARAM_DESCRIPTION, Default: T
#' @examples
#' \dontrun{
#' if(interactive()){
#'  getClusterColors(obj = combined.obj, ident = GetClusteringRuns()[2] )
#'  }
#' }
#' @seealso
#'  \code{\link[scales]{hue_pal}}
#' @export
#' @importFrom scales hue_pal
getClusterColors <- function(obj = combined.obj
                             , use_new_palettes = TRUE
                             , palette =  c("alphabet", "alphabet2", "glasbey", "polychrome", "stepped")[3]
                             , ident = GetClusteringRuns()[1]
                             , show = T) {
  (identities <- levels(as.factor(obj[[ident]][,1])))
  n.clusters <- length(unique(obj[[ident]][,1]))
  color_palette <- if (use_new_palettes) {
    DiscretePalette(n = n.clusters, palette = palette)
  } else {
    scales::hue_pal()(length(identities))
  }
  # color_check(color_palette)
  # names(color_palette) <- sort(as.factor(identities))
  names(color_palette) <- (identities)
  identvec <- obj[[ident]][,1]
  colz <- color_palette[identvec]
  names(colz) <- identvec
  if (show) color_check(unique(colz)) # MarkdownReports
  colz
}


# _________________________________________________________________________________________________
#' @title get.clustercomposition
#' @description Get cluster composition: which datasets contribute to each cluster?
#' @param obj Seurat object, Default: combined.obj
#' @param x Bars along the X axis, Default: 'integrated_snn_res.0.3'
#' @param y Vertical split of each bar, Default: 'project'
#' @param color Color, Default: y
#' @param plot  Show plot, Default: T
#' @param ScaleTo100pc Scale the Y Axis, Default: T
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @examples get.clustercomposition(); get.clustercomposition()
#' @export


get.clustercomposition <- function(obj = combined.obj, ident = 'integrated_snn_res.0.3', splitby = 'ShortNames'
                                   , color = y
                                   , plot = TRUE, ScaleTo100pc = TRUE
                                   , ...) {
  setwd(OutDir)
  clUMAP(obj = obj, ident = x, save.plot = T, suffix = "as.in.barplot")

  (df.meta <- obj@meta.data[, c(ident, splitby)])

  df.meta %>%
    dplyr::group_by_(splitby) %>%
    summarise()

  categ.per.cluster <- ggbarplot(obj@meta.data
                                 , x = x
                                 , y = y
                                 , color = y
                                 , ...
  )
  if (ScaleTo100pc) categ.per.cluster <- categ.per.cluster + scale_y_discrete(labels = scales::percent_format())
  if (plot) categ.per.cluster

  ggExpress::qqSave(categ.per.cluster, ...)
}



# _________________________________________________________________________________________________
#' remove.residual.small.clusters
#' @description E.g.: after subsetting often some residual cells remain in clusters originally denfined in the full dataset.
#' @param identitites Identities to scan for residual clusters
#' @param obj Seurat object, Default: combined.obj
#' @param max.cells Max number of cells in cluster to be removed. Default: 0.5% of the dataset
#' @export

remove.residual.small.clusters <- function(identitites = GetOrderedClusteringRuns()
                                           , obj = combined.obj
                                           , max.cells = round((ncol(obj))/2000)
                                           ) {
  META <- obj@meta.data
  all.cells <- rownames(META)

  iprint("max.cells:",max.cells,"| Scanning over these", length(identitites), "identities:", identitites)
  small.clusters <- cells.to.remove <- list.fromNames(identitites)

  for (i in 1:length(identitites)) {
    colX <- identitites[i]
    tbl <- table(META[[colX]])

    small.clusters[[i]] <- which_names(tbl <=max.cells )
    cells.to.remove[[i]] <- all.cells[which(META[[colX]] %in% small.clusters)]
    if (length(cells.to.remove[[i]])) iprint(length(cells.to.remove[[i]]), "cells in small clusters:", small.clusters[[i]]) # , head(cells.to.remove[[i]])
  }

  all.cells.2.remove <- unique(unlist(cells.to.remove))
  if (length(all.cells.2.remove)) {
    iprint(">>> a total of", length(all.cells.2.remove), "cells are removed which belonged to a small cluster in any of the identities.")
  } else { iprint(">>> No cells are removed because belonging to small cluster.")}

  cells.2.keep <- setdiff(all.cells, all.cells.2.remove)
  obj <- subset(x = obj, cells = cells.2.keep)

  return(obj)
}


# _________________________________________________________________________________________________
#' drop.levels.Seurat
#' @description Drop levels in clustering vectors in metadata (e.g. after subsetting)
#' @param obj Seurat object, Default: combined.obj
#' @export

drop.levels.Seurat <- function(obj = combined.obj) {
  META <- obj@meta.data
  colclasses <- sapply(META, class)
  drop_in_these <- names(colclasses[ colclasses == 'factor'])
  iprint("Dropping levels in", length(drop_in_these), drop_in_these)
  for (i in 1:length(drop_in_these)) {
    colX <- drop_in_these[i]
    META[[colX]] <- droplevels(META[[colX]])
  }

  obj@meta.data <- META
  return(obj)
}


# _________________________________________________________________________________________________
# Monocle.Utils.R
# ____________________________________________________________________ ----
# source('~/GitHub/Packages/Seurat.utils/Functions/Monocle.Utils.R')
# rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)

# Functions ------------------------
# source('https://raw.githubusercontent.com/vertesy/Seurat.Pipeline/main/elements/Load.packages.local.R')
# try(source("~/GitHub/Packages/Seurat.multicore/00.Load.Seurat3.Multicore.LOCAL.R"));


# _________________________________________________________________________________________________
# _________________________________________________________________________________________________


# mplotGene ------------------------
#' @title mplotGene
#' @description Plot genes in Monocle.
#' @param gene gene of interest, Default: 'PGK1'
#' @param reduction UMAP, tSNE, or PCA (Dim. reduction to use), Default: 'UMAP'
#' @param obj Seurat object, Default: cds_from_seurat
#' @examples
#' \dontrun{
#' if(interactive()){
#'  mplotGene(); mplotGene("GAPDH")
#'  }
#' }
#' @export
mplotGene <- function(gene = "PGK1", reduction = "UMAP", obj = cds_from_seurat) {
  pl1 <- plot_cells(cds = obj,
                    # color_cells_by = 'clusters_low',
                    genes = gene,
                    reduction_method =reduction,
                    # color_cells_by = 'clusters_superlow',
                    label_groups_by_cluster = F,
                    label_leaves = F,
                    label_branch_points = F, label_cell_groups = F
                    , group_label_size = 10
                    , cell_size = 1, alpha = .5)
  ggExpress::qqSave(pl1, fname = ppp(reduction, gene, ".png"), w = 14, h = 7)
}



# mplotManyGenes ------------------------
#' @title mplotManyGenes
#' @description Plot many genes in Monocle.
#' @param ls.genes PARAM_DESCRIPTION, Default: c(`S-phase` = "TOP2A",
#'    `G2M-phase` = "HIST1H4C", oRG = "ID4",
#'    oRG = "HOPX", `Intermediate progenitor` = "EOMES", `Intermediate progenitor1` = "TAC3",
#'    Astroglia = "GFAP", Astrocyte = "S100B", `Immature neurons` = "SLA",
#'    Interneurons = "DLX6-AS1", `Hypoxia/Stress` = "DDIT4", Glycolytic = "PDK1",
#'    `Low-Quality` = "POLR2A", PGC = "DCN", `dl-EN` = "KAZN",
#'    `ul-EN` = "SATB2")
#' @examples
#' \dontrun{
#' if(interactive()){
#'  mplotManyGenes()
#'  }
#' }
#' @export
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



#' @title m3DplotGene
#' @description Plot a gene in 3D in Monocle.
#' @param gene gene of interest, Default: 'PGK1'
#' @param reduction UMAP, tSNE, or PCA (Dim. reduction to use), Default: 'UMAP'
#' @param obj Seurat object, Default: cds.10pc
#' @param ttl.suffix PARAM_DESCRIPTION, Default: 'expression'
#' @param suffix A suffix added to the filename, Default: ''
#' @param cex Point size, Default: 5
#' @examples
#' \dontrun{
#' if(interactive()){
#'  m3DplotGene("DDIT4")
#'  }
#' }
#' @export
m3DplotGene <- function(gene = "PGK1", reduction = "UMAP", obj = cds.10pc, ttl.suffix = "expression"
                        , suffix = "", cex = 5) {
  gene <- intersect(rownames(obj), gene)
  if (length(gene)) {
    pl1 <- plot_cells_3d(cds = obj
                         # color_cells_by = 'clusters_low',
                         , genes = gene
                         , reduction_method = reduction
                         # color_cells_by = 'clusters_superlow',
                         # label_groups_by_cluster = F,
                         # label_leaves = F,
                         # label_branch_points = F, label_cell_groups = F
                         # , group_label_size = 10
                         , cell_size = cex
                         , alpha = .5) %>%
      plotly::layout(title = paste0(gene, ttl.suffix))
    SavePlotlyAsHtml(pl1, category. = gene, suffix. = suffix)
  } else { iprint("gene not found") }
}



# m3DplotKeyGenes ------------------------
#' @title m3DplotKeyGenes
#' @description Plot many genes in 3D in Monocle.
#' @param obj Seurat object, Default: cds.10pc
#' @param cex Point size, Default: iround(log10(idim(obj)[2]))
#' @param ls.genes PARAM_DESCRIPTION, Default: c(`S-phase` = "TOP2A", `G2M-phase` = "HIST1H4C", oRG = "ID4",
#'    oRG = "HOPX", `Intermediate progenitor` = "EOMES", `Intermediate progenitor1` = "TAC3",
#'    Astroglia = "GFAP", Astrocyte = "S100B", `Immature neurons` = "SLA",
#'    Interneurons = "DLX6-AS1", `Hypoxia/Stress` = "DDIT4", Glycolytic = "PDK1",
#'    `Low-Quality` = "POLR2A", PGC = "DCN", `dl-EN` = "KAZN",
#'    `ul-EN` = "SATB2")
#' @param reduction UMAP, tSNE, or PCA (Dim. reduction to use), Default: 'UMAP'
#' @param suffix A suffix added to the filename, Default: ''
#' @examples
#' \dontrun{
#' if(interactive()){
#'  m3DplotKeyGenes()
#'  }
#' }
#' @export
m3DplotKeyGenes <- function(obj = cds.10pc, cex = iround(log10(idim(obj)[2]))
                            , ls.genes =  c(
                              `S-phase` = "TOP2A", `G2M-phase` = "HIST1H4C"
                              , `oRG` = "ID4", `oRG` = "HOPX"
                              , `Intermediate progenitor` = "EOMES"
                              ,  `Intermediate progenitor1` = "TAC3"
                              , Astroglia = "GFAP", Astrocyte = "S100B"
                              , `Immature neurons` = "SLA", Interneurons = "DLX6-AS1"
                              , `Hypoxia/Stress` = "DDIT4", Glycolytic = "PDK1"
                              , `Low-Quality` = "POLR2A", `PGC` = "DCN"
                              , `dl-EN` = "KAZN", `ul-EN` = "SATB2")
                            , reduction = "UMAP", suffix = "") {

  ls.genes <- intersect(rownames(obj), ls.genes)
  iprint(length(ls.genes), "genes found.")
  create_set_SubDir(ppp("3D.gex.plots", substitute(obj)))
  for (g in ls.genes) {
    m3DplotGene(gene = g, obj = obj, cex = cex)
  }
  create_set_OutDir(ParentDir)
}


# subsetMonocleObject ------------------------
#' @title subsetMonocleObject
#' @description Subset a compressed Seurat Obj and save it in wd.
#' @param obj Seurat object, Default: cds_from_seurat
#' @param fraction_ PARAM_DESCRIPTION, Default: 0.1
#' @param nCells PARAM_DESCRIPTION, Default: F
#' @param seed_ PARAM_DESCRIPTION, Default: 1989
#' @examples
#' \dontrun{
#' if(interactive()){
#'  cds.10pc <- subsetMonocleObject(cds_from_seurat, fraction_ = .1);
#'  cds.25pc <- subsetMonocleObject(cds_from_seurat, fraction_ = .25)
#'  }
#' }
#' @export
#' @importFrom Stringendo percentage_formatter
subsetMonocleObject <- function(obj = cds_from_seurat, fraction_ = 0.1, nCells = F, seed_ = 1989 ) { # Subset a compressed Seurat Obj and save it in wd.
  set.seed(seed_)
  if (isFALSE(nCells)) {
    all.cells <- colnames(obj)
    n.keep <- floor(length(all.cells) * fraction_)
    cellIDs.keep <- sample(all.cells, size = n.keep, replace = FALSE)
    iprint(length(cellIDs.keep), "or",Stringendo::percentage_formatter(fraction_),"of the cells are kept. Seed:", head(cellIDs.keep), seed_)
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




# m3.get.umap ------------------------
#' @title m3.get.umap
#' @description Fetch the umap coordinates from obj@int_colData@listData$reducedDims[[slot]]
#' @param obj Seurat object, Default: cds_from_seurat
#' @param slot slot in the Seurat object. Default: 'UMAP'
#' @param dim Numer of dimensions used, Default: (2:3)[2]
#' @examples
#' \dontrun{
#' if(interactive()){
#'  m3.get.umap(obj = cds_from_seurat, slot = 'UMAP', dim = (2:3)[2])
#'  }
#' }
#' @export
m3.get.umap <- function(obj = cds_from_seurat, slot = 'UMAP', dim = (2:3)[2]) {
  iprint(dim, 'dimensional', slot)
  df.reduction <- obj@int_colData@listData$reducedDims[[slot]]
  stopif(condition = is.null(df.reduction), message = paste('slot not found', slot))
  stopifnot(ncol(df.reduction) == dim)
  return(df.reduction)
}


# _________________________________________________________________________________________________
#' @title m3.backup.umap
#' @description Backup umap coordinates to obj@int_colData@listData$reducedDims[[new.slot]]
#' @param obj Seurat object, Default: cds_from_seurat
#' @param slot slot in the Seurat object. Default: 'UMAP'
#' @param dim Numer of dimensions used, Default: (2:3)[2]
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @examples
#' \dontrun{
#' if(interactive()){
#'  cds_from_seurat <- m3.backup.umap(obj = cds_from_seurat, slot = 'UMAP', dim = (2:3)[2])
#'  }
#' }
#' @export
m3.backup.umap <- function(obj = cds_from_seurat, slot = 'UMAP', dim = (2:3)[2], ...) {
  new.slot <- paste0(slot,'.',dim, 'D')
  obj@int_colData@listData$reducedDims[[new.slot]] <- m3.get.umap(obj = obj, slot = slot, dim = dim, ... )
  iprint('obj@int_colData@listData$reducedDims$', new.slot)
  return(obj)
}



# _________________________________________________________________________________________________
#' @title m3.recall.umap
#' @description Fetch UMAP coordinates.
#' @param obj Seurat object, Default: cds_from_seurat
#' @param slot slot in the Seurat object. Default: 'UMAP'
#' @param dim Numer of dimensions used, Default: (2:3)[2]
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @examples
#' \dontrun{
#' if(interactive()){
#'  cds_from_seurat <- m3.recall.umap(obj = cds_from_seurat, slot = 'UMAP', dim = (2:3)[2]); cds_from_seurat.bac <- cds_from_seurat
#'  }
#' }
#' @export
m3.recall.umap <- function(obj = cds_from_seurat, slot = 'UMAP', dim = (2:3)[2], ...) {
  backup.slot <- paste0(slot,'.',dim, 'D')
  old.dim <- ncol(obj@int_colData@listData$reducedDims[[slot]])
  obj@int_colData@listData$reducedDims[[slot]] <- obj@int_colData@listData$reducedDims[[backup.slot]]
  iprint(old.dim, 'dimensional', slot, 'replaced by', dim, slot, 'reduction.')
  return(obj)
}



# _________________________________________________________________________________________________
#' @title m3.export.umap.2.Seurat
#' @description Export umap coordinates.
#' @param mobj PARAM_DESCRIPTION, Default: cds_from_seurat
#' @param sobj PARAM_DESCRIPTION, Default: combined.obj
#' @param def.dim PARAM_DESCRIPTION, Default: 2
#' @examples
#' \dontrun{
#' if(interactive()){
#'  combined.obj.bac <- combined.obj; # combined.obj <- combined.obj.bac
#'  combined.obj <- m3.export.umap.2.Seurat(mobj = cds_from_seurat, sobj = combined.obj); qUMAP("DDIT4")
#'  }
#' }
#' @export
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

# _________________________________________________________________________________________________






# _________________________________________________________________________________________________
# MULTI-seq.functions.R ----
# ____________________________________________________________________ ----
# source('~/GitHub/Packages/Seurat.utils/Functions/MULTI-seq.functions.R')
# try(source('https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/MULTI-seq.functions.R'))

# Requirements __________________________________________
# try(require(MarkdownReports),  silent = T)
# try(require(pheatmap),  silent = T)
# May also require


# BarTableSweepList -----------------------------------------------------------------------
#' @title BarTableSweepList
#' @description BarTableSweepList
#' @param min PARAM_DESCRIPTION, Default: 0.01
#' @param max PARAM_DESCRIPTION, Default: 0.99
#' @param step PARAM_DESCRIPTION, Default: 0.02
#' @param bar_table PARAM_DESCRIPTION, Default: bar.table
#' @examples
#' \dontrun{
#' if(interactive()){
#'  bar.table_sweep.list <- BarTableSweepList(bar_table = bar.table.solo) # Quantile Sweep List
#'  }
#' }
#' @export
BarTableSweepList <- function(min = 0.01, max = 0.99, step = 0.02, bar_table =bar.table) {
  bar.table_sweep.list <- list()
  n <- 0
  Quantiles = seq(from = min, to = max, by = step)
  for (n in 1:length(Quantiles)) { # print(q)
    bar.table_sweep.list[[n]] <- classifyCells(bar_table, q = Quantiles[n])
    names(bar.table_sweep.list)[n] <- paste("q=",Quantiles[n], sep="")
  }
  return(bar.table_sweep.list)
}





# _________________________________________________________________________________________________
#' @title mSeq.map.all96.BCs
#' @description mSeq.map.all96.BCs
#' @param readTable PARAM_DESCRIPTION, Default: readTable
#' @param CellIDs PARAM_DESCRIPTION, Default: CellIDs
#' @param path2allBCs PARAM_DESCRIPTION, Default: '~/Google_Drive/Science/IMBA/MULTI.seq/from.US/All.MULTI-seq_barcodes.Mar2019.tsv'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  bar.table <-mSeq.map.all96.BCs(readTable = readTable, CellIDs = CellIDs)
#'  }
#' }
#' @export
mSeq.map.all96.BCs <- function(readTable = readTable, CellIDs = CellIDs
                               , path2allBCs = '~/Google_Drive/Science/IMBA/MULTI.seq/from.US/All.MULTI-seq_barcodes.Mar2019.tsv'
) {
  (bar.ref <- read_tsv(path2allBCs)[[1]]) # Vector of reference all MULTI-seq sample barcode sequences.
  MULTIseq.align(readTable = readTable, cellIDs = CellIDs, ref = bar.ref)
}



# _________________________________________________________________________________________________

#' @title aux_plotAllMseqBCs
#' @description aux_plotAllMseqBCs
#' @param bar.table PARAM_DESCRIPTION, Default: bar.table[, 1:96]
#' @param barcodes.used PARAM_DESCRIPTION, Default: BCs.used
#' @param plotname Title of the plot, Default: 'Barcode seq depth'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  aux_plotAllMseqBCs(bar.table = bar.table[,1:96], barcodes.used = BCs.used, plotname = "Barcode seq depth")
#'  }
#' }
#' @export
aux_plotAllMseqBCs <- function(bar.table = bar.table[,1:96], barcodes.used = BCs.used
                               , plotname = "Barcode seq depth") {
  stopifnot(is.numeric(BCs.used))
  BC.depth <- colSums(bar.table)[1:96]
  if (min(BC.depth) < 1) { BC.depth <- BC.depth+1 }
  log.depth <- log10(BC.depth); range(log.depth)
  ccc <- colnames(bar.table) %in% BCs.used

  wbarplot(log.depth, col = ccc, lwd = 1, lcol = 1, lty = 2, plotname = plotname
           , vline = (range(BCs.used)+c(-1,1))
           , hline = quantile(log.depth[setdiff(1:69, BCs.used)], .95)
           , ylab = "log10(total reads / BC)", main =  plotname)
  wlegend.label(
    "    Horiz. line at 95%
    of unused BC's.
    Vertical line: range
    of used BCs.",poz = 2, cex = 1)
}


# _________________________________________________________________________________________________
# bar.table.log <- t(log10(bar.table[,BCs.used]+1))
# bar.table.log <- CodeAndRoll2::clip.outliers.at.percentile(bar.table.log)



# ____________________________________________________________________ ----
# plotting.dim.reduction.2D.R ----
# _________________________________________________________________________________________________
# source('~/GitHub/Packages/Seurat.utils/Functions/plotting.dim.reduction.2D.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Plotting.dim.reduction.2D.R"))
# Source: self + web

# Requirements __________________________________________
# library(plotly)
# try(source("~/GitHub/Packages/ggExpressDev/ggExpress.functions.R"), silent = T)
# try(source("https://raw.githubusercontent.com/vertesy/ggExpressDev/main/ggExpress.functions.R"), silent = T)

# May also require
# try (source('/GitHub/Packages/CodeAndRoll/CodeAndRoll.R'),silent= F) # generic utilities funtions
# require('MarkdownReports') # require("devtools") # plotting related utilities functions # devtools::install_github(repo = "vertesy/MarkdownReports")


# _________________________________________________________________________________________________
#' @title qUMAP
#' @description The quickest way to draw a gene expression UMAP #
#' @param feature PARAM_DESCRIPTION, Default: 'TOP2A'
#' @param obj Seurat object, Default: combined.obj
#' @param title Title of the plot, Default: feature
#' @param sub Subtitle of the plot, Default: NULL
#' @param reduction UMAP, tSNE, or PCA (Dim. reduction to use), Default: 'umap'
#' @param splitby PARAM_DESCRIPTION, Default: NULL
#' @param prefix A prefix added before the filename, Default: NULL
#' @param suffix A suffix added to the end of the filename, Default: subtitle
#' @param save.plot Save the plot into a file?, Default: T
#' @param PNG Save file as .png?, Default: T
#' @param h height of the plot, Default: 7
#' @param w width of the plot, Default: NULL
#' @param nr.cols Number of columns to combine multiple feature plots to, ignored if split.by is not NULL, Default: NULL
#' @param assay RNA or integrated assay, Default: c("RNA", "integrated")[1]
#' @param HGNC.lookup PARAM_DESCRIPTION, Default: TRUE
#' @param axes Show axes? Default: T
#' @param aspect.ratio Fix the aspect ratio?  Default: Yes, at 0.6
#' @param make.uppercase PARAM_DESCRIPTION, Default: TRUE
#' @param qlow PARAM_DESCRIPTION, Default: 'q10'
#' @param qhigh PARAM_DESCRIPTION, Default: 'q90'
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @examples
#' \dontrun{
#' if(interactive()){
#'  qUMAP('nFeature_RNA'); qUMAP('TOP2A')
#'  }
#' }
#' @export

qUMAP <- function( feature= 'TOP2A', obj =  combined.obj  # The quickest way to draw a gene expression UMAP
                   , title = feature, sub =NULL
                   , reduction ="umap", splitby = NULL
                   , prefix = NULL
                   , suffix = make.names(sub)
                   , save.plot = MarkdownHelpers::TRUE.unless('b.save.wplots')
                   , PNG = T
                   , h = 7, w = NULL, nr.cols = NULL
                   , assay = c("RNA","integrated")[1]
                   , axes = T
                   , aspect.ratio = c(FALSE, 0.6)[2]
                   , HGNC.lookup= TRUE
                   , make.uppercase = FALSE
                   , qlow = "q10", qhigh = "q90", ...) {

  if ( !(feature %in% colnames(obj@meta.data) | feature %in% rownames(obj) ) ) {
    feature <- check.genes(list.of.genes = feature, obj = obj, verbose = F
                           , HGNC.lookup = HGNC.lookup, makeuppercase = make.uppercase)
  }

  DefaultAssay(obj) <- assay
  ggplot.obj <- Seurat::FeaturePlot(obj, features = feature
                            , reduction = reduction
                            , min.cutoff = qlow, max.cutoff = qhigh
                            # , plotname = ppp(toupper(reduction), feature)
                            , ncol = nr.cols
                            , split.by = splitby
                            , ...) +
    ggtitle(label = title, subtitle = sub) +
    if (!axes) NoAxes() else NULL

  ggplot.obj <- ggplot.obj + if (aspect.ratio) ggplot2::coord_fixed(ratio = aspect.ratio) else NULL

  if (save.plot) {
    fname = ww.FnP_parser(Stringendo::sppp(prefix, toupper(reduction), feature, assay, suffix), if (PNG) "png" else "pdf")
    try(save_plot(filename = fname, plot = ggplot.obj, base_height = h, base_width = w)) #, ncol = 1, nrow = 1
  }
  return(ggplot.obj)
}



# _________________________________________________________________________________________________
#' @title clUMAP
#' @description The quickest way to draw a clustering result  UMAP #
#' @param ident identity used, Default: 'integrated_snn_res.0.5'
#' @param obj Seurat object, Default: combined.obj
#' @param reduction UMAP, tSNE, or PCA (Dim. reduction to use), Default: 'umap'
#' @param splitby PARAM_DESCRIPTION, Default: NULL
#' @param title Title of the plot, Default: ident
#' @param sub Subtitle of the plot, Default: NULL
#' @param prefix A prefix added before the filename, Default: NULL
#' @param suffix A suffix added to the filename, Default: sub
#' @param label.cex PARAM_DESCRIPTION, Default: 7
#' @param h height of the plot, Default: 7
#' @param w width of the plot, Default: NULL
#' @param nr.cols Number of columns to combine multiple feature plots to, ignored if split.by is not NULL, Default: NULL
#' @param plotname Title of the plot, Default: ppp(toupper(reduction), ident)
#' @param cols Colors used, Default: getDiscretePalette(ident.used = ident, show.colors = F)
#' @param palette Color palette.
#' @param highlight.clusters PARAM_DESCRIPTION, Default: NULL
#' @param cells.highlight PARAM_DESCRIPTION, Default: NULL
#' @param label PARAM_DESCRIPTION, Default: T
#' @param repel PARAM_DESCRIPTION, Default: T
#' @param legend PARAM_DESCRIPTION, Default: !label
#' @param axes Show axes? Default: T
#' @param aspect.ratio Fix the aspect ratio? Default: Yes, at 0.6
#' @param MaxCategThrHP PARAM_DESCRIPTION, Default: 200
#' @param save.plot Save the plot into a file?, Default: T
#' @param PNG Save file as .png?, Default: T
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @examples
#' \dontrun{
#' if(interactive()){
#'  clUMAP('cl.names.KnownMarkers.0.5' ); clUMAP('cl.names.KnownMarkers.0.5', cols = NULL)
#'  }
#' }
#' @export
# #' @param save.object PARAM_DESCRIPTION, Default: F

clUMAP <- function(ident = "integrated_snn_res.0.5", obj =  combined.obj   # The quickest way to draw a clustering result  UMAP
                   , reduction ="umap", splitby = NULL
                   , title = ident, sub =NULL
                   , prefix = NULL
                   , suffix = make.names(sub)
                   , label.cex = 7
                   , h = 7, w = NULL, nr.cols = NULL
                   , plotname = ppp(toupper(reduction), ident)
                   , cols = NULL, palette =  c("alphabet", "alphabet2", "glasbey", "polychrome", "stepped")[3]
                   , highlight.clusters = NULL, cells.highlight = NULL
                   , label = T, repel = T, legend = !label, MaxCategThrHP = 200
                   , axes = T
                   , aspect.ratio = c(FALSE, 0.6)[2]
                   , save.plot = MarkdownHelpers::TRUE.unless('b.save.wplots')
                   , PNG = T
                   # , save.object = F
                   , ...) {

  IdentFound <- (ident %in%  colnames(obj@meta.data))
  if (!IdentFound) {
    ident <- GetClusteringRuns(obj = obj, pat = "_res.*[0,1]\\.[0-9]$")[1]
    iprint("Identity not found. Plotting", ident)
  }
  identity <- obj[[ident]]
  NtCategs <- length(unique(identity[,1]))

  if ( !missing(highlight.clusters)) {
    idx.ok <- identity[,1] %in% highlight.clusters
    highlight.these <- rownames(identity)[idx.ok]
  } else { highlight.these <- NULL}
  if ( !missing(cells.highlight)) { highlight.these <- cells.highlight } # overwrite, if directly defined


  if (is.null(cols)) {
    cols = if (NtCategs > 5) getDiscretePalette(ident.used = ident, palette.used = palette, obj = obj, show.colors = F)
  }

  if( NtCategs > MaxCategThrHP ) {
    iprint("Too many categories (",NtCategs,") in ", ident, "- use qUMAP for continous variables.")
  } else {
    if( length(unique(identity)) < MaxCategThrHP )
      ggplot.obj <-
        Seurat::DimPlot(object = obj, group.by = ident
                , cols = cols
                , reduction = reduction, split.by = splitby
                , ncol = nr.cols, cells.highlight = highlight.these
                , label = label, repel = repel, label.size = label.cex, ...) +
        ggtitle(label = title, subtitle = sub) +
        if (!legend) NoLegend() else NULL

    ggplot.obj <- ggplot.obj + if (!axes) NoAxes() else NULL
    ggplot.obj <- ggplot.obj + if (aspect.ratio) ggplot2::coord_fixed(ratio = aspect.ratio) else NULL

    if (save.plot) {
      pname = Stringendo::sppp(prefix, plotname, suffix, sppp(highlight.clusters))
      fname = ww.FnP_parser(pname, if (PNG) "png" else "pdf")
      try(save_plot(filename = fname, plot = ggplot.obj, base_height = h, base_width = w)) #, ncol = 1, nrow = 1
    }
    # if(save.object) saveRDS(object = ggplot.obj, file = ppp(fname, 'ggobj.RDS'))
    return(ggplot.obj)
  } # if not too many categories
}




# _________________________________________________________________________________________________
#' @title gg_color_hue
#' @description reproduce the ggplot2 default color palette #
#' @param n PARAM_DESCRIPTION
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
gg_color_hue <- function(n) { # reproduce the ggplot2 default color palette
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
# https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette

# _________________________________________________________________________________________________
#' @title save2umaps.A4
#' @description Save 2 umaps on 1 A4
#' @param plot_list PARAM_DESCRIPTION
#' @param pname PARAM_DESCRIPTION, Default: F
#' @param suffix A suffix added to the filename, Default: NULL
#' @param scale PARAM_DESCRIPTION, Default: 1
#' @param nrow PARAM_DESCRIPTION, Default: 2
#' @param ncol PARAM_DESCRIPTION, Default: 1
#' @param h height of the plot, Default: hA4 * scale
#' @param w width of the plot, Default: wA4 * scale
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
save2umaps.A4 <- function(plot_list, pname = F, suffix = NULL, scale = 1
                          , nrow = 2, ncol = 1
                          , h = hA4 * scale, w = wA4 * scale, ...) { # Save 2 umaps on an A4 page.
  if (pname ==F) pname = Stringendo::sppp(substitute(plot_list), suffix)
  p1 = plot_grid(plotlist = plot_list, nrow = nrow, ncol = ncol, labels = LETTERS[1:length(plot_list)], ...  )
  save_plot(plot = p1, filename = extPNG(pname), base_height = h, base_width = w)
}

# _________________________________________________________________________________________________
#' @title save4umaps.A4
#' @description Save 4 umaps on 1 A4
#' @param plot_list A list of ggplot objects, each of which is one panel.
#' @param pname Plotname, Default: F
#' @param suffix A suffix added to the filename, Default: NULL
#' @param scale Scaling factor of the canvas, Default: 1
#' @param nrow number of rows for panelson the page, Default: 2
#' @param ncol number of columns for panelson the page, Default: 2
#' @param h height of the plot, Default: wA4 * scale
#' @param w width of the plot, Default: hA4 * scale
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
save4umaps.A4 <- function(plot_list, pname = F, suffix = NULL, scale = 1
                          , nrow = 2, ncol = 2
                          , h = wA4 * scale, w = hA4 * scale, ...) { # Save 4 umaps on an A4 page.
  if (pname==F) pname =  Stringendo::sppp(substitute(plot_list), suffix)
  p1 = plot_grid(plotlist = plot_list, nrow = nrow, ncol = ncol, labels = LETTERS[1:length(plot_list)], ...  )
  save_plot(plot = p1, filename = extPNG(pname), base_height = h, base_width = w)
}


# _________________________________________________________________________________________________
#' @title umapNamedClusters
#' @description Plot and save umap based on a metadata column. #
#' @param obj Seurat object, Default: combined.obj
#' @param metaD.colname PARAM_DESCRIPTION, Default: metaD.colname.labeled
#' @param ext File extension for saving, Default: 'png'
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @examples
#' \dontrun{
#' if(interactive()){
#'  umapNamedClusters(obj = combined.obj, metaD.colname = metaD.colname.labeled)
#'  }
#' }
#' @export
umapNamedClusters <- function(obj = combined.obj, metaD.colname = metaD.colname.labeled, ext = "png", ...) { # Plot and save umap based on a metadata column.
  fname = ppp("Named.clusters", metaD.colname, ext)
  p.named =
    Seurat::DimPlot(obj, reduction = "umap", group.by = metaD.colname, label = T, ...) +
    NoLegend() +
    ggtitle(metaD.colname)
  save_plot(p.named, filename = fname); p.named
}




# _________________________________________________________________________________________________

# _________________________________________________________________________________________________
#' @title qqSaveGridA4
#' @description Save 2 or 4 ggplot objects using plot_grid() on an A4 page #
#' @param plotlist PARAM_DESCRIPTION, Default: pl
#' @param plots PARAM_DESCRIPTION, Default: 1:2
#' @param NrPlots PARAM_DESCRIPTION, Default: length(plots)
#' @param height PARAM_DESCRIPTION, Default: hA4
#' @param width PARAM_DESCRIPTION, Default: wA4
#' @param fname File name, Default: 'Fractions.Organoid-to-organoid variation.png'
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @examples
#' \dontrun{
#' if(interactive()){
#'  qqSaveGridA4(plotlist= pl, plots = 1:2, fname = "Fractions.per.Cl.png"); qqSaveGridA4(plotlist= pl, plots = 1:4, fname = "Fractions.per.Cl.4.png")
#'  }
#' }
#' @export
qqSaveGridA4 <- function(plotlist= pl # Save 2 or 4 ggplot objects using plot_grid() on an A4 page
                         , plots = 1:2, NrPlots = length(plots), height = hA4, width = wA4
                         , fname = "Fractions.Organoid-to-organoid variation.png", ...) {
  stopifnot(NrPlots %in% c(2,4))
  iprint(NrPlots,"plots found,", plots,"are saved.")
  pg.cf = plot_grid(plotlist = plotlist[plots], nrow = 2, ncol = NrPlots/2, labels = LETTERS[1:NrPlots], ...  )
  if (NrPlots == 4) list2env(list(height = width, width = height), envir = as.environment(environment()))
  save_plot(filename = fname,
            plot = pg.cf, base_height = height, base_width = width)
  ww.FnP_parser(fname)
}





# _________________________________________________________________________________________________

# umapHiLightSel highlight a set of cells based on clusterIDs provided---------------
#' @title umapHiLightSel
#' @description Highlight a set of cells based on clusterIDs provided. #
#' @param obj Seurat object, Default: combined.obj
#' @param COI PARAM_DESCRIPTION, Default: c("0", "2", "4", "5", "11")
#' @param res.cl PARAM_DESCRIPTION, Default: 'integrated_snn_res.0.3'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
umapHiLightSel <- function(obj = combined.obj, # Highlight a set of cells based on clusterIDs provided.
                           COI =  c("0", "2", "4", "5",  "11"), res.cl = 'integrated_snn_res.0.3') {
  cellsSel = getCellIDs.from.meta(obj, values = COI, ColName.meta = res.cl)
  Seurat::DimPlot(obj, reduction = "umap", group.by = res.cl,
          label = T, cells.highlight = cellsSel)
  ggsave(filename = extPNG(kollapse("cells",COI, collapseby = '.')))
}




# Save multiple FeaturePlot from a list of genes on A4 jpeg ------------------------
#' @title multiFeaturePlot.A4
#' @description Save multiple FeaturePlots, as jpeg, on A4 for each gene, which are stored as a list of gene names. #
#' @param list.of.genes PARAM_DESCRIPTION
#' @param obj Seurat object, Default: combined.obj
#' @param foldername PARAM_DESCRIPTION, Default: substitute(list.of.genes)
#' @param plot.reduction PARAM_DESCRIPTION, Default: 'umap'
#' @param intersectionAssay PARAM_DESCRIPTION, Default: c("RNA", "integrated")[1]
#' @param layout PARAM_DESCRIPTION, Default: c("tall", "wide", FALSE)[2]
#' @param colors PARAM_DESCRIPTION, Default: c("grey", "red")
#' @param nr.Col PARAM_DESCRIPTION, Default: 2
#' @param nr.Row PARAM_DESCRIPTION, Default: 4
#' @param cex Point size, Default: round(0.1/(nr.Col * nr.Row), digits = 2)
#' @param gene.min.exp PARAM_DESCRIPTION, Default: 'q01'
#' @param gene.max.exp PARAM_DESCRIPTION, Default: 'q99'
#' @param subdir PARAM_DESCRIPTION, Default: T
#' @param prefix PARAM_DESCRIPTION, Default: NULL
#' @param suffix A suffix added to the filename, Default: NULL
#' @param background_col background color def: White
#' @param saveGeneList PARAM_DESCRIPTION, Default: FALSE
#' @param w width of the plot, Default: wA4
#' @param h height of the plot, Default: hA4
#' @param scaling PARAM_DESCRIPTION, Default: 1
#' @param aspect.ratio Fix the aspect ratio? Default: Yes, at 0.6
#' @param format PARAM_DESCRIPTION, Default: c("jpg", "pdf", "png")[1]
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[tictoc]{tic}}
#'  \code{\link[cowplot]{plot_grid}}
#' @export
#' @importFrom tictoc tic toc
#' @importFrom cowplot plot_grid


multiFeaturePlot.A4 <- function(list.of.genes # Save multiple FeaturePlots, as jpeg, on A4 for each gene, which are stored as a list of gene names.
                                , obj = combined.obj
                                , foldername = substitute(list.of.genes), plot.reduction='umap'
                                , intersectionAssay = c('RNA', 'integrated')[1]
                                , layout = c('tall', 'wide', FALSE )[2]
                                , colors = c("grey", "red"), nr.Col = 2, nr.Row =4, cex = round(0.1/(nr.Col*nr.Row), digits = 2)
                                , gene.min.exp = 'q01', gene.max.exp = 'q99', subdir =T
                                , prefix = NULL , suffix = NULL
                                , background_col = "white"
                                , aspect.ratio = c(FALSE, 0.6)[2]
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

  lsG = iterBy.over(1:length(list.of.genes.found), by = nr.Row * nr.Col)
  for (i in 1:length(lsG)) {
    genes = list.of.genes.found[lsG[[i]]]
    iprint(i,genes )
    plotname = kpp(c(prefix, plot.reduction,i, genes, suffix, format ))

    plot.list = Seurat::FeaturePlot(object = obj, features = genes, reduction = plot.reduction, combine = F
                            , ncol = nr.Col, cols = colors
                            , min.cutoff = gene.min.exp, max.cutoff = gene.max.exp
                            , pt.size = cex, ...)

    for (i in 1:length(plot.list)) {
      plot.list[[i]] <- plot.list[[i]] + NoLegend() + NoAxes()
      if (aspect.ratio) plot.list[[i]] <- plot.list[[i]] + ggplot2::coord_fixed(ratio = aspect.ratio)
    }

    pltGrid <- cowplot::plot_grid(plotlist = plot.list, ncol = nr.Col, nrow = nr.Row )
    ggsave(filename = plotname, width = w, height = h, bg = background_col, plot = pltGrid)
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

#' @title multiFeatureHeatmap.A4
#' @description Save multiple FeatureHeatmaps from a list of genes on A4 jpeg #
#' @param obj Seurat object, Default: combined.obj
#' @param list.of.genes PARAM_DESCRIPTION
#' @param gene.per.page PARAM_DESCRIPTION, Default: 5
#' @param group.cells.by PARAM_DESCRIPTION, Default: 'batch'
#' @param plot.reduction PARAM_DESCRIPTION, Default: 'umap'
#' @param cex Point size, Default: iround(3/gene.per.page)
#' @param sep_scale PARAM_DESCRIPTION, Default: F
#' @param gene.min.exp PARAM_DESCRIPTION, Default: 'q5'
#' @param gene.max.exp PARAM_DESCRIPTION, Default: 'q95'
#' @param jpeg.res PARAM_DESCRIPTION, Default: 225
#' @param jpeg.q PARAM_DESCRIPTION, Default: 90
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[tictoc]{tic}}
#' @export
#' @importFrom tictoc tic toc
multiFeatureHeatmap.A4 <- function(obj = combined.obj # Save multiple FeatureHeatmaps from a list of genes on A4 jpeg
                                   , list.of.genes, gene.per.page = 5
                                   , group.cells.by= "batch", plot.reduction='umap'
                                   , cex = iround(3/gene.per.page), sep_scale = F
                                   , gene.min.exp = 'q5', gene.max.exp = 'q95'
                                   , jpeg.res = 225, jpeg.q = 90, ...) {

  tictoc::tic()
  list.of.genes = check.genes(list.of.genes, obj = obj)

  lsG = iterBy.over(1:length(list.of.genes), by = gene.per.page)
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


# __________________________________________

#' @title plot.UMAP.tSNE.sidebyside
#' @description Plot a UMAP and tSNE sidebyside #
#' @param obj Seurat object, Default: combined.obj
#' @param grouping PARAM_DESCRIPTION, Default: 'res.0.6'
#' @param no_legend PARAM_DESCRIPTION, Default: F
#' @param do_return PARAM_DESCRIPTION, Default: TRUE
#' @param do_label PARAM_DESCRIPTION, Default: T
#' @param label_size PARAM_DESCRIPTION, Default: 10
#' @param vector_friendly PARAM_DESCRIPTION, Default: TRUE
#' @param cells_use PARAM_DESCRIPTION, Default: NULL
#' @param no_axes PARAM_DESCRIPTION, Default: T
#' @param pt_size PARAM_DESCRIPTION, Default: 0.5
#' @param name.suffix PARAM_DESCRIPTION, Default: NULL
#' @param width PARAM_DESCRIPTION, Default: hA4
#' @param heigth PARAM_DESCRIPTION, Default: 1.75 * wA4
#' @param filetype PARAM_DESCRIPTION, Default: 'pdf'
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[cowplot]{save_plot}}
#' @export plot.UMAP.tSNE.sidebyside
#' @importFrom cowplot save_plot
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

  p1 <- Seurat::DimPlot(object = obj, reduction.use = "tsne", no.axes = no_axes, cells.use = cells_use
                , no.legend = no_legend, do.return = do_return, do.label = do_label, label.size = label_size
                , group.by = grouping, vector.friendly = vector_friendly, pt.size = pt_size, ...) +
    ggtitle("tSNE") + theme(plot.title = element_text(hjust = 0.5))

  p2 <- Seurat::DimPlot(object = obj, reduction.use = "umap", no.axes = no_axes, cells.use = cells_use
                , no.legend = T, do.return = do_return, do.label = do_label, label.size = label_size
                , group.by = grouping, vector.friendly = vector_friendly, pt.size = pt_size, ...) +
    ggtitle("UMAP") + theme(plot.title = element_text(hjust = 0.5))

  plots = plot_grid(p1, p2, labels = c("A", "B"), ncol = 2)
  plotname = kpp( 'UMAP.tSNE', grouping, name.suffix, filetype)

  cowplot::save_plot(filename = plotname, plot = plots
                     , ncol = 2 # we're saving a grid plot of 2 columns
                     , nrow = 1 # and 2 rows
                     , base_width = width
                     , base_height = heigth
                     # each individual subplot should have an aspect ratio of 1.3
                     # , base_aspect_ratio = 1.5
  )
}

# _________________________________________________________________________________________________
#' @title PlotTopGenesPerCluster
#' @description Plot the top N diff. exp. genes in each cluster
#' @param obj Seurat object, Default: combined.obj
#' @param cl_res PARAM_DESCRIPTION, Default: res
#' @param nrGenes PARAM_DESCRIPTION, Default: p$n.markers
#' @param order.by Sort output tibble by which column, Default: c("combined.score", "avg_log2FC", "p_val_adj")[1]
#' @param df_markers PARAM_DESCRIPTION, Default: combined.obj@misc$df.markers[[paste0("res.", cl_res)]]
#' @examples
#' \dontrun{
#' if(interactive()){
#'  PlotTopGenesPerCluster(obj = combined.obj, cl_res = 0.5, nrGenes = p$'n.markers')
#'  }
#' }
#' @export
PlotTopGenesPerCluster <- function(obj = combined.obj, cl_res = res, nrGenes = p$'n.markers'
                                   , order.by = c("combined.score","avg_log2FC", "p_val_adj")[1]
                                   , df_markers = obj@misc$"df.markers"[[paste0("res.",cl_res)]]) {
  topX.markers <- GetTopMarkers(df = df_markers,  n= nrGenes
                                , order.by = order.by )
  ls.topMarkers <-  splitbyitsnames(topX.markers)
  for (i in 1:length(ls.topMarkers)) {
    multiFeaturePlot.A4(list.of.genes = ls.topMarkers[[i]], obj = obj, subdir = F
                        , prefix = ppp("DEG.markers.res",cl_res,"cluster",names(ls.topMarkers)[i]))
  }
}




# _________________________________________________________________________________________________

#' @title qFeatureScatter
#' @description Quickly plot and save a FeatureScatter plot.
#' @param feature1 PARAM_DESCRIPTION, Default: 'TOP2A'
#' @param feature2 PARAM_DESCRIPTION, Default: 'ID2'
#' @param obj Seurat object, Default: combined.obj
#' @param ext File extension for saving, Default: 'png'
#' @param logX logX
#' @param logY logY
#' @param plot PARAM_DESCRIPTION, Default: TRUE
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
qFeatureScatter <- function(feature1 = "TOP2A", feature2 = "ID2", obj = combined.obj
                            , ext ="png", plot = TRUE
                            , logX = F, logY = F
                            , ...) {
  plotname <- kpp(feature1,"VS", feature2)
  p <- FeatureScatter(object = obj, feature1 = feature1, feature2 = feature2, ...) +
    ggtitle(paste("Correlation", plotname)) +
    theme_linedraw()

  if (logX) p <- p + scale_x_log10()
  if (logY) p <- p + scale_y_log10()

  fname = kpp("FeatureScatter", plotname)
  ggExpress::qqSave(ggobj = p, title = plotname, ext = ext, w = 8, h = 5)
  if (plot) p
}

# _________________________________________________________________________________________________
#' @title qQC.plots.BrainOrg
#' @description Quickly plot key QC markers in brain organoids
#' @param QC.Features Any numeric metadata columns
#' @param obj Seurat object, Default: combined.obj
#'
#' @examples qQC.plots.BrainOrg.RV()
#' @export

qQC.plots.BrainOrg <- function(obj = combined.obj, title = "Top 4 QC markers on UMAP"
                               , nrow = 2, ncol = 2
                               , QC.Features = c('nFeature_RNA', 'percent.ribo', 'percent.mito', 'log10.HGA_Markers')) {
  print(QC.Features)
  n.found <- setdiff(QC.Features, colnames(obj@meta.data))
  stopif(length(n.found),message = paste("n.found:", n.found))
  px <- list(
    'A' = qUMAP(QC.Features[1], save.plot = F, obj = obj) + NoAxes(),
    'B' = qUMAP(QC.Features[2], save.plot = F, obj = obj) + NoAxes(),
    'C' = qUMAP(QC.Features[3], save.plot = F, obj = obj) + NoAxes(),
    'D' = qUMAP(QC.Features[4], save.plot = F, obj = obj) + NoAxes()
  )
  qA4_grid_plot(plot_list = px
                , plotname = title
                , w = hA4, h = wA4
                , nrow = nrow, ncol = ncol
  )
}

# _________________________________________________________________________________________________
#' @title qMarkerCheck.BrainOrg
#' @description Quickly plot key markers in brain organoids
#' @param obj Seurat object, Default: combined.obj
#' @param custom.genes PARAM_DESCRIPTION, Default: F
#' @param suffix Folder name suffix, Default: ""
#' @examples
#' \dontrun{
#' if(interactive()){
#'  qMarkerCheck.BrainOrg(combined.obj)
#'  }
#' }
#' @export
qMarkerCheck.BrainOrg <- function(obj = combined.obj, custom.genes = F, suffix = "") {
  Signature.Genes.Top16 <- if (!isFALSE(custom.genes)) custom.genes else
  {
    Signature.Genes.Top16  <- c(
      `dl-EN` = "KAZN", `ul-EN` = "SATB2" # dl-EN = deep layer excitatory neuron
      , `Immature neurons` = "SLA", Interneurons = "DLX6-AS1"
      , `Intermediate progenitor` = "EOMES",  `Intermediate progenitor1` = "TAC3"
      , `S-phase` = "TOP2A", `G2M-phase` = "HIST1H4C"
      , `oRG` = "ID4", `oRG` = "HOPX" # oRG outer radial glia
      , Astroglia = "GFAP", Astrocyte = "S100B"
      , `Hypoxia/Stress` = "DDIT4", Glycolytic = "PDK1"
      , `Low-Quality` = "POLR2A", `Choroid.Plexus` = "DCN"
      # , `Choroid.Plexus` = "OTX2", `Choroid.Plexus` = "BMP4"
    )
    print(Signature.Genes.Top16)
  }
  print(as_tibble_from_namedVec(Signature.Genes.Top16))
  multiFeaturePlot.A4(obj = obj, list.of.genes = Signature.Genes.Top16, layout = "tall"
                      , foldername = sppp('Signature.Genes.Top16', suffix))
}





# _________________________________________________________________________________________________
#' @title getDiscretePalette
#' @description Generate a Discrete color Palette.
#' @param ident.used PARAM_DESCRIPTION, Default: GetClusteringRuns()[1]
#' @param obj Seurat object, Default: combined.obj
#' @param palette.used PARAM_DESCRIPTION, Default: c("alphabet", "alphabet2", "glasbey", "polychrome", "stepped")[1]
#' @param show.colors PARAM_DESCRIPTION, Default: F
#' @examples
#' \dontrun{
#' if(interactive()){
#'  getDiscretePalette()
#'  }
#' }
#' @export
getDiscretePalette <- function(ident.used = GetClusteringRuns()[1]
                               , obj = combined.obj
                               , palette.used = c("alphabet", "alphabet2", "glasbey", "polychrome", "stepped")[1]
                               , show.colors = F) {
  n.clusters <-  nrow(unique(obj[[ident.used]]))
  colz <- DiscretePalette(n = n.clusters, palette = palette.used)
  if (anyNA(colz)) {
    colzOK <- na.omit.strip(colz)
    repNeeded <- ceiling(length(colz) / length(colzOK) )
    colzFixed <- rep(colzOK,  repNeeded)[1:length(colz)]
    stopif(anyNA(colzFixed))
    colz <- colzFixed
  }
  if (show.colors) Color_Check(colz)
  return(colz)
}


# _________________________________________________________________________________________________
# plotting.dim.reduction.3D.R
# ____________________________________________________________________ ----
# source('~/GitHub/Packages/Seurat.utils/Functions/plotting.dim.reduction.3D.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Plotting.dim.reduction.3D.R"))
# Source: self + https://github.com/Dragonmasterx87/Interactive-3D-Plotting-in-Seurat-3.0.0

# Requirements __________________________________________
# try(library(plotly), silent = T)
# try(library(MarkdownReports), silent = T)
# try(library(htmlwidgets), silent = T)

# May also require
# try (source('~/GitHub/Packages/CodeAndRoll/CodeAndRoll.R'),silent= T) # generic utilities funtions
# require('MarkdownReports') # require("devtools") # plotting related utilities functions # devtools::install_github(repo = "vertesy/MarkdownReports")


# _________________________________________________________________________________________________
#' @title ww.check.if.3D.reduction.exist
#' @description ww.check.if.3D.reduction.exist in backup slot #
#' @param obj Seurat object, Default: obj
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
ww.check.if.3D.reduction.exist <- function(obj = obj) { # ww.check.if.3D.reduction.exist in backup slot
  if( !("UMAP_3" %in% colnames(obj@reductions$'umap'))) {
    stopif2( is.null(obj@misc$reductions.backup$'umap3d')
             , "No 3D umap found in backup slot, @misc$reductions.backup. Run SetupReductionsNtoKdimensions() first.")
    RecallReduction(obj = obj, dim = 3, reduction = "umap")
  } else { # Reduction found in normal UMAP slot
    obj
  }
}

# _________________________________________________________________________________________________
#' @title ww.check.quantile.cutoff.and.clip.outliers
#' @description Helper function.
#' @param expr.vec PARAM_DESCRIPTION, Default: plotting.data[, gene]
#' @param quantileCutoffX PARAM_DESCRIPTION, Default: quantileCutoff
#' @param min.cells.expressing PARAM_DESCRIPTION, Default: 10
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
ww.check.quantile.cutoff.and.clip.outliers <- function(expr.vec = plotting.data[,gene], quantileCutoffX = quantileCutoff, min.cells.expressing = 10) {
  expr.vec.clipped <- CodeAndRoll2::clip.outliers.at.percentile(expr.vec, probs = c(1 - quantileCutoffX, quantileCutoffX))
  if( sum(expr.vec.clipped > 0) > min.cells.expressing ){
    expr.vec <- expr.vec.clipped
  } else {
    iprint("WARNING: quantile.cutoff too stringent, would leave <", min.cells.expressing, "cells. It is NOT applied.")
  }
  return(expr.vec)
}

# _________________________________________________________________________________________________
#' @title plot3D.umap.gene
#' @description Plot a 3D umap with gene expression. Uses plotly. Based on github.com/Dragonmasterx87. #
#' @param gene gene of interest, Default: 'TOP2A'
#' @param obj Seurat object, Default: combined.obj
#' @param quantileCutoff PARAM_DESCRIPTION, Default: 0.99
#' @param def.assay PARAM_DESCRIPTION, Default: c("integrated", "RNA")[2]
#' @param suffix A suffix added to the filename, Default: NULL
#' @param AutoAnnotBy PARAM_DESCRIPTION, Default: GetNamedClusteringRuns(obj)[1]
#' @param alpha PARAM_DESCRIPTION, Default: 0.5
#' @param dotsize PARAM_DESCRIPTION, Default: 1.25
#' @examples
#' \dontrun{
#' if(interactive()){
#'  plot3D.umap.gene(obj = combined.obj, gene = "DDIT4", quantileCutoff = .95)
#'  plot3D.umap.gene(obj = combined.obj, gene = "percent.mito", quantileCutoff = .95) # for continous meta variables
#'  plot3D.umap.gene(obj = combined.obj, gene = "nFeature_RNA", quantileCutoff = .95) # for continous meta variables
#'  }
#' }
#' @export
plot3D.umap.gene <- function(gene="TOP2A", obj = combined.obj # Plot a 3D umap with gene expression. Uses plotly. Based on github.com/Dragonmasterx87.
                             , quantileCutoff = .99, def.assay = c("integrated", "RNA")[2]
                             , suffix = NULL, AutoAnnotBy = GetNamedClusteringRuns(obj)[1]
                             , alpha = .5, dotsize = 1.25 ){
  # stopifnot(AutoAnnotBy %in% colnames(obj@meta.data) | AutoAnnotBy = FALSE)

  obj <- ww.check.if.3D.reduction.exist(obj = obj)
  stopifnot((gene %in% rownames(obj) | gene %in% colnames(obj@meta.data)))
  DefaultAssay(object = obj) <- def.assay; iprint(DefaultAssay(object = obj), "assay")

  plotting.data <- FetchData(object = obj, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "Expression" = gene), slot = 'data')

  plotting.data$'Expression' <- ww.check.quantile.cutoff.and.clip.outliers(expr.vec = plotting.data[,gene], quantileCutoffX = quantileCutoff, min.cells.expressing = 10)
  CodeAndRoll2::clip.outliers.at.percentile(plotting.data[,gene], probs = c(1 - quantileCutoff, quantileCutoff))
  plotting.data$'label' <- paste(rownames(plotting.data), " - ", plotting.data[,gene], sep = "")

  ls.ann.auto <- if (AutoAnnotBy != FALSE) {
    Annotate4Plotly3D(obj = obj, plotting.data. = plotting.data, AnnotCateg = AutoAnnotBy)
  } else { NULL }

  plt <- plotly::plot_ly(data = plotting.data
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
  ) %>% plotly::layout(title = gene, scene = list(annotations = ls.ann.auto))
  SavePlotlyAsHtml(plt, category. = gene, suffix. = suffix)
  return(plt)
}




# _________________________________________________________________________________________________
#' @title plot3D.umap
#' @description Plot a 3D umap based on one of the metadata columns. Uses plotly. Based on github.com/Dragonmasterx87. #
#' @param category PARAM_DESCRIPTION, Default: 'v.project'
#' @param obj Seurat object, Default: combined.obj
#' @param suffix A suffix added to the filename, Default: NULL
#' @param AutoAnnotBy PARAM_DESCRIPTION, Default: GetNamedClusteringRuns(obj)[1]
#' @param dotsize PARAM_DESCRIPTION, Default: 1.25
#' @examples
#' \dontrun{
#' if(interactive()){
#'  plot3D.umap(combined.obj, category = "Phase")
#'  }
#' }
#' @export
plot3D.umap <- function(category="v.project", obj = combined.obj # Plot a 3D umap based on one of the metadata columns. Uses plotly. Based on github.com/Dragonmasterx87.
                        , suffix = NULL, AutoAnnotBy = GetNamedClusteringRuns(obj)[1]
                        , dotsize = 1.25) {

  stopifnot(category %in% colnames(obj@meta.data))
  obj <- ww.check.if.3D.reduction.exist(obj = obj)

  plotting.data <- FetchData(object = obj, vars = c("UMAP_1", "UMAP_2", "UMAP_3", category))
  colnames(plotting.data)[4] = "category"
  plotting.data$label <- paste(rownames(plotting.data))   # Make a column of row name identities (these will be your cell/barcode names)

  ls.ann.auto <- if (AutoAnnotBy != FALSE) {
    Annotate4Plotly3D(obj = obj, plotting.data. = plotting.data, AnnotCateg = AutoAnnotBy)
  } else { NULL }

  plt <- plotly::plot_ly(data = plotting.data
                         , x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3
                         , type = "scatter3d"
                         , mode = "markers"
                         , marker = list(size = dotsize)
                         , text = ~label
                         , color = ~category
                         , colors = gg_color_hue(length(unique(plotting.data$'category')))
                         # , hoverinfo="text"
  ) %>% plotly::layout(title = category, scene = list(annotations = ls.ann.auto))
  SavePlotlyAsHtml(plt, category. = category, suffix. = suffix)
  return(plt)
}


# _________________________________________________________________________________________________
#' @title SavePlotlyAsHtml
#' @description Save Plotly 3D scatterplot as an html file. #
#' @param plotly_obj PARAM_DESCRIPTION
#' @param category. PARAM_DESCRIPTION, Default: category
#' @param suffix. PARAM_DESCRIPTION, Default: NULL
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[htmlwidgets]{saveWidget}}
#' @export
#' @importFrom htmlwidgets saveWidget
SavePlotlyAsHtml <- function(plotly_obj, category.=category, suffix. = NULL) { # Save Plotly 3D scatterplot as an html file.
  OutputDir <- if (exists("OutDir")) OutDir else getwd()
  name.trunk <- kpp("umap.3D", category., suffix., idate(), "html")
  fname <- kpps(OutputDir, name.trunk)
  iprint("Plot saved as:", fname)
  htmlwidgets::saveWidget(plotly_obj, file = fname, selfcontained = TRUE, title = category.)
}


# _________________________________________________________________________________________________
#' @title BackupReduction
#' @description Backup UMAP to `obj@misc$reductions.backup` from `obj@reductions$umap`. #
#' @param obj Seurat object, Default: combined.obj
#' @param dim Numer of dimensions used, Default: 2
#' @param reduction UMAP, tSNE, or PCA (Dim. reduction to use), Default: 'umap'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  obj <- BackupReduction(obj = obj, dim = 2, reduction = umap"")
#'  }
#' }
#' @export
BackupReduction <- function(obj = combined.obj, dim = 2, reduction="umap") { # Backup UMAP to `obj@misc$reductions.backup` from `obj@reductions$umap`.
  if (is.null(obj@misc$"reductions.backup")) obj@misc$"reductions.backup" <- list()
  dslot = paste0(reduction,dim,"d")
  obj@misc$reductions.backup[[dslot]] <- obj@reductions[[reduction]]
  return(obj)
}


# _________________________________________________________________________________________________
#' @title SetupReductionsNtoKdimensions
#' @description Calculate N-to-K dimensional umaps (default = 2:3); and back them up UMAP to `obj@misc$reductions.backup` from @reductions$umap #
#' @param obj Seurat object, Default: combined.obj
#' @param nPCs PARAM_DESCRIPTION, Default: p$n.PC
#' @param dimensions PARAM_DESCRIPTION, Default: 3:2
#' @param reduction UMAP, tSNE, or PCA (Dim. reduction to use), Default: 'umap'
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @examples
#' \dontrun{
#' if(interactive()){
#'  combined.obj <- SetupReductionsNtoKdimensions(obj = combined.obj, nPCs = p$'n.PC', dimensions = 2:3, reduction="umap"); qUMAP()
#'  }
#' }
#' @export
SetupReductionsNtoKdimensions <- function(obj = combined.obj, nPCs = p$'n.PC', dimensions = 3:2, reduction="umap", ...) { # Calculate N-to-K dimensional umaps (default = 2:3); and back them up UMAP to `obj@misc$reductions.backup` from @reductions$umap
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


# _________________________________________________________________________________________________
#' @title RecallReduction
#' @description Set active UMAP to `obj@reductions$umap` from `obj@misc$reductions.backup`. #
#' @param obj Seurat object, Default: combined.obj
#' @param dim Numer of dimensions used, Default: 2
#' @param reduction UMAP, tSNE, or PCA (Dim. reduction to use), Default: 'umap'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  combined.obj <- RecallReduction(obj = combined.obj, dim = 2, reduction="umap"); qUMAP()
#'  combined.obj <- RecallReduction(obj = combined.obj, dim = 3, reduction="umap"); qUMAP()
#'  }
#' }
#' @export
RecallReduction <- function(obj = combined.obj, dim = 2, reduction="umap") { # Set active UMAP to `obj@reductions$umap` from `obj@misc$reductions.backup`.
  dslot = paste0(reduction,dim,"d")
  reduction.backup <- obj@misc$reductions.backup[[dslot]]
  msg <-  paste(dim, "dimensional", reduction, "from obj@misc$reductions.backup" )
  stopif(is.null(reduction.backup), message = paste0(msg," is NOT FOUND")); iprint(msg, "is set active. " )
  stopifnot(dim == ncol(reduction.backup))
  obj@reductions[[reduction]] <- reduction.backup
  return(obj)
}



# _________________________________________________________________________________________________
#' @title Annotate4Plotly3D
#' @description Create annotation labels for 3D plots. Source https://plot.ly/r/text-and-annotations/#3d-annotations #
#' @param obj Seurat object, Default: combined.obj
#' @param plotting.data. PARAM_DESCRIPTION, Default: plotting.data
#' @param AnnotCateg PARAM_DESCRIPTION, Default: AutoAnnotBy
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
Annotate4Plotly3D <- function(obj = combined.obj # Create annotation labels for 3D plots. Source https://plot.ly/r/text-and-annotations/#3d-annotations
                              , plotting.data. = plotting.data
                              , AnnotCateg = AutoAnnotBy) {
  stopifnot(AnnotCateg %in% colnames(obj@meta.data))

  plotting.data.$'annot' <- FetchData(object = obj, vars = c(AnnotCateg))[,1]
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

# _________________________________________________________________________________________________
#' @title Plot3D.ListOfGenes
#' @description Plot and save list of 3D UMAP ot tSNE plots using plotly. #
#' @param obj Seurat object, Default: combined.obj
#' @param annotate.by PARAM_DESCRIPTION, Default: 'integrated_snn_res.0.7'
#' @param opacity PARAM_DESCRIPTION, Default: 0.5
#' @param cex Point size, Default: 1.25
#' @param default.assay PARAM_DESCRIPTION, Default: c("integrated", "RNA")[2]
#' @param ListOfGenes PARAM_DESCRIPTION, Default: c("BCL11B", "FEZF2", "EOMES", "DLX6-AS1", "HOPX", "DDIT4")
#' @param SubFolderName PARAM_DESCRIPTION, Default: ppp("plot3D", substitute(ListOfGenes))
#' @examples
#' \dontrun{
#' if(interactive()){
#'  CellTypeMarkers <- c(  "PGK1", "CTIP2" = "BCL11B" , "FEZF2", "EOMES", "DLX6-AS1", "HOPX", "DDIT4","TOP2A", "PTGDS", "EDNRB", "EGFR", "SCGN", "NR2F2", "EMX2", "GAD2", "DLX2", "SATB2")
#'  Plot3D.ListOfGenes(obj = combined.obj, ListOfGenes = CellTypeMarkers)
#'  }
#' }
#' @export
Plot3D.ListOfGenes <- function(obj = combined.obj # Plot and save list of 3D UMAP ot tSNE plots using plotly.
                               , annotate.by = "integrated_snn_res.0.7", opacity = 0.5, cex = 1.25, default.assay = c("integrated", "RNA")[2]
                               , ListOfGenes = c("BCL11B" , "FEZF2", "EOMES", "DLX6-AS1", "HOPX", "DDIT4")
                               , SubFolderName = ppp("plot3D", substitute(ListOfGenes))) {


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


# _________________________________________________________________________________________________
# _________________________________________________________________________________________________
#' @title Plot3D.ListOfCategories
#' @description Plot and save list of 3D UMAP ot tSNE plots using plotly. #
#' @param obj Seurat object, Default: combined.obj
#' @param annotate.by PARAM_DESCRIPTION, Default: 'integrated_snn_res.0.7'
#' @param cex Point size, Default: 1.25
#' @param default.assay PARAM_DESCRIPTION, Default: c("integrated", "RNA")[2]
#' @param ListOfCategories PARAM_DESCRIPTION, Default: c("v.project", "experiment", "Phase", "integrated_snn_res.0.7")
#' @param SubFolderName PARAM_DESCRIPTION, Default: ppp("plot3D", substitute(ListOfCategories))
#' @examples
#' \dontrun{
#' if(interactive()){
#'  categ3Dplots <- c("v.project","experiment", "Phase", "integrated_snn_res.0.7", "Area", "Individual", "Type")
#'  Plot3D.ListOfCategories(obj = combined.obj, ListOfCategories = categ3Dplots)
#'  }
#' }
#' @export
Plot3D.ListOfCategories <- function(obj = combined.obj # Plot and save list of 3D UMAP ot tSNE plots using plotly.
                                    , annotate.by = "integrated_snn_res.0.7", cex = 1.25, default.assay = c("integrated", "RNA")[2]
                                    , ListOfCategories = c("v.project","experiment", "Phase", "integrated_snn_res.0.7")
                                    , SubFolderName = ppp("plot3D", substitute(ListOfCategories))) {

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


# _________________________________________________________________________________________________
# _________________________________________________________________________________________________


# _________________________________________________________________________________________________
# plotting.filtering.R
# ____________________________________________________________________ ----
# source('~/GitHub/Packages/Seurat.utils/Functions/plotting.filtering.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Plotting.filtering.R"))


# _________________________________________________________________________________________________
#' @title PlotFilters
#' @description Plot filtering threshold and distributions, using four panels to highlight the relation between Gene- and UMI-count, ribosomal- and mitochondrial-content. #
#' @param ls.obj List of Seurat objects, Default: ls.Seurat
#' @param parentdir PARAM_DESCRIPTION, Default: OutDirOrig
#' @param suffices PARAM_DESCRIPTION, Default: names(ls.obj)
#' @param filetype PARAM_DESCRIPTION, Default: '.png'
#' @param below.mito PARAM_DESCRIPTION, Default: p$thr.lp.mito
#' @param above.mito PARAM_DESCRIPTION, Default: p$thr.hp.mito
#' @param below.ribo PARAM_DESCRIPTION, Default: p$thr.lp.ribo
#' @param above.ribo PARAM_DESCRIPTION, Default: p$thr.hp.ribo
#' @param below.nFeature_RNA PARAM_DESCRIPTION, Default: p$thr.lp.nFeature_RNA
#' @param above.nFeature_RNA PARAM_DESCRIPTION, Default: p$thr.hp.nFeature_RNA
#' @param subdir PARAM_DESCRIPTION, Default: kpp("Filtering.plots", "mito", p$thr.hp.mito, p$thr.lp.mito,
#'    "ribo", p$thr.hp.ribo, p$thr.lp.ribo, "nFeature", p$thr.hp.nFeature_RNA,
#'    p$thr.lp.nFeature_RNA, "/")
#' @param transparency PARAM_DESCRIPTION, Default: 0.25
#' @param cex Point size, Default: 0.75
#' @param theme.used PARAM_DESCRIPTION, Default: theme_bw(base_size = 18)
#' @param LabelDistFromTop PARAM_DESCRIPTION, Default: 200
#' @examples
#' \dontrun{
#' if(interactive()){
#'  PlotFilters(ls.Seurat)
#'  }
#' }
#' @seealso
#'  \code{\link[ggplot2]{ggplot}}, \code{\link[ggplot2]{labs}}, \code{\link[ggplot2]{geom_point}}
#' @export
#' @importFrom ggplot2 ggplot ggtitle geom_point
#' @importFrom Stringendo percentage_formatter
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
    , "] and the fraction of mitochondrial [", Stringendo::percentage_formatter(above.mito), ";" ,Stringendo::percentage_formatter(below.mito)
    , "] and ribosomal [",Stringendo::percentage_formatter(above.ribo), ";" ,Stringendo::percentage_formatter(below.ribo), "] reads."
  )


  theme_set(theme.used)
  create_set_OutDir(parentdir, subdir)
  # require(ggplot2)
  if (suffices == length(ls.obj)) print("ls.obj elements have no names (required).")

  for (i in 1:length(ls.obj)) {
    print(suffices[i])
    mm =  ls.obj[[i]]@meta.data

    AllMetaColumnsPresent <- all(c('nFeature_RNA', 'percent.mito', 'percent.ribo') %in% colnames(mm))
    if (!AllMetaColumnsPresent) {
      print(c('nFeature_RNA', 'percent.mito', 'percent.ribo'))
      print(c('nFeature_RNA', 'percent.mito', 'percent.ribo') %in% colnames(mm))
      print("Try to run:")
      print('objX <- add.meta.fraction(obj = objX, col.name = "percent.mito", gene.symbol.pattern =  "^MT\\.|^MT-")')
      print('objX <- add.meta.fraction(obj = objX, col.name = "percent.ribo", gene.symbol.pattern =  "^RPL|^RPS")')
      stop()
    }



    filt.nFeature_RNA = (mm$'nFeature_RNA' < below.nFeature_RNA & mm$'nFeature_RNA' > above.nFeature_RNA)
    filt.below.mito = (mm$'percent.mito' < below.mito & mm$'percent.mito' > above.mito)

    # filt.below.mito = (mm$'percent.mito' < below.mito)
    filt.below.ribo = (mm$'percent.ribo' < below.ribo & mm$'percent.ribo' > above.ribo)

    mm =  cbind(mm, filt.nFeature_RNA, filt.below.mito, filt.below.ribo)

    mm$colour.thr.nFeature <- cut(mm$'nFeature_RNA',
                                  breaks = c(-Inf, above.nFeature_RNA, below.nFeature_RNA, Inf),
                                  labels = c(paste0("LQ (<", above.nFeature_RNA,")"),
                                             paste0("HQ (", above.nFeature_RNA,"< X <", below.nFeature_RNA,")"),
                                             paste0("Dbl/Outlier (>", below.nFeature_RNA,")")
                                  )
    )

    A = ggplot(data = mm, aes(x = nFeature_RNA, fill = colour.thr.nFeature)) +
      geom_histogram(binwidth = 100) +
      ggtitle(paste("Cells between", above.nFeature_RNA,"and",below.nFeature_RNA, " UMIs are selected (", pc_TRUE(filt.nFeature_RNA), ")")) +
      geom_vline(xintercept = below.nFeature_RNA) +
      geom_vline(xintercept = above.nFeature_RNA);
    # A

    B = ggplot2::ggplot(mm, aes(x = nFeature_RNA, y = percent.mito)) +
      ggplot2::ggtitle(paste("Cells below", Stringendo::percentage_formatter(below.mito),
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
      ggtitle(paste("Cells below", Stringendo::percentage_formatter(below.ribo),
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

  # _________________________________________________________________________________________________
  create_set_Original_OutDir()
}



# _________________________________________________________________________________________________
# plotting.statistics.and.QC.R
# ____________________________________________________________________ ----
# source('~/GitHub/Packages/Seurat.utils/Functions/Plotting.statistics.and.QC.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Plotting.statistics.and.QC.R"))

# Source: self + web

# Requirements __________________________________________
# require(Seurat)
# require(ggplot2)
# tools for tools::toTitleCase

# May also require
# try (source('/GitHub/Packages/CodeAndRoll/CodeAndRoll.R'),silent= F) # generic utilities funtions
# require('MarkdownReports') # require("devtools") # plotting related utilities functions # devtools::install_github(repo = "vertesy/MarkdownReports")

# PCA percent of variation associated with each PC ------------------------------------------------------------
#' @title seu.PC.var.explained
#' @description Determine percent of variation associated with each PC. For normal prcomp objects, see: PCA.percent.var.explained().
#' @param obj Seurat object, Default: combined.obj
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
seu.PC.var.explained <- function(obj =  combined.obj) { # Determine percent of variation associated with each PC.
  pct <- obj@reductions$pca@stdev / sum(obj@reductions$pca@stdev) * 100
  names(pct) =1:length(obj@reductions$pca@stdev)
  return(pct)
}

# plot percent of variation associated with each PC ------------------------------------------------------------
#' @title seu.plot.PC.var.explained
#' @description Plot the percent of variation associated with each PC. #
#' @param obj Seurat object, Default: combined.obj
#' @param use.MDrep PARAM_DESCRIPTION, Default: F
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
seu.plot.PC.var.explained <- function(obj =  combined.obj, use.MDrep = F) { # Plot the percent of variation associated with each PC.
  pct <- seu.PC.var.explained(obj)
  if (use.MDrep) {
    wbarplot(pct , xlab = "Principal Components", ylab = "% of variation explained")
    barplot_label(round(pct, digits = 2), barplotted_variable = pct, cex = .5 )
  } else {
    ggExpress::qbarplot(vec = pct, xlab = "Principal Components", ylab =  "% of variation explained", w = 10, h = 5, hline = 1 )
  }
}


# scBarplot.CellFractions ------------------------------------------------------------
#' @title scBarplot.CellFractions
#' @description Barplot the Fraction of cells per cluster.
#' @param obj Seurat object, Default: combined.obj
#' @param group.by PARAM_DESCRIPTION, Default: 'integrated_snn_res.0.5.ordered'
#' @param fill.by PARAM_DESCRIPTION, Default: 'age'
#' @param downsample PARAM_DESCRIPTION, Default: T
#' @param plotname Title of the plot, Default: paste(tools::toTitleCase(fill.by), "proportions")
#' @param hlines PARAM_DESCRIPTION, Default: c(0.25, 0.5, 0.75)
#' @param seedNr PARAM_DESCRIPTION, Default: 1989
#' @param return_table return contingency table instead of a barplot, Default: F
#' @examples
#' \dontrun{
#' if(interactive()){
#'  scBarplot.CellFractions(obj = combined.obj, group.by = "integrated_snn_res.0.1", fill.by = "Phase", downsample = T)
#'  scBarplot.CellFractions(obj = combined.obj, group.by = "integrated_snn_res.0.1", fill.by = "Phase", downsample = F)
#'  }
#' }
#' @seealso
#'  \code{\link[tools]{toTitleCase}}
#' @export
#' @importFrom tools toTitleCase
scBarplot.CellFractions <- function(obj = combined.obj
                                    , group.by = "integrated_snn_res.0.5.ordered", fill.by = "age", downsample = T
                                    , plotname = paste(tools::toTitleCase(fill.by), "proportions"), hlines = c(.25, .5, .75)
                                    , return_table = F
                                    , seedNr = 1989) {
  set.seed(seedNr)
  pname.suffix <- capt.suffix <- NULL
  if (downsample) {
    downsample <- min (table(obj@meta.data[[fill.by]]))
    pname.suffix <- "(downsampled)"
    capt.suffix <- paste("Downsampled to", downsample, "cells in the smallest", fill.by, "group.")
  }
  caption_ <- paste("Numbers denote # cells.", capt.suffix)
  pname_ <- paste(plotname, pname.suffix)

  contingency.table <- table(obj@meta.data[ ,group.by], obj@meta.data[, fill.by])
  print(contingency.table)
  if (return_table) {
    list(
      'values' = contingency.table,
      'percentages' = CodeAndRoll2::colDivide(mat = contingency.table, vec = colSums(contingency.table))
    )
  } else {
    obj@meta.data %>%
      # group_by( (!!as.name(fill.by)) ) %>%
      { if (downsample) sample_n(., downsample) else . } %>%
      group_by( (!!as.name(group.by)) ) %>%

      ggplot( aes(fill = (!!(as.name(fill.by))),  x = (!!(as.name(group.by)))) ) +
      geom_hline( yintercept = hlines, lwd = 1.5)  +
      geom_bar( position = "fill" ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      geom_text(aes(label = ..count..), stat='count',position = position_fill(vjust = 0.5)) +
      labs(title = pname_,  x = "Clusters", y = "Fraction", caption = caption_) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  } # else barplot
}




# _________________________________________________________________________________________________
#' @title scBarplot.CellsPerCluster
#' @description Barplot the Fraction of cells per cluster. (dupl?)
#'
#' @param obj Seurat object, Default: combined.obj
#' @param ident identity used, Default: 'cl.names.KnownMarkers.0.5'
#' @param label True: displays cell count, but you can provide anything in a vector.
#' @param palette Color palette. Default: glasbey.
#' @param return_table Should it return the plotting data instead of the plot?
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @param sort Sort by cluster size? Default: F
#' @param suffix File name suffix
#'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  scBarplot.CellsPerCluster(); scBarplot.CellsPerCluster(sort = T)
#'  }
#' }
#' @export scBarplot.CellsPerCluster

scBarplot.CellsPerCluster <- function(ident =  GetOrderedClusteringRuns()[1]
                                      , sort = F
                                      , label = list(T, 'percent')[[1]]
                                      , suffix = if (label == 'percent') 'percent' else NULL
                                      , palette = c("alphabet", "alphabet2", "glasbey", "polychrome", "stepped")[3]
                                      , obj = combined.obj
                                      , return_table = F
                                      , ...) {
  cell.per.cl <- obj[[ident]][,1]
  cell.per.cluster <- (table(cell.per.cl))
  if (sort) cell.per.cluster <- sort(cell.per.cluster)
  lbl <- if (isFALSE(label)) { NULL
    } else if (label == 'percent') { percentage_formatter(cell.per.cluster/sum(cell.per.cluster))
      } else if (label == 'T') { cell.per.cluster
        } else label

  n.clusters <- length(cell.per.cluster)
  if (return_table) {
    cell.per.cluster
  } else {
    ggExpress::qbarplot(cell.per.cluster, subtitle = ident, suffix = kpp(ident, suffix)
             , col = 1:n.clusters
             , xlab.angle = 45
             , label = lbl
             # , col = getClusterColors(ident = ident, show = T)
             , palette_use = DiscretePalette(n = n.clusters, palette = palette)
             , ...)
  }

}



# scBarplot.CellsPerObject ------------------------------------------------------------

#' @title scBarplot.CellsPerObject
#' @description Take a List of Seurat objects and draw a barplot for the number of cells per object. #
#' @param ls.Seu PARAM_DESCRIPTION, Default: ls.Seurat
#' @param plotname Title of the plot, Default: 'Nr.Cells.After.Filtering'
#' @param names PARAM_DESCRIPTION, Default: F
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
scBarplot.CellsPerObject <- function(ls.Seu = ls.Seurat
                                     , plotname="Nr.Cells.After.Filtering", xlab.angle = 45
                                     , names = F, ...) {
  cellCounts = unlapply(ls.Seu, ncol)
  names(cellCounts) = if (length(names) == length(ls.Seurat)) names else names(ls.Seurat)
  qbarplot(cellCounts, plotname = plotname
           , subtitle = paste(sum(cellCounts), "cells in total")
           , label = cellCounts
           , xlab.angle =  xlab.angle, ylab="Cells"
           , ...)

}

# _________________________________________________________________________________________________
#' @title BulkGEScatterPlot
#' @description Plot bulk scatterplots to identify differential expressed genes across conditions #
#' @param obj Seurat object, Default: combined.obj
#' @param clusters PARAM_DESCRIPTION, Default: 'cl.names.KnownMarkers.0.2'
#' @param TwoCategIdent PARAM_DESCRIPTION, Default: 'age'
#' @param genes.from.bulk.DE PARAM_DESCRIPTION, Default: rownames(df.markers.per.AGE)
#' @examples
#' \dontrun{
#' if(interactive()){
#'  BulkGEScatterPlot(obj = combined.obj, clusters = "cl.names.KnownMarkers.0.2", TwoCategIdent = 'age', genes.from.bulk.DE = rownames(df.markers.per.AGE))
#'  }
#' }
#' @export
BulkGEScatterPlot <- function(obj = combined.obj # Plot bulk scatterplots to identify differential expressed genes across conditions
                              , clusters = "cl.names.KnownMarkers.0.2", TwoCategIdent = 'age', genes.from.bulk.DE = rownames(df.markers.per.AGE)) {

  (SplitIdents <- unique(obj[[TwoCategIdent]][,1]))
  stopifnot(length(SplitIdents) == 2)

  Idents(obj) <- clusters
  IdentsUsed <- sort.natural(as.character(unique(Idents(obj))))
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
    p.clAv.AutoLabel[[i]] <- LabelPoints(plot = p.clAv[[i]], points = genes.to.label[[i]], xnudge = 0, ynudge = 0, repel = TRUE, size = 2);
    p.clAv.AutoLabel[[i]]

    "Pre-identified genes"
    p.clAv[[i]] <- LabelPoints(plot = p.clAv[[i]], points = genes.from.bulk.DE, repel = TRUE, size = 2);
  }

  PlotIter <- iterBy.over(1:NrPlots, by = 4)
  for (i in 1:length(PlotIter)) {
    plotLS = p.clAv.AutoLabel[PlotIter[[i]]]
    qqSaveGridA4(plotlist = plotLS, plots = 1:4, fname = ppp("BulkGEScatterPlot.AutoGenes",kpp(PlotIter[[i]]), "png"))

    plotLS = p.clAv[PlotIter[[i]]]
    qqSaveGridA4(plotlist = plotLS, plots = 1:4, fname= ppp("BulkGEScatterPlot.BulkGenes",kpp(PlotIter[[i]]), "png"))
  }
}







# _________________________________________________________________________________________________

#' @title sparse.cor
#' @description Sparse, fast correlation.
#' @param smat PARAM_DESCRIPTION
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#'
#'
sparse.cor <- function(smat){
  n <- nrow(smat)
  cMeans <- colMeans(smat)
  covmat <- (as.matrix(crossprod(smat)) - n * tcrossprod(cMeans))/(n - 1)
  sdvec <- sqrt(diag(covmat))
  cormat <- covmat / tcrossprod(sdvec)
  list(cov = covmat, cor = cormat)
}


sparse.cor4 <- function(x){
  n <- nrow(x)
  cMeans <- colMeans(x)
  covmat <- (as.matrix(crossprod(x)) - n*tcrossprod(cMeans))/(n-1)
  sdvec <- sqrt(diag(covmat))
  cormat <- covmat/tcrossprod(sdvec)
  list(cov=covmat,cor=cormat)
}

# Calc.Cor.Seurat ------------------------------------------------
#' @title Calc.Cor.Seurat
#' @description Calculate gene correlation on a Seurat object.
#' @param assay.use PARAM_DESCRIPTION, Default: 'RNA'
#' @param slot.use PARAM_DESCRIPTION, Default: 'data'
#' @param quantileX Quantile level, Default: 0.95
#' @param max.cells PARAM_DESCRIPTION, Default: 10000
#' @param seed random seed used, Default: p$seed
#' @param seed random seed used, Default: p$seed
#' @param digits PARAM_DESCRIPTION, Default: 2
#' @param obj Seurat object, Default: combined.obj
#' @examples
#' \dontrun{
#' if(interactive()){
#'  combined.obj <- calc.q99.Expression.and.set.all.genes(combined.obj, quantileX = 0.99, max.cells =  400000, set.all.genes = F)
#'  combined.obj <- Calc.Cor.Seurat(assay.use = "RNA", slot.use = "data", digits = 2, obj = combined.obj, quantile = 0.99, max.cells = 40000)
#'  }
#' }
#' @export
#' @importFrom tictoc tic toc
Calc.Cor.Seurat <- function(assay.use = "RNA", slot.use = "data"
                            , quantileX = 0.95, max.cells =  40000, seed = p$"seed"
                            , digits = 2, obj = combined.obj) {
  expr.mat <- GetAssayData(slot = slot.use, assay = assay.use, object = obj)
  if (ncol(expr.mat) > max.cells) {
    set.seed(seed = seed)
    cells.use <- sample(x = colnames(expr.mat), size = max.cells)
  } else {
    cells.use <- ncol(expr.mat)
  }

  qname = paste0("q", quantileX * 100)
  quantile_name = kpp("expr", qname)

  if (is.null(obj@misc[[quantile_name]])) iprint("Call: combined.obj <- calc.q99.Expression.and.set.all.genes(combined.obj, quantileX =",quantileX," first )")
  genes.HE = which_names(obj@misc[[quantile_name]] > 0)
  iprint("Pearson correlation is calculated for", length(genes.HE), "HE genes with expr.",qname,": > 0.")
  tictoc::tic(); ls.cor <- sparse.cor(smat = t(expr.mat[genes.HE, cells.use])); toc()
  ls.cor <- lapply(ls.cor, round, digits = 2)

  slot__name <- kpp(slot.use, assay.use, quantile_name)
  obj@misc[[kpp('cor', slot__name)]] <- ls.cor$'cor'
  obj@misc[[kpp('cov', slot__name)]] <- ls.cor$'cov'
  iprint("Stored under obj@misc$", kpp('cor', slot.use, assay.use), "or cov... ." )
  return(obj)
}


# _________________________________________________________________________________________________
#' @title plot.Metadata.Cor.Heatmap
#' @description Plot a heatmap with Metadata correlation values.
#' @param columns PARAM_DESCRIPTION, Default: c("nCount_RNA", "nFeature_RNA", "percent.mito", "percent.ribo")
#' @param cormethod PARAM_DESCRIPTION, Default: c("pearson", "spearman")[1]
#' @param main PARAM_DESCRIPTION, Default: paste("Metadata", cormethod, "correlations")
#' @param obj Seurat object, Default: combined.obj
#' @param w width of the plot, Default: 10
#' @param h height of the plot, Default: w
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[ggcorrplot]{ggcorrplot}}
#' @export plot.Metadata.Cor.Heatmap
#' @importFrom ggcorrplot ggcorrplot
plot.Metadata.Cor.Heatmap <- function(
  columns = c( "nCount_RNA", "nFeature_RNA", "percent.mito", "percent.ribo")
  , cormethod = c('pearson', 'spearman')[1]
  , main =paste( "Metadata", cormethod,"correlations")
  , obj = combined.obj
  , w = 10, h = w
  , ...){

  meta.data <- obj@meta.data
  columns.found <- intersect(colnames(meta.data), columns)

  corX <- cor(meta.data[ , columns.found], method = cormethod)
  pl <- ggcorrplot::ggcorrplot(corX, hc.order = TRUE, title = main
                               , type = "full", lab = T)
  ggExpress::qqSave(pl, fname = ppp(make.names(main),'pdf'), w = w, h = h)
  pl
}




# _________________________________________________________________________________________________
#' @title plot.Metadata.median.fraction.barplot
#' @description Barplot Metadata median values
#' @param columns PARAM_DESCRIPTION, Default: c("percent.mito", "percent.ribo")
#' @param suffix A suffix added to the filename, Default: NULL
#' @param group.by PARAM_DESCRIPTION, Default: GetClusteringRuns(obj = obj)[2]
#' @param method PARAM_DESCRIPTION, Default: c("median", "mean")[1]
#' @param min.thr PARAM_DESCRIPTION, Default: 2.5
#' @param return.matrix PARAM_DESCRIPTION, Default: F
#' @param main PARAM_DESCRIPTION, Default: paste(method, "read fractions per transcript class and cluster",
#'    suffix)
#' @param ylab PARAM_DESCRIPTION, Default: 'Fraction of transcriptome (%)'
#' @param percentify PARAM_DESCRIPTION, Default: T
#' @param subt PARAM_DESCRIPTION, Default: NULL
#' @param position PARAM_DESCRIPTION, Default: position_stack()
#' @param w width of the plot, Default: 10
#' @param h height of the plot, Default: 6
#' @param obj Seurat object, Default: combined.obj
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[dplyr]{summarise_all}}
#'  \code{\link[reshape2]{melt}}
#' @export plot.Metadata.median.fraction.barplot
#' @importFrom dplyr summarize_all
#' @importFrom reshape2 melt
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
  ggExpress::qqSave(pl, fname = ppp(make.names(main),'pdf'), w = w, h = h)
  pl
  if (return.matrix) mat.cluster.medians1 else pl
}

# plot.Metadata.median.fraction.barplot()



# _________________________________________________________________________________________________
#' @title plot.Gene.Cor.Heatmap
#' @description Plot a gene correlation heatmap.
#' @param genes Genes of iinterest, Default: WU.2017.139.IEGsf
#' @param assay.use PARAM_DESCRIPTION, Default: 'RNA'
#' @param slot.use PARAM_DESCRIPTION, Default: c("data", "scale.data", "data.imputed")[1]
#' @param quantileX Quantile level, Default: 0.95
#' @param min.g.cor PARAM_DESCRIPTION, Default: 0.3
#' @param calc.COR PARAM_DESCRIPTION, Default: FALSE
#' @param cutRows PARAM_DESCRIPTION, Default: NULL
#' @param cutCols PARAM_DESCRIPTION, Default: cutRows
#' @param obj Seurat object, Default: combined.obj
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export plot.Gene.Cor.Heatmap
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

  qname = paste0("expr.q", quantileX * 100)
  slotname_cor.mat <- kpp('cor', slot.use, assay.use, qname)
  cor.mat <- obj@misc[[slotname_cor.mat]]

  if (is.null(cor.mat)) {
    iprint(slotname_cor.mat, " not found in @misc.")
    iprint("Correlation slots present in @misc:", CodeAndRoll2::grepv(names(obj@misc), pattern = "^cor"))

    # Calculate ------------------------------------
    if (calc.COR) {
      print("Calculating correlation now.")
      genes.found <- check.genes(list.of.genes = genes)
      iprint(length(genes.found), "genes are found in the object.")
      if (length(genes.found) > 200) iprint("Too many genes found in data, cor will be slow: ", length(genes.found))
      ls.cor <- sparse.cor(t(expr.mat[genes.found,]))
      cor.mat <- ls.cor$cor
    } else { stop() }
  } else {
    print("Correlation is pre-calculated")
    genes.found <- intersect(genes, rownames(cor.mat))
    iprint(length(genes.found), "genes are found in the correlation matrix.")
    cor.mat <- cor.mat[genes.found, genes.found]
  }


  # Filter ------------------------------------
  diag(cor.mat) <- NaN
  corgene.names <- union(
    which_names(rowMax(cor.mat) >= min.g.cor),
    which_names(rowMin(cor.mat) <= -min.g.cor)
  )
  iprint(length(corgene.names), "genes are more (anti-)correlated than +/-:", min.g.cor)

  pname = paste0("Pearson correlations of ", substitute(genes),"\n min.cor:", min.g.cor, " | ",  assay.use ,'.', slot.use )
  o.heatmap <- pheatmap(cor.mat[corgene.names,corgene.names],main = pname, cutree_rows = cutRows, cutree_cols = cutCols, ...)
  wplot_save_pheatmap(o.heatmap, filename = make.names(pname))

  # return values
  maxCorrz <- rowMax(cor.mat)[corgene.names]; names(maxCorrz) <- corgene.names
  dput(maxCorrz)
}



# plot.clust.size.distr ------------------------------------------------
#' @title plot.clust.size.distr
#' @description Barplot of Histogram of cluster size distribution
#' @param obj Seurat object, Default: combined.obj
#' @param ident identity used, Default: GetClusteringRuns()[2]
#' @param plot PARAM_DESCRIPTION, Default: T
#' @param thr.hist PARAM_DESCRIPTION, Default: 30
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @examples
#' \dontrun{
#' if(interactive()){
#'  plot.clust.size.distr()
#'  }
#' }
#' @export plot.clust.size.distr
#' @importFrom Stringendo percentage_formatter
plot.clust.size.distr <- function(obj = combined.obj, ident = GetClusteringRuns()[2]
                                  , plot = T, thr.hist = 30, ...) {
  clust.size.distr <- table(obj@meta.data[,ident])
  print(clust.size.distr)
  resX <- gsub(pattern = ".*res\\.", replacement = '',x = ident)
  ptitle <- ppp('clust.size.distr', ident)
  psubtitle <- paste("Nr.clusters:", length(clust.size.distr)
                     , "| median:", median(clust.size.distr)
                     , "| CV:", Stringendo::percentage_formatter(cv(clust.size.distr))
  )
  xlb = "Cluster size (cells)"
  xlim = c(0, max(clust.size.distr))

  if (plot) {
    if (length(clust.size.distr) < thr.hist) {
      ggExpress::qbarplot(clust.size.distr, plotname = ptitle, subtitle = psubtitle, xlab = xlb, ...)
    } else {
      ggExpress::qhistogram(vec = clust.size.distr, plotname = ptitle, subtitle = psubtitle, xlab = xlb, xlim = xlim, ...)
    }
  } else {    "return vector"
    clust.size.distr
  }

}




# _________________________________________________________________________________________________
#' @title gene.expression.level.plots
#' @description Histogram of gene expression levels.
#' @param gene gene of interest, Default: 'TOP2A'
#' @param obj Seurat object, Default: ls.Seurat[[1]]
#' @param slot slot in the Seurat object. Default: c("counts", "data")[2]
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
gene.expression.level.plots <- function(gene = 'TOP2A', obj = ls.Seurat[[1]], slot = c('counts', 'data')[2]
                                        , ... ) {
  slot = 'data'
  print(gene)
  if (gene %in% rownames(obj)) {
    GEX.Counts <- GetAssayData(object = obj, assay = 'RNA', slot = slot)

    GEX.Counts.total <- rowSums(GEX.Counts)
    genes.expression <- GEX.Counts.total[gene]
    mean.expr <- iround(mean(GEX.Counts[gene,]))

    suffx = if (slot == 'counts') 'raw' else 'normalised, logtransformed'
    (pname = paste(gene, 'and the', suffx, 'transcript count distribution'))

    ggExpress::qhistogram(GEX.Counts.total, vline = genes.expression, logX = T, w = 7, h = 4
               , subtitle = paste('It belong to the top', pc_TRUE(GEX.Counts.total > genes.expression), 'of genes (black line). Mean expr:', mean.expr)
               , plotname = pname, xlab = 'Total Transcripts in Dataset', ylab = 'Number of Genes'
               , ...)
  } else { print("     !!! Gene not found in object!")}
}

# _________________________________________________________________________________________________
#' @title PrctCellExpringGene
#' @description From Github/Ryan-Zhu https://github.com/satijalab/seurat/issues/371 #
#' @param genes Genes of iinterest
#' @param group.by PARAM_DESCRIPTION, Default: 'all'
#' @param obj Seurat object, Default: combined.obj
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
PrctCellExpringGene <- function(genes, group.by = "all", obj = combined.obj){ # From Github/Ryan-Zhu https://github.com/satijalab/seurat/issues/371
  if(group.by == "all"){
    prct = unlist(lapply(genes, ww.calc_helper, object = obj))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }

  else{
    list = SplitObject(obj, group.by)
    factors = names(list)
    results = lapply(list, PrctCellExpringGene, genes = genes)
    for (i in 1:length(factors)) {
      results[[i]]$Feature = factors[i]
    }
    combined = do.call("rbind", results)
    return(combined)
  }
}


# _________________________________________________________________________________________________
#' @title ww.calc_helper
#' @description From Github/Ryan-Zhu https://github.com/satijalab/seurat/issues/371 #
#' @param obj Seurat object
#' @param genes Genes of iinterest
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
ww.calc_helper <- function(obj, genes){ # From Github/Ryan-Zhu https://github.com/satijalab/seurat/issues/371
  counts = obj[['RNA']]@counts
  ncells = ncol(counts)
  if (genes %in% row.names(counts)) {
    sum(counts[genes, ] > 0) / ncells
  } else{
    return(NA)
  }
}

# _________________________________________________________________________________________________

#' @title scBarplot.FractionAboveThr
#' @description Barplot the fraction of cell above a threshold value (based on a meta.data column), in each cluster.
#' @param thrX PARAM_DESCRIPTION, Default: 0
#' @param value.col PARAM_DESCRIPTION, Default: 'percent.ribo'
#' @param id.col PARAM_DESCRIPTION, Default: 'cl.names.top.gene.res.0.3'
#' @param obj Seurat object, Default: combined.obj
#' @param suffix Suffix for filename
#' @param return.df PARAM_DESCRIPTION, Default: F
#' @param label label barplot
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @examples
#' \dontrun{
#' if(interactive()){
#'  scBarplot.FractionAboveThr(id.col =  'cl.names.top.gene.res.0.3', value.col = 'percent.ribo', thrX = 0)
#'  }
#' }
#' @seealso
#'  \code{\link[dplyr]{select}}, \code{\link[dplyr]{se-deprecated}}
#' @export
#' @importFrom dplyr select group_by_
scBarplot.FractionAboveThr <- function(thrX = 0.3, suffix= NULL, value.col = 'percent.ribo', id.col =  'cl.names.top.gene.res.0.3'
                                      , obj = combined.obj, return.df = F, label = F
                                      , ...) { # Calculat the fraction of cells per cluster above a certain threhold
  meta = obj@meta.data
  (df_cells_above <- meta %>%
      dplyr::select(c(id.col, value.col))  %>%
      dplyr::group_by_(id.col)  %>%
      summarize(n_cells = n(),
                n_cells_above = sum(!!as.name(value.col) > thrX),
                fr_n_cells_above = n_cells_above / n_cells)
  )

  # (v.fr_n_cells_above <- 100* as.named.vector(df_cells_above[3]))
  df_2vec <- df_cells_above[,c(1,4)]
  (v.fr_n_cells_above <- 100* deframe(df_2vec))
  if (label == TRUE) lab =  percentage_formatter(deframe(df_2vec)) else lab = NULL

  pname <- make.names(paste('Cells with', value.col, '>', thrX, id.col))
  ggobj <- ggExpress::qbarplot(v.fr_n_cells_above, suffix = suffix
                    , xlab = 'Clusters', ylab = '% Cells'
                    , plotname = pname, label = lab
                    , subtitle = id.col, xlab.angle = 45)
  if (return.df) return(df_cells_above) else ggobj
}



# _________________________________________________________________________________________________
#' @title scBarplot.FractionBelowThr
#' @description Barplot the fraction of cell below a threshold value (based on a meta.data column), in each cluster.
#' @param thrX PARAM_DESCRIPTION, Default: 0.01
#' @param value.col PARAM_DESCRIPTION, Default: 'percent.ribo'
#' @param id.col PARAM_DESCRIPTION, Default: 'cl.names.top.gene.res.0.3'
#' @param obj Seurat object, Default: combined.obj
#' @param return.df PARAM_DESCRIPTION, Default: F
#' @examples
#' \dontrun{
#' if(interactive()){
#'  scBarplot.FractionBelowThr(id.col =  'cl.names.top.gene.res.0.3', value.col = 'percent.ribo', thrX = 0.01, return.df = T)
#'  }
#' }
#' @seealso
#'  \code{\link[dplyr]{select}}, \code{\link[dplyr]{se-deprecated}}
#' @export
#' @importFrom dplyr select group_by_
scBarplot.FractionBelowThr <- function(thrX = 0.01, value.col = 'percent.ribo', id.col =  'cl.names.top.gene.res.0.3'
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

  pname <- make.names(paste('Cells with', value.col, '<', thrX, id.col))
  ggobj <- ggExpress::qbarplot(v.fr_n_cells_below, xlab = 'Clusters', ylab = '% Cells'
                    , plotname = pname
                    , subtitle = id.col, xlab.angle = 45)
  if (return.df) return(df_cells_below) else ggobj
}





# _________________________________________________________________________________________________

# _________________________________________________________________________________________________

# _________________________________________________________________________________________________
# Read.Write.Save.Load.functions.R
# ____________________________________________________________________ ----
# source('~/GitHub/Packages/Seurat.utils/Functions/Read.Write.Save.Load.functions.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Read.Write.Save.Load.functions.R"))

"Multicore read / write (I/O) functions are https://github.com/vertesy/Seurat.multicore"
"Single core read / write (I/O) functions are in https://github.com/vertesy/Seurat.utils/"


# _________________________________________________________________________________________________
#' @title Convert10Xfolders
#' @description Take a parent directory with a number of subfolders, each containing the standard output of 10X Cell Ranger. (1.) It loads the filtered data matrices; (2.) converts them to Seurat objects, and (3.) saves them as *.RDS files. #
#' @param InputDir Input directory
#' @param regex PARAM_DESCRIPTION, Default: F
#' @param folderPattern PARAM_DESCRIPTION, Default: c("filtered_feature", "SoupX_decont")[1]
#' @param min.cells PARAM_DESCRIPTION, Default: 5
#' @param min.features PARAM_DESCRIPTION, Default: 200
#' @param updateHGNC PARAM_DESCRIPTION, Default: T
#' @param ShowStats PARAM_DESCRIPTION, Default: T
#' @param writeCBCtable write out a list of cell barcodes (CBC) as tsv, Default: T
#' @param depth Depth of scan (How many levels below InputDir). Def 2
#' @param sample.barcoding Cell Ranger run with sample barcoding. The folder structure is different.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  Convert10Xfolders(InputDir)
#'  }
#' }
#' @export
Convert10Xfolders <- function(InputDir # Take a parent directory with a number of subfolders, each containing the standard output of 10X Cell Ranger. (1.) It loads the filtered data matrices; (2.) converts them to Seurat objects, and (3.) saves them as *.RDS files.
                               , regex = F, folderPattern = c("filtered_feature", "SoupX_decont")[1]
                               , min.cells = 5, min.features = 200
                               , updateHGNC = T, ShowStats = T
                               , writeCBCtable = TRUE
                               , sample.barcoding = F
                               , depth = 2) {

  if (sample.barcoding) depth = 3

  finOrig <- list.dirs.depth.n(InputDir, depth = depth)
  fin <- CodeAndRoll2::grepv(x = finOrig, pattern = folderPattern, perl = regex)

  iprint(length(fin), "samples found.")

  if (sample.barcoding) {
    samples <- basename(list.dirs(InputDir, recursive = F))
    iprint("Samples:", samples)
  }

  if (length(fin)) {
    for (i in 1:length(fin)) { print(i)
      pathIN = fin[i]; print(pathIN)

      # sample.barcoding----
      fnameIN = if (sample.barcoding) {
        samples[i]
      } else {
        strsplit(basename(dirname(pathIN)),split = "_")[[1]][1]
      }

      print(fnameIN)

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


      ncells = ncol(seu)
      fnameOUT = ppp(paste0(InputDir, '/', fnameIN), 'min.cells', min.cells, 'min.features', min.features, 'cells', ncells, "Rds")
      print(fnameOUT)


      # update----
      if (updateHGNC) seu <- UpdateGenesSeurat(seu, EnforceUnique = T, ShowStats = T)
      saveRDS(seu, file = fnameOUT)

      # write cellIDs ----
      if (writeCBCtable) {
        fnameCBC <- ppp(fnameOUT, "CBC.tsv")
        CBCs <- t(t(colnames(seu)))
        write.simple.tsv(CBCs, ManualName = fnameCBC)

      }

    }
  } else { iprint("No subfolders found with pattern", folderPattern, "in dirs like: ", finOrig[1:3]) }
}



# _________________________________________________________________________________________________
#' @title Convert10Xfolders.old
#' @description Take a parent directory with a number of subfolders, each containing the standard output of 10X Cell Ranger. (1.) It loads the filtered data matrices; (2.) converts them to Seurat objects, and (3.) saves them as *.RDS files. #
#' @param InputDir Input directory
#' @param folderPattern PARAM_DESCRIPTION, Default: c("filtered", "SoupX_decont")[1]
#' @param min.cells PARAM_DESCRIPTION, Default: 10
#' @param min.features PARAM_DESCRIPTION, Default: 200
#' @param updateHGNC PARAM_DESCRIPTION, Default: T
#' @param ShowStats PARAM_DESCRIPTION, Default: T
#' @examples
#' \dontrun{
#' if(interactive()){
#'  Convert10Xfolders(InputDir = InputDir)
#'  }
#' }
#' @export
Convert10Xfolders.old <- function(InputDir # Take a parent directory with a number of subfolders, each containing the standard output of 10X Cell Ranger. (1.) It loads the filtered data matrices; (2.) converts them to Seurat objects, and (3.) saves them as *.RDS files.
                                  , folderPattern = c("filtered", "SoupX_decont")[1]
                                  , min.cells = 10, min.features = 200, updateHGNC = T, ShowStats = T) {
  fin <- list.dirs(InputDir, recursive = F)
  fin <- CodeAndRoll2::grepv(x = fin, pattern = folderPattern, perl = F)

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


# _________________________________________________________________________________________________
#' @title ConvertDropSeqfolders
#' @description Take a parent directory with a number of subfolders, each containing the standard output of 10X Cell Ranger. (1.) It loads the filtered data matrices; (2.) converts them to Seurat objects, and (3.) saves them as *.RDS files. #
#' @param InputDir Input directory
#' @param folderPattern PARAM_DESCRIPTION, Default: 'SRR*'
#' @param filePattern PARAM_DESCRIPTION, Default: 'expression.tsv.gz'
#' @param useVroom PARAM_DESCRIPTION, Default: T
#' @param col_types.vroom PARAM_DESCRIPTION, Default: list(GENE = "c", .default = "d")
#' @param min.cells PARAM_DESCRIPTION, Default: 10
#' @param min.features PARAM_DESCRIPTION, Default: 200
#' @param updateHGNC PARAM_DESCRIPTION, Default: T
#' @param ShowStats PARAM_DESCRIPTION, Default: T
#' @param minDimension PARAM_DESCRIPTION, Default: 10
#' @param overwrite PARAM_DESCRIPTION, Default: FALSE
#' @examples
#' \dontrun{
#' if(interactive()){
#'  ConvertDropSeqfolders(InputDir)
#'  }
#' }
#' @seealso
#'  \code{\link[vroom]{vroom}}
#'  \code{\link[readr]{read_delim}}
#' @export
#' @importFrom vroom vroom
#' @importFrom readr read_tsv
ConvertDropSeqfolders <- function(InputDir # Take a parent directory with a number of subfolders, each containing the standard output of 10X Cell Ranger. (1.) It loads the filtered data matrices; (2.) converts them to Seurat objects, and (3.) saves them as *.RDS files.
                                  , folderPattern = "SRR*", filePattern = "expression.tsv.gz"
                                  , useVroom = T, col_types.vroom = list("GENE" = "c", .default = "d")
                                  , min.cells = 10, min.features = 200, updateHGNC = T, ShowStats = T, minDimension = 10, overwrite = FALSE) {
  InputDir <- FixPath(InputDir)
  fin <- list.dirs(InputDir, recursive = F)
  fin <- CodeAndRoll2::grepv(x = fin, pattern = folderPattern, perl = F)

  for (i in 1:length(fin)) { print(i)
    pathIN <- FixPath(fin[i]); print(pathIN)
    fnameIN <- basename(fin[i])
    subdir <- paste0(InputDir, fnameIN)
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


# _________________________________________________________________________________________________
#' @title LoadAllSeurats
#' @description Load all Seurat objects found in a directory. Also works with symbolic links (but not with aliases). #
#' @param InputDir Input directory
#' @param file.pattern PARAM_DESCRIPTION, Default: '^filtered.+Rds$'
#' @param string.remove1 PARAM_DESCRIPTION, Default: c(F, "filtered_feature_bc_matrix.", "raw_feature_bc_matrix.")[2]
#' @param string.remove2 PARAM_DESCRIPTION, Default: c(F, ".min.cells.10.min.features.200.Rds")[2]
#' @examples
#' \dontrun{
#' if(interactive()){
#'  ls.Seurat <- LoadAllSeurats(InputDir)
#'  }
#' }
#' @export
#' @importFrom tictoc tic toc
LoadAllSeurats <- function(InputDir # Load all Seurat objects found in a directory. Also works with symbolic links (but not with aliases).
                           , file.pattern = "^filtered.+Rds$"
                           , string.remove1 = c(F, "filtered_feature_bc_matrix.", "raw_feature_bc_matrix." )[2]
                           , string.remove2 = c(F, ".min.cells.10.min.features.200.Rds")[2]) {
  tictoc::tic()
  InputDir <- AddTrailingSlash(InputDir) # add '/' if necessary

  fin.orig <- list.files(InputDir, include.dirs = F, pattern = file.pattern)
  print(fin.orig)
  fin <- if (!isFALSE(string.remove1)) sapply(fin.orig, gsub, pattern = string.remove1, replacement = "") else fin.orig
  fin <- if (!isFALSE(string.remove2)) sapply(fin, gsub, pattern = string.remove2, replacement = "") else fin

  ls.Seu <- list.fromNames(fin)
  for (i in 1:length(fin)) {print(fin[i]); ls.Seu[[i]] <- readRDS(paste0(InputDir, fin.orig[i]))}
  print(tictoc::toc())
  return(ls.Seu)
}



# _________________________________________________________________________________________________
#' @title read10x
#' @description read10x from gzipped matrix.mtx, features.tsv and barcodes.tsv #
#' @param dir PARAM_DESCRIPTION
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[tictoc]{tic}}
#'  \code{\link[R.utils]{compressFile}}
#'  \code{\link[Seurat]{Read10X}}
#' @export
#' @importFrom tictoc tic toc
#' @importFrom R.utils gunzip gzip
#' @importFrom Seurat Read10X
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


# _________________________________________________________________________________________________
#' @title saveRDS.compress.in.BG
#' @description Save and RDS object and compress resulting file in the background using system(gzip). OS X or unix.
#' @param obj Seurat object.
#' @param compr PARAM_DESCRIPTION, Default: FALSE
#' @param fname File name
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[tictoc]{tic}}
#' @export
#' @importFrom tictoc tic toc
saveRDS.compress.in.BG <- function(obj, compr = FALSE, fname) {
  try(tictoc::tic(), silent = T)
  saveRDS(object = obj, compress = compr, file = fname)
  try(tictoc::toc(), silent = T)
  print(paste("Saved, being compressed", fname))
  system(command = paste0("gzip '", fname,"'"),  wait = FALSE) # execute in the background
  try(say(), silent = T)
}



#' getProject  -----------------------------------------------
#'
#' @description Try to get the project name you are wokring on in Rstudio.
#' @returns The final subfolder of your project, or NULL, if you are not running one
#' @export
#'
#' @examples getProject()
getProject <- function() {
  tryCatch(basename(rstudioapi::getActiveProject()), error=function(e){})
}


# Save an object -----------------------------------------------
#' @title isave.RDS
#' @description Save and RDS object.
#' @param obj Seurat object
#' @param project project code appended to the saved file name. Default: try(basename(rstudioapi::getActiveProject()), silent=T) using  getProject().
#' @param prefix PARAM_DESCRIPTION, Default: NULL
#' @param suffix A suffix added to the filename, Default: NULL
#' @param inOutDir PARAM_DESCRIPTION, Default: F
#' @param alternative_path_rdata PARAM_DESCRIPTION, Default: paste0("~/Dropbox/Abel.IMBA/AnalysisD/_RDS.files/", basename(OutDir))
#' @param homepath homepath to replace '~', Default: paste0("~/Dropbox/Abel.IMBA/AnalysisD/_RDS.files/", basename(OutDir))
#' @param showMemObject PARAM_DESCRIPTION, Default: T
#' @param saveParams PARAM_DESCRIPTION, Default: T
#' @examples
#' \dontrun{
#' if(interactive()){
#'  isave.RDS(my.R.object)
#'  }
#' }
#' @export

isave.RDS <- function(obj, prefix =NULL, suffix = NULL, inOutDir = F
                      , project = getProject()
                      , alternative_path_rdata = paste0("~/Dropbox (VBC)/Abel.IMBA/AnalysisD/_RDS.files/", basename(OutDir))
                      , homepath = '/Users/abel.vertesy/'
                      , showMemObject = T, saveParams =T){ # Faster saving of workspace, and compression outside R, when it can run in the background. Seemingly quite CPU hungry and not very efficient compression.
  path_rdata = if (inOutDir) OutDir else alternative_path_rdata
  dir.create(path_rdata)

  if (showMemObject) { try(memory.biggest.objects(), silent = T) }
  if ( "seurat" %in% is(obj) & saveParams) {
    try(obj@misc$p <- p, silent = T)
    try(obj@misc$all.genes  <- all.genes, silent = T)
  }
  fnameBase = kppu(prefix, substitute(obj), project, suffix, idate(Format = "%Y.%m.%d_%H.%M"))
  fnameBase = trimws(fnameBase, whitespace = '_')
  FNN <- paste0(path_rdata, fnameBase , ".Rds")
  FNN <- gsub(pattern = '~/', replacement = homepath, x = FNN)
  print(FNN)
  saveRDS.compress.in.BG(obj = obj, fname =  FNN)
}




# Save workspace -----------------------------------------------
# requires MarkdownReports (github) and defining OutDir
# requires github/vertesy/CodeAndRoll.r

#' @title isave.image
#' @description Save and RData image.
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @param path_rdata PARAM_DESCRIPTION, Default: paste0("~/Dropbox/Abel.IMBA/AnalysisD/_Rdata.files/", basename(OutDir))
#' @param showMemObject PARAM_DESCRIPTION, Default: T
#' @param options PARAM_DESCRIPTION, Default: c("--force", NULL)[1]
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[Stringendo]{kollapse}}, \code{\link[MarkdownReports]{iprint}}
#' @export
#' @importFrom Stringendo kollapse iprint
isave.image <- function(..., path_rdata = paste0("~/Dropbox/Abel.IMBA/AnalysisD/_Rdata.files/", basename(OutDir))
                        , showMemObject = T, options = c("--force", NULL)[1]
){ # Faster saving of workspace, and compression outside R, when it can run in the background. Seemingly quite CPU hungry and not very efficient compression.

  dir.create(path_rdata)

  if (showMemObject) { try(memory.biggest.objects(), silent = T) }
  fname = Stringendo::kollapse(path_rdata, "/",idate(),...,".Rdata")
  print(fname)
  if (nchar(fname) > 2000) stop()
  save.image(file = fname, compress = F)
  iprint("Saved, being compressed", fname)
  system(paste("gzip", options, fname),  wait = FALSE) # execute in the background
}


# Save workspace -----------------------------------------------
# requires MarkdownReports (github) and defining OutDir
# requires github/vertesy/CodeAndRoll.r

#' @title qsave.image
#' @description Faster saving of workspace, and compression outside R, when it can run in the background. Seemingly quite CPU hungry and not very efficient compression. #
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @param showMemObject PARAM_DESCRIPTION, Default: T
#' @param options PARAM_DESCRIPTION, Default: c("--force", NULL)[1]
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[Stringendo]{kollapse}}, \code{\link[function]{iprint}}
#' @export
#' @importFrom Stringendo kollapse iprint
#' @importFrom tictoc tic toc
qsave.image <- function(..., showMemObject = T, options = c("--force", NULL)[1]){ # Faster saving of workspace, and compression outside R, when it can run in the background. Seemingly quite CPU hungry and not very efficient compression.
  fname = Stringendo::kollapse(getwd(), "/",basename(OutDir),idate(),...,".Rdata")
  print(fname)
  if (nchar(fname) > 2000) stop()
  tictoc::tic()
  save.image(file = fname, compress = F)
  iprint("Saved, being compressed", fname)
  system(paste("gzip", options, fname),  wait = FALSE) # execute in the background
  cat(toc)
}


# subsetSeuObj -----------------------------------------------------------------------
#' @title subsetSeuObj
#' @description Subset a compressed Seurat Obj and save it in wd. #
#' @param obj Seurat object, Default: ls.Seurat[[i]]
#' @param fraction_ PARAM_DESCRIPTION, Default: 0.25
#' @param nCells PARAM_DESCRIPTION, Default: F
#' @param seed_ PARAM_DESCRIPTION, Default: 1989
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
#' @importFrom Stringendo percentage_formatter
subsetSeuObj <- function(obj = ls.Seurat[[i]], fraction_ = 0.25, nCells = F, seed_ = 1989 ) { # Subset a compressed Seurat Obj and save it in wd.
  set.seed(seed_)
  if (isFALSE(nCells)) {
    cellIDs.keep = sampleNpc(metaDF = obj@meta.data, pc = fraction_)
    iprint(length(cellIDs.keep), "or",Stringendo::percentage_formatter(fraction_),"of the cells are kept. Seed:", seed_)
  } else if (nCells > 1) {
    nKeep = min(ncol(obj), nCells)
    # print(nKeep)
    cellIDs.keep = sample(colnames(obj), size = nKeep, replace = F)
    if (nKeep < nCells) iprint("Only",nCells,"cells were found in the object, so downsampling is not possible.")
  }
  obj <- subset(x = obj, cells = cellIDs.keep) # downsample
  return(obj)
}

# _________________________________________________________________________________________________
#' @title subsetSeuObj.and.Save
#' @description Subset a compressed Seurat Obj and save it in wd. #
#' @param obj Seurat object, Default: ORC
#' @param fraction PARAM_DESCRIPTION, Default: 0.25
#' @param seed random seed used, Default: 1989
#' @param min.features Minimum features
#' @param dir PARAM_DESCRIPTION, Default: OutDir
#' @param suffix A suffix added to the filename, Default: ''
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
subsetSeuObj.and.Save <- function(obj = ORC, fraction = 0.25, seed = 1989, dir = OutDir
                                  , min.features = p$'min.features', suffix = '') { # Subset a compressed Seurat Obj and save it in wd.
  obj_Xpc <- subsetSeuObj(obj = obj, fraction_ =  fraction, seed_ = seed)
  nr.cells.kept <- ncol(obj_Xpc)
  saveRDS.compress.in.BG(obj = obj_Xpc, fname = ppp(paste0(dir, substitute(obj)),suffix, nr.cells.kept, 'cells.with.min.features', min.features,"Rds" ) )
}


# _________________________________________________________________________________________________
#' @title Downsample.Seurat.Objects
#' @description Downsample a list of Seurat objects
#' @param ls.obj List of Seurat objects, Default: ls.Seurat
#' @param NrCells PARAM_DESCRIPTION, Default: p$dSample.Organoids
#' @examples
#' \dontrun{
#' if(interactive()){
#'  Downsample.Seurat.Objects(NrCells = 2000); Downsample.Seurat.Objects(NrCells = 200)
#'  }
#' }
#' @export
#' @importFrom tictoc tic toc
#' @importFrom Stringendo percentage_formatter
Downsample.Seurat.Objects <- function(ls.obj = ls.Seurat, NrCells = p$"dSample.Organoids") {
  names.ls = names(ls.obj)
  n.datasets = length(ls.obj)
  iprint(NrCells, "cells")
  tictoc::tic()
  if (foreach::getDoParRegistered() ) {
    ls.obj.downsampled <- foreach(i = 1:n.datasets ) %dopar% {
      iprint(names(ls.obj)[i], Stringendo::percentage_formatter(i/n.datasets, digitz = 2))
      subsetSeuObj(obj = ls.obj[[i]], nCells = NrCells)
    }; names(ls.obj.downsampled)  <- names.ls
  } else {
    ls.obj.downsampled <- list.fromNames(names.ls)
    for (i in 1:n.datasets ) {
      iprint(names(ls.obj)[i], Stringendo::percentage_formatter(i/n.datasets, digitz = 2))
      ls.obj.downsampled[[i]] <- subsetSeuObj(obj = ls.obj[[i]], nCells = NrCells)
    };
  } # else
  toc();

  print(head(unlapply(ls.obj, ncol)))
  print(head(unlapply(ls.obj.downsampled, ncol)))

  isave.RDS(obj = ls.obj.downsampled, suffix = ppp(NrCells, "cells"), inOutDir = T)

}



# _________________________________________________________________________________________________
#' @title Downsample.Seurat.Objects.PC
#' @description Downsample a list of Seurat objects, by fraction
#' @param ls.obj List of Seurat objects, Default: ls.Seurat
#' @param NrCells PARAM_DESCRIPTION, Default: p$dSample.Organoids
#' @examples
#' \dontrun{
#' if(interactive()){
#'  Downsample.Seurat.Objects.PC()
#'  }
#' }
#' @export
#' @importFrom tictoc tic toc
#' @importFrom Stringendo percentage_formatter

Downsample.Seurat.Objects.PC <- function(ls.obj = ls.Seurat, fraction = 0.1) {
  names.ls = names(ls.obj)
  n.datasets = length(ls.obj)
  iprint(fraction, "fraction")
  tictoc::tic()
  if (foreach::getDoParRegistered() ) {
    ls.obj.downsampled <- foreach(i = 1:n.datasets ) %dopar% {
      subsetSeuObj(obj = ls.obj[[i]], fraction_ = fraction)
    }; names(ls.obj.downsampled)  <- names.ls
  } else {
    ls.obj.downsampled <- list.fromNames(names.ls)
    for (i in 1:n.datasets ) {
      cells = round(ncol(ls.obj[[1]]) * fraction)
      iprint(names(ls.obj)[i], cells, "cells=", Stringendo::percentage_formatter(i/n.datasets, digitz = 2))
      ls.obj.downsampled[[i]] <- subsetSeuObj(obj = ls.obj[[i]], fraction_ = fraction)
    };
  }; toc(); # else

  NrCells <- sum(unlapply(ls.obj, ncol))

  print(head(unlapply(ls.obj, ncol)))
  print(head(unlapply(ls.obj.downsampled, ncol)))
  isave.RDS(obj = ls.obj.downsampled, suffix = ppp(NrCells, "cells"), inOutDir = T)

}


# _________________________________________________________________________________________________
# Seurat.object.manipulations.etc.R
# ____________________________________________________________________ ----
# source('~/GitHub/Packages/Seurat.utils/Functions/Seurat.object.manipulations.etc.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Seurat.object.manipulations.etc.R"))

# _________________________________________________________________________________________________
#' @title clip10Xcellname
#' @description Clip all suffices after underscore (10X adds it per chip-lane, Seurat adds in during integration). #
#' @param cellnames PARAM_DESCRIPTION
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
clip10Xcellname <- function(cellnames) stringr::str_split_fixed(cellnames, "_", n = 2)[,1] # Clip all suffices after underscore (10X adds it per chip-lane, Seurat adds in during integration).

# _________________________________________________________________________________________________
#' @title make10Xcellname
#' @description Add a suffix to cell names, so that it mimics the lane-suffix, e.g.: "_1". #
#' @param cellnames PARAM_DESCRIPTION
#' @param suffix A suffix added to the filename, Default: '_1'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
make10Xcellname <- function(cellnames, suffix="_1") paste0(cellnames, suffix) # Add a suffix to cell names, so that it mimics the lane-suffix, e.g.: "_1".


# _________________________________________________________________________________________________
#' @title seu.Make.Cl.Label.per.cell
#' @description Take a named vector (of e.g. values ="gene names", names = clusterID), and a vector of cell-IDs and make a vector of "GeneName.ClusterID". #
#' @param TopGenes PARAM_DESCRIPTION
#' @param clID.per.cell PARAM_DESCRIPTION
#' @examples
#' \dontrun{
#' if(interactive()){
#'  seu.Make.Cl.Label.per.cell(TopGenes = TopGenes.Classic, clID.per.cell = getMetadataColumn(ColName.metadata = metaD.CL.colname)  )
#'  }
#' }
#' @export
seu.Make.Cl.Label.per.cell <- function(TopGenes, clID.per.cell) { # Take a named vector (of e.g. values ="gene names", names = clusterID), and a vector of cell-IDs and make a vector of "GeneName.ClusterID".
  Cl.names_class= TopGenes[ clID.per.cell ]
  Cl.names_wNr = paste0(Cl.names_class,' (',names(Cl.names_class),')')
  return(Cl.names_wNr)
}


# FeaturePlot with different defaults ------------------------------------------------------------------
#' @title GetMostVarGenes
#' @description Get the most variable rGenes #
#' @param obj Seurat object, Default: org
#' @param nGenes PARAM_DESCRIPTION, Default: p$nVarGenes
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
GetMostVarGenes <- function(obj = org, nGenes = p$nVarGenes) { # Get the most variable rGenes
  head(rownames(slot(object = obj, name = "hvg.info")), n = nGenes)
}

# gene.name.check for read .mtx /write .rds script ---------------------------------------
#' @title gene.name.check
#' @description Check gene names in a seurat object, for naming conventions (e.g.: mitochondrial reads have - or .). Use for reading .mtx & writing .rds files. #
#' @param Seu.obj PARAM_DESCRIPTION, Default: ls.Seurat[[1]]
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
gene.name.check <- function(Seu.obj = ls.Seurat[[1]] ) { # Check gene names in a seurat object, for naming conventions (e.g.: mitochondrial reads have - or .). Use for reading .mtx & writing .rds files.
  rn = rownames(GetAssayData(object = Seu.obj, slot = "counts"))
  llprint("### Gene name pattern")

  llogit('`rn = rownames(GetAssayData(object = ls.Seurat[[1]], slot = "counts"))`')
  llogit('`head(CodeAndRoll2::grepv(rn, pattern = "-"), 10)`')
  print('pattern = -')
  llprint(head(CodeAndRoll2::grepv(rn, pattern = "-"), 10))

  llogit('`head(CodeAndRoll2::grepv(rn, pattern = "_"), 10)`')
  print('pattern = _')
  llprint(head(CodeAndRoll2::grepv(rn, pattern = "_"), 10))

  llogit('`head(CodeAndRoll2::grepv(rn, pattern = "\\."), 10)`')
  print('pattern = \\.')
  llprint(head(CodeAndRoll2::grepv(rn, pattern = "\\."), 10))

  llogit('`head(CodeAndRoll2::grepv(rn, pattern = "\\.AS[1-9]"), 10)`')
  print('pattern = \\.AS[1-9]')
  llprint(head(CodeAndRoll2::grepv(rn, pattern = "\\.AS[1-9]"), 10))
}


# _________________________________________________________________________________________________
#' @title check.genes
#' @description Check if a gene name exists in a Seurat object, or in HGNC?
#' @param list.of.genes PARAM_DESCRIPTION, Default: ClassicMarkers
#' @param makeuppercase PARAM_DESCRIPTION, Default: FALSE
#' @param verbose PARAM_DESCRIPTION, Default: TRUE
#' @param HGNC.lookup PARAM_DESCRIPTION, Default: FALSE
#' @param obj Seurat object, Default: combined.obj
#' @param assay.slot PARAM_DESCRIPTION, Default: c("RNA", "integrated")[1]
#' @param dataslot PARAM_DESCRIPTION, Default: c("counts", "data")[2]
#' @examples
#' \dontrun{
#' if(interactive()){
#'  check.genes("top2a", makeuppercase = TRUE); check.genes("VGLUT2", verbose = F, HGNC.lookup = T)
#'  }
#' }
#' @export
#' @importFrom Stringendo percentage_formatter
check.genes <- function(list.of.genes = ClassicMarkers, makeuppercase = FALSE, verbose  = TRUE, HGNC.lookup = FALSE
                        , obj = combined.obj
                        , assay.slot = c('RNA', 'integrated')[1]
                        , dataslot = c("counts", "data")[2]) { # Check if genes exist in your dataset.
  if (makeuppercase) list.of.genes <- toupper(list.of.genes)
  all_genes = rownames(GetAssayData(object = obj, assay = assay.slot, slot = dataslot)); length(all_genes)
  missingGenes = setdiff(list.of.genes, all_genes)
  if (length(missingGenes) > 0) {
    if (verbose) { iprint(length(missingGenes), "or", Stringendo::percentage_formatter(length(missingGenes) / length(list.of.genes)), "genes not found in the data, e.g:", head(missingGenes, n = 10))  }
    if (HGNC.lookup) {
      if (exists('qHGNC', mode='function')) { try(DatabaseLinke.R::qHGNC(missingGenes)) } else { print("load qHGNC() function, see database.linker")}
    }
  }
  intersect(list.of.genes, all_genes)
}



# _________________________________________________________________________________________________
#' @title fixZeroIndexing.seurat
#' @description Fix zero indexing in seurat clustering, to 1-based indexing. replace zero indexed clusternames.
#' @param ColName.metadata PARAM_DESCRIPTION, Default: 'res.0.6'
#' @param obj Seurat object, Default: org
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
fixZeroIndexing.seurat <- function(ColName.metadata = 'res.0.6', obj = org) { # Fix zero indexing in seurat clustering, to 1-based indexing
  obj@meta.data[ ,ColName.metadata] =  as.numeric(obj@meta.data[ ,ColName.metadata])+1
  print(obj@meta.data[ ,ColName.metadata])
  return(obj)
}


# _________________________________________________________________________________________________
#' @title CalculateFractionInTrome
#' @description Calculate the fraction of a set of genes within the full Transcriptome of each cell. #
#' @param geneset PARAM_DESCRIPTION, Default: c("MALAT1")
#' @param obj Seurat object, Default: combined.obj
#' @param dataslot PARAM_DESCRIPTION, Default: c("counts", "data")[2]
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
CalculateFractionInTrome <- function(genesCalc.Cor.Seuratet = c("MALAT1") # Calculate the fraction of a set of genes within the full Transcriptome of each cell.
                                     , obj = combined.obj
                                     , dataslot = c("counts", "data")[2]
) {
  print("    >>>> Use add.meta.fraction() <<<<")
  geneset <- check.genes(list.of.genes = geneset)
  stopifnot(length(geneset)>0)

  mat <- as.matrix(slot(obj@assays$RNA, name = dataslot))
  mat.sub <- mat[geneset,,drop = F]
  RC.per.cell.geneset <- colSums(mat.sub)

  RC.per.cell <- colSums(mat)
  gene.fraction.per.cell <- 100*RC.per.cell.geneset / RC.per.cell
  return(gene.fraction.per.cell)
}

# _________________________________________________________________________________________________
#' @title AddNewAnnotation
#' @description Create a new metadata column based on an exisiting metadata column and a list of mappings (name <- IDs). #
#' @param obj Seurat object, Default: obj
#' @param source PARAM_DESCRIPTION, Default: 'RNA_snn_res.0.5'
#' @param named.list.of.identities PARAM_DESCRIPTION, Default: ls.Subset.ClusterLists
#' @examples
#' \dontrun{
#' if(interactive()){
#'  ls.Subset.ClusterLists = list( "hESC.h9" = c("4", "10", "14"), "hESC.176" = c("0", "1", "2")); AddNewAnnotation()
#'  }
#' }
#' @export
AddNewAnnotation <- function(obj = obj # Create a new metadata column based on an exisiting metadata column and a list of mappings (name <- IDs).
                             , source = "RNA_snn_res.0.5", named.list.of.identities = ls.Subset.ClusterLists) {
  NewID <- as.named.vector(obj[[source]])

  for (i in 1:length(named.list.of.identities)) {
    lx <- as.character(named.list.of.identities[[i]])
    name.lx <- names(named.list.of.identities)[i]
    NewID <- translate(vec = NewID, oldvalues = lx, newvalues = name.lx)
  }
  print(table(NewID))
  return(NewID)
}


# _________________________________________________________________________________________________
#' @title whitelist.subset.ls.Seurat
#' @description Subset cells in a (list of) Seurat objects, based on an externally provided list of cell IDs.
#' @param ls.obj List of Seurat objects, Default: ls.Seurat
#' @param metadir PARAM_DESCRIPTION, Default: p$cellWhiteList
#' @param whitelist.file PARAM_DESCRIPTION, Default: 'NonStressedCellIDs.2020.10.21_18h.tsv'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
whitelist.subset.ls.Seurat <- function(ls.obj = ls.Seurat
                                       , metadir = p$'cellWhiteList' #  '~/Dropbox/Abel.IMBA/MetadataD/POL.meta/cell.lists/'
                                       , whitelist.file = "NonStressedCellIDs.2020.10.21_18h.tsv"
) {
  cells.before <- unlapply(ls.obj, ncol)
  # Find file
  df.cell.whitelist <- ReadWriter::read.simple.tsv(metadir, whitelist.file)
  dsets <- table(df.cell.whitelist[,1])

  ls.orig.idents <- lapply(lapply(ls.Seurat, getMetadataColumn, ColName.metadata = "orig.ident"), unique)
  stopif(any(unlapply(ls.orig.idents, l) == length(ls.Seurat)), message = "Some ls.Seurat objects have 1+ orig identity.")

  dsets.in.lsSeu <- unlist(ls.orig.idents)
  isMathced <- all(dsets.in.lsSeu == names(dsets)) # Stop if either ls.Seurat OR the metadata has identities not found in the other, in the same order.
  stopif(!isMathced, message = paste("either ls.Seurat OR the metadata has identities not found in the other, or they are not in same order."
                                     , kpps(dsets.in.lsSeu),"vs.", kpps(names(dsets) ) )
  )

  # identX <- ls.orig.idents[[1]]
  for (i in 1:length(ls.orig.idents)) {
    identX <- ls.orig.idents[[i]]; print(identX)

    # Extract and process cellIDs
    idx.match <- which(df.cell.whitelist[,1] == identX)
    cell.whitelist <- rownames(df.cell.whitelist)[idx.match]
    cell.whitelist <- substr(x = cell.whitelist
                             , start = 1 ,stop = nchar(cell.whitelist)-2)

    # Extract and process cellIDs
    ls.obj[[i]] <- subset(x = ls.obj[[i]], cells = cell.whitelist)
  }
  cells.after <- unlapply(ls.obj, ncol)
  iprint("cells.before",cells.before,"cells.after",cells.after)
  return(ls.obj)
}

# _________________________________________________________________________________________________
#' @title FindCorrelatedGenes
#' @description Find correlated genes in a Seurat object
#' @param gene gene of interest, Default: 'TOP2A'
#' @param obj Seurat object, Default: combined.obj
#' @param assay RNA or integrated assay, Default: 'RNA'
#' @param slot slot in the Seurat object. Default: 'data'
#' @param HEonly PARAM_DESCRIPTION, Default: F
#' @param minExpr PARAM_DESCRIPTION, Default: 1
#' @param minCells PARAM_DESCRIPTION, Default: 1000
#' @param trailingNgenes PARAM_DESCRIPTION, Default: 1000
#' @examples
#' \dontrun{
#' if(interactive()){
#'  FindCorrelatedGenes(gene ="TOP2A", obj = combined.obj); write_clip(names(head(topGenes[-(1:6)], n = 50)))
#'  }
#' }
#' @seealso
#'  \code{\link[matrixStats]{rowSums2}}
#' @export
#' @importFrom matrixStats rowSums2
#' @importFrom tictoc tic toc
FindCorrelatedGenes <- function(gene ="TOP2A", obj = combined.obj, assay = "RNA", slot = "data"
                                , HEonly =F , minExpr = 1, minCells = 1000
                                , trailingNgenes = 1000) {
  tictoc::tic()
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
  tictoc::toc()
  wbarplot(head(topGenes, n =25))
  topGenes
}





# _________________________________________________________________________________________________
# _________________________________________________________________________________________________

# _________________________________________________________________________________________________
# Seurat.update.gene.symbols.HGNC.R
# ____________________________________________________________________ ----
# source('~/GitHub/Packages/Seurat.utils/Functions/Seurat.update.gene.symbols.HGNC.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Seurat.update.gene.symbols.HGNC.R"))
# require(HGNChelper)

# _________________________________________________________________________________________________
#' @title UpdateGenesSeurat
#' @description Update genes symbols that are stored in a Seurat object. It returns a data frame. The last column are the updated gene names. #
#' @param obj Seurat object, Default: ls.Seurat[[i]]
#' @param species_ PARAM_DESCRIPTION, Default: 'human'
#' @param EnforceUnique PARAM_DESCRIPTION, Default: T
#' @param ShowStats PARAM_DESCRIPTION, Default: F
#' @examples
#' \dontrun{
#' if(interactive()){
#'  UpdateGenesSeurat()
#'  }
#' }
#' @seealso
#'  \code{\link[HGNChelper]{checkGeneSymbols}}
#' @export
#' @importFrom HGNChelper checkGeneSymbols
UpdateGenesSeurat <- function(obj = ls.Seurat[[i]], species_="human", EnforceUnique = T, ShowStats = F ) { # Update genes symbols that are stored in a Seurat object. It returns a data frame. The last column are the updated gene names.
  HGNC.updated <- HGNChelper::checkGeneSymbols(rownames(obj), unmapped.as.na = FALSE, map = NULL, species = species_)
  if (EnforceUnique) HGNC.updated <- HGNC.EnforceUnique(HGNC.updated)
  if (ShowStats) print(GetUpdateStats(HGNC.updated))
  obj <- RenameGenesSeurat(obj, newnames = HGNC.updated$Suggested.Symbol)
  return(obj)
}


# _________________________________________________________________________________________________
#' @title RenameGenesSeurat
#' @description Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data. #
#' @param obj Seurat object, Default: ls.Seurat[[i]]
#' @param newnames PARAM_DESCRIPTION, Default: HGNC.updated[[i]]$Suggested.Symbol
#' @examples
#' \dontrun{
#' if(interactive()){
#'  RenameGenesSeurat(obj = SeuratObj, newnames = HGNC.updated.genes$Suggested.Symbol)
#'  }
#' }
#' @export
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



# _________________________________________________________________________________________________
#' @title RemoveGenesSeurat
#' @description Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes metadata; obj@assays$RNA@counts, @data and @scale.data. #
#' @param obj Seurat object, Default: ls.Seurat[[i]]
#' @param symbols2remove PARAM_DESCRIPTION, Default: c("TOP2A")
#' @examples
#' \dontrun{
#' if(interactive()){
#'  RemoveGenesSeurat(obj = SeuratObj, symbols2remove = "TOP2A")
#'  }
#' }
#' @export
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



# _________________________________________________________________________________________________
#' @title HGNC.EnforceUnique
#' @description Enforce Unique names after HGNC symbol update. updatedSymbols is the output of HGNChelper::checkGeneSymbols. #
#' @param updatedSymbols PARAM_DESCRIPTION
#' @examples
#' \dontrun{
#' if(interactive()){
#'  x <- HGNC.EnforceUnique(updatedSymbols = SymUpd)
#'  # While "make.unique" is not the ideal solution, because it generates mismatched, in my integration example it does reduce the mismatching genes from ~800 to 4
#'  }
#' }
#' @export
HGNC.EnforceUnique <- function(updatedSymbols) { # Enforce Unique names after HGNC symbol update. updatedSymbols is the output of HGNChelper::checkGeneSymbols.
  NGL <- updatedSymbols[,3]
  if (any.duplicated(NGL)) {
    updatedSymbols[,3] <- make.unique(NGL); "Unique names are enforced by suffixing .1, .2, etc."
  }
  return(updatedSymbols)
}




# _________________________________________________________________________________________________
#' @title GetUpdateStats
#' @description Plot the Symbol-update statistics. Works on the data frame returned by `UpdateGenesSeurat()`. #
#' @param genes Genes of iinterest, Default: HGNC.updated[[i]]
#' @examples
#' \dontrun{
#' if(interactive()){
#'  GetUpdateStats(genes = HGNC.updated.genes)
#'  }
#' }
#' @export
#' @importFrom Stringendo percentage_formatter
GetUpdateStats <- function(genes = HGNC.updated[[i]]) { # Plot the Symbol-update statistics. Works on the data frame returned by `UpdateGenesSeurat()`.
  (MarkedAsUpdated <- genes[genes$Approved == FALSE, ])
  (AcutallyUpdated <- sum(MarkedAsUpdated[,1] != MarkedAsUpdated[,3]))
  (UpdateStats = c("Updated (%)"=Stringendo::percentage_formatter(AcutallyUpdated / nrow(genes)), "Updated Genes"=floor(AcutallyUpdated), "Total Genes"=floor(nrow(genes))))
  return(UpdateStats)
}


# _________________________________________________________________________________________________
#' @title PlotUpdateStats
#' @description Scatter plot of update stats. #
#' @param mat PARAM_DESCRIPTION, Default: UpdateStatMat
#' @param column.names PARAM_DESCRIPTION, Default: c("Updated (%)", "Updated (Nr.)")
#' @examples
#' \dontrun{
#' if(interactive()){
#'  PlotUpdateStats(mat = result.of.GetUpdateStats)
#'  }
#' }
#' @export
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


# _________________________________________________________________________________________________
# Soup.Analysis.of.ambient.RNA.R
# ____________________________________________________________________ ----
# source('~/GitHub/Packages/Seurat.utils/Functions/Soup.Analysis.of.ambient.RNA.R')
# try (source('https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Soup.Analysis.of.ambient.RNA.R'))
# Source: self + web


# require(tibble)

# _________________________________________________________________________________________________
#' @title plotTheSoup
#' @description Plot stats about the ambient RNA content in a 10X experiment.
#' @param CellRangerOutputDir PARAM_DESCRIPTION, Default: '~/Data/114593/114593'
#' @param SeqRun PARAM_DESCRIPTION, Default: gsub("*([0-9]+).*", "\\1", x = basename(CellRangerOutputDir))
#' @examples
#' \dontrun{
#' if(interactive()){
#'  plotTheSoup(CellRangerOutputDir = "~/Data/114593/114593" , SeqRun = gsub('*([0-9]+).*','\\1', x = basename(CellRangerOutputDir)))
#'  }
#' }
#' @seealso
#'  \code{\link[Matrix]{colSums}}
#'  \code{\link[tibble]{rownames}}
#'  \code{\link[ggrepel]{geom_label_repel}}
#' @export
#' @importFrom Matrix rowSums
#' @importFrom tibble rownames_to_column
#' @importFrom ggrepel geom_text_repel
#' @importFrom Stringendo percentage_formatter
plotTheSoup <- function(CellRangerOutputDir = "~/Data/114593/114593"
                        , SeqRun = gsub('*([0-9]+).*','\\1', x = basename(CellRangerOutputDir))) { # Plot the ambient RNA content of droplets without a cell (background droplets).

  ls.Alpha = 1
  # Setup ___________________________________
  # require(Matrix); require(ggrepel)

  dirz <- list.dirs(CellRangerOutputDir, full.names = F, recursive = F)
  path.raw <- file.path(CellRangerOutputDir, grep(x = dirz, pattern = "^raw_*", value = T))
  path.filt <- file.path(CellRangerOutputDir, grep(x = dirz, pattern = "^filt_*", value = T))
  CR.matrices <- list.fromNames(c("raw", "filt"))

  # Adapter for Markdownreports background variable "OutDir ___________________________________
  if (exists('OutDir')) OutDirBac <- OutDir
  OutDir <- file.path(CellRangerOutputDir,paste0(kpp("SoupStatistics", SeqRun)))
  try(dir.create(OutDir))
  ww.assign_to_global("OutDir", OutDir, 1)

  # Read In ___________________________________
  print("Reading raw CellRanger output matrices")
  CR.matrices$'raw' <- Read10X(path.raw)
  if (length(CR.matrices$'raw') == 2 ) { CR.matrices$'raw' <- CR.matrices$'raw'[[1]] } # Maybe AB table is present too at slot 2!

  print("Reading filtered CellRanger output matrices")
  CR.matrices$'filt' <- Read10X(path.filt)
  if (length(CR.matrices$'filt') == 2 ) { CR.matrices$'filt' <- CR.matrices$'filt'[[1]] } # Maybe AB table is present too at slot 2!

  # Profiling the soup ___________________________________
  print("Profiling the soup")
  GEMs.all <- CR.matrices$'raw'@Dimnames[[2]]
  GEMs.cells <- CR.matrices$'filt'@Dimnames[[2]]
  iprint("There are", length(GEMs.all), "GEMs sequenced, and",l(GEMs.cells), "are cells among those." )

  GEMs.soup <- setdiff(GEMs.all, GEMs.cells)
  CR.matrices$'soup' <- CR.matrices$'raw'[,GEMs.soup]
  CR.matrices$'soup.total.RC' <- Matrix::rowSums(CR.matrices$'soup')
  CR.matrices$'soup.total.sum' <- sum(CR.matrices$'soup')
  CR.matrices$'cells.total.sum' <- sum(CR.matrices$'filt')

  CR.matrices$'soup.rel.RC'  <- CR.matrices$'soup.total.RC' / CR.matrices$'soup.total.sum'

  # Diff Exp ___________________________________
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

  # ggplot prepare ___________________________________
  Soup.VS.Cells.Av.Exp.gg <- tibble::rownames_to_column(as.data.frame(Soup.VS.Cells.Av.Exp.log10), "gene")
  (Soup.VS.Cells.Av.Exp.gg <- as_tibble(Soup.VS.Cells.Av.Exp.gg))
  soup.rate <- Soup.VS.Cells.Av.Exp.gg$Soup / (Soup.VS.Cells.Av.Exp.gg$Cells + Soup.VS.Cells.Av.Exp.gg$Soup)
  cell.rate <- Soup.VS.Cells.Av.Exp.gg$Cells / (Soup.VS.Cells.Av.Exp.gg$Cells + Soup.VS.Cells.Av.Exp.gg$Soup)

  axl.pfx <- "Total Expression in"
  axl.sfx <- "[log10(mRNA+1)]"



  HGNC <- Soup.VS.Cells.Av.Exp.gg$gene
  Class <- rep("Other", times = nrow(Soup.VS.Cells.Av.Exp.gg))
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
           , aes(x= Soup, y= Cells, label = gene, col= Class))  +
    geom_abline(slope = 1, col='darkgrey') + geom_point()+
    scale_alpha_manual(guide='none', values = ls.Alpha) +
    xlab(paste(axl.pfx, "Soup", axl.sfx)) + ylab(paste(axl.pfx, "Cells", axl.sfx)) +
    ggtitle("Soup VS. Cells | gene classes")

  ggsave(pgg, filename = file.path(OutDir, fname))

  # ggplot ___________________________________
  quantiles <- c(0.025, 0.01, 0.0025)

  i = 1
  for (i in 1:length(quantiles)) {
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
      ggplot(Soup.VS.Cells.Av.Exp.gg, aes(x= Soup, y= Cells, label = gene,
                                          col= Outlier))  +
      geom_point() + theme(legend.position = "none") +
      xlab(paste(axl.pfx, "Soup", axl.sfx)) + ylab(paste(axl.pfx, "Cells", axl.sfx)) +
      ggtitle("Soup VS. Cells", subtitle = pr) +
      ggrepel::geom_text_repel(aes(label= ifelse(Outlier
                                                 , as.character(gene),'')))
    ggsave(pgg, filename = file.path(OutDir, fname))
  }


  # Per Gene ___________________________________
  PC.mRNA.in.Soup <- sum(CR.matrices$'soup')/sum(CR.matrices$'raw')
  PC.mRNA.in.Cells <- 100*sum(CR.matrices$'filt')/sum(CR.matrices$'raw')
  wbarplot(variable = PC.mRNA.in.Cells, col ="seagreen", plotname = kppd("PC.mRNA.in.Cells", SeqRun)
           , ylim = c(0,100), ylab = "% mRNA in cells"
           , sub = "% mRNA is more meaningful than % reads reported by CR")
  barplot_label(barplotted_variable = PC.mRNA.in.Cells
                , labels = Stringendo::percentage_formatter(PC.mRNA.in.Cells/100, digitz = 2)
                , TopOffset = 10)


  # Plot top gene's expression ___________________________________
  Soup.GEMs.top.Genes = 100*head(sort(CR.matrices$'soup.rel.RC', decreasing = T), n = 20)

  wbarplot(Soup.GEMs.top.Genes, plotname = kppd("Soup.GEMs.top.Genes", SeqRun)
           , ylab="% mRNA in the Soup"
           , sub = paste("Within the", SeqRun, "dataset")
           , tilted_text = T
           , ylim = c(0, max(Soup.GEMs.top.Genes)*1.5))
  barplot_label(barplotted_variable = Soup.GEMs.top.Genes
                , labels = Stringendo::percentage_formatter(Soup.GEMs.top.Genes/100, digitz = 2)
                , TopOffset = -.5, srt = 90, cex=.75)

  # Plot summarize expression ___________________________________
  soupProfile <- CR.matrices$'soup.total.RC'
  {
    soup.RP.sum   <- sum(soupProfile[grep('^RPL|^RPS', names(soupProfile))])
    soup.RPL.sum   <- sum(soupProfile[grep('^RPL', names(soupProfile))])
    soup.RPS.sum   <- sum(soupProfile[grep('^RPS', names(soupProfile))])
    soup.mito.sum <- sum(soupProfile[grep('^MT-', names(soupProfile))])
    soup.LINC.sum <- sum(soupProfile[grep('^LINC', names(soupProfile))])
    soup.AC.sum <- sum(soupProfile[grep('^AC', names(soupProfile))])
    soup.AL.sum <- sum(soupProfile[grep('^AL', names(soupProfile))])
    genes.non.Above <- soupProfile[CodeAndRoll2::grepv('^RPL|^RPS|^MT-|^LINC|^AC|^AL', names(soupProfile), invert = T)]
  }
  head(sort(genes.non.Above), n = 50)


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
  ccc <- c("#FF4E00","#778B04","#8ea604","#8ea604","#F5BB00","#F5BB00","#EC9F05",rep(x = "#BF3100", times = NrColumns2Show-6)) # ,"#"


  Soup.GEMs.top.Genes.summarized = 100 * soupProfile.summarized[1:NrColumns2Show] / CR.matrices$'soup.total.sum'
  maxx <- max(Soup.GEMs.top.Genes.summarized)
  wbarplot(Soup.GEMs.top.Genes.summarized, plotname = kppd("Soup.GEMs.top.Genes.summarized", SeqRun)
           , ylab="% mRNA in the Soup", ylim = c(0, maxx+3)
           , sub = paste("Within the", SeqRun, "dataset")
           , tilted_text = T, col = ccc)
  barplot_label(barplotted_variable = Soup.GEMs.top.Genes.summarized
                , srt = 45, labels = Stringendo::percentage_formatter(Soup.GEMs.top.Genes.summarized/100, digitz = 2)
                , TopOffset = -1.5)

  # Absolute.fraction ___________________________________
  Absolute.fraction.soupProfile.summarized <- Soup.GEMs.top.Genes.summarized * PC.mRNA.in.Soup

  maxx <- max(Absolute.fraction.soupProfile.summarized)
  wbarplot(Absolute.fraction.soupProfile.summarized, plotname = kppd("Absolute.fraction.soupProfile.summarized", SeqRun)
           , ylab="% of mRNA in cells", ylim = c(0, maxx*1.33)
           , sub = paste(Stringendo::percentage_formatter(PC.mRNA.in.Soup), "of mRNA counts are in the Soup, in the dataset ", SeqRun)
           , tilted_text = T, col = ccc)
  barplot_label(barplotted_variable = Absolute.fraction.soupProfile.summarized
                , srt = 45, labels = Stringendo::percentage_formatter(Absolute.fraction.soupProfile.summarized/100, digitz = 2)
                # formatC(Absolute.fraction.soupProfile.summarized, format="f", big.mark = " ", digits = 0)
                , TopOffset = -maxx*0.15)

  # ___________________________________
  Soup.GEMs.top.Genes.non.summarized <- 100* sort(genes.non.Above, decreasing = T)[1:20]/ CR.matrices$'soup.total.sum'
  maxx <- max(Soup.GEMs.top.Genes.non.summarized)
  wbarplot(Soup.GEMs.top.Genes.non.summarized, plotname = kppd("Soup.GEMs.top.Genes.non.summarized", SeqRun)
           , ylab="% mRNA in the Soup"
           , sub = paste("Within the", SeqRun, "dataset")
           , tilted_text = T, col = "#BF3100"
           , ylim = c(0, maxx*1.5))
  barplot_label(barplotted_variable = Soup.GEMs.top.Genes.non.summarized
                , labels = Stringendo::percentage_formatter(Soup.GEMs.top.Genes.non.summarized/100, digitz = 2)
                # , labels = paste0(round(1e6 * Soup.GEMs.top.Genes.non.summarized), " ppm")
                , TopOffset = -maxx*0.2, srt = 90, cex=.75)
  if (exists('OutDirBac'))  ww.assign_to_global("OutDir", OutDirBac, 1)
} # plotTheSoup




#' @title load10Xv3
#' @description Load 10X output folders.
#' @param dataDir PARAM_DESCRIPTION
#' @param cellIDs PARAM_DESCRIPTION, Default: NULL
#' @param channelName PARAM_DESCRIPTION, Default: NULL
#' @param readArgs PARAM_DESCRIPTION, Default: list()
#' @param includeFeatures PARAM_DESCRIPTION, Default: c("Gene Expression")
#' @param verbose PARAM_DESCRIPTION, Default: TRUE
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @seealso
#'  \code{\link[SoupX]{SoupChannel}}
#' @export
#' @importFrom SoupX SoupChannel
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

# _________________________________________________________________________________________________


#' @title parallel.computing.by.future
#' @description Run gc(), load multi-session computing and extend memory limits.
#'
#' @param workers_ cores
#' @param maxMemSize memory limit
#'
#' @export
parallel.computing.by.future <- function(workers_ = 6, maxMemSize = 4000 * 1024^2) { # Run gc(), load multi-session computing and extend memory limits.
  # https://satijalab.org/seurat/v3.0/future_vignette.html
  cat(
    "1. If you load futures before you finished using foreach loops,
    NormalizeData inside a foreach loop fails (Error writing to connection)
    -> I assume 'future' and 'doMC' are not compatible

    2. If you setup computing on e.g. six cores, it runs 6 instances of R with the entire memory space copied.
    If you run out of memory, the system starts using the SSD as memory, and it slows you down extremely extremely extremely.
    -> Therefore it is important to clean up the memory space before setting up multicore computation.

    Loaded: library(future), workers set to 6 (def),set Max mem size to 2GB (def)."   )

  gc(full = T)
  try(memory.biggest.objects(), silent = T)

  library(future)
  plan("multiprocess", workers = workers_)
  # plan("multisession", workers = workers_)
  # So to set Max mem size to 2GB, you would run :
  options(future.globals.maxSize = maxMemSize)
}



# _________________________________________________________________________________________________
#' getClusterNames
#'
#' @param obj Seurat object
#' @param ident ident
#' @export

getClusterNames <- function(obj = combined.obj, ident = GetClusteringRuns(obj)[2]) {
  print(GetClusteringRuns(obj))
  clz <- as.character(sort(deframe(unique(obj[[ident]]))))
  cat(dput(clz))
}




# _________________________________________________________________________________________________
#' RenameClustering
#'
#' @param namedVector named vector, where values = new, names(vec) = old
#' @param orig.ident meta.data colname original
#' @param suffix.new.ident How to name (suffix) the new identity. Default: "ManualNames"
#' @param new.ident meta.data colname new
#' @param obj Seurat object
#' @export

RenameClustering <- function(namedVector = ManualNames
                             , orig.ident =  "RNA_snn_res.0.3"
                             , suffix.new.ident = "ManualNames"
                             , new.ident = ppp(orig.ident, suffix.new.ident)
                             , obj = combined.obj) {
  NewX <- translate(vec = as.character(obj@meta.data[ ,orig.ident])
                    , oldvalues = names(namedVector)
                    , newvalues = namedVector)
  obj@meta.data[[new.ident]] <- NewX
  clUMAP(orig.ident)
  clUMAP(new.ident)
  return(obj)
}

# _________________________________________________________________________________________________
#' IntersectWithExpressed
#'
#' @param genes A vector of gene names.
#' @param obj Seurat object
#' @param genes.shown Number of genes printed (head).
#' @export

IntersectWithExpressed <- function(genes, obj=combined.obj, genes.shown = 10) { # Intersect a set of genes with genes in the Seurat object.
  print('IntersectWithExpressed()')
  # print(head(genes, n=15))
  diff = setdiff(genes, rownames(obj))
  Stringendo::iprint(length(diff),"genes (of",length(genes), ") are MISSING from the Seurat object:",head(diff, genes.shown))
  return(intersect(rownames(obj), genes))
}


# _________________________________________________________________________________________________
#' seu.RemoveMetadata
#'
#' @param obj Seurat object, Default: combined.obj
#' @param cols_remove columns to remove
#' @export

seu.RemoveMetadata <- function(obj = combined.obj
                               , cols_remove = grepv(colnames(obj@meta.data), pattern = "^integr|^cl.names", perl = T)
) {

  CNN <- colnames(obj@meta.data)
  iprint('cols_remove:', cols_remove); print('')
  (cols_keep <- setdiff(CNN, cols_remove))
  obj@meta.data <- obj@meta.data[, cols_keep]
  iprint('meta.data colnames kept:', colnames(obj@meta.data))

  return(obj)
}



# _________________________________________________________________________________________________
#' Percent.in.Trome
#' @description Gene expression as fraction of all UMI's
#' @param obj Seurat object
#' @param n.genes.barplot number of top genes shows
#' @param width.barplot barplot width
#' @return Seurat object
#' @export
#' @examples # combined.obj <- Percent.in.Trome()

Percent.in.Trome <- function(obj = combined.obj, n.genes.barplot = 25, width.barplot = round(n.genes.barplot/4)) {
  m.expr <- combined.obj@assays$RNA@counts
  total.Expr <- sort(rowSums(m.expr), decreasing = T)
  relative.total.Expr <- total.Expr / sum(total.Expr)
  print(head(iround(100*relative.total.Expr), n = n.genes.barplot))

  qhistogram(relative.total.Expr*100, logX = F, logY = T
             , plotname = "Gene expression as fraction of all UMI's"
             , subtitle = "Percentage in RNA-counts"
             , xlab = "Percent in Transcriptome (total per gene)"
             , ylab = "Number of genes"
             , xlab.angle = 45
  ) # + geom_hline(yintercept = 10)

  Highest.Expressed.Genes <- head(iround(100*relative.total.Expr), n = n.genes.barplot)
  qbarplot(Highest.Expressed.Genes, w = width.barplot
           , plotname = "Percentage of highest expressed genes"
           , subtitle = "Total, in RNA-counts"
           , xlab = ""
           , ylab = "Gene expression as percent of all UMI's"
           , xlab.angle = 45
  )
  print("!!!")
  print("TotalReadFraction is stored under combined.obj@misc$'TotalReadFraction'  ")
  print("!!!")
  combined.obj@misc$'TotalReadFraction' <- relative.total.Expr
  return(combined.obj)

}
