######################################################################
# Cluster.Auto-naming.DE.R
######################################################################
# source ('~/GitHub/Seurat.utils/Cluster.Auto-naming.DE.R')

# ------------------------------------------------------------------------
require(princurve)

# ------------------------------------------------------------------------------------
StoreTop25Markers <- function(obj = combined.obj # Save the top 25 makers based on `avg_logFC` output table of `FindAllMarkers()` (df_markers) under `@misc$df.markers$res...`. By default, it rounds up insignificant digits up to 3.
                              , df_markers = df.markers, res = 0.5) {
  top25.markers <-
    df_markers %>%
    group_by(cluster) %>%
    top_n(n = 25, wt = avg_logFC) %>%
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


# GetTopMarkers ------------------------------------------------------------------------------------
GetTopMarkers <- function(df = df.markers # Get the vector of N most diff. exp. genes.
                          , n = p$'n.markers', order.by = c("avg_logFC", "p_val_adj")[1]) {
  'Works on active Idents() -> thus we call cluster'

  TopMarkers <- df %>%
    arrange(desc(!!as.name(order.by))) %>%
    group_by(cluster) %>%
    dplyr::slice(1:n) %>%
    # group_by(cluster) %>% # OLD WRONG SOLUTION
    # top_n(n = n, wt = (!!as.name(order.by))) %>% # OLD WRONG SOLUTION
    dplyr::select(gene) %>%
    col2named.vec.tbl()

  return(TopMarkers)
}
# GetTopMarkers(df = df.markers, n=3 )

# ------------------------------------------------------------------------------------
AutoLabelTop.logFC <- function(obj = combined.obj # Create a new "named identity" column in the metadata of a Seurat object, with `Ident` set to a clustering output matching the `res` parameter of the function. It requires the output table of `FindAllMarkers()`. If you used `StoreAllMarkers()` is stored under `@misc$df.markers$res...`, which location is assumed by default.
                               , res = 0.5, plot.top.genes = T
                               , df_markers = combined.obj@misc$"df.markers"[[p0("res.",res)]] ) {
  top.markers <- GetTopMarkers(df = df_markers, n=1)
  # df_markers %>%
  # group_by(cluster) %>%
  # top_n(n = 1, wt = avg_logFC) %>%
  # dplyr::select(gene) %>%
  # col2named.vec.tbl()
  stopifnot(length(unique(Idents(object = obj))) == length(top.markers))

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
AutoLabel.KnownMarkers <- function(obj = combined.obj # Create a new "named identity" column in the metadata of a Seurat object, with `Ident` set to a clustering output matching the `res` parameter of the function. It requires the output table of `FindAllMarkers()`. If you used `StoreAllMarkers()` is stored under `@misc$df.markers$res...`, which location is assumed by default.
                                   , KnownMarkers=c("TOP2A", "EOMES", "SLA", "HOPX", "S100B", "DLX6-AS1", "POU5F1","SALL4","DDIT4", "PDK1", "SATB2", "FEZF2")
                                   , res = 0.5 , df_markers = combined.obj@misc$"df.markers"[[p0("res.",res)]] ) {
  matching.clusters <-
    df_markers %>%
    arrange(desc(avg_logFC)) %>%
    filter(gene %in%  KnownMarkers) %>%
    dplyr::select(avg_logFC, p_val_adj, cluster, gene ) %>%
    group_by(gene) %>%
    top_n(n=1, wt=avg_logFC) %>% # Select the top cluster for each gene
    arrange(cluster)

  print(matching.clusters)

  unique.matches <-
    matching.clusters %>%
    group_by(cluster) %>% # Select rows with unique values based on column "cluster"
    distinct(cluster,.keep_all = T)  %>%
    dplyr::select(gene)

  print("Best matches:")
  print(unique.matches)

  top.markers.df <-
    df_markers %>%
    group_by(cluster) %>%
    top_n(n = 1, wt = (avg_logFC)) %>% # Find the top expressed gene
    dplyr::select(gene)

  top.markers <-
    top.markers.df %>%
    col2named.vec.tbl()

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
                                  , ident = "cl.names.top.gene.res.0.5", reduct ="umap",title = ident) {
  DimPlot(object = obj, reduction = reduct, group.by=ident, label=T, repel=T) + NoLegend() + ggtitle(title)
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
  NewLabel <- as.character(0:(length(MedianClusterCoordinate)-1))
  NewMeta <- translate(vec = obj[[res]], oldvalues = OldLabel, newvalues = NewLabel)
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

  dim_name <- ppu(toupper(reduction),dimension)
  coord.umap <- FetchData(object = obj, vars = dim_name)
  fit <- principal_curve(x = as.matrix(coord.umap))
  if (plotit) {
    plot(fit, xlim=range(coord.umap[,1]), ylim=range(coord.umap[,2])
         , main = "principal_curve")
    # points(fit)
    points(coord.umap, pch=18, cex=.25)
    whiskers(coord.umap, fit$s, lwd=.1)
    wplot_save_this(plotname = "principal_curve")
  }

  ls.perCl <- split(swap * fit$lambda, f = obj[[res]])
  MedianClusterCoordinate <- unlapply(ls.perCl, median)
  OldLabel <- names(sort(MedianClusterCoordinate))
  NewLabel <- as.character(0:(length(MedianClusterCoordinate)-1))
  NewMeta <- translate(vec = obj[[res]], oldvalues = OldLabel, newvalues = NewLabel)
  NewMetaCol <- kpp(res,"prin.curve")
  iprint("NewMetaCol:",NewMetaCol)
  obj[[NewMetaCol]] <- NewMeta
  return(obj)
}
# DimPlot.ClusterNames(ident = "integrated_snn_res.0.5")
# combined.obj <- AutoNumber.by.PrinCurve(obj = combined.obj, dimension=1:2, reduction="umap", plotit=T, swap= -1, res = "integrated_snn_res.0.5" )
# DimPlot.ClusterNames(ident = "integrated_snn_res.0.5.prin.curve")
