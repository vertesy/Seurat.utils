######################################################################
# Cluster.Auto-naming.DE.R
######################################################################
# source ('~/GitHub/Seurat.utils/Cluster.Auto-naming.DE.R')

# ------------------------------------------------------------------------


# ------------------------------------------------------------------------------------
StoreAllMarkers <- function(obj =  combined.obj # Save the output table of `FindAllMarkers()` (df_markers) under `@misc$df.markers$res...`. By default, it rounds up insignificant digits up to 3.
                            , df_markers = df.markers, res = 0.5, digit=c(0,3)[2]) {
  if (digit) df_markers[,1:5] <- signif(df_markers[,1:5], digits = digit)
  combined.obj@misc$'df.markers'[[ppp('res',res)]] <- df_markers
}
# combined.obj <- StoreAllMarkers(df_markers = df.markers, res = 0.5)

# ------------------------------------------------------------------------------------
AutoLabelTop.logFC <- function(obj =  combined.obj # Create a new "named identity" column in the metadata of a Seurat object, with `Ident` set to a clustering output matching the `res` parameter of the function. It requires the output table of `FindAllMarkers()`. If you used `StoreAllMarkers()` is stored under `@misc$df.markers$res...`, which location is assumed by default.  
                               , res = 0.5 , df_markers = combined.obj@misc$"df.markers"[[p0("res.",res)]] ) {
  top.markers <- 
    df_markers %>%
    group_by(cluster) %>%
    top_n(n = 1, wt = avg_logFC) %>%
    select(gene) %>%
    col2named.vec.tbl()
  stopifnot(length(unique(Idents(object = combined.obj))) == length(top.markers))
  
  (top.markers.ID <- ppp(names(top.markers), top.markers))
  names(top.markers.ID) <- names(top.markers)
  named.ident <- top.markers.ID[Idents(object = combined.obj)]
  
  namedIDslot <- ppp('cl.names.top.gene.res', res )
  combined.obj[[namedIDslot]] = named.ident
  return(combined.obj)
}
# combined.obj <- AutoLabelTop.logFC(); combined.obj$"cl.names.top.gene.res.0.5"


# ------------------------------------------------------------------------------------
AutoLabel.KnownMarkers <- function(obj =  combined.obj # Create a new "named identity" column in the metadata of a Seurat object, with `Ident` set to a clustering output matching the `res` parameter of the function. It requires the output table of `FindAllMarkers()`. If you used `StoreAllMarkers()` is stored under `@misc$df.markers$res...`, which location is assumed by default.  
                                   , KnownMarkers=c("TOP2A", "EOMES", "SLA", "HOPX", "S100B", "DLX6-AS1", "POU5F1","SALL4","DDIT4", "PDK1", "SATB2", "FEZF2") 
                                   , res = 0.5 , df_markers = combined.obj@misc$"df.markers"[[p0("res.",res)]] ) {
  matching.clusters <- 
    df_markers %>%
    arrange(desc(avg_logFC)) %>%
    filter(gene %in%  KnownMarkers) %>%
    select(avg_logFC, p_val_adj, cluster, gene ) %>% 
    group_by(gene) %>% 
    top_n(n=1, wt=avg_logFC) %>% # Select the top cluster for each gene
    arrange(cluster) 
  
  print(matching.clusters)
  
  unique.matches <- 
    matching.clusters %>%
    group_by(cluster) %>% # Select rows with unique values based on column "cluster"
    distinct(cluster,.keep_all = T)  %>%
    select(gene)
  
  print("Best matches:")
  print(unique.matches)
  
  top.markers.df <- 
    df_markers %>%
    group_by(cluster) %>%
    top_n(n = 1, wt = (avg_logFC)) %>% # Find the top expressed gene
    select(gene)
  
  missing.annotations <- 
    top.markers.df %>% 
    filter(!cluster %in%  unique.matches$cluster) # filter for clusters that do not have a unique label already
  
  named.annotations <- 
    rbind(unique.matches, missing.annotations) %>%  # merge the 2 df's
    arrange(cluster) %>%
    col2named.vec.tbl() # requires github.com/vertesy/CodeAndRoll

  (top.markers.ID <- ppp(names(named.annotations), named.annotations))
  names(top.markers.ID) <- names(top.markers)
  named.ident <- top.markers.ID[Idents(object = combined.obj)]
  
  namedIDslot <- ppp('cl.names.KnownMarkers', res )
  combined.obj[[namedIDslot]] = named.ident
  return(combined.obj)
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


