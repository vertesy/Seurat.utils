######################################################################
# plotting.dim.reduction.3D.R
######################################################################
# source ('~/GitHub/Seurat.utils/plotting.dim.reduction.3D.R')
# Source: self + https://github.com/Dragonmasterx87/Interactive-3D-Plotting-in-Seurat-3.0.0

# Requirements ------------------------
library(plotly)
library(MarkdownReportsDev)
# May also require
# try (source ('~/GitHub/TheCorvinas/R/CodeAndRoll.R'),silent= F) # generic utilities funtions
# require('MarkdownReportsDev') # require("devtools") # plotting related utilities functions # devtools::install_github(repo = "vertesy/MarkdownReportsDev")


# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
plot3D.umap <- function(obj=combined.obj, # Plot a 3D umap based on one of the metadata columns. Uses plotly. Based on github.com/Dragonmasterx87.
  category="v.project", AutoAnnotByCluster=c(FALSE, category, "integrated_snn_res.0.7")[3]) {
  stopifnot(category %in% colnames(obj@meta.data))
  stopifnot("UMAP_3" %in% colnames(obj@reductions$umap))
  plotting.data <- FetchData(object = obj, vars = c("UMAP_1", "UMAP_2", "UMAP_3", category))
  colnames(plotting.data)[4] = "category"
  plotting.data$label <- paste(rownames(plotting.data))   # Make a column of row name identities (these will be your cell/barcode names)

  if (AutoAnnotByCluster != FALSE) {
    # https://plot.ly/r/text-and-annotations/#3d-annotations
    stopifnot(AutoAnnotByCluster %in% colnames(obj@meta.data))

    plotting.data$'annot' <- FetchData(object = obj, vars = c(AutoAnnotByCluster))[,1]

    auto_annot <-
      plotting.data %>%
      group_by(annot)%>%
      summarise(showarrow=F
                , xanchor = "left"
                , xshift = 10
                , opacity = 0.7
                ,"x" = mean(UMAP_1)
                , "y" = mean(UMAP_2)
                , "z" = mean(UMAP_3)

      )
    names(auto_annot)[1]="text"
    ls.ann.auto = apply(auto_annot, 1, as.list)
  } else {ls.ann.auto <- NULL}

  plot_ly(data = plotting.data
          , x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3
          , color = ~category
          , colors = gg_color_hue(length(unique(plotting.data$'category')))
          , type = "scatter3d"
          , mode = "markers"
          , marker = list(size = 1)
          , text=~label
          # , hoverinfo="text"
  ) %>%  layout(scene = list(title=category, annotations=ls.ann.auto))

}


# ------------------------------------------------------------------------
BackupReduction <- function(obj = combined.obj, dim=2, reduction="umap") { # Backup UMAP to `obj@misc$reductions.backup` from `obj@reductions$umap`.
  dslot=p0(reduction,dim,"d")
  obj@misc$reductions.backup[[dslot]] <- obj@reductions[[reduction]]
  return(obj)
}
# Example
# obj <- BackupReduction(obj = obj, dim=2, reduction=umap"")

# ------------------------------------------------------------------------
SetupReductionsNtoKdimensions <- function(obj = combined.obj, nPCs = p$'n.PC', dimensions=3:2, reduction="umap") { # Calculate N-to-K dimensional umaps (default = 2:3); and back them up UMAP to `obj@misc$reductions.backup` from @reductions$umap
  red <- reduction
  for (d in dimensions) {
    iprint(d, "dimensional", red, "is calculated")
    obj <- if (reduction == "umap") {
      RunUMAP(obj, dims = 1:nPCs, n.components = d)
    } else if (reduction == "tsne") {
      RunTSNE(obj, dims = 1:nPCs, n.components = d)
    } else if (reduction == "pca") {
      RunPCA(obj, dims = 1:nPCs, n.components = d)
    }
    obj <- BackupReduction(obj = obj, dim=d, reduction=red)
  }
  return(obj)
}

# Example
# combined.obj <- SetupReductionsNtoKdimensions(obj = combined.obj, nPCs = p$'n.PC', dimensions=2:3, reduction="umap")
# qUMAP()

# ------------------------------------------------------------------------
RecallReduction <- function(obj = combined.obj, dim=2, reduction="umap") { # Set active UMAP to `obj@reductions$umap` from `obj@misc$reductions.backup`.
  dslot=p0(reduction,dim,"d")
  iprint(dim, "dimensional", reduction, "is set active. Source")
  obj@reductions[[reduction]] <- obj@misc$reductions.backup[[dslot]]
  return(obj)
}

# Example
# combined.obj <- RecallReduction(obj = combined.obj, dim=2, reduction="umap")
# qUMAP()
# combined.obj <- RecallReduction(obj = combined.obj, dim=3, reduction="umap")
# qUMAP()


# ------------------------------------------------------------------------
