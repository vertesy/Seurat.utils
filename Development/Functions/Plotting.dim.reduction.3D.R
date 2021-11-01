######################################################################
# plotting.dim.reduction.3D.R
######################################################################
# source('~/GitHub/Packages/Seurat.utils/Functions/plotting.dim.reduction.3D.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Plotting.dim.reduction.3D.R"))
# Source: self + https://github.com/Dragonmasterx87/Interactive-3D-Plotting-in-Seurat-3.0.0

# Requirements ------------------------
try(library(plotly), silent = T)
try(library(MarkdownReportsDev), silent = T)
try(library(htmlwidgets), silent = T)

# May also require
# try (source('~/GitHub/Packages/CodeAndRoll/CodeAndRoll.R'),silent= T) # generic utilities funtions
# require('MarkdownReportsDev') # require("devtools") # plotting related utilities functions # devtools::install_github(repo = "vertesy/MarkdownReportsDev")


# ------------------------------------------------------------------------
ww.check.if.3D.reduction.exist <- function(obj = obj) { # ww.check.if.3D.reduction.exist in backup slot
  if( !("UMAP_3" %in% colnames(obj@reductions$'umap'))) {
    stopif2( is.null(combined.obj@misc$reductions.backup$'umap3d')
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
