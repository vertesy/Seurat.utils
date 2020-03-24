######################################################################
# plotting.dim.reduction.3D.R
######################################################################
# source ('~/GitHub/Seurat.utils/plotting.dim.reduction.3D.R')
# Source: self + https://github.com/Dragonmasterx87/Interactive-3D-Plotting-in-Seurat-3.0.0

# Requirements ------------------------
library(plotly)
library(MarkdownReportsDev)
library(htmlwidgets)

# May also require
# try (source ('~/GitHub/CodeAndRoll/CodeAndRoll.R'),silent= F) # generic utilities funtions
# require('MarkdownReportsDev') # require("devtools") # plotting related utilities functions # devtools::install_github(repo = "vertesy/MarkdownReportsDev")


# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
plot3D.umap.gene <- function(obj=combined.obj # Plot a 3D umap with gene expression. Uses plotly. Based on github.com/Dragonmasterx87.
                             , gene="TOP2A", quantileCutoff = .99, alpha = .5, dotsize=1.25, def.assay = c("integrated", "RNA")[2]
                             , AutoAnnotBy=c(FALSE, category="v.project", "integrated_snn_res.0.7")[3]) {
  stopifnot(category %in% colnames(obj@meta.data))
  stopifnot("UMAP_3" %in% colnames(obj@reductions$umap))
  stopifnot(gene %in% rownames(obj))

  DefaultAssay(object = obj) <- def.assay

  plotting.data <- FetchData(object = obj, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "Expression"=gene), slot = 'data')
  # Cutoff <- quantile(plotting.data[,gene], probs = quantileCutoff)
  # plotting.data$'Expression' <- ifelse(test = plotting.data[,gene] < Cutoff, yes = plotting.data[,gene], no = Cutoff)
  plotting.data$'Expression' <- clip.outliers(plotting.data[,gene], probs = c(1-quantileCutoff, quantileCutoff))
  plotting.data$'label' <- paste(rownames(plotting.data)," - ", plotting.data[,gene], sep="")

  ls.ann.auto <- if (AutoAnnotBy != FALSE) {
    Annotate4Plotly3D(obj. = obj, plotting.data. = plotting.data, AnnotCateg = AutoAnnotBy)
  } else { NULL }

  plt <- plot_ly(data = plotting.data
                 , x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3
                 , type = "scatter3d"
                 , mode = "markers"
                 , marker = list(size = dotsize)
                 , text=~label
                 , color = ~Expression
                 , opacity = alpha
                 , colors = c('darkgrey', 'red')
                 #, hoverinfo="text"
  ) %>% layout(title=gene, scene = list(annotations=ls.ann.auto))
  SavePlotlyAsHtml(plt, category. = gene)
  return(plt)
}
# plot3D.umap.gene(obj = combined.obj, gene = "DDIT4", quantileCutoff = .95)

# ------------------------------------------------------------------------
plot3D.umap <- function(obj=combined.obj, # Plot a 3D umap based on one of the metadata columns. Uses plotly. Based on github.com/Dragonmasterx87.
  category="v.project", AutoAnnotBy=c(FALSE, category, "integrated_snn_res.0.7")[3]) {
  stopifnot(category %in% colnames(obj@meta.data))
  stopifnot("UMAP_3" %in% colnames(obj@reductions$umap))
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
          , marker = list(size = 1)
          , text=~label
          , color = ~category
          , colors = gg_color_hue(length(unique(plotting.data$'category')))
          # , hoverinfo="text"
  ) %>% layout(title=category, scene = list(annotations=ls.ann.auto))
  SavePlotlyAsHtml(plt, category. = category)
  return(plt)
}
# plot3D.umap(combined.obj, category = "Phase")

# ------------------------------------------------------------------------
SavePlotlyAsHtml <- function(plotly_obj, category.=category) { # Save Plotly 3D scatterplot as an html file.
  OutputDir <- if(exists("OutDir")) OutDir else getwd()
  fname <- kpp(OutputDir,"/umap.3D",category.,idate(),"html"); iprint("Plot saved as:",fname)
  htmlwidgets::saveWidget(plotly_obj, file = fname, selfcontained = TRUE, title = category.)
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
Annotate4Plotly3D <- function(obj. = combined.obj # Create annotation labels for 3D plots. Source https://plot.ly/r/text-and-annotations/#3d-annotations
                              , plotting.data. = plotting.data
                              , AnnotCateg = AutoAnnotBy) {
  stopifnot(AnnotCateg %in% colnames(obj@meta.data))

  plotting.data.$'annot' <- FetchData(object = obj, vars = c(AnnotCateg))[,1]
  auto_annot <-
    plotting.data. %>%
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
  return(ls.ann.auto)
}

# ------------------------------------------------------------------------
Plot3D.ListOfGenes <- function(obj = combined.obj # Plot and save list of 3D UMAP ot tSNE plots using plotly.
                               , annotate.by = "integrated_snn_res.0.7", opacity = 0.5, default.assay = c("integrated", "RNA")[2]
                               , ListOfGenes=c( "BCL11B" , "FEZF2", "EOMES", "DLX6-AS1", "HOPX", "DDIT4") ) {
  obj. <- obj; rm("obj")
  DefaultAssay(object = obj.) <- default.assay

  MissingGenes <- setdiff(ListOfGenes, rownames(obj.))
  if ( length(MissingGenes)) iprint("These genes are not found, and omitted:", MissingGenes, ". Try to change default assay.")
  ListOfGenes <- intersect(ListOfGenes, rownames(obj.))

  for (i in 1:length(ListOfGenes)) {
    g <- ListOfGenes[i]; print(g)
    plot3D.umap.gene(obj = obj., gene = g, AutoAnnotBy = annotate.by, alpha = opacity, def.assay = default.assay)
  }
  try(oo())
}
# gois <- c(  "PGK1", "CTIP2" = "BCL11B" , "FEZF2", "EOMES", "DLX6-AS1", "HOPX", "DDIT4","TOP2A", "PTGDS", "EDNRB", "EGFR", "SCGN", "NR2F2", "EMX2", "GAD2", "DLX2", "SATB2")
# Plot3D.ListOfGenes(obj = combined.obj, ListOfGenes = gois)


# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
