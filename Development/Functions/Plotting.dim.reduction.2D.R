######################################################################
# plotting.dim.reduction.2D.R
######################################################################
# source('~/GitHub/Packages/Seurat.utils/Functions/plotting.dim.reduction.2D.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Plotting.dim.reduction.2D.R"))
# Source: self + web

# Requirements ------------------------
library(plotly)
# try(source("~/GitHub/Packages/ggExpressDev/ggExpress.functions.R"), silent = T)
try(source("https://raw.githubusercontent.com/vertesy/ggExpressDev/main/ggExpress.functions.R"), silent = T)

# May also require
# try (source('/GitHub/Packages/CodeAndRoll/CodeAndRoll.R'),silent= F) # generic utilities funtions
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

