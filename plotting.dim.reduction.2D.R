######################################################################
# plotting.dim.reduction.2D.R
######################################################################
# source ('~/GitHub/Seurat.utils/plotting.dim.reduction.2D.R')
# Source: self + web

# Requirements ------------------------
library(plotly)
# May also require
# try (source ('~/GitHub/CodeAndRoll/CodeAndRoll.R'),silent= F) # generic utilities funtions
# require('MarkdownReportsDev') # require("devtools") # plotting related utilities functions # devtools::install_github(repo = "vertesy/MarkdownReportsDev")


# Quick umaps# ------------------------------------------------------------------------
qUMAP <- function(f= 'TOP2A', obj =  combined.obj, splitby = NULL, qlow = "q10", qhigh = "q90") { # The quickest way to a draw a UMAP
  FeaturePlot(combined.obj, reduction = 'umap'
              , min.cutoff = qlow, max.cutoff = qhigh
              , split.by = splitby
              , features = f)
}
# qUMAP(  )

# ------------------------------------------------------------------------
gg_color_hue <- function(n) { # reproduce the ggplot2 default color palette
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
# https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette

# ------------------------------------------------------------------------
save2umaps.A4 <- function(plot_list, pname = F) { # Save 2 umaps on an A4 page.
  if (pname ==F) pname = substitute(plot_list)
  p1 = plot_grid(plotlist = plot_list, nrow = 2, ncol = 1, labels = LETTERS[1:length(plot_list)]  )
  save_plot(plot = p1, filename = extPNG(pname), base_height = hA4, base_width = wA4)
}

# ------------------------------------------------------------------------
save4umaps.A4 <- function(plot_list, pname = F) { # Save 4 umaps on an A4 page.
  if (pname==F) pname = substitute(plot_list)
  p1 = plot_grid(plotlist = plot_list, nrow = 2, ncol = 2, labels = LETTERS[1:length(plot_list)]  )
  save_plot(plot = p1, filename = extPNG(pname), base_height = wA4, base_width = hA4)
}

# ------------------------------------------------------------------------
umapNamedClusters <- function(obj = combined.obj, metaD.colname = metaD.colname.labeled, ext = "png") { # Plot and save umap based on a metadata column.
  fname = ppp("Named.clusters", metaD.colname, ext)
  p.named =
    DimPlot(obj, reduction = "umap", group.by = metaD.colname, label = T) +
    NoLegend() +
    ggtitle(metaD.colname)
  save_plot(p.named, filename = fname); p.named
}
# umapNamedClusters(obj = combined.obj, metaD.colname = metaD.colname.labeled)


# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------

# umapHiLightSel highlight a set of cells based on clusterIDs provided---------------
umapHiLightSel <- function(obj = combined.obj, # Highlight a set of cells based on clusterIDs provided.
                          COI =  c("0", "2", "4", "5",  "11"), res.cl = 'integrated_snn_res.0.3') {
  cellsSel = getCellIDs.from.meta(obj, values = COI, ColName.meta = res.cl)
  DimPlot(combined.obj, reduction = "umap", group.by = res.cl,
          label = T, cells.highlight = cellsSel)
  ggsave(filename = extPNG(kollapse("cells",COI, collapseby = '.')))
}



# Save multiple FeaturePlot from a list of genes on A4 jpeg ------------------------
multiFeaturePlot.A4 <- function(obj = combined.obj # Save multiple FeaturePlots, as jpeg, on A4 for each gene, which are stored as a list of gene names.
                                , list.of.genes, plot.reduction='umap', intersectionAssay = c('RNA', 'integrated')[1]
                                , colors=c("grey", "red"), nr.Col=2, nr.Row =4, cex = round(0.1/(nr.Col*nr.Row), digits = 2)
                                , gene.min.exp = 'q01', gene.max.exp = 'q99', subdir =T
                                , jpeg.res = 225, jpeg.q = 90) {
  tictoc::tic()
  ParentDir = OutDir
  if (subdir) create_set_SubDir(... = p0(substitute(list.of.genes),'.', plot.reduction),'/')

  list.of.genes = check.genes(list.of.genes = list.of.genes, obj = obj, assay.slot = intersectionAssay)
  lsG = iterBy.over(1:length(list.of.genes), by=nr.Row*nr.Col)
  for (i in 1:length(lsG)) {
    genes = list.of.genes[lsG[[i]]]
    iprint(i,genes )
    plotname = kpp(c(plot.reduction,i, genes, 'jpg' ))

    plot.list = FeaturePlot(object = obj, features =genes, reduction = plot.reduction, combine = F
                            , ncol = nr.Col, cols = colors
                            , min.cutoff = gene.min.exp, max.cutoff = gene.max.exp
                            , pt.size = cex)

    for(i in 1:length(plot.list)) {
      plot.list[[i]] <- plot.list[[i]] + NoLegend() + NoAxes()
    }

    ggsave(filename = plotname, width = wA4, height = hA4,
           plot = cowplot::plot_grid(plotlist = plot.list, ncol = nr.Col, nrow = nr.Row)
    )
  }

  if (subdir) create_set_OutDir(... = ParentDir)
  tictoc::toc()
};








# Save multiple FeatureHeatmaps from a list of genes on A4 jpeg -----------------------
# code for quantile: https://github.com/satijalab/seurat/blob/master/R/plotting_internal.R

multiFeatureHeatmap.A4 <- function(obj = combined.obj # Save multiple FeatureHeatmaps from a list of genes on A4 jpeg
                                   , list.of.genes, gene.per.page=5
                                   , group.cells.by= "batch", plot.reduction='umap'
                                   , cex = iround(3/gene.per.page), sep_scale = F
                                   , gene.min.exp = 'q5', gene.max.exp = 'q95'
                                   , jpeg.res = 225, jpeg.q = 90) {

  tictoc::tic()
  list.of.genes = check.genes(list.of.genes, obj = obj)

  lsG = iterBy.over(1:length(list.of.genes), by=gene.per.page)
  for (i in 1:length(lsG)) { print(i )
    genes = list.of.genes[lsG[[i]]]
    plotname = kpp(c("FeatureHeatmap",plot.reduction,i, genes, 'jpg' ))
    print(plotname)
    jjpegA4(plotname, r = jpeg.res, q = jpeg.q)
    try(
      FeatureHeatmap(obj, features.plot =genes , group.by = group.cells.by
                     , reduction.use = plot.reduction, do.return = F
                     , sep.scale = sep_scale, min.exp = gene.min.exp, max.exp = gene.max.exp
                     , pt.size = cex, key.position = "top")
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
                                      width = hA4, heigth = 1.75*wA4, filetype = "pdf") {

  p1 <- DimPlot(object = obj, reduction.use = "tsne", no.axes = no_axes, cells.use = cells_use
                , no.legend = no_legend, do.return = do_return, do.label = do_label, label.size = label_size
                , group.by = grouping, vector.friendly = vector_friendly, pt.size = pt_size) +
    ggtitle("tSNE") + theme(plot.title = element_text(hjust = 0.5))

  p2 <- DimPlot(object = obj, reduction.use = "umap", no.axes = no_axes, cells.use = cells_use
                , no.legend = T, do.return = do_return, do.label = do_label, label.size = label_size
                , group.by = grouping, vector.friendly = vector_friendly, pt.size = pt_size) +
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
