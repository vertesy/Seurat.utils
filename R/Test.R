

plot.qUMAPs.in.a.folder <- function(genes, obj = combined.obj, foldername = NULL, intersectionAssay = 'RNA', plot.reduction = 'umap', ...) {

  ParentDir = OutDir
  if (is.null(foldername)) foldername = substitute(genes)
  create_set_SubDir( paste0(foldername,'-', plot.reduction),'/')
  list.of.genes.found = check.genes(list.of.genes = genes, obj = obj, assay.slot = intersectionAssay, makeuppercase = F)


  for (g in list.of.genes.found) {
    qUMAP(g, reduction = plot.reduction, ...)
  }
  MarkdownReports::create_set_OutDir(... = ParentDir)
}



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
                                , solo = MarkdownHelpers::FALSE.unless('b.solo.plot')
                                , raster = MarkdownHelpers::FALSE.unless('b.raster')
                                , raster.dpi = c(512, 512)/4
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

    if (solo) {


    } else {
      plot.list = Seurat::FeaturePlot(object = obj, features = genes, reduction = plot.reduction, combine = F
                                      , ncol = nr.Col, cols = colors
                                      , min.cutoff = gene.min.exp, max.cutoff = gene.max.exp
                                      , raster = TRUE, raster.dpi = raster.dpi
                                      , pt.size = cex, ...)

      for (i in 1:length(plot.list)) {
        plot.list[[i]] <- plot.list[[i]] + NoLegend() + NoAxes()
        if (aspect.ratio) plot.list[[i]] <- plot.list[[i]] + ggplot2::coord_fixed(ratio = aspect.ratio)
      }

      pltGrid <- cowplot::plot_grid(plotlist = plot.list, ncol = nr.Col, nrow = nr.Row )
      ggsave(filename = plotname, width = w, height = h, bg = background_col, plot = pltGrid)
    }
  }


  if (subdir) MarkdownReports::create_set_OutDir(... = ParentDir)
  if (saveGeneList) {
    if (is.null(obj@misc$gene.lists)) obj@misc$gene.lists <- list()
    obj@misc$gene.lists[[substitute(list.of.genes)]] <- list.of.genes.found
    print("Genes saved under: obj@misc$gene.lists")
    return(obj)
  }
  tictoc::toc()
};

