######################################################################
# plotting.statistics.and.QC.R
######################################################################
# source('~/GitHub/Packages/Seurat.utils/plotting.statistics.and.QC.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Plotting.statistics.and.QC.R"))

# Source: self + web

# Requirements ------------------------
require(Seurat)
require(ggplot2)
# tools for tools::toTitleCase

# May also require
# try (source('/GitHub/Packages/CodeAndRoll/CodeAndRoll.R'),silent= F) # generic utilities funtions
# require('MarkdownReportsDev') # require("devtools") # plotting related utilities functions # devtools::install_github(repo = "vertesy/MarkdownReportsDev")

# PCA percent of variation associated with each PC ------------------------------------------------------------
seu.PC.var.explained <- function(obj =  combined.obj) { # Determine percent of variation associated with each PC.
  pct <- obj@reductions$pca@stdev / sum(obj@reductions$pca@stdev) * 100
  names(pct) =1:length(obj@reductions$pca@stdev)
  return(pct)
}

# plot percent of variation associated with each PC ------------------------------------------------------------
seu.plot.PC.var.explained <- function(obj =  combined.obj, use.MDrep = F) { # Plot the percent of variation associated with each PC.
  pct <- seu.PC.var.explained(obj)
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

    ggplot( aes(fill=(!!(as.name(fill.by))), x = (!!(as.name(group.by)))) ) +
    geom_hline( yintercept = hlines, lwd=1.5)  +
    geom_bar( position="fill" ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_text(aes(label=..count..), stat='count',position=position_fill(vjust=0.5)) +
    labs(title =pname_, x = "Clusters", y = "Fraction", caption = caption_) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
}
# CellFractionsBarplot2(obj = combined.obj, group.by = "integrated_snn_res.0.1", fill.by = "Phase", downsample = T)
# CellFractionsBarplot2(obj = combined.obj, group.by = "integrated_snn_res.0.1", fill.by = "Phase", downsample = F)

# BulkGEScatterPlot ------------------------------------------------------------------------
BulkGEScatterPlot <- function(obj = combined.obj # Plot bulk scatterplots to identify differential expressed genes across conditions
                              , clusters = "cl.names.KnownMarkers.0.2", TwoCategIdent = 'age', genes.from.bulk.DE = rownames(df.markers.per.AGE)) {

  (SplitIdents <- unique(obj[[TwoCategIdent]][,1]))
  stopifnot(length(SplitIdents) == 2)

  Idents(obj) <- clusters
  (IdentsUsed <- sort.natural(as.character(unique(Idents(obj)))))
  NrPlots <- length(IdentsUsed)
  p.clAv <- p.clAv.AutoLabel <- genes.to.label <- list.fromNames(IdentsUsed)

  i=1
  for (i in 1:NrPlots) {
    print(IdentsUsed[i])
    ClX <- subset(obj, idents = IdentsUsed[i])
    Idents(ClX) <- TwoCategIdent
    avg.ClX.cells <- log2(AverageExpression(ClX, verbose = FALSE)$RNA + 1)
    avg.ClX.cells$gene <- rownames(avg.ClX.cells)

    # plot ----
    p.clAv[[i]] <- p.clAv.AutoLabel[[i]] <-
      ggplot(avg.ClX.cells, aes(x = !!as.name(SplitIdents[1]), y = !!as.name(SplitIdents[2]) )) +
      geom_point(data = avg.ClX.cells, color=rgb(0,.5,0,0.25), size=1) +
      FontSize(x.title = 8, x.text = 8, y.title = 8, y.text = 8)+
      geom_abline(slope = 1, intercept = 0, color='grey') +
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

  PlotIter <- iterBy.over(1:NrPlots, by = 4)
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


# plot.Gene.Cor.Heatmap ------------------------------------------------------------------------

plot.Gene.Cor.Heatmap <- function(genes = WU.2017.139.IEGsf
                                  , assay.use = "RNA", slot.use = "data"
                                  , min.g.cor =  0.3
                                  , obj = combined.obj, ...) {
  expr.mat <- GetAssayData(slot = slot.use, assay = assay.use, object = obj)

  slotname_cor.mat <- kpp('cor', slot.use, assay.use)
  cor.mat <- obj@misc[[slotname_cor.mat]]
  stopif(is_null(cor.mat), message = kollapse(slotname_cor.mat, " not found in @misc."))

  intersect(genes, rownames(cor.mat))
  genes.found <- check.genes(genes)
  if (l(genes.found) > 200) iprint("Too many genes found in data, cor will be slow: ", l(genes.found))
  ls.cor <- sparse.cor(t(expr.mat[genes.found,]))

  # Filter
  corrz <- ls.cor$cor
  diag(corrz) <- NaN
  corgene.names <- union(
    which_names(rowMax(corrz) >= min.g.cor),
    which_names(rowMin(corrz) <= -min.g.cor)
  )

  pname = p0("Pearson correlations of ", substitute(genes),"\n min.cor:", min.g.cor, " | ",  assay.use ,'.', slot.use )
  o.heatmap <- pheatmap(corrz[corgene.names,corgene.names],main = pname,..., annota)
  wplot_save_pheatmap(o.heatmap, filename = make.names(pname))
}

# Calc.Cor.Seurat ------------------------------------------------
Calc.Cor.Seurat <- function(assay.use = "RNA", slot.use = "data"
                            , quantileX = 0.99, max.cells =  40000, seed = p$"seed"
                            , digits = 2, obj = combined.obj) {
  expr.mat <- GetAssayData(slot = slot.use, assay = assay.use, object = obj)
  if (ncol(expr.mat) > max.cells) {
    set.seed(seed = seed)
    cells.use <- sample(x = colnames(expr.mat), size = max.cells)
  }

  qname = p0("q", quantileX * 100)
  slot_name = kpp("expr", qname)

  if (is.null(obj@misc[[slot_name]])) iprint("Call: combined.obj <- Calcq90Expression(combined.obj, quantileX =",quantileX," first )")
  genes.HE = which_names(obj@misc[[slot_name]] > 0)
  iprint("Pearson correlation is calculated for", l(genes.HE), "HE genes with expr.q99 > 0.")
  tic(); ls.cor <- sparse.cor(smat = t(expr.mat[genes.HE, cells.use])); toc()
  ls.cor <- lapply(ls.cor, round, digits = 2)

  obj@misc[[kpp('cor', slot.use, assay.use)]] <- ls.cor$'cor'
  obj@misc[[kpp('cov', slot.use, assay.use)]] <- ls.cor$'cov'
  iprint("Stored under obj@misc$", kpp('cor', slot.use, assay.use), "or cov... ." )
  return(obj)
}
# combined.obj <- Calcq90Expression(combined.obj, quantileX = 0.99, max.cells =  400000, set.all.genes = F)
# combined.obj <- Calc.Cor.Seurat(assay.use = "RNA", slot.use = "data", digits = 2, obj = combined.obj, quantile = 0.99, max.cells = 40000)

