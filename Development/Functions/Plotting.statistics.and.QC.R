######################################################################
# plotting.statistics.and.QC.R
######################################################################
# source('~/GitHub/Packages/Seurat.utils/Functions/Plotting.statistics.and.QC.R')
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
scCalcPCAVarExplained <- function(obj =  combined.obj) { # Determine percent of variation associated with each PC.
  pct <- obj@reductions$pca@stdev / sum(obj@reductions$pca@stdev) * 100
  names(pct) =1:length(obj@reductions$pca@stdev)
  return(pct)
}

# plot percent of variation associated with each PC ------------------------------------------------------------
scPlotPCAvarExplained <- function(obj =  combined.obj, use.MDrep = F) { # Plot the percent of variation associated with each PC.
  pct <- scCalcPCAVarExplained(obj)
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

    ggplot( aes(fill = (!!(as.name(fill.by))),  x = (!!(as.name(group.by)))) ) +
    geom_hline( yintercept = hlines, lwd=1.5)  +
    geom_bar( position = "fill" ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_text(aes(label = ..count..), stat='count',position = position_fill(vjust = 0.5)) +
    labs(title = pname_,  x = "Clusters", y = "Fraction", caption = caption_) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
}
# CellFractionsBarplot2(obj = combined.obj, group.by = "integrated_snn_res.0.1", fill.by = "Phase", downsample = T)
# CellFractionsBarplot2(obj = combined.obj, group.by = "integrated_snn_res.0.1", fill.by = "Phase", downsample = F)



#  ------------------------------------------------
barplot.cells.per.cluster <- function(obj = combined.obj, ident =  "cl.names.KnownMarkers.0.5", sort = F) {
  cell.per.cluster <- (table(obj[[ident]][,1]))
  if (sort) cell.per.cluster <- sort(cell.per.cluster)
  qbarplot(cell.per.cluster, subtitle = ident, suffix = ident
           , col = rainbow(l(cell.per.cluster))
           , xlab.angle = 45
           # , col = getClusterColors(ident = ident, show = T)
           , palette_use = NULL, )
}
# barplot.cells.per.cluster()
# barplot.cells.per.cluster(sort=T)



# BulkGEScatterPlot ------------------------------------------------------------------------
BulkGEScatterPlot <- function(obj = combined.obj # Plot bulk scatterplots to identify differential expressed genes across conditions
                              , clusters = "cl.names.KnownMarkers.0.2", TwoCategIdent = 'age', genes.from.bulk.DE = rownames(df.markers.per.AGE)) {

  (SplitIdents <- unique(obj[[TwoCategIdent]][,1]))
  stopifnot(length(SplitIdents) == 2)

  Idents(obj) <- clusters
  IdentsUsed <- gtools::mixedsort(as.character(unique(Idents(obj))))
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
    p.clAv.AutoLabel[[i]] <- LabelPoints(plot = p.clAv[[i]], points = genes.to.label[[i]], xnudge = 0, ynudge = 0, repel = TRUE, size=2);
    p.clAv.AutoLabel[[i]]

    "Pre-identified genes"
    p.clAv[[i]] <- LabelPoints(plot = p.clAv[[i]], points = genes.from.bulk.DE, repel = TRUE, size=2);
  }

  PlotIter <- CodeAndRoll2::split_vec_to_list_by_N(1:NrPlots, by = 4)
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


# plot.Metadata.Cor.Heatmap ------------------------------------------------------------------------
plot.Metadata.Cor.Heatmap <- function(
  columns = c( "nCount_RNA", "nFeature_RNA", "percent.mito", "percent.ribo")
  , cormethod = c('pearson', 'spearman')[1]
  , main =paste( "Metadata", cormethod,"correlations")
  , obj = combined.obj
  , w = 10, h = w
  , ...){
  library(ggcorrplot)


  meta.data <- obj@meta.data
  columns.found <- intersect(colnames(meta.data), columns)

  corX <- cor(meta.data[ , columns.found], method = cormethod)
  pl <- ggcorrplot(corX, hc.order = TRUE, title = main
                   , type = "full", lab = T)
  qqSave(pl, fname = ppp(make.names(main),'pdf'), w = w, h = h)
  pl
}




# plot.Metadata.median.fraction.barplot ------------------------------------------------------------------------
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
  qqSave(pl, fname = ppp(make.names(main),'pdf'), w = w, h = h)
  pl
  if (return.matrix) mat.cluster.medians1 else pl
}

# plot.Metadata.median.fraction.barplot()



# plot.Gene.Cor.Heatmap ------------------------------------------------------------------------
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

  qname = p0("expr.q", quantileX * 100)
  slotname_cor.mat <- kpp('cor', slot.use, assay.use, qname)
  cor.mat <- obj@misc[[slotname_cor.mat]]

  if (is_null(cor.mat)) {
    iprint(slotname_cor.mat, " not found in @misc.")
    iprint("Correlation slots present in @misc:",grepv(names(combined.obj@misc), pattern = "^cor"))

    # Calculate ------------------------------------
    if (calc.COR) {
      print("Calculating correlation now.")
      genes.found <- check.genes(genes)
      iprint(l(genes.found), "genes are found in the object.")
      if (l(genes.found) > 200) iprint("Too many genes found in data, cor will be slow: ", l(genes.found))
      ls.cor <- sparse.cor(t(expr.mat[genes.found,]))
      cor.mat <- ls.cor$cor
    } else { stop() }
  } else {
    print("Correlation is pre-calculated")
    genes.found <- intersect(genes, rownames(cor.mat))
    iprint(l(genes.found), "genes are found in the correlation matrix.")
    cor.mat <- cor.mat[genes.found, genes.found]
  }


  # Filter ------------------------------------
  diag(cor.mat) <- NaN
  corgene.names <- union(
    which_names(rowMax(cor.mat) >= min.g.cor),
    which_names(rowMin(cor.mat) <= -min.g.cor)
  )
  iprint(l(corgene.names), "genes are more (anti-)correlated than +/-:", min.g.cor)

  pname = p0("Pearson correlations of ", substitute(genes),"\n min.cor:", min.g.cor, " | ",  assay.use ,'.', slot.use )
  o.heatmap <- pheatmap(cor.mat[corgene.names,corgene.names],main = pname, cutree_rows = cutRows, cutree_cols = cutCols, ...)
  wplot_save_pheatmap(o.heatmap, filename = make.names(pname))

  # return values
  maxCorrz <- rowMax(cor.mat)[corgene.names]; names(maxCorrz) <- corgene.names
  dput(maxCorrz)
}

# Calc.Cor.Seurat ------------------------------------------------
Calc.Cor.Seurat <- function(assay.use = "RNA", slot.use = "data"
                            , quantileX = 0.95, max.cells =  40000, seed = p$"seed"
                            , digits = 2, obj = combined.obj) {
  expr.mat <- GetAssayData(slot = slot.use, assay = assay.use, object = obj)
  if (ncol(expr.mat) > max.cells) {
    set.seed(seed = seed)
    cells.use <- sample(x = colnames(expr.mat), size = max.cells)
  }

  qname = p0("q", quantileX * 100)
  quantile_name = kpp("expr", qname)

  if (is.null(obj@misc[[quantile_name]])) iprint("Call: combined.obj <- calc.q90.Expression.and.set.all.genes(combined.obj, quantileX =",quantileX," first )")
  genes.HE = which_names(obj@misc[[quantile_name]] > 0)
  iprint("Pearson correlation is calculated for", l(genes.HE), "HE genes with expr.",qname,": > 0.")
  tic(); ls.cor <- sparse.cor(smat = t(expr.mat[genes.HE, cells.use])); toc()
  ls.cor <- lapply(ls.cor, round, digits = 2)

  slot__name <- kpp(slot.use, assay.use, quantile_name)
  obj@misc[[kpp('cor', slot__name)]] <- ls.cor$'cor'
  obj@misc[[kpp('cov', slot__name)]] <- ls.cor$'cov'
  iprint("Stored under obj@misc$", kpp('cor', slot.use, assay.use), "or cov... ." )
  return(obj)
}
# combined.obj <- calc.q90.Expression.and.set.all.genes(combined.obj, quantileX = 0.99, max.cells =  400000, set.all.genes = F)
# combined.obj <- Calc.Cor.Seurat(assay.use = "RNA", slot.use = "data", digits = 2, obj = combined.obj, quantile = 0.99, max.cells = 40000)

# plot.clust.size.distr ------------------------------------------------
plot.clust.size.distr <- function(obj = combined.obj, ident = GetClusteringRuns()[2]
                                  , plot = T, thr.hist = 30, ...) {
  clust.size.distr <- table(obj@meta.data[,ident])
  print(clust.size.distr)
  resX <- gsub(pattern = ".*res\\.", replacement = '',x = ident)
  ptitle <- ppp('clust.size.distr', ident)
  psubtitle <- paste("Nr.clusters:", l(clust.size.distr)
                     , "| median:", median(clust.size.distr)
                     , "| CV:", percentage_formatter(cv(clust.size.distr))
  )
  xlb = "Cluster size (cells)"
  xlim = c(0, max(clust.size.distr))

  if (plot) {
    if (l(clust.size.distr) < thr.hist) {
      qbarplot(clust.size.distr, plotname = ptitle, subtitle = psubtitle, xlab = xlb, ...)
    } else {
      qhistogram(vec = clust.size.distr, plotname = ptitle, subtitle = psubtitle, xlab = xlb, xlim = xlim, ...)
    }
  } else {    "return vector"
    clust.size.distr
  }

}
# plot.clust.size.distr()



#  ------------------------------------------------
geneExpressionLevelPlots <- function(gene = 'TOP2A', obj = ls.Seurat[[1]], slot = c('counts', 'data')[2] ) {
  slot = 'data'
  print(gene)
  if (gene %in% rownames(obj)) {
    GEX.Counts <- GetAssayData(object = obj, assay = 'RNA', slot = slot)

    GEX.Counts.total <- rowSums(GEX.Counts)
    genes.expression <- GEX.Counts.total[gene]
    mean.expr <- iround(mean(GEX.Counts[gene,]))

    suffx = if (slot == 'counts') 'raw' else 'normalised, logtransformed'
    (pname = paste(gene, 'and the', suffx,'transcript count distribution'))

    qhistogram(GEX.Counts.total, vline = genes.expression, logX = T, w = 6, h = 4
               , subtitle = paste('It belong to the top', pc_TRUE(GEX.Counts.total > genes.expression), 'of genes (black line). Mean expr:', mean.expr)
               , plotname = pname, xlab = 'Total Transcripts in Dataset', ylab = 'Number of Genes')
  } else { print("     !!! Gene not found in object!")}
}

#  ------------------------------------------------
PrctCellExpringGene <- function(genes, group.by = "all", object = combined.obj){ # From Github/Ryan-Zhu https://github.com/satijalab/seurat/issues/371
  if(group.by == "all"){
    prct = unlist(lapply(genes, ww.calc_helper, object = object))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }

  else{
    list = SplitObject(object, group.by)
    factors = names(list)
    results = lapply(list, PrctCellExpringGene, genes = genes)
    for (i in 1:length(factors)) {
      results[[i]]$Feature = factors[i]
    }
    combined = do.call("rbind", results)
    return(combined)
  }
}


#  ------------------------------------------------
ww.calc_helper <- function(object, genes){ # From Github/Ryan-Zhu https://github.com/satijalab/seurat/issues/371
  counts = object[['RNA']]@counts
  ncells = ncol(counts)
  if (genes %in% row.names(counts)) {
    sum(counts[genes, ] > 0) / ncells
  } else{
    return(NA)
  }
}

#  ------------------------------------------------

# scBarplotFractionAboveThr <- function(thrX = 0.01, columns.used = c('cl.names.top.gene.res.0.3', 'percent.ribo')
#                                       , obj = combined.obj, ) { # Calculat the fraction of cells per cluster above a certain threhold
#   meta = obj@meta.data
#   (fr_ribo_low_cells <- meta %>%
#       dplyr::select(columns.used)  %>%
#       dplyr::group_by_(columns.used[1])  %>%
#       summarize(n_cells = n(),
#                 n_ribo_low_cells = sum(!!as.name(columns.used[2]) < thrX),
#                 fr_ribo_low_cells = n_ribo_low_cells / n_cells) %>%
#       FirstCol2RowNames())
#   print(fr_ribo_low_cells)
#
#   (v.fr_ribo_low_cells <- 100* as.named.vector(fr_ribo_low_cells[3]))
#   qbarplot(v.fr_ribo_low_cells, xlab.angle = 45, xlab = 'Clusters', ylab = '% Cells')
# }




scBarplotFractionAboveThr <- function(thrX = 0., value.col = 'TVA', id.col =  'cl.names.top.gene.res.0.3'
                                      , obj = combined.obj, return.df = F) { # Calculat the fraction of cells per cluster above a certain threhold
  meta = obj@meta.data
  (df_cells_above <- meta %>%
      dplyr::select(c(id.col, value.col))  %>%
      dplyr::group_by_(id.col)  %>%
      summarize(n_cells = n(),
                n_cells_above = sum(!!as.name(value.col) > thrX),
                fr_n_cells_above = n_cells_above / n_cells) %>%
      FirstCol2RowNames())

  (v.fr_n_cells_above <- 100* as.named.vector(df_cells_above[3]))
  ggobj <- qbarplot(v.fr_n_cells_above, xlab = 'Clusters', ylab = '% Cells'
                    , plotname = paste('Cells with', value.col, '>', thrX)
                    , subtitle = id.col, xlab.angle = 45)
  if (return.df) return(df_cells_above) else ggobj
}

# combined.obj$'TVA' = combined.obj@assays$RNA@data['AAV.TVA.p2a.G.p2a.GFP.WPRE.bGH',]
# scBarplotFractionAboveThr(id.col =  'cl.names.top.gene.res.0.3', value.col = 'TVA', thrX = 0)


#  ------------------------------------------------
scBarplotFractionBelowThr <- function(thrX = 0.01, value.col = 'percent.ribo', id.col =  'cl.names.top.gene.res.0.3'
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
  ggobj <- qbarplot(v.fr_n_cells_below, xlab = 'Clusters', ylab = '% Cells'
                    , plotname = make.names(paste('Cells with', value.col, '<', thrX))
                    , subtitle = id.col, xlab.angle = 45)
  if (return.df) return(df_cells_below) else ggobj
}

# scBarplotFractionBelowThr(id.col =  'cl.names.top.gene.res.0.3', value.col = 'percent.ribo', thrX = 0.01, return.df = T)



#  ------------------------------------------------

#  ------------------------------------------------
