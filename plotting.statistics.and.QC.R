######################################################################
# plotting.statistics.and.QC.R
######################################################################
# source ('~/GitHub/Seurat.utils/plotting.statistics.and.QC.R')
# Source: self + web

# Requirements ------------------------
require(Seurat)
require(ggplot2)
try(require(SoupX))
# May also require
# try (source ('~/GitHub/CodeAndRoll/CodeAndRoll.R'),silent= F) # generic utilities funtions
# require('MarkdownReportsDev') # require("devtools") # plotting related utilities functions # devtools::install_github(repo = "vertesy/MarkdownReportsDev")


# PCA percent of variation associated with each PC ------------------------------------------------------------
seu.PC.var.explained <- function(obj =  combined.obj) { # Determine percent of variation associated with each PC.
  pct <- obj@reductions$pca@stdev / sum(obj@reductions$pca@stdev) * 100
  names(pct) =1:length(obj@reductions$pca@stdev)
  return(pct)
}

# plot percent of variation associated with each PC ------------------------------------------------------------
seu.plot.PC.var.explained <- function(obj =  combined.obj) { # Plot the percent of variation associated with each PC.
  pct <- seu.PC.var.explained(obj)
  wbarplot(pct, ylab= "% of variation explained" , xlab="Principal Components")
  barplot_label(round(pct, digits = 2), barplotted_variable = pct, cex=.5 )
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
                                  , plotname = paste(TitleCase(fill.by), "proportions"), hlines = c(.25, .5, .75), seedNr = 1989) {
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
    { if(downsample) sample_n(., downsample) else . } %>%
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
    avg.ClX.cells <- log2(AverageExpression(ClX, verbose = FALSE)$RNA+1)
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
    qqSaveGridA4(plotlist = plotLS, plots = 1:4, fname = ppp("BulkGEScatterPlot.BulkGenes",kpp(PlotIter[[i]]), "png"))
  }
}
# BulkGEScatterPlot(obj = combined.obj, clusters = "cl.names.KnownMarkers.0.2", TwoCategIdent = 'age', genes.from.bulk.DE = rownames(df.markers.per.AGE))


# plotTheSoup ------------------------------------------------------------------------
plotTheSoup <- function(CellR.OutputDir = "~/Dropbox/Abel.IMBA/Data/SoupX_pbmc4k_demo/") { # Plot the ambient RNA content of droplets without a cell (background droplets).
  # Read In ------------------------
  sc = load10X(CellR.OutputDir, keepDroplets = TRUE)
  # Profiling the soup ------------------------
  sc = estimateSoup(sc)

  # Plot top gene's expression ----------------------------------------------------------------
  soupProfile = head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20)
  soupX.cellfree.RNA.profile = 100 * col2named.vector(soupProfile[,1,drop=F])
  wbarplot(soupX.cellfree.RNA.profile
           , ylab="% Reads in the Soup"
           , sub = paste("Within the", basename(CellR.OutputDir), "dataset")
           , tilted_text = T)
  barplot_label(barplotted_variable = soupX.cellfree.RNA.profile
                , labels = percentage_formatter(soupX.cellfree.RNA.profile/100, digitz = 2)
                , TopOffset = .4, srt = 90, cex=.75)

  # Plot summarize expression ----------------------------------------------------------------
  soup.RP.sum   <- colSums(soupProfile[grep('^RPL|^RPS', rownames(soupProfile)),])
  soup.RPL.sum   <- colSums(soupProfile[grep('^RPL', rownames(soupProfile)),])
  soup.RPS.sum   <- colSums(soupProfile[grep('^RPS', rownames(soupProfile)),])
  soup.mito.sum <- colSums(soupProfile[grep('^MT-', rownames(soupProfile)),])

  soupProfile.summarized <- rbind(
    'Ribosomal' = soup.RP.sum,
    'Ribosomal.L' = soup.RPL.sum,
    'Ribosomal.S' = soup.RPS.sum,
    'Mitochondial' = soup.mito.sum,
    soupProfile[grep('^RPL|^RPS|^MT-', rownames(soupProfile), invert = T),]
  )

  NrColumns2Show  = min(10, nrow(soupProfile.summarized))
  soupX.cellfree.RNA.profile.summarized = 100 * col2named.vector(soupProfile.summarized[1:NrColumns2Show,1,drop=F])
  wbarplot(soupX.cellfree.RNA.profile.summarized
           , ylab="% Reads in the Soup"
           , sub = paste("Within the", basename(CellR.OutputDir), "dataset")
           , tilted_text = T)
  barplot_label(barplotted_variable = soupX.cellfree.RNA.profile.summarized
                , labels = percentage_formatter(soupX.cellfree.RNA.profile.summarized/100, digitz = 2)
                , TopOffset = .5)
  remove("sc")
  detach(SoupX)
} # plotTheSoup
