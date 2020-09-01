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
plotTheSoup <- function(CellRangerOutputDir = "~/Data/114593/114593"
                        , SeqRun = gsub('*([0-9]+).*','\\1', x = basename(CellRangerOutputDir))) { # Plot the ambient RNA content of droplets without a cell (background droplets).

  # Setup ------------------------
  require(Matrix); require(ggrepel)

  dirz <- list.dirs(CellRangerOutputDir, full.names = F, recursive = F)
  path.raw <- file.path(CellRangerOutputDir, grep(x = dirz, pattern = "^raw_*", value = T))
  path.filt <- file.path(CellRangerOutputDir, grep(x = dirz, pattern = "^filt_*", value = T))
  CR.matrices <- list.fromNames(c("raw", "filt"))

  # Adapter for Markdownreports background variable "OutDir" ----------------------------------------------------------------
  if(exists('OutDir')) OutDirBac <- OutDir
  OutDir <- file.path(CellRangerOutputDir,p0(kpp("SoupStatistics", SeqRun)))
  try(dir.create(OutDir))
  ww.assign_to_global("OutDir", OutDir, 1)


  # Read In ------------------------
  print("Reading raw CellRanger output matrices")
  CR.matrices$'raw' <- Read10X(path.raw)
  print("Reading filtered CellRanger output matrices")
  CR.matrices$'filt' <- Read10X(path.filt)

  # Profiling the soup ------------------------
  print("Profiling the soup")
  GEMs.all <- CR.matrices$'raw'@Dimnames[[2]]
  GEMs.cells <- CR.matrices$'filt'@Dimnames[[2]]
  iprint("There are", l(GEMs.all), "GEMs sequenced, and",l(GEMs.cells), "are cells among those." )

  GEMs.soup <- setdiff(GEMs.all, GEMs.cells)
  CR.matrices$'soup' <- CR.matrices$'raw'[,GEMs.soup]
  CR.matrices$'soup.total.RC' <- Matrix::rowSums(CR.matrices$'soup')
  CR.matrices$'soup.total.sum' <- sum(CR.matrices$'soup')
  CR.matrices$'cells.total.sum' <- sum(CR.matrices$'filt')

  CR.matrices$'soup.rel.RC'  <- CR.matrices$'soup.total.RC' / CR.matrices$'soup.total.sum'

  # Diff Exp ----------------------------------------------------------------
  Soup.VS.Cells.Av.Exp <- cbind(
    'Soup' = Matrix::rowSums(CR.matrices$'soup'),
    'Cells' = Matrix::rowSums(CR.matrices$'filt')
  )
  colnames(Soup.VS.Cells.Av.Exp)
  idx.HE <- rowSums(Soup.VS.Cells.Av.Exp)>10; pc_TRUE(idx.HE)
  Soup.VS.Cells.Av.Exp <- Soup.VS.Cells.Av.Exp[idx.HE,]; idim(Soup.VS.Cells.Av.Exp)
  Soup.VS.Cells.Av.Exp.log10 <- log10(Soup.VS.Cells.Av.Exp+1)
  wplot(Soup.VS.Cells.Av.Exp.log10, col = rgb(0,0,0,.25)
        , xlab = "Total Expression in Soup [log10(mRNA+1)]"
        , ylab = "Total Expression in Cells [log10(mRNA+1)]"
  )


  # ggplot prepare ----------------------------------------------------------------
  Soup.VS.Cells.Av.Exp.gg <- tibble::rownames_to_column(as.data.frame(Soup.VS.Cells.Av.Exp.log10), "gene")
  (Soup.VS.Cells.Av.Exp.gg <- as_tibble(Soup.VS.Cells.Av.Exp.gg))
  soup.rate <- Soup.VS.Cells.Av.Exp.gg$Soup / (Soup.VS.Cells.Av.Exp.gg$Cells + Soup.VS.Cells.Av.Exp.gg$Soup)
  cell.rate <- Soup.VS.Cells.Av.Exp.gg$Cells / (Soup.VS.Cells.Av.Exp.gg$Cells + Soup.VS.Cells.Av.Exp.gg$Soup)

  axl.pfx <- "Total Expression in"
  axl.sfx <- "[log10(mRNA+1)]"

  # ggplot ----------------------------------------------------------------
  quantiles <- c(0.025, 0.01, 0.0025)

  i=3
  for (i in 1:l(quantiles)) {
    pr <- quantiles[i]; print(pr)
    HP.thr <- 2*pr/quantiles[2]
    idx.HE2 <- rowSums(Soup.VS.Cells.Av.Exp) > HP.thr

    fname <- kpp("Soup.VS.Cells.Av.Exp.quantile",pr,"pdf")

    Outlier <- idx.HE2 &
      (cell.rate < quantile(cell.rate, probs = pr) |
         soup.rate < quantile(soup.rate, probs = pr))

    pgg <-
      ggplot(Soup.VS.Cells.Av.Exp.gg, aes(x= Soup, y= Cells, label=gene,
                                          col= Outlier))  +
      geom_point() + theme(legend.position = "none") +
      xlab(paste(axl.pfx, "Soup", axl.sfx)) + ylab(paste(axl.pfx, "Cells", axl.sfx)) +
      ggtitle("Soup VS. Cells", subtitle = pr) +
      geom_text_repel(aes(label= ifelse(Outlier
                                        , as.character(gene),'')))
    ggsave(pgg, filename = file.path(OutDir, fname))
  }


  # Per Gene ----------------------------------------------------------------
  PC.mRNA.in.Soup <- sum(CR.matrices$'soup')/sum(CR.matrices$'raw')
  PC.mRNA.in.Cells <- 100*sum(CR.matrices$'filt')/sum(CR.matrices$'raw')
  wbarplot(variable = PC.mRNA.in.Cells, col ="seagreen", plotname = kppd("PC.mRNA.in.Cells", SeqRun)
           , ylim = c(0,100), ylab = "% mRNA in cells"
           , sub = "% mRNA is more meaningful than % reads reported by CR")
  barplot_label(barplotted_variable = PC.mRNA.in.Cells
                , labels = percentage_formatter(PC.mRNA.in.Cells/100, digitz = 2)
                , TopOffset = 10)


  # Plot top gene's expression ----------------------------------------------------------------
  Soup.GEMs.top.Genes = 100*head(sort(CR.matrices$'soup.rel.RC', decreasing = T), n = 20)

  wbarplot(Soup.GEMs.top.Genes, plotname = kppd("Soup.GEMs.top.Genes", SeqRun)
           , ylab="% mRNA in the Soup"
           , sub = paste("Within the", SeqRun, "dataset")
           , tilted_text = T
           , ylim = c(0, max(Soup.GEMs.top.Genes)*1.5))
  barplot_label(barplotted_variable = Soup.GEMs.top.Genes
                , labels = percentage_formatter(Soup.GEMs.top.Genes/100, digitz = 2)
                , TopOffset = -.5, srt = 90, cex=.75)

  # Plot summarize expression ----------------------------------------------------------------
  soupProfile <- CR.matrices$'soup.total.RC'
  {
    soup.RP.sum   <- sum(soupProfile[grep('^RPL|^RPS', names(soupProfile))])
    soup.RPL.sum   <- sum(soupProfile[grep('^RPL', names(soupProfile))])
    soup.RPS.sum   <- sum(soupProfile[grep('^RPS', names(soupProfile))])
    soup.mito.sum <- sum(soupProfile[grep('^MT-', names(soupProfile))])
    soup.LINC.sum <- sum(soupProfile[grep('^LINC', names(soupProfile))])
    soup.AC.sum <- sum(soupProfile[grep('^AC', names(soupProfile))])
    soup.AL.sum <- sum(soupProfile[grep('^AL', names(soupProfile))])
    genes.non.Above <- soupProfile[grepv('^RPL|^RPS|^MT-|^LINC|^AC|^AL', names(soupProfile), invert = T)]
  }
  head(sort(genes.non.Above), n=50)


  soupProfile.summarized <- c(
    'Mitochondial' = soup.mito.sum,
    'Ribosomal' = soup.RP.sum,
    'Ribosomal.L' = soup.RPL.sum,
    'Ribosomal.S' = soup.RPS.sum,
    'GenBank (AC)' = soup.AC.sum,
    'EMBL (AL)' = soup.AL.sum,
    'LINC' = soup.LINC.sum,
    sort(genes.non.Above, decreasing = T)
  )
  NrColumns2Show  = min(10, nrow(soupProfile.summarized))
  ccc <- c("#FF4E00","#778B04","#8ea604","#8ea604","#F5BB00","#F5BB00","#EC9F05",rep(x = "#BF3100", times=NrColumns2Show-6)) # ,"#"


  Soup.GEMs.top.Genes.summarized = 100 * soupProfile.summarized[1:NrColumns2Show] / CR.matrices$'soup.total.sum'
  maxx <- max(Soup.GEMs.top.Genes.summarized)
  wbarplot(Soup.GEMs.top.Genes.summarized, plotname = kppd("Soup.GEMs.top.Genes.summarized", SeqRun)
           , ylab="% mRNA in the Soup", ylim = c(0, maxx+3)
           , sub = paste("Within the", SeqRun, "dataset")
           , tilted_text = T, col = ccc)
  barplot_label(barplotted_variable = Soup.GEMs.top.Genes.summarized
                , srt = 45, labels = percentage_formatter(Soup.GEMs.top.Genes.summarized/100, digitz = 2)
                , TopOffset = -1.5)

  # Absolute.fraction ---------------------------
  Absolute.fraction.soupProfile.summarized <- Soup.GEMs.top.Genes.summarized * PC.mRNA.in.Soup

  maxx <- max(Absolute.fraction.soupProfile.summarized)
  wbarplot(Absolute.fraction.soupProfile.summarized, plotname = kppd("Absolute.fraction.soupProfile.summarized", SeqRun)
           , ylab="% of mRNA in cells", ylim = c(0, maxx*1.33)
           , sub = paste(percentage_formatter(PC.mRNA.in.Soup), "of mRNA counts are in the Soup, in the dataset ", SeqRun)
           , tilted_text = T, col = ccc)
  barplot_label(barplotted_variable = Absolute.fraction.soupProfile.summarized
                , srt = 45, labels = percentage_formatter(Absolute.fraction.soupProfile.summarized/100, digitz = 2)
                # formatC(Absolute.fraction.soupProfile.summarized, format="f", big.mark = " ", digits=0)
                , TopOffset = -maxx*0.15)

  # -----
  Soup.GEMs.top.Genes.non.summarized <- 100* sort(genes.non.Above, decreasing = T)[1:20]/ CR.matrices$'soup.total.sum'
  maxx <- max(Soup.GEMs.top.Genes.non.summarized)
  wbarplot(Soup.GEMs.top.Genes.non.summarized, plotname = kppd("Soup.GEMs.top.Genes.non.summarized", SeqRun)
           , ylab="% mRNA in the Soup"
           , sub = paste("Within the", SeqRun, "dataset")
           , tilted_text = T, col = "#BF3100"
           , ylim = c(0, maxx*1.5))
  barplot_label(barplotted_variable = Soup.GEMs.top.Genes.non.summarized
                , labels = percentage_formatter(Soup.GEMs.top.Genes.non.summarized/100, digitz = 2)
                # , labels = p0(round(1e6 * Soup.GEMs.top.Genes.non.summarized), " ppm")
                , TopOffset = -maxx*0.2, srt = 90, cex=.75)

  # Diff Exp ----------------------------------------------------------------


  # Adapter for Markdownreports background variable "OutDir" ----------------------------------------------------------------
  if (exists('OutDirBac'))  ww.assign_to_global("OutDir", OutDirBac, 1)


} # plotTheSoup
# plotTheSoup()





# # plotTheSoup ------------------------------------------------------------------------
# plotTheSoup.old <- function(CellR.OutputDir = "~/Dropbox/Abel.IMBA/Data/SoupX_pbmc4k_demo/") { # Plot the ambient RNA content of droplets without a cell (background droplets).
#   # Read In ------------------------
#   sc = load10X(CellR.OutputDir, keepDroplets = TRUE)
#   # Profiling the soup ------------------------
#   sc = estimateSoup(sc)
#
#   # Plot top gene's expression ----------------------------------------------------------------
#   soupProfile = head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20)
#   soupX.cellfree.RNA.profile = 100 * col2named.vector(soupProfile[,1,drop=F])
#   wbarplot(soupX.cellfree.RNA.profile
#            , ylab="% Reads in the Soup"
#            , sub = paste("Within the", basename(CellR.OutputDir), "dataset")
#            , tilted_text = T)
#   barplot_label(barplotted_variable = soupX.cellfree.RNA.profile
#                 , labels = percentage_formatter(soupX.cellfree.RNA.profile/100, digitz = 2)
#                 , TopOffset = .4, srt = 90, cex=.75)
#
#   # Plot summarize expression ----------------------------------------------------------------
#   soup.RP.sum   <- colSums(soupProfile[grep('^RPL|^RPS', rownames(soupProfile)),])
#   soup.RPL.sum   <- colSums(soupProfile[grep('^RPL', rownames(soupProfile)),])
#   soup.RPS.sum   <- colSums(soupProfile[grep('^RPS', rownames(soupProfile)),])
#   soup.mito.sum <- colSums(soupProfile[grep('^MT-', rownames(soupProfile)),])
#
#   soupProfile.summarized <- rbind(
#     'Ribosomal' = soup.RP.sum,
#     'Ribosomal.L' = soup.RPL.sum,
#     'Ribosomal.S' = soup.RPS.sum,
#     'Mitochondial' = soup.mito.sum,
#     soupProfile[grep('^RPL|^RPS|^MT-', rownames(soupProfile), invert = T),]
#   )
#
#   NrColumns2Show  = min(10, nrow(soupProfile.summarized))
#   soupX.cellfree.RNA.profile.summarized = 100 * col2named.vector(soupProfile.summarized[1:NrColumns2Show,1,drop=F])
#   wbarplot(soupX.cellfree.RNA.profile.summarized
#            , ylab="% Reads in the Soup"
#            , sub = paste("Within the", basename(CellR.OutputDir), "dataset")
#            , tilted_text = T)
#   barplot_label(barplotted_variable = soupX.cellfree.RNA.profile.summarized
#                 , labels = percentage_formatter(soupX.cellfree.RNA.profile.summarized/100, digitz = 2)
#                 , TopOffset = .5)
#   remove("sc")
#   detach(SoupX)
# } # plotTheSoup
