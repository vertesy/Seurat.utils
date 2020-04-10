######################################################################
# plotting.statistics.and.QC.R
######################################################################
# source ('~/GitHub/Seurat.utils/plotting.statistics.and.QC.R')
# Source: self + web

# Requirements ------------------------
require(Seurat)
require(ggplot2)
require(SoupX)
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
                                  , plotname = paste(TitleCase(fill.by), "proportions"), seedNr = 1989) {
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

    geom_hline(yintercept=c(.25, .5, .75))  +
    geom_bar( position="fill" ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_text(aes(label=..count..), stat='count',position=position_fill(vjust=0.5)) +
    labs(title =pname_, x = "Clusters", y = "Fraction", caption = caption_)
}
# CellFractionsBarplot2(obj = combined.obj, group.by = "integrated_snn_res.0.1", fill.by = "Phase", downsample = T)
# CellFractionsBarplot2(obj = combined.obj, group.by = "integrated_snn_res.0.1", fill.by = "Phase", downsample = F)


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
