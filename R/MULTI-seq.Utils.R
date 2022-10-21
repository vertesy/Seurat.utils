# _________________________________________________________________________________________________
# MULTI-seq.Utils.R ______________________________ ----
# _________________________________________________________________________________________________
# source('~/GitHub/Packages/Seurat.utils/R/MULTI-seq.Utils.R')
# try(source('https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/MULTI-seq.functions.R'))

# Requirements __________________________________________
# try(require(MarkdownReports),  silent = T)
# try(require(pheatmap),  silent = T)
# May also require


# _________________________________________________________________________________________________
# BarTableSweepList
#' @title BarTableSweepList
#' @description BarTableSweepList
#' @param min PARAM_DESCRIPTION, Default: 0.01
#' @param max PARAM_DESCRIPTION, Default: 0.99
#' @param step PARAM_DESCRIPTION, Default: 0.02
#' @param bar_table PARAM_DESCRIPTION, Default: bar.table
#' @examples
#' \dontrun{
#' if(interactive()){
#'  bar.table_sweep.list <- BarTableSweepList(bar_table = bar.table.solo) # Quantile Sweep List
#'  }
#' }
#' @export
BarTableSweepList <- function(min = 0.01, max = 0.99, step = 0.02, bar_table =bar.table) {
  bar.table_sweep.list <- list()
  n <- 0
  Quantiles = seq(from = min, to = max, by = step)
  for (n in 1:length(Quantiles)) { # print(q)
    bar.table_sweep.list[[n]] <- classifyCells(bar_table, q = Quantiles[n])
    names(bar.table_sweep.list)[n] <- paste("q=",Quantiles[n], sep="")
  }
  return(bar.table_sweep.list)
}





# _________________________________________________________________________________________________
#' @title mSeq.map.all96.BCs
#' @description mSeq.map.all96.BCs
#' @param readTable PARAM_DESCRIPTION, Default: readTable
#' @param CellIDs PARAM_DESCRIPTION, Default: CellIDs
#' @param path2allBCs PARAM_DESCRIPTION, Default: '~/Google_Drive/Science/IMBA/MULTI.seq/from.US/All.MULTI-seq_barcodes.Mar2019.tsv'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  bar.table <-mSeq.map.all96.BCs(readTable = readTable, CellIDs = CellIDs)
#'  }
#' }
#' @export
mSeq.map.all96.BCs <- function(readTable = readTable, CellIDs = CellIDs
                               , path2allBCs = '~/Google_Drive/Science/IMBA/MULTI.seq/from.US/All.MULTI-seq_barcodes.Mar2019.tsv'
) {
  (bar.ref <- read_tsv(path2allBCs)[[1]]) # Vector of reference all MULTI-seq sample barcode sequences.
  MULTIseq.align(readTable = readTable, cellIDs = CellIDs, ref = bar.ref)
}



# _________________________________________________________________________________________________

#' @title aux_plotAllMseqBCs
#' @description aux_plotAllMseqBCs
#' @param bar.table PARAM_DESCRIPTION, Default: bar.table[, 1:96]
#' @param barcodes.used PARAM_DESCRIPTION, Default: BCs.used
#' @param plotname Title of the plot, Default: 'Barcode seq depth'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  aux_plotAllMseqBCs(bar.table = bar.table[,1:96], barcodes.used = BCs.used, plotname = "Barcode seq depth")
#'  }
#' }
#' @export
aux_plotAllMseqBCs <- function(bar.table = bar.table[,1:96], barcodes.used = BCs.used
                               , plotname = "Barcode seq depth") {
  stopifnot(is.numeric(BCs.used))
  BC.depth <- colSums(bar.table)[1:96]
  if (min(BC.depth) < 1) { BC.depth <- BC.depth+1 }
  log.depth <- log10(BC.depth); range(log.depth)
  ccc <- colnames(bar.table) %in% BCs.used

  wbarplot(log.depth, col = ccc, lwd = 1, lcol = 1, lty = 2, plotname = plotname
           , vline = (range(BCs.used)+c(-1,1))
           , hline = quantile(log.depth[setdiff(1:69, BCs.used)], .95)
           , ylab = "log10(total reads / BC)", main =  plotname)
  wlegend.label(
    "    Horiz. line at 95%
    of unused BC's.
    Vertical line: range
    of used BCs.",poz = 2, cex = 1)
}


# _________________________________________________________________________________________________
# bar.table.log <- t(log10(bar.table[,BCs.used]+1))
# bar.table.log <- CodeAndRoll2::clip.outliers.at.percentile(bar.table.log)

