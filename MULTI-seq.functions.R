######################################################################
# MULTI-seq.functions.R
######################################################################
# source ('~/GitHub/Seurat.utils/MULTI-seq.functions.R')

# Requirements ------------------------
require(MarkdownReportsDev)
require(pheatmap)
# May also require


# BarTableSweepList -----------------------------------------------------------------------
BarTableSweepList <- function(min=0.01, max=0.99, step=0.02, bar_table =bar.table) {
  bar.table_sweep.list <- list()
  n <- 0
  Quantiles = seq(from = min, to = max, by=step)
  for (n in 1:length(Quantiles)) { # print(q)
    bar.table_sweep.list[[n]] <- classifyCells(bar_table, q=Quantiles[n])
    names(bar.table_sweep.list)[n] <- paste("q=",Quantiles[n], sep="")
  }
  return(bar.table_sweep.list)
}
# bar.table_sweep.list <- BarTableSweepList(bar_table = bar.table.solo) # Quantile Sweep List




# mSeq.map.all96.BCs ------------------------------------------------------------------------
mSeq.map.all96.BCs <- function(readTable = readTable, CellIDs = CellIDs
                               , path2allBCs = '~/Google_Drive/Science/IMBA/MULTI.seq/from.US/All.MULTI-seq_barcodes.Mar2019.tsv'
) {
  (bar.ref <- read_tsv(path2allBCs)[[1]]) # Vector of reference all MULTI-seq sample barcode sequences.
   MULTIseq.align(readTable = readTable, cellIDs = CellIDs, ref = bar.ref)
}
# bar.table <-mSeq.map.all96.BCs(readTable = readTable, CellIDs = CellIDs)


# aux.plotAllMseqBCs ------------------------------------------------------------------------

aux.plotAllMseqBCs <- function(bar.table = bar.table[,1:96], barcodes.used = BCs.used
                               , plotname = "Barcode seq depth") {
  stopifnot(is.numeric(BCs.used))
  BC.depth <- colSums(bar.table)[1:96]
  if (min(BC.depth) < 1) { BC.depth <- BC.depth+1 }
  log.depth <- log10(BC.depth); range(log.depth)
  ccc <- colnames(bar.table) %in% BCs.used

  wbarplot(log.depth, col = ccc, lwd=1, lcol=1, lty=2, plotname = plotname
           , vline = (range(BCs.used)+c(-1,1))
           , hline = quantile(log.depth[setdiff(1:69, BCs.used)], .95)
           , ylab = "log10(total reads / BC)", main =  plotname)
  wlegend.label(
    "    Horiz. line at 95%
    of unused BC's.
    Vertical line: range
    of used BCs.",poz = 2, cex=1)
}
# aux.plotAllMseqBCs(bar.table = bar.table[,1:96], barcodes.used = BCs.used, plotname = "Barcode seq depth")

# ------------------------------------------------------------------------
# bar.table.log <- t(log10(bar.table[,BCs.used]+1))
# bar.table.log <- clip.outliers(bar.table.log)

# p$'dist' <- c("correlation", "manhattan") # 'canberra', 'binary', 'minkowski', "euclidean")
# for (i in 1:l(p$'dist')) {
#   distX <- p$'dist'[i]; print(distX)
#   pheatmap(bar.table.log, clustering_distance_cols = distX# , cluster_rows = FALSE
#            , show_colnames = F, cutree_cols = nrow(bar.table.log), treeheight_col = 0
#            , filename = kpp("Barcode.log10.RC", distX,"heatmap", idate(Format = "%Y.%m.%d_%H.%M.%S"), "png"))
# }
# oo()


# subset.ReadTable <- function(variables) {
#
# }
# ------------------------------------------------------------------------
