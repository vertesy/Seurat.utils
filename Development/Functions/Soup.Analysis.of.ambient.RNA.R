######################################################################
# Soup.Analysis.of.ambient.RNA.R
######################################################################
# source('~/GitHub/Packages/Seurat.utils/Functions/Soup.Analysis.of.ambient.RNA.R')
# try (source('https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Soup.Analysis.of.ambient.RNA.R'))
# Source: self + web

# Requirements ------------------------
require(tibble)

# plotTheSoup ------------------------------------------------------------------------
plotTheSoup <- function(CellRangerOutputDir = "~/Data/114593/114593"
                        , SeqRun = gsub('*([0-9]+).*','\\1', x = basename(CellRangerOutputDir))) { # Plot the ambient RNA content of droplets without a cell (background droplets).

  ls.Alpha=1
  # Setup ------------------------
  require(Matrix); require(ggrepel)

  dirz <- list.dirs(CellRangerOutputDir, full.names = F, recursive = F)
  path.raw <- file.path(CellRangerOutputDir, grep(x = dirz, pattern = "^raw_*", value = T))
  path.filt <- file.path(CellRangerOutputDir, grep(x = dirz, pattern = "^filt_*", value = T))
  CR.matrices <- list.fromNames(c("raw", "filt"))

  # Adapter for Markdownreports background variable "OutDir" ----------------------------------------------------------------
  if (exists('OutDir')) OutDirBac <- OutDir
  OutDir <- file.path(CellRangerOutputDir,p0(kpp("SoupStatistics", SeqRun)))
  try(dir.create(OutDir))
  ww.assign_to_global("OutDir", OutDir, 1)

  # Read In ------------------------
  print("Reading raw CellRanger output matrices")
  CR.matrices$'raw' <- Read10X(path.raw)
  if (length(CR.matrices$'raw') == 2 ) { CR.matrices$'raw' <- CR.matrices$'raw'[[1]] } # Maybe AB table is present too at slot 2!

  print("Reading filtered CellRanger output matrices")
  CR.matrices$'filt' <- Read10X(path.filt)
  if (length(CR.matrices$'filt') == 2 ) { CR.matrices$'filt' <- CR.matrices$'filt'[[1]] } # Maybe AB table is present too at slot 2!

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
  # wplot(Soup.VS.Cells.Av.Exp.log10, col = rgb(0,0,0,.25), PNG = T
  #       , xlab = "Total Expression in Soup [log10(mRNA+1)]"
  #       , ylab = "Total Expression in Cells [log10(mRNA+1)]"
  # )

  # ggplot prepare ----------------------------------------------------------------
  Soup.VS.Cells.Av.Exp.gg <- tibble::rownames_to_column(as.data.frame(Soup.VS.Cells.Av.Exp.log10), "gene")
  (Soup.VS.Cells.Av.Exp.gg <- as_tibble(Soup.VS.Cells.Av.Exp.gg))
  soup.rate <- Soup.VS.Cells.Av.Exp.gg$Soup / (Soup.VS.Cells.Av.Exp.gg$Cells + Soup.VS.Cells.Av.Exp.gg$Soup)
  cell.rate <- Soup.VS.Cells.Av.Exp.gg$Cells / (Soup.VS.Cells.Av.Exp.gg$Cells + Soup.VS.Cells.Av.Exp.gg$Soup)

  axl.pfx <- "Total Expression in"
  axl.sfx <- "[log10(mRNA+1)]"



  HGNC <- Soup.VS.Cells.Av.Exp.gg$gene
  Class <- rep("Other", times=nrow(Soup.VS.Cells.Av.Exp.gg))
  Class[grep('^RPL|^RPS', HGNC)]  <- "RP"
  Class[grep('^MT-', HGNC)] <- "MT"
  Class[grep('^LINC', HGNC)]  <- "LINC"
  Class[grep('^AC', HGNC)]  <- "AC"
  Class[grep('^AL', HGNC)]  <- "AL"
  wpie(table(Class))
  Soup.VS.Cells.Av.Exp.gg$Class <- Class

  fname <- kpp("Soup.VS.Cells.Av.Exp.GeneClasses",SeqRun,"pdf")
  pgg <-
    ggplot(Soup.VS.Cells.Av.Exp.gg %>% arrange(-nchar(Class) )
           , aes(x= Soup, y= Cells, label=gene, col= Class))  +
    geom_abline(slope=1, col='darkgrey') + geom_point()+
    scale_alpha_manual(guide='none', values = ls.Alpha) +
    xlab(paste(axl.pfx, "Soup", axl.sfx)) + ylab(paste(axl.pfx, "Cells", axl.sfx)) +
    ggtitle("Soup VS. Cells | gene classes")

  ggsave(pgg, filename = file.path(OutDir, fname))

  # ggplot ----------------------------------------------------------------
  quantiles <- c(0.025, 0.01, 0.0025)

  i=1
  for (i in 1:l(quantiles)) {
    pr <- quantiles[i]; print(pr)
    HP.thr <- 200*pr/quantiles[2]
    idx.HE2 <- rowSums(Soup.VS.Cells.Av.Exp) > HP.thr
    pc_TRUE(idx.HE2)

    fname <- kpp("Soup.VS.Cells.Av.Exp.quantile",pr,SeqRun,"pdf")

    Outlier <- idx.HE2 &
      (cell.rate < quantile(cell.rate, probs = pr) |
         soup.rate < quantile(soup.rate, probs = pr))

    pc_TRUE(Outlier); sum(Outlier)
    HP.thr.mod <- HP.thr
    while (sum(Outlier) > 40) {
      HP.thr.mod <- HP.thr.mod *2
      Outlier <- Outlier &  rowSums(Soup.VS.Cells.Av.Exp) > HP.thr.mod
    }
    sum(Outlier)


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
# plotTheSoup(CellRangerOutputDir = "~/Data/114593/114593" , SeqRun = gsub('*([0-9]+).*','\\1', x = basename(CellRangerOutputDir)))



load10Xv3 <- function(dataDir, cellIDs = NULL, channelName = NULL, readArgs = list(),
          includeFeatures = c("Gene Expression"), verbose = TRUE,
          ...)
{

  # include
  dirz <- list.dirs(dataDir, full.names = F, recursive = F)
  path.raw <- file.path(dataDir, grep(x = dirz, pattern = "^raw_*", value = T))
  path.filt <- file.path(dataDir, grep(x = dirz, pattern = "^filt_*", value = T))
  CR.matrices <- list.fromNames(c("raw", "filt"))


  (isV3 = any(grepl(x = dirz, pattern = "^raw_feature_bc*")))
  tgt = path.raw

  if (!isV3)
    tgt = file.path(tgt, list.files(tgt))
  if (verbose)
    message(sprintf("Loading raw count data"))
  dat = do.call(Read10X, c(list(data.dir = tgt), readArgs))
  if (verbose)
    message(sprintf("Loading cell-only count data"))
  if (!is.null(cellIDs)) {
    if (all(grepl("\\-1$", cellIDs)))
      cellIDs = gsub("\\-1$", "", cellIDs)
    if (!all(cellIDs %in% colnames(dat)))
      stop("Not all supplied cellIDs found in raw data.")
    datCells = dat[, match(cellIDs, colnames(dat))]
  }
  else {
    tgt = path.filt
    if (!isV3)
      tgt = file.path(tgt, list.files(tgt))
    datCells = do.call(Read10X, c(list(data.dir = tgt),
                                  readArgs))
    if (is.list(dat)) {
      dat = do.call(rbind, dat[includeFeatures])
      datCells = do.call(rbind, datCells[includeFeatures])
    }
  }
  if (verbose)
    message(sprintf("Loading extra analysis data where available"))
  mDat = NULL
  tgt = file.path(dataDir, "analysis", "clustering", "graphclust",
                  "clusters.csv")
  if (file.exists(tgt)) {
    clusters = read.csv(tgt)
    mDat = data.frame(clusters = clusters$Cluster, row.names = clusters$Barcode)
  }
  tgt = file.path(dataDir, "analysis", "clustering", "kmeans_10_clusters",
                  "clusters.csv")
  if (file.exists(tgt)) {
    clusters = read.csv(tgt)
    mDat$clustersFine = clusters$Cluster
  }
  tgt = file.path(dataDir, "analysis", "tsne", "2_components",
                  "projection.csv")
  if (file.exists(tgt)) {
    tsne = read.csv(tgt)
    if (is.null(mDat)) {
      mDat = data.frame(tSNE1 = tsne$TSNE.1, tSNE2 = tsne$TSNE.2,
                        row.names = tsne$Barcode)
    }
    else {
      mDat$tSNE1 = tsne$TSNE.1[match(rownames(mDat), tsne$Barcode)]
      mDat$tSNE2 = tsne$TSNE.2[match(rownames(mDat), tsne$Barcode)]
    }
    DR = c("tSNE1", "tSNE2")
  }
  else {
    DR = NULL
  }
  if (!is.null(mDat) && any(rownames(mDat) != colnames(datCells))) {
    rownames(mDat) = gsub("-1$", "", rownames(mDat))
    if (any(rownames(mDat) != colnames(datCells)))
      stop("Error matching meta-data to cell names.")
  }
  if (is.null(channelName))
    channelName = ifelse(is.null(names(dataDir)), dataDir,
                         names(dataDir))
  channel = SoupChannel(tod = dat, toc = datCells, metaData = mDat,
                        channelName = channelName, dataDir = dataDir, dataType = "10X",
                        isV3 = isV3, DR = DR, ...)
  return(channel)
}




# dataDir="/Volumes/copy.your.own.data.here/A.Vertesy/SEO/HNV73DRXX_R10015/HNV73DRXX_R10015/aligned_rna/124851_rnacount"

