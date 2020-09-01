######################################################################
# Soup.Analysis.of.ambient.RNA.R
######################################################################
# source ('~/GitHub/Seurat.utils/Soup.Analysis.of.ambient.RNA.R')
# Source: self + web

# Requirements ------------------------

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
  CR.matrices$'raw' <- Read10X(path.raw, )
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
