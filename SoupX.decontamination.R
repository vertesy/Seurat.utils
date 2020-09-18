######################################################################
# SoupX.decontamination.R
######################################################################
source("~/GitHub/Seurat.utils/SoupX.decontamination.R")


# Libraries ------------------------
require(SoupX)
require(Seurat)
require(MarkdownReportsDev)
require(DropletUtils)
require(cowplot)

try(source("~/GitHub/Seurat.utils/Seurat.Utils.Load.R"));
source ('~/GitHub/Seurat.utils/Soup.Analysis.of.ambient.RNA.R')

# Parameters ------------------------
v.parentfolder = c(
  # SEO -----
  # "/Volumes/single.cell.RNA.seq/A.Vertesy/SEO/HNV73DRXX_R10015/HNV73DRXX_R10015/aligned_rna/124851_rnacount",
  # "/Volumes/single.cell.RNA.seq/A.Vertesy/SEO/HNV73DRXX_R10015/HNV73DRXX_R10015/aligned_rna/124719_rnacount"
  # # POL -----
  '123062' = "/Volumes/HN3V3DRXX_R9836/demultiplexed/123062/123062_premRNA_POL/outs",
  '123063' = "/Volumes/HN3V3DRXX_R9836/demultiplexed/123063/123063_premRNA_POL/outs"
    # "/Users/abel.vertesy/Data/POL/filtered_feature_bc_matrix.123062"
    # "/Volumes/single.cell.RNA.seq/C.Bosone/2020.06/HN3V3DRXX_R9836/aligned_rna/123063_rnacount",
  # TSC -----
  # "/Volumes/single.cell.RNA.seq/O.Eichmueller/diffmedia.d110/101146/101146_premRNA_POL2/outs/",
  # "/Volumes/single.cell.RNA.seq/O.Eichmueller/diffmedia.d110/101147/101147_premRNA_POL2/outs/"
  # CON -----
  # "/Volumes/single.cell.RNA.seq/A.Vertesy/CONN/Connectome.pilot.1/03.Fulldepth/HKJHKDRXX_R9539/114593/114593_premRNA_local/outs/",
  # INM -----
  # # "/Volumes/single.cell.RNA.seq/S.Bajaj/HJGHFDRXX_all/aligned/104566",
  # # "/Volumes/single.cell.RNA.seq/S.Bajaj/HJGHFDRXX_all/aligned/104567",
  # # "/Volumes/single.cell.RNA.seq/S.Bajaj/HJGHFDRXX_all/aligned/104568",
  # "/Volumes/single.cell.RNA.seq/S.Bajaj/HJGHFDRXX_all/aligned/104569"
)

i=2
for (i in 1:l(v.parentfolder)) {

  # Setup ------------------------
  InDir <- v.parentfolder[i]
  DataSet <- names(v.parentfolder)[i]
  OutDir <- OutDirDecont <- pps(InDir,ppp("SoupX_decont_filt_matrix",DataSet))
  dir.create(OutDirDecont)

  # Estimate ------------------------
  sc = load10Xv3(InDir, verbose = T)
  sc = autoEstCont(sc, verbose = T);  wplot_save_this(plotname =  "SoupX.contamination.fraciton.rho" )
  out = adjustCounts(sc, verbose = T)

  # Analyze ------------------------

  p.DCX <- plotChangeMap(sc, out, "DCX") + ggtitle("Change in DCX expression due to soup correction")
  p.SATB2 <- plotChangeMap(sc, out, "SATB2") + ggtitle("Change in SATB2 expression due to soup correction")
  p.VIM <- plotChangeMap(sc, out, "VIM") + ggtitle("Change in VIM expression due to soup correction")
  p.DLX6.AS1 <- plotChangeMap(sc, out, "DLX6-AS1") + ggtitle("Change in DLX6-AS1 expression due to soup correction")

  plgr <- plot_grid(plotlist = list(p.DCX, p.SATB2, p.VIM, p.DLX6.AS1), nrow = 2)
  save_plot(filename =pps(OutDirDecont, "Change.in.gene.expression.due.to.soup.correction.png")
              , plot = plgr, base_height = wA4, base_width = wA4)

  # Write out ------------------------
  WriteMtxLocal = T
  OutDirDecont_local = "~/Data/Soup.stats"
  if (WriteMtxLocal) {
    DropletUtils:::write10xCounts(path = OutDirDecont_local, x = out, overwrite = TRUE)
    system(paste("gzip", pps(OutDirLocal,"barcodes.tsv.gz")),  wait = FALSE) # execute in the background
    system(paste("gzip", pps(OutDirLocal,"features.tsv.gz")),  wait = FALSE) # execute in the background
    system(paste("gzip", pps(OutDirLocal,"matrix.mtx.gz")),  wait = FALSE) # execute in the background
  } else {
    DropletUtils:::write10xCounts(path = OutDirDecont, x = out, overwrite = TRUE)
    # system(paste("gzip", pathDecontMtx),  wait = FALSE) # execute in the background
  }

  # plotMarkerDistribution(sc); wplot_save_this(plotname =  "SoupX.Auto.MarkerDistribution" )
  say()
}


# plotTheSoup ------------------------------------------------------------------------------------------------
if (F) {

  for (i in 1:l(DataSets)) {
    DataSet <- names(v.parentfolder)[i]
    # Setup ------------------------
    InDir <- v.parentfolder[i]; print(InDir)
    DataSet <- names(v.parentfolder)[i]

    # Plotting ------------------------
    plotTheSoup(CellRangerOutputDir = InDir, SeqRun = DataSet)
    system(paste("say Soup Plot ",i,"Ready"))
  }

}


# zacc ----------------------------------------------------

# # DataSets <- c('123062', '123063')
# DataSets <- c('123062', '123063', '101146', '101147')
# DataSets <- c('101146', '101147')
