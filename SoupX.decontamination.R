######################################################################
# SoupX.decontamination.R
######################################################################
source("~/GitHub/Seurat.utils/SoupX.decontamination.R")


# Libraries ------------------------
require(SoupX)
require(Seurat)
require(MarkdownReportsDev)
try(source("~/GitHub/Seurat.utils/Seurat.Utils.Load.R"));


# Parameters ------------------------
v.parentfolder = c(
  # SEO -----
  # "/Volumes/copy.your.own.data.here/A.Vertesy/SEO/HNV73DRXX_R10015/HNV73DRXX_R10015/aligned_rna/124851_rnacount",
  # "/Volumes/single.cell.RNA.seq/C.Bosone/2020.06/HN3V3DRXX_R9836/aligned_rna/123062_rnacount",
  # # POL -----
  # "/Volumes/HN3V3DRXX_R9836/demultiplexed/123062/123062_premRNA_POL/outs/",
  "/Volumes/HN3V3DRXX_R9836/demultiplexed/123063/123063_premRNA_POL/outs/"
  # "/Users/abel.vertesy/Data/POL/filtered_feature_bc_matrix.123062"
  # "/Volumes/single.cell.RNA.seq/C.Bosone/2020.06/HN3V3DRXX_R9836/aligned_rna/123063_rnacount",
  # CON -----
  # "/Volumes/single.cell.RNA.seq/A.Vertesy/CONN/Connectome.pilot.1/03.Fulldepth/HKJHKDRXX_R9539/114593/114593_premRNA_local/outs/",
  # INM -----
  # # "/Volumes/single.cell.RNA.seq/S.Bajaj/HJGHFDRXX_all/aligned/104566",
  # # "/Volumes/single.cell.RNA.seq/S.Bajaj/HJGHFDRXX_all/aligned/104567",
  # # "/Volumes/single.cell.RNA.seq/S.Bajaj/HJGHFDRXX_all/aligned/104568",
  # "/Volumes/single.cell.RNA.seq/S.Bajaj/HJGHFDRXX_all/aligned/104569"
)

i=1
for (i in 1:l(v.parentfolder)) {

  # Setup ------------------------
  InDir <- v.parentfolder[i]
  OutDir <- OutDirDecont <- file.path(InDir,"SoupX_decont_filt_feature_bc_matrix")
  dir.create(OutDirDecont)

  # Estimate ------------------------
  sc = load10Xv3(InDir, verbose = T)
  sc = autoEstCont(sc, verbose = T);  wplot_save_this(plotname =  "SoupX.contamination.fraciton.rho" )
  out = adjustCounts(sc, verbose = T)

  # Write out ------------------------
  (pathDecontMtx <- file.path(OutDirDecont,"matrix.mtx"))
  Matrix::writeMM(obj = out, file=pathDecontMtx )
  # system(paste("gzip", pathDecontMtx),  wait = FALSE) # execute in the background
  plotMarkerDistribution(sc); wplot_save_this(plotname =  "SoupX.Auto.MarkerDistribution" )
  say()
}

strsplit(x = v.parentfolder,split =  '/')

DataSets <- c('123062', '123063')
i=2
for (i in 1:l(v.parentfolder)) {
  # Plotting ------------------------
  plotTheSoup(CellRangerOutputDir = InDir, SeqRun = DataSets[i])
  say()
}

