######################################################################
# Visualize Single-cell RNA-seq data using a Seurat object
######################################################################
# source('Visualize.Gene.expression.v2.R')
# rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)

"This is an example usage script for first time users of Seurat & my associated function libraries."

# 1. Install function libraries ------------------------
if (FALSE) {
  "Run these only once"
  install.packages(Seurat)
  install.packages(dplyr)
  install.packages(cowplot)
  install.packages(stringr)
  install.packages(tictoc)
  install.packages(tibble)
  install.packages(ggplot2)

  install.packages("devtools"); # If you don't have it
  require("devtools")
  devtools::install_github(repo = "vertesy/MarkdownReportsDev")
}

# 2. Load function libraries ------------------------
try (source('https://raw.githubusercontent.com/vertesy/CodeAndRoll/master/CodeAndRoll.R'),silent= F)
try(source('https://raw.githubusercontent.com/vertesy/Seurat.utils/master/00.Load.Seurat.Utils.WEB.R'), silent =   T)
try(source('https://raw.githubusercontent.com/vertesy/Seurat.multicore/master/00.Load.Seurat3.Multicore.R'), silent=T)
# If it fails, see alternative at the end

require('MarkdownReportsDev')
# If installation of MarkdownReportsDev did not work:
# try (source('https://raw.githubusercontent.com/vertesy/MarkdownReportsDev/master/R/MarkdownReportsDev.R'),silent= F)

require(Seurat)
require(dplyr)
require(cowplot)
require(stringr)
require(tictoc)
require(tibble)
require(ggplot2)


# 2b Setup ------------------------
"Where do you want to save the output files?"
OutDir = "~/Downloads/scAnalysis"
setup_MarkdownReports(OutDir = OutDir, scriptname = "Visualize.Gene.expression.v2.R")
OutDirOrig = OutDir


# 3. load Seurat object (sc data) -----------------------------------------------------------------
"Where did you save combined.obj.RDS.gz ?"
combined.obj <- readRDS("path/to/combined.obj.RDS.gz")

# 4. create list of genes (all.genes -----------------------------------------------------------------
recallAllGenes()

# 5. search for genes -----------------------------------------------------------------
"Start typing after $ and hit TAB for autocomplete"
all.genes$SATB2

# 5b. plot genes -----------------------------------------------------------------
qUMAP(feature = "SATB2", obj = combined.obj)

# 6. find.clustering -----------------------------------------------------------------
"Start typing after $ and hit TAB for autocomplete"
combined.obj$integrated_snn_res.0.5

# 6b. plot clusters -----------------------------------------------------------------
clUMAP(ident = "integrated_snn_res.0.5", obj = combined.obj)




# Custom gene sets  --------------------------------------------------------------------------------

SEGA.old = c("HES5", "ID3", "STMN1", "HMGN2", "SLC6A9", "VCAM1", "S1PR1",
         "NOTCH2", "S100A6", "BCAN", "NES", "ATP1A2", "PRDX6", "AGT",
         "ID2", "MYCN", "TMSB10", "DBI", "GAD1", "DLX1", "DLX2", "NEUROD1",
         "MAP2", "IGFBP5", "ACSL3", "PTMA", "SLC6A11", "SLC6A1", "RPL32",
         "RPSA", "SYNPR", "ALDH1L1", "MCM2", "HES1", "FGFR3", "PROM1",
         "HOPX", "SPARCL1", "BMPR1B", "H2AFZ", "CCNA2", "CPE", "SLC1A3",
         "THBS4", "LIX1", "NREP", "ACSL6", "SPARC", "NPM1", "HNRNPAB",
         "ID4", "TFEB", "NR2E1", "GJA1", "FABP7", "MOXD1", "DLL1", "LFNG",
         "SP8", "PPIA", "EGFR", "PTN", "RARRES2", "TMSB4X", "TSPAN7",
         "DCX", "BMP1", "CLU", "SFRP1", "PABPC1", "MYC", "NTRK2", "TNC",
         "RPL12", "CD81", "DKK3", "PAX6", "SLC1A2", "CCND1", "NCAM1",
         "PARD3", "CDK1", "HTRA1", "MKI67", "FOXM1", "CCND2", "CD9", "CDCA3",
         "ENO2", "SLCO1C1", "LDHB", "LIMA1", "CDK2", "RPL41", "HMGA2",
         "ASCL1", "RAN", "TPT1", "ITM2B", "CLDN10", "FOS", "CKB", "NUPR1",
         "TUBB3", "TP53", "AURKB", "ALDOC", "GFAP", "LUC7L3", "TOB1",
         "SOX9", "RBFOX3", "FASN", "MIB1", "AQP4", "PRNP", "PCNA", "JAG1",
         "CST3", "E2F1", "RPL18A", "FXYD1", "DLL3", "APOE", "SMARCB1",
         "MIAT", "OLIG2", "GSTM1", "GPR37L1", "BMP6", "CSPG4", "MT1A",
         "ETNPPL")

# Markers ----------------------------
ClassicMarkers = c(
  "Apical precursor" = 						"SOX2",
  "Stem cells" = 						  		"ID4",
  "Cycling cells" = 							"TOP2A",
  "IPC" = 												"EOMES",
  "IPC.late" = 		 		 		 		 	  "TAC3",

  "Astro-like" = 									"S100B",
  "Interneurons, (OPC, MSN)"   = 			"DLX6-AS1",
  "GABA-ergic inter neuron" = 		"GAD2",

  "General neuronal" = 						"NEUROD6",
  "CTIP2 " = 	  									"BCL11B",
  "Upper layer " =   										"SATB2",
  "Deep layer neurons" = 					"FEZF2",

  "Weird sc gene" = 										"MALAT1",
  "RPL34" = 										"RPL34",
  "Mito" =   										"MT-ATP6",
  "PDK" =   										"PDK1",
  "SLA" = 									    "SLA",
  "IGFBP5" = 										"IGFBP5"
)


# Alternative for SEURAT.UTILS - You can try to load each script individually ------------------------
if (FALSE) {
  "SEURAT.UTILS components"
  "You probably don't need all for your work"
  try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Seurat.update.gene.symbols.HGNC.R"))
  try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/metadata.manipulation.R"))
  try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/plotting.dim.reduction.2D.R"))
  try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/plotting.dim.reduction.3D.R"))
  try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/plotting.filtering.R"))
  try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/plotting.statistics.and.QC.R"))
  try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Read.Write.Save.Load.functions.R"))
  try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Seurat.object.manipulations.etc.R"))
  try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Cluster.Auto-naming.DE.R"))
  try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Seurat.gene.sets.and.GO.terms.R"))
  try (source('https://raw.githubusercontent.com/vertesy/Seurat.utils/master/MULTI-seq.functions.R'))

}

# Alternative for SEURAT.MULTICORE - You can try to load each script individually ------------------------
if (FALSE) {
  " SEURAT.MULTICORE components"
  "You probably don't need these at all for your work"
  try(source("https://raw.githubusercontent.com/vertesy/Seurat.multicore/master/Seurat3.Multicore.Read.Write.R"))
  try(source("https://raw.githubusercontent.com/vertesy/Seurat.multicore/master/Seurat3.Multicore.Generic.Functions.R"))

}


# Custom functions you need--------------------------------------------------------------------------------
if (FALSE) {
  "probably not needed"
  mmeta <- function(ColName.metadata = 'batch', obj = org, as_numeric =F) { # get a metadata column as a named vector
    x = as.named.vector(obj@meta.data[ ,ColName.metadata, drop=F])
    if (as_numeric) {
      as.numeric.wNames(x)+1
    } else {x}
  }

  jjpegA4 <- function(filename, r = 225, q = 90) {
    jpeg(file=filename,width=8.27, height=11.69, units = 'in', quality = q,res = r)
  }

}
