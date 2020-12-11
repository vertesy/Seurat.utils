######################################################################
# SoupX.Marker.Analysis.R
######################################################################
# source('SoupX.Marker.Analysis.R')
# rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)

# Functions ------------------------
try (source('/GitHub/Packages/CodeAndRoll/CodeAndRoll.R'),silent= F)
try(require('MarkdownReportsDev'), silent = T)
# source('~/Github/TheCorvinas/R/DatabaseLinke.r')

# Setup ------------------------
setup_MarkdownReports(OutDir = "~/Data/POL.GFP/OLD/mtx.broken/SoupX.Marker.Analysis", scriptname = 'SoupX.Marker.Analysis.R')

# Metadata ------------------------

# Parameters ------------------------


# Read In ------------------------


SoupX.Marker.Analysis = FALSE
if (SoupX.Marker.Analysis) {

  "Required POL object to be loaded"
    MarkersFound.123062 <- c('NEAT1', 'FBXL7', 'ZBTB20', 'DACH1', 'AKAP12', 'SEZ6L', 'SLC1A3', 'ZEB1', 'EGR1', 'RORA', 'SAT1'
                             , 'SVIL', 'FOS', 'GMDS', 'DBI', 'IFITM3', 'LGALS1', 'FABP5', 'CDH4', 'AL589740.1')
    MarkersFound.123063 <- c('LGALS1', 'IFITM3', 'SPARC', 'TFF3', 'TIMP1', 'FN1', 'IGFBP4', 'ZFP36L1', 'DCN', 'LTBP1'
                             , 'ANXA2', 'IGFBP3', 'WWTR1', 'ID4', 'CDH11', 'S100A10', 'S100A11', 'MT2A', 'BNC2', 'ID3' )

# ------------------------------------------------------------------------------------------------
  MarkersFound.SoupX.QuickMarkers <- list(
    MarkersFound.123062 = MarkersFound.123062,
    MarkersFound.123063 = MarkersFound.123063
  )
  wvenn(MarkersFound.SoupX.QuickMarkers)


  multiFeaturePlot.A4(list.of.genes = MarkersFound.123062, obj = combined.obj, subdir =T)
  multiFeaturePlot.A4(list.of.genes = MarkersFound.123063, obj = combined.obj, subdir =T)


}

# QC ------------------------

# ------------------------
# ------------------------
# ------------------------
# ------------------------
# ------------------------


