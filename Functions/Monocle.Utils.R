######################################################################
# Monocle.Utils.R
######################################################################
# source('~/GitHub/Packages/Seurat.utils/Functions/Monocle.Utils.R')
# rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)

# Functions ------------------------
try(source('~/GitHub/Packages/CodeAndRoll/CodeAndRoll.R'),silent= F)
require('MarkdownReportsDev')
# source('~/Github/TheCorvinas/R/DatabaseLinke.r')

# source('https://raw.githubusercontent.com/vertesy/Seurat.Pipeline/main/elements/Load.packages.local.R')
# try(source("~/GitHub/Packages/Seurat.multicore/00.Load.Seurat3.Multicore.LOCAL.R"));

# Setup ------------------------
OutDir = "~/Dropbox/Abel.IMBA/AnalysisD/"
setup_MarkdownReports(OutDir = OutDir, scriptname = "")
OutDirOrig = OutDir

# Metadata ------------------------

# Parameters ------------------------


# ------------------------
# ------------------------
# ------------------------
# ------------------------

# mplotGene ------------------------
mplotGene <- function(gene = "PGK1", reduc = "UMAP", obj = cds_from_seurat) {
  pl1 <- plot_cells(cds = obj,
                    # color_cells_by = 'clusters_low',
                    genes = gene,
                    reduction_method =reduc,
                    # color_cells_by = 'clusters_superlow',
                    label_groups_by_cluster = F,
                    label_leaves = F,
                    label_branch_points = F, label_cell_groups = F
                    , group_label_size = 10
                    , cell_size = 1, alpha = .5)
  qqSave(pl1, fname = ppp(reduc, gene, ".png"), w = 14, h=7)
}
# mplotGene()
# mplotGene("GAPDH")


# mplotManyGenes ------------------------
mplotManyGenes <- function(ls.genes = c(
  `S-phase` = "TOP2A", `G2M-phase` = "HIST1H4C"
  , `oRG` = "ID4", `oRG` = "HOPX"
  , `Intermediate progenitor` = "EOMES",  `Intermediate progenitor1` = "TAC3"
  , Astroglia = "GFAP", Astrocyte = "S100B"
  , `Immature neurons` = "SLA", Interneurons = "DLX6-AS1"
  , `Hypoxia/Stress` = "DDIT4", Glycolytic = "PDK1"
  , `Low-Quality` = "POLR2A", `PGC` = "DCN"
  , `dl-EN` = "KAZN", `ul-EN` = "SATB2")
) {
  create_set_SubDir("gene.expression")
  for (g in ls.genes) {
    mplotGene(g)
  }
  create_set_OutDir(ParentDir)
}
# mplotManyGenes()

# m3DplotGene ------------------------
m3DplotGene <- function(gene = "PGK1", reduc = "UMAP", obj = cds.10pc, ttl.suffix = "expression", suffix = "", cex = 5) {
  pl1 <- plot_cells_3d(cds = obj
                       # color_cells_by = 'clusters_low',
                       , genes = gene
                       , reduction_method = reduc
                       # color_cells_by = 'clusters_superlow',
                       # label_groups_by_cluster = F,
                       # label_leaves = F,
                       # label_branch_points = F, label_cell_groups = F
                       # , group_label_size = 10
                       , cell_size = cex
                       , alpha = .5) %>% layout(title=p0(gene, ttl.suffix))
  SavePlotlyAsHtml(pl1, category. = gene, suffix. = suffix)
}
m3DplotGene("DDIT4")


# m3DplotKeyGenes ------------------------
m3DplotKeyGenes <- function(obj = cds.10pc, cex = iround(log10(idim(obj)[2]))
                            , ls.genes =  c(
                                    `S-phase` = "TOP2A", `G2M-phase` = "HIST1H4C"
                                    , `oRG` = "ID4", `oRG` = "HOPX"
                                    , `Intermediate progenitor` = "EOMES",  `Intermediate progenitor1` = "TAC3"
                                    , Astroglia = "GFAP", Astrocyte = "S100B"
                                    , `Immature neurons` = "SLA", Interneurons = "DLX6-AS1"
                                    , `Hypoxia/Stress` = "DDIT4", Glycolytic = "PDK1"
                                    , `Low-Quality` = "POLR2A", `PGC` = "DCN"
                                    , `dl-EN` = "KAZN", `ul-EN` = "SATB2")
                            , reduc = "UMAP", suffix = "") {

  create_set_SubDir("3D.gex.plots")
  for (g in ls.genes) {
    m3DplotGene(gene = g, obj = obj, cex = cex)
  }
  create_set_OutDir(ParentDir)
}
# m3DplotKeyGenes()

# subsetMonocleObject ------------------------
subsetMonocleObject <-function(obj = cds_from_seurat, fraction_ = 0.1, nCells = F, seed_ = 1989 ) { # Subset a compressed Seurat Obj and save it in wd.
  set.seed(seed_)
  if (isFALSE(nCells)) {
    all.cells <- colnames(obj)
    n.keep <- floor(l(all.cells) * fraction_)
    cellIDs.keep <- sample(all.cells, size = n.keep, replace = FALSE)
    iprint(length(cellIDs.keep), "or",percentage_formatter(fraction_),"of the cells are kept. Seed:", head(cellIDs.keep), seed_)
    # cellIDs.keep
  } else if (nCells > 1) {
    nKeep = min(ncol(obj), nCells)
    # print(nKeep)
    cellIDs.keep = sample(colnames(obj), size = nKeep, replace = F)
    if (nKeep < nCells) iprint("Only",nCells,"cells were found in the object, so downsampling is not possible.")
  }
  obj <- obj[,cellIDs.keep] # downsample
  return(obj)
}
# cds.10pc <- subsetMonocleObject(cds_from_seurat, fraction_ = .1)
# cds.25pc <- subsetMonocleObject(cds_from_seurat, fraction_ = .25)



# ------------------------
# ------------------------
# ------------------------


