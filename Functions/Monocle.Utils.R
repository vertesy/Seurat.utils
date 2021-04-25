######################################################################
# Monocle.Utils.R
######################################################################
# source('~/GitHub/Packages/Seurat.utils/Functions/Monocle.Utils.R')
# rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)

# Functions ------------------------
# source('https://raw.githubusercontent.com/vertesy/Seurat.Pipeline/main/elements/Load.packages.local.R')
# try(source("~/GitHub/Packages/Seurat.multicore/00.Load.Seurat3.Multicore.LOCAL.R"));


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
  gene <- intersect(rownames(obj), gene)
  if (l(gene)) {
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
  } else { iprint("gene not found") }
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

  ls.genes <- intersect(rownames(obj), ls.genes)
  iprint(l(ls.genes), "genes found.")
  create_set_SubDir(ppp("3D.gex.plots", substitute(obj)))
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



# m3.get.umap ------------------------
m3.get.umap <- function(obj = cds_from_seurat, slot = 'UMAP', dim = (2:3)[2]) {
  iprint(dim, 'dimensional', slot)
  df.reduc <- obj@int_colData@listData$reducedDims[[slot]]
  stopif(condition = is_null(df.reduc), message = paste('slot not found', slot))
  stopifnot(ncol(df.reduc) == dim)
  return(df.reduc)
}
# m3.get.umap(obj = cds_from_seurat, slot = 'UMAP', dim = (2:3)[2])

# ------------------------
m3.backup.umap <- function(obj = cds_from_seurat, slot = 'UMAP', dim = (2:3)[2], ...) {
  new.slot <- p0(slot,'.',dim, 'D')
  obj@int_colData@listData$reducedDims[[new.slot]] <- m3.get.umap(obj = obj, slot = slot, dim = dim, ... )
  iprint('obj@int_colData@listData$reducedDims$', new.slot)
  return(obj)
}
# cds_from_seurat <- m3.backup.umap(obj = cds_from_seurat, slot = 'UMAP', dim = (2:3)[2])


# ------------------------
m3.recall.umap <- function(obj = cds_from_seurat, slot = 'UMAP', dim = (2:3)[2], ...) {
  backup.slot <- p0(slot,'.',dim, 'D')
  old.dim <- ncol(obj@int_colData@listData$reducedDims[[slot]])
  obj@int_colData@listData$reducedDims[[slot]] <- obj@int_colData@listData$reducedDims[[backup.slot]]
  iprint(old.dim, 'dimensional', slot, 'replaced by', dim, slot, 'reduction.')
  return(obj)
}
cds_from_seurat <- m3.recall.umap(obj = cds_from_seurat, slot = 'UMAP', dim = (2:3)[2])
# cds_from_seurat.bac <- cds_from_seurat


# ------------------------
m3.export.umap.2.Seurat <- function(mobj = cds_from_seurat, sobj = combined.obj, def.dim = 2) {

  reduc_bac <- sobj@misc$reductions.backup
  reduc_bac$'umap2d.Seurat' <- reduc_bac$'umap2d'
  reduc_bac$'umap2d'@cell.embeddings <- m3.get.umap(obj = cds_from_seurat, slot = 'UMAP.2D', dim = 2)
  colnames(reduc_bac$'umap2d'@cell.embeddings) <- paste('UMAP', 1:2, sep = '_')

  reduc_bac$'umap3d.Seurat' <- reduc_bac$'umap3d'
  reduc_bac$'umap3d'@cell.embeddings <- m3.get.umap(obj = cds_from_seurat, slot = 'UMAP.3D', dim = 3)
  colnames(reduc_bac$'umap3d'@cell.embeddings) <- paste('UMAP', 1:3, sep = '_')

  sobj@misc$reductions.backup <- reduc_bac
  sobj@reductions$umap <-  reduc_bac$'umap2d'  # if (def.dim == 2) { reduc_bac$'umap2d' } else if (def.dim == 3) { reduc_bac$'umap3d' }

  return(sobj)
}
# combined.obj.bac <- combined.obj
# combined.obj <- combined.obj.bac
combined.obj <- m3.export.umap.2.Seurat(mobj = cds_from_seurat, sobj = combined.obj); qUMAP("DDIT4")

# ------------------------






