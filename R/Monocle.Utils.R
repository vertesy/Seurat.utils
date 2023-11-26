# _________________________________________________________________________________________________
# Monocle.Utils.R ______________________________ ----
# _________________________________________________________________________________________________
# source("~/GitHub/Packages/Seurat.utils/R/Monocle.Utils.R")


# _________________________________________________________________________________________________
# _________________________________________________________________________________________________


# ________________________________________________________________________________
#' @title mplotGene
#'
#' @description Plots a specific gene in Monocle using a chosen dimensionality reduction method.
#' @param gene A string representing the gene of interest. Default: 'PGK1'.
#' @param reduction A string representing the dimensionality reduction method to use. Can be 'UMAP', 'tSNE', or 'PCA'. Default: 'UMAP'.
#' @param obj A Seurat object that will be plotted. Default: cds_from_seurat.
#' @return A plot showing the gene of interest.
#' @details This function takes a specific gene and plots it in Monocle using a specified dimensionality reduction method. The plot is then saved using the ggExpress::qqSave function.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   mplotGene()
#'   mplotGene("GAPDH")
#' }
#' }
#' @seealso
#' \code{\link[plot_cells]{plot_cells}}, \code{\link[qqSave]{qqSave}}
#' @export
mplotGene <- function(gene = "PGK1", reduction = "UMAP", obj = cds_from_seurat) {
  pl1 <- plot_cells(
    cds = obj,
    # color_cells_by = 'clusters_low',
    genes = gene,
    reduction_method = reduction,
    # color_cells_by = 'clusters_superlow',
    label_groups_by_cluster = F,
    label_leaves = F,
    label_branch_points = F, label_cell_groups = F,
    group_label_size = 10,
    cell_size = 1, alpha = .5
  )
  ggExpress::qqSave(pl1, fname = ppp(reduction, gene, ".png"), w = 14, h = 7)
}



# ________________________________________________________________________________
#' @title mplotManyGenes
#'
#' @description Plot many genes in Monocle.
#' @param ls.genes PARAM_DESCRIPTION, Default: c(`S-phase` = "TOP2A",
#'    `G2M-phase` = "HIST1H4C", oRG = "ID4",
#'    oRG = "HOPX", `Intermediate progenitor` = "EOMES", `Intermediate progenitor1` = "TAC3",
#'    Astroglia = "GFAP", Astrocyte = "S100B", `Immature neurons` = "SLA",
#'    Interneurons = "DLX6-AS1", `Hypoxia/Stress` = "DDIT4", Glycolytic = "PDK1",
#'    `Low-Quality` = "POLR2A", PGC = "DCN", `dl-EN` = "KAZN",
#'    `ul-EN` = "SATB2")
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   mplotManyGenes()
#' }
#' }
#' @export
mplotManyGenes <- function(ls.genes = c(
                             `S-phase` = "TOP2A", `G2M-phase` = "HIST1H4C",
                             `oRG` = "ID4", `oRG` = "HOPX",
                             `Intermediate progenitor` = "EOMES", `Intermediate progenitor1` = "TAC3",
                             Astroglia = "GFAP", Astrocyte = "S100B",
                             `Immature neurons` = "SLA", Interneurons = "DLX6-AS1",
                             `Hypoxia/Stress` = "DDIT4", Glycolytic = "PDK1",
                             `Low-Quality` = "POLR2A", `PGC` = "DCN",
                             `dl-EN` = "KAZN", `ul-EN` = "SATB2"
                           )) {
  create_set_SubDir("gene.expression")
  for (g in ls.genes) {
    mplotGene(g)
  }
  create_set_OutDir(ParentDir)
}



# ________________________________________________________________________________
#' @title m3DplotGene
#'
#' @description Plot a gene in 3D in Monocle.
#' @param gene gene of interest, Default: 'PGK1'
#' @param reduction UMAP, tSNE, or PCA (Dim. reduction to use), Default: 'UMAP'
#' @param obj Seurat object, Default: cds.10pc
#' @param ttl.suffix PARAM_DESCRIPTION, Default: 'expression'
#' @param suffix A suffix added to the filename, Default: ''
#' @param cex Point size, Default: 5
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   m3DplotGene("DDIT4")
#' }
#' }
#' @export
m3DplotGene <- function(
    gene = "PGK1", reduction = "UMAP", obj = cds.10pc, ttl.suffix = "expression",
    suffix = "", cex = 5) {
  gene <- intersect(rownames(obj), gene)
  if (length(gene)) {
    pl1 <- plot_cells_3d(
      cds = obj
      # color_cells_by = 'clusters_low',
      , genes = gene,
      reduction_method = reduction
      # color_cells_by = 'clusters_superlow',
      # label_groups_by_cluster = F,
      # label_leaves = F,
      # label_branch_points = F, label_cell_groups = F
      # , group_label_size = 10
      , cell_size = cex,
      alpha = .5
    ) %>%
      plotly::layout(title = paste0(gene, ttl.suffix))
    SavePlotlyAsHtml(pl1, category. = gene, suffix. = suffix)
  } else {
    iprint("gene not found")
  }
}



# ________________________________________________________________________________
#' @title m3DplotKeyGenes
#'
#' @description Plot many genes in 3D in Monocle.
#' @param obj Seurat object, Default: cds.10pc
#' @param cex Point size, Default: iround(log10(idim(obj)[2]))
#' @param ls.genes PARAM_DESCRIPTION, Default: c(`S-phase` = "TOP2A", `G2M-phase` = "HIST1H4C", oRG = "ID4",
#'    oRG = "HOPX", `Intermediate progenitor` = "EOMES", `Intermediate progenitor1` = "TAC3",
#'    Astroglia = "GFAP", Astrocyte = "S100B", `Immature neurons` = "SLA",
#'    Interneurons = "DLX6-AS1", `Hypoxia/Stress` = "DDIT4", Glycolytic = "PDK1",
#'    `Low-Quality` = "POLR2A", PGC = "DCN", `dl-EN` = "KAZN",
#'    `ul-EN` = "SATB2")
#' @param reduction UMAP, tSNE, or PCA (Dim. reduction to use), Default: 'UMAP'
#' @param suffix A suffix added to the filename, Default: ''
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   m3DplotKeyGenes()
#' }
#' }
#' @export
m3DplotKeyGenes <- function(
    obj = cds.10pc, cex = iround(log10(idim(obj)[2])),
    ls.genes = c(
      `S-phase` = "TOP2A", `G2M-phase` = "HIST1H4C",
      `oRG` = "ID4", `oRG` = "HOPX",
      `Intermediate progenitor` = "EOMES",
      `Intermediate progenitor1` = "TAC3",
      Astroglia = "GFAP", Astrocyte = "S100B",
      `Immature neurons` = "SLA", Interneurons = "DLX6-AS1",
      `Hypoxia/Stress` = "DDIT4", Glycolytic = "PDK1",
      `Low-Quality` = "POLR2A", `PGC` = "DCN",
      `dl-EN` = "KAZN", `ul-EN` = "SATB2"
    ),
    reduction = "UMAP", suffix = "") {
  ls.genes <- intersect(rownames(obj), ls.genes)
  iprint(length(ls.genes), "genes found.")
  create_set_SubDir(ppp("3D.gex.plots", substitute(obj)))
  for (g in ls.genes) {
    m3DplotGene(gene = g, obj = obj, cex = cex)
  }
  create_set_OutDir(ParentDir)
}


# subsetMonocleObject ------------------------
#' @title subsetMonocleObject
#'
#' @description Subsets a compressed Seurat object and saves it in the working directory.
#' @param obj A Seurat object to be subsetted. Default: cds_from_seurat.
#' @param fraction_ A numeric value indicating the fraction of the object to be kept. Default: 0.1.
#' @param nCells A numeric or boolean value indicating the number of cells to keep. If FALSE, the function calculates the number of cells to keep based on the fraction parameter. Default: FALSE.
#' @param seed_ A numeric value specifying the seed for random number generation. Default: 1989.
#' @return A subset of the original Seurat object.
#' @details This function subsets a Seurat object by either a specified fraction or a specific number of cells. The subset object is then returned by the function.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   cds.10pc <- subsetMonocleObject(cds_from_seurat, fraction_ = .1)
#'   cds.25pc <- subsetMonocleObject(cds_from_seurat, fraction_ = .25)
#' }
#' }
#' @seealso
#' \code{\link[sample]{sample}}, \code{\link[percentage_formatter]{percentage_formatter}}
#' @export
#' @importFrom Stringendo percentage_formatter

subsetMonocleObject <- function(obj = cds_from_seurat, fraction_ = 0.1, nCells = F, seed_ = 1989) { # Subset a compressed Seurat Obj and save it in wd.
  set.seed(seed_)
  if (isFALSE(nCells)) {
    all.cells <- colnames(obj)
    n.keep <- floor(length(all.cells) * fraction_)
    cellIDs.keep <- sample(all.cells, size = n.keep, replace = FALSE)
    iprint(length(cellIDs.keep), "or", Stringendo::percentage_formatter(fraction_), "of the cells are kept. Seed:", head(cellIDs.keep), seed_)
    # cellIDs.keep
  } else if (nCells > 1) {
    nKeep <- min(ncol(obj), nCells)
    # print(nKeep)
    cellIDs.keep <- sample(colnames(obj), size = nKeep, replace = F)
    if (nKeep < nCells) iprint("Only", nCells, "cells were found in the object, so downsampling is not possible.")
  }
  obj <- obj[, cellIDs.keep] # downsample
  return(obj)
}




# m3.get.umap ------------------------
#' @title m3.get.umap
#' @description Fetch the umap coordinates from obj@int_colData@listData$reducedDims[[slot]]
#' @param obj Seurat object, Default: cds_from_seurat
#' @param slot slot in the Seurat object. Default: 'UMAP'
#' @param dim Numer of dimensions used, Default: (2:3)[2]
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   m3.get.umap(obj = cds_from_seurat, slot = "UMAP", dim = (2:3)[2])
#' }
#' }
#' @export
m3.get.umap <- function(obj = cds_from_seurat, slot = "UMAP", dim = (2:3)[2]) {
  iprint(dim, "dimensional", slot)
  df.reduction <- obj@int_colData@listData$reducedDims[[slot]]
  stopif(condition = is.null(df.reduction), message = paste("slot not found", slot))
  stopifnot(ncol(df.reduction) == dim)
  return(df.reduction)
}


# _________________________________________________________________________________________________
#' @title m3.backup.umap
#' @description Backup umap coordinates to obj@int_colData@listData$reducedDims[[new.slot]]
#' @param obj Seurat object, Default: cds_from_seurat
#' @param slot slot in the Seurat object. Default: 'UMAP'
#' @param dim Numer of dimensions used, Default: (2:3)[2]
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   cds_from_seurat <- m3.backup.umap(obj = cds_from_seurat, slot = "UMAP", dim = (2:3)[2])
#' }
#' }
#' @export
m3.backup.umap <- function(obj = cds_from_seurat, slot = "UMAP", dim = (2:3)[2], ...) {
  new.slot <- paste0(slot, ".", dim, "D")
  obj@int_colData@listData$reducedDims[[new.slot]] <- m3.get.umap(obj = obj, slot = slot, dim = dim, ...)
  iprint("obj@int_colData@listData$reducedDims$", new.slot)
  return(obj)
}



# _________________________________________________________________________________________________
#' @title m3.recall.umap
#'
#' @description Fetch UMAP coordinates.
#' @param obj Seurat object, Default: cds_from_seurat
#' @param slot slot in the Seurat object. Default: 'UMAP'
#' @param dim Numer of dimensions used, Default: (2:3)[2]
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   cds_from_seurat <- m3.recall.umap(obj = cds_from_seurat, slot = "UMAP", dim = (2:3)[2])
#'   cds_from_seurat.bac <- cds_from_seurat
#' }
#' }
#' @export
m3.recall.umap <- function(obj = cds_from_seurat, slot = "UMAP", dim = (2:3)[2], ...) {
  backup.slot <- paste0(slot, ".", dim, "D")
  old.dim <- ncol(obj@int_colData@listData$reducedDims[[slot]])
  obj@int_colData@listData$reducedDims[[slot]] <- obj@int_colData@listData$reducedDims[[backup.slot]]
  iprint(old.dim, "dimensional", slot, "replaced by", dim, slot, "reduction.")
  return(obj)
}



# _________________________________________________________________________________________________
#' @title m3.export.umap.2.Seurat
#' @description Exports UMAP coordinates from a monocle object to a Seurat object.
#' @param mobj The Monocle object from which to extract UMAP coordinates. Default: cds_from_seurat.
#' @param sobj The Seurat object to which UMAP coordinates will be added. Default: combined.obj.
#' @param def.dim The dimension of the UMAP coordinates to extract (2D or 3D). Default: 2.
#' @return A Seurat object with the UMAP coordinates added.
#' @details This function extracts UMAP coordinates from a Monocle object and adds them to a Seurat object. The dimensions of the UMAP coordinates can be specified.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   combined.obj.bac <- combined.obj # combined.obj <- combined.obj.bac
#'   combined.obj <- m3.export.umap.2.Seurat(mobj = cds_from_seurat, sobj = combined.obj)
#'   qUMAP("DDIT4")
#' }
#' }
#' @seealso
#' \code{\link[m3.get.umap]{m3.get.umap}}
#' @export
m3.export.umap.2.Seurat <- function(mobj = cds_from_seurat, sobj = combined.obj, def.dim = 2) {
  # ... function body ...
}

m3.export.umap.2.Seurat <- function(mobj = cds_from_seurat, sobj = combined.obj, def.dim = 2) {
  reduc_bac <- sobj@misc$reductions.backup
  reduc_bac$"umap2d.Seurat" <- reduc_bac$"umap2d"
  reduc_bac$"umap2d"@cell.embeddings <- m3.get.umap(obj = cds_from_seurat, slot = "UMAP.2D", dim = 2)
  colnames(reduc_bac$"umap2d"@cell.embeddings) <- paste("UMAP", 1:2, sep = "_")

  reduc_bac$"umap3d.Seurat" <- reduc_bac$"umap3d"
  reduc_bac$"umap3d"@cell.embeddings <- m3.get.umap(obj = cds_from_seurat, slot = "UMAP.3D", dim = 3)
  colnames(reduc_bac$"umap3d"@cell.embeddings) <- paste("UMAP", 1:3, sep = "_")

  sobj@misc$reductions.backup <- reduc_bac
  sobj@reductions$umap <- reduc_bac$"umap2d" # if (def.dim == 2) { reduc_bac$'umap2d' } else if (def.dim == 3) { reduc_bac$'umap3d' }

  return(sobj)
}

# _________________________________________________________________________________________________
