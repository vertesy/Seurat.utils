# ____________________________________________________________________
# Seurat.utils.less.used.R ----
# ____________________________________________________________________
# source("~/GitHub/Packages/Seurat.utils/R/Seurat.utils.less.used.R")



# ____________________________________________________________________________________
#' @title plot.UMAP.tSNE.sidebyside
#'
#' @description Plot a UMAP and tSNE side by side.
#' @param obj Seurat object. Default: combined.obj
#' @param grouping Variable to group cells by. Default: 'res.0.6'
#' @param no_legend Logical, whether to display legend. Default: FALSE
#' @param do_return Logical, whether to return plot object. Default: TRUE
#' @param do_label Logical, whether to display labels. Default: TRUE
#' @param label_size Size of labels. Default: 10
#' @param vector_friendly Logical, whether to optimize for vector outputs. Default: TRUE
#' @param cells_use A vector of cell names to use for the plot. Default: NULL
#' @param no_axes Logical, whether to hide axes. Default: TRUE
#' @param pt_size Size of points. Default: 0.5
#' @param name.suffix Suffix to append to the plot's name. Default: NULL
#' @param width Width of the plot. Default: hA4
#' @param heigth Height of the plot. Default: 1.75 * wA4
#' @param filetype Filetype to save plot as. Default: 'pdf'
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @seealso
#'  \code{\link[cowplot]{save_plot}}
#' @importFrom cowplot save_plot plot_grid
#'
#' @export plot.UMAP.tSNE.sidebyside
plot.UMAP.tSNE.sidebyside <- function(obj = combined.obj, grouping = "res.0.6", # Plot a UMAP and tSNE sidebyside
                                      no_legend = FALSE,
                                      do_return = TRUE,
                                      do_label = TRUE,
                                      label_size = 10,
                                      vector_friendly = TRUE,
                                      cells_use = NULL,
                                      no_axes = TRUE,
                                      pt_size = 0.5,
                                      name.suffix = NULL,
                                      width = hA4, heigth = 1.75 * wA4, filetype = "pdf", ...) {
  p1 <- Seurat::DimPlot(
    object = obj, reduction.use = "tsne", no.axes = no_axes, cells.use = cells_use,
    no.legend = no_legend, do.return = do_return, do.label = do_label, label.size = label_size,
    group.by = grouping, vector.friendly = vector_friendly, pt.size = pt_size, ...
  ) +
    ggtitle("tSNE") + theme(plot.title = element_text(hjust = 0.5))

  p2 <- Seurat::DimPlot(
    object = obj, reduction.use = "umap", no.axes = no_axes, cells.use = cells_use,
    no.legend = TRUE, do.return = do_return, do.label = do_label, label.size = label_size,
    group.by = grouping, vector.friendly = vector_friendly, pt.size = pt_size, ...
  ) +
    ggtitle("UMAP") + theme(plot.title = element_text(hjust = 0.5))

  plots <- cowplot::plot_grid(p1, p2, labels = c("A", "B"), ncol = 2)
  plotname <- kpp("UMAP.tSNE", grouping, name.suffix, filetype)

  cowplot::save_plot(
    filename = plotname, plot = plots,
    ncol = 2 # we're saving a grid plot of 2 columns
    , nrow = 1 # and 2 rows
    , base_width = width,
    base_height = heigth
    # each individual subplot should have an aspect ratio of 1.3
    # , base_aspect_ratio = 1.5
  )
}

# _________________________________________________________________________________________________
#' @title Plot and Save UMAP without legend
#'
#' @description Generates a UMAP plot colored by a specified metadata column and saves the plot to a file.
#'
#' @param obj Seurat object to be visualized; Default: `combined.obj`.
#' @param metaD.colname Metadata column name to color the UMAP by; Default: 'metaD.colname.labeled'.
#' @param ext File extension for the saved plot, supports 'png', 'pdf', etc.; Default: 'png'.
#' @param ... Additional arguments passed to Seurat's `DimPlot`.
#'
#' @return Displays a UMAP plot and saves it to the current working directory.
#'
#' @examples
#' \dontrun{
#'   if (interactive()) {
#'     umapNamedClusters(obj = combined.obj, metaD.colname = "metaD.colname.labeled")
#'   }
#' }
#'
#' @export
#' @importFrom Seurat DimPlot
#' @importFrom ggplot2 ggtitle
#' @importFrom cowplot save_plot

umapNamedClusters <- function(obj = combined.obj,
                              metaD.colname = metaD.colname.labeled,
                              ext = "png", ...) {
  warning("This function is deprecated. No support.")
  fname <- ppp("Named.clusters", metaD.colname, ext)
  p.named <-
    Seurat::DimPlot(obj, reduction = "umap", group.by = metaD.colname, label = TRUE, ...) +
    NoLegend() +
    ggtitle(metaD.colname)
  save_plot(p.named, filename = fname)
  p.named
}



# _________________________________________________________________________________________________
#' @title Convert10Xfolders.old
#'
#' @description This function takes a parent directory with a number of subfolders, each containing the standard output of 10X Cell Ranger. It (1) loads the filtered data matrices, (2) converts them to Seurat objects, and (3) saves them as .RDS files.
#' @param InputDir A character string specifying the input directory.
#' @param folderPattern A character vector specifying the pattern of folder names to be searched. Default is 'filtered'.
#' @param min.cells An integer value specifying the minimum number of cells. Default is 10.
#' @param min.features An integer value specifying the minimum number of features. Default is 200.
#' @param updateHGNC A logical value indicating whether to update the HGNC. Default is TRUE.
#' @param ShowStats A logical value indicating whether to show statistics. Default is TRUE.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   Convert10Xfolders.old(InputDir = InputDir)
#' }
#' }
#' @export
Convert10Xfolders.old <- function(
    InputDir,
    folderPattern = c("filtered", "SoupX_decont")[1],
    min.cells = 10, min.features = 200,
    updateHGNC = TRUE, ShowStats = TRUE) {
  # ... function body ...
}

Convert10Xfolders.old <- function(
    InputDir # Take a parent directory with a number of subfolders, each containing the standard output of 10X Cell Ranger. (1.) It loads the filtered data matrices; (2.) converts them to Seurat objects, and (3.) saves them as *.RDS files.
    , folderPattern = c("filtered", "SoupX_decont")[1],
    min.cells = 10, min.features = 200, updateHGNC = TRUE, ShowStats = TRUE) {
  fin <- list.dirs(InputDir, recursive = FALSE)
  fin <- CodeAndRoll2::grepv(x = fin, pattern = folderPattern, perl = FALSE)

  for (i in 1:length(fin)) {
    pathIN <- fin[i]
    print(pathIN)
    fnameIN <- basename(fin[i])
    fnameOUT <- ppp(paste0(InputDir, "/", fnameIN), "min.cells", min.cells, "min.features", min.features, "Rds")
    count_matrix <- Read10X(pathIN)

    if (!is.list(count_matrix) | length(count_matrix) == 1) {
      seu <- CreateSeuratObject(
        counts = count_matrix, project = fnameIN,
        min.cells = min.cells, min.features = min.features
      )
    } else if (is.list(count_matrix) & length(count_matrix) == 2) {
      seu <- CreateSeuratObject(
        counts = count_matrix[[1]], project = fnameIN,
        min.cells = min.cells, min.features = min.features
      )

      # LSB, Lipid Sample barcode (Multi-seq) --- --- --- --- --- ---
      LSB <- CreateSeuratObject(counts = count_matrix[[2]], project = fnameIN)
      LSBnameOUT <- ppp(paste0(InputDir, "/LSB.", fnameIN), "Rds")
      saveRDS(LSB, file = LSBnameOUT)
    } else {
      print("More than 2 elements in the list of matrices")
    }
    # update --- --- --- ---
    if (updateHGNC) seu <- UpdateGenesSeurat(seu, EnforceUnique = TRUE, ShowStats = TRUE)
    saveRDS(seu, file = fnameOUT)
  }
}


# _________________________________________________________________________________________________
#' @title set.all.genes
#'
#' @description It is just a reminder to use calc.q99.Expression.and.set.all.genes to create the all.genes variable
#' @param obj Seurat object, Default: combined.obj
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   set.all.genes()
#'   all.genes
#' }
#' }
#' @export
set.all.genes <- function(obj = combined.obj) iprint("Use calc.q99.Expression.and.set.all.genes()")

