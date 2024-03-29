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
#' if (interactive()) {
#'   umapNamedClusters(obj = combined.obj, metaD.colname = "metaD.colname.labeled")
#' }
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
#' @title AutoNumber.by.PrinCurve
#'
#' @description Relabel cluster numbers along the principal curve of 2 UMAP (or tSNE) dimensions. #
#' @param obj Seurat object, Default: combined.obj
#' @param dim Dimensions to use, Default: 1:2
#' @param plotit Plot results (& show it), Default: TRUE
#' @param swap Swap Lambda paramter (multiplied with this) , Default: -1
#' @param reduction UMAP, tSNE, or PCA (Dim. reduction to use), Default: 'umap'
#' @param res Clustering resoluton to use, Default: 'integrated_snn_res.0.5'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   DimPlot.ClusterNames(ident = "integrated_snn_res.0.5")
#'   combined.obj <- AutoNumber.by.PrinCurve(
#'     obj = combined.obj, dim = 1:2, reduction = "umap", plotit = TRUE,
#'     swap = -1, res = "integrated_snn_res.0.5"
#'   )
#'   DimPlot.ClusterNames(ident = "integrated_snn_res.0.5.prin.curve")
#' }
#' }
#' @seealso
#'  \code{\link[princurve]{principal_curve}}
#' @importFrom princurve principal_curve whiskers
#' @importFrom MarkdownReports wplot_save_this
#' @importFrom Seurat FetchData
#'
#' @export
AutoNumber.by.PrinCurve <- function(
    obj = combined.obj # Relabel cluster numbers along the principal curve of 2 UMAP (or tSNE) dimensions.
    , dim = 1:2, plotit = TRUE, swap = -1,
    reduction = "umap", res = "integrated_snn_res.0.5") {
  # require(princurve)
  dim_name <- ppu(toupper(reduction), dim)
  coord.umap <- FetchData(object = obj, vars = dim_name)
  fit <- princurve::principal_curve(x = as.matrix(coord.umap))
  if (plotit) {
    plot(fit,
      xlim = range(coord.umap[, 1]), ylim = range(coord.umap[, 2]),
      main = "principal_curve"
    )
    # points(fit)
    points(coord.umap, pch = 18, cex = .25)
    princurve::whiskers(coord.umap, fit$s, lwd = .1)
    MarkdownReports::wplot_save_this(plotname = "principal_curve")
  }

  ls.perCl <- split(swap * fit$lambda, f = obj[[res]])
  MedianClusterCoordinate <- unlapply(ls.perCl, median)
  OldLabel <- names(sort(MedianClusterCoordinate))
  NewLabel <- as.character(0:(length(MedianClusterCoordinate) - 1))
  NewMeta <- translate(vec = obj[[res]], oldvalues = OldLabel, newvalues = NewLabel)
  NewMetaCol <- kpp(res, "prin.curve")
  iprint("NewMetaCol:", NewMetaCol)
  obj[[NewMetaCol]] <- NewMeta
  return(obj)
}



# _________________________________________________________________________________________________
# Read and Write Seurat Objects ----
# _________________________________________________________________________________________________



#' @title Load 10X Genomics Version 3 Data
#'
#' @description Loads 10X Genomics data from a specified directory containing output folders for raw and filtered data.
#' This function is designed to handle data from 10X Genomics Chromium Single Cell technologies (version 3).
#'
#' @param dataDir A string specifying the directory that contains the 10X Genomics output folders.
#' This directory should include subdirectories for raw and filtered data, typically named starting with
#' `raw_` and `filt_`, respectively.
#' @param cellIDs An optional vector of cell IDs to include in the loaded data. Default is `NULL`,
#' indicating that all available cells will be included. This is useful for subsetting the data based
#' on specific cell IDs.
#' @param channelName An optional string specifying the channel name for the data being loaded.
#' This can be used to label the data according to the experimental condition or sample name. Default is `NULL`.
#' @param readArgs A list of additional arguments to pass to the internal `Read10X` function used for
#' loading the data. Default is an empty list.
#' @param includeFeatures A character vector specifying which features to include in the loaded data.
#' Common values include "Gene Expression", "Antibody Capture", and "CRISPR Guide Capture".
#' Default is `c("Gene Expression")`.
#' @param verbose A logical flag indicating whether to print progress messages and status updates as the
#' data is loaded. Default is `TRUE`.
#' @param ... Additional arguments passed to other internally called functions, if applicable.
#' @return An object of class `SoupChannel`, representing the loaded 10X data. This object includes
#' raw counts, filtered counts, and optionally, additional metadata and dimensionality reduction coordinates
#' (e.g., t-SNE).
#'
#' @details This function provides a comprehensive approach to loading and organizing 10X Genomics data
#' for further analysis. It accommodates the data structure commonly found in 10X Genomics version 3 outputs
#' and allows for the inclusion of various types of molecular data as well as optional cell and channel
#' specifications.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Assuming `dataDir` is the path to your 10X data directory
#'   channel <- load10Xv3(dataDir = "path/to/10X/data")
#'   # Now `channel` contains the loaded 10X data as a `SoupChannel` object
#' }
#' }
#'
#' @seealso \code{\link[SoupX]{SoupChannel}} for the structure and utilities of the `SoupChannel` class.
#'
#' @export
#' @importFrom SoupX SoupChannel
load10Xv3 <- function(dataDir, cellIDs = NULL, channelName = NULL, readArgs = list(),
                      includeFeatures = c("Gene Expression"), verbose = TRUE,
                      ...) {
  # include
  dirz <- list.dirs(dataDir, full.names = FALSE, recursive = FALSE)
  path.raw <- file.path(dataDir, grep(x = dirz, pattern = "^raw_*", value = TRUE))
  path.filt <- file.path(dataDir, grep(x = dirz, pattern = "^filt_*", value = TRUE))
  CR.matrices <- list.fromNames(c("raw", "filt"))


  (isV3 <- any(grepl(x = dirz, pattern = "^raw_feature_bc*")))
  tgt <- path.raw

  if (!isV3) {
    tgt <- file.path(tgt, list.files(tgt))
  }
  if (verbose) {
    message(sprintf("Loading raw count data"))
  }
  dat <- do.call(Read10X, c(list(data.dir = tgt), readArgs))
  if (verbose) {
    message(sprintf("Loading cell-only count data"))
  }
  if (!is.null(cellIDs)) {
    if (all(grepl("\\-1$", cellIDs))) {
      cellIDs <- gsub("\\-1$", "", cellIDs)
    }
    if (!all(cellIDs %in% colnames(dat))) {
      stop("Not all supplied cellIDs found in raw data.")
    }
    datCells <- dat[, match(cellIDs, colnames(dat))]
  } else {
    tgt <- path.filt
    if (!isV3) {
      tgt <- file.path(tgt, list.files(tgt))
    }
    datCells <- do.call(Read10X, c(
      list(data.dir = tgt),
      readArgs
    ))
    if (is.list(dat)) {
      dat <- do.call(rbind, dat[includeFeatures])
      datCells <- do.call(rbind, datCells[includeFeatures])
    }
  }
  if (verbose) {
    message(sprintf("Loading extra analysis data where available"))
  }
  mDat <- NULL
  tgt <- file.path(
    dataDir, "analysis", "clustering", "graphclust",
    "clusters.csv"
  )
  if (file.exists(tgt)) {
    clusters <- read.csv(tgt)
    mDat <- data.frame(clusters = clusters$Cluster, row.names = clusters$Barcode)
  }
  tgt <- file.path(
    dataDir, "analysis", "clustering", "kmeans_10_clusters",
    "clusters.csv"
  )
  if (file.exists(tgt)) {
    clusters <- read.csv(tgt)
    mDat$clustersFine <- clusters$Cluster
  }
  tgt <- file.path(
    dataDir, "analysis", "tsne", "2_components",
    "projection.csv"
  )

  if (file.exists(tgt)) {
    tsne <- read.csv(tgt)
    if (is.null(mDat)) {
      mDat <- data.frame(
        tSNE1 = tsne$TSNE.1, tSNE2 = tsne$TSNE.2,
        row.names = tsne$Barcode
      )
    } else {
      mDat$tSNE1 <- tsne$TSNE.1[match(rownames(mDat), tsne$Barcode)]
      mDat$tSNE2 <- tsne$TSNE.2[match(rownames(mDat), tsne$Barcode)]
    }
    DR <- c("tSNE1", "tSNE2")
  } else {
    DR <- NULL
  }
  if (!is.null(mDat) && any(rownames(mDat) != colnames(datCells))) {
    rownames(mDat) <- gsub("-1$", "", rownames(mDat))
    if (any(rownames(mDat) != colnames(datCells))) {
      stop("Error matching meta-data to cell names.")
    }
  }
  if (is.null(channelName)) {
    channelName <- ifelse(is.null(names(dataDir)), dataDir,
      names(dataDir)
    )
  }

  "Maybe the one below should be within the above if statement?"
  channel <- SoupX::SoupChannel(
    tod = dat, toc = datCells, metaData = mDat,
    channelName = channelName, dataDir = dataDir, dataType = "10X",
    isV3 = isV3, DR = DR, ...
  )
  return(channel)
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




# _________________________________________________________________________________________________
# Deprecated ----
# _________________________________________________________________________________________________
save2umaps.A4 <- function(...) .Deprecated("save2plots.A4()")

save4umaps.A4 <- function(...) .Deprecated("save4plots.A4()")



# -----------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------
# Main script / functions


# will it be used?
cellID_to_cellType_v1 <- function(cellIDs, ident, obj = aaa) {
  celltypes <- as.named.vector.df(obj@meta.data[, ident], verbose = FALSE)
  celltypes[cellIDs]
}

cellID_to_cellType <- function(cellIDs, ident_w_names) {
  ident_w_names[cellIDs]
}


# _________________________________________________________________________________________________
#' @title Create.MiscSlot
#'
#' @description Create a new slot in the 'misc' slot of a Seurat object.
#' @param obj Seurat object
#' @param NewSlotName Name of the new element inside obj@misc.
#' @export

Create.MiscSlot <- function(obj, NewSlotName = "UVI.tables", SubSlotName = NULL) {
  .Deprecated("addToMiscOrToolsSlot")
  # if (is.null(obj@misc[[NewSlotName]])) obj@misc[[NewSlotName]] <- list() else iprint(NewSlotName, "already exists in @misc.")
  # if (is.null(obj@misc[[NewSlotName]][[SubSlotName]])) obj@misc[[NewSlotName]][[SubSlotName]] <- list() else iprint(SubSlotName, "subslot already exists in @misc$NewSlot.")
  return(obj)
}



# # _________________________________________________________________________________________________
# sparse.cor4 <- function(x){
#   n <- nrow(x)
#   cMeans <- colMeans(x)
#   covmat <- (as.matrix(crossprod(x)) - n*tcrossprod(cMeans))/(n-1)
#   sdvec <- sqrt(diag(covmat))
#   cormat <- covmat/tcrossprod(sdvec)
#   list(cov=covmat,cor=cormat)
# }


# # _________________________________________________________________________________________________
# #' @title Find Specific Files in Specified Subdirectories
# #'
# #' @description This function searches through specified subdirectories within a root directory
# #' to find files that match a specified pattern and returns a character vector with their full paths.
# #' The printed output excludes the root directory part from the paths.
# #'
# #' @param root_dir The root directory.
# #' @param subdir A character vector of subdirectory names within the root directory to be scanned.
# #' @param file_name_pattern The pattern of the file name to search for.
# #' @param recursive Boolean indicating whether to search recursively within subdirectories.
# #' @return A character vector containing the full paths to the located files.
# # #' @importFrom fs dir_ls
# #' @export

# findBamFilesInSubdirs <- function(root_dir, subdir, file_name_pattern = "possorted_genome_bam.bam", recursive = TRUE) {
#   stopifnot(is.character(root_dir), length(root_dir) == 1, dir.exists(root_dir),
#             is.character(subdir), all(dir.exists(file.path(root_dir, subdir))),
#             is.character(file_name_pattern), length(file_name_pattern) == 1,
#             is.logical(recursive))

#   pattern <- paste0("**/", file_name_pattern)
#   paths_to_search <- file.path(root_dir, subdir)
#   bams <- c()

#   for (path in paths_to_search) {
#     iprint("Searching in:", path)
#     found_files <- fs::dir_ls(path, recurse = recursive, glob = pattern, type = "file")
#     iprint(length(found_files), "files found.")
#     bams <- c(bams, found_files)
#   }

#   # Replace root_dir in the paths with an empty string for printing
#   bams_print <- gsub(paste0("^", root_dir, "/?"), "", bams)
#   iprint(length(bams), bams_print)

#   return(bams)
# }


# -----------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------
# VISUALIZATION

# panelCorPearson <- function(x, y, digits = 2, prefix = "", cex.cor = 2, method = "pearson") {
#   usr <- par("usr"); on.exit(par(usr))
#   par(usr = c(0, 1, 0, 1))
#   r <- abs(cor(x, y, method = method, use = "complete.obs"))
#   txt <- format(c(r, 0.123456789), digits = digits)[1]
#   txt <- paste(prefix, txt, sep = "")
#   if (missing(cex.cor)) cex <- 0.8/strwidth(txt)
#
#   test <- cor.test(x, y)
#   Signif <- symnum(test$p.value, corr = FALSE, na = FALSE,
#                    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
#                    symbols = c("***", "**", "*", ".", " "))
#
#   text(0.5, 0.5, txt, cex = cex * r)
#   text(.8, .8, Signif, cex = cex,  col = 2)
# }

# getDiscretePalette <- function(
    #     ident.used = GetClusteringRuns()[1],
#     obj = combined.obj,
#     palette.used = c("alphabet", "alphabet2", "glasbey", "polychrome", "stepped")[1],
#     show.colors = FALSE, seed = 1989) {
#
#   n.clusters <- nrow(unique(obj[[ident.used]]))
#
#   colorz <- Seurat::DiscretePalette(n = n.clusters, palette = palette.used)
#
#   if (anyNA(colorz)) {
#
#     colorsOK <- colorz[!is.na(colorz)] # Extract non-NA values
#     n.colz <- length(colorsOK)
#
#     msg <- paste("More categories then present in the palette", n.clusters, "vs."
#                  , n.colz, "in", palette.used, "-> recycling.")
#     warning(msg, immediate. = TRUE)
#
#     # Resample non-NA values and replace NA values
#     set.seed(seed)
#
#     if (n.clusters > 10 * n.colz) {
#       colorz <- sample(gplots::rich.colors(n.clusters))
#     } else {
#       colorz <- sample(x = colorsOK, size = n.clusters, replace = TRUE)
#     }
#
#     stopif(anyNA(colorz))
#
#   }
#   if (show.colors) MarkdownHelpers::color_check(colorz)
#   return(colorz)
# }


# -----------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------
# META

# transferMetadataV1 <- function(from, to, colname_from, colname_to = colname_from, verbose = TRUE, overwrite = FALSE) {
#
#   stopifnot(
#     is(from, "Seurat"), is(to, "Seurat"),
#     is.character(colname_from), is.character(colname_to),
#     "Column not found" = colname_from %in% colnames(from@meta.data),
#     "Column already exists" = !(colname_to %in% colnames(to@meta.data)) | overwrite
#   )
#
#   # Extract the metadata column to transfer
#   data.to.transfer <- data.frame(new.metadata = from[[colname_from]])
#
#   # Check cell overlaps
#   cells_in_both <- intersect(colnames(from), colnames(to))
#   cells_only_in_from <- setdiff(colnames(from), colnames(to))
#   cells_only_in_to <- setdiff(colnames(to), colnames(from))
#
#   if (verbose) {
#     cat("Number and % of cells matching between objects:", length(cells_in_both),
#         "(", sprintf("%.2f%%", length(cells_in_both) / length(colnames(from)) * 100), "of from and",
#         sprintf("%.2f%%", length(cells_in_both) / length(colnames(to)) * 100), "of to)\n")
#     cat("Number and % of cells only in obj1 (from):", length(cells_only_in_from),
#         "(", sprintf("%.2f%%", length(cells_only_in_from) / length(colnames(from)) * 100), ")\n")
#     cat("Number and % of cells only in obj2 (to):", length(cells_only_in_to),
#         "(", sprintf("%.2f%%", length(cells_only_in_to) / length(colnames(to)) * 100), ")\n")
#   }
#
#   # Add the metadata to the 2nd obj
#   to <- Seurat::AddMetaData(object = to, metadata = data.to.transfer, col.name = colname_to )
#
#   return(to)
# }
