# ____________________________________________________________________
# Seurat.utils.less.used.R ----
# ____________________________________________________________________
# file.edit("~/GitHub/Packages/Seurat.utils/R/Seurat.utils.less.used.R")
# setwd("~/GitHub/Packages/Seurat.utils")



# _________________________________________________________________________________________________
#' @title Convert10Xfolders - legacy version
#'
#' @description This function takes a parent directory with a number of subfolders, each
#' containing the standard output of 10X Cell Ranger. It (1) loads the filtered data matrices,
#' (2) converts them to Seurat objects, and (3) saves them as .RDS files.
#' @param InputDir A character string specifying the input directory.
#' @param regex A logical value. If TRUE, the folderPattern is treated as a regular expression. Default: `FALSE`.
#' @param folderPattern A character vector specifying the pattern of folder names to be searched. Default: 'filtered_feature'.
#' @param min.cells An integer value specifying the minimum number of cells. Default: 5.
#' @param min.features An integer value specifying the minimum number of features. Default: 200.
#' @param updateHGNC A logical value indicating whether to update the HGNC. Default: `TRUE`.
#' @param ShowStats A logical value indicating whether to show statistics. Default: `TRUE`.
#' @param writeCBCtable A logical value indicating whether to write out a list of cell barcodes (CBC) as a tsv file. Default: `TRUE`.
#' @param depth An integer value specifying the depth of scan (i.e., how many levels below the InputDir). Default: 2.
#' @param sample.barcoding A logical value indicating whether Cell Ranger was run with sample barcoding. Default: `FALSE`.
#' @param sort_alphanumeric Sort files alphanumerically? Default: `TRUE`.
#' @examples
#' \dontrun{
#' if (interactive()) Convert10Xfolders(InputDir)
#' }
#' @export
Convert10Xfolders_v1 <- function(
    InputDir,
    regex = FALSE,
    folderPattern = c("filtered_feature", "raw_feature", "SoupX_decont")[1],
    depth = 4,
    min.cells = 5, min.features = 200,
    updateHGNC = TRUE, ShowStats = TRUE,
    writeCBCtable = TRUE,
    sample.barcoding = FALSE,
    nthreads = .getNrCores(),
    preset = "high",
    ext = "qs",
    sort_alphanumeric = TRUE,
    ...) {
  warning("Since v2.5.0, the output is saved in the more efficient qs format! See qs package.", immediate. = TRUE)

  finOrig <- ReplaceRepeatedSlashes(list.dirs.depth.n(InputDir, depth = depth))
  fin <- CodeAndRoll2::grepv(x = finOrig, pattern = folderPattern, perl = regex)

  iprint(length(fin), "samples found.")

  samples <- basename(list.dirs(InputDir, recursive = FALSE))
  if (sort_alphanumeric) samples <- gtools::mixedsort(samples)
  iprint("Samples:", samples)

  if (!length(fin) > 0) {
    stop(paste("No subfolders found with pattern", folderPattern, "in dirs like: ", finOrig[1:3]))
  }

  for (i in 1:length(fin)) {
    print(i)
    pathIN <- Stringendo::FixPath(fin[i])
    print(pathIN)

    # sample.barcoding --- --- ---
    fnameIN <- if (sample.barcoding) {
      samples[i]
    } else {
      basename(dirname(dirname(pathIN)))
    }
    print("")
    print(fnameIN)

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

      LSBnameOUT <- ppp(paste0(InputDir, "/LSB.", fnameIN), "qs")
      qs::qsave(x = LSB, file = LSBnameOUT)
    } else {
      print("More than 2 elements in the list of matrices")
    }

    ncells <- ncol(seu)
    fname_X <- Stringendo::sppp(
      fnameIN, "min.cells", min.cells, "min.features", min.features,
      "cells", ncells
    )
    print(fname_X)

    f.path.out <- Stringendo::ParseFullFilePath(path = InputDir, file_name = fname_X, extension = ext)
    message(f.path.out)

    # update --- --- ---
    if (updateHGNC) seu <- UpdateGenesSeurat(seu, EnforceUnique = TRUE, ShowStats = TRUE)

    # write out --- --- ---
    qs::qsave(x = seu, file = f.path.out, nthreads = nthreads, preset = preset)

    # write cellIDs ---  --- ---
    if (writeCBCtable) {
      CBCs <- t(t(colnames(seu)))
      colnames(CBCs) <- "CBC"
      ReadWriter::write.simple.tsv(input_df = CBCs, manual_file_name = sppp(fnameIN, "CBC"), manual_directory = InputDir)
    }
  } # for
}




# ____________________________________________________________________________________
#' @title plot.UMAP.tSNE.sidebyside
#'
#' @description Plot a UMAP and tSNE side by side.
#' @param obj Seurat object. Default: combined.obj
#' @param grouping Variable to group cells by. Default: 'res.0.6'
#' @param no_legend Logical, whether to display legend. Default: `FALSE`.
#' @param do_return Logical, whether to return plot object. Default: `TRUE`.
#' @param do_label Logical, whether to display labels. Default: `TRUE`.
#' @param label_size Size of labels. Default: 10
#' @param vector_friendly Logical, whether to optimize for vector outputs. Default: `TRUE`.
#' @param cells_use A vector of cell names to use for the plot. Default: NULL
#' @param no_axes Logical, whether to hide axes. Default: `TRUE`.
#' @param pt_size Size of points. Default: 0.5
#' @param name.suffix Suffix to append to the plot's name. Default: NULL
#' @param width Width of the plot. Default: hA4
#' @param height Height of the plot. Default: 1.75 * wA4
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
                                      width = hA4, height = 1.75 * wA4, filetype = "pdf", ...) {
  p1 <- Seurat::DimPlot(
    object = obj, reduction = "tsne", cells = cells_use,
    group.by = grouping, label = do_label, label.size = label_size,
    pt.size = pt_size, raster = !vector_friendly, ...
  ) +
    ggtitle("tSNE") + theme(plot.title = element_text(hjust = 0.5))

  if (no_axes) p1 <- p1 + Seurat::NoAxes()
  if (no_legend) p1 <- p1 + Seurat::NoLegend()

  p2 <- Seurat::DimPlot(
    object = obj, reduction = "umap", cells = cells_use,
    group.by = grouping, label = do_label, label.size = label_size,
    pt.size = pt_size, raster = !vector_friendly, ...
  ) +
    ggtitle("UMAP") + theme(plot.title = element_text(hjust = 0.5))

  if (no_axes) p2 <- p2 + Seurat::NoAxes()
  p2 <- p2 + Seurat::NoLegend()

  plots <- cowplot::plot_grid(p1, p2, labels = c("A", "B"), ncol = 2)
  plotname <- kpp("UMAP.tSNE", grouping, name.suffix, filetype)

  cowplot::save_plot(
    filename = plotname, plot = plots,
    ncol = 2 # we're saving a grid plot of 2 columns
    , nrow = 1 # and 2 rows
    , base_width = width,
    base_height = height
    # each individual subplot should have an aspect ratio of 1.3
    # , base_aspect_ratio = 1.5
  )

  if (isTRUE(do_return)) return(plots)

  invisible(NULL)
}


# _________________________________________________________________________________________________
#' @title Plot multiple categorical variables in combined UMAPs
#'
#' @description Generates and saves multiple UMAP plots for clustering results, adjusting the
#' layout and plot dimensions. Supports the generation of plots in different
#' formats and customization of the visual appearance.
#'
#' @param idents A vector of cluster identities to plot. Default: `GetClusteringRuns()[1:4]`.
#' @param obj The Seurat object containing clustering information. Default: `combined.obj`.
#' @param foldername The name of the folder to save plots. Default: `substitute(ident)`.
#' @param plot.reduction The dimensionality reduction technique to use for plotting. Default: "umap".
#' @param intersectionAssay The assay to use for intersection. Default: "RNA".
#' @param layout The layout orientation, either "tall", "wide", or `FALSE` to disable. Default: "wide".
#' @param nr.Col Number of columns in the plot grid. Default: 2.
#' @param nr.Row Number of rows in the plot grid. Default: 4.
#' @param cex The character expansion size for plot text, automatically adjusted. Default: `round(0.1 / (nr.Col * nr.Row), digits = 2)`.
#' @param label Logical indicating if labels should be displayed on the plots. Default: `FALSE`.
#' @param legend Logical indicating if a legend should be included in the plots. Default: `!label`.
#' @param subdir Logical indicating if a subdirectory should be created for saving plots. Default: `TRUE`.
#' @param prefix Optional prefix for plot filenames. Default: `NULL`.
#' @param suffix Optional suffix for plot filenames. Default: `NULL`.
#' @param background_col The background color of the plot. Default: "white".
#' @param aspect.ratio The aspect ratio of the plot, `FALSE` to disable fixed ratio. Default: 0.6.
#' @param saveGeneList Logical indicating if a list of genes should be saved. Default: `FALSE`.
#' @param w The width of the plot in inches. Default: `8.27`.
#' @param h The height of the plot in inches. Default: `11.69`.
#' @param scaling The scaling factor to apply to plot dimensions. Default: 1.
#' @param format The file format for saving plots. Default: "jpg".
#' @param ... Additional arguments passed to plotting functions.
#'
#' @return Invisible `NULL`. Plots are saved to files.
#' @examples
#' \dontrun{
#' multi_clUMAP.A4(idents = c("S1", "S2"), obj = YourSeuratObject)
#' }
#' @export

multi_clUMAP.A4 <- function(
    obj = combined.obj,
    idents = GetClusteringRuns(obj)[1:4],
    foldername = "clUMAPs_multi",
    plot.reduction = "umap",
    intersectionAssay = c("RNA", "integrated")[1],
    layout = c("tall", "wide", FALSE)[2],
    # colors = c("grey", "red"),
    nr.Col = 2, nr.Row = 4,
    cex = round(0.1 / (nr.Col * nr.Row), digits = 2),
    label = FALSE, # can be a vector of length idents
    legend = !label,
    subdir = TRUE,
    prefix = NULL, suffix = NULL,
    background_col = "white",
    aspect.ratio = c(FALSE, 0.6)[2],
    saveGeneList = FALSE,
    w = 8.27, h = 11.69, scaling = 1,
    format = c("jpg", "pdf", "png")[1],
    ...) {
  .Deprecated("qClusteringUMAPS")
  message("multi_clUMAP.A4() is kept because it can plot more than 4 resolutions, inti a subfolder.")

  tictoc::tic()
  ParentDir <- OutDir
  if (is.null(foldername)) foldername <- "clusters"
  if (subdir) create_set_SubDir(paste0(foldername, "-", plot.reduction), "/")

  DefaultAssay(obj) <- intersectionAssay

  # Adjust plot dimensions and grid layout based on specified layout
  .adjustLayout(layout, scaling, wA4 = 8.27, hA4 = 11.69, environment())


  # Split clusters into lists for plotting
  ls.idents <- CodeAndRoll2::split_vec_to_list_by_N(1:length(idents), by = nr.Row * nr.Col)
  for (i in 1:length(ls.idents)) {
    idents_on_this_page <- idents[ls.idents[[i]]]
    iprint("page:", i, "| idents", kppc(idents_on_this_page))
    (plotname <- kpp(c(prefix, plot.reduction, i, "idents", ls.idents[[i]], suffix, format)))

    plot.list <- list()
    for (i in seq(idents_on_this_page)) {
      # browser()
      if (length(label) == 1) {
        label_X <- label
        legend_X <- legend
      } else {
        label_X <- label[i]
        legend_X <- legend[i]
      }

      ident_X <- idents_on_this_page[i]
      imessage("plotting:", ident_X)
      plot.list[[i]] <- clUMAP(
        ident = ident_X, obj = obj, plotname = label_X,
        label = label_X, legend = legend_X, save.plot = FALSE, h = h, w = w, ...
      )
    }

    # Customize plot appearance
    for (i in 1:length(plot.list)) {
      plot.list[[i]] <- plot.list[[i]] + NoAxes()
      if (aspect.ratio) plot.list[[i]] <- plot.list[[i]]
      ggplot2::coord_fixed(ratio = aspect.ratio)
    }

    # Save plots
    pltGrid <- cowplot::plot_grid(plotlist = plot.list, ncol = nr.Col, nrow = nr.Row)
    cowplot::ggsave2(filename = plotname, width = w, height = h, bg = background_col, plot = pltGrid)
  } # for ls.idents

  if (subdir) MarkdownReports::create_set_OutDir(ParentDir)
  tictoc::toc()
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
#' @description Relabel cluster numbers based on their position along the principal curve fitted to the specified dimensionality reduction (UMAP, tSNE, or PCA).
#' @param obj Seurat object, Default: combined.obj
#' @param dim Dimensions to use, Default: 1:2
#' @param plotit Plot results (& show it), Default: `TRUE`.
#' @param swap Swap Lambda parameter (multiplied with this) , Default: -1
#' @param reduction UMAP, tSNE, or PCA (Dim. reduction to use), Default: 'umap'
#' @param res Clustering resolution to use, Default: 'integrated_snn_res.0.5'
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
  NewMeta <- translate(vec = obj[[res]], old = OldLabel, new = NewLabel)
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
#' @param cellIDs An optional vector of cell IDs to include in the loaded data. Default: `NULL`,
#' indicating that all available cells will be included. This is useful for subsetting the data based
#' on specific cell IDs.
#' @param channelName An optional string specifying the channel name for the data being loaded.
#' This can be used to label the data according to the experimental condition or sample name. Default: `NULL`.
#' @param readArgs A list of additional arguments to pass to the internal `Read10X` function used for
#' loading the data. Default: an empty list.
#' @param includeFeatures A character vector specifying which features to include in the loaded data.
#' Common values include "Gene Expression", "Antibody Capture", and "CRISPR Guide Capture".
#' Default: `c("Gene Expression")`.
#' @param verbose A logical flag indicating whether to print progress messages and status updates as the
#' data is loaded. Default: `TRUE`.
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

  stopifnot("Package 'SoupX' must be installed to use this function." = require("SoupX"))
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
#' @param folderPattern A character vector specifying the pattern of folder names to be searched. Default: 'filtered'.
#' @param min.cells An integer value specifying the minimum number of cells. Default: 10.
#' @param min.features An integer value specifying the minimum number of features. Default: 200.
#' @param updateHGNC A logical value indicating whether to update the HGNC. Default: `TRUE`.
#' @param ShowStats A logical value indicating whether to show statistics. Default: `TRUE`.
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
# Layer Removal ----
# _________________________________________________________________________________________________

# _________________________________________________________________________________________________
#' @title Remove Scale Data from Seurat Objects
#'
#' @param ls.obj A list of Seurat objects.
#' @return A list of Seurat objects with `scale.data` slot removed from RNA assays.
#' @examples
#' # Assuming `seuratList` is a list of Seurat objects
#' seuratList <- removeScaleData(seuratList)
#' @export
removeScaleData <- function(ls.obj) {
  lapply(ls.obj, function(x) {
    x@assays$RNA@layers$scale.data <- NULL
    x
  })
}


# _________________________________________________________________________________________________
#' @title Remove Layers from Seurat Object by Pattern
#'
#' @description This function removes layers from a Seurat object's RNA assay based on a specified regular expression pattern.
#' It first backs up the object before removing layers that match the pattern.
#'
#' @param obj A Seurat object.
#' @param pattern A regular expression pattern to match layer names.
#' @param perl A logical value indicating whether to use Perl-compatible regular expressions.
#' Default: `TRUE`.
#'
#' @importFrom CodeAndRoll2 grepv
#' @return A Seurat object with specified layers removed.
#' @export
removeLayersByPattern <- function(obj, pattern = "sc[0-9][0-9]_", perl = TRUE) {
  message(paste("pattern: ", pattern))
  stopifnot("obj must be a Seurat object" = inherits(obj, "Seurat"))

  layerNames <- Layers(obj)
  layersToRemove <- CodeAndRoll2::grepv(pattern, x = layerNames, perl = perl)
  message(paste(length(layersToRemove), "form", length(layerNames), "layers are removed."))
  obj@assays$RNA@layers[layersToRemove] <- NULL
  return(obj)
}

# _________________________________________________________________________________________________
# Deprecated ----
# _________________________________________________________________________________________________
set.all.genes <- function() .Deprecated("calc.q99.Expression.and.set.all.genes()")
save2umaps.A4 <- function(...) .Deprecated("save2plots.A4()")
save4umaps.A4 <- function(...) .Deprecated("save4plots.A4()")
plotGeneExpHist <- function(...) .Deprecated("plotGeneExprHistAcrossCells()")
geneExpressionLevelPlots <- function(...) .Deprecated("plotGeneExpressionInBackgroundHist()")
get.clustercomposition <- function(...) .Deprecated("No longer provided.")
multiFeatureHeatmap.A4 <- function(...) .Deprecated("No longer provided.")
Annotate4Plotly3D <- function(...) .Deprecated(".Annotate4Plotly3D() - with dot/invisible.")
Percent.in.Trome <- function(...) .Deprecated("PercentInTranscriptome()")
.parseRegressionVariablesForScaleData <- function(...) .Deprecated(".getRegressionVariablesForScaleData()")
seu.add.meta.from.vector <- function(...) .Deprecated("addMetaDataSafe()")

# _________________________________________________________________________________________________
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




# _________________________________________________________________________________________________
# Archived ----
# _________________________________________________________________________________________________

".Deprecated"

#' @title Regress Out and Recalculate Seurat
#'
#' @description The function performs a series of calculations and manipulations on a Seurat object,
#' including identifying variable features, scaling data, running PCA, setting up reductions, finding neighbors,
#' and finding clusters. It optionally performs t-SNE and saves the object.
#'
#' @param obj The Seurat object.
#' @param n.var.features The number of variable features to use. Default: the 'n.var.genes' element from a list 'p'.
#' @param features.scale A logical value indicating whether to scale the data. Default: `TRUE`.
#' @param vars.to.regress A vector of variable names to be regressed out.
#' @param suffix A character string to be used as a suffix when saving the object.
#' @param nPCs The number of principal components to use. Default: the 'n.PC' element from a list 'p'.
#' @param clust_resolutions The resolution for clustering. Default: the 'snn_res' element from a list 'p'.
#' @param calc_tSNE Logical, if TRUE, t-SNE will be performed. Default: `FALSE`.
#' @param plot_umaps Logical, if TRUE, UMAP plots will be generated. Default: `TRUE`.
#' @param save_obj Logical, if TRUE, the object will be saved. Default: `TRUE`.
#' @param assayX The assay to be used in scaling data. Default: 'RNA'.
#' @return Seurat object after calculations and manipulations.
#' @importFrom Seurat FindVariableFeatures ScaleData RunPCA FindNeighbors FindClusters RunTSNE
#' @importFrom MarkdownReports create_set_OutDir
#' @examples
#' \dontrun{
#' # Assuming 'seurat_obj' is a valid Seurat object and 'vars' is a vector of variable names to be regressed out.
#' result <- regress_out_and_recalculate_seurat(seurat_obj, vars, suffix = "_regressed")
#' }
#' @importFrom tictoc tic toc
#'
#' @export
regress_out_and_recalculate_seurat <- function(
    obj,
    n.var.features = p$"n.var.genes", # p is a list of parameters, 2000
    features.scale = n.var.features,
    vars.to.regress,
    suffix,
    nPCs = p$"n.PC",
    clust_resolutions = p$"snn_res",
    calc_tSNE = FALSE,
    plot_umaps = TRUE,
    save_obj = TRUE,
    assayX = "RNA") {
  .Deprecated("processSeuratObject")

  tictoc::tic()
  print("FindVariableFeatures")
  obj <- FindVariableFeatures(obj, mean.function = "FastExpMean", dispersion.function = "FastLogVMR", nfeatures = n.var.features)
  tictoc::toc()

  tictoc::tic()
  print("calc.q99.Expression.and.set.all.genes")
  obj <- calc.q99.Expression.and.set.all.genes(obj = obj, quantileX = .99)
  tictoc::toc()

  tictoc::tic()
  print("ScaleData")
  obj <- ScaleData(obj, assay = assayX, verbose = TRUE, vars.to.regress = vars.to.regress, features = features.scale)
  tictoc::toc()

  tictoc::tic()
  print("RunPCA")
  obj <- RunPCA(obj, npcs = nPCs, verbose = TRUE)
  tictoc::toc()

  tictoc::tic()
  print("SetupReductionsNtoKdimensions")
  obj <- SetupReductionsNtoKdimensions(obj = obj, nPCs = nPCs, dimensions = 3:2, reduction = "umap")
  tictoc::toc()

  tictoc::tic()
  print("FindNeighbors")
  obj <- FindNeighbors(obj, reduction = "pca", dims = 1:nPCs)
  tictoc::toc()

  tictoc::tic()
  print("FindClusters")
  obj <- FindClusters(obj, resolution = clust_resolutions)
  tictoc::toc()

  if (calc_tSNE) {
    tictoc::tic()
    print("RunTSNE")
    obj <- RunTSNE(obj, reduction = "pca", dims = 1:nPCs)
    tictoc::toc()
  }

  # orig.dir <- getwd()
  # new_path <- FixPath(orig.dir, suffix)
  # MarkdownReports::create_set_OutDir(new_path)

  clz <- GetClusteringRuns(obj, pat = "*snn_res.*[0-9]$")

  if (plot_umaps) {
    print("Plotting umaps")
    for (v in clz) clUMAP(ident = v, obj = obj, sub = suffix)

    # MarkdownReports::create_set_OutDir(new_path, 'UMAP_stats')
    for (v in vars.to.regress) qUMAP(feature = v, obj = obj, sub = suffix)
    # MarkdownReports::create_set_OutDir(new_path)
  }


  if (save_obj) {
    print("Save RDS")
    isave.RDS(obj, suffix = suffix, inOutDir = TRUE)
  }

  return(obj)
}


# _________________________________________________________________________________________________
#' @title Proportion of Cells Expressing Given Genes
#'
#' @description Calculates the proportion of cells expressing one or more specified genes.
#'
#' @param genes Character vector of gene names of interest.
#' @param group.by Optional grouping variable for analysis (e.g., cell type). Default: 'all'.
#' @param obj Seurat object to analyze. Default: `combined.obj`.
#' @param ... Additional arguments.
#'
#' @return Data frame with genes and their cell expression proportion, optionally grouped.
#'
#' @examples
#' \dontrun{
#' PrctCellExpringGene(genes = c("LTB", "GNLY"), obj = combined.obj)
#' }
#'
#' @source Adapted from Ryan-Zhu on GitHub.
#'
#' @export
PrctCellExpringGene <- function(genes, group.by = "all", obj = combined.obj,
                                ...) {
  .Deprecated("PctCellsExpressingGenes")
  #
  nf <- setdiff(genes, c(Features(obj, assay = "RNA"), colnames(obj@m@data)))

  if (length(nf) > 0) message("Some genes/ features not found: ", nf)

  stopifnot("Some genes not foun!." = all(genes %in% Features(obj)))

  if (group.by == "all") {
    prct <- 1:length(genes)
    for (i in seq(prct)) prct[i] <- ww.calc_helper(genes = genes[1], obj = obj)
    result <- data.frame("Markers" = genes, "Cell_proportion" = prct)
    return(result)
  } else {
    ls.Seurat <- Seurat::SplitObject(object = obj, split.by = group.by)
    factors <- names(ls.Seurat)

    # This is a self referencing function, how does this supposed to even work??
    results <- lapply(ls.Seurat, PrctCellExpringGene, genes = genes)
    for (i in 1:length(factors)) {
      results[[i]]$Feature <- factors[i]
    }
    combined <- do.call("rbind", results)
    return(combined)
  }
}


# _________________________________________________________________________________________________
#' @title Helper to calculate Cell Expression Proportion for Gene
#'
#' @description Computes the proportion of cells expressing a specific gene within a Seurat object.
#'
#' @param obj Seurat object with cell data.
#' @param genes Single gene name as a character string.
#' @param slot Slot to use for the analysis. Default: 'RNA'.
#'
#' @return Proportion of cells expressing the gene. Returns `NA` if the gene is not found.
#'
#' @examples
#' \dontrun{
#' ww.calc_helper(obj = seurat_object, genes = "Gene1")
#' }
#'
#' @source Adapted from Ryan-Zhu on GitHub.
#'
#' @export
ww.calc_helper <- function(obj, genes, slot = "RNA") {
  .Deprecated("Unused function.")
  # stopifnot("Some genes not found!." = all(genes %in% row.names(obj)))
  counts <- obj[[slot]]@counts
  ncells <- ncol(counts)
  if (genes %in% row.names(counts)) {
    sum(counts[genes, ] > 0) / ncells
  } else {
    return(NA)
  }
}



#' @title Cluster Size Distribution Plot (Barplot or Histogram)
#'
#' @description Generates a bar plot or histogram to visualize the size distribution of clusters
#' within a Seurat object, based on the specified clustering identity.
#'
#' @param obj Seurat object for analysis. Default: `combined.obj`.
#' @param ident Clustering identity to base the plot on.
#' Default: The second entry from `GetClusteringRuns()`.
#' @param plot Whether to display the plot (TRUE) or return cluster sizes (FALSE). Default: `TRUE`.
#' @param thr.hist Threshold for switching from a bar plot to a histogram based on the number of
#' clusters. Default: 30.
#' @param ... Extra parameters for the plot.
#'
#' @examples
#' \dontrun{
#' plotClustSizeDistr()
#' }
#'
#' @importFrom ggExpress qbarplot qhistogram
#'
# #' @export
plotClustSizeDistr <- function(
    obj = combined.obj, ident,
    plot = TRUE, thr.hist = 30,
    ...) {
  .Deprecated("plotClusterSizeDistribution()")

  stopifnot(ident %in% colnames(obj@meta.data))

  clust.size.distr <- table(obj@meta.data[, ident])
  print(clust.size.distr)
  resX <- gsub(pattern = ".*res\\.", replacement = "", x = ident)
  ptitle <- paste("Cluster sizes at ", ident)
  psubtitle <- paste(
    "Nr.clusters:", length(clust.size.distr),
    "| median size:", median(clust.size.distr),
    "| CV:", percentage_formatter(cv(clust.size.distr))
  )
  xlb <- "Cluster size (cells)"
  ylb <- "Nr of Clusters"
  xlim <- c(0, max(clust.size.distr))

  if (plot) {
    if (length(clust.size.distr) < thr.hist) {
      ggExpress::qbarplot(clust.size.distr,
                          plotname = ptitle, subtitle = psubtitle,
                          label = clust.size.distr, xlab = "Clusters", ylab = xlb, ...
      )
    } else {
      ggExpress::qhistogram(
        vec = clust.size.distr, plotname = ptitle, subtitle = psubtitle,
        xlab = xlb, ylab = ylb, xlim = xlim, ...
      )
    }
  } else {
    "return vector"
    clust.size.distr
  }
}




#' # _________________________________________________________________________________________________
#' #' @title seu.add.meta.from.vector
#' #'
#' #' @description Adds a new metadata column to a Seurat object.
#' #' @param obj A Seurat object to which the new metadata column will be added. Default: combined.obj.
#' #' @param metaD.colname A string specifying the name of the new metadata column. Default: metaD.colname.labeled.
#' #' @param Label.per.cell A vector of labels for each cell, to be added as new metadata. Default: Cl.Label.per.cell.
#' #' @return A Seurat object with the new metadata column added.
#' #' @examples
#' #' \dontrun{
#' #' if (interactive()) {
#' #'   # Example usage:
#' #'   combined.obj <- seu.add.meta.from.vector(
#' #'     obj = combined.obj,
#' #'     metaD.colname = metaD.colname.labeled,
#' #'     Label.per.cell = Cl.Label.per.cell
#' #'   )
#' #' }
#' #' }
#' #' @export
#' seu.add.meta.from.vector <- function(obj = combined.obj, metaD.colname, Label.per.cell = Cl.Label.per.cell) {
#'   .Deprecated("addMetaDataSafe")
#'   obj@meta.data[, metaD.colname] <- Label.per.cell
#'   iprint(metaD.colname, "contains the named identities. Use Idents(combined.obj) = '...'. The names are:", unique(Label.per.cell))
#'   return(obj)
#' }
#'
#'

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


# _________________________________________________________________________________________________
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


# _________________________________________________________________________________________________
# _________________________________________________________________________________________________
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

# _________________________________________________________________________________________________

# #' @title Cluster Composition Analysis
# #'
# #' @description Analyzes and visualizes the composition of clusters in a Seurat object, indicating
# #' the contribution of different datasets to each cluster.
# #'
# #' @param obj Seurat object to analyze. Default: `combined.obj`.
# #' @param ident Cluster identity resolution to use. Default: 'integrated_snn_res.0.3'.
# #' @param splitby Variable to split the data by, typically a project or dataset identifier.
# #' Default: 'ShortNames'.
# #' @param color Bar color. Default: as defined by `splitby`.
# #' @param plot Whether to display the plot. Default: `TRUE`.
# #' @param ScaleTo100pc Whether to scale Y axis to 100%. Default: `TRUE`.
# #' @param ... Additional parameters for plotting functions.
# #'
# #' @return If `plot` is TRUE, displays a bar plot showing the composition of each cluster. Otherwise,
# #' performs the analysis without plotting.
# #'
# #' @examples
# #' get.clustercomposition()
# #'
# #' @export
# #' @importFrom dplyr group_by_ summarise
# #' @importFrom scales percent_format
# get.clustercomposition <- function(
#     obj = combined.obj,
#     ident = GetClusteringRuns()[1],
#     splitby = "orig.ident",
#     color = splitby,
#     plot = TRUE, ScaleTo100pc = TRUE,
#     ...) {
#
#   stopifnot(ident %in% colnames(obj@meta.data),
#             splitby %in% colnames(obj@meta.data)
#   )
#   (df.meta <- obj@meta.data[, c(ident, splitby)])
#
#   try(setwd(OutDir), silent = TRUE)
#
#   df.meta |>
#     dplyr::group_by_(splitby) |>
#     summarise()
#
#   categ.per.cluster <- ggbarplot(obj@meta.data,
#                                  x = ident,
#                                  y = splitby,
#                                  color = splitby,
#                                  ...
#   )
#   if (ScaleTo100pc) categ.per.cluster <- categ.per.cluster + scale_y_discrete(labels = scales::percent_format())
#   if (plot) categ.per.cluster
#
#   # ggExpress::qqSave(categ.per.cluster, ...)
# }



# _________________________________________________________________________________________________
# #' @title Cluster Composition Analysis
# #'
# #' @description Analyzes and visualizes the composition of clusters in a Seurat object, indicating
# #' the contribution of different datasets to each cluster.
# #'
# #' @param obj Seurat object to analyze. Default: `combined.obj`.
# #' @param ident Cluster identity resolution to use. Default: 'integrated_snn_res.0.3'.
# #' @param splitby Variable to split the data by, typically a project or dataset identifier.
# #' Default: 'ShortNames'.
# #' @param color Bar color. Default: as defined by `splitby`.
# #' @param plot Whether to display the plot. Default: `TRUE`.
# #' @param ScaleTo100pc Whether to scale Y axis to 100%. Default: `TRUE`.
# #' @param ... Additional parameters for plotting functions.
# #'
# #' @return If `plot` is TRUE, displays a bar plot showing the composition of each cluster. Otherwise,
# #' performs the analysis without plotting.
# #'
# #' @examples
# #' get.clustercomposition()
# #'
# #' @export
# #' @importFrom dplyr group_by_ summarise
# #' @importFrom scales percent_format
# get.clustercomposition <- function(
#     obj = combined.obj,
#     ident = GetClusteringRuns()[1],
#     splitby = "orig.ident",
#     color = splitby,
#     plot = TRUE, ScaleTo100pc = TRUE,
#     ...) {
#
#   stopifnot(ident %in% colnames(obj@meta.data),
#             splitby %in% colnames(obj@meta.data)
#             )
#   (df.meta <- obj@meta.data[, c(ident, splitby)])
#
#   try(setwd(OutDir), silent = TRUE)
#
#   df.meta |>
#     dplyr::group_by_(splitby) |>
#     summarise()
#
#   categ.per.cluster <- ggbarplot(obj@meta.data,
#                                  x = ident,
#                                  y = splitby,
#                                  color = splitby,
#                                  ...
#   )
#   if (ScaleTo100pc) categ.per.cluster <- categ.per.cluster + scale_y_discrete(labels = scales::percent_format())
#   if (plot) categ.per.cluster
#
#   # ggExpress::qqSave(categ.per.cluster, ...)
# }


# # _________________________________________________________________________________________________
# # Save multiple FeatureHeatmaps from a list of genes on A4 jpeg
# # code for quantile: https://github.com/satijalab/seurat/blob/master/R/plotting_internal.R
#
# #' @title multiFeatureHeatmap.A4
# #'
# #' @description Save multiple FeatureHeatmaps from a list of genes on A4 jpeg.
# #' @param obj Seurat object, Default: combined.obj
# #' @param list.of.genes A list of genes to plot. No default.
# #' @param gene.per.page Number of genes to plot per page. Default: 5
# #' @param group.cells.by Cell grouping variable for the heatmap. Default: 'batch'
# #' @param plot.reduction Dimension reduction technique to use for plots. Default: 'umap'
# #' @param cex Point size in the plot. Default: iround(3/gene.per.page)
# #' @param sep_scale Logical, whether to scale the features separately. Default: `FALSE`.
# #' @param gene.min.exp Minimum gene expression level for plotting. Default: 'q5'
# #' @param gene.max.exp Maximum gene expression level for plotting. Default: 'q95'
# #' @param jpeg.res Resolution of the jpeg output. Default: 225
# #' @param jpeg.q Quality of the jpeg output. Default: 90
# #' @param ... Pass any other parameter to the internally called functions (most of them should work).
# #' @seealso
# #'  \code{\link[tictoc]{tic}}
# #' @importFrom tictoc tic toc
# #'
# #' @export
# multiFeatureHeatmap.A4 <- function(
#     obj = combined.obj,
#     list.of.genes, gene.per.page = 5,
#     group.cells.by = "batch", plot.reduction = "umap",
#     cex = iround(3 / gene.per.page), sep_scale = FALSE,
#     gene.min.exp = "q5", gene.max.exp = "q95",
#     jpeg.res = 225, jpeg.q = 90,
#     ...) {
#   tictoc::tic()
#   list.of.genes <- check.genes(list.of.genes, obj = obj)
#
#   lsG <- CodeAndRoll2::split_vec_to_list_by_N(1:length(list.of.genes), by = gene.per.page)
#   for (i in 1:length(lsG)) {
#     print(i)
#     genes <- list.of.genes[lsG[[i]]]
#     plotname <- kpp(c("FeatureHeatmap", plot.reduction, i, genes, "jpg"))
#     print(plotname)
#     jjpegA4(plotname, r = jpeg.res, q = jpeg.q)
#     try(
#       FeatureHeatmap(obj,
#                      features.plot = genes, group.by = group.cells.by,
#                      reduction.use = plot.reduction, do.return = FALSE,
#                      sep.scale = sep_scale, min.exp = gene.min.exp, max.exp = gene.max.exp,
#                      pt.size = cex, key.position = "top", ...
#       ),
#       silent = FALSE
#     )
#     try.dev.off()
#   }
#   tictoc::toc()
# }
#


# # _________________________________________________________________________________________________
# #' @title ww.check.if.3D.reduction.exist
# #'
# #' @description ww.check.if.3D.reduction.exist in backup slot #
# #' @param obj Seurat object, Default: obj
# #' @export
# ww.check.if.3D.reduction.exist <- function(obj = obj) {
#   if (!("UMAP_3" %in% colnames(obj@reductions$"umap"))) {
#     stopif(
#       is.null(obj@misc$reductions.backup$"umap3d"),
#       "No 3D umap found in backup slot, @misc$reductions.backup. Run SetupReductionsNtoKdimensions() first."
#     )
#     RecallReduction(obj = obj, dim = 3, reduction = "umap")
#   } else { # Reduction found in normal UMAP slot
#     obj
#   }
# }
