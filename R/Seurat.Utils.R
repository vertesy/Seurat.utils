# ____________________________________________________________________
# Seurat.utils ----
# ____________________________________________________________________
# source("~/GitHub/Packages/Seurat.utils/R/Seurat.Utils.R")
# source("~/GitHub/Packages/Seurat.utils/R/Seurat.Utils.Visualization.R")
# source("~/GitHub/Packages/Seurat.utils/R/Seurat.Utils.Metadata.R")

# devtools::load_all("~/GitHub/Packages/Seurat.utils")
# devtools::document("~/GitHub/Packages/Seurat.utils"); devtools::load_all("~/GitHub/Packages/Seurat.utils")

# _________________________________________________________________________________________________
# Cluster.Auto-naming.DE.R
# _________________________________________________________________________________________________
# source('~/GitHub/Packages/Seurat.utils/Functions/Cluster.Auto-naming.DE.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Cluster.Auto-naming.DE.R"))

# _________________________________________________________________________________________________
# require(princurve) # only for AutoNumber.by.PrinCurve


# _________________________________________________________________________________________________
# General ______________________________ ----
# _________________________________________________________________________________________________


#' @title parallel.computing.by.future
#'
#' @description Run gc(), load multi-session computing and extend memory limits.
#' @param cores Number of cores
#' @param maxMemSize memory limit
#'
#' @export
parallel.computing.by.future <- function(cores = 4, maxMemSize = 4000 * 1024^2) { # Run gc(), load multi-session computing and extend memory limits.
  # https://satijalab.org/seurat/v3.0/future_vignette.html
  cat(
    "1. If you load futures before you finished using foreach loops,
    NormalizeData inside a foreach loop fails (Error writing to connection)
    -> I assume 'future' and 'doMC' are not compatible

    2. If you setup computing on e.g. six cores, it runs 6 instances of R with the entire memory space copied.
    If you run out of memory, the system starts using the SSD as memory, and it slows you down extremely extremely extremely.
    -> Therefore it is important to clean up the memory space before setting up multicore computation.

    Loaded: library(future), workers set to 6 (def),set Max mem size to 2GB (def)."
  )

  gc(full = TRUE)
  try(memory.biggest.objects(), silent = TRUE)
  user_input <- readline(prompt = "Are you sure that memory should not be cleaned before paralellizng? (y/n)")

  if (user_input == "y") {
    iprint("N. cores", cores)
    library(future)
    # plan("multiprocess", workers = cores)
    plan("multisession", workers = cores)
    # So to set Max mem size to 2GB, you would run :
    options(future.globals.maxSize = maxMemSize)
  } else {
    print("No parallelization")
  }
}


# _________________________________________________________________________________________________
#' @title Intersect Genes with Seurat Object
#'
#' @description Intersects a set of gene names with those found in a Seurat object.
#' @param genes A vector of gene names to be intersected with the Seurat object.
#' @param obj A Seurat object containing gene expression data.
#' @param n_genes_shown Number of missing genes to be printed. Default: 10.
#' @param strict All genes to be present in the Seurat object?  Default: TRUE.
#' @return A vector of gene names that are found both in the input 'genes' vector and the
#'         Seurat object.
#'
#' @export
IntersectGeneLsWithObject <- function(genes, obj = combined.obj, n_genes_shown = 10, strict = TRUE) {
  message(">>> Running IntersectGeneLsWithObject(), formerly IntersectWithExpressed()")

  stopifnot(
    is.character(genes),
    is(obj, "Seurat"),
    is.numeric(n_genes_shown) && n_genes_shown > 0,
    is.logical(strict)
  )
  stopifnot(length(genes) > 0, length(rownames(obj)) > 0)

  # Strict mode: Ensure all genes are present in the Seurat object
  if (strict) stopifnot(all(genes %in% rownames(obj)))

  # Finding genes that are missing in the Seurat object
  missing_in_obj <- setdiff(genes, rownames(obj))
  Stringendo::iprint(length(missing_in_obj), " (of ", length(genes),
                     ") genes are MISSING from the Seurat object with (", length(rownames(obj)),
                     ") genes. E.g.:", head(missing_in_obj, n_genes_shown))

  # Finding genes that are found in both the input list and the Seurat object
  g_found <- intersect(rownames(obj), genes)

  # Output argument assertion
  stopifnot(length(g_found) > 0)

  return(g_found)
}


# _________________________________________________________________________________________________
#' @title SmallestNonAboveX
#'
#' @description replace small values with the next smallest value found, which is >X.
#' @param vec Numeric input vector
#' @param X Threshold, Default: 0
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   SmallestNonZero(vec = df.markers$"p_val")
#' }
#' }
#' @export
SmallestNonAboveX <- function(vec, X = 0) { # replace small values with the next smallest value found, which is >X.
  newmin <- min(vec[vec > X])
  vec[vec <= X] <- newmin
  vec
}


# _________________________________________________________________________________________________
#' @title AreTheseCellNamesTheSame
#'
#' @description Compare two character vectors (e.g.: cell IDs) how much they overlap and plot a Venn Diagram.
#' @param vec1 Character vector, eg. with cell names
#' @param vec2 Character vector, eg. with cell names
#' @param names Names for plotting
#' @param min.overlap Threhold below there is no there is no meaninful overlap between the tow vectors. The function aborts with an error.
#'
#' @export
#' @examples # reTheseCellNamesTheSame()
AreTheseCellNamesTheSame <- function(
    vec1 = names(UVI.annot),
    vec2 = names(nr_UVI),
    names = c("Cells in Targ.Ampl", "Cells in GEX"),
    min.overlap = 0.33) {
  Cellname.Overlap <- list(vec1, vec2)
  names(Cellname.Overlap) <- if (!isFALSE(names)) names else c(substitute(vec1), substitute(vec2))

  cells.in.both <- intersect(vec1, vec2)
  sbb <- percentage_formatter(length(cells.in.both) / length(vec2), suffix = "of cells (GEX) in have a UVI assigned")
  ggExpress::qvenn(Cellname.Overlap, subt = sbb)
  iprint("Venn Diagramm saved.")
  iprint(sbb)

  Nr.overlapping <- length(intersect(vec1, vec2))
  Nr.total <- length(union(vec1, vec2))
  Percent_Overlapping <- Nr.overlapping / Nr.total
  print("")
  report <- percentage_formatter(Percent_Overlapping,
    prefix = "In total,",
    suffix = paste("of the cellIDs overlap across", names(Cellname.Overlap)[1], "and", names(Cellname.Overlap)[2])
  )
  print(report[1])
  stopifnot(Percent_Overlapping > min.overlap)
}



# _________________________________________________________________________________________________
#' Create.MiscSlot
#'
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
#' @title Add to Misc or Tools Slot
#'
#' @description This function adds a sub-slot to either the 'misc' or 'tools' slot of a Seurat object,
#' allowing for flexible data storage within the object structure.
#' If the sub-slot already exists, it can either be overwritten or a warning will be issued.
#'
#' @param obj A Seurat object.
#' @param slot_name The name of the slot to which the sub-slot should be added ('misc' or 'tools').
#' @param sub_slot_value The value to be assigned to the sub-slot.
#' @param sub_slot_name The name of the sub-slot. Automatically derived from 'sub_slot_value' if not provided.
#' @param overwrite A boolean indicating whether to overwrite an existing sub-slot with the same name.
#'
#' @return The modified Seurat object with the new or updated sub-slot.
#' @export
addToMiscOrToolsSlot <- function(obj, slot_name = 'misc', sub_slot_value = NULL,
                                 sub_slot_name = deparse(substitute(sub_slot_value)),
                                 overwrite = FALSE) {

  stopifnot(is(obj, "Seurat"), is.character(slot_name), length(slot_name) == 1)
  stopifnot(is.null(sub_slot_value) || !is.null(sub_slot_name),
            is.character(sub_slot_name), length(sub_slot_name) == 1)

  # Accessing the specified slot
  slot_orig <- slot(object = obj, name = slot_name)

  # Creating new slot or reporting if it exists
  if (sub_slot_name %in% names(slot_orig) && !overwrite) {
    warning(paste(sub_slot_name, "in", slot_name, "already exists. Not overwritten."), immediate. = TRUE)
  } else {
    slot_orig[[sub_slot_name]] <- sub_slot_value
  }

  # Assigning the modified slot back to the object
  slot(object = obj, name = slot_name) <- slot_orig

  return(obj)
}


# _________________________________________________________________________________________________
#' @title calc.q99.Expression.and.set.all.genes
#'
#' @description Calculate the gene expression of the e.g.: 90th quantile (expression in the top 10% cells). #
#' @param obj Seurat object, Default: combined.obj
#' @param quantileX Quantile level, Default: 0.9
#' @param max.cells Max number of cells to do the calculation on. Downsample if excdeeded. Default: 1e+05
#' @param slot slot in the Seurat object. Default: 'data'
#' @param assay RNA or integrated assay, Default: c("RNA", "integrated")[1]
#' @param set.misc Create the "all.genes" variable in @misc? Default: TRUE
#' @param assign_to_global_env Create the "all.genes" variable in the global env?, Default: TRUE
#' @param show Show plot? Default: TRUE
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   combined.obj <- calc.q99.Expression.and.set.all.genes(
#'     obj = combined.obj, quantileX = 0.9,
#'     max.cells = 25000
#'   )
#'   head(sort(as.numeric.wNames(obj@misc$expr.q90), decreasing = TRUE))
#'   combined.obj <- calc.q99.Expression.and.set.all.genes(
#'     obj = combined.obj, quantileX = 0.95,
#'     max.cells = 25000, set.all.genes = FALSE
#'   )
#' }
#' }
#' @seealso
#'  \code{\link[sparseMatrixStats]{character(0)}}
#' @importFrom tictoc tic toc
#' @importFrom sparseMatrixStats rowQuantiles
#'
#' @export

calc.q99.Expression.and.set.all.genes <- function(
    obj = combined.obj # Calculate the gene expression of the e.g.: 90th quantile (expression in the top 10% cells).
    , quantileX = 0.99, max.cells = 1e5
    , slot = "data"
    , assay = c("RNA", "integrated")[1]
    , set.misc = TRUE,
    assign_to_global_env = TRUE,
    show = TRUE) {
  tictoc::tic()
  x <- GetAssayData(object = obj, assay = assay, slot = slot) # , assay = 'RNA'
  if (ncol(x) > max.cells) {
    dsampled <- sample(x = 1:ncol(x), size = max.cells)
    x <- x[, dsampled]
  }
  qname <- paste0("q", quantileX * 100)
  slot_name <- kpp("expr", qname)

  # expr.q99 = iround(apply(x, 1, quantile, probs = quantileX) )
  print("Calculating Gene Quantiles")
  expr.q99.df <- sparseMatrixStats::rowQuantiles(x, probs = quantileX)
  expr.q99 <- iround(expr.q99.df)

  log2.gene.expr.of.the.90th.quantile <- as.numeric(log2(expr.q99 + 1)) # strip names
  n.cells <- floor(ncol(obj) * (1 - quantileX))
  qnameP <- paste0(100 * quantileX, "th quantile")
  try(
    ggExpress::qhistogram(log2.gene.expr.of.the.90th.quantile,
      ext = "pdf", breaks = 30,
      plotname = paste("Gene expression in the", qnameP),
      subtitle = kollapse(pc_TRUE(expr.q99 > 0, NumberAndPC = TRUE), " genes have ", qname, " expr. > 0."),
      caption = paste(n.cells, "cells in", qnameP),
      xlab = paste0("log2(expr. in the ", qnameP, "quantile+1) [UMI]"),
      ylab = "Nr. of genes",
      plot = show, save = TRUE,
      vline = .15,
      filtercol = TRUE,
      palette_use = "npg"
    ),
    silent = TRUE
  )

  all.genes <- percent_rank(expr.q99)
  names(all.genes) <- names(expr.q99)
  all.genes <- as.list(sort(all.genes, decreasing = TRUE))

  if (assign_to_global_env) assign("all.genes", all.genes, envir = as.environment(1))

  # if (set.all.genes) obj@misc$'all.genes' = all.genes
  if (set.misc) obj@misc[[slot_name]] <- expr.q99

  iprint("Quantile", quantileX, "is now stored under obj@misc$all.genes and $", slot_name, " Please execute all.genes <- obj@misc$all.genes.")
  return(obj)
}


# _________________________________________________________________________________________________
# Clustering ______________________________ ----
# _________________________________________________________________________________________________



# _________________________________________________________________________________________________
#' @title RenameClustering
#'
#' @description Rename clustering in a Seurat object.
#' @param namedVector named vector, where values = new, names(vec) = old
#' @param orig.ident meta.data colname original
#' @param suffix.new.ident How to name (suffix) the new identity. Default: "ManualNames"
#' @param new.ident meta.data colname new
#' @param suffix.plot Suffix description (short string) to be added to the umap plots.
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @param obj Seurat object
#'
#' @export

RenameClustering <- function(
    namedVector = ManualNames,
    orig.ident = "RNA_snn_res.0.3",
    suffix.new.ident = "ManualNames",
    new.ident = ppp(orig.ident, suffix.new.ident),
    obj = combined.obj,
    suffix.plot = "",
    plot_umaps = TRUE,
    ...) {
  NewX <- translate(
    vec = as.character(obj@meta.data[, orig.ident]),
    oldvalues = names(namedVector),
    newvalues = namedVector
  )

  obj <- AddMetaData(object = obj, metadata = NewX, col.name = new.ident)

  iprint("new.ident is", new.ident, "created from", orig.ident)
  print("")

  if (plot_umaps) {
    stopifnot(is.character(suffix.plot))
    suffix.plot <- make.names(suffix.plot)
    print(clUMAP(orig.ident, suffix = suffix.plot, sub = suffix.plot, obj = obj, ...))
    print(clUMAP(new.ident, suffix = suffix.plot, sub = suffix.plot, obj = obj, ...))
  } else {
    iprint("New ident:", new.ident)
  }

  return(obj)
}


# _________________________________________________________________________________________________
#' @title Shorten Clustering Names
#'
#' @description This function takes in a string representing a clustering name,
#' and shortens it according to specific rules. It replaces "snn_res." with "",
#' "best.matching.names" with "bmatch", "ordered" with "ord",
#' "ManualNames" with "mNames", and ".long" at the end of the string with ".L".
#'
#' @param str A character string representing the clustering name to be shortened.
#'
#' @return A character string representing the shortened clustering name.
#'
#' @examples
#' \dontrun{
#' shorten_clustering_names("RNA_snn_res.0.5.ordered.ManualNames") # Returns 'RNA.0.5.ord.mNames'
#' shorten_clustering_names("RNA_snn_res.0.3.best.matching.names.ManualNames.long") # Returns 'RNA.0.3.bmatch.mNames.L'
#' shorten_clustering_names("RNA_snn_res.1.7.ordered.ManualNames.Simplest") # Returns 'RNA.1.7.ord.mNames.Simplest'
#' shorten_clustering_names("RNA_snn_res.0.5.ordered.ManualNames.Simpler") # Returns 'RNA.0.5.ord.mNames.Simpler'
#' }
#'
#' @export
shorten_clustering_names <- function(str) {
  # Replace 'snn_res' with nothing
  str <- gsub("snn_res.", "", str)
  # Replace 'best.matching.names' with 'bmatch'
  str <- gsub("best.matching.names", "bmatch", str)
  # Replace 'ordered' with 'ord'
  str <- gsub("ordered", "ord", str)
  # Replace 'ManualNames' with 'mNames'
  str <- gsub("ManualNames", "mNames", str)
  # Replace 'long' with 'L'
  str <- gsub(".long$", ".L", str)
  return(str)
}



# _________________________________________________________________________________________________
#' getClusterNames
#'
#' @param obj Seurat object
#' @param ident ident
#' @export

getClusterNames <- function(obj = combined.obj, ident = GetClusteringRuns(obj)[2]) {
  print(GetClusteringRuns(obj))
  clz <- as.character(sort(deframe(unique(obj[[ident]]))))
  cat(dput(clz))
}




# _________________________________________________________________________________________________
#' @title GetClusteringRuns
#'
#' @description Get Clustering Runs: metadata column names #
#' @param obj Seurat object, Default: combined.obj
#' @param res Clustering resoluton to use, Default: FALSE
#' @param pat Pettern to match, Default: '*snn_res.*[0-9]$'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   GetClusteringRuns()
#' }
#' }
#' @export
GetClusteringRuns <- function(obj = combined.obj, res = FALSE, pat = "*snn_res.*[0-9]$") { # Get Clustering Runs: metadata column names
  if (res) pat <- gsub(x = pat, pattern = "\\[.*\\]", replacement = res)
  clustering.results <- CodeAndRoll2::grepv(x = colnames(obj@meta.data), pattern = pat)
  if (identical(clustering.results, character(0))) warning("No matching column found!", immediate. = TRUE)
  return(clustering.results)
}



# _________________________________________________________________________________________________
#' @title GetNamedClusteringRuns
#'
#' @description Get Clustering Runs: metadata column names
#' @param obj Seurat object, Default: combined.obj
#' @param res Clustering resoluton to use, Default: c(FALSE, 0.5)[1]
#' @param topgene Match clustering named after top expressed gene (see vertesy/Seurat.pipeline/~Diff gene expr.), Default: FALSE
#' @param pat Pettern to match, Default: '^cl.names.Known.*[0,1]\.[0-9]$'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   GetNamedClusteringRuns()
#' }
#' }
#' @export
GetNamedClusteringRuns <- function(
    obj = combined.obj
    , res = list(FALSE, 0.5)[[1]], topgene = FALSE,
    pat = c("^cl.names.Known.*[0,1]\\.[0-9]$", "Name|name")[2]) {
  if (res) pat <- gsub(x = pat, pattern = "\\[.*\\]", replacement = res)
  if (topgene) pat <- gsub(x = pat, pattern = "Known", replacement = "top")
  clustering.results <- CodeAndRoll2::grepv(x = colnames(obj@meta.data), pattern = pat)
  if (identical(clustering.results, character(0))) {
    warning("No matching column found! Trying GetClusteringRuns(..., pat = '*_res.*[0,1]\\.[0-9]$)"
            , immediate. = TRUE)
    clustering.results <- GetClusteringRuns(obj = obj, res = FALSE, pat = "*_res.*[0,1]\\.[0-9]$")
  }
  return(clustering.results)
}



# _________________________________________________________________________________________________
#' @title GetOrderedClusteringRuns
#'
#' @description Get Clustering Runs: metadata column names #
#' @param obj Seurat object, Default: combined.obj
#' @param res Clustering resoluton to use, Default: FALSE
#' @param pat Pettern to match, Default: '*snn_res.*[0,1]\.[0-9]\.ordered$'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   GetOrderedClusteringRuns()
#'   GetOrderedClusteringRuns(res = 0.5)
#' }
#' }
#' @export
GetOrderedClusteringRuns <- function(obj = combined.obj, res = FALSE, pat = "*snn_res.*[0,1]\\.[0-9]\\.ordered$") { # Get Clustering Runs: metadata column names
  if (res) pat <- gsub(x = pat, pattern = "\\[.*\\]", replacement = res)
  clustering.results <- CodeAndRoll2::grepv(x = colnames(obj@meta.data), pattern = pat)
  if (identical(clustering.results, character(0))) warning("No matching column found!", immediate. = TRUE)
  return(clustering.results)
}



# _________________________________________________________________________________________________
#' @title GetNumberOfClusters
#'
#' @description Get Number Of Clusters #
#' @param obj Seurat object, Default: combined.obj
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   GetNumberOfClusters()
#' }
#' }
#' @export
GetNumberOfClusters <- function(obj = combined.obj) { # Get Number Of Clusters
  clustering.results <- GetClusteringRuns(obj)
  print("## Number of clusters: ---------")
  for (cc in clustering.results) {
    NrCl <- length(unique(obj@meta.data[[cc]]))
    iprint(cc, "   ", NrCl)
  }
}


# _________________________________________________________________________________________________
#' @title calc.cluster.averages
#'
#' @description Calculates the average of a metadata column (numeric) per cluster.
#' @param col_name The name of the column for which the average is calculated. Default: 'Score.GO.0006096'.
#' @param plot.UMAP.too Whether to plot a UMAP as well. Default: TRUE.
#' @param return.plot Whether to return the plot. Default: FALSE.
#' @param obj The main Seurat object used for calculations. Default: combined.obj.
#' @param split_by Cluster to split by. Default: First entry of GetNamedClusteringRuns().
#' @param scale.zscore Whether to scale z-scores. Default: FALSE.
#' @param simplify Whether to simplify the result. Default: TRUE.
#' @param plotit Whether to plot the results. Default: TRUE.
#' @param histogram Whether to produce a histogram. Default: FALSE.
#' @param nbins The number of bins for the histogram. Default: 50.
#' @param suffix Suffix added to the filename. Default: NULL.
#' @param stat Statistical method applied, "mean" or "median". Default: "median".
#' @param quantile.thr The threshold for quantiles. Default: 0.9.
#' @param absolute.thr Absolute threshold used in computations. Default: FALSE.
#' @param filter The filter mode: 'above', 'below', or FALSE. Default: FALSE.
#' @param ylab.text Text for the y-axis label. Default: "Cluster" followed by the statistical method and "score".
#' @param title Title for the plot. Default: "Cluster" followed by the statistical method and column name.
#' @param subtitle The subtitle for the plot. Default: NULL.
#' @param width The width of the plot. Default: 8.
#' @param height The height of the plot. Default: 6.
#' @param ... Additional parameters passed to the internally called functions.
#' @param xlb The label for the x-axis. Default depends on the 'absolute.thr' parameter.
#' @param fname The filename for the plot. Default is based on column name and split_by value.
#' @export
#' @importFrom Stringendo percentage_formatter

calc.cluster.averages <- function(
    col_name = "Score.GO.0006096",
    plot.UMAP.too = TRUE,
    return.plot = FALSE,
    obj = combined.obj,
    split_by = GetNamedClusteringRuns()[1],
    scale.zscore = FALSE,
    simplify = TRUE, plotit = TRUE,
    histogram = FALSE, nbins = 50,
    suffix = NULL,
    stat = c("mean", "median")[2],
    quantile.thr = 0.9,
    absolute.thr = FALSE,
    filter = c(FALSE, "above", "below")[1],
    ylab.text = paste("Cluster", stat, "score"),
    title = paste("Cluster", stat, col_name),
    prefix.cl.names = FALSE,
    report = TRUE,
    subtitle = NULL,
    width = 8, height = 6,
    ...
    # , ylb = paste(ylab.text, col_name)
    # , xlb = paste("Clusters >",Stringendo::percentage_formatter(quantile.thr),"quantile are highlighted. |", split_by)
    , xlb = if (absolute.thr) {
      paste("Threshold at", absolute.thr)
    } else {
      paste(
        "Black lines: ", kppd(Stringendo::percentage_formatter(c(1 - quantile.thr, quantile.thr))), "quantiles |",
        "Cl. >", Stringendo::percentage_formatter(quantile.thr), "are highlighted. |", split_by
      )
    },
    fname = ppp(col_name, split_by, "cluster.average.barplot.pdf", ...)) { # calc.cluster.averages of a m
  iprint(substitute(obj), "split by", split_by)
  if (absolute.thr) iprint("In case of the absolute threshold, only the returned values are correct, the plot annotations are not!")

  if (plot.UMAP.too) qUMAP(obj = obj, feature = col_name)

  df.summary <-
    obj@meta.data %>%
    select_at(c(col_name, split_by)) %>%
    group_by_at(split_by) %>%
    summarize(
      "nr.cells" = n(),
      "mean" = mean(!!sym(col_name), na.rm = TRUE),
      "SEM" = sem(!!sym(col_name), na.rm = TRUE),
      "median" = median(!!sym(col_name), na.rm = TRUE),
      "SE.median" = 1.2533 * sem(!!sym(col_name), na.rm = TRUE)
    )

  if (simplify) {
    av.score <- df.summary[[stat]]
    names(av.score) <- if (!isFALSE(prefix.cl.names)) ppp("cl", df.summary[[1]]) else df.summary[[1]]
    av.score <- sortbyitsnames(av.score)
    if (scale.zscore) av.score <- (scale(av.score)[, 1])

    cutoff <- if (absolute.thr) absolute.thr else quantile(av.score, quantile.thr)
    cutoff.low <- if (absolute.thr) NULL else quantile(av.score, (1 - quantile.thr))

    iprint("quantile.thr:", quantile.thr)
    if (plotit) {
      if (histogram) {
        p <- ggExpress::qhistogram(
          vec = as.numeric(av.score), save = FALSE,
          vline = cutoff,
          plotname = ppp(title, quantile.thr),
          bins = nbins,
          subtitle = paste(subtitle, "| median in blue/dashed"),
          ylab = ylab.text,
          xlab = xlb # Abused
          , xlab.angle = 45
          # , ylim = c(-1,1)
          , ...
          # , ext = "png", w = 7, h = 5
        ) + geom_vline(xintercept = cutoff.low, lty = 2)
        print(p)
        title_ <- ppp(title, suffix, flag.nameiftrue(scale.zscore))
        ggExpress::qqSave(ggobj = p, title = title_, ext = "png", w = width, h = height)
      } else {
        p <- ggExpress::qbarplot(
          vec = av.score, save = FALSE,
          hline = cutoff,
          plotname = title,
          suffix = quantile.thr,
          subtitle = subtitle,
          ylab = ylab.text,
          xlab = xlb # Abused
          , xlab.angle = 45
          # , ylim = c(-1,1)
          , ...
          # , ext = "png", w = 7, h = 5
        ) + geom_hline(yintercept = cutoff.low, lty = 2)

        print(p)
        title_ <- ppp(title, suffix, flag.nameiftrue(scale.zscore))
        qqSave(ggobj = p, title = title_, fname = ppp(title_, split_by, "png"), w = width, h = height)
      }
    }

    if (report) print(paste0(col_name, ": ", paste(iround(av.score), collapse = " vs. ")))
    if (filter == "below") {
      return(filter_LP(av.score, threshold = cutoff, plot.hist = FALSE))
    } else if (filter == "above") {
      return(filter_HP(av.score, threshold = cutoff, plot.hist = FALSE))
    } else {
      return(av.score)
    }
  } else if (return.plot) { # if /not/ simplify
    return(p)
  } else {
    return(df.summary)
  }
}






# _________________________________________________________________________________________________
#' @title plot.expression.rank.q90
#'
#' @description Plot gene expression based on the expression at the 90th quantile (so you will not lose genes expressed in few cells).
#' @param obj Seurat object, Default: combined.obj
#' @param gene gene of interest, Default: 'ACTB'
#' @param filterZero Remove genes whose quantile-90 expression in 0? Default: TRUE
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   plot.expression.rank.q90(gene = "SATB2")
#' }
#' }
#' @importFrom Stringendo percentage_formatter
#' @importFrom MarkdownReports whist
#'
#' @export plot.expression.rank.q90
plot.expression.rank.q90 <- function(obj = combined.obj, gene = "ACTB", filterZero = TRUE) {
  expr.GOI <- obj@misc$expr.q90[gene]
  expr.all <- unlist(obj@misc$expr.q90)
  gene.found <- gene %in% names(expr.all)
  stopifnot(gene.found)

  if (expr.GOI == 0) iprint(gene, "is not expressed. q90-av.exp:", expr.GOI) else if (expr.GOI < 0.05) iprint(gene, "is lowly expressed. q90-av.exp:", expr.GOI)
  if (filterZero) {
    iprint("Zero 'q90 expression' genes (", pc_TRUE(expr.all == 0), ") are removed.")
    expr.all <- expr.all[expr.all > 0]
  }
  counts <- sum(obj@assays$RNA@counts[gene, ])
  if (expr.GOI == 0) {
    quantile.GOI <- 0
    title <- paste(gene, "is too lowly expressed: q90-av.exp is zero. \n There are", counts, "counts.")
  } else {
    pos.GOI <- which(names(expr.all) == gene)
    quantile.GOI <- ecdf(expr.all)(expr.all)[pos.GOI]
    title <- paste(gene, "is in the", Stringendo::percentage_formatter(quantile.GOI), "quantile of 'q90-av' expression. \n There are", counts, "counts")
  }
  suppressWarnings(
    MarkdownReports::whist(expr.all,
      vline = expr.GOI, breaks = 100, main = title, plotname = make.names(title),
      ylab = "Genes",
      xlab = "Av. mRNA in the 10% top expressing cells (q90 av.exp.)"
    )
  )
}


# _________________________________________________________________________________________________
# Interacting with the environment ______________________________ ----
# _________________________________________________________________________________________________
# Subsetting, downsampling and manipulating the Seurat object



# _________________________________________________________________________________________________
#' @title set.mm
#'
#' @description Helps to find metadata columns. It creates a list with the names of of 'obj@meta.data'.
#' @param obj Seurat object, Default: combined.obj
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   set.mm()
#'   mm
#' }
#' }
#' @export
set.mm <- function(obj = combined.obj) {
  mm <- CodeAndRoll2::list.fromNames(colnames(obj@meta.data))
  assign(x = "mm", value = mm, envir = as.environment(1))
}


# _________________________________________________________________________________________________
#' @title Get the First Seurat Object from a List of Seurat Objects
#'
#' @description #' If provided with a list of Seurat objects, this function returns the first
#' Seurat object in the list. If the input is a single Seurat object, it returns
#' the object itself. It is assumed that all elements of the list are Seurat
#' objects if the input is a list.
#'
#' @param obj A Seurat object, a list of Seurat objects, or any other list.
#'
#' @return The first Seurat object from the list or the Seurat object itself.
#' If the input is not a Seurat object or a list containing at least one Seurat
#' object, the function will throw an error.
#' @export
ww.get.1st.Seur.element <- function(obj) {
  if (is(obj)[1] == "list") {
    iprint("A list of objects is provided, taking the 1st from", length(obj), "elements.")
    obj <- obj[[1]]
  }
  stopifnot(is(obj) == "Seurat")
  return(obj)
}
# ww.get.1st.Seur.element(ls.Seurat[[1]])


# _________________________________________________________________________________________________
#' @title recall.all.genes
#'
#' @description all.genes set by calc.q99.Expression.and.set.all.genes() #
#' @param obj Seurat object, Default: combined.obj
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   recall.all.genes()
#'   all.genes
#' }
#' }
#' @importFrom MarkdownHelpers ww.assign_to_global
#'
#' @export
recall.all.genes <- function(obj = combined.obj) { # all.genes set by calc.q99.Expression.and.set.all.genes()
  obj <- ww.get.1st.Seur.element(obj)

  if ("all.genes" %in% names(obj@misc)) {
    if (!exists("all.genes")) {
      all.genes <- obj@misc$all.genes
      print(head(unlist(all.genes)))
      MarkdownHelpers::ww.assign_to_global(name = "all.genes", value = all.genes, verbose = FALSE)
    } else {
      print("  ->   Variable 'all.genes' exits in the global namespace.")
    }
  } else {
    print("  ->   Slot 'all.genes' does not exist in obj@misc.")
  }
}


# _________________________________________________________________________________________________
#' @title recall.meta.tags.n.datasets
#'
#' @description Recall  meta.tags from obj@misc to "meta.tags" in the global environment.
#' @param obj Seurat object, Default: combined.obj
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   recall.n.datasets()
#'   n.datasets
#' }
#' }
#' @importFrom MarkdownHelpers ww.assign_to_global
#'
#' @export
recall.meta.tags.n.datasets <- function(obj = combined.obj) {
  obj <- ww.get.1st.Seur.element(obj)

  if ("n.datasets" %in% names(obj@misc)) {
    if (!exists("n.datasets")) {
      n.datasets <- obj@misc$n.datasets
      print(head(unlist(n.datasets)))
      MarkdownHelpers::ww.assign_to_global(name = "n.datasets", value = n.datasets)
    } else {
      print("  ->   Variable 'n.datasets' already exists in the global namespace.")
    }
  } else {
    print("  ->   Slot 'n.datasets' does not exist in obj@misc.")
  }


  if ("meta.tags" %in% names(obj@misc)) {
    if (!exists("meta.tags")) {
      meta.tags <- obj@misc$meta.tags
      print(head(unlist(meta.tags)))
      MarkdownHelpers::ww.assign_to_global(name = "meta.tags", value = meta.tags)
    } else {
      iprint("  ->   Variable 'meta.tags' already exists in the global namespace.")
    }
  } else {
    print("  ->   Slot 'meta.tags' does not exist in obj@misc.")
  }
}


# _________________________________________________________________________________________________
#' @title recall.parameters
#'
#' @description Recall parameters from obj@misc to "p" in the global environment.
#' @param obj Seurat object, Default: combined.obj
#' @param overwrite Overwrite already existing in environment? Default: FALSE
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   recall.parameters()
#'   p
#' }
#' }
#' @importFrom MarkdownHelpers ww.assign_to_global
#'
#' @export
recall.parameters <- function(obj = combined.obj, overwrite = FALSE) {
  obj <- ww.get.1st.Seur.element(obj)

  if ("p" %in% names(obj@misc)) {
    if (!exists("p")) iprint("variable 'p' exits in the global namespace:", head(p))

    if (!exists("p") | (exists("p") & overwrite == TRUE)) {
      MarkdownHelpers::ww.assign_to_global(name = "p", value = obj@misc$"p")
      print("Overwritten.")
    } else {
      print("Not overwritten.")
    }
  } else {
    print("  ->   Slot 'p' does not exist in obj@misc.")
  }
}



# _________________________________________________________________________________________________
#' @title recall.genes.ls
#'
#' @description Recall genes.ls from obj@misc to "genes.ls" in the global environment.
#' @param obj Seurat object, Default: combined.obj
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   recall.genes.ls()
#'   genes.ls
#' }
#' }
#' @importFrom MarkdownHelpers ww.assign_to_global
#'
#' @export

recall.genes.ls <- function(obj = combined.obj, overwrite = FALSE) { # genes.ls
  obj <- ww.get.1st.Seur.element(obj)

  if ("genes.ls" %in% names(obj@misc)) {
    if (!exists("genes.ls")) iprint("variable 'genes.ls' exits in the global namespace:", head(p))

    if (!exists("genes.ls") | (exists("genes.ls") & overwrite == TRUE)) {
      MarkdownHelpers::ww.assign_to_global(name = "genes.ls", value = obj@misc$"genes.ls")
      print("Overwritten.")
    } else {
      print("Not overwritten.")
    }
  } else {
    print("  ->   Slot 'genes.ls' does not exist in obj@misc.")
  }
}



# _________________________________________________________________________________________________
#' @title save.parameters
#'
#' @description Save parameters to obj@misc$p
#' @param obj Seurat object, Default: combined.obj
#' @param params List of parameters, Default: p
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   save.parameters(obj = combined.obj, params = p)
#' }
#' }
#' @export
save.parameters <- function(obj = combined.obj, params = p, overwrite = TRUE) {
  obj <- ww.get.1st.Seur.element(obj)

  if (!is.null(obj@misc$"p") && overwrite) {
    print("Overwriting already existing obj@misc$p. Old version:")
    print(head(unlist(obj@misc$"p")))
    obj@misc$p <- params
  } else if (is.null(obj@misc$"p")) {
    obj@misc$p <- params
  }
}



# _________________________________________________________________________________________________
# List level metadata for ______________________________ ----
# _________________________________________________________________________________________________


#' @title Create Single-Cell Metadata Object for a collection of Seurat Objects
#'
#' @description This function creates a metadata object to correspond to a list of
#'   single-cell experiments, for storing parent level information.
#'   It initializes the object with the experiment and project name, and the
#'   creation date. The created object is of class 'scMetadata_class'.
#' @param experiment The name of the experiment for which metadata is being created.
#' @param project_ The project information to be associated with the metadata object.
#'   This defaults to the current project obtained using Seurat.utils::getProject().
#'
#' @return An 'scCollectionMetadata_class' object containing the metadata for a collection of experiment.
#' @export
#'
#' @examples
#' sc_meta <- create_scCombinedMeta(experiment = "Experiment1")
create_scCombinedMeta <- function(experiment, project_ = getProject()) {
  x <- list(
    experiment.corresponding = experiment,
    initialized = format(Sys.time(), format = "%Y.%m.%d | %H:%M:%S"),
    project = project_
  )
  class(x) <- "scCollectionMetadata_class"
  print(x)
  return(x)
}





# _________________________________________________________________________________________________
# Subsetting the Seurat object ______________________________ ----
# _________________________________________________________________________________________________
# Subsetting, downsampling and manipulating the Seurat object


# _________________________________________________________________________________________________
#' @title subsetSeuObj
#'
#' @description Subset a compressed Seurat object and save it in the working directory.
#' @param obj A Seurat object to subset. Default: the i-th element of the list 'ls.Seurat'.
#' @param fraction_ The fraction of the object's data to keep. Default: 0.25.
#' @param nCells If set to a number greater than 1, indicates the absolute number of cells to keep. If FALSE, the function uses 'fraction_' to determine the number of cells. Default: FALSE.
#' @param seed_ A seed for random number generation to ensure reproducible results. Default: 1989.
#' @export
#' @importFrom Stringendo percentage_formatter

subsetSeuObj <- function(obj = ls.Seurat[[i]], fraction_ = 0.25, nCells = FALSE, seed_ = 1989) { # Subset a compressed Seurat Obj and save it in wd.
  set.seed(seed_)
  if (isFALSE(nCells)) {
    cellIDs.keep <- sampleNpc(metaDF = obj@meta.data, pc = fraction_)
    iprint(length(cellIDs.keep), "or", Stringendo::percentage_formatter(fraction_), "of the cells are kept. Seed:", seed_)
  } else if (nCells > 1) {
    nKeep <- min(ncol(obj), nCells)
    # print(nKeep)
    cellIDs.keep <- sample(colnames(obj), size = nKeep, replace = FALSE)
    if (nKeep < nCells) iprint("Only", nCells, "cells were found in the object, so downsampling is not possible.")
  }
  obj <- subset(x = obj, cells = cellIDs.keep) # downsample
  return(obj)
}

# _________________________________________________________________________________________________
#' @title subsetSeuObj.and.Save
#'
#' @description Subset a compressed Seurat Obj and save it in wd. #
#' @param obj Seurat object, Default: ORC
#' @param fraction Fractional size to downsample to. Default: 0.25
#' @param seed random seed used, Default: 1989
#' @param min.features Minimum features
#' @param dir Directory to save to. Default: OutDir
#' @param suffix A suffix added to the filename, Default: ''
#' @export
subsetSeuObj.and.Save <- function(
    obj = ORC, fraction = 0.25, seed = 1989, dir = OutDir,
    min.features = p$"min.features", suffix = "") { # Subset a compressed Seurat Obj and save it in wd.
  obj_Xpc <- subsetSeuObj(obj = obj, fraction_ = fraction, seed_ = seed) # , nthreads = 6
  nr.cells.kept <- ncol(obj_Xpc)
  # Seurat.utils:::.saveRDS.compress.in.BG(obj = obj_Xpc, fname = ppp(paste0(dir, substitute(obj)), suffix, nr.cells.kept, 'cells.with.min.features', min.features,"Rds" ) )
  xsave(obj_Xpc,
    suffix = ppp(suffix, nr.cells.kept, "cells.with.min.features", min.features),
    nthreads = nthreads, project = getProject(), showMemObject = TRUE, saveParams = FALSE
  )
}


# _________________________________________________________________________________________________
#' @title subsetSeuObj.ident.class
#'
#' @description Subset a Seurat Obj to a given column
#' @param obj Seurat object, Default: ORC
#' @param ident identity
#' @param clusters which value to match
#' @param invert invert selecion
#' @export

subsetSeuObj.ident.class <- function(
    obj = combined.obj, ident = "RNA_snn_res.0.5.ordered.ManualNames",
    clusters = "Neuron, unclear", invert = FALSE) {
  Idents(obj) <- ident
  cellz <- WhichCells(obj, idents = clusters, invert = invert)
  iprint(length(cellz), "cells are selected from", ncol(obj), "using", ident)
  subset(x = obj, cells = cellz)
}


# _________________________________________________________________________________________________
#' @title Downsample.Seurat.Objects
#'
#' @description Downsample a list of Seurat objects
#' @param ls.obj List of Seurat objects. Default: ls.Seurat
#' @param NrCells Number of cells to downsample to.
#' @param save_object save object by isaveRDS, otherwise return it. Default: FALSE
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   Downsample.Seurat.Objects(NrCells = 2000)
#'   Downsample.Seurat.Objects(NrCells = 200)
#' }
#' }
#' @export
#' @importFrom tictoc tic toc
#' @importFrom Stringendo percentage_formatter
#' @importFrom foreach foreach getDoParRegistered
Downsample.Seurat.Objects <- function(
    ls.obj = ls.Seurat, NrCells = p$"dSample.Organoids",
    save_object = FALSE) {
  # Check if 'ls_obj' is a list of Seurat objects and 'obj_IDs' is a character vector of the same length
  if (!is.list(ls.obj) & inherits(ls.obj, "Seurat")) ls.obj <- list(ls.obj)
  stopifnot(is.list(ls.obj) & all(sapply(ls.obj, function(x) inherits(x, "Seurat"))))

  names.ls <- names(ls.obj)
  n.datasets <- length(ls.obj)
  iprint(NrCells, "cells")
  tictoc::tic()
  if (foreach::getDoParRegistered()) {
    ls.obj.downsampled <- foreach::foreach(i = 1:n.datasets) %dopar% {
      iprint(names(ls.obj)[i], Stringendo::percentage_formatter(i / n.datasets, digitz = 2))
      subsetSeuObj(obj = ls.obj[[i]], nCells = NrCells)
    }
    names(ls.obj.downsampled) <- names.ls
  } else {
    ls.obj.downsampled <- list.fromNames(names.ls)
    for (i in 1:n.datasets) {
      iprint(names(ls.obj)[i], Stringendo::percentage_formatter(i / n.datasets, digitz = 2))
      ls.obj.downsampled[[i]] <- subsetSeuObj(obj = ls.obj[[i]], nCells = NrCells)
    }
  } # else
  tictoc::toc()

  print(head(unlapply(ls.obj, ncol)))
  print(head(unlapply(ls.obj.downsampled, ncol)))

  if (save_object) {
    isave.RDS(obj = ls.obj.downsampled, suffix = ppp(NrCells, "cells"), inOutDir = TRUE)
  } else {
    return(ls.obj.downsampled)
  }
}



# _________________________________________________________________________________________________
#' @title Downsample.Seurat.Objects.PC
#'
#' @description Downsample a list of Seurat objects, by fraction
#' @param ls.obj List of Seurat objects, Default: ls.Seurat
#' @param NrCells Number of cells to downsample to. Default: p$dSample.Organoids
#' @param save_object save object by isaveRDS, otherwise return it. Default: FALSE
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   Downsample.Seurat.Objects.PC()
#' }
#' }
#' @importFrom tictoc tic toc
#' @importFrom Stringendo percentage_formatter
#' @importFrom foreach getDoParRegistered foreach
#'
#' @export
Downsample.Seurat.Objects.PC <- function(
    ls.obj = ls.Seurat, fraction = 0.1,
    save_object = FALSE) {
  # Check if 'ls_obj' is a list of Seurat objects and 'obj_IDs' is a character vector of the same length
  if (!is.list(ls.obj) & inherits(ls.obj, "Seurat")) ls.obj <- list(ls.obj)
  stopifnot(is.list(ls.obj) & all(sapply(ls.obj, function(x) inherits(x, "Seurat"))))

  names.ls <- names(ls.obj)
  n.datasets <- length(ls.obj)
  iprint(fraction, "fraction")

  tictoc::tic()
  if (foreach::getDoParRegistered()) {
    ls.obj.downsampled <- foreach::foreach(i = 1:n.datasets) %dopar% {
      subsetSeuObj(obj = ls.obj[[i]], fraction_ = fraction)
    }
    names(ls.obj.downsampled) <- names.ls
  } else {
    ls.obj.downsampled <- list.fromNames(names.ls)
    for (i in 1:n.datasets) {
      cells <- round(ncol(ls.obj[[1]]) * fraction)
      iprint(names(ls.obj)[i], cells, "cells=", Stringendo::percentage_formatter(i / n.datasets, digitz = 2))
      ls.obj.downsampled[[i]] <- subsetSeuObj(obj = ls.obj[[i]], fraction_ = fraction)
    }
  }
  tictoc::toc() # else

  NrCells <- sum(unlapply(ls.obj, ncol))

  print(head(unlapply(ls.obj, ncol)))
  print(head(unlapply(ls.obj.downsampled, ncol)))
  if (save_object) {
    isave.RDS(obj = ls.obj.downsampled, suffix = ppp(NrCells, "cells"), inOutDir = TRUE)
  } else {
    return(ls.obj.downsampled)
  }
}


# _________________________________________________________________________________________________
#' @title remove.residual.small.clusters
#'
#' @description E.g.: after subsetting often some residual cells remain in clusters originally defined in the full dataset.
#' @param identitites Identities to scan for residual clusters
#' @param obj Seurat object, Default: combined.obj
#' @param max.cells Max number of cells in cluster to be removed. Default: 0.5% of the dataset, or 5 cells.
#' @export
remove.residual.small.clusters <- function(
    obj = combined.obj,
    identitites = GetClusteringRuns(obj),
    max.cells = max(round((ncol(obj)) / 2000), 5)) {
  META <- obj@meta.data
  all.cells <- rownames(META)

  iprint("max.cells:", max.cells, "| Scanning over these", length(identitites), "identities:", identitites)
  small.clusters <- cells.to.remove <- list.fromNames(identitites)

  for (i in 1:length(identitites)) {
    colX <- identitites[i]
    tbl <- table(META[[colX]])

    small.clusters[[i]] <- which_names(tbl <= max.cells)
    cells.to.remove[[i]] <- all.cells[which(META[[colX]] %in% small.clusters[[i]])]
    if (length(cells.to.remove[[i]])) {
      iprint(
        length(cells.to.remove[[i]]), "cells in small clusters:", small.clusters[[i]],
        "| Cell counts:", tbl[small.clusters[[i]]]
      )
    }
  }

  all.cells.2.remove <- unique(unlist(cells.to.remove))
  if (length(all.cells.2.remove)) {
    iprint(">>> a total of", length(all.cells.2.remove), "cells are removed which belonged to a small cluster in any of the identities.")
  } else {
    iprint(">>> No cells are removed because belonging to small cluster.")
  }

  cells.2.keep <- setdiff(all.cells, all.cells.2.remove)
  obj <- subset(x = obj, cells = cells.2.keep)

  return(obj)
}


# _________________________________________________________________________________________________
#' @title dropLevelsSeurat
#'
#' @description Drop unused levels from factor variables in a Seurat object.
#' @param obj A Seurat object.
#' @param verbose Logical. Whether to print a message indicating which levels are being dropped.
#' @export
dropLevelsSeurat <- function(obj = combined.obj, verbose = TRUE) {
  META <- obj@meta.data
  colclasses <- sapply(META, class)
  drop_in_these <- names(colclasses[colclasses == "factor"])
  if (verbose) iprint("Dropping levels in", length(drop_in_these), drop_in_these)
  for (i in 1:length(drop_in_these)) {
    colX <- drop_in_these[i]
    META[[colX]] <- droplevels(META[[colX]])
  }

  obj@meta.data <- META
  return(obj)
}



# ____________________________________________________________________
#' @title Remove Clusters and Drop Levels
#'
#' @description This function removes residual small clusters from specified Seurat objects and drops levels in factor-like metadata.
#' @param ls_obj A list of Seurat objects.
#' @param object_names A character vector containing the names of the Seurat objects to process. Default is names of all objects in the `ls_obj`.
#' @param indices A numeric vector indicating which datasets to process by their position in the `object_names` vector. By default, it processes the second and third datasets.
#' @param ... Additional parameters passed to the `remove.residual.small.clusters` function.
#'
#' @details This function applies `remove.residual.small.clusters` and `dropLevelsSeurat` to the Seurat objects specified by the `indices` in the `object_names`.
#' It operates in place, modifying the input `ls_obj` list.
#'
#' @return The function returns the modified list of Seurat objects.
#' @examples
#' \dontrun{
#' # Process the 2nd and 3rd datasets
#' remove_clusters_and_drop_levels(ls_obj, indices = c(2, 3))
#' }
#'
#' @export
remove_clusters_and_drop_levels <- function(
    ls_obj, object_names = names(ls_obj),
    indices = 2:3, ...) {
  for (index in indices) {
    dataset_name <- object_names[index]
    obj <- ls_obj[[dataset_name]]
    obj <- remove.residual.small.clusters(obj = obj, identitites = GetClusteringRuns(obj), ...)
    obj <- dropLevelsSeurat(obj)
    ls_obj[[dataset_name]] <- obj
  }
  return(ls_obj)
}




# _________________________________________________________________________________________________
#' @title Remove Cells by Dimension Reduction
#'
#' @description This function applies a cutoff in the specified dimension of a given dimension reduction (UMAP, PCA, or t-SNE) to remove cells.
#' @param reduction A string specifying the dimension reduction technique to be used ('umap', 'pca', or 'tsne'). Default is 'umap'.
#' @param umap_dim An integer specifying which dimension (axis) to apply the cutoff. Default is 1.
#' @param obj A Seurat object. Default is 'combined.obj'.
#' @param cutoff A numerical value indicating the cutoff value for the specified dimension. Default is 0.
#' @param cut_below A logical value indicating whether to remove cells below (TRUE) or above (FALSE) the cutoff line. Default is TRUE.
#' @param only_plot_cutoff Simulate and plot cutoff only.
#' @param ... Any other parameters to be passed to internally called functions.
#' @return A Seurat object with cells removed according to the specified cutoff.
#' @export
remove.cells.by.UMAP <- function(
    reduction = "umap",
    umap_dim = 1,
    obj = combined.obj,
    cutoff = 0,
    cut_below = TRUE,
    only_plot_cutoff = FALSE,
    ...) {
  # Plot cells
  sfx <- if (cut_below) "below" else "above"
  p <- clUMAP(obj = obj, save.plot = FALSE, sub = paste0("cutoff ", sfx, ": ", cutoff), ...)

  # Add cutoff line to plot
  if (umap_dim == 1) {
    p <- p + geom_vline(xintercept = cutoff)
  } else if (umap_dim == 2) {
    p <- p + geom_hline(yintercept = cutoff)
  }
  print(p)
  qqSave(p, fname = kpp("UMAP.with.cutoff", umap_dim, sfx, cutoff, "png"), h = 7, w = 7)

  if (!only_plot_cutoff) {
    # Retrieve cell embeddings
    cell_embedding <- obj@reductions[[reduction]]@cell.embeddings
    all_cells <- rownames(cell_embedding)
    stopifnot(ncol(cell_embedding) > 0)
    embedding_dim_x <- cell_embedding[, umap_dim]

    # Determine cells to remove based on cutoff
    cells_to_remove <- if (cut_below) which_names(embedding_dim_x < cutoff) else which_names(embedding_dim_x >= cutoff)

    # Report on cells removed
    if (length(cells_to_remove)) {
      iprint(">>> A total of", length(cells_to_remove), "cells are removed which fell on UMAP aside cutoff:", cutoff)
    } else {
      iprint(">>> No cells are removed because of the UMAP dimension cutoff.")
    }

    # Subset object to include only cells not removed
    cells_to_keep <- setdiff(all_cells, cells_to_remove)
    obj <- subset(x = obj, cells = cells_to_keep)
  } # only_plot_cutoff
  return(obj)
}



# _________________________________________________________________________________________________
# DGEA ______________________________ ----
# _________________________________________________________________________________________________


# _________________________________________________________________________________________________
#' @title Add.DE.combined.score
#'
#' @description Add a combined score to differential expression (DE) results. The score is calculated as log-fold change (LFC) times negative logarithm of scaled p-value (LFC * -log10( p_cutoff / pval_scaling )).
#' @param df A data frame that holds the result of a differential gene expression analysis, typically obtained via the 'FindAllMarkers' function. Default: df.markers.
#' @param p_val_min The minimum p-value considered. All values below this threshold are set to this value. Default: 1e-25.
#' @param pval_scaling The value to scale p-values by in the calculation of the combined score. Default: 0.001.
#' @param colP The name of the column in the input data frame that holds p-values. Default: 'p_val'.
#' @param colLFC The name of the column in the input data frame that holds log-fold change values. By default, it selects the first column not named "avg_logFC" or "avg_log2FC".
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   df.markers <- Add.DE.combined.score(df.markers)
#' }
#' }
#' @export
Add.DE.combined.score <- function(
    df = df.markers, p_val_min = 1e-25, pval_scaling = 0.001, colP = "p_val",
    colLFC = CodeAndRoll2::grepv(pattern = c("avg_logFC|avg_log2FC"), x = colnames(df), perl = TRUE)
    # , colLFC = "avg_log2FC"
    ) { # Score = -LOG10(p_val) * avg_log2FC
  p_cutoff <- SmallestNonAboveX(vec = df[[colP]], X = p_val_min)
  df$"combined.score" <- round(df[[colLFC]] * -log10(p_cutoff / pval_scaling))
  return(df)
}




# _________________________________________________________________________________________________
#' @title StoreTop25Markers
#'
#' @description Save the top 25 makers based on `avg_log2FC` output table of `FindAllMarkers()` (df_markers) under `@misc$df.markers$res...`. By default, it rounds up insignificant digits up to 3. #
#' @param obj Seurat object, Default: combined.obj
#' @param df_markers Data frame, result of DGEA analysis (FindAllMarkers), Default: df.markers
#' @param res Clustering resoluton to use, Default: 0.5
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   combined.obj <- StoreTop25Markers(df_markers = df.markers, res = 0.)
#' }
#' }
#' @seealso
#'  \code{\link[dplyr]{select}}
#' @export
#' @importFrom dplyr select
StoreTop25Markers <- function(obj = combined.obj # Save the top 25 makers based on `avg_log2FC` output table of `FindAllMarkers()` (df_markers) under `@misc$df.markers$res...`. By default, it rounds up insignificant digits up to 3.
                              , df_markers = df.markers, res = 0.5) {
  top25.markers <-
    df_markers %>%
    group_by(cluster) %>%
    top_n(n = 25, wt = avg_2logFC) %>%
    dplyr::select(gene) %>%
    col2named.vec.tbl() %>%
    splitbyitsnames()

  obj@misc$"top25.markers"[[ppp("res", res)]] <- top25.markers
  return(obj)
}


# _________________________________________________________________________________________________
#' @title StoreAllMarkers
#'
#' @description Save the output table of `FindAllMarkers()` (df_markers) under
#' `@misc$df.markers$res...`. By default, it rounds up insignificant digits up to 3.
#' @param obj Seurat object, Default: combined.obj
#' @param df_markers Data frame, result of DGEA analysis (FindAllMarkers), Default: df.markers
#' @param res Clustering resoluton to use, Default: 0.5
#' @param digit Number of digits to keep, Default: c(0, 3)[2]
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   combined.obj <- StoreAllMarkers(df_markers = df.markers, res = 0.5)
#' }
#' }
#' @export
StoreAllMarkers <- function(
    obj = combined.obj,
    df_markers = df.markers, res = 0.5, digit = c(0, 3)[2]) {
  if (digit) df_markers[, 1:5] <- signif(df_markers[, 1:5], digits = digit)
  obj@misc$"df.markers"[[ppp("res", res)]] <- df_markers
  iprint("DF markers are stored under:", "obj@misc$df.markers$", ppp("res", res))
  return(obj)
}


# _________________________________________________________________________________________________
#' @title GetTopMarkersDF
#'
#' @description Get the vector of N most diff. exp. genes. #
#' @param dfDE Data frame, result of DGEA analysis (FindAllMarkers), Default: df.markers
#' @param n Number of markers to return, Default: p$n.markers
#' @param order.by Sort output tibble by which column, Default: c("avg_log2FC", "p_val_adj")[1]
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   GetTopMarkers(df = df.markers, n = 3)
#' }
#' }
#' @seealso
#'  \code{\link[dplyr]{slice}}, \code{\link[dplyr]{select}}
#' @export
#' @importFrom dplyr slice select
GetTopMarkersDF <- function(dfDE = df.markers # Get the vector of N most diff. exp. genes.
                            , n = p$"n.markers", order.by = c("avg_log2FC", "p_val_adj")[1]
                            , exclude = c("^AL*|^AC*|^LINC*")
                            ) {
  "Works on active Idents() -> thus we call cluster"

  library(dplyr)

  # browser()
  TopMarkers <- dfDE |>
    filter(!grepl(exclude, gene, perl = TRUE)) |>
    arrange(desc(!!as.name(order.by))) |>
    dplyr::group_by(cluster) |>
    dplyr::slice(1:n) |>
    dplyr::select(cluster, gene, avg_log2FC)

  # TopMarkers <- dfDE %>%
  #   arrange(desc(!!as.name(order.by))) %>%
  #   group_by(cluster) %>%
  #   dplyr::slice(1:n) %>%
  #   dplyr::select(gene)

  return(TopMarkers)
}


# _________________________________________________________________________________________________
#' @title GetTopMarkers
#'
#' @description Get the vector of N most diff. exp. genes. #
#' @param dfDE Data frame, result of DGEA analysis (FindAllMarkers), Default: df.markers
#' @param n Number of markers to return, Default: p$n.markers
#' @param order.by Sort output tibble by which column, Default: c("combined.score", "avg_log2FC", "p_val_adj")[2]
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   GetTopMarkers(df = df.markers, n = 3)
#' }
#' }
#' @seealso
#'  \code{\link[dplyr]{slice}}, \code{\link[dplyr]{select}}
#' @export
#' @importFrom dplyr slice select
GetTopMarkers <- function(dfDE = df.markers # Get the vector of N most diff. exp. genes.
                          , n = p$"n.markers", order.by = c("combined.score", "avg_log2FC", "p_val_adj")[2]) {
  "Works on active Idents() -> thus we call cluster"
  TopMarkers <- dfDE %>%
    arrange(desc(!!as.name(order.by))) %>%
    group_by(cluster) %>%
    dplyr::slice(1:n) %>%
    dplyr::select(gene) %>%
    col2named.vec.tbl()

  return(TopMarkers)
}




# _________________________________________________________________________________________________
#' @title AutoLabelTop.logFC
#'
#' @description Create a new "named identity" column in the metadata of a Seurat object,
#' with `Ident` set to a clustering output matching the `res` parameter of the function.
#' It requires the output table of `FindAllMarkers()`.  If you used `StoreAllMarkers()`
#' is stored under `@misc$df.markers$res...`, which location is assumed by default.
#' @param obj A Seurat object, with default value `combined.obj`.
#' @param group.by The clustering group to be used, defaults to the first entry by
#' `GetClusteringRuns()`.
#' @param res Clustering resolution to use, with a default value of 0.1.
#' @param plot.top.genes Logical indicating whether to show a plot, default is `TRUE`.
#' @param suffix Suffix for the naming, defaults to the value of `res`.
#' @param order.by Sorting criterion for the output tibble, defaults to the second element
#' of `c("combined.score", "avg_log2FC", "p_val_adj")`.
#' @param exclude A vector of regular expressions to specify genes to exclude, with
#' default value `c("^AL*|^AC*|^LINC*")`.
#' @param df_markers Data frame resulting from DGEA analysis (`FindAllMarkers`). The default
#' is `combined.obj@misc$df.markers[[paste0("res.", res)]]`.
#' @param plotEnrichment Logical indicating whether to plot enrichment, default is `TRUE`.
#'
#' @examples
#' \dontrun{
#'   if (interactive()) {
#'     combined.obj <- AutoLabelTop.logFC()
#'     combined.obj$"cl.names.top.gene.res.0.5"
#'   }
#' }

#' @export
AutoLabelTop.logFC <- function(
    obj = combined.obj,
    group.by = GetClusteringRuns(obj)[1],
    res = 0.1, plot.top.genes = TRUE,
    suffix = res,
    order.by = c("combined.score", "avg_log2FC", "p_val_adj")[2],
    exclude = c("^AL*|^AC*|^LINC*"),
    df_markers = obj@misc$"df.markers"[[paste0("res.", res)]],
    plotEnrichment = TRUE) {

  stopifnot(!is.null("df_markers"),
            order.by %in% colnames(df_markers)
  )

  df.top.markers <- GetTopMarkersDF(dfDE = df_markers, order.by = order.by, n = 1, exclude = exclude)

  if (plotEnrichment) {
    top_log2FC <- df.top.markers$"avg_log2FC"
    names(top_log2FC) <- ppp(df.top.markers$"cluster", df.top.markers$"gene")
    ggExpress::qbarplot(top_log2FC, label = iround(top_log2FC)
                        , subtitle = suffix
                        , ylab = "avg_log2FC", xlab = "clusters"
                        , suffix = suffix )
  }

  top.markers <- col2named.vec.tbl(df.top.markers[, 1:2])

  obj@misc[[ppp("top.markers.res", res)]] <- top.markers

  ids_CBC <- deframe(obj[[group.by]])
  ids <- unique(ids_CBC)

  if (length(ids) != length(top.markers)) {
    warning("Not all clusters returned DE-genes!", immediate. = TRUE)
    missing <- setdiff(ids, names(top.markers))
    names(missing) <- missing
    iprint("missing:", missing)
    top.markers <- sortbyitsnames(c(top.markers, missing))
  }

  top.markers.ID <- ppp(names(top.markers), top.markers)
  names(top.markers.ID) <- names(top.markers)
  named.group.by <- top.markers.ID[ids_CBC]

  namedIDslot <- ppp("cl.names.top.gene.res", res)
  obj <- addMetaDataSafe(obj, metadata = as.character(named.group.by), col.name = namedIDslot, overwrite = TRUE)
  if (plot.top.genes) multiFeaturePlot.A4(list.of.genes = top.markers, suffix = suffix)

  return(obj)
}







# _________________________________________________________________________________________________
#' @title AutoLabel.KnownMarkers
#'
#' @description Creates a new "named identity" column in the metadata of a Seurat object,
#'  setting 'Ident' to a clustering output matching the 'res' parameter.
#'  This function requires the output table of `FindAllMarkers()`.
#' If you used `StoreAllMarkers()`, the output is stored under `@misc$df.markers$res...`,
#' which is the default location.
#' @param obj A Seurat object to work with. Default: combined.obj.
#' @param topN The top 'N' genes to consider. Default: 1.
#' @param res The clustering resolution to use. Default: 0.5.
#' @param KnownMarkers A character vector containing known marker genes to be used for annotation.
#' Default: c("TOP2A", "EOMES", "SLA", "HOPX", "S100B", "DLX6-AS1", "POU5F1", "SALL4", "DDIT4",
#' "PDK1", "SATB2", "FEZF2").
#' @param order.by Specifies the column to sort the output tibble by.
#' Default: 'combined.score' (First among "combined.score", "avg_log2FC", "p_val_adj").
#' @param df_markers The data frame of markers. By default, it is stored under
#'  `@misc$df.markers$res...` in the provided Seurat object.
#'  Default: combined.obj@misc$df.markers[[paste0("res.", res)]].
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   combined.obj <- AutoLabel.KnownMarkers()
#'   DimPlot.ClusterNames(ident = "cl.names.KnownMarkers.0.5")
#' }
#' }
#' @seealso
#'  \code{\link[dplyr]{select}}, \code{\link[dplyr]{slice}}
#' @export
#' @importFrom dplyr select slice
AutoLabel.KnownMarkers <- function(
    obj = combined.obj, topN = 1, res = 0.5,
    KnownMarkers = c(
      "TOP2A", "EOMES", "SLA", "HOPX", "S100B",
      "DLX6-AS1", "POU5F1", "SALL4", "DDIT4", "PDK1",
      "SATB2", "FEZF2"
    ),
    order.by = c("combined.score", "avg_log2FC", "p_val_adj")[1],
    df_markers = obj@misc$"df.markers"[[paste0("res.", res)]]) {
  stopifnot(!is.null("df_markers"))

  lfcCOL <- CodeAndRoll2::grepv(pattern = c("avg_logFC|avg_log2FC"), x = colnames(df_markers), perl = TRUE)
  keep <- unique(c(lfcCOL, "p_val_adj", "cluster", order.by, "gene"))


  matching.clusters <-
    df_markers %>%
    dplyr::select(keep) %>%
    arrange(desc(!!as.name(order.by))) %>%
    filter(gene %in% KnownMarkers) %>%
    group_by(gene) %>%
    dplyr::slice(1:topN) %>%
    arrange(desc(!!as.name(order.by))) %>%
    # top_n(n = 1, wt = avg_log2FC) %>% # Select the top cluster for each gene
    arrange(cluster)

  print(matching.clusters)

  unique.matches <-
    matching.clusters %>%
    group_by(cluster) %>% # Select rows with unique values based on column "cluster"
    distinct(cluster, .keep_all = TRUE) %>%
    dplyr::select(gene)

  print("Best matches:")
  print(unique.matches)

  top.markers.df <- GetTopMarkersDF(dfDE = df_markers, order.by = lfcCOL, n = 1)
  top.markers <- top.markers.df %>% col2named.vec.tbl()

  missing.annotations <-
    top.markers.df %>%
    filter(!cluster %in% unique.matches$cluster) # filter for clusters that do not have a unique label already

  named.annotations <-
    rbind(unique.matches, missing.annotations) %>% # merge the 2 df's
    arrange(cluster) %>%
    col2named.vec.tbl() # requires github.com/vertesy/CodeAndRoll

  (top.markers.ID <- ppp(names(named.annotations), named.annotations))
  names(top.markers.ID) <- names(top.markers)
  named.ident <- top.markers.ID[Idents(object = obj)]

  namedIDslot <- ppp("cl.names.KnownMarkers", res)
  obj[[namedIDslot]] <- named.ident
  return(obj)
}



# ________________________________________________________________________
#' @title scEnhancedVolcano
#'
#' @description This function creates an enhanced volcano plot.
#'
#' @param toptable A data frame with the results of differential gene expression analysis.
#' @param lab A vector of gene symbols to label on the plot.
#' @param suffix A string to append to the title of the plot.
#' @param title The title of the plot.
#' @param subtitle The subtitle of the plot.
#' @param x The x-axis, which is typically the average log2 fold change.
#' @param y The y-axis, which is typically the adjusted p-value.
#' @param selectLab A vector of gene symbols to select for labeling.
#' @param h The height of the plot.
#' @param w The width of the plot.
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @return A ggplot object.
#' @importFrom EnhancedVolcano EnhancedVolcano
#'
#' @export

scEnhancedVolcano <- function(
    toptable = df.markers.X, lab = rownames(toptable),
    suffix = "",
    title = kppd("DGEA", suffix),
    subtitle = paste("min FC:", iround(2^abs(min(toptable$"avg_log2FC")))),
    x = "avg_log2FC", y = "p_val_adj",
    selectLab = trail(lab, 5),
    pCutoffCol = "p_val_adj",
    h = 8, w = h, ...) {
  # Create an enhanced volcano plot.
  pobj <- EnhancedVolcano::EnhancedVolcano(
    toptable = toptable, ,
    title = title, subtitle = subtitle,
    lab = lab, selectLab = selectLab,
    pCutoffCol = pCutoffCol,
    x = x, y = y,
    ...
  )

  # Save the plot.
  qqSave(ggobj = pobj, title = title, h = h, w = w)
  return(pobj)
}



# _________________________________________________________________________________________________
#' @title BulkGEScatterPlot
#'
#' @description Plots scatterplots of bulk gene expression to identify differentially expressed genes across conditions.
#' @param obj The Seurat object to use for the analysis. Default: combined.obj.
#' @param clusters A string specifying the identity class in the Seurat object to use for cluster assignment. Default: 'cl.names.KnownMarkers.0.2'.
#' @param TwoCategIdent A string specifying the binary categorical identity to split the data for comparison. Default: 'age'.
#' @param genes.from.bulk.DE A character vector specifying the genes obtained from bulk differential expression analysis to be highlighted in the scatterplots. Default: rownames(df.markers.per.AGE).
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   BulkGEScatterPlot(obj = combined.obj, clusters = "cl.names.KnownMarkers.0.2", TwoCategIdent = "age", genes.from.bulk.DE = rownames(df.markers.per.AGE))
#' }
#' }
#' @export
BulkGEScatterPlot <- function(obj = combined.obj # Plot bulk scatterplots to identify differential expressed genes across conditions
                              , clusters = "cl.names.KnownMarkers.0.2", TwoCategIdent = "age", genes.from.bulk.DE = rownames(df.markers.per.AGE)) {
  (SplitIdents <- unique(obj[[TwoCategIdent]][, 1]))
  stopifnot(length(SplitIdents) == 2)

  Idents(obj) <- clusters
  IdentsUsed <- gtools::mixedsort(as.character(unique(Idents(obj))))
  NrPlots <- length(IdentsUsed)
  p.clAv <- p.clAv.AutoLabel <- genes.to.label <- list.fromNames(IdentsUsed)

  # i = 1
  for (i in 1:NrPlots) {
    print(IdentsUsed[i])
    ClX <- subset(obj, idents = IdentsUsed[i])
    Idents(ClX) <- TwoCategIdent
    avg.ClX.cells <- log2(AverageExpression(ClX, verbose = FALSE)$RNA + 1)
    avg.ClX.cells$gene <- rownames(avg.ClX.cells)

    # plot --- --- ---
    p.clAv[[i]] <- p.clAv.AutoLabel[[i]] <-
      ggplot(avg.ClX.cells, aes(x = !!as.name(SplitIdents[1]), y = !!as.name(SplitIdents[2]))) +
      geom_point(data = avg.ClX.cells, color = rgb(0, .5, 0, 0.25), size = 1) +
      FontSize(x.title = 8, x.text = 8, y.title = 8, y.text = 8) +
      geom_abline(slope = 1, intercept = 0, color = "grey") +
      ggtitle(paste("Cluster", IdentsUsed[i])) +
      # ggtitle(paste0("Cluster ", i) ) +
      scale_x_log10() +
      scale_y_log10() +
      annotation_logticks()
    # p.clAv[[i]]

    "Auto identify divergent genes"
    dist.from.axis <- eucl.dist.pairwise(avg.ClX.cells[, 1:2])
    genes.to.label[[i]] <- names(head(sort(dist.from.axis, decreasing = TRUE), n = 20))
    p.clAv.AutoLabel[[i]] <- LabelPoints(plot = p.clAv[[i]], points = genes.to.label[[i]], xnudge = 0, ynudge = 0, repel = TRUE, size = 2)
    p.clAv.AutoLabel[[i]]

    "Pre-identified genes"
    p.clAv[[i]] <- LabelPoints(plot = p.clAv[[i]], points = genes.from.bulk.DE, repel = TRUE, size = 2)
  }

  PlotIter <- CodeAndRoll2::split_vec_to_list_by_N(1:NrPlots, by = 4)
  for (i in 1:length(PlotIter)) {
    plotLS <- p.clAv.AutoLabel[PlotIter[[i]]]
    qqSaveGridA4(plotlist = plotLS, plots = 1:4, fname = ppp("BulkGEScatterPlot.AutoGenes", kpp(PlotIter[[i]]), "png"))

    plotLS <- p.clAv[PlotIter[[i]]]
    qqSaveGridA4(plotlist = plotLS, plots = 1:4, fname = ppp("BulkGEScatterPlot.BulkGenes", kpp(PlotIter[[i]]), "png"))
  }
}





# _________________________________________________________________________________________________
# Compositional analysis ______________________________ ----
# _________________________________________________________________________________________________

#' @title get.clustercomposition
#'
#' @description Get cluster composition: which datasets contribute to each cluster?
#' @param obj Seurat object, Default: combined.obj
#' @param x Bars along the X axis, Default: 'integrated_snn_res.0.3'
#' @param y Vertical split of each bar, Default: 'project'
#' @param color Color, Default: y
#' @param plot  Show plot, Default: TRUE
#' @param ScaleTo100pc Scale the Y Axis, Default: TRUE
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @examples get.clustercomposition()
#' get.clustercomposition()
#' @export
#' @importFrom dplyr group_by_
#' @importFrom scales percent_format
get.clustercomposition <- function(
    obj = combined.obj, ident = "integrated_snn_res.0.3", splitby = "ShortNames",
    color = y,
    plot = TRUE, ScaleTo100pc = TRUE,
    ...) {
  setwd(OutDir)
  clUMAP(obj = obj, ident = x, save.plot = TRUE, suffix = "as.in.barplot")

  (df.meta <- obj@meta.data[, c(ident, splitby)])

  df.meta %>%
    dplyr::group_by_(splitby) %>%
    summarise()

  categ.per.cluster <- ggbarplot(obj@meta.data,
    x = x,
    y = y,
    color = y,
    ...
  )
  if (ScaleTo100pc) categ.per.cluster <- categ.per.cluster + scale_y_discrete(labels = scales::percent_format())
  if (plot) categ.per.cluster

  ggExpress::qqSave(categ.per.cluster, ...)
}


# _________________________________________________________________________________________________
#' @title Generate Barplot of Cell Fractions
#'
#' @description This function generates a bar plot of cell fractions per cluster
#' from a Seurat object. It offers the option to downsample data, which equalizes
#' the number of cells in each group to the number in the smallest group.
#' The plot's bars are grouped by one variable and filled by another.
#'
#' @param obj A Seurat object. The default is 'combined.obj'.
#' @param group.by The variable to group by for the bar plot. The default is 'integrated_snn_res.0.5.ordered'.
#' @param fill.by The variable to fill by for the bar plot. The default is 'age'.
#' @param downsample Logical indicating whether to downsample data. The default is TRUE.
#' @param plotname The title of the plot. The default is 'paste(toTitleCase(fill.by), "proportions")'.
#' @param hlines A numeric vector for the y-intercepts of horizontal lines on the plot. The default is c(0.25, 0.5, 0.75).
#' @param return_table Logical flag indicating whether to return a contingency table instead of a bar plot. Default: FALSE.
#' @param save_plot Logical flag indicating whether to save the plot. Default is TRUE.
#' @param seedNr Seed for random number generation. The default is 1989.
#' @param w The width of the plot. The default is 10.
#' @param h The height of the plot. The default is calculated as 'ceiling(0.5 * w)'.
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @param show_numbers Show numbers on bar
#' @param draw_plot Show plot
#' @param min_frequency Smallest fraction to show individually. Default: 0.025
#' @param custom_col_palette Use custom color palette? Default: c("Standard", "glasbey")[1]
#' @param suffix Plot suffix
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   scBarplot.CellFractions(obj = combined.obj, group.by = "integrated_snn_res.0.1", fill.by = "Phase", downsample = TRUE)
#'   scBarplot.CellFractions(obj = combined.obj, group.by = "integrated_snn_res.0.1", fill.by = "Phase", downsample = FALSE)
#' }
#' }
#' @seealso
#'  \code{\link[tools]{toTitleCase}}
#' @importFrom tools toTitleCase
#' @importFrom dplyr sample_n
#'
#' @export
scBarplot.CellFractions <- function(
    obj = combined.obj,
    group.by = "integrated_snn_res.0.5.ordered", fill.by = "age",
    downsample = TRUE,
    plotname = paste(tools::toTitleCase(fill.by), "proportions"),
    suffix = NULL,
    sub_title = suffix,
    hlines = c(.25, .5, .75),
    return_table = FALSE,
    save_plot = TRUE,
    seedNr = 1989,
    w = 10, h = ceiling(0.5 * w),
    draw_plot = TRUE,
    show_numbers = TRUE,
    min_frequency = 0.025,
    custom_col_palette = c("Standard", "glasbey")[1],
    ...) {
  set.seed(seedNr)
  pname.suffix <- capt.suffix <- NULL
  if (downsample) {
    tbl_X <- table(obj@meta.data[[fill.by]])
    downsample <- min(tbl_X)
    largest_grp <- max(tbl_X)
    pname.suffix <- "(downsampled)"
    capt.suffix <- paste("Downsampled from (max)", largest_grp, "\nto", downsample, "cells in the smallest", fill.by, "group.")
  }
  caption_ <- paste("Numbers denote # cells.", capt.suffix)
  if (min_frequency > 0) caption_ <- paste(caption_, "\nCategories <", percentage_formatter(min_frequency), "are shown together as 'Other'")
  pname_ <- paste(plotname, pname.suffix)

  contingency.table <- table(obj@meta.data[, group.by], obj@meta.data[, fill.by])
  print(contingency.table)

  if (draw_plot) {
    # calculate the proportions and add up small fractions
    prop_table <- obj@meta.data %>%
      group_by(!!as.name(fill.by)) %>%
      summarise(proportion = n() / nrow(obj@meta.data)) %>%
      mutate("category" = ifelse(proportion < min_frequency, "Other", as.character(!!as.name(fill.by))))

    print(unique(prop_table$category))
    palette_x <- color_scale[unique(prop_table$category)]
    print(palette_x)

    # join the proportions back to the original data
    obj@meta.data <- left_join(obj@meta.data, prop_table, by = fill.by)

    subtt <- FixPlotName(group.by, sub_title)
    pl <- obj@meta.data %>%
      {
        if (downsample) dplyr::sample_n(., downsample) else .
      } %>%
      group_by(group_by = !!sym(group.by)) %>%
      ggplot(aes(fill = category, x = !!sym(group.by))) +
      geom_hline(yintercept = hlines, lwd = 1.5) +
      geom_bar(position = "fill") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        title = pname_, subtitle = subtt,
        x = "Clusters", y = "Fraction", caption = caption_
      ) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

    if (custom_col_palette != "Standard") {
      palette_x <- color_scale
      palette_x <- palette_x[unique(prop_table$category)]
      pl <- pl + scale_fill_manual(values = palette_x)
    }

    if (show_numbers) {
      pl <- pl + geom_text(aes(label = ..count..), stat = "count", position = position_fill(vjust = 0.5))
    }

    if (save_plot) {
      sfx <- shorten_clustering_names(group.by)
      if (!is.null(suffix)) sfx <- sppp(sfx, suffix)
      if (min_frequency) sfx <- sppp(sfx, min_frequency)
      qqSave(
        ggobj = pl, title = plotname, also.pdf = TRUE, w = w, h = h,
        suffix = sfx, ...
      )
    } # save_plot
  } # draw_plot

  if (return_table) {
    ls.tables <- list(
      "values" = contingency.table,
      "percentages" = CodeAndRoll2::rowDivide(mat = contingency.table, vec = rowSums(contingency.table))
    )
    return(ls.tables)
  } else {
    return(pl)
  } # else barplot
}





# _________________________________________________________________________________________________
#' @title scBarplot.CellsPerCluster
#'
#' @description Barplot the Fraction of cells per cluster. (dupl?)
#' @param obj Seurat object, Default: combined.obj
#' @param ident identity used, Default: 'cl.names.KnownMarkers.0.5'
#' @param label True: displays cell count, but you can provide anything in a vector.
#' @param palette Color palette. Default: glasbey.
#' @param return_table Should it return the plotting data instead of the plot?
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @param sort Sort by cluster size? Default: FALSE
#' @param suffix File name suffix
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   scBarplot.CellsPerCluster()
#'   scBarplot.CellsPerCluster(sort = TRUE)
#' }
#' }
#' @export scBarplot.CellsPerCluster

scBarplot.CellsPerCluster <- function(
    ident = GetOrderedClusteringRuns()[1],
    sort = FALSE,
    label = list(TRUE, "percent")[[1]],
    suffix = if (label == "percent") "percent" else NULL,
    palette = c("alphabet", "alphabet2", "glasbey", "polychrome", "stepped")[3],
    obj = combined.obj,
    return_table = FALSE,
    ylab_adj = 1.1,
    ...) {
  cell.per.cl <- obj[[ident]][, 1]
  cell.per.cluster <- (table(cell.per.cl))
  if (sort) cell.per.cluster <- sort(cell.per.cluster)
  lbl <- if (isFALSE(label)) {
    NULL
  } else if (label == "percent") {
    percentage_formatter(cell.per.cluster / sum(cell.per.cluster))
  } else if (isTRUE(label)) {
    cell.per.cluster
  } else {
    label
  }

  n.clusters <- length(cell.per.cluster)
  pl <- ggExpress::qbarplot(cell.per.cluster,
    subtitle = ident,
    suffix = kpp(ident, suffix),
    col = 1:n.clusters,
    xlab.angle = 45,
    ylim = c(0, ylab_adj * max(cell.per.cluster)),
    label = lbl,
    ylab = "Cells"
    # , col = getClusterColors(ident = ident, show = TRUE)
    , palette_use = DiscretePalette(n = n.clusters, palette = palette),
    ...
  )

  if (return_table) {
    return(cell.per.cluster)
  } else {
    return(pl)
  }
}



# _________________________________________________________________________________________________
#' @title scBarplot.CellsPerObject
#'
#' @description Creates a bar plot for the number of cells per object from a list of Seurat objects.
#' @param ls.Seu A list of Seurat objects. Default: ls.Seurat.
#' @param plotname A string specifying the title of the plot. Default: 'Nr.Cells.After.Filtering'.
#' @param xlab.angle The angle at which the x-axis labels should be displayed. Default: 45.
#' @param names A logical value indicating whether to use the provided names as labels on the x-axis. If FALSE, the names of the Seurat objects will be used. Default: FALSE.
#' @param ... Additional arguments to be passed to the qbarplot function.
#' @export
scBarplot.CellsPerObject <- function(
    ls.Seu = ls.Seurat,
    plotname = "Nr.Cells.After.Filtering", xlab.angle = 45,
    names = FALSE, ...) {
  cellCounts <- unlapply(ls.Seu, ncol)
  names(cellCounts) <- if (length(names) == length(ls.Seurat)) names else names(ls.Seurat)
  qbarplot(cellCounts,
    plotname = plotname,
    subtitle = paste(sum(cellCounts), "cells in total"),
    label = cellCounts,
    xlab.angle = xlab.angle,
    ylab = "Cells",
    ...
  )
}


# _________________________________________________________________________________________________
#' @title plot.clust.size.distr
#'
#' @description Creates a bar plot or histogram of the cluster size distribution from a given Seurat object.
#' @param obj The Seurat object to be used for the plot. Default: combined.obj.
#' @param ident The identity or clustering to be used for the plot. Default: The second result from GetClusteringRuns().
#' @param plot A logical value indicating whether to plot the data. If FALSE, a vector of the cluster size distribution will be returned. Default: TRUE.
#' @param thr.hist A threshold for the number of clusters above which a histogram will be plotted instead of a bar plot. Default: 30.
#' @param ... Additional arguments to be passed to the internally called plotting functions.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   plot.clust.size.distr()
#' }
#' }
#' @export plot.clust.size.distr
#' @importFrom Stringendo percentage_formatter
plot.clust.size.distr <- function(
    obj = combined.obj, ident = GetClusteringRuns(obj)[2],
    plot = TRUE, thr.hist = 30, ...) {
  clust.size.distr <- table(obj@meta.data[, ident])
  print(clust.size.distr)
  resX <- gsub(pattern = ".*res\\.", replacement = "", x = ident)
  ptitle <- ppp("clust.size.distr", ident)
  psubtitle <- paste(
    "Nr.clusters:", length(clust.size.distr),
    "| median:", median(clust.size.distr),
    "| CV:", Stringendo::percentage_formatter(cv(clust.size.distr))
  )
  xlb <- "Clusters"
  ylb <- "Cluster size (cells)"
  xlim <- c(0, max(clust.size.distr))

  if (plot) {
    if (length(clust.size.distr) < thr.hist) {
      ggExpress::qbarplot(clust.size.distr,
        plotname = ptitle, subtitle = psubtitle,
        label = clust.size.distr, xlab = xlb, ylab = ylb, ...
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

# _________________________________________________________________________________________________
# Correlations _________________________ ----
# _________________________________________________________________________________________________

#' @title sparse.cor
#'
#' @description Calculate a sparse correlation matrix.
#' @param smat A sparse matrix.
#' @return A list with two elements:
#'   * `cov`: The covariance matrix.
#'   * `cor`: The correlation matrix.
#'
#' @export

sparse.cor <- function(smat) {
  n <- nrow(smat)
  cMeans <- colMeans(smat)
  covmat <- (as.matrix(crossprod(smat)) - n * tcrossprod(cMeans)) / (n - 1)
  sdvec <- sqrt(diag(covmat))
  cormat <- covmat / tcrossprod(sdvec)
  list(cov = covmat, cor = cormat)
}

#
# sparse.cor4 <- function(x){
#   n <- nrow(x)
#   cMeans <- colMeans(x)
#   covmat <- (as.matrix(crossprod(x)) - n*tcrossprod(cMeans))/(n-1)
#   sdvec <- sqrt(diag(covmat))
#   cormat <- covmat/tcrossprod(sdvec)
#   list(cov=covmat,cor=cormat)
# }

# _________________________________________________________________________________________________
#' @title Calc.Cor.Seurat
#'
#' @description Calculate gene correlation on a Seurat object.
#' @param assay.use The assay to use from the Seurat object. Default: 'RNA'
#' @param slot.use The slot to use from the assay in the Seurat object. Default: 'data'
#' @param quantileX The quantile level for the calculation. Default: 0.95
#' @param max.cells Maximum number of cells to be used in the calculation. Default: 40000
#' @param seed The random seed used for the calculation. Default: p$seed
#' @param digits The number of decimal places to round the correlation and covariance values. Default: 2
#' @param obj The Seurat object to perform calculations on. Default: combined.obj
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   combined.obj <- calc.q99.Expression.and.set.all.genes(combined.obj, quantileX = 0.99, max.cells = 400000, set.all.genes = FALSE)
#'   combined.obj <- Calc.Cor.Seurat(assay.use = "RNA", slot.use = "data", digits = 2, obj = combined.obj, quantile = 0.99, max.cells = 40000)
#' }
#' }
#' @importFrom tictoc tic toc
#'
#' @export
Calc.Cor.Seurat <- function(
    assay.use = "RNA", slot.use = "data",
    quantileX = 0.95, max.cells = 40000, seed = p$"seed",
    digits = 2, obj = combined.obj) {
  expr.mat <- GetAssayData(slot = slot.use, assay = assay.use, object = obj)
  if (ncol(expr.mat) > max.cells) {
    set.seed(seed = seed)
    cells.use <- sample(x = colnames(expr.mat), size = max.cells)
  } else {
    cells.use <- ncol(expr.mat)
  }

  qname <- paste0("q", quantileX * 100)
  quantile_name <- kpp("expr", qname)

  if (is.null(obj@misc[[quantile_name]])) iprint("Call: combined.obj <- calc.q99.Expression.and.set.all.genes(combined.obj, quantileX =", quantileX, " first )")
  genes.HE <- which_names(obj@misc[[quantile_name]] > 0)
  iprint("Pearson correlation is calculated for", length(genes.HE), "HE genes with expr.", qname, ": > 0.")
  tictoc::tic()
  ls.cor <- sparse.cor(smat = t(expr.mat[genes.HE, cells.use]))
  tictoc::toc()
  ls.cor <- lapply(ls.cor, round, digits = 2)

  slot__name <- kpp(slot.use, assay.use, quantile_name)
  obj@misc[[kpp("cor", slot__name)]] <- ls.cor$"cor"
  obj@misc[[kpp("cov", slot__name)]] <- ls.cor$"cov"
  iprint("Stored under obj@misc$", kpp("cor", slot.use, assay.use), "or cov... .")
  return(obj)
}


# _________________________________________________________________________________________________
#' @title plot.Gene.Cor.Heatmap
#'
#' @description Plot a gene correlation heatmap.
#' @param genes The list of genes of interest. Default: WU.2017.139.IEGsf
#' @param assay.use The assay to use from the Seurat object. Default: 'RNA'
#' @param slot.use The slot to use from the assay in the Seurat object. Default: First item in the vector c("data", "scale.data", "data.imputed")
#' @param quantileX The quantile level for the calculation. Default: 0.95
#' @param min.g.cor The minimum gene correlation. Genes with correlation above this value or below its negative will be included. Default: 0.3
#' @param calc.COR A logical flag, if TRUE, the correlation will be calculated if not found. Default: FALSE
#' @param cutRows A number, the dendrogram will be cut at this height into clusters. Default: NULL
#' @param cutCols A number, the dendrogram will be cut at this height into clusters. If NULL, cutCols will be the same as cutRows. Default: cutRows
#' @param obj The Seurat object to perform calculations on. Default: combined.obj
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @importFrom pheatmap pheatmap
#' @importFrom MarkdownReports wplot_save_pheatmap
#'
#' @export plot.Gene.Cor.Heatmap
plot.Gene.Cor.Heatmap <- function(
    genes = WU.2017.139.IEGsf,
    assay.use = "RNA", slot.use = c("data", "scale.data", "data.imputed")[1], quantileX = 0.95,
    min.g.cor = 0.3, calc.COR = FALSE,
    cutRows = NULL, cutCols = cutRows,
    obj = combined.obj, ...) {
  expr.mat <- GetAssayData(slot = slot.use, assay = assay.use, object = obj)
  if (slot.use == c("data.imputed")) {
    "WIP"
  }
  expr.mat <- GetAssayData(slot = slot.use, assay = assay.use, object = obj)

  qname <- paste0("expr.q", quantileX * 100)
  slotname_cor.mat <- kpp("cor", slot.use, assay.use, qname)
  cor.mat <- obj@misc[[slotname_cor.mat]]

  if (is.null(cor.mat)) {
    iprint(slotname_cor.mat, " not found in @misc.")
    iprint("Correlation slots present in @misc:", CodeAndRoll2::grepv(names(obj@misc), pattern = "^cor"))

    # Calculate --- --- --- --- ---
    if (calc.COR) {
      print("Calculating correlation now.")
      genes.found <- check.genes(list.of.genes = genes)
      iprint(length(genes.found), "genes are found in the object.")
      if (length(genes.found) > 200) iprint("Too many genes found in data, cor will be slow: ", length(genes.found))
      ls.cor <- sparse.cor(t(expr.mat[genes.found, ]))
      cor.mat <- ls.cor$cor
    } else {
      stop()
    }
  } else {
    print("Correlation is pre-calculated")
    genes.found <- intersect(genes, rownames(cor.mat))
    iprint(length(genes.found), "genes are found in the correlation matrix.")
    cor.mat <- cor.mat[genes.found, genes.found]
  }


  # Filter --- --- --- --- --- ---
  diag(cor.mat) <- NaN
  corgene.names <- union(
    which_names(rowMax(cor.mat) >= min.g.cor),
    which_names(rowMin(cor.mat) <= -min.g.cor)
  )
  iprint(length(corgene.names), "genes are more (anti-)correlated than +/-:", min.g.cor)

  pname <- paste0("Pearson correlations of ", substitute(genes), "\n min.cor:", min.g.cor, " | ", assay.use, ".", slot.use)
  o.heatmap <- pheatmap::pheatmap(cor.mat[corgene.names, corgene.names], main = pname, cutree_rows = cutRows, cutree_cols = cutCols, ...)
  MarkdownReports::wplot_save_pheatmap( o.heatmap, plotname = make.names(pname))

  # return values
  maxCorrz <- rowMax(cor.mat)[corgene.names]
  names(maxCorrz) <- corgene.names
  dput(maxCorrz)
}





# _________________________________________________________________________________________________
# Seurat.object.manipulations.etc.R ______________________________ ----
# _________________________________________________________________________________________________
# source('~/GitHub/Packages/Seurat.utils/Functions/Seurat.object.manipulations.etc.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Seurat.object.manipulations.etc.R"))



# _________________________________________________________________________________________________
#' @title prefix_cells_seurat
#'
#' @description This function adds prefixes from 'obj_IDs' to cell names in Seurat S4 objects from 'ls_obj'
#'
#' @param ls_obj List. A list of Seurat S4 objects
#' @param obj_IDs Character vector. A vector of sample IDs corresponding to each Seurat S4 object
#'
#' @examples
#' # ls_obj <- list(seurat_obj1, seurat_obj2)
#' # obj_IDs <- c("sample1", "sample2")
#' # ls_obj_prefixed <- prefix_cells_seurat(ls_obj = ls_obj, obj_IDs = obj_IDs)
#' @export
prefix_cells_seurat <- function(ls_obj, obj_IDs) {
  # Check if 'ls_obj' is a list of Seurat objects and 'obj_IDs' is a character vector of the same length
  if (!is.list(ls_obj) & inherits(ls_obj, "Seurat")) ls_obj <- list(ls_obj)
  stopifnot(is.list(ls_obj) & all(sapply(ls_obj, function(x) inherits(x, "Seurat"))))
  stopifnot(is.character(obj_IDs) & length(ls_obj) == length(obj_IDs))

  names_orig <- names(ls_obj)

  # Iterate over Seurat objects
  ls_obj_prefixed <- lapply(seq_along(ls_obj), function(i) {
    # Get the Seurat object and corresponding prefix
    obj <- ls_obj[[i]]
    prefix <- obj_IDs[i]

    # Add prefix to cell names
    new_cell_names <- paste0(prefix, "_", colnames(obj))

    # Rename cells in the Seurat object
    obj <- RenameCells(obj, new.names = new_cell_names)

    return(obj)
  })
  print(lapply(lapply(ls_obj_prefixed, colnames), head))

  names(ls_obj_prefixed) <- names_orig
  return(ls_obj_prefixed)
}

# _________________________________________________________________________________________________
#' @title Check Prefix in Seurat Object Cell IDs
#'
#' @description This function checks if a prefix has been added to the standard cell-IDs (16 characters of A,T,C,G)
#' in a Seurat object. If so, it prints the number of unique prefixes found,
#' issues a warning if more than one unique prefix is found, and returns the identified prefix(es).
#' @param obj A Seurat object with cell IDs possibly prefixed.
#' @param cell_ID_pattern Pattern to match cellIDs (with any suffix).
#' @return A character vector of the identified prefix(es).
#'
#' @examples
#' # Assuming 'obj' is your Seurat object
#' # prefix <- find_prefix_in_cell_IDs(obj)
#'
#' @export

find_prefix_in_cell_IDs <- function(obj, cell_ID_pattern = "[ATCG]{16}.*$") {
  stopifnot(inherits(obj, "Seurat"))

  # Extract cell IDs
  cell_IDs <- colnames(obj)

  # Remove the standard 16-character cell-IDs
  potential_prefixes <- gsub(pattern = cell_ID_pattern, replacement = "", x = cell_IDs)

  # Check if there is no prefix
  if (all(potential_prefixes == "")) {
    print("No prefix found in cell IDs.")
    return(NULL)
  }

  # Identify unique prefixes
  unique_prefixes <- unique(potential_prefixes)

  # Print the number of unique prefixes
  print(paste(length(unique_prefixes), "unique prefix(es) found:", head(unique_prefixes)))

  # Issue a warning if more than one unique prefix is found
  if (length(unique_prefixes) > 1) {
    warning("Multiple unique prefixes identified in cell IDs:", head(unique_prefixes), immediate. = TRUE)
  }

  # Return the identified prefix(es)
  return(unique_prefixes)
}




# _________________________________________________________________________________________________
#' @title seu.Make.Cl.Label.per.cell
#'
#' @description Take a named vector (of e.g. values ="gene names", names = clusterID), and a vector of cell-IDs and make a vector of "GeneName.ClusterID".
#' @param TopGenes A named vector, where values are gene names and names are cluster IDs.
#' @param clID.per.cell A vector of cell-IDs used to create the output vector.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   seu.Make.Cl.Label.per.cell(TopGenes = TopGenes.Classic, clID.per.cell = getMetadataColumn(ColName.metadata = metaD.CL.colname))
#' }
#' }
#' @export

seu.Make.Cl.Label.per.cell <- function(TopGenes, clID.per.cell) { # Take a named vector (of e.g. values ="gene names", names = clusterID), and a vector of cell-IDs and make a vector of "GeneName.ClusterID".
  Cl.names_class <- TopGenes[clID.per.cell]
  Cl.names_wNr <- paste0(Cl.names_class, " (", names(Cl.names_class), ")")
  return(Cl.names_wNr)
}


# _________________________________________________________________________________________________
#' @title GetMostVarGenes
#'
#' @description Get the N most variable Genes
#' @param obj A Seurat object.
#' @param nGenes Number of genes, Default: p$nVarGenes
#' @export
GetMostVarGenes <- function(obj, nGenes = p$nVarGenes) { # Get the most variable rGenes
  head(rownames(slot(object = obj, name = "hvg.info")), n = nGenes)
}

# _________________________________________________________________________________________________
#' @title gene.name.check
#'
#' @description Check gene names in a seurat object, for naming conventions (e.g.: mitochondrial reads have - or .). Use for reading .mtx & writing .rds files. #
#' @param Seu.obj A Seurat object.
#' @importFrom MarkdownHelpers llprint llogit
#'
#' @export
gene.name.check <- function(Seu.obj) { # Check gene names in a seurat object, for naming conventions (e.g.: mitochondrial reads have - or .). Use for reading .mtx & writing .rds files.
  rn <- rownames(GetAssayData(object = Seu.obj, slot = "counts"))
  MarkdownHelpers::llprint("### Gene name pattern")

  MarkdownHelpers::llogit('`rn = rownames(GetAssayData(object = ls.Seurat[[1]], slot = "counts"))`')
  MarkdownHelpers::llogit('`head(CodeAndRoll2::grepv(rn, pattern = "-"), 10)`')
  print("pattern = -")
  MarkdownHelpers::llprint(head(CodeAndRoll2::grepv(rn, pattern = "-"), 10))

  MarkdownHelpers::llogit('`head(CodeAndRoll2::grepv(rn, pattern = "_"), 10)`')
  print("pattern = _")
  MarkdownHelpers::llprint(head(CodeAndRoll2::grepv(rn, pattern = "_"), 10))

  MarkdownHelpers::llogit('`head(CodeAndRoll2::grepv(rn, pattern = "\\."), 10)`')
  print("pattern = \\.")
  MarkdownHelpers::llprint(head(CodeAndRoll2::grepv(rn, pattern = "\\."), 10))

  MarkdownHelpers::llogit('`head(CodeAndRoll2::grepv(rn, pattern = "\\.AS[1-9]"), 10)`')
  print("pattern = \\.AS[1-9]")
  MarkdownHelpers::llprint(head(CodeAndRoll2::grepv(rn, pattern = "\\.AS[1-9]"), 10))
}


# _________________________________________________________________________________________________
#' @title check.genes
#'
#' @description Check if a gene name exists in a Seurat object, or in HGNC?
#' @param list.of.genes A vector of gene names to check for existence in the Seurat object or HGNC. Default: ClassicMarkers
#' @param makeuppercase Logical, if TRUE, transforms all gene names to uppercase. Default: FALSE
#' @param verbose Logical, if TRUE, prints out information about missing genes. Default: TRUE
#' @param HGNC.lookup Logical, if TRUE, looks up missing genes in HGNC database. Default: FALSE
#' @param obj Seurat object that the function uses to check the gene names. Default: combined.obj
#' @param assay.slot The assay slot of the Seurat object to consider. Default: 'RNA'
#' @param dataslot The data slot of the Seurat object to consider. Default: 'data'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   check.genes("top2a", makeuppercase = TRUE)
#'   check.genes("VGLUT2", verbose = FALSE, HGNC.lookup = TRUE)
#' }
#' }
#' @importFrom DatabaseLinke.R qHGNC
#' @importFrom Stringendo percentage_formatter
#'
#' @export
check.genes <- function(
    list.of.genes = ClassicMarkers, makeuppercase = FALSE, verbose = TRUE, HGNC.lookup = FALSE,
    obj = combined.obj,
    assay.slot = c("RNA", "integrated")[1],
    dataslot = c("counts", "data")[2]) { # Check if genes exist in your dataset.
  if (makeuppercase) list.of.genes <- toupper(list.of.genes)
  all_genes <- rownames(GetAssayData(object = obj, assay = assay.slot, slot = dataslot))
  length(all_genes)
  missingGenes <- setdiff(list.of.genes, all_genes)
  if (length(missingGenes) > 0) {
    if (verbose) {
      iprint(length(missingGenes), "or", Stringendo::percentage_formatter(length(missingGenes) / length(list.of.genes)), "genes not found in the data, e.g:", head(missingGenes, n = 10))
    }
    if (HGNC.lookup) {
      if (exists("qHGNC", mode = "function")) {
        try(DatabaseLinke.R::qHGNC(missingGenes))
      } else {
        print("load qHGNC() function, see database.linker")
      }
    }
  }
  intersect(list.of.genes, all_genes)
}



# _________________________________________________________________________________________________
#' @title fixZeroIndexing.seurat
#'
#' @description Fix zero indexing in seurat clustering, to 1-based indexing. replace zero indexed clusternames.
#' @param ColName.metadata Metadata column name to use, Default: 'res.0.6'
#' @param obj Seurat object, Default: org
#' @export
fixZeroIndexing.seurat <- function(ColName.metadata = "res.0.6", obj = org) { # Fix zero indexing in seurat clustering, to 1-based indexing
  obj@meta.data[, ColName.metadata] <- as.numeric(obj@meta.data[, ColName.metadata]) + 1
  print(obj@meta.data[, ColName.metadata])
  return(obj)
}


# _________________________________________________________________________________________________
#' @title CalculateFractionInTranscriptome
#'
#' @description This function calculates the fraction of a set of genes within the full transcriptome of each cell.
#' @param geneset A character vector specifying the set of genes for which the fraction in the transcriptome is to be calculated. Default is c("MALAT1").
#' @param obj A Seurat object from which the gene data is extracted. Default is 'combined.obj'.
#' @param dataslot A character vector specifying the data slot to be used in the calculation. Default is 'data' (second element in the vector c("counts", "data")).
#' @return A numeric vector containing the fraction of the specified genes in the transcriptome of each cell.
#' @export
CalculateFractionInTrome <- function(
    genesCalc.Cor.Seuratet = c("MALAT1") # Calculate the fraction of a set of genes within the full Transcriptome of each cell.
    , obj = combined.obj,
    dataslot = c("counts", "data")[2]) {
  warning("    >>>> Use addMetaFraction() <<<<", immediate. = TRUE)
  geneset <- check.genes(list.of.genes = geneset)
  stopifnot(length(geneset) > 0)

  mat <- as.matrix(slot(obj@assays$RNA, name = dataslot))
  mat.sub <- mat[geneset, , drop = FALSE]
  RC.per.cell.geneset <- colSums(mat.sub)

  RC.per.cell <- colSums(mat)
  gene.fraction.per.cell <- 100 * RC.per.cell.geneset / RC.per.cell
  return(gene.fraction.per.cell)
}

# _________________________________________________________________________________________________
#' @title AddNewAnnotation
#'
#' @description This function creates a new metadata column based on an existing metadata column and a list of mappings (name <- IDs).
#' @param obj A Seurat object for which the new annotation is to be created. Default is 'obj'.
#' @param source A character string specifying the existing metadata column to be used as the basis for the new annotation. Default is 'RNA_snn_res.0.5'.
#' @param named.list.of.identities A named list providing the mappings for the new annotation. Default is 'ls.Subset.ClusterLists'.
#' @return A character vector representing the new metadata column.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   ls.Subset.ClusterLists <- list("hESC.h9" = c("4", "10", "14"), "hESC.176" = c("0", "1", "2"))
#'   AddNewAnnotation()
#' }
#' }
#' @export
AddNewAnnotation <- function(obj = obj # Create a new metadata column based on an exisiting metadata column and a list of mappings (name <- IDs).
                             , source = "RNA_snn_res.0.5", named.list.of.identities = ls.Subset.ClusterLists) {
  NewID <- as.named.vector.df(obj[[source]])

  for (i in 1:length(named.list.of.identities)) {
    lx <- as.character(named.list.of.identities[[i]])
    name.lx <- names(named.list.of.identities)[i]
    NewID <- translate(vec = NewID, oldvalues = lx, newvalues = name.lx)
  }
  print(table(NewID))
  return(NewID)
}


# _________________________________________________________________________________________________
#' @title whitelist.subset.ls.Seurat
#'
#' @description Subsets cells in a list of Seurat objects based on an externally provided list of cell IDs.
#' @param ls.obj A list of Seurat objects. Default: ls.Seurat.
#' @param metadir Directory for the metadata. Default: p$cellWhiteList.
#' @param whitelist.file Filename of the whitelist containing cell IDs. Default: "NonStressedCellIDs.2020.10.21_18h.tsv".
#' @return A list of Seurat objects containing only the cells specified in the whitelist.
#' @details The function first validates the presence of all identities from the metadata in the Seurat objects. If all identities are present, the function subsets each Seurat object based on the whitelist of cell IDs.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   ls.Seurat.subset <- whitelist.subset.ls.Seurat(ls.obj = ls.Seurat, metadir = p$"cellWhiteList", whitelist.file = "NonStressedCellIDs.2020.10.21_18h.tsv")
#' }
#' }
#' @seealso
#' \code{\link[Seurat]{subset}}
#' @importFrom ReadWriter read.simple.tsv
#'
#' @export
whitelist.subset.ls.Seurat <- function(
    ls.obj = ls.Seurat,
    metadir = p$"cellWhiteList" #  '~/Dropbox/Abel.IMBA/MetadataD/POL.meta/cell.lists/'
    , whitelist.file = "NonStressedCellIDs.2020.10.21_18h.tsv") {
  cells.before <- unlapply(ls.obj, ncol)
  # Find file
  df.cell.whitelist <- ReadWriter::read.simple.tsv(metadir, whitelist.file)
  dsets <- table(df.cell.whitelist[, 1])

  ls.orig.idents <- lapply(lapply(ls.Seurat, getMetadataColumn, ColName.metadata = "orig.ident"), unique)
  stopif(any(unlapply(ls.orig.idents, l) == length(ls.Seurat)), message = "Some ls.Seurat objects have 1+ orig identity.")

  dsets.in.lsSeu <- unlist(ls.orig.idents)
  isMathced <- all(dsets.in.lsSeu == names(dsets)) # Stop if either ls.Seurat OR the metadata has identities not found in the other, in the same order.
  stopif(!isMathced, message = paste(
    "either ls.Seurat OR the metadata has identities not found in the other, or they are not in same order.",
    kpps(dsets.in.lsSeu), "vs.", kpps(names(dsets))
  ))

  # identX <- ls.orig.idents[[1]]
  for (i in 1:length(ls.orig.idents)) {
    identX <- ls.orig.idents[[i]]
    print(identX)

    # Extract and process cellIDs
    idx.match <- which(df.cell.whitelist[, 1] == identX)
    cell.whitelist <- rownames(df.cell.whitelist)[idx.match]
    cell.whitelist <- substr(
      x = cell.whitelist,
      start = 1, stop = nchar(cell.whitelist) - 2
    )

    # Extract and process cellIDs
    ls.obj[[i]] <- subset(x = ls.obj[[i]], cells = cell.whitelist)
  }
  cells.after <- unlapply(ls.obj, ncol)
  iprint("cells.before", cells.before, "cells.after", cells.after)
  return(ls.obj)
}

# _________________________________________________________________________________________________
#' @title FindCorrelatedGenes
#'
#' @description Find correlated genes in a Seurat object
#' @param gene Gene of interest. Default: 'TOP2A'
#' @param obj Seurat object to find the correlated genes from. Default: combined.obj
#' @param assay Assay to be used from the Seurat object. Default: 'RNA'
#' @param slot Slot to be used from the specified assay in the Seurat object. Default: 'data'
#' @param HEonly Logical, if TRUE, filters matrix to high-expressing genes only. Default: FALSE
#' @param minExpr Minimum expression level for a gene to be considered. Default: 1
#' @param minCells Minimum number of cells expressing a gene for the gene to be considered. Default: 1000
#' @param trailingNgenes Number of top genes to consider based on their correlation. Default: 1000
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   FindCorrelatedGenes(gene = "TOP2A", obj = combined.obj)
#'   write_clip(names(head(topGenes[-(1:6)], n = 50)))
#' }
#' }
#' @seealso
#'  \code{\link[matrixStats]{rowSums2}}
#' @importFrom matrixStats rowSums2
#' @importFrom tictoc tic toc
#' @importFrom MarkdownReports wbarplot
#'
#' @export
FindCorrelatedGenes <- function(
    gene = "TOP2A", obj = combined.obj, assay = "RNA", slot = "data",
    HEonly = FALSE, minExpr = 1, minCells = 1000,
    trailingNgenes = 1000) {
  tictoc::tic()
  AssayData <- GetAssayData(object = obj, assay = assay, slot = slot)
  matrix_mod <- iround(as.matrix(AssayData))
  if (HEonly) {
    idx.pass <- (matrixStats::rowSums2(matrix_mod > minExpr) > minCells)
    pc_TRUE(idx.pass)
    matrix_mod <- matrix_mod[which(idx.pass), ]
  }
  geneExpr <- as.numeric(matrix_mod[gene, ])
  correlations <- apply(matrix_mod, 1, cor, geneExpr)
  topGenes <- trail(sort(correlations, decreasing = TRUE), N = trailingNgenes)
  tictoc::toc()
  MarkdownReports::wbarplot(head(topGenes, n = 25))
  topGenes
}



# _________________________________________________________________________________________________
# _________________________________________________________________________________________________



# _________________________________________________________________________________________________
# Seurat.update.gene.symbols.HGNC.R ______________________________ ----
# _________________________________________________________________________________________________
# source('~/GitHub/Packages/Seurat.utils/Functions/Seurat.update.gene.symbols.HGNC.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Seurat.update.gene.symbols.HGNC.R"))
# require(HGNChelper)

# _________________________________________________________________________________________________
#' @title UpdateGenesSeurat
#'
#' @description Update genes symbols that are stored in a Seurat object. It returns a data frame. The last column are the updated gene names.
#' @param obj Seurat object to update gene symbols in. Default: ls.Seurat[[i]]
#' @param species_ Species to which the gene symbols correspond, used to check gene symbols. Default: 'human'
#' @param EnforceUnique Logical, if TRUE, enforces unique gene symbols. Default: TRUE
#' @param ShowStats Logical, if TRUE, displays statistics of gene symbols update. Default: FALSE
#' @seealso
#'  \code{\link[HGNChelper]{checkGeneSymbols}}
#'
#' @export
#' @importFrom HGNChelper checkGeneSymbols

UpdateGenesSeurat <- function(obj = ls.Seurat[[i]], species_ = "human", EnforceUnique = TRUE, ShowStats = FALSE) { # Update genes symbols that are stored in a Seurat object. It returns a data frame. The last column are the updated gene names.
  HGNC.updated <- HGNChelper::checkGeneSymbols(rownames(obj), unmapped.as.na = FALSE, map = NULL, species = species_)
  if (EnforceUnique) HGNC.updated <- HGNC.EnforceUnique(HGNC.updated)
  if (ShowStats) print(GetUpdateStats(HGNC.updated))
  obj <- RenameGenesSeurat(obj, newnames = HGNC.updated$"Suggested.Symbol")
  return(obj)
}


# _________________________________________________________________________________________________
#' @title RenameGenesSeurat
#'
#' @description Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data. #
#' @param obj Seurat object, Default: ls.Seurat[[i]]
#' @param assay Which Seurat assay to replace. Default: RNA. Disclaimer: Intended use on simple objects that ONLY contain an RNA object. I highly advise against selectively replacing name in other assays that may have slots that cannot be updated by this function.
#' @param newnames A vector of new gene names. Default: HGNC.updated[[i]]$Suggested.Symbol
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   RenameGenesSeurat(obj = SeuratObj, newnames = HGNC.updated.genes$Suggested.Symbol)
#' }
#' }
#' @export
RenameGenesSeurat <- function(obj = ls.Seurat[[i]],
                              newnames = HGNC.updated[[i]]$Suggested.Symbol,
                              assay = "RNA",
                              slots = c("data", "counts", "scale.data", "meta.features")
                              ) {
  warning("Run this before integration and downstream processing. It only attempts to change
          @counts, @data, @scale.data and @meta.features in obj@assays$YOUR_ASSAY.", immediate. = TRUE)

  if (nrow(obj) == length(newnames)) {
    print(paste("Present:", SeuratObject::Layers(obj@assays[[assay]])))
    for (s in slots) {
      # browser()
      nrO <- nrow(SeuratObject::GetAssayData(object = obj, assay = assay, layer = s))
      obj <- .check_and_rename(obj, assay, newnames = newnames, layer.name = s)
      nrN <- nrow(SeuratObject::GetAssayData(object = obj, assay = assay, layer = s))
      stopifnot(nrN == nrO)
    }
  } else {
    warning("Unequal gene sets: nrow(assayobj) != nrow(newnames). No renaming performed!", immediate. = TRUE)
  }
  return(obj)
}


# _________________________________________________________________________________________________
#' @title Check and Rename Gene Names in Seurat Assay Object
#'
#' @description This function renames rows (genes) in a specified slot of a Seurat assay object.
#' It supports slots storing data as either a dense or a sparse matrix (dgCMatrix) or data.frame.
#'
#' @param obj A Seurat object.
#' @param assay An Assay name in a Seurat object.
#' @param newnames A character vector of new gene names to be assigned.
#' @param layer.name A string specifying the slot in the Assay object to be updated.
#'                 Valid options typically include 'counts', 'data', or 'scale.data'.
#'
#' @return An Assay object with updated gene names in the specified slot.
#' @examples
#' \dontrun{
#'   # Assuming 'seurat_obj' is a Seurat object and 'new_gene_names' is a vector of gene names
#'   updated_assay <- check_and_rename(assayobj = seurat_obj[["RNA"]],
#'                                     newnames = new_gene_names,
#'                                     layer.name = "counts")
#' }

.check_and_rename <- function(obj, assay, newnames, layer.name) {
  cat(layer.name, fill = TRUE)

  stopifnot(
    is(obj, "Seurat"),
    is.character(assay),
    is.character(layer.name),
    is.character(newnames),
    nrow(obj) == length(newnames)
    )

  assayobj <- obj@assays[[assay]]
  feature.list <- rownames(assayobj@features@.Data)

  if (length(feature.list) == length(newnames)) {
    rownames(assayobj@features@.Data) <- newnames
    nrX <- length(rownames(assayobj@features@.Data))
  } else {
    iprint("length feature.list", length(feature.list), "length newnames", length(newnames))
    stop()
  }

  if(layer.name %in% SeuratObject::Layers(assayobj)) {

    matrix_n <- SeuratObject::LayerData(assayobj, layer = layer.name)
    nr1 <- nrow(matrix_n)

    if (all(dim(matrix_n)) > 0) {
      # browser()
      stopifnot(nrow(matrix_n) == length(newnames))

      if ("dgCMatrix" %in% class(matrix_n)) {
        message(assay, "@", layer.name, " is of type dgeCMatrix!")
        matrix_n@Dimnames[[1]] <- newnames

      } else if ("matrix" %in% class(matrix_n)) {
        message(assay, "@", layer.name, " is of type Matrix!")
        rownames(matrix_n) <- newnames

      } else if ("data.frame" %in% class(matrix_n)) {
        message(assay, "@", layer.name, " is of type data.frame!")
        rownames(matrix_n) <- newnames

      } else {
        warning(">>> No renaming: ", assay, "@", layer.name,
                " not of type dgeCMatrix / Matrix / data.frame.", immediate. = TRUE)
      }
      stopifnot(nr1 == nrow(matrix_n))

      SeuratObject::LayerData(assayobj, layer = layer.name) <- matrix_n
      nr3 <- nrow(SeuratObject::LayerData(assayobj, layer = layer.name))
      stopifnot(nr3 == nrX)
    }

  } else {
    warning(paste(">>>", assay, "@", layer.name, "does not exist!"), immediate. = TRUE)

  }
  # obj <- SetAssayData(obj, layer = layer.name, new.data = matrix_n)
  obj@assays[[assay]] <- assayobj
  return(obj)
}

# _________________________________________________________________________________________________
#' @title RemoveGenesSeurat
#'
#' @description Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes metadata; obj@assays$RNA@counts, @data and @scale.data. #
#' @param obj Seurat object, Default: ls.Seurat[[i]]
#' @param symbols2remove Genes to remove from a Seurat object. Default: c("TOP2A")
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   RemoveGenesSeurat(obj = SeuratObj, symbols2remove = "TOP2A")
#' }
#' }
#' @export
RemoveGenesSeurat <- function(obj = ls.Seurat[[i]], symbols2remove = c("TOP2A")) { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes metadata; obj@assays$RNA@counts, @data and @scale.data.
  print("Run this as the first thing after creating the Seurat object. It only removes genes from: metadata; obj@assays$RNA@counts, @data and @scale.data.")
  RNA <- obj@assays$RNA

  if (length(RNA@counts)) {
    NotFound <- setdiff(symbols2remove, RNA@counts@Dimnames[[1]])
    if (length(NotFound) == 0) {
      RNA@counts@Dimnames[[1]] <- symbols2remove
      print("Genes removed from RNA@counts")
    } else {
      print("Not All Genes Found in RNA@counts. Missing:")
      print(NotFound)
    }
  }
  if (length(RNA@data)) {
    if (length(setdiff(symbols2remove, RNA@data@Dimnames[[1]])) == 0) {
      RNA@data@Dimnames[[1]] <- symbols2remove
      print("Genes removed from RNA@data.")
    } else {
      print("Not All Genes Found in RNA@data")
    }
  }
  if (length(RNA@scale.data)) {
    if (length(setdiff(symbols2remove, RNA@scale.data@Dimnames[[1]])) == 0) {
      RNA@scale.data@Dimnames[[1]] <- symbols2remove
      print("Genes removed from RNA@scale.data.")
    } else {
      print("Not All Genes Found in RNA@scale.data")
    }
  }
  if (length(obj@meta.data)) {
    if (length(setdiff(symbols2remove, rownames(obj@meta.data))) == 0) {
      rownames(obj@meta.data) <- symbols2remove
      print("Genes removed from @meta.data.")
    } else {
      print("Not All Genes Found in @metadata")
    }
  }
  obj@assays$RNA <- RNA
  return(obj)
}



# _________________________________________________________________________________________________
#' @title HGNC.EnforceUnique
#'
#' @description Enforce Unique names after HGNC symbol update. While "make.unique" is not the ideal
#' solution, because it generates mismatches, in my integration example it does reduce the
#' mismatching genes from ~800 to 4
#' @param updatedSymbols Gene symbols, it is the output of HGNChelper's checkGeneSymbols)_.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   x <- HGNC.EnforceUnique(updatedSymbols = SymUpd)
#' }
#' }
#' @export
HGNC.EnforceUnique <- function(updatedSymbols) {
  NGL <- updatedSymbols[, 3]
  if (any.duplicated(NGL)) {
    updatedSymbols[, 3] <- make.unique(NGL)
    "Unique names are enforced by suffixing .1, .2, etc."
  }
  return(updatedSymbols)
}




# _________________________________________________________________________________________________
#' @title GetUpdateStats
#'
#' @description Plot the Symbol-update statistics. Works on the data frame returned by `UpdateGenesSeurat()`. #
#' @param genes Genes of iinterest, Default: HGNC.updated[[i]]
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   GetUpdateStats(genes = HGNC.updated.genes)
#' }
#' }
#' @importFrom Stringendo percentage_formatter
#' @export
GetUpdateStats <- function(genes = HGNC.updated[[i]]) { # Plot the Symbol-update statistics. Works on the data frame returned by `UpdateGenesSeurat()`.
  MarkedAsUpdated <- genes[genes$Approved == FALSE, ]
  AcutallyUpdated <- sum(MarkedAsUpdated[, 1] != MarkedAsUpdated[, 3])
  UpdateStats <- c("Updated (%)" = Stringendo::percentage_formatter(AcutallyUpdated / nrow(genes)), "Updated Genes" = floor(AcutallyUpdated), "Total Genes" = floor(nrow(genes)))
  return(UpdateStats)
}


# _________________________________________________________________________________________________
#' @title PlotUpdateStats
#'
#' @description Creates a scatter plot of update statistics.
#' @param mat A matrix containing update statistics. Default: UpdateStatMat.
#' @param column.names A character vector of column names in the mat parameter. Default: c("Updated (%)", "Updated (Nr.)").
#' @return A scatter plot displaying update statistics.
#' @details This function takes a matrix containing update statistics and column names to plot the corresponding statistics. It colorizes the genes and plots the percentage of total genes updated against the number of genes updated.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   PlotUpdateStats(mat = result.of.GetUpdateStats)
#' }
#' }
#' @seealso
#' \code{\link[wplot]{wplot}}, \code{\link[wcolorize]{wcolorize}}
#' @importFrom MarkdownReports wplot wlegend
#'
#' @export
PlotUpdateStats <- function(mat = UpdateStatMat, column.names = c("Updated (%)", "Updated (Nr.)")) { # Scatter plot of update stats.
  stopifnot(column.names %in% colnames(UpdateStatMat))
  HGNC.UpdateStatistics <- mat[, column.names]
  HGNC.UpdateStatistics[, "Updated (%)"] <- 100 * HGNC.UpdateStatistics[, "Updated (%)"]
  colnames(HGNC.UpdateStatistics) <- c("Gene Symbols updated (% of Total Genes)", "Number of Gene Symbols updated")
  lll <- wcolorize(vector = rownames(HGNC.UpdateStatistics))
  MarkdownReports::wplot(HGNC.UpdateStatistics,
    col = lll,
    xlim = c(0, max(HGNC.UpdateStatistics[, 1])),
    ylim = c(0, max(HGNC.UpdateStatistics[, 2]))
  )
  MarkdownReports::wlegend(NamedColorVec = lll, poz = 1)
}



# _________________________________________________________________________________________________
# Handling SNP demux table results coming from SoupOrCell ______________________________ ----
# _________________________________________________________________________________________________






# _________________________________________________________________________________________________


# _________________________________________________________________________________________________
# Read.Write.Save.Load.functions.R ______________________________ ----
# _________________________________________________________________________________________________
# source('~/GitHub/Packages/Seurat.utils/Functions/Read.Write.Save.Load.functions.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Read.Write.Save.Load.functions.R"))

"Multicore read / write (I/O) functions are https://github.com/vertesy/Seurat.multicore"
"Single core read / write (I/O) functions are in https://github.com/vertesy/Seurat.utils/"


# _________________________________________________________________________________________________
#' @title Convert10Xfolders
#'
#' @description This function takes a parent directory with a number of subfolders, each containing the standard output of 10X Cell Ranger. It (1) loads the filtered data matrices, (2) converts them to Seurat objects, and (3) saves them as .RDS files.
#' @param InputDir A character string specifying the input directory.
#' @param regex A logical value. If TRUE, the folderPattern is treated as a regular expression. Default is FALSE.
#' @param folderPattern A character vector specifying the pattern of folder names to be searched. Default is 'filtered_feature'.
#' @param min.cells An integer value specifying the minimum number of cells. Default is 5.
#' @param min.features An integer value specifying the minimum number of features. Default is 200.
#' @param updateHGNC A logical value indicating whether to update the HGNC. Default is TRUE.
#' @param ShowStats A logical value indicating whether to show statistics. Default is TRUE.
#' @param writeCBCtable A logical value indicating whether to write out a list of cell barcodes (CBC) as a tsv file. Default is TRUE.
#' @param depth An integer value specifying the depth of scan (i.e., how many levels below the InputDir). Default is 2.
#' @param sample.barcoding A logical value indicating whether Cell Ranger was run with sample barcoding. Default is FALSE.
#' @param sort_alphanumeric sort files alphanumeric? Default: TRUE.
#' @examples
#' \dontrun{
#' if (interactive()) Convert10Xfolders(InputDir)
#' }
#' @export
Convert10Xfolders <- function(
    InputDir,
    regex = FALSE,
    folderPattern = c("filtered_feature", "SoupX_decont")[1],
    depth = 4,
    min.cells = 5, min.features = 200,
    updateHGNC = TRUE, ShowStats = TRUE,
    writeCBCtable = TRUE,
    sample.barcoding = FALSE,
    nthreads = 12,
    preset = "high",
    ext = "qs",
    sort_alphanumeric = TRUE,
    ...) {

  warning("Since v2.5.0, the output is saved in the more effcient qs format! See qs package.", immediate. = TRUE)

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
    print(""); print(fnameIN)

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
      qs::qsave(x = LSB, file = LSBnameOUT)

    } else {
      print("More than 2 elements in the list of matrices")
    }

    ncells <- ncol(seu)
    fname_X <- Stringendo::sppp(fnameIN, "min.cells", min.cells, "min.features", min.features,
                                "cells", ncells)
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



# _________________________________________________________________________________________________
#' @title ConvertDropSeqfolders
#'
#' @description This function takes a parent directory with a number of subfolders, each containing the standard output of 10X Cell Ranger. It (1) loads the filtered data matrices, (2) converts them to Seurat objects, and (3) saves them as .RDS files.
#' @param InputDir A character string specifying the input directory.
#' @param folderPattern A character string specifying the pattern of folder names to be searched. Default is 'SRR*'.
#' @param filePattern A character string specifying the pattern of file names to be searched. Default is 'expression.tsv.gz'.
#' @param useVroom A logical value indicating whether to use vroom. Default is TRUE.
#' @param col_types.vroom A list defining column types for vroom. Default is list("GENE" = "c", .default = "d").
#' @param min.cells An integer value specifying the minimum number of cells. Default is 10.
#' @param min.features An integer value specifying the minimum number of features. Default is 200.
#' @param updateHGNC A logical value indicating whether to update the HGNC. Default is TRUE.
#' @param ShowStats A logical value indicating whether to show statistics. Default is TRUE.
#' @param minDimension An integer value specifying the minimum dimension. Default is 10.
#' @param overwrite A logical value indicating whether to overwrite files. Default is FALSE.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   ConvertDropSeqfolders(InputDir = InputDir)
#' }
#' }
#' @seealso
#'  \code{\link[vroom]{vroom}}
#'  \code{\link[readr]{read_delim}}
#' @export
#' @importFrom vroom vroom
#' @importFrom readr read_tsv
ConvertDropSeqfolders <- function(
    InputDir # Take a parent directory with a number of subfolders, each containing the standard output of 10X Cell Ranger. (1.) It loads the filtered data matrices; (2.) converts them to Seurat objects, and (3.) saves them as *.RDS files.
    , folderPattern = "SRR*", filePattern = "expression.tsv.gz",
    useVroom = TRUE, col_types.vroom = list("GENE" = "c", .default = "d"),
    min.cells = 10, min.features = 200, updateHGNC = TRUE, ShowStats = TRUE, minDimension = 10, overwrite = FALSE) {
  InputDir <- FixPath(InputDir)
  fin <- list.dirs(InputDir, recursive = FALSE)
  fin <- CodeAndRoll2::grepv(x = fin, pattern = folderPattern, perl = FALSE)

  for (i in 1:length(fin)) {
    print(i)
    pathIN <- FixPath(fin[i])
    print(pathIN)
    fnameIN <- basename(fin[i])
    subdir <- paste0(InputDir, fnameIN)
    fnameOUT <- ppp(subdir, "min.cells", min.cells, "min.features", min.features, "Rds")
    print(fnameOUT)
    if (!overwrite) {
      OutFile <- list.files(InputDir, pattern = basename(fnameOUT), recursive = TRUE)
      if (length(OutFile) > 0) {
        if (grepl(pattern = ".Rds$", OutFile, perl = TRUE)) {
          iprint("      RDS OBJECT ALREADY EXISTS.")
          next
        }
      } # if length
    }
    CountTable <- list.files(subdir, pattern = filePattern, recursive = FALSE)
    stopifnot(length(CountTable) == 1)
    count_matrix <- if (useVroom) {
      vroom::vroom(file = kpps(subdir, CountTable), col_types = col_types.vroom)
    } else {
      readr::read_tsv(file = kpps(subdir, CountTable))
    }

    if (nrow(count_matrix) < minDimension | ncol(count_matrix) < minDimension) {
      iprint("")
      iprint("      EXPRESSION MATRIX TOO SMALL.", nrow(count_matrix), "x", ncol(count_matrix), ". Not processed.")
    } else {
      count_matrix <- FirstCol2RowNames(count_matrix)[, -1] # remove 1st "Cell column" # https://github.com/vertesy/SEO/issues/63
      seu <- CreateSeuratObject(
        counts = count_matrix, project = fnameIN,
        min.cells = min.cells, min.features = min.features
      )
      if (ncol(seu) < 1000) print("Only", ncol(seu), "cells survived filtering in the Seurat obj!")
      if (nrow(seu) < 1000) print("Only", nrow(seu), "genes found in the Seurat obj!")

      # update HGNC --- --- --- --- ---
      Sys.setenv("R_MAX_VSIZE" = 32000000000)
      if (updateHGNC) seu <- UpdateGenesSeurat(seu, EnforceUnique = TRUE, ShowStats = TRUE)
      saveRDS(seu, file = fnameOUT)
    }
  }
}


# _________________________________________________________________________________________________
#' @title LoadAllSeurats
#'
#' @description This function loads all Seurat objects found in a directory. It also works with symbolic links (but not with aliases).
#' @param InputDir A character string specifying the input directory.
#' @param file.pattern A character string specifying the pattern of file names to be searched. Default is '^filtered.+Rds$'.
#' @param string.remove1 A character string or FALSE. If a string is provided, it is removed from file names. Default is "filtered_feature_bc_matrix.".
#' @param string.replace1 A character string of the new text instead of "string.remove1".
#' @param string.remove2 A character string or FALSE. If a string is provided, it is removed from file names. Default is ".min.cells.10.min.features.200.Rds".
#' @param sort_alphanumeric sort files alphanumeric? Default: TRUE.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   ls.Seurat <- LoadAllSeurats(InputDir = InputDir)
#' }
#' }
#' @export
#' @importFrom tictoc tic toc
LoadAllSeurats <- function(
    InputDir,
    file.pattern = "^filtered.+Rds$",
    string.remove1 = list(FALSE, "filtered_feature_bc_matrix.", "raw_feature_bc_matrix.")[[2]],
    string.replace1 = "",
    string.remove2 = list(FALSE, ".min.cells.10.min.features.200.Rds")[[2]],
    sort_alphanumeric = TRUE) {

  tictoc::tic()
  InputDir <- FixPath(InputDir)

  print(file.pattern)
  use_rds <- grepl(pattern = "Rds", x = file.pattern) && !grepl(pattern = "qs", x = file.pattern)
  print(use_rds)

  fin.orig <- list.files(InputDir, include.dirs = FALSE, pattern = file.pattern)
  print(fin.orig)
  print(length(fin.orig))
  stopifnot(length(fin.orig)>0)
  fin <- if (!isFALSE(string.remove1)) sapply(fin.orig, gsub, pattern = string.remove1, replacement = string.replace1) else fin.orig
  fin <- if (!isFALSE(string.remove2)) sapply(fin, gsub, pattern = string.remove2, replacement = "") else fin
  if (sort_alphanumeric) fin <- gtools::mixedsort(fin)


  ls.Seu <- list.fromNames(fin)
  for (i in 1:length(fin)) {
    print(fin[i])
    FNP <- paste0(InputDir, fin.orig[i])
    # print(paste("Attempting to load file:", FNP))  # Debug print

    if (use_rds) {
      ls.Seu[[i]] <-readRDS(FNP)
    } else if (!use_rds) {
      ls.Seu[[i]] <-qs::qread(file = FNP)
    } else {
      warning("File pattern ambigous. Use either qs or rds:", file.pattern, immediate. = TRUE)
    }
  } # for
  print(tictoc::toc())
  return(ls.Seu)
}




# _________________________________________________________________________________________________
#' @title read10x
#'
#' @description This function reads a 10X dataset from gzipped matrix.mtx, features.tsv and barcodes.tsv files.
#' @param dir A character string specifying the directory where the gzipped files are located.
#' @return A Seurat object constructed from the 10X dataset.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   seuratObject <- read10x(dir = dir)
#' }
#' }
#' @export
#' @seealso
#'  \code{\link[tictoc]{tic}}
#'  \code{\link[R.utils]{compressFile}}
#'  \code{\link[Seurat]{Read10X}}
#' @importFrom tictoc tic toc
#' @importFrom R.utils gunzip gzip
#' @importFrom Seurat Read10X
read10x <- function(dir) { # read10x from gzipped matrix.mtx, features.tsv and barcodes.tsv
  tictoc::tic()
  names <- c("barcodes.tsv", "features.tsv", "matrix.mtx")
  for (i in 1:length(names)) {
    R.utils::gunzip(paste0(dir, "/", names[i], ".gz"))
  }
  file.copy(paste0(dir, "/features.tsv"), paste0(dir, "/genes.tsv"))
  mat <- Seurat::Read10X(dir)
  file.remove(paste0(dir, "/genes.tsv"))
  for (i in 1:length(names)) {
    R.utils::gzip(paste0(dir, "/", names[i]))
  }
  tictoc::toc()
  mat
}



# _________________________________________________________________________________________________
#' @title load10Xv3
#'
#' @description Load 10X output folders.
#' @param dataDir A character string specifying the directory containing the 10X output folders.
#' @param cellIDs A vector specifying the cell IDs. Default: NULL.
#' @param channelName A character string specifying the channel name. Default: NULL.
#' @param readArgs A list of arguments to pass to the internal Read10X function. Default: list().
#' @param includeFeatures A character vector specifying which features to include. Default: c("Gene Expression").
#' @param verbose A logical flag indicating whether to print status messages. Default: TRUE.
#' @param ... Additional arguments to pass to the internally called functions.
#' @return An object of class "SoupChannel" representing the loaded 10X data.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   channel <- load10Xv3(dataDir = dataDir)
#' }
#' }
#' @seealso
#'  \code{\link[SoupX]{SoupChannel}}
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
#' @title .saveRDS.compress.in.BG
#'
#' @description Save and RDS object and compress resulting file in the background using system(gzip). OS X or unix.
#' @param obj Seurat object.
#' @param compress_internally Compress by R? Default: FALSE (still compressed in background via CLI).
#' @param compr Compress at all? Default: TRUE
#' @param fname File name
#' @param ... Additional parameters passed to saveRDS() function.
#' @seealso
#'  \code{\link[tictoc]{tic}}
#' @importFrom tictoc tic toc
.saveRDS.compress.in.BG <- function(obj, compr = FALSE, fname, compress_internally = FALSE, ...) {
  try(tictoc::tic(), silent = TRUE)
  saveRDS(object = obj, compress = compress_internally, file = fname, ...)
  try(tictoc::toc(), silent = TRUE)
  if (compr) system(command = paste0("gzip '", fname, "'"), wait = FALSE) # execute in the background
  print(paste("Saved, optionally being .gz compressed", fname))
  try(say(), silent = TRUE)
}




# _________________________________________________________________________________________________
#' @title isave.RDS
#'
#' @description Save an RDS object, using a faster and efficient compression method that runs in the background.
#' @param obj The object to be saved, typically a Seurat object.
#' @param prefix A string prefix added to the filename. Default is NULL.
#' @param suffix A string suffix added to the filename. Default is NULL.
#' @param inOutDir A boolean flag, if TRUE the OutDir is used as save directory, if FALSE the alternative_path_rdata is used. Default is TRUE
#' @param project A string representing the project code. This is appended to the saved file name. Default is the active project determined by getProject().
#' @param alternative_path_rdata A string that specifies the alternative path for storing the RDS file if inOutDir is FALSE. Default is "~/Dropbox (VBC)/Abel.IMBA/AnalysisD/_RDS.files/" appended with the basename of OutDir.
#' @param homepath A string representing the homepath. Will be replaced by '~' in the file path. Default is '~/'.
#' @param showMemObject A boolean flag, if TRUE the function will print out the memory size of the largest objects in the workspace. Default is TRUE.
#' @param saveParams A boolean flag, if TRUE the parameters 'p' and 'all.genes' are added to the 'misc' slot of the Seurat object if the object is of class Seurat. Default is TRUE.
#' @param compress Compress .Rds file after writing? Default is TRUE.
#' @param test_read Provide command to test validity by reading in the object just written.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   isave.RDS(my.R.object)
#' }
#' }
#' @export
isave.RDS <- function(
    obj, prefix = NULL, suffix = NULL, inOutDir = TRUE,
    project = getProject(),
    alternative_path_rdata = paste0("~/Dropbox (VBC)/Abel.IMBA/AnalysisD/_RDS.files/", basename(OutDir)),
    homepath = if (Sys.info()[1] == "Darwin") "~/" else "/users/abel.vertesy/",
    showMemObject = TRUE, saveParams = TRUE,
    compress = TRUE,
    test_read = FALSE) {
  warning("isave.RDS() is deprecated. Use xsave() to save in .qs format.", immediate. = TRUE)
  path_rdata <- if (inOutDir) OutDir else alternative_path_rdata
  dir.create(path_rdata)

  if (showMemObject) {
    try(memory.biggest.objects(), silent = TRUE)
  }
  if ("seurat" %in% is(obj) & saveParams) {
    try(obj@misc$p <- p, silent = TRUE)
    try(obj@misc$all.genes <- all.genes, silent = TRUE)
  }
  fnameBase <- kppu(prefix, substitute(obj), project, suffix, idate(Format = "%Y.%m.%d_%H.%M"))
  fnameBase <- trimws(fnameBase, whitespace = "_")
  FNN <- paste0(path_rdata, fnameBase, ".Rds")
  FNN <- gsub(pattern = "~/", replacement = homepath, x = FNN)
  print(FNN)
  if (test_read) {
    print(paste0('xx5 <- read_rds(\\"', FNN, '\\")'))
  } else {
    Seurat.utils:::.saveRDS.compress.in.BG(obj = obj, fname = FNN, compr = compress, compress_internally = FALSE)
  }
}

# _________________________________________________________________________________________________
#' @title Save an R Object Using 'qs' Package for Fast Compressed Saving
#'
#' @description This function saves an R object to a file in a quick and efficient format using the 'qs' package.
#' It constructs the file name based on various inputs and stores additional metadata if the object is a Seurat object.
#' The saving path can be adjusted by the presence of 'OutDir' in the global environment or defaults to the working directory.
#'
#' @param obj The R object to be saved.
#' @param prefix Optional; a prefix to add to the filename.
#' @param suffix Optional; a suffix to add to the filename.
#' @param nthreads Number of threads to use when saving, defaults to 12.
#' @param preset Compression preset, defaults to 'high'.
#' @param project The project name to be included in the filename, defaults to the result of `getProject()`.
#' @param out_dir Output Directory
#' @param background_job Logical; if TRUE save runs as "background job"
#' @param showMemObject Logical; if TRUE, displays the memory size of the largest objects.
#' @param saveParams Logical; if TRUE and if the object is a Seurat object, additional parameters are saved within it.
#'
#' @return Invisible; The function is called for its side effects (saving a file) and does not return anything.
#'
#' @note The function uses the 'qs' package for quick and efficient serialization of objects and includes a timing feature from the 'tictoc' package.
#' @seealso \code{\link[qs]{qsave}} for the underlying save function used.
#' @importFrom qs qsave
#' @importFrom tictoc tic toc
#' @importFrom job job
#' @importFrom rstudioapi isAvailable
#'
#' @export
xsave <- function(
    obj, prefix = NULL,
    suffix = NULL,
    nthreads = 12,
    preset = "high",
    project = getProject(),
    out_dir = if (exists("OutDir")) OutDir else getwd(),
    background_job = FALSE,
    showMemObject = TRUE, saveParams = TRUE) {


  try(tictoc::tic(), silent = TRUE)
  if (showMemObject) {
    try(memory.biggest.objects(), silent = TRUE)
  }

  if ("seurat" %in% is(obj) & saveParams) {
    try(obj@misc$p <- p, silent = TRUE)
    try(obj@misc$all.genes <- all.genes, silent = TRUE)
  }

  fnameBase <- trimws(kppu(prefix, substitute(obj), suffix, project, preset, "compr", idate(Format = "%Y.%m.%d_%H.%M")), whitespace = "_")
  FNN <- paste0(out_dir, fnameBase, ".qs")
  iprint(substitute(obj), '<- xread("', FNN, '")')

  if (background_job & rstudioapi::isAvailable()) {
    "This part is not debugged yet!"

    message("Started saving as background job.")
    job::job(
      {
        qs::qsave(x = obj, file = FNN, nthreads = nthreads, preset = preset)
      },
      import = c("obj", "FNN", "nthreads", "preset")
    )
  } else {
    qs::qsave(x = obj, file = FNN, nthreads = nthreads, preset = preset)
  }

  try(tictoc::toc(), silent = TRUE)
}

# _________________________________________________________________________________________________
#' @title Read an R Object Using 'qs' Package for Fast Decompression
#'
#' @description This function reads an R object from a file saved in a format specific to the 'qs' package,
#' which is designed for quick and efficient compression and decompression of R objects.
#' It also times the read operation, providing feedback on the duration of the operation.
#'
#' @param file A character string specifying the path to the file where the R object is saved.
#' @param nthreads The number of threads to use when reading the object, defaults to 4.
#' @param ... Further arguments passed on to the 'qs::qread' function.
#' @return The R object that was saved in the specified file.
#' @note The function uses the 'qs' package for fast and efficient deserialization of objects
#' and includes a timing feature from the 'tictoc' package.
#' @seealso \code{\link[qs]{qread}} for the underlying read function used.
#' @importFrom qs qread
#' @importFrom tictoc tic toc
#' @importFrom job job
#' @importFrom rstudioapi isAvailable
#'
#' @export
xread <- function(file, nthreads = 4, ...) {
  stopifnot(file.exists(file))
  try(tictoc::tic(), silent = TRUE)

  # if (background_job & rstudioapi::isAvailable()) {
  #   "This part is not debugged yet!"
  #   "This part is not debugged yet!"
  #
  #   message("Started reading in as background job.")
  #   job::job(
  #     {
  #       qs::qread(file = file, nthreads = nthreads, ...)
  #     },
  #     import = c("file", "nthreads")
  #   )
  # } else {
    x <- qs::qread(file = file, nthreads = nthreads, ...)
  # }

  iprint(is(x)[1], "of length:", length(x))
  try(tictoc::toc(), silent = TRUE)
  invisible(x)
}



# _________________________________________________________________________________________________
# Save workspace
# requires MarkdownReports (github) and defining OutDir
# requires github/vertesy/CodeAndRoll.r

#' @title isave.image
#'
#' @description Save an image of the current workspace using a faster and efficient compression method that runs in the background.
#' @param ... Additional parameters passed to the idate() function in the creation of the file name.
#' @param path_rdata A string that specifies the path for storing the image of the workspace. Default is "~/Dropbox/Abel.IMBA/AnalysisD/_Rdata.files/" appended with the basename of OutDir.
#' @param showMemObject A boolean flag, if TRUE the function will print out the memory size of the largest objects in the workspace. Default is TRUE.
#' @param options A string for gzip options. Default is "--force".
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   isave.image(my.R.image)
#' }
#' }
#' @export
#' @importFrom Stringendo kollapse iprint
isave.image <- function(
    ..., path_rdata = paste0("~/Dropbox/Abel.IMBA/AnalysisD/_Rdata.files/", basename(OutDir)),
    showMemObject = TRUE, options = c("--force", NULL)[1]) { # Faster saving of workspace, and compression outside R, when it can run in the background. Seemingly quite CPU hungry and not very efficient compression.

  dir.create(path_rdata)

  if (showMemObject) {
    try(memory.biggest.objects(), silent = TRUE)
  }
  fname <- Stringendo::kollapse(path_rdata, "/", idate(), ..., ".Rdata")
  print(fname)
  if (nchar(fname) > 2000) stop()
  save.image(file = fname, compress = FALSE)
  iprint("Saved, being compressed", fname)
  system(paste("gzip", options, fname), wait = FALSE) # execute in the background
}


# _________________________________________________________________________________________________
#' @title Save workspace - qsave.image
#'
#' @description Faster saving of workspace, and compression outside R, when it can run in the background. Seemingly quite CPU hungry and not very efficient compression. #
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @param options Options passed on to gzip, via CLI. Default: c("--force", NULL)[1]
#' @seealso
#'  \code{\link[Stringendo]{kollapse}}, \code{\link[function]{iprint}}
#' @export
#' @importFrom Stringendo kollapse iprint
#' @importFrom tictoc tic toc
qsave.image <- function(..., showMemObject = TRUE, options = c("--force", NULL)[1]) { # Faster saving of workspace, and compression outside R, when it can run in the background. Seemingly quite CPU hungry and not very efficient compression.
  fname <- Stringendo::kollapse(getwd(), "/", basename(OutDir), idate(), ..., ".Rdata")
  print(fname)
  if (nchar(fname) > 2000) stop()
  tictoc::tic()
  save.image(file = fname, compress = FALSE)
  iprint("Saved, being compressed", fname)
  system(paste("gzip", options, fname), wait = FALSE) # execute in the background
  cat(tictoc::toc)
}


# _________________________________________________________________________________________________
#' @title Find 'Outs' Subdirectories in Specified Subdirectories
#'
#' @description This function searches through specified subdirectories within a root directory
#' to find all subdirectories named 'outs' and returns a character vector with their full paths.
#'
#' @param root_dir The root directory.
#' @param subdir A character vector of subdirectory names within the root directory to be scanned.
#' @param recursive Boolean indicating whether to search recursively within subdirectories.
#' @return A character vector containing the full paths to the 'outs' subdirectories.
#' @importFrom fs dir_ls
#' @export
find10XoutputFolders <- function(root_dir, subdir, recursive = TRUE) {
  stopifnot(is.character(root_dir), length(root_dir) == 1, dir.exists(root_dir),
            is.character(subdir), all(dir.exists(file.path(root_dir, subdir))),
            is.logical(recursive))

  outs_dirs <- c()
  for (i in seq_along(subdir)) {
    path <- file.path(root_dir, subdir[i])
    printProgress(i, length(subdir), "Processing subdirectory")

    iprint("Searching in:", path)
    found_dirs <- fs::dir_ls(path, recurse = recursive, glob = "*/outs", type = "directory")
    iprint(length(found_dirs), "output folders found.")
    outs_dirs <- c(outs_dirs, found_dirs)
  }

  # Replace root_dir in the paths with an empty string for printing
  outs_print <- gsub(paste0("^", root_dir, "/?"), "", outs_dirs)
  iprint(length(outs_dirs), outs_print)

  return(outs_dirs)
}

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
#' @title clip10Xcellname
#'
#' @description Clip all suffices after underscore (10X adds it per chip-lane, Seurat adds in during integration). #
#' @param cellnames Character vector containing the cellIDs (with numeric suffixes).
#' @export
#' @importFrom stringr str_split_fixed
clip10Xcellname <- function(cellnames) {
  stringr::str_split_fixed(cellnames, "_", n = 2)[, 1]
}

# _________________________________________________________________________________________________
#' @title make10Xcellname
#'
#' @description Add a suffix to cell names, so that it mimics the lane-suffix, e.g.: "_1". #
#' @param cellnames Character vector containing the cellIDs (WITHOUT numeric suffixes).
#' @param suffix A suffix added to the filename, Default: '_1'
#' @export
make10Xcellname <- function(cellnames, suffix = "_1") {
  paste0(cellnames, suffix)
}




# _________________________________________________________________________________________________
# Soup.Analysis.of.ambient.RNA.R ______________________________ ----
# _________________________________________________________________________________________________
# source('~/GitHub/Packages/Seurat.utils/Functions/Soup.Analysis.of.ambient.RNA.R')
# try (source('https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Soup.Analysis.of.ambient.RNA.R'))
# Source: self + web



# _________________________________________________________________________________________________
#' @title plotTheSoup
#'
#' @description Plot stats about the ambient RNA content in a 10X experiment.
#' @param CellRanger_outs_Dir CellRanger 'outs' (output) directory, Default: '~/Data/114593/114593'
#' @param SeqRun Aka SampleName (the folder above 'outs;). Default: str_extract(CellRanger_outs_Dir, "[[:alnum:]_]+(?=/outs/)").
#' @seealso
#'  \code{\link[Matrix]{colSums}}
#'  \code{\link[tibble]{rownames}}
#'  \code{\link[ggrepel]{geom_label_repel}}
#' @importFrom Matrix rowSums
#' @importFrom tibble rownames_to_column
#' @importFrom ggrepel geom_text_repel
#' @importFrom Stringendo percentage_formatter
#' @importFrom MarkdownReports wbarplot create_set_OutDir
#' @importFrom MarkdownHelpers ww.assign_to_global
#' @importFrom dplyr as_tibble
#'
#' @export
plotTheSoup <- function(CellRanger_outs_Dir = "~/Data/114593/114593",
                        SeqRun = str_extract(CellRanger_outs_Dir, "[[:alnum:]_]+(?=/outs/)"),
                        ls.Alpha = 1) {
  # The regular expression `[[:alnum:]_]+(?=/outs/)` matches one or more alphanumeric characters or
  # underscores that are followed by the `/outs/` portion in the string. It ensures that the desired
  # substring is captured, but it does not include the `/outs/` in the matched result.
  # `[[:alnum:]_]+` matches one or more alphanumeric characters or underscores. The `[:alnum:]`
  # character class represents all alphabetic characters (both uppercase and lowercase) and digits.
  # The underscore character `_` is included as well. The `+` quantifier specifies that there should
  # be one or more occurrences of these characters in a row.
  # `(?=/outs/)` This is a positive lookahead assertion. It matches a position in the string where
  # `/outs/` is present immediately after. It doesn't consume any characters from the string; it
  # just checks for the presence of `/outs/` after the matched substring.

  # Check the directory exists and SeqRun ID is longer than 4
  stopifnot(dir.exists(CellRanger_outs_Dir))
  iprint("SeqRun", SeqRun)
  stopifnot(nchar(SeqRun) > 4)
  Subfolders_10X_outs <- list.dirs(CellRanger_outs_Dir, full.names = FALSE, recursive = FALSE)
  stopifnot(length(Subfolders_10X_outs) > 0)

  # Identify raw and filtered files ___________________________________
  path.raw <- file.path(CellRanger_outs_Dir, grep(x = Subfolders_10X_outs, pattern = "^raw_*", value = TRUE))
  path.filt <- file.path(CellRanger_outs_Dir, grep(x = Subfolders_10X_outs, pattern = "^filt_*", value = TRUE))
  CR.matrices <- list.fromNames(c("raw", "filt"))

  # Adapter for Markdownreports background variable "OutDir"
  OutDirBac <- if (exists("OutDir")) OutDir else getwd()
  OutDir <- file.path(CellRanger_outs_Dir, paste0(kpp("SoupStatistics", SeqRun)))
  MarkdownReports::create_set_OutDir(OutDir)
  MarkdownHelpers::ww.assign_to_global("OutDir", OutDir, 1)

  # Read raw and filtered data ___________________________________
  print("Reading raw CellRanger output matrices")
  CR.matrices$"raw" <- Seurat::Read10X(path.raw)
  if (length(CR.matrices$"raw") == 2) {
    CR.matrices$"raw" <- CR.matrices$"raw"[[1]]
  } # Maybe AB table is present too at slot 2!
  print("Reading filtered CellRanger output matrices")
  CR.matrices$"filt" <- Seurat::Read10X(path.filt)
  if (length(CR.matrices$"filt") == 2) {
    CR.matrices$"filt" <- CR.matrices$"filt"[[1]]
  } # Maybe AB table is present too at slot 2!

  # Profiling the soup ___________________________________
  print("Profiling the soup")
  GEMs.all <- CR.matrices$"raw"@Dimnames[[2]]
  GEMs.cells <- CR.matrices$"filt"@Dimnames[[2]]
  iprint("There are", length(GEMs.all), "GEMs sequenced, and", length(GEMs.cells), "are cells among those.")
  EmptyDroplets.and.Cells <- c("EmptyDroplets" = length(GEMs.all) - length(GEMs.cells), "Cells" = length(GEMs.cells))
  ggExpress::qbarplot(EmptyDroplets.and.Cells, label = EmptyDroplets.and.Cells, palette_use = "npg", col = 1:2, ylab = "GEMs")

  GEMs.soup <- setdiff(GEMs.all, GEMs.cells)
  CR.matrices$"soup" <- CR.matrices$"raw"[, GEMs.soup]
  CR.matrices$"soup.total.RC" <- Matrix::rowSums(CR.matrices$"soup")
  CR.matrices$"soup.total.sum" <- sum(CR.matrices$"soup")
  CR.matrices$"cells.total.sum" <- sum(CR.matrices$"filt")

  CR.matrices$"soup.rel.RC" <- CR.matrices$"soup.total.RC" / CR.matrices$"soup.total.sum"

  # Diff Exp ___________________________________
  Soup.VS.Cells.Av.Exp <- cbind(
    "Soup" = Matrix::rowSums(CR.matrices$"soup"),
    "Cells" = Matrix::rowSums(CR.matrices$"filt")
  )
  colnames(Soup.VS.Cells.Av.Exp)
  idx.HE <- rowSums(Soup.VS.Cells.Av.Exp) > 10
  pc_TRUE(idx.HE)
  Soup.VS.Cells.Av.Exp <- Soup.VS.Cells.Av.Exp[idx.HE, ]
  idim(Soup.VS.Cells.Av.Exp)
  Soup.VS.Cells.Av.Exp.log10 <- log10(Soup.VS.Cells.Av.Exp + 1)

  # ggplot prepare ___________________________________
  Soup.VS.Cells.Av.Exp.gg <- tibble::rownames_to_column(as.data.frame(Soup.VS.Cells.Av.Exp.log10), "gene")
  (Soup.VS.Cells.Av.Exp.gg <- dplyr::as_tibble(Soup.VS.Cells.Av.Exp.gg))
  soup.rate <- Soup.VS.Cells.Av.Exp.gg$Soup / (Soup.VS.Cells.Av.Exp.gg$Cells + Soup.VS.Cells.Av.Exp.gg$Soup)
  cell.rate <- Soup.VS.Cells.Av.Exp.gg$Cells / (Soup.VS.Cells.Av.Exp.gg$Cells + Soup.VS.Cells.Av.Exp.gg$Soup)

  axl.pfx <- "Total Expression in"
  axl.sfx <- "[log10(mRNA+1)]"

  HGNC <- Soup.VS.Cells.Av.Exp.gg$gene
  Class <- rep("Other", times = nrow(Soup.VS.Cells.Av.Exp.gg))
  Class[grep("^RPL|^RPS", HGNC)] <- "RP"
  Class[grep("^MT-", HGNC)] <- "MT"
  Class[grep("^LINC", HGNC)] <- "LINC"
  Class[grep("^AC", HGNC)] <- "AC"
  Class[grep("^AL", HGNC)] <- "AL"
  Fraction.of.Geneclasses <- table(Class)
  ggExpress::qpie(Fraction.of.Geneclasses)
  # wpie(Fraction.of.Geneclasses)
  Soup.VS.Cells.Av.Exp.gg$Class <- Class

  fname <- kpp("Soup.VS.Cells.Av.Exp.GeneClasses", SeqRun, "pdf")
  pgg <-
    ggplot(
      Soup.VS.Cells.Av.Exp.gg %>% arrange(-nchar(Class)),
      aes(x = Soup, y = Cells, label = gene, col = Class)
    ) +
    geom_abline(slope = 1, col = "darkgrey") +
    geom_point() +
    scale_alpha_manual(guide = "none", values = ls.Alpha) +
    xlab(paste(axl.pfx, "Soup", axl.sfx)) +
    ylab(paste(axl.pfx, "Cells", axl.sfx)) +
    ggtitle("Soup VS. Cells | gene classes")

  ggsave(pgg, filename = file.path(OutDir, fname))

  # ggplot ___________________________________
  quantiles <- c(0.025, 0.01, 0.0025)

  i <- 1
  for (i in 1:length(quantiles)) {
    pr <- quantiles[i]
    print(pr)
    HP.thr <- 200 * pr / quantiles[2]
    idx.HE2 <- rowSums(Soup.VS.Cells.Av.Exp) > HP.thr
    pc_TRUE(idx.HE2)

    fname <- kpp("Soup.VS.Cells.Av.Exp.quantile", pr, SeqRun, "pdf")

    Outlier <- idx.HE2 &
      (cell.rate < quantile(cell.rate, probs = pr) |
        soup.rate < quantile(soup.rate, probs = pr))

    pc_TRUE(Outlier)
    sum(Outlier)
    HP.thr.mod <- HP.thr
    while (sum(Outlier) > 40) {
      HP.thr.mod <- HP.thr.mod * 2
      Outlier <- Outlier & rowSums(Soup.VS.Cells.Av.Exp) > HP.thr.mod
    }



    pgg <-
      ggplot(Soup.VS.Cells.Av.Exp.gg, aes(
        x = Soup, y = Cells, label = gene,
        col = Outlier
      )) +
      geom_point() +
      theme(legend.position = "none") +
      xlab(paste(axl.pfx, "Soup", axl.sfx)) +
      ylab(paste(axl.pfx, "Cells", axl.sfx)) +
      ggtitle("Soup VS. Cells", subtitle = pr) +
      ggrepel::geom_text_repel(aes(label = ifelse(Outlier,
        as.character(gene), ""
      )))
    ggsave(pgg, filename = file.path(OutDir, fname))
  }


  # Per Gene ___________________________________
  PC.mRNA.in.Soup <- sum(CR.matrices$"soup") / sum(CR.matrices$"raw")
  PC.mRNA.in.Cells <- 100 * sum(CR.matrices$"filt") / sum(CR.matrices$"raw")
  MarkdownReports::wbarplot(
    variable = PC.mRNA.in.Cells, col = "seagreen", plotname = kppd("PC.mRNA.in.Cells", SeqRun),
    ylim = c(0, 100), ylab = "% mRNA in cells",
    sub = "% mRNA is more meaningful than % reads reported by CR"
  )
  barplot_label(
    barplotted_variable = PC.mRNA.in.Cells,
    labels = Stringendo::percentage_formatter(PC.mRNA.in.Cells / 100, digitz = 2),
    TopOffset = 10
  )


  # Plot top gene's expression ___________________________________
  Soup.GEMs.top.Genes <- 100 * head(sort(CR.matrices$"soup.rel.RC", decreasing = TRUE), n = 20)

  MarkdownReports::wbarplot(Soup.GEMs.top.Genes,
    plotname = kppd("Soup.GEMs.top.Genes", SeqRun),
    ylab = "% mRNA in the Soup",
    sub = paste("Within the", SeqRun, "dataset"),
    tilted_text = TRUE,
    ylim = c(0, max(Soup.GEMs.top.Genes) * 1.5)
  )
  barplot_label(
    barplotted_variable = Soup.GEMs.top.Genes,
    labels = Stringendo::percentage_formatter(Soup.GEMs.top.Genes / 100, digitz = 2),
    TopOffset = -.5, srt = 90, cex = .75
  )

  # Plot summarize expression ___________________________________
  soupProfile <- CR.matrices$"soup.total.RC"
  {
    soup.RP.sum <- sum(soupProfile[grep("^RPL|^RPS", names(soupProfile))])
    soup.RPL.sum <- sum(soupProfile[grep("^RPL", names(soupProfile))])
    soup.RPS.sum <- sum(soupProfile[grep("^RPS", names(soupProfile))])
    soup.mito.sum <- sum(soupProfile[grep("^MT-", names(soupProfile))])
    soup.LINC.sum <- sum(soupProfile[grep("^LINC", names(soupProfile))])
    soup.AC.sum <- sum(soupProfile[grep("^AC", names(soupProfile))])
    soup.AL.sum <- sum(soupProfile[grep("^AL", names(soupProfile))])
    genes.non.Above <- soupProfile[CodeAndRoll2::grepv("^RPL|^RPS|^MT-|^LINC|^AC|^AL", names(soupProfile), invert = TRUE)]
  }
  head(sort(genes.non.Above), n = 50)


  soupProfile.summarized <- c(
    "Mitochondial" = soup.mito.sum,
    "Ribosomal" = soup.RP.sum,
    "Ribosomal.L" = soup.RPL.sum,
    "Ribosomal.S" = soup.RPS.sum,
    "GenBank (AC)" = soup.AC.sum,
    "EMBL (AL)" = soup.AL.sum,
    "LINC" = soup.LINC.sum,
    sort(genes.non.Above, decreasing = TRUE)
  )
  NrColumns2Show <- min(10, nrow(soupProfile.summarized))
  ccc <- c("#FF4E00", "#778B04", "#8ea604", "#8ea604", "#F5BB00", "#F5BB00", "#EC9F05", rep(x = "#BF3100", times = NrColumns2Show - 6)) # ,"#"


  Soup.GEMs.top.Genes.summarized <- 100 * soupProfile.summarized[1:NrColumns2Show] / CR.matrices$"soup.total.sum"
  maxx <- max(Soup.GEMs.top.Genes.summarized)
  MarkdownReports::wbarplot(Soup.GEMs.top.Genes.summarized,
    plotname = kppd("Soup.GEMs.top.Genes.summarized", SeqRun),
    ylab = "% mRNA in the Soup", ylim = c(0, maxx + 3),
    sub = paste("Within the", SeqRun, "dataset"),
    tilted_text = TRUE, col = ccc
  )
  barplot_label(
    barplotted_variable = Soup.GEMs.top.Genes.summarized,
    srt = 45, labels = Stringendo::percentage_formatter(Soup.GEMs.top.Genes.summarized / 100, digitz = 2),
    TopOffset = -1.5
  )

  # Absolute.fraction ___________________________________
  Absolute.fraction.soupProfile.summarized <- Soup.GEMs.top.Genes.summarized * PC.mRNA.in.Soup

  maxx <- max(Absolute.fraction.soupProfile.summarized)
  MarkdownReports::wbarplot(Absolute.fraction.soupProfile.summarized,
    plotname = kppd("Absolute.fraction.soupProfile.summarized", SeqRun),
    ylab = "% of mRNA in cells", ylim = c(0, maxx * 1.33),
    sub = paste(Stringendo::percentage_formatter(PC.mRNA.in.Soup), "of mRNA counts are in the Soup, in the dataset ", SeqRun),
    tilted_text = TRUE, col = ccc
  )
  barplot_label(
    barplotted_variable = Absolute.fraction.soupProfile.summarized,
    srt = 45, labels = Stringendo::percentage_formatter(Absolute.fraction.soupProfile.summarized / 100, digitz = 2)
    # formatC(Absolute.fraction.soupProfile.summarized, format="f", big.mark = " ", digits = 0)
    , TopOffset = -maxx * 0.15
  )

  # ___________________________________
  Soup.GEMs.top.Genes.non.summarized <- 100 * sort(genes.non.Above, decreasing = TRUE)[1:20] / CR.matrices$"soup.total.sum"
  maxx <- max(Soup.GEMs.top.Genes.non.summarized)
  MarkdownReports::wbarplot(Soup.GEMs.top.Genes.non.summarized,
    plotname = kppd("Soup.GEMs.top.Genes.non.summarized", SeqRun),
    ylab = "% mRNA in the Soup",
    sub = paste("Within the", SeqRun, "dataset"),
    tilted_text = TRUE, col = "#BF3100",
    ylim = c(0, maxx * 1.5)
  )
  barplot_label(
    barplotted_variable = Soup.GEMs.top.Genes.non.summarized,
    labels = Stringendo::percentage_formatter(Soup.GEMs.top.Genes.non.summarized / 100, digitz = 2),
    TopOffset = -maxx * 0.2, srt = 90, cex = .75
  )

  if (exists("OutDirBac")) MarkdownHelpers::ww.assign_to_global("OutDir", OutDirBac, 1)
} # plotTheSoup





# _________________________________________________________________________________________________
# Jaccard.toolkit _____________________________ ----
# _________________________________________________________________________________________________
# try(source('~/GitHub/Packages/Seurat.utils/Functions/Jaccard.toolkit.R'))
# try(source('https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Jaccard.toolkit.R'))


#  __________________________________________
# Fast direct calculation from a list


# _________________________________________________________________________________________________
#' @title jJaccardIndexVec
#'
#' @description Calculate jaccard similarity for 2 vecotrs. Helper to jPairwiseJaccardIndexList.
#' @param A Set A, Default: 1:3
#' @param B Set B, Default: 2:4
#' @export
jJaccardIndexVec <- function(A = 1:3, B = 2:4) length(intersect(A, B)) / length(union(A, B))

# _________________________________________________________________________________________________
#' @title jPairwiseJaccardIndexList
#'
#' @description Create a pairwise jaccard similarity matrix across all combinations of columns in binary.presence.matrix. Modified from: https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/ #
#' @param lsG List of genes, Default: ls_genes
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   jPairwiseJaccardIndexList(lsG = ls_genes)
#' }
#' }
#' @export
#' @importFrom Stringendo percentage_formatter
jPairwiseJaccardIndexList <- function(lsG = ls_genes) { # Create a pairwise jaccard similarity matrix across all combinations of columns in binary.presence.matrix. Modified from: https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/
  if (length(names(lsG)) < length(lsG)) {
    iprint("Gene lists were not (all) named, now renamed as:")
    names(lsG) <- ppp("dataset", 1:length(lsG))
    print(names(lsG))
  }
  m <- matrix.fromNames(rowname_vec = names(lsG), colname_vec = names(lsG))
  n.sets <- length(lsG)
  for (r in 1:n.sets) {
    # print(Stringendo::percentage_formatter(r/n.sets))
    for (c in 1:n.sets) {
      if (c == r) {
        m[r, c] <- 1
      } else {
        m[r, c] <- signif(jJaccardIndexVec(lsG[[r]], lsG[[c]]), digits = 2)
      }
    }
  }
  return(m)
}


# Much slower Indirect calculation via PresenceMatrix
# _________________________________________________________________________________________________

# _________________________________________________________________________________________________
#' @title jPresenceMatrix
#'
#' @description Make a binary presence matrix from a list. Source: https://stackoverflow.com/questions/56155707/r-how-to-create-a-binary-relation-matrix-from-a-list-of-strings #
#' @param string_list List of strings to compare overlapping entries. Default: lst(a = 1:3, b = 2:5, c = 4:9, d = -1:4)
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   df.presence <- jPresenceMatrix(string_list = lst(a = 1:3, b = 2:5, c = 4:9, d = -1:4))
#' }
#' }
#' @export
jPresenceMatrix <- function(string_list = lst(a = 1:3, b = 2:5, c = 4:9, d = -1:4)) { # Make a binary presence matrix from a list. Source: https://stackoverflow.com/questions/56155707/r-how-to-create-a-binary-relation-matrix-from-a-list-of-strings
  df.presence <- string_list %>%
    enframe() %>%
    unnest(cols = "value") %>%
    count(name, value) %>%
    spread(value, n, fill = 0)
  df.presence2 <- FirstCol2RowNames(df.presence)
  return(t(df.presence2))
}


# _________________________________________________________________________________________________
#' @title jJaccardIndexBinary
#'
#' @description Calculate Jaccard Index. Modified from: https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/ #
#' @param x Set X
#' @param y Set Y
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   JaccardSimilarity <- jJaccardIndexBinary(
#'     x = sample(x = 0:1, size = 100, replace = TRUE),
#'     y = sample(x = 0:1, size = 100, replace = TRUE)
#'   )
#' }
#' }
#' @export
jJaccardIndexBinary <- function(x, y) { # Calculate Jaccard Index. Modified from: https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/
  elements.found <- sort(unique(union(x, y)))
  stopifnot(length(elements.found) == 2) # check if you only have [0,1]
  stopifnot(as.numeric(elements.found) == 0:1) # check if you only have [0,1]

  M.11 <- sum(x == 1 & y == 1)
  M.10 <- sum(x == 1 & y == 0)
  M.01 <- sum(x == 0 & y == 1)
  return(M.11 / (M.11 + M.10 + M.01))
}




# _________________________________________________________________________________________________
#' @title jPairwiseJaccardIndex
#'
#' @description Create a pairwise jaccard similarity matrix across all combinations of columns in binary.presence.matrix. Modified from: https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/ #
#' @param binary.presence.matrix A boolean matrix. Default: df.presence
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   PairwiseJaccardIndices <- jPairwiseJaccardIndex(binary.presence.matrix = df.presence)
#' }
#' }
#' @export
#' @importFrom Stringendo percentage_formatter
jPairwiseJaccardIndex <- function(binary.presence.matrix = df.presence) { # Create a pairwise jaccard similarity matrix across all combinations of columns in binary.presence.matrix. Modified from: https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/
  m <- matrix.fromNames(rowname_vec = colnames(binary.presence.matrix), colname_vec = colnames(binary.presence.matrix))
  n.sets <- ncol(binary.presence.matrix)
  for (r in 1:n.sets) {
    print(Stringendo::percentage_formatter(r / n.sets))
    for (c in 1:n.sets) {
      if (c == r) {
        m[r, c] <- 1
      } else {
        m[r, c] <- signif(jJaccardIndexBinary(binary.presence.matrix[, r], binary.presence.matrix[, c]), digits = 2)
      }
    }
  }
  return(m)
}


# _________________________________________________________________________________________________
# New additions,  categorized _____________________________ ------
# _________________________________________________________________________________________________



#' @title Regress Out and Recalculate Seurat
#'
#' @description The function performs a series of calculations and manipulations on a Seurat object,
#' including identifying variable features, scaling data, running PCA, setting up reductions, finding neighbors,
#' and finding clusters. It optionally performs t-SNE and saves the object.
#'
#' @param obj The Seurat object.
#' @param vars.to.regress A vector of variable names to be regressed out.
#' @param suffix A character string to be used as a suffix when saving the object.
#' @param nPCs The number of principal components to use. Default is the 'n.PC' element from a list 'p'.
#' @param clust_resolutions The resolution for clustering. Default is the 'snn_res' element from a list 'p'.
#' @param calc_tSNE Logical, if TRUE, t-SNE will be performed. Default is FALSE.
#' @param plot_umaps Logical, if TRUE, UMAP plots will be generated. Default is TRUE.
#' @param save_obj Logical, if TRUE, the object will be saved. Default is TRUE.
#' @param assayX The assay to be used in scaling data. Default is 'RNA'.
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
    vars.to.regress,
    suffix,
    nPCs = p$"n.PC",
    clust_resolutions = p$"snn_res",
    calc_tSNE = FALSE,
    plot_umaps = TRUE,
    save_obj = TRUE,
    assayX = "RNA") {
  tictoc::tic()
  print("FindVariableFeatures")
  obj <- FindVariableFeatures(obj, mean.function = "FastExpMean", dispersion.function = "FastLogVMR", nfeatures = 10000)
  tictoc::toc()

  tictoc::tic()
  print("calc.q99.Expression.and.set.all.genes")
  obj <- calc.q99.Expression.and.set.all.genes(obj = obj, quantileX = .99)
  tictoc::toc()

  tictoc::tic()
  print("ScaleData")
  obj <- ScaleData(obj, assay = assayX, verbose = TRUE, vars.to.regress = vars.to.regress)
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
# Temp _____________________________ ------

"THIS SHOULD BE MOVED"
# _________________________________________________________________________________________________
#' @title getProject
#'
#' @description Try to get the project name you are wokring on in Rstudio.
#' @returns The final subfolder of your project, or NULL, if you are not running one
#' @importFrom rstudioapi getActiveProject
#' @export
#'
#' @examples getProject()
getProject <- function() {
  tryCatch(basename(rstudioapi::getActiveProject()), error = function(e) {})
}


# _________________________________________________________________________________________________
# Temp _____________________________ ------
# _________________________________________________________________________________________________

# will it be used?
cellID_to_cellType_v1 <- function(cellIDs, ident, obj = aaa) {
  celltypes <- as.named.vector.df(obj@meta.data[, ident], verbose = FALSE)
  celltypes[cellIDs]
}

cellID_to_cellType <- function(cellIDs, ident_w_names) {
  ident_w_names[cellIDs]
}

