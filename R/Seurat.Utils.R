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
parallel.computing.by.future <- function(cores = 4, maxMemSize = 4000 * 1024^2) {
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
#' @param verbose verbose
#' @return A vector of gene names that are found both in the input 'genes' vector and the
#'         Seurat object.
#'
#' @export
IntersectGeneLsWithObject <- function(genes, obj = combined.obj, n_genes_shown = 10,
                                      species_ = "human", EnforceUnique = TRUE, ShowStats = TRUE,
                                      strict = TRUE, verbose = TRUE) {
  message(">>> Running IntersectGeneLsWithObject()")
  # "formerly IntersectWithExpressed(), which still exist in gruffi."

  stopifnot(
    is.character(genes),
    is(obj, "Seurat"),
    is.numeric(n_genes_shown) && n_genes_shown > 0,
    is.logical(strict)
  )
  stopifnot(length(genes) > 0, length(rownames(obj)) > 0)

  # Strict mode: Ensure all genes are present in the Seurat object
  all.genes.found <- all(genes %in% rownames(obj))
  if (!all.genes.found) {
    symbols.missing <- setdiff(genes, rownames(obj))
    iprint(length(symbols.missing), "symbols.missing:", symbols.missing)
    message("running HGNChelper::checkGeneSymbols() to update symbols")

    HGNC.updated <- HGNChelper::checkGeneSymbols(genes, unmapped.as.na = FALSE, map = NULL, species = species_)
    if (ShowStats) {
      HGNC.updated
      print(GetUpdateStats(HGNC.updated))
    }

    if (EnforceUnique) HGNC.updated <- HGNC.EnforceUnique(HGNC.updated)
    genes <- HGNC.updated$Suggested.Symbol

    # UpdateSymbolList(symbols.missing) # Does not catch CTIP2 !!!
    if (strict) stopifnot(all(genes %in% rownames(obj)))
  }

  # Finding genes that are missing in the Seurat object
  missing_in_obj <- setdiff(genes, rownames(obj))
  if (verbose) {
    Stringendo::iprint(
      length(missing_in_obj), " (of ", length(genes),
      ") genes are MISSING from the Seurat object with (", length(rownames(obj)),
      ") genes. E.g.:", head(missing_in_obj, n_genes_shown)
    )
  }

  # Finding genes that are found in both the input list and the Seurat object
  g_found <- intersect(genes, rownames(obj))

  # Output argument assertion
  stopifnot(length(g_found) > 0)

  return(g_found)
}

# _________________________________________________________________________________________________

#' @title Intersect Genes with the List of Noticeably Expressed Genes
#'
#' @description Intersects a vector of gene names with a Seurat object to find genes that are both
#' in the input list and have expression levels in the top quantiles as defined by the object's
#' q99 expression data. It aims to filter genes based on their expression levels being above a
#' specified threshold. Additionally, it offers an option to sort the genes by their expression
#' levels in decreasing order.
#'
#' @param genes A vector of gene names to be intersected with the Seurat object.
#' @param obj A Seurat object containing gene expression data. Default: `combined.obj`.
#' @param above The expression level threshold above which genes are considered noticeably
#' expressed. Default: 0.
#' @param sort A logical flag indicating whether to sort the filtered genes by their expression
#' levels in decreasing order. Default: FALSE.
#' @return A vector of gene names that are found both in the input 'genes' vector and the Seurat
#' object, and have expression levels above the specified 'above' threshold. If `sort` is TRUE,
#' these genes are returned in decreasing order of their expression levels.
#'
#' @examples
#' # Assuming `genes` is a vector of gene names and `
#'
#' @export
SelectHighlyExpressedGenesq99 <- function(genes, obj = combined.obj,
                                          above = 0, sort = FALSE, strict = FALSE) {
  message("Running SelectHighlyExpressedGenesq99()...")
  stopifnot(is.character(genes), is(obj, "Seurat"), is.numeric(above))

  genes.expr <- IntersectGeneLsWithObject(genes = genes, obj = obj, verbose = FALSE, strict = strict)
  if (length(genes.expr) < length(genes)) message("Some genes not expressed. Recommend to IntersectGeneLsWithObject() first.")

  q99.expression <- obj@misc$expr.q99
  print(pc_TRUE(q99.expression == 0, suffix = "of genes at q99.expression are zero"))
  genes.expr.high <- q99.expression[genes.expr]
  if (sort) genes.expr.high <- sort.decreasing(genes.expr.high)
  print(genes.expr.high)
  genes.filt <- names(genes.expr.high)[genes.expr.high > above]

  SFX <- kppws("of the genes are above min. q99 expression of:", above)
  print(pc_TRUE(genes.expr %in% genes.filt, suffix = SFX))

  return(genes.filt)
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
SmallestNonAboveX <- function(vec, X = 0) {
  newmin <- min(vec[vec > X])
  vec[vec <= X] <- newmin
  vec
}


# _________________________________________________________________________________________________
#' @title AreTheseCellNamesTheSame
#'
#' @description Assert and compare two character vectors (e.g.: cell IDs) how much they overlap and
#' plot a Venn Diagram. The function aborts with an error if overlap is too small.
#' @param vec1 Character vector, eg. with cell names
#' @param vec2 Character vector, eg. with cell names
#' @param names Names for plotting
#' @param min.overlap Threshold below there is no there is no meaningful overlap between the two vectors.
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
#' @title Add to Misc or Tools Slot
#'
#' @description This function creates and adds a sub-slot to either the 'misc' or 'tools' slot of a
#' Seurat object. If the sub-slot already exists, it can either be overwritten or a warning will be issued.
#'
#' @param obj A Seurat object.
#' @param pocket_name Which main pocket to use: 'misc' or 'tools'. Default: 'misc'.
#' @param slot_value The value to be assigned to the sub-slot.
#' @param slot_name The name of the sub-slot. Automatically derived from 'sub_slot_value' if not provided.
#' @param sub_slot_value The value to be assigned to the sub-slot.
#' @param sub_slot_name The name of the sub-slot. Automatically derived from 'sub_slot_value' if not provided.
#' @param overwrite A boolean indicating whether to overwrite an existing sub-slot with the same name.
#'
#' @return The modified Seurat object with the new or updated sub-slot.
#'
#' @export
addToMiscOrToolsSlot <- function(obj, pocket_name = "misc",
                                 slot_value = NULL,
                                 slot_name = deparse(substitute(slot_value)),
                                 sub_slot_value = NULL,
                                 sub_slot_name = deparse(substitute(sub_slot_value)),
                                 overwrite = FALSE) {
  stopifnot(is(obj, "Seurat"),
    pocket_name %in% c("misc", "tools"),
    is.character(slot_name), length(slot_name) == 1,
    is.character(sub_slot_name), length(sub_slot_name) == 1,
    "slot name or value is provided" = is.null(slot_value) || !is.null(slot_name),
    "sub_slot name or value is provided" = is.null(sub_slot_value) || !is.null(sub_slot_name)
  )

  # Accessing the specified slot
  pocket <- slot(object = obj, name = pocket_name)

  # Creating new sub_slot or reporting if it exists
  if (slot_name %in% names(pocket) && !overwrite) {
    warning(paste(slot_name, "in", pocket_name, "already exists. Not overwritten."), immediate. = TRUE)
  } else {
    pocket[[slot_name]] <- slot_value
  }

  # Creating new sub_sub_slot or reporting if it exists
  if (sub_slot_name %in% names(pocket[[slot_name]]) && !overwrite) {
    warning(paste(sub_slot_name, "in", pocket_name, "@", slot_name, "already exists. Not overwritten."), immediate. = TRUE)
  } else {
    pocket[[slot_name]][[sub_slot_name]] <- sub_slot_value
  }


  # Assigning the modified slot back to the object
  slot(object = obj, name = pocket_name) <- pocket

  return(obj)
}

# _________________________________________________________________________________________________
#' @title Display Slots in the @tools of an Seurat Object
#'
#' @description
#' `showToolsSlots` prints the names of slots in the `@tools` of a given object.
#' It specifically targets list elements, skipping over data frames and other non-list objects.
#'
#' @param obj An object whose `@tools` slot needs to be examined.
#'
#' @details
#' The function iterates over the slots in the `@tools` of `obj`. If a slot is a list
#' (and not a data frame), it prints the names of elements within this list. If the slot
#' is not a list or is a data frame, it skips printing the names. The function currently
#' does not use the `indent` parameter but it could be incorporated in future enhancements
#' to control the formatting of the output.
#'
#' @examples showToolsSlots(obj)
#'
#' @export
showToolsSlots <- function(obj, max.level = 1, subslot = NULL, ...) {
  slotX <- if (is.null(subslot)) obj@tools else obj@tools[[subslot]]
  str(slotX, max.level = max.level, ...)

  # tools_slot <- names(obj@tools)
  # # i=4
  # for (i in seq(tools_slot)) {
  #   cat("", fill = TRUE)
  #   message("obj@tools$", tools_slot[i])
  #
  #   x <- obj@tools[[tools_slot[i]]]
  #   if (!is.data.frame(x) & is.list(x)) {
  #     print(paste("  ", names(x)), width = 12)
  #   } else {
  #     return(idim(x))
  #   }
  # }
}


# _________________________________________________________________________________________________
#' @title Display Slots in the @misc of an Seurat Object
#'
#' @description See `showToolsSlots` for details. Prints the names of slots in the `@misc` of a given object.
#' It specifically targets list elements, skipping over data frames and other non-list objects.
#'
#' @param obj An object whose `@misc` slot needs to be examined.
#' @param max.level Max depth to dive into sub-elements.
#' @param subslot A subslot within `@misc`.
#' @param ... ...
#'
#' @examples showToolsSlots(obj)
#'
#' @export
showMiscSlots <- function(obj, max.level = 1, subslot = NULL,
                          ...) {
  slotX <- if (is.null(subslot)) obj@misc else obj@misc[[subslot]]
  str(slotX, max.level = max.level, ...)

  # Path to slot
  msg <- paste0(substitute(obj), "@misc")
  if (!is.null(subslot)) msg <- paste0(msg, "$", substitute(subslot))
  message(msg)
}




# _________________________________________________________________________________________________
#' @title calc.q99.Expression.and.set.all.genes

#' @description Calculate the gene expression of the e.g.: 99th quantile (expression in the top 1% cells).
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
    obj = combined.obj,
    quantileX = 0.99, max.cells = 1e5,
    slot = "data",
    assay = c("RNA", "integrated")[1],
    set.misc = TRUE,
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

  iprint(
    "Quantile", quantileX, "is now stored under obj@misc$all.genes and $", slot_name,
    " Please execute all.genes <- obj@misc$all.genes."
  )
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
    suffix.plot <- if (nchar(suffix.plot)) make.names(suffix.plot)
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
#' @title Retrieve Cluster Names
#'
#' @description Extracts cluster names based on a specified identity class from a Seurat object.
#'
#' @param obj A Seurat object. Default: `combined.obj`.
#' @param ident The identity class from which to retrieve cluster names.
#' Default uses the second clustering run from `GetClusteringRuns(obj)`.
#' @examples
#' \dontrun{
#'  getClusterNames(obj = combined.obj, ident = GetClusteringRuns(obj)[2])
#'  }
#' @return Prints and returns the sorted unique cluster names as a character vector.
#' @export
getClusterNames <- function(obj = combined.obj, ident = GetClusteringRuns(obj)[2]) {
  print(GetClusteringRuns(obj))
  clz <- as.character(sort(deframe(unique(obj[[ident]]))))
  cat(dput(clz))
}



# _________________________________________________________________________________________________
#' @title GetClusteringRuns
#'
#' @description Get Clustering Runs: metadata column names.
#' @param obj Seurat object, Default: combined.obj
#' @param res Clustering resoluton to use, Default: FALSE
#' @param pat Pettern to match, Default: `*snn_res.*[0-9]$`
#' @examples
#' \dontrun{
#' if (interactive()) {
  #'   GetClusteringRuns(obj = combined.obj, pat = '*snn_res.*[0-9]$')
#' }
#' }
#' @export
GetClusteringRuns <- function(obj = combined.obj, res = FALSE, pat = '*snn_res.*[0-9]$') {
  if (res) pat <- gsub(x = pat, pattern = "\\[.*\\]", replacement = res)
  clustering.results <- CodeAndRoll2::grepv(x = colnames(obj@meta.data), pattern = pat)
  if (identical(clustering.results, character(0))) warning("No matching column found!", immediate. = TRUE)
  dput(clustering.results)
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
    obj = combined.obj,
    res = list(FALSE, 0.5)[[1]], topgene = FALSE,
    pat = c("^cl.names.Known.*[0,1]\\.[0-9]$", "Name|name")[2]) {
  if (res) pat <- gsub(x = pat, pattern = "\\[.*\\]", replacement = res)
  if (topgene) pat <- gsub(x = pat, pattern = "Known", replacement = "top")
  clustering.results <- CodeAndRoll2::grepv(x = colnames(obj@meta.data), pattern = pat)
  if (identical(clustering.results, character(0))) {
    warning("No matching column found! Trying GetClusteringRuns(..., pat = '*_res.*[0,1]\\.[0-9]$)",
      immediate. = TRUE
    )
    clustering.results <- GetClusteringRuns(obj = obj, res = FALSE, pat = "*_res.*[0,1]\\.[0-9]$")
  }
  dput(clustering.results)
  return(clustering.results)
}



# _________________________________________________________________________________________________
#' @title GetOrderedClusteringRuns
#'
#' @description Get Clustering Runs: metadata column names.
#' @param obj Seurat object, Default: combined.obj.
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
GetOrderedClusteringRuns <- function(obj = combined.obj, res = FALSE, pat = "*snn_res.*[0,1]\\.[0-9]\\.ordered$") {
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
#' @description Plot gene expression based on the expression at the 90th quantile
#' (so you will not lose genes expressed in few cells).
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
    title <- paste(
      gene, "is in the", Stringendo::percentage_formatter(quantile.GOI),
      "quantile of 'q90-av' expression. \n There are", counts, "counts"
    )
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
# 3D umaps ______________________________ ----
# _________________________________________________________________________________________________



# _________________________________________________________________________________________________
#' @title BackupReduction
#'
#' @description Backup UMAP to `obj@misc$reductions.backup` from `obj@reductions$umap`. #
#' @param obj Seurat object, Default: combined.obj
#' @param dim Numer of dimensions used, Default: 2
#' @param reduction UMAP, tSNE, or PCA (Dim. reduction to use), Default: 'umap'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   obj <- BackupReduction(obj = obj, dim = 2, reduction = "umap")
#' }
#' }
#' @export
BackupReduction <- function(obj = combined.obj, dim = 2, reduction = "umap") {
  if (is.null(obj@misc$"reductions.backup")) obj@misc$"reductions.backup" <- list()
  dslot <- paste0(reduction, dim, "d")
  obj@misc$reductions.backup[[dslot]] <- obj@reductions[[reduction]]
  return(obj)
}


# _________________________________________________________________________________________________
#' @title SetupReductionsNtoKdimensions
#'
#' @description Function to calculate N-to-K dimensional umaps (default = 2:3); and back them up to
#' slots `obj@misc$reductions.backup` from @reductions$umap
#' @param obj A Seurat object. Default: combined.obj
#' @param nPCs A numeric value representing the number of principal components to use. Default: p$n.PC
#' @param dimensions A numeric vector specifying the dimensions to use for the dimensionality reductions. Default: 3:2
#' @param reduction A character string specifying the type of dimensionality reduction to perform.
#' Can be "umap", "tsne", or "pca". Default: 'umap'
#' @return The input Seurat object with computed dimensionality reductions and backups of these reductions.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   combined.obj <- SetupReductionsNtoKdimensions(obj = combined.obj, nPCs = 10, dimensions = 2:3, reduction = "umap")
#' }
#' }
#' @export
SetupReductionsNtoKdimensions <- function(obj = combined.obj, nPCs = p$"n.PC", dimensions = 3:2,
                                          reduction = "umap", ...) {
  red <- reduction
  for (d in dimensions) {
    iprint(d, "dimensional", red, "is calculated")
    obj <- if (reduction == "umap") {
      RunUMAP(obj, dims = 1:nPCs, n.components = d, ...)
    } else if (reduction == "tsne") {
      RunTSNE(obj, dims = 1:nPCs, n.components = d, ...)
    } else if (reduction == "pca") {
      RunPCA(obj, dims = 1:nPCs, n.components = d, ...)
    }
    obj <- BackupReduction(obj = obj, dim = d, reduction = red)
  }
  return(obj)
}


# _________________________________________________________________________________________________
#' @title RecallReduction
#'
#' @description Set active UMAP to `obj@reductions$umap` from `obj@misc$reductions.backup`.
#' @param obj Seurat object, Default: combined.obj
#' @param dim Numer of dimensions used, Default: 2
#' @param reduction UMAP, tSNE, or PCA (Dim. reduction to use), Default: 'umap'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   combined.obj <- RecallReduction(obj = combined.obj, dim = 2, reduction = "umap")
#'   qUMAP()
#'   combined.obj <- RecallReduction(obj = combined.obj, dim = 3, reduction = "umap")
#'   qUMAP()
#' }
#' }
#' @export
RecallReduction <- function(obj = combined.obj, dim = 2, reduction = "umap") {
  dslot <- paste0(reduction, dim, "d")
  reduction.backup <- obj@misc$reductions.backup[[dslot]]
  msg <- paste(dim, "dimensional", reduction, "from obj@misc$reductions.backup")
  stopif(is.null(reduction.backup), message = paste0(msg, " is NOT FOUND"))
  iprint(msg, "is set active. ")
  stopifnot(dim == ncol(reduction.backup))
  obj@reductions[[reduction]] <- reduction.backup
  return(obj)
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
#' @title Recall all.genes global variable from a Seurat object
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
recall.all.genes <- function(obj = combined.obj, overwrite = FALSE) {
  obj <- ww.get.1st.Seur.element(obj)

  if ("all.genes" %in% names(obj@misc)) {

    if (!exists("all.genes") | overwrite) {
      all.genes <- obj@misc$all.genes
      print(head(unlist(all.genes)))
      MarkdownHelpers::ww.assign_to_global(name = "all.genes", value = all.genes, verbose = FALSE)
      message("all.genes is now (re)defined in the global environment.")
    } else {
      message("  ->   Variable 'all.genes' exits in the global namespace, and overwrite is: FALSE")
    }
  } else {
    message("  ->   Slot 'all.genes' does not exist in obj@misc.")
    hits <- grepv(pattern = "expr.", names(obj@misc))
    if (!is.null(hits)) {
      message("Found instead (", hits, "). Returning 1st element: ", hits[1])
      all.genes <- obj@misc[[hits[1]]]
      MarkdownHelpers::ww.assign_to_global(name = "all.genes", value = as.list(all.genes), verbose = FALSE)
    }
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
    p_found <- exists("p", envir = .GlobalEnv)
    if (p_found) message("  ->   Variable 'p' exits in the global namespace.")

    if (!p_found | (p_found & overwrite == TRUE)) {
      MarkdownHelpers::ww.assign_to_global(name = "p", value = obj@misc$"p", verbose = F)
      message("p is now (re)defined in the global environment.")
    } else {
      message("p not overwritten.")
    }
  } else {
    message("  ->   Slot 'p' does not exist in obj@misc.")
  }
}



# _________________________________________________________________________________________________
#' @title recall.genes.ls
#'
#' @description Recall genes.ls from obj@misc to "genes.ls" in the global environment.
#' @param obj Seurat object, Default: combined.obj
#' @param overwrite Overwrite already existing in environment? Default: FALSE
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
    if (!exists("genes.ls")) message("variable 'genes.ls' exits in the global namespace: ", head(p))

    if (!exists("genes.ls") | (exists("genes.ls") & overwrite == TRUE)) {
      MarkdownHelpers::ww.assign_to_global(name = "genes.ls", value = obj@misc$"genes.ls")
      message("Overwritten.")
    } else {
      message("Not overwritten.")
    }
  } else {
    message("  ->   Slot 'genes.ls' does not exist in obj@misc.")
  }
}


# _________________________________________________________________________________________________
#' @title Save Parameters to Seurat Object
#'
#' @description Stores a list of parameters within the `@misc$p` slot of a Seurat object,
#' allowing for easy reference and tracking of analysis parameters used.
#'
#' @param obj Seurat object to update; Default: `combined.obj`.
#' @param params List of parameters to save; Default: `p`.
#' @param overwrite Logical indicating if existing parameters should be overwritten; Default: TRUE.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   save.parameters(obj = combined.obj, params = p)
#' }
#' }
#'
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
# Merging objects and @misc ______________________________ ----
# _________________________________________________________________________________________________


#' @title Copy Specified Elements from One Seurat Object's @misc to Another's
#'
#' @description Copies specified elements from the `@misc` slot of one Seurat object to the `@misc` slot
#' of another. It warns if some specified elements are missing in the source object or if elements are
#' overwritten in the destination object, depending on the `overwrite` argument.
#'
#' @param obj.from The source Seurat object from which elements in the `@misc` slot are to be copied.
#' @param obj.to The destination Seurat object to which elements in the `@misc` slot are to be copied.
#' @param elements.needed A vector of strings specifying the names of the elements in the `@misc` slot of
#' `obj.from` that should be copied to `obj.to`.
#' @param overwrite Logical indicating whether to overwrite elements in `obj.to` that already exist.
#' If `TRUE`, existing elements will be overwritten with a warning; if `FALSE`, the function will
#' stop with an error if it tries to copy an element that already exists in `obj.to`.
#' @return Returns the modified destination Seurat object (`obj.to`) with the specified elements
#' added to or updated in its `@misc` slot.
#' @examples
#' # Assuming `obj1` and `obj2` are Seurat objects and you wish to copy specific elements
#' # from obj1 to obj2, possibly overwriting existing elements in obj2
#' obj2 <- copyMiscElements(obj1, obj2, c("element1", "element2"), overwrite = TRUE)
#'
#' @export
copyMiscElements <- function(obj.from, obj.to, elements.needed, overwrite = TRUE) {
  obj.from <- ww.get.1st.Seur.element(obj.from)

  stopifnot(
    inherits(obj.from, "Seurat"),
    inherits(obj.to, "Seurat")
  )

  # Check for missing elements in obj.to@misc
  elements.from <- names(obj.from@misc)
  missing <- setdiff(elements.needed, elements.from)
  if (length(missing) > 0) {
    warning("Missing elements in obj.from@misc: ", paste(missing, collapse = ", "), immediate. = TRUE)
  }

  # Check for existing elements in obj.to@misc
  elements.already.exisiting <- intersect(elements.needed, names(obj.to@misc))
  if (length(elements.already.exisiting) > 0) {
    if (!overwrite) {
      stop(
        "The following elements already exist in obj.to@misc and 'overwrite' is FALSE: ",
        paste(elements.already.exisiting, collapse = ", ")
      )
    } else {
      warning("Overwriting the following elements in obj.to@misc: ",
        paste(elements.already.exisiting, collapse = ", "),
        immediate. = TRUE
      )
    }
  }

  # Copy specified elements from obj.from to obj.to
  existingElementsFrom <- intersect(elements.needed, names(obj.from@misc))
  for (element in existingElementsFrom) {
    obj.to@misc[[element]] <- obj.from@misc[[element]]
  }
  iprint("@misc contains: ", names(obj.to@misc))

  return(obj.to)
}


# _________________________________________________________________________________________________
#' @title Copy Tools Slots from Multiple Seurat Objects
#'
#' @description This function copies the `@tools` slots from a list of Seurat objects into a new slot
#' of a target Seurat object. This allows for the aggregation of tools information from multiple
#' experiments or datasets into a single, consolidated Seurat object.
#'
#' @param ls.obj A list of Seurat objects from which the `@tools` slots will be copied.
#' @param obj.to The target Seurat object to which the `@tools` slots will be added.
#' @param overwrite A logical parameter that is kept for compatibility but not used in this version.
#' Its presence does not affect the function's behavior.
#' @param new.slot The name of the new slot within `obj.to@tools` where the copied `@tools` information
#' will be stored. This allows for the organization of copied tools under a specific label, facilitating
#' easy access and interpretation.
#' @return Returns the modified target Seurat object (`obj.to`) with a new `@tools` slot containing
#' the copied information from the list of Seurat objects.
#' @examples
#' # Assuming `ls.obj` is a list of Seurat objects and `obj.to` is a target Seurat object
#' obj.to <- copyCompleteToolsSlots(ls.obj, obj.to, overwrite = TRUE, new.slot = "per.experiment")
#' @export
copyCompleteToolsSlots <- function(ls.obj, obj.to, overwrite = TRUE, new.slot = "per.experiment") {
  stopifnot(
    inherits(obj.to, "Seurat"),
    all(sapply(ls.obj, inherits, "Seurat"))
  )

  ls.tools <- lapply(ls.obj, function(x) x@tools)
  obj.to@tools[[new.slot]] <- ls.tools

  return(obj.to)
}




# _________________________________________________________________________________________________
# Subsetting the Seurat object ______________________________ ----
# _________________________________________________________________________________________________


#' @title Subset a Seurat Object by Identity
#'
#' @description Subsets a Seurat object based on a specified identity column and values. It allows
#'   for an optional inversion of the selection.
#'
#' @param obj A Seurat object. Default: `NULL`.
#' @param ident The name of the identity column to use for subsetting. It is recommended to
#'   specify this explicitly. Default: First entry from the result of `GetClusteringRuns()`.
#' @param clusters A vector of cluster values for which cells should be matched and retained.
#'   This parameter does not have a default value and must be specified.
#' @param invert A logical indicating whether to invert the selection, keeping cells that do
#'   not match the specified clusters. Default: `FALSE`.
#'
#' @return A Seurat object subsetted based on the specified identity and clusters.
#'
#' @examples
#' # Assuming `seurat_obj` is your Seurat object and you want to subset based on cluster 1
#' subsetted_obj <- subsetSeuObjByIdent(
#'   obj = seurat_obj, ident = "your_ident_column",
#'   clusters = c(1), invert = FALSE
#' )
#'
#' @export
subsetSeuObjByIdent <- function(
    obj = combined.obj, ident = GetClusteringRuns()[1],
    clusters,
    invert = FALSE) {

  # Input checks
  stopifnot(
    "obj must be a Seurat object" = inherits(obj, "Seurat"),
    "ident must be a character and exist in obj@meta.data" = is.character(ident) && ident %in% colnames(obj@meta.data),
    "clusters must exist in ident" = all(clusters %in% unique(combined.obj[[ident]][ ,1] ))
  )

  Idents(obj) <- ident
  cellz <- WhichCells(obj, idents = clusters, invert = invert)
  message(length(cellz), " cells are selected from ", ncol(obj),
          ", using values: ", clusters,
          ", from ", ident, ".")
  subset(x = obj, cells = cellz)
}


# _________________________________________________________________________________________________
#' @title downsampleSeuObj
#'
#' @description Subset a compressed Seurat object and save it in the working directory.
#' @param obj A Seurat object to subset. Default: the i-th element of the list 'ls.Seurat'.
#' @param fractionCells The fraction of the object's data to keep. Default: 0.25.
#' @param nCells If set to a number greater than 1, indicates the absolute number of cells to keep.
#' If FALSE, the function uses 'fractionCells' to determine the number of cells. Default: FALSE.
#' @param seed A seed for random number generation to ensure reproducible results. Default: 1989.
#' @export
#' @importFrom Stringendo percentage_formatter

downsampleSeuObj <- function(obj = ls.Seurat[[i]], fractionCells = 0.25, nCells = FALSE,
                             seed = 1989) {
  set.seed(seed)
  if (isFALSE(nCells)) {
    cellIDs.keep <- sampleNpc(metaDF = obj@meta.data, pc = fractionCells)
    iprint(
      length(cellIDs.keep), "or", Stringendo::percentage_formatter(fractionCells),
      "of the cells are kept. Seed:", seed
    )
  } else if (nCells > 1) {
    nKeep <- min(ncol(obj), nCells)
    # print(nKeep)
    cellIDs.keep <- sample(colnames(obj), size = nKeep, replace = FALSE)
    if (nKeep < nCells) {
      iprint(
        "Only", nCells,
        "cells were found in the object, so downsampling is not possible."
      )
    }
  }
  obj <- subset(x = obj, cells = cellIDs.keep) # downsample
  return(obj)
}

# _________________________________________________________________________________________________
#' @title downsampleSeuObj.and.Save
#'
#' @description Subset a compressed Seurat Obj and save it in wd. #
#' @param obj Seurat object, Default: ORC
#' @param fraction Fractional size to downsample to. Default: 0.25
#' @param seed random seed used, Default: 1989
#' @param min.features Minimum features
#' @param dir Directory to save to. Default: OutDir
#' @param suffix A suffix added to the filename, Default: ''
#' @export
downsampleSeuObj.and.Save <- function(
    obj = ORC, fraction = 0.25, seed = 1989, dir = OutDir,
    min.features = p$"min.features", suffix = fraction,
    nthreads = .getNrCores()
    ) {
  obj_Xpc <- downsampleSeuObj(obj = obj, fractionCells = fraction, seed = seed)
  nr.cells.kept <- ncol(obj_Xpc)

  # Seurat.utils:::.saveRDS.compress.in.BG(obj = obj_Xpc, fname = ppp(paste0(dir, substitute(obj)),
  # suffix, nr.cells.kept, 'cells.with.min.features', min.features,"Rds" ) )
  xsave(obj_Xpc,
    suffix = ppp(suffix, nr.cells.kept, "cells.with.min.features", min.features),
    nthreads = nthreads, project = getProject(), showMemObject = TRUE, saveParams = FALSE
  )
}



# _________________________________________________________________________________________________
#' @title Sample Cells From Identifiers in Seurat Object
#'
#' @description This function samples a specified maximum number of cells from each identity class
#' in a Seurat object, in the meta.data. It ensures that the sampling does not exceed the total
#' number of cells available per identity.
#'
#' @param obj A Seurat object from which cells are to be sampled.
#' @param ident A character vector specifying the identity class from which cells are to be sampled.
#' @param max.cells A positive integer indicating the maximum number of cells to sample from each identity class.
#' @param verbose Logical indicating if messages about the sampling process should be printed to the console. Defaults to TRUE.
#'
#' @return Returns a Seurat object containing only the sampled cells.
#'
#' @details This function checks for the presence of the specified identity class within the object's metadata.
#' If the number of cells within any identity class is less than or equal to the `max.cells` parameter,
#' all cells from that class are retained. Otherwise, a random sample of `max.cells` is taken from the class.
#' The function updates the identity of the cells in the returned Seurat object to reflect the sampled cells.
#' If `verbose` is TRUE, it prints the total number of cells sampled and provides a visual summary of the fraction
#' of cells retained per identity class.
#'
#' @examples
#' # Assuming `seuratObj` is a Seurat object with identities stored in its metadata
#' sampledSeuratObj <- downsampleSeuObjByIdentAndMaxcells(obj = seuratObj, ident = "cellType", max.cells = 100)
#'
#' @importFrom CodeAndRoll2 as.named.vector.df
#'
#' @export
#'
downsampleSeuObjByIdentAndMaxcells <- function(obj,
                                               ident = GetNamedClusteringRuns()[1],
                                               max.cells = min(table(combined.obj[[ident]])),
                                               verbose = TRUE,
                                               seed = 1989) {
  stopifnot(
    "obj must be a Seurat object" = inherits(obj, "Seurat"),
    "ident must be a character and exist in obj@meta.data" = is.character(ident) && ident %in% colnames(obj@meta.data),
    "max.cells must be a positive integer" = is.numeric(max.cells) && max.cells > 0,
    max.cells < ncol(obj)
  )

  data <- CodeAndRoll2::as.named.vector.df(obj[[ident]])
  uniqueCategories <- unique(data)

  set.seed(seed)
  sampledNames <- lapply(uniqueCategories, function(category) {
    namesInCategory <- names(data[data == category])
    if (length(namesInCategory) <= max.cells) {
      return(namesInCategory)
    } else {
      return(sample(namesInCategory, max.cells))
    }
  })

  sampledCells <- unlist(sampledNames)

  Idents(obj) <- ident
  obj2 <- subset(x = obj, cells = sampledCells)

  subb <- paste0("From ", ncol(obj), " reduced to ", ncol(obj2), " cells.")
  message(subb)

  if (verbose) {
    cat("Total cells sampled:", length(sampledCells), "\n")
    nr_remaining_cells <- orig_cells <- table(data)
    nr_remaining_cells[nr_remaining_cells > max.cells] <- max.cells
    fr_remaining_per_cluster <- iround(nr_remaining_cells / orig_cells)
    print(fr_remaining_per_cluster)
    pobj <- qbarplot(
      vec = fr_remaining_per_cluster, subtitle = subb, label = fr_remaining_per_cluster,
      ylab = "fr. of cells", save = FALSE
    )
    print(pobj)
  }
  return(obj2)
}


# _________________________________________________________________________________________________
#' @title Remove Residual Small Clusters from Seurat Object
#'
#' @description Removes clusters containing fewer cells than specified by `max.cells`
#' from a Seurat object. This function is particularly useful after subsetting a dataset,
#' where small, possibly unrepresentative clusters may remain.
#'
#' @param obj Seurat object from which small clusters will be removed; Default: `combined.obj`.
#' @param identitites Vector of clustering identities to examine for small clusters;
#' Default: `GetClusteringRuns(obj)`.
#' @param max.cells Maximum number of cells a cluster can contain to still be considered for removal.
#' Default: The lesser of 0.5% of the dataset or 5 cells.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   combined.obj <- removeResidualSmallClusters(obj = combined.obj)
#' }
#' }
#'
#' @export
removeResidualSmallClusters <- function(
    obj = combined.obj,
    identitites = GetClusteringRuns(obj),
    max.cells = max(round((ncol(obj)) / 2000), 5)) {
  META <- obj@meta.data
  all.cells <- rownames(META)

  iprint("max.cells:", max.cells, "| Scanning over these", length(identitites), "identities:", identitites)
  small.clusters <- cells.to.remove <- list.fromNames(identitites)

  for (i in 1:length(identitites)) {
    # colX <- identitites[i]
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
    iprint(
      ">>> a total of", length(all.cells.2.remove),
      "cells are removed which belonged to a small cluster in any of the identities."
    )
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
#' @description This function removes residual small clusters from specified Seurat objects and
#' drops levels in factor-like metadata.
#' @param ls_obj A list of Seurat objects.
#' @param object_names A character vector containing the names of the Seurat objects to process.
#' Default is names of all objects in the `ls_obj`.
#' @param indices A numeric vector indicating which datasets to process by their position in
#' the `object_names` vector. By default, it processes the second and third datasets.
#' @param ... Additional parameters passed to the `removeResidualSmallClusters` function.
#'
#' @details This function applies `removeResidualSmallClusters` and `dropLevelsSeurat` to
#' the Seurat objects specified by the `indices` in the `object_names`.
#' It operates in place, modifying the input `ls_obj` list.
#'
#' @return The function returns the modified list of Seurat objects.
#' @examples
#' \dontrun{
#' # Process the 2nd and 3rd datasets
#' removeClustersAndDropLevels(ls_obj, indices = c(2, 3))
#' }
#'
#' @export
removeClustersAndDropLevels <- function(
    ls_obj, object_names = names(ls_obj),
    indices = 2:3, ...) {
  for (index in indices) {
    dataset_name <- object_names[index]
    obj <- ls_obj[[dataset_name]]
    obj <- removeResidualSmallClusters(obj = obj, identitites = GetClusteringRuns(obj), ...)
    obj <- dropLevelsSeurat(obj)
    ls_obj[[dataset_name]] <- obj
  }
  return(ls_obj)
}




# _________________________________________________________________________________________________
#' @title Remove Cells by Dimension Reduction
#'
#' @description This function applies a cutoff in the specified dimension of a given
#' dimension reduction (UMAP, PCA, or t-SNE) to remove cells.
#' @param reduction A string specifying the dimension reduction technique to be used
#' ('umap', 'pca', or 'tsne'). Default is 'umap'.
#' @param umap_dim An integer specifying which dimension (axis) to apply the cutoff. Default is 1.
#' @param obj A Seurat object. Default is 'combined.obj'.
#' @param cutoff A numerical value indicating the cutoff value for the specified dimension. Default is 0.
#' @param cut_below A logical value indicating whether to remove cells below (TRUE) or
#' above (FALSE) the cutoff line. Default is TRUE.
#' @param only_plot_cutoff Simulate and plot cutoff only.
#' @param ... Any other parameters to be passed to internally called functions.
#' @return A Seurat object with cells removed according to the specified cutoff.
#' @export
removeCellsByUmap <- function(
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
# Downsampling Lists of Seurat objects ______________________________ ----
# _________________________________________________________________________________________________


#' @title Downsample a List of Seurat Objects to a Specific Number of Cells
#'
#' @description Downsampling each Seurat object in a list to a specified number of cells. This function is
#' particularly useful for creating smaller, more manageable subsets of large single-cell datasets for
#' preliminary analyses or testing.
#'
#' @param ls.obj List of Seurat objects to be downsampled; Default: `ls.Seurat`.
#' @param NrCells Target number of cells to downsample each Seurat object to.
#' @param save_object Logical indicating whether to save the downsampled Seurat objects using `isaveRDS`
#' or to return them; Default: FALSE.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   downsampledSeuratList <- downsampleListSeuObjsNCells(
#'     ls.obj =
#'       list(yourSeuratObj1, yourSeuratObj2), NrCells = 2000
#'   )
#'   downsampledSeuratList <- downsampleListSeuObjsNCells(NrCells = 200)
#' }
#' }
#'
#' @export
#' @importFrom tictoc tic toc
#' @importFrom Stringendo percentage_formatter
#' @importFrom foreach foreach %dopar% getDoParRegistered

downsampleListSeuObjsNCells <- function(
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
      downsampleSeuObj(obj = ls.obj[[i]], nCells = NrCells)
    }
    names(ls.obj.downsampled) <- names.ls
  } else {
    ls.obj.downsampled <- list.fromNames(names.ls)
    for (i in 1:n.datasets) {
      iprint(names(ls.obj)[i], Stringendo::percentage_formatter(i / n.datasets, digitz = 2))
      ls.obj.downsampled[[i]] <- downsampleSeuObj(obj = ls.obj[[i]], nCells = NrCells)
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
#' @title Downsample a List of Seurat Objects to a Fraction
#'
#' @description Downsampling a list of Seurat objects to a specified fraction of their original size.
#' This is useful for reducing dataset size for quicker processing or testing workflows.
#'
#' @param ls.obj List of Seurat objects to be downsampled; Default: `ls.Seurat`.
#' @param fraction Fraction of cells to retain in each Seurat object; Default: 0.1.
#' @param save_object Logical indicating whether to save the downsampled Seurat objects using
#' `isaveRDS` or return them; Default: FALSE.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   downsampled_objs <- downsampleListSeuObjsPercent(ls.obj = yourListOfSeuratObjects, fraction = 0.1)
#' }
#' }
#'
#' @export
#' @importFrom tictoc tic toc
#' @importFrom Stringendo percentage_formatter
#' @importFrom foreach foreach %dopar% getDoParRegistered
#'
downsampleListSeuObjsPercent <- function(
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
      downsampleSeuObj(obj = ls.obj[[i]], fractionCells = fraction)
    }
    names(ls.obj.downsampled) <- names.ls
  } else {
    ls.obj.downsampled <- list.fromNames(names.ls)
    for (i in 1:n.datasets) {
      cells <- round(ncol(ls.obj[[1]]) * fraction)
      iprint(names(ls.obj)[i], cells, "cells=", Stringendo::percentage_formatter(i / n.datasets, digitz = 2))
      ls.obj.downsampled[[i]] <- downsampleSeuObj(obj = ls.obj[[i]], fractionCells = fraction)
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
# DGEA ______________________________ ----
# _________________________________________________________________________________________________


# _________________________________________________________________________________________________
#' @title Add.DE.combined.score
#'
#' @description Add a combined score to differential expression (DE) results. The score is
#' calculated as log-fold change (LFC) times negative logarithm of scaled
#' p-value (LFC * -log10( p_cutoff / pval_scaling )).
#' @param df A data frame that holds the result of a differential gene expression analysis,
#' typically obtained via the 'FindAllMarkers' function. Default: df.markers.
#' @param p_val_min The minimum p-value considered. All values below this threshold are set to
#' this value. Default: 1e-25.
#' @param pval_scaling The value to scale p-values by in the calculation of the combined score. Default: 0.001.
#' @param colP The name of the column in the input data frame that holds p-values. Default: 'p_val'.
#' @param colLFC The name of the column in the input data frame that holds log-fold change values.
#' By default, it selects the first column not named "avg_logFC" or "avg_log2FC".
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
#' @title Save Top 25 Markers per Cluster
#'
#' @description Stores the top 25 markers for each cluster identified in a Seurat object, based on
#' the `avg_log2FC` from the output table of `FindAllMarkers()`. The result is saved under `@misc$df.markers$res...`,
#' rounding insignificant digits to three decimal places.
#'
#' @param obj Seurat object to update with top 25 markers information; Default: `combined.obj`.
#' @param df_markers Data frame containing results from differential gene expression analysis
#' via `FindAllMarkers()`, specifying significant markers across clusters; Default: `df.markers`.
#' @param res Clustering resolution at which the markers were identified; Default: 0.5.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   combined.obj <- StoreTop25Markers(obj = combined.obj, df_markers = df.markers, res = 0.5)
#' }
#' }
#'
#' @seealso \code{\link[Seurat]{FindAllMarkers}}, \code{\link[dplyr]{top_n}}
#'
#' @export
#' @importFrom Seurat FindAllMarkers
#' @importFrom dplyr group_by top_n select arrange
#' @importFrom magrittr `%>%`

StoreTop25Markers <- function(
    obj = combined.obj,
    df_markers = df.markers, res = 0.5) {
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
#' @title Store All Differential Expression Markers
#'
#' @description Saves the complete output table from `FindAllMarkers()` to a Seurat object, facilitating
#' easy access to differential expression analysis results. This function rounds numerical values to a
#' specified number of digits to maintain readability and manage file sizes.
#'
#' @param obj Seurat object to update with differential expression markers; Default: `combined.obj`.
#' @param df_markers Data frame containing the results from differential gene expression analysis
#' (`FindAllMarkers()` output); Default: `df.markers`.
#' @param res Clustering resolution identifier for storing and referencing the markers; Default: 0.5.
#' @param digit Number of significant digits to retain in numerical values; Default: 3.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   combined.obj <- StoreAllMarkers(obj = combined.obj, df_markers = df.markers, res = 0.5)
#' }
#' }
#'
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
#' @title Get Top Differential Expression Genes Data Frame
#'
#' @description Retrieves a data frame of the top N differentially expressed genes from
#' differential gene expression analysis results, offering an option to exclude certain genes
#' based on patterns.
#'
#' @param dfDE Data frame containing the results of differential gene expression analysis
#' (e.g., output from `FindAllMarkers()`); Default: `df.markers`.
#' @param n Number of top markers to retrieve per cluster; Default: `p$n.markers`.
#' @param order.by Priority column for sorting markers before selection, such as `"avg_log2FC"`;
#' Default: `"avg_log2FC"`.
#' @param exclude Vector of regex patterns to exclude genes from the top markers list;
#' Default: `c("^AL*|^AC*|^LINC*")`.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   topMarkersDF <- GetTopMarkersDF(dfDE = df.markers, n = 3)
#' }
#' }
#'
#' @seealso \code{\link[Seurat]{FindAllMarkers}}, \code{\link[dplyr]{arrange}},
#'  \code{\link[dplyr]{filter}}, \code{\link[dplyr]{group_by}}
#'
#' @export
#' @importFrom dplyr arrange group_by slice select filter
# #' @importFrom magrittr `%>%`
GetTopMarkersDF <- function(
    dfDE = df.markers # Get the vector of N most diff. exp. genes.
    , n = p$"n.markers", order.by = c("avg_log2FC", "p_val_adj")[1],
    exclude = c("^AL*|^AC*|^LINC*")) {
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
#' @title Get Top Differential Expression Markers from DGEA Results
#'
#' @description Retrieves the top N differentially expressed genes from the results of a differential
#' gene expression analysis, such as that provided by `FindAllMarkers()`.
#'
#' @param dfDE Data frame containing differential expression analysis results; Default: `df.markers`.
#' @param n Number of top markers to retrieve for each cluster; Default: `p$n.markers`.
#' @param order.by Column by which to sort the markers before selection, typically prioritizing
#' markers by significance or effect size; Default: `"avg_log2FC"`.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   topMarkers <- GetTopMarkers(df = df.markers, n = 3)
#' }
#' }
#'
#' @seealso \code{\link[Seurat]{FindAllMarkers}}, \code{\link[dplyr]{arrange}}, \code{\link[dplyr]{group_by}}
#'
#' @export
#' @importFrom dplyr arrange group_by slice select
# #' @importFrom magrittr `%>%`
GetTopMarkers <- function(dfDE = df.markers,
                          n = p$"n.markers",
                          order.by = c("combined.score", "avg_log2FC", "p_val_adj")[2]) {
  message("Works on active Idents()") # thus we call cluster
  TopMarkers <- dfDE |>
    arrange(desc(!!as.name(order.by))) |>
    group_by(cluster) |>
    dplyr::slice(1:n) |>
    dplyr::select(gene) |>
    col2named.vec.tbl()

  # TopMarkers <- dfDE %>%
  #   arrange(desc(!!as.name(order.by))) %>%
  #   group_by(cluster) %>%
  #   dplyr::slice(1:n) %>%
  #   dplyr::select(gene) %>%
  #   col2named.vec.tbl()

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
#' if (interactive()) {
#'   combined.obj <- AutoLabelTop.logFC()
#'   combined.obj$"cl.names.top.gene.res.0.5"
#' }
#' }

#' @export
AutoLabelTop.logFC <- function(
    obj = combined.obj,
    group.by = GetClusteringRuns(obj)[1],
    res = 0.1, plot.top.genes = TRUE,
    suffix = res,
    order.by = c("combined.score", "avg_log2FC", "p_val_adj")[2],
    exclude = c("^AL*|^AC*|^LINC*|^C[0-9]orf*"),
    df_markers = obj@misc$"df.markers"[[paste0("res.", res)]],
    plotEnrichment = TRUE) {
  stopifnot(
    !is.null("df_markers"),
    order.by %in% colnames(df_markers)
  )

  df.top.markers <- GetTopMarkersDF(dfDE = df_markers, order.by = order.by, n = 1, exclude = exclude)

  if (plotEnrichment) {
    top_log2FC <- df.top.markers$"avg_log2FC"
    names(top_log2FC) <- ppp(df.top.markers$"cluster", df.top.markers$"gene")
    ggExpress::qbarplot(top_log2FC,
      label = iround(top_log2FC),
      subtitle = suffix,
      ylab = "avg_log2FC", xlab = "clusters",
      suffix = suffix
    )
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

  "Error Here"
  "Error Here"
  "Error Here"
  "Error Here"
  "Error Here"
  "Error Here"

  top.markers.df <- GetTopMarkersDF(dfDE = df_markers, order.by = lfcCOL, n = 1)
  top.markers <- top.markers.df %>% col2named.vec.tbl()

  missing.annotations <-
    top.markers.df %>%
    filter(!cluster %in% unique.matches$cluster) # filter for clusters that do not have a unique label already

  named.annotations <-
    rbind(unique.matches, missing.annotations) %>% # merge the 2 df's
    arrange(cluster) %>%
    CodeAndRoll2::col2named.vec.tbl()

  (top.markers.ID <- ppp(names(named.annotations), named.annotations))
  names(top.markers.ID) <- names(top.markers)
  named.ident <- top.markers.ID[Idents(object = obj)]

  namedIDslot <- ppp("cl.names.KnownMarkers", res)
  obj[[namedIDslot]] <- named.ident
  return(obj)
}



# _________________________________________________________________________________________________
# Correlations _________________________ ----
# _________________________________________________________________________________________________

#' @title Calculate Sparse Correlation Matrix
#'
#' @description Computes a sparse correlation matrix from a given sparse matrix input. This function is
#' useful for efficiently handling large datasets where most values are zero, facilitating the calculation
#' of both covariance and correlation matrices without converting to a dense format.
#'
#' @param smat A sparse matrix object, typically of class Matrix from the Matrix package.
#' @return A list with two elements:
#'   * `cov`: The covariance matrix derived from the input sparse matrix.
#'   * `cor`: The correlation matrix derived from the covariance matrix.
#'
#' @examples
#' \dontrun{
#' library(Matrix)
#' smat <- Matrix(rnorm(1000), nrow = 100, sparse = TRUE)
#' cor_res <- sparse.cor(smat)
#' print(cor_res$cor)
#' }
#'
#' @export
#' @importFrom Matrix colMeans crossprod tcrossprod
#' @importFrom stats sd
sparse.cor <- function(smat) {
  n <- nrow(smat)
  cMeans <- colMeans(smat)
  covmat <- (as.matrix(crossprod(smat)) - n * tcrossprod(cMeans)) / (n - 1)
  sdvec <- sqrt(diag(covmat))
  cormat <- covmat / tcrossprod(sdvec)
  list(cov = covmat, cor = cormat)
}


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
    assay.use = "RNA",
    slot.use = "data",
    quantileX = 0.95,
    max.cells = 40000,
    seed = p$"seed",
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

  if (is.null(obj@misc[[quantile_name]])) {
    iprint(
      "Call: combined.obj <- calc.q99.Expression.and.set.all.genes(combined.obj, quantileX =",
      quantileX, " first )"
    )
  }
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
#' @title Plot Gene Correlation Heatmap
#'
#' @description Generates a heatmap visualization of gene correlations based on expression data.
#' Useful for identifying groups of genes that exhibit similar expression patterns across different conditions
#' or cell types in a Seurat object.
#'
#' @param genes Vector of gene symbols to include in the correlation analysi.
#' @param assay.use Assay from which to retrieve expression data within the Seurat object; Default: 'RNA'.
#' @param slot.use Specifies which slot of the assay to use for expression data ('data', 'scale.data', 'data.imputed');
#' Default: first item ('data').
#' @param quantileX Quantile level for calculating expression thresholds; Default: 0.95.
#' @param min.g.cor Minimum absolute gene correlation value for inclusion in the heatmap; Default: 0.3.
#' @param calc.COR Logical flag to calculate correlation matrix if not found in `@misc`; Default: FALSE.
#' @param cutRows Height at which to cut the dendrogram for rows, determining cluster formation; Default: NULL.
#' @param cutCols Height at which to cut the dendrogram for columns, determining cluster formation;
#' Default: same as `cutRows`.
#' @param obj Seurat object containing the data; Default: `combined.obj`.
#' @param ... Additional parameters passed to the internally called functions.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   plot.Gene.Cor.Heatmap(genes = c("Gene1", "Gene2", "Gene3"), obj = combined.obj)
#' }
#' }
#'
#' @importFrom Seurat GetAssayData
#' @importFrom pheatmap pheatmap
#' @importFrom MarkdownReports wplot_save_pheatmap
#'
#' @export plot.Gene.Cor.Heatmap
#'
plot.Gene.Cor.Heatmap <- function(
    genes,
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
  MarkdownReports::wplot_save_pheatmap(o.heatmap, plotname = make.names(pname))

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
#' @title Add Prefixes to Cell Names in Seurat Objects
#'
#' @description Adds prefixes derived from a vector of identifiers to cell names in a list of Seurat objects.
#' This is useful for ensuring unique cell names across multiple samples or conditions when combining or comparing datasets.
#'
#' @param ls_obj List of Seurat S4 objects to which prefixes will be added. Each object should correspond
#' to a different sample or condition.
#' @param obj_IDs Character vector of identifiers that will be used as prefixes. Each identifier in the vector
#' corresponds to a Seurat object in `ls_obj`. The length of `obj_IDs` must match the length of `ls_obj`.
#'
#' @examples
#' \dontrun{
#' # Assuming seurat_obj1 and seurat_obj2 are Seurat objects
#' ls_obj <- list(seurat_obj1, seurat_obj2)
#' obj_IDs <- c("sample1", "sample2")
#' ls_obj_prefixed <- prefix_cells_seurat(ls_obj = ls_obj, obj_IDs = obj_IDs)
#' # Now each cell name in seurat_obj1 and seurat_obj2 will be prefixed with 'sample1_' and 'sample2_', respectively.
#' }
#'
#' @return A list of Seurat objects with updated cell names, incorporating the specified prefixes.
#'
#' @export
#' @importFrom Seurat RenameCells
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
#' @description This function checks if a prefix has been added to the standard
#' cell-IDs (16 characters of A,TRUE,C,G) in a Seurat object. If so, it prints the number of unique prefixes found,
#' issues a warning if more than one unique prefix is found, and returns the identified prefix(es).
#'
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
#' @title Create Cluster Labels for Each Cell
#'
#' @description Generates labels for each cell by combining gene names and cluster IDs. This function
#' takes a named vector, typically representing top genes for clusters (values) and their corresponding
#' cluster IDs (names), along with a vector of cell IDs. It then creates a new vector where each cell
#' is labeled with its top gene and cluster ID in the format "GeneName.ClusterID".
#'
#' @param TopGenes A named vector with gene names as values and cluster IDs as names,
#' representing the top or defining gene for each cluster.
#' @param clID.per.cell A vector of cluster IDs for each cell, used to match each cell with its
#' corresponding top gene from `TopGenes`.
#'
#' @return A vector where each element corresponds to a cell labeled with both its defining gene
#' name and cluster ID, in the format "GeneName.ClusterID".
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Assuming `TopGenes.Classic` is a named vector of top genes and cluster IDs,
#'   # and `metaD.CL.colname` is a column in metadata with cluster IDs per cell
#'   cellLabels <- seu.Make.Cl.Label.per.cell(
#'     TopGenes = TopGenes.Classic,
#'     clID.per.cell = getMetadataColumn(ColName.metadata = metaD.CL.colname)
#'   )
#'   # `cellLabels` now contains labels for each cell in the format "GeneName.ClusterID"
#' }
#' }
#'
#' @export
seu.Make.Cl.Label.per.cell <- function(TopGenes, clID.per.cell) {
  Cl.names_class <- TopGenes[clID.per.cell]
  Cl.names_wNr <- paste0(Cl.names_class, " (", names(Cl.names_class), ")")
  return(Cl.names_wNr)
}


# _________________________________________________________________________________________________
#' @title Retrieve the Top Variable Genes from a Seurat Object
#'
#' @description Retrieves the names of the most variable genes from a Seurat object,
#' typically used to focus subsequent analyses on genes with the greatest variation across cells.
#'
#' @param obj A Seurat object containing gene expression data and,
#' pre-computed highly variable gene information.
#' @param nGenes The number of most variable genes to retrieve; Default: `p$nVarGenes`.
#'
#' @return A vector containing the names of the most variable genes.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Assuming `combined.obj` is a Seurat object with computed variable genes
#'   varGenes <- GetMostVarGenes(obj = combined.obj, nGenes = 100)
#' }
#' }
#'
#' @export
#' @importFrom Seurat FindVariableFeatures
#'
GetMostVarGenes <- function(obj, nGenes = p$nVarGenes) {
  head(rownames(slot(object = obj, name = "hvg.info")), n = nGenes)
}

# _________________________________________________________________________________________________
#' @title Check Gene Names in Seurat Object
#'
#' @description Examines gene names in a Seurat object for specific naming conventions,
#' such as the presence of hyphens (-) or dots (.) often found in mitochondrial gene names.
#' This function is useful for ensuring gene names conform to expected patterns,
#' especially when preparing data for compatibility with other tools or databases.
#'
#' @param Seu.obj A Seurat object containing gene expression data.
#'
#' @details This function prints out examples of gene names that contain specific characters
#' of interest (e.g., '-', '_', '.', '.AS[1-9]'). It is primarily used for data inspection
#' and cleaning before further analysis or data export.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Assuming `combined.obj` is your Seurat object
#'   gene.name.check(Seu.obj = combined.obj)
#'   # This will print examples of gene names containing '-', '_', '.', and '.AS[1-9]'
#' }
#' }
#'
#' @seealso \code{\link[Seurat]{GetAssayData}}
#'
#' @importFrom Seurat GetAssayData
#' @importFrom CodeAndRoll2 grepv
#' @importFrom MarkdownHelpers llprint llogit
#'
#' @export
gene.name.check <- function(Seu.obj) {
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
#' @title Check if Gene Names exist in Seurat Object or HGNC Database
#'
#' @description Verifies the presence of specified gene names within a Seurat object or
#' queries them against the HGNC database. This function is useful for ensuring gene names are
#' correctly formatted and exist within the dataset or are recognized gene symbols.
#'
#' @param list.of.genes A vector of gene names to be checked; Default: `ClassicMarkers`.
#' @param makeuppercase If `TRUE`, converts all gene names to uppercase before checking; Default: `FALSE`.
#' @param verbose If `TRUE`, prints information about any missing genes; Default: `TRUE`.
#' @param HGNC.lookup If `TRUE`, attempts to look up any missing genes in the HGNC database to
#' verify their existence; Default: `FALSE`.
#' @param obj The Seurat object against which the gene names will be checked; Default: `combined.obj`.
#' @param assay.slot Assay slot of the Seurat object to check for gene names; Default: `'RNA'`.
#' @param dataslot Data slot of the assay to check for gene names; Default: `'data'`.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Check for the presence of a gene name in uppercase
#'   check.genes(list.of.genes = "top2a", makeuppercase = TRUE, obj = combined.obj)
#'
#'   # Check for a gene name with verbose output and HGNC lookup
#'   check.genes(list.of.genes = "VGLUT2", verbose = TRUE, HGNC.lookup = TRUE, obj = combined.obj)
#' }
#' }
#'
#' @seealso \code{\link[Seurat]{GetAssayData}}, \code{\link[DatabaseLinke.R]{qHGNC}}
#'
#' @export
#' @importFrom DatabaseLinke.R qHGNC
#' @importFrom Seurat GetAssayData
#' @importFrom Stringendo percentage_formatter
#'
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
      iprint(
        length(missingGenes), "or", Stringendo::percentage_formatter(length(missingGenes) / length(list.of.genes)),
        "genes not found in the data, e.g:", head(missingGenes, n = 10)
      )
    }
    if (HGNC.lookup) {
      if (exists("qHGNC", mode = "function")) {
        try(DatabaseLinke.R::qHGNC(missingGenes))
      } else {
        warning("DatabaseLinke.R's qHGNC() function is needed, please install from github.", immediate. = TRUE)
      }
    }
  }
  intersect(list.of.genes, all_genes)
}



# _________________________________________________________________________________________________
#' @title Fix Zero Indexing in Seurat Clustering
#'
#' @description Adjusts Seurat object metadata to fix zero-based cluster indexing, converting it to one-based indexing.
#' This function modifies a specified metadata column in the Seurat object to replace zero-indexed cluster names with one-based indexing.
#'
#' @param ColName.metadata The name of the metadata column containing zero-based cluster indices; Default: `'res.0.6'`.
#' @param obj The Seurat object to be modified; Default: `org`.
#'
#' @return The Seurat object with the specified metadata column's cluster indices adjusted to one-based indexing.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Assuming `org` is a Seurat object with zero-based cluster indexing
#'   org <- fixZeroIndexing.seurat(ColName.metadata = "res.0.6", obj = org)
#'   # Now, `org` has its cluster indices in the 'res.0.6' metadata column adjusted to one-based indexing
#' }
#' }
#'
#' @export
fixZeroIndexing.seurat <- function(ColName.metadata = "res.0.6", obj = org) {
  obj@meta.data[, ColName.metadata] <- as.numeric(obj@meta.data[, ColName.metadata]) + 1
  print(obj@meta.data[, ColName.metadata])
  return(obj)
}


# _________________________________________________________________________________________________
#' @title Calculate Fraction of Genes in Transcriptome
#'
#' @description Calculates the fraction of specified genes within the entire transcriptome of
#' each cell in a Seurat object.
#' This function is useful for assessing the relative abundance of a set of genes across cells,
#' such as identifying cells with high expression of marker genes.
#'
#' @param geneset A character vector of gene symbols for which the fraction in the transcriptome will be calculated.
#' Default: `c("MALAT1")`. The function will check for the existence of these genes in the Seurat object.
#' @param obj A Seurat object containing gene expression data; Default: `combined.obj`.
#' The function extracts gene expression data from this object to calculate fractions.
#' @param dataslot The data slot from which to extract expression data. This can be `"counts"`
#' for raw counts or `"data"` for normalized data; Default: second element (`"data"`).
#'
#' @return A numeric vector where each element represents the fraction of the specified geneset's expression
#' relative to the total transcriptome of a cell, expressed as a percentage. The names of the vector correspond to cell IDs.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Assuming `combined.obj` is your Seurat object
#'   fractionInTranscriptome <- CalculateFractionInTranscriptome(geneset = c("MALAT1", "GAPDH"), obj = combined.obj)
#'   # This will return the fraction of MALAT1 and GAPDH in the transcriptome of each cell
#' }
#' }
#'
#' @note This function calls `check.genes` to verify the existence of the specified genes within the Seurat object.
#' If genes are not found, it will return a warning.
#'
#' @seealso \code{\link[Seurat]{GetAssayData}} for retrieving expression data from a Seurat object.
#'
#' @export
#'
CalculateFractionInTrome <- function(
    genesCalc.Cor.Seuratet = c("MALAT1"),
    obj = combined.obj,
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
#' @description This function creates a new metadata column based on an existing metadata column
#' and a list of mappings (name <- IDs).
#' @param obj A Seurat object for which the new annotation is to be created. Default is 'obj'.
#' @param source A character string specifying the existing metadata column to be used as the
#' basis for the new annotation. Default is 'RNA_snn_res.0.5'.
#' @param named.list.of.identities A named list providing the mappings for the new annotation.
#' Default is 'ls.Subset.ClusterLists'.
#' @return A character vector representing the new metadata column.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   ls.Subset.ClusterLists <- list("hESC.h9" = c("4", "10", "14"), "hESC.176" = c("0", "1", "2"))
#'   AddNewAnnotation()
#' }
#' }
#' @export
AddNewAnnotation <- function(
    obj = obj,
    source = "RNA_snn_res.0.5", named.list.of.identities = ls.Subset.ClusterLists) {
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
#' @details The function first validates the presence of all identities from the metadata in the
#' Seurat objects. If all identities are present, the function subsets each Seurat object based on
#' the whitelist of cell IDs.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   ls.Seurat.subset <- whitelist.subset.ls.Seurat(
#'     ls.obj = ls.Seurat, metadir = p$"cellWhiteList",
#'     whitelist.file = "NonStressedCellIDs.2020.10.21_18h.tsv"
#'   )
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



#' @title Update Gene Symbols in a Seurat Object
#'
#' @description This function updates gene symbols in a Seurat object based on current gene
#' nomenclature guidelines, using HGNChelper(). It checks and updates gene symbols to their
#' latest approved versions,ensuring that gene annotations are current and consistent.
#' The function optionally enforces unique gene symbols and provides statistics on the update process.
#'
#' @param obj A Seurat object containing gene expression data; Default: `ls.Seurat[[i]]`
#' (ensure to replace `i` with the actual index or variable referencing your Seurat object).
#' @param species_ The species for which the gene symbols are checked and updated,
#' used to ensure the correct gene nomenclature is applied; Default: `'human'`.
#' Supports `'human'`, `'mouse'`, etc., as specified in the `HGNChelper` package.
#' @param EnforceUnique Logical flag indicating whether to enforce unique gene symbols
#' within the Seurat object. When set to `TRUE`, it resolves issues with duplicated gene symbols
#' by appending unique identifiers; Default: `TRUE`.
#' @param ShowStats Logical flag indicating whether to display statistics about the gene
#' symbol update process. When set to `TRUE`, it prints detailed information on the console
#' about the changes made; Default: `FALSE`.
#'
#' @return A modified Seurat object with updated gene symbols. The function directly modifies
#' the input Seurat object, ensuring that gene symbols adhere to the latest nomenclature.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Assuming `mySeuratObject` is your Seurat object
#'   updatedSeuratObject <- UpdateGenesSeurat(
#'     obj = mySeuratObject, species_ = "human",
#'     EnforceUnique = TRUE, ShowStats = TRUE
#'   )
#'   # `updatedSeuratObject` now has updated gene symbols
#' }
#' }
#'
#' @seealso
#' \code{\link[HGNChelper]{checkGeneSymbols}} for details on checking and updating gene symbols.
#'
#' @export
#' @importFrom HGNChelper checkGeneSymbols
#'
UpdateGenesSeurat <- function(obj = ls.Seurat[[i]], species_ = "human", EnforceUnique = TRUE, ShowStats = FALSE) {
  HGNC.updated <- HGNChelper::checkGeneSymbols(rownames(obj), unmapped.as.na = FALSE, map = NULL, species = species_)
  if (EnforceUnique) HGNC.updated <- HGNC.EnforceUnique(HGNC.updated)
  if (ShowStats) {
    print(HGNC.updated)
    print(GetUpdateStats(HGNC.updated))
  }
  obj <- RenameGenesSeurat(obj, newnames = HGNC.updated$"Suggested.Symbol")
  return(obj)
}


# _________________________________________________________________________________________________
#' @title Rename Gene Symbols in a Seurat Object
#'
#' @description This function replaces gene names across various slots within a specified assay
#' of a Seurat object. It is designed to be run prior to any data integration or downstream analysis
#' processes. The function targets the `@counts`, `@data`, and `@meta.features` slots within
#' the specified assay, ensuring consistency in gene nomenclature across the object.
#'
#' @param obj A Seurat object containing the assay and slots to be updated; Default: `ls.Seurat[[i]]`
#' (replace `i` with the appropriate index).
#' @param newnames A character vector containing the new gene names intended to replace the
#' existing ones; Default: `HGNC.updated[[i]]$Suggested.Symbol`. Ensure this matches the order
#' and length of the genes in the specified assay.
#' @param assay The name of the assay within the Seurat object where gene names will be updated;
#' Default: `"RNA"`. This function assumes simple objects containing only an RNA assay.
#' @param slots A character vector specifying which slots within the assay to update. Possible
#' values include `"data"`, `"counts"`, and `"meta.features"`; other layers can be specified if present.
#'
#' @details It is crucial to run this function before any data integration or further analysis
#' to ensure gene symbol consistency. The function does not support complex objects with multiple
#' assays where dependencies between assays might lead to inconsistencies. Use with caution and
#' verify the results.
#'
#' @note This function modifies the Seurat object in place, changing gene symbols directly within
#' the specified slots. Be sure to have a backup of your Seurat object if needed before applying
#' this function.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Assuming `SeuratObj` is your Seurat object
#'   # and `HGNC.updated.genes` contains the updated gene symbols
#'   SeuratObj <- RenameGenesSeurat(
#'     obj = SeuratObj,
#'     newnames = HGNC.updated.genes$Suggested.Symbol
#'   )
#'   # `SeuratObj` now has updated gene symbols in the specified assay and slots
#' }
#' }
#'
#' @export
RenameGenesSeurat <- function(obj = ls.Seurat[[i]],
                              newnames = HGNC.updated[[i]]$Suggested.Symbol,
                              assay = "RNA",
                              slots = c("data", "counts", "meta.features")) {
  message(assay)
  warning("Run this before integration and downstream processing. It only attempts to change
          @counts, @data, and @meta.features in obj@assays$YOUR_ASSAY.", immediate. = TRUE)

  stopifnot(
    "Unequal gene name sets: nrow(assayobj) != nrow(newnames):" =
      nrow(obj) == length(newnames)
  )

  if (obj@version < 5) warning("obj@version < 5. Old versions are not supported. Update the obj!", immediate. = TRUE)

  if ("scale.data" %in% slots) {
    n_genes_sc_dta <- nrow(obj@assays[[assay]]$"scale.data")
    stopifnot(
      "scale.data does has different number of genes than newnames!" =
        n_genes_sc_dta == length(newnames)
    )
  }

  LayersFound <- SeuratObject::Layers(obj@assays[[assay]])
  iprint("Present: ", sort(LayersFound))

  slots <- sort(intersect(slots, LayersFound))
  iprint("Replaced: ", slots)

  for (slotX in slots) {
    print(slotX)
    if (slotX == "scale.data") browser()
    nrO <- nrow(SeuratObject::GetAssayData(object = obj, assay = assay, layer = slotX))
    obj <- .check_and_rename(obj, assay, newnames = newnames, layer.name = slotX)
    nrN <- nrow(SeuratObject::GetAssayData(object = obj, assay = assay, layer = slotX))
    stopifnot(nrN == nrO)
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
#' # Assuming 'seurat_obj' is a Seurat object and 'new_gene_names' is a vector of gene names
#' updated_assay <- check_and_rename(
#'   assayobj = seurat_obj[["RNA"]],
#'   newnames = new_gene_names,
#'   layer.name = "counts"
#' )
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

  if (layer.name %in% SeuratObject::Layers(assayobj)) {
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
          " not of type dgeCMatrix / Matrix / data.frame.",
          immediate. = TRUE
        )
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
#' @title Remove Specific Genes from a Seurat Object
#'
#' @description Removes specified genes from the metadata, counts, data, and scale.data slots of a Seurat object.
#' This operation is typically performed prior to data integration to ensure that gene sets are consistent
#' across multiple datasets. The function modifies the Seurat object in place.
#'
#' @param obj A Seurat object; Default: `ls.Seurat[[i]]` (please ensure to replace `i` with the actual index or variable).
#' @param symbols2remove A character vector specifying the genes to be removed from the Seurat object;
#' Default: `c("TOP2A")`.
#'
#' @details This function directly modifies the `@counts`, `@data`, and `@scale.data` slots within
#' the RNA assay of the provided Seurat object, as well as the `@meta.data` slot. It's important to run
#' this function as one of the initial steps after creating the Seurat object and before proceeding
#' with downstream analyses or integration processes.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Assuming `SeuratObj` is your Seurat object and you want to remove the gene "TOP2A"
#'   updatedSeuratObj <- RemoveGenesSeurat(obj = SeuratObj, symbols2remove = "TOP2A")
#'   # Now `updatedSeuratObj` does not contain "TOP2A" in the specified slots
#' }
#' }
#'
#' @return A Seurat object with the specified genes removed from the mentioned slots.
#'
#' @export
RemoveGenesSeurat <- function(obj = ls.Seurat[[i]], symbols2remove = c("TOP2A")) {
  print("Run this as the first thing after creating the Seurat object.
        It only removes genes from: metadata; obj@assays$RNA@counts, @data and @scale.data.")
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
#' @title Enforce Unique HGNC Gene Symbols
#'
#' @description Ensures that gene symbols are unique after being updated with HGNC symbols. This function
#' applies a suffix to duplicate gene symbols to enforce uniqueness. While using `make.unique` might not
#' be the ideal solution due to potential mismatches, it significantly reduces the number of mismatching
#' genes in certain scenarios, making it a practical approach for data integration tasks.
#'
#' @param updatedSymbols A data frame or matrix containing gene symbols updated via `HGNChelper::checkGeneSymbols()`.
#' The third column should contain the updated gene symbols that are to be made unique.
#'
#' @return A modified version of the input data frame or matrix with unique gene symbols in the third column.
#' If duplicates were found, they are made unique by appending `.1`, `.2`, etc., to the repeated symbols.
#'
#' @details The function specifically targets the issue of duplicate gene symbols which can occur after
#' updating gene symbols to their latest HGNC-approved versions. Duplicate symbols can introduce
#' ambiguity in gene expression datasets, affecting downstream analyses like differential expression or
#' data integration. By ensuring each gene symbol is unique, this function helps maintain the integrity
#' of the dataset.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Assuming `SymUpd` is your data frame of updated symbols from HGNChelper::checkGeneSymbols()
#'   uniqueSymbols <- HGNC.EnforceUnique(updatedSymbols = SymUpd)
#'   # `uniqueSymbols` now contains unique gene symbols in its third column
#' }
#' }
#'
#' @note This function is a workaround for ensuring unique gene symbols and might not be suitable
#' for all datasets or analyses. It's important to review the results and ensure that the gene
#' symbols accurately represent your data.
#'
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
#' @title Gene Symbol Update Statistics
#'
#' @description Generates statistics on the gene symbol updates performed by `UpdateGenesSeurat()`.
#' This function analyzes the data frame of gene symbols before and after the update process,
#' providing insights into the proportion and total number of genes that were updated.
#'
#' @param genes A data frame of gene symbols before and after update, typically the output of
#' `UpdateGenesSeurat()`. Default: `HGNC.updated[[i]]` (where `i` is the index of the desired
#' Seurat object in a list).
#'
#' @return A named vector with statistics on gene updates, including the percentage of updated genes,
#' the absolute number of updated genes, and the total number of genes processed.
#'
#' @details The function examines the `Approved` column of the input data frame to identify
#' gene symbols marked for update and compares the original and suggested symbols to determine
#' actual updates. The statistics highlight the efficiency and impact of the gene symbol
#' updating process, aiding in the assessment of data preprocessing steps.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Assuming `HGNC.updated.genes` is your data frame containing the original and
#'   # suggested gene symbols, as returned by `UpdateGenesSeurat()`
#'   updateStats <- GetUpdateStats(genes = HGNC.updated.genes)
#'   # `updateStats` now contains the update statistics, including percentage and count of updated genes
#' }
#' }
#'
#' @note The function requires the input data frame to have specific columns as produced by
#' `HGNChelper::checkGeneSymbols()` and subsequently processed by `UpdateGenesSeurat()`.
#' Ensure that the input adheres to this format for accurate statistics.
#'
#' @seealso \code{\link{UpdateGenesSeurat}}, for the function that updates gene symbols and produces
#' the input data frame for this function.
#'
#' @importFrom Stringendo percentage_formatter
#'
#' @export
GetUpdateStats <- function(genes = HGNC.updated[[i]]) {
  MarkedAsUpdated <- genes[genes$Approved == FALSE, ]
  AcutallyUpdated <- sum(MarkedAsUpdated[, 1] != MarkedAsUpdated[, 3])
  UpdateStats <- c(
    "Updated (%)" = Stringendo::percentage_formatter(AcutallyUpdated / nrow(genes)),
    "Updated Genes" = floor(AcutallyUpdated), "Total Genes" = floor(nrow(genes))
  )
  return(UpdateStats)
}


# _________________________________________________________________________________________________
#' @title PlotUpdateStats
#'
#' @description Creates a scatter plot of update statistics.
#' @param mat A matrix containing update statistics. Default: UpdateStatMat.
#' @param column.names A character vector of column names in the mat parameter. Default: c("Updated (%)", "Updated (Nr.)").
#' @return A scatter plot displaying update statistics.
#' @details This function takes a matrix containing update statistics and column names to plot
#' the corresponding statistics. It colorizes the genes and plots the percentage of total genes
#' updated against the number of genes updated.
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
#' @description This function takes a parent directory with a number of subfolders, each
#' containing the standard output of 10X Cell Ranger. It (1) loads the filtered data matrices,
#' (2) converts them to Seurat objects, and (3) saves them as .RDS files.
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
    nthreads = .getNrCores(),
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
      LSBnameOUT <- ppp(paste0(InputDir, "/LSB.", fnameIN), "Rds")
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



# _________________________________________________________________________________________________
#' @title ConvertDropSeqfolders
#'
#' @description This function takes a parent directory with a number of subfolders, each
#' containing the standard output of 10X Cell Ranger. It (1) loads the filtered data matrices,
#' (2) converts them to Seurat objects, and (3) saves them as .RDS files.
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
    InputDir,
    folderPattern = "SRR*", filePattern = "expression.tsv.gz",
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
#' @description This function loads all Seurat objects found in a directory. It also works with
#' symbolic links (but not with aliases).
#' @param InputDir A character string specifying the input directory.
#' @param file.pattern A character string specifying the pattern of file names to be searched.
#' Default is '^filtered.+Rds$'.
#' @param string.remove1 A character string or FALSE. If a string is provided, it is removed from
#' file names. Default is "filtered_feature_bc_matrix.".
#' @param string.replace1 A character string of the new text instead of "string.remove1".
#' @param string.remove2 A character string or FALSE. If a string is provided, it is removed from
#' file names. Default is ".min.cells.10.min.features.200.Rds".
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
  stopifnot(length(fin.orig) > 0)
  fin <- if (!isFALSE(string.remove1)) sapply(fin.orig, gsub, pattern = string.remove1, replacement = string.replace1) else fin.orig
  fin <- if (!isFALSE(string.remove2)) sapply(fin, gsub, pattern = string.remove2, replacement = "") else fin
  if (sort_alphanumeric) fin <- gtools::mixedsort(fin)


  ls.Seu <- list.fromNames(fin)
  for (i in 1:length(fin)) {
    print(fin[i])
    FNP <- paste0(InputDir, fin.orig[i])
    # print(paste("Attempting to load file:", FNP))  # Debug print

    if (use_rds) {
      ls.Seu[[i]] <- readRDS(FNP)
    } else if (!use_rds) {
      ls.Seu[[i]] <- qs::qread(file = FNP)
    } else {
      warning("File pattern ambigous. Use either qs or rds:", file.pattern, immediate. = TRUE)
    }
  } # for
  print(tictoc::toc())
  return(ls.Seu)
}




# _________________________________________________________________________________________________
#' @title Load 10X Genomics Data as Seurat Object
#'
#' @description Reads 10X Genomics dataset files (gzipped) including matrix, features, and barcodes,
#' to a single expression matrix. This function handles the unzipping of these files, reads the data,
#' and re-compresses the files back to their original gzipped format.
#'
#' @param dir A character string specifying the path to the directory containing the 10X dataset files.
#' This directory should contain `matrix.mtx.gz`, `features.tsv.gz`, and `barcodes.tsv.gz` files.
#'
#' @return A Seurat object containing the single-cell RNA-seq data extracted from the provided 10X
#' Genomics dataset.
#'
#' @details This function facilitates the loading of 10X Genomics datasets into R for analysis with
#' the Seurat package. It specifically caters to gzipped versions of the `matrix.mtx`, `features.tsv`,
#' and `barcodes.tsv` files, automating their decompression, reading, and subsequent recompression.
#' The function relies on Seurat's `Read10X` function for data reading and object construction.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Replace `path_to_10x_data` with the path to your 10X data directory
#'   seuratObject <- read10x(dir = "path_to_10x_data")
#'   # `seuratObject` is now a Seurat object containing the loaded 10X data
#' }
#' }
#'
#' @note Ensure that the specified directory contains the required gzipped files.
#' If the `features.tsv.gz` file is named differently (e.g., `genes.tsv.gz`), please rename it
#' accordingly before running this function.
#'
#' @seealso \code{\link[Seurat]{Read10X}} for the underlying function used to read the 10X data.
#'
#' @importFrom tictoc tic toc
#' @importFrom R.utils gunzip gzip
#' @importFrom Seurat Read10X
#'
#' @export
read10x <- function(dir) {
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
#' @param inOutDir A boolean flag, if TRUE the OutDir is used as save directory, if FALSE the
#' alternative_path_rdata is used. Default is TRUE
#' @param project A string representing the project code. This is appended to the saved file name.
#' Default is the active project determined by getProject().
#' @param alternative_path_rdata A string that specifies the alternative path for storing the
#' RDS file if inOutDir is FALSE. Default is "~/Dropbox (VBC)/Abel.IMBA/AnalysisD/_RDS.files/"
#' appended with the basename of OutDir.
#' @param homepath A string representing the homepath. Will be replaced by '~' in the file path. Default is '~/'.
#' @param showMemObject A boolean flag, if TRUE the function will print out the memory size of the
#' largest objects in the workspace. Default is TRUE.
#' @param saveParams A boolean flag, if TRUE the parameters 'p' and 'all.genes' are added to the
#' 'misc' slot of the Seurat object if the object is of class Seurat. Default is TRUE.
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
  if ("Seurat" %in% is(obj) & saveParams) {
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
#' @param suffix Optional; a suffix to add to the filename.
#' @param prefix Optional; a prefix to add to the filename.
#' @param nthreads Number of threads to use when saving, defaults to 12.
#' @param preset Compression preset, defaults to 'high'.
#' @param project The project name to be included in the filename, defaults to the result of `getProject()`.
#' @param out_dir Output Directory
#' @param background_job Logical; if TRUE save runs as "background job"
#' @param showMemObject Logical; if TRUE, displays the memory size of the largest objects.
#' @param saveParams Logical; if TRUE and if the object is a Seurat object, additional parameters
#' are saved within it.
#' @param saveLocation Logical; if TRUE and if the object is a Seurat object, file location is saved
#' into misc slot.
#'
#' @return Invisible; The function is called for its side effects (saving a file) and does not return anything.
#'
#' @note The function uses the 'qs' package for quick and efficient serialization of objects and
#' includes a timing feature from the 'tictoc' package.
#' @seealso \code{\link[qs]{qsave}} for the underlying save function used.
#' @importFrom qs qsave
#' @importFrom tictoc tic toc
#' @importFrom job job
#' @importFrom rstudioapi isAvailable
#'
#' @export
xsave <- function(
    obj,
    suffix = NULL,
    prefix = NULL,
    nthreads = .getNrCores(12),
    preset = "high",
    project = getProject(),
    out_dir = if (exists("OutDir")) OutDir else getwd(),
    background_job = FALSE,
    showMemObject = TRUE, saveParams = TRUE,
    saveLocation = TRUE) {
  message(nthreads, " threads.")

  try(tictoc::tic(), silent = TRUE)
  if (showMemObject) {
    try(memory.biggest.objects(), silent = TRUE)
  }

  annot.suffix <- if (inherits(obj, "Seurat")) kpp(ncol(obj), "cells") else if (is.list(obj)) kppd("ls", length(obj)) else NULL
  fnameBase <- trimws(kppu(
    prefix, substitute(obj), annot.suffix, suffix, project,
    idate(Format = "%Y.%m.%d_%H.%M")
  ), whitespace = "_")

  FNN <- paste0(out_dir, fnameBase, ".qs")
  print(paste0(substitute(obj), " <- xread('", FNN, "')"))

  if ("Seurat" %in% is(obj)) {
    if (saveParams) {
      try(obj@misc$"p" <- p, silent = TRUE)
      try(obj@misc$"all.genes" <- all.genes, silent = TRUE)
    }
    if (saveLocation) {
      loc <- 1
      try(obj@misc$"file.location" <- loc, silent = TRUE)
    }
  }


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
xread <- function(file, nthreads = 4,
                  loadParamsAndAllGenes = TRUE,
                  overwriteParams = FALSE,
                  overwriteAllGenes = FALSE,
                  ...) {
  stopifnot(file.exists(file))

  message(nthreads, " threads.")
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
  obj <- qs::qread(file = file, nthreads = nthreads, ...)
  # }

  report <- if (is(obj, "Seurat")) {
    kppws("with", ncol(obj), "cells &", ncol(obj@meta.data), "meta colums.")
  } else if (is.list(obj)) {
    kppws("is a list of:", length(obj))
  } else {
    kppws("of length:", length(obj))
  }


  if ("Seurat" %in% is(obj)) {
    if (loadParamsAndAllGenes) {
      p_local <- obj@misc$"p"
      all.genes_local <- obj@misc$"all.genes"

      if (is.null(p_local)) {
        message("No parameter list 'p' found in object@misc.")
      } else {
          recall.parameters(obj = obj, overwrite = overwriteParams )
      }

      if (is.null(all.genes_local)) {
        message("No gene list 'all.genes' found in object@misc.")
      } else {
        recall.all.genes(obj = obj, overwrite = overwriteAllGenes)
      }
    }
  }


  iprint(is(obj)[1], report)
  try(tictoc::toc(), silent = TRUE)
  invisible(obj)
}


# _________________________________________________________________________________________________
#' @title Get the number of CPUs to use for CBE processing
#'
#' @description This function checks for the presence of a global `CBE.params` list and,
#' if found and contains a `cpus` entry, returns the number of CPUs specified by `cpus` minus one.
#' Otherwise, it returns a default number of CPUs.
#'
#' @param n.cpus.def The default number of CPUs to return if `CBE.params` does not exist
#' or does not contain a `cpus` entry. Defaults to 8.
#'
#' @return The number of CPUs to use for CBE processing. If `CBE.params$cpus` is set,
#' returns `CBE.params$cpus - 1`, ensuring at least 1 CPU is returned. Otherwise, returns `n.cpus.def`.
#'
#' @examples
#' # Assuming CBE.params does not exist or does not have a `cpus` entry
#' getCPUsCBE() # returns 8 by default
#'
#' # Assuming CBE.params exists and has a `cpus` entry of 4
#' getCPUsCBE() # returns 3
.getNrCores <- function(n.cpus.def = 8) {
  # Check if 'CBE.params' exists and contains 'cpus'
  if (exists("CBE.params") && is.list(CBE.params) &&
      is.numeric(CBE.params$"cpus") && CBE.params$"cpus" > 0) {
    max(CBE.params$"cpus" - 1, 1)
  } else {
    return(n.cpus.def)
  }
}



# _________________________________________________________________________________________________
# Save workspace
# requires MarkdownReports (github) and defining OutDir
# requires github/vertesy/CodeAndRoll.r

#' @title isave.image
#'
#' @description Save an image of the current workspace using a faster and efficient compression
#' method that runs in the background.
#' @param ... Additional parameters passed to the idate() function in the creation of the file name.
#' @param path_rdata A string that specifies the path for storing the image of the workspace.
#' Default is "~/Dropbox/Abel.IMBA/AnalysisD/_Rdata.files/" appended with the basename of OutDir.
#' @param showMemObject A boolean flag, if TRUE the function will print out the memory size of the
#' largest objects in the workspace. Default is TRUE.
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
    showMemObject = TRUE, options = c("--force", NULL)[1]) {
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
#' @description Faster saving of workspace, and compression outside R, when it can run in the background.
#' Seemingly quite CPU hungry and not very efficient compression. #
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @param options Options passed on to gzip, via CLI. Default: c("--force", NULL)[1]
#' @seealso
#'  \code{\link[Stringendo]{kollapse}}, \code{\link[function]{iprint}}
#' @export
#' @importFrom Stringendo kollapse iprint
#' @importFrom tictoc tic toc
qsave.image <- function(..., showMemObject = TRUE, options = c("--force", NULL)[1]) {
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
  stopifnot(
    is.character(root_dir), length(root_dir) == 1, dir.exists(root_dir),
    is.character(subdir), all(dir.exists(file.path(root_dir, subdir))),
    is.logical(recursive)
  )

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


# _________________________________________________________________________________________________
#' @title Clip Suffixes from 10X Cell Names
#'
#' @description Removes suffixes from cell names that are added by 10X technology and Seurat during data processing.
#'
#' @param cellnames A vector of cell names with potential suffixes.
#' @return A vector of cell names with suffixes removed.
#' @examples
#' cellnames <- c("cell1_1", "cell2_2")
#' clip10Xcellname(cellnames)
#' @export
#' @importFrom stringr str_split_fixed
clip10Xcellname <- function(cellnames) {
  stringr::str_split_fixed(cellnames, "_", n = 2)[, 1]
}

# _________________________________________________________________________________________________
#' @title Add Suffix to Cell Names (e.g. lane suffix: _1)
#'
#' @description Appends a specified suffix to cell names to mimic lane suffixes used in 10X datasets.
#'
#' @param cellnames A vector of cell names without numeric suffixes.
#' @param suffix The suffix to add to each cell name. Default is '_1'.
#' @return A vector of cell names with the specified suffix appended.
#' @examples
#' cellnames <- c("cell1", "cell2")
#' make10Xcellname(cellnames)
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
      Soup.VS.Cells.Av.Exp.gg |>
        arrange(-nchar(Class)), aes(x = Soup, y = Cells, label = gene, col = Class)
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
  ccc <- c("#FF4E00", "#778B04", "#8ea604", "#8ea604", "#F5BB00", "#F5BB00", "#EC9F05", rep(x = "#BF3100", times = NrColumns2Show - 6))


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
#' @description Create a pairwise jaccard similarity matrix across all combinations of columns in
#' binary.presence.matrix. Modified from:
#' https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/
#' @param lsG List of genes, Default: ls_genes
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   jPairwiseJaccardIndexList(lsG = ls_genes)
#' }
#' }
#' @export
#' @importFrom Stringendo percentage_formatter
jPairwiseJaccardIndexList <- function(lsG = ls_genes) {
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
#' @description Make a binary presence matrix from a list. Source:
#' https://stackoverflow.com/questions/56155707/r-how-to-create-a-binary-relation-matrix-from-a-list-of-strings #
#' @param string_list List of strings to compare overlapping entries.
#' Default: lst(a = 1:3, b = 2:5, c = 4:9, d = -1:4)
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   df.presence <- jPresenceMatrix(string_list = lst(a = 1:3, b = 2:5, c = 4:9, d = -1:4))
#' }
#' }
#' @export
jPresenceMatrix <- function(string_list = lst(a = 1:3, b = 2:5, c = 4:9, d = -1:4)) {
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
#' @description Calculate Jaccard Index. Modified from:
#' https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/ #
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
jJaccardIndexBinary <- function(x, y) {
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
#' @description Create a pairwise jaccard similarity matrix across all combinations of columns in
#' binary.presence.matrix. Modified from:
#' https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/
#' @param binary.presence.matrix A boolean matrix. Default: df.presence
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   PairwiseJaccardIndices <- jPairwiseJaccardIndex(binary.presence.matrix = df.presence)
#' }
#' }
#' @export
#' @importFrom Stringendo percentage_formatter
jPairwiseJaccardIndex <- function(binary.presence.matrix = df.presence) {
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
# Variable Features _____________________________ ------
# _________________________________________________________________________________________________

#' @title Compare variable features and their ranks in two Seurat objects.
#'
#' @description This function compares variable features (genes) between two Seurat objects,
#'   reporting the number of genes in each, the percentage of common genes, the percentage
#'   of unique genes in each object, and the similarity in the ranking of overlapping genes
#'   using Spearman's rank correlation coefficient. Optionally, it can also generate a scatterplot
#'   of the ranks of common genes using ggpubr's ggscatter. The function returns the common genes
#'   and the Spearman's rank correlation coefficient.
#'
#' @param obj1 The first Seurat object for comparison. Default: NULL.
#' @param obj2 The second Seurat object for comparison. Default: NULL.
#' @param cor.plot An optional boolean indicating whether to generate a scatterplot of the ranks
#'   of common genes. Default: FALSE.
#' @param plot_venn plot_venn
#' @param suffix suffix
#' @param save.plot save.plot
#' @return A list containing the common genes and Spearman's rank correlation coefficient.
#'   If cor.plot is TRUE, a scatterplot is also generated.
#' @importFrom Seurat VariableFeatures
#' @importFrom stats cor
#' @importFrom ggExpress qvenn qscatter
#' @examples
#' # Assuming obj1 and obj2 are Seurat objects
#' result <- compareVarFeaturesAndRanks(obj1, obj2, cor.plot = TRUE)
#' @export
compareVarFeaturesAndRanks <- function(
    obj1 = NULL, obj2 = NULL, cor.plot = TRUE, save.plot = TRUE,
    plot_venn = TRUE,
    suffix = NULL,
    ...) {
  stopifnot(!is.null(obj1), !is.null(obj2))
  stopifnot(is(obj1, "Seurat"), is(obj2, "Seurat"))

  name1 <- deparse(substitute(obj1))
  name2 <- deparse(substitute(obj2))

  var.genes1 <- Seurat::VariableFeatures(obj1)
  var.genes2 <- Seurat::VariableFeatures(obj2)

  if (plot_venn) {
    variable.genes.overlap <- list(var.genes1, var.genes2)
    names(variable.genes.overlap) <- c(name1, name2)
    ggExpress::qvenn(list = variable.genes.overlap, suffix = sppp(suffix, c(name1, name2)))
  }

  nr_genes1 <- length(var.genes1)
  nr_genes2 <- length(var.genes2)
  common_genes <- intersect(var.genes1, var.genes2)
  percent_common <- length(common_genes) / max(nr_genes1, nr_genes2) * 100
  percent_uniq1 <- (nr_genes1 - length(common_genes)) / nr_genes1 * 100
  percent_uniq2 <- (nr_genes2 - length(common_genes)) / nr_genes2 * 100

  ranks1 <- match(common_genes, var.genes1)
  ranks2 <- match(common_genes, var.genes2)

  spearman_correlation <- cor(ranks1, ranks2, method = "spearman")

  stopifnot(is.numeric(spearman_correlation))

  cat(sprintf("Nr of genes in obj1: %d\n", nr_genes1))
  cat(sprintf("Nr of genes in obj2: %d\n", nr_genes2))
  cat(sprintf("%% Common genes: %.2f%%\n", percent_common))
  cat(sprintf("%% Unique genes in obj1: %.2f%%\n", percent_uniq1))
  cat(sprintf("%% Unique genes in obj2: %.2f%%\n", percent_uniq2))
  cat(sprintf("Spearman's rank correlation: %.2f\n", spearman_correlation))

  if (cor.plot) {
    plot_data <- data.frame(ranks1, ranks2)
    colnames(plot_data) <- paste("Rank in", c(name1, name2))
    TTL <- paste("Spearman Rank Correlation of Shared Variable Genes")

    SUB <- paste(
      "between objects:", name1, "&", name2, "\n",
      length(common_genes), "or", percent_common, "% overlap from objects:",
      nr_genes1, "&", nr_genes2, "genes."
    )
    CPT <- paste("median ranks:", median(ranks1), "/", median(ranks2))
    file_name <- paste0(
      "Spearman_Rank_Correlation_of_",
      name1, "_and_", name2,
      "_", sprintf("%.2f", spearman_correlation), ".png"
    )
    print(head(plot_data))
    plt <- ggExpress::qscatter(
      df_XYcol = plot_data,
      plotname = TTL,
      subtitle = SUB,
      caption = CPT,
      # abline = c(0,1),
      save = save.plot,
      filename = file_name,
      correlation_r2 = TRUE,
      also.pdf = FALSE,
      cor.coef = TRUE, cor.method = "spearman",
      ...
    )
    print(plt)
  }

  unique.genes <- symdiff(var.genes1, var.genes2)
  names(unique.genes) <- paste0("Unique.", c(name1, name2))
  return(list(
    "common_genes" = common_genes,
    "unique.genes" = unique.genes,
    "spearman_correlation" = spearman_correlation
  ))
}


# _________________________________________________________________________________________________
# New additions,  categorized _____________________________ ------
# _________________________________________________________________________________________________

#' @title Process Seurat Objects in Parallel
#'
#' @description Applies a series of Seurat processing steps to each Seurat object in a list.
#'              The operations include scaling data, running PCA, UMAP, finding neighbors, and finding clusters.
#'              This is done in parallel using multiple cores.
#'
#' @param obj A Seurat object to be processed.
#' @param param.list A list of parameters used in the processing steps.
#' @return A Seurat object after applying scaling, PCA, UMAP, neighbor finding, and clustering.
#' @examples
#' # Assuming ls.Seurat is a list of Seurat objects and params is a list of parameters
#' # results <- mclapply(ls.Seurat, processSeuratObject, params, mc.cores = 4)
#' @importFrom Seurat ScaleData RunPCA RunUMAP FindNeighbors FindClusters
#' @export
processSeuratObject <- function(obj, param.list = p, compute = TRUE,
                                save = TRUE, plot = TRUE,
                                nfeatures = param.list$"n.var.genes") {
  warning("Make sure you cleaned up the memory!", immediate. = TRUE)
  stopifnot(require(tictoc))
  message("nfeatures: ", nfeatures)

  # Assertions to check input types
  stopifnot(
    "Seurat" %in% class(obj),
    is.list(param.list),
    all(c("n.PC", "snn_res") %in% names(param.list)),
    is.numeric(param.list$"n.PC"),
    is.numeric(param.list$"snn_res"),
    is.character(param.list$"variables.2.regress") | is.null(param.list$"variables.2.regress")
  )
  .checkListElements(param_list = p, elements = c("variables.2.regress.combined", "n.PC", "snn_res"))


  gc()
  if (compute) {
    message("------------------- FindVariableFeatures -------------------")
    tic()
    obj <- FindVariableFeatures(obj, mean.function = "FastExpMean", dispersion.function = "FastLogVMR", nfeatures = nfeatures)
    toc()
    obj <- calc.q99.Expression.and.set.all.genes(obj = obj, quantileX = .99)
    toc()
    message("------------------- ScaleData -------------------")
    tic()
    obj <- ScaleData(obj, assay = "RNA", verbose = TRUE, vars.to.regress = param.list$"variables.2.regress.combined")
    message("------------------- PCA /UMAP -------------------")
    tic()
    obj <- RunPCA(obj, npcs = param.list$"n.PC", verbose = TRUE)
    toc()
    tic()
    obj <- RunUMAP(obj, reduction = "pca", dims = 1:param.list$"n.PC")
    toc()
    message("------------------- FindNeighbors & Clusters -------------------")
    tic()
    obj <- FindNeighbors(obj, reduction = "pca", dims = 1:param.list$"n.PC")
    toc()
    tic()
    obj <- FindClusters(obj, resolution = param.list$"snn_res")
    toc()
  }

  message("------------------- Save -------------------")
  if (save) xsave(obj, suffix = "reprocessed")

  if (plot) {
    message("scPlotPCAvarExplained")
    scPlotPCAvarExplained(obj)

    message("qQC.plots.BrainOrg")
    qQC.plots.BrainOrg(obj = obj)

    message("multi_clUMAP.A4")
    multi_clUMAP.A4(obj = obj)

    message("qClusteringUMAPS")
    qClusteringUMAPS(obj = obj)

    message("suPlotVariableFeatures")
    suPlotVariableFeatures(obj = obj)

    if (ncol(obj) < 50000) { # TEMP
      message("qMarkerCheck.BrainOrg")
      try(qMarkerCheck.BrainOrg(obj = obj), silent = TRUE)
    }

    Signature.Genes.Top20 <- c(
      `dl-EN` = "KAZN", `ul-EN` = "SATB2" # dl-EN = deep layer excitatory neuron
      , `Immature neurons` = "SLA", Interneurons = "DLX6-AS1",
      Interneurons = "ERBB4", Interneurons = "SCGN",
      `Intermediate progenitor` = "EOMES" # ,  `Intermediate progenitor1` = "TAC3"
      , `S-phase` = "TOP2A", `G2M-phase` = "H4C3" # formerly: HIST1H4C
      , `oRG` = "HOPX", `oRG` = "ID4" # oRG outer radial glia
      , Astroglia = "GFAP",
      Astrocyte = "S100B", `Hypoxia/Stress` = "DDIT4",
      `Choroid.Plexus` = "TTR", `Low-Quality` = "POLR2A",
      `Mesenchyme` = "DCN", Glycolytic = "PDK1",
      `Choroid.Plexus` = "OTX2", `Mesenchyme` = "DCN"
    )
    message("plotQUMAPsInAFolder")
    try(plotQUMAPsInAFolder(genes = Signature.Genes.Top20, obj = obj), silent = TRUE)
  }
  tictoc::toc()

  return(obj)
}



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
  .Deprecated("processSeuratObject")

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
#' @title Check List Elements
#'
#' @description Tests if list elements are defined and reports the value or issues a warning.
#'
#' @param param_list A list containing variables to be checked. Default: `NULL`.
#' @param elements A character vector of element names in `param_list` to check.
#' Default: `character(0)`.
#'
#' @return A message for each element that is defined, and a warning for elements that are not.
#' @examples
#' param_list <- list(a = 1, b = NULL)
#' elements <- c("a", "b", "c")
#' .checkListElements(param_list, elements)
.checkListElements <- function(param_list = NULL, elements = character(0)) {
  stopifnot(
    is.list(param_list),
    is.character(elements)
  )

  sapply(elements, function(element) {
    if (is.null(param_list[[element]])) {
      warning(sprintf("`%s` is not defined", element), immediate. = TRUE, call. = FALSE)
    } else {
      message(sprintf("`%s` is: %s", element, param_list[[element]]))
    }
  }, USE.NAMES = FALSE)

  invisible(NULL)
}



# _________________________________________________________________________________________________



# _________________________________________________________________________________________________
#' @title Get number of scaled features
#'
#' @param obj A Seurat object containing scaled data in  `obj@assays$RNA@scale.data`.
#' @return Integer representing the number of scaled features
.getNrScaledFeatures <- function(obj, assay = Seurat::DefaultAssay(obj)) {

  message("Seurat version: ", obj@version)
  message("Assay searched: ", assay)

  # !!! Below may have been necessary bc of a bug in version 5.0.0
  if (obj@version >= 5) {
    if ("scale.data" %in% names(obj@assays[[assay]])) {
      nrow(obj@assays[[assay]]@"scale.data")
    } else {
      nrow(obj@assays[[assay]]@layers$"scale.data")
    }
  } else {
    nrow(obj@assays[[assay]]@"scale.data")
  }

}



# _________________________________________________________________________________________________
#' @title Get number of principal components
#'
#' @param obj A Seurat object containing PCA cell embeddings in `reductions$pca@cell.embeddings`
#' @return Integer representing the number of principal components
.getNrPCs <- function(obj) {
  ncol(obj@reductions$pca@"cell.embeddings")
}

# _________________________________________________________________________________________________
#' @title Parse regression variables for name tagging
#'
#' @param p list of parameters
#' @return Integer representing the number of principal components
.parseRegressionVariablesForScaleData <- function(element = "variables.2.regress.combined", par.list = p) {
  (regV <- par.list[[element]])
  txt <- if (is.null(regV)) "No.Regr" else kpp("Regr", regV)
  return(txt)
}


# _________________________________________________________________________________________________
#' @title Parse key parameters from an object and format as a string
#'
#' @description This function extracts the number of scaled features, the number of principal components,
#' and formats additional information including regression variables.
#' @param obj An object to extract information from.
#' @param regressionVariables A list or vector containing variables for regression.
#' @param nrVarFeatures You can provide this number manually. Default: NULL.
#' @param suffix A suffix string to add.
#' @return A character string summarizing the key parameters.
#'
.parseKeyParams <- function(obj,
                            regressionVariables = obj@misc$p$"variables.2.regress.combined",
                            nrVarFeatures = NULL,
                            return.as.name = FALSE,
                            assay = Seurat::DefaultAssay(obj),
                            suffix = NULL) {
  scaledFeatures <- .getNrScaledFeatures(obj, assay = assay)

  if (!is.null(nrVarFeatures)) {
    if (nrVarFeatures != scaledFeatures) {
      warning("nrVarFeatures !=  scaledFeatures. Reporting nrVarFeatures: ", nrVarFeatures, immediate. = TRUE)
    }
    scaledFeatures <- nrVarFeatures
  } # else use scaledFeatures

  pcs <- .getNrPCs(obj)
  regressionInfo <- kppc(regressionVariables)
  reg <- if (!is.null(regressionVariables)) paste0(" regress ", regressionInfo) else "no regression"
  if (return.as.name) {
    reg <- ReplaceSpecialCharacters(RemoveWhitespaces(reg, replacement = "."))
    tag <- kpp(scaledFeatures, "ScaledFeatures", pcs, "PCs", reg, suffix)
  } else {
    tag <- paste0(scaledFeatures, " ScaledFeatures | ", pcs, " PCs | ", reg, " ", suffix)
  }

  return(tag)
}



# _________________________________________________________________________________________________
#' @title Parse basic obj stats
#'
#' @description Parse cell and feature count to a string.
#' @param obj An object to extract information from.
#' @return A character string summarizing the key parameters.
#'
.parseBasicObjStats <- function(obj, sep = " ", simple = FALSE, suffix = NULL) {
  n.cells <- format(ncol(obj), big.mark = sep, scientific = FALSE)
  n.feat <- format(nrow(obj), big.mark = sep, scientific = FALSE)
    if (simple) {
      return(paste(n.cells, "cells.", suffix))
    } else {
      return(paste(n.cells, "cells,", n.feat, "features.", suffix))
    }
}




# _________________________________________________________________________________________________
# Temp _____________________________ ------
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
#' @param seuratObj A Seurat object.
#' @param pattern A regular expression pattern to match layer names.
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
