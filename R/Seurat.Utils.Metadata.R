# ____________________________________________________________________
# Seurat.Utils.Metadata ----
# ____________________________________________________________________
# source("~/GitHub/Packages/Seurat.utils/R/Seurat.Utils.Metadata.R")
# devtools::load_all(path = '~/GitHub/Packages/Seurat.utils');
# devtools::document("~/GitHub/Packages/Seurat.utils"); devtools::load_all("~/GitHub/Packages/Seurat.utils")
# source("~/GitHub/Packages/Seurat.utils/R/Seurat.Utils.R")
# source('~/.pack.R')



# _________________________________________________________________________________________________
# Extract and check metadata columns  ______________________________ ----
# _________________________________________________________________________________________________



#' @title Get Metadata Column Names Matching Pattern
#'
#' @description Retrieves column names from an object's metadata that match a specified pattern.
#'
#' @param pattern A character string containing a regular expression to match against the
#' column names. Default: "RNA".
#' @param obj An object containing a `meta.data` slot, typically from combined datasets.
#' Default: `combined.obj`.
#'
#' @return A character vector of column names matching the pattern.
#'
#' @examples
#' # Assuming `combined.obj` is an object with a meta.data slot
#' getMetaColnames()
#'
#' @export
getMetaColnames <- function(obj = combined.obj,
                            pattern = "RNA") {
  stopifnot(inherits(obj, "Seurat"))

  # Retrieve column names matching the pattern
  matchedColnames <- grep(pattern = pattern, x = colnames(obj@meta.data), value = TRUE)

  # Output assertion
  if (is.null(matchedColnames)) {
    warning("No matching meta data!", immediate. = TRUE)
  } else {
    message(length(matchedColnames), " columns matching pattern '", pattern, "'.")
  }

  dput(matchedColnames)
  return(matchedColnames)
}


# _________________________________________________________________________________________________
#' @title Check if a Column Exists in the Metadata of an S4 Object
#'
#' @description This function checks whether a given column exists in the meta.data of a Seurat object.
#' @param obj A Seurat object.
#' @param col_name A character string specifying the name of the column.
#'
#' @return A logical value indicating whether the column exists (TRUE) or not (FALSE).
#' @export
metaColnameExists <- function(col_name, obj = combined.obj) {
  col_name %in% colnames(obj@meta.data)
}

# _________________________________________________________________________________________________
#' @title getMetadataColumn
#'
#' @description Retrieves a specified metadata column from a Seurat object and returns it as a named vector.
#' @param ColName.metadata A string specifying the name of the metadata column to be retrieved. Default: 'batch'.
#' @param obj A Seurat object from which the metadata column will be retrieved. Default: combined.obj.
#' @param as_numeric A logical flag indicating whether the returned values should be converted to numeric format. Default: FALSE (FALSE).
#' @return A named vector containing the values from the specified metadata column. If 'as_numeric' is TRUE, the values are converted to numeric format.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Example usage:
#'   batch_metadata <- getMetadataColumn(ColName.metadata = "batch", obj = combined.obj, as_numeric = TRUE)
#' }
#' }
#' @export
getMetadataColumn <- function(ColName.metadata = "batch", obj = combined.obj, as_numeric = FALSE) {
  stopifnot(ColName.metadata %in% colnames(obj@meta.data))

  x <- as.named.vector.df(obj@meta.data[, ColName.metadata, drop = FALSE])
  if (as_numeric) {
    as.numeric.wNames(x) + 1
  } else {
    x
  }
}

# mmeta <- getMetadataColumn

# _________________________________________________________________________________________________
#' @title Get Unique Levels of a Seurat Object Ident Slot
#'
#' @description This function extracts the unique levels present in the 'ident' slot of a Seurat object.
#' The function throws an error if the number of levels exceeds 'max_levels'.
#' The function optionally prints the R code to recreate the 'Levels' vector using 'dput'.
#'
#' @param obj A Seurat object.
#' @param ident A character string representing the name of the slot in the Seurat object.
#' @param max_levels An integer that sets the maximum number of levels allowed. Default is 100.
#' @param dput A logical that decides whether to print the R code to recreate the 'Levels' vector. Default is TRUE.
#'
#' @return A vector of unique levels present in the 'ident' slot of the Seurat object.
#'
#' @importFrom tibble deframe
#' @import Seurat
#' @export

get_levels_seu <- function(obj, ident, max_levels = 100, dput = TRUE) {
  Levels <- unique(deframe(obj[[ident]]))
  stopifnot(length(Levels) < max_levels)
  if (dput) {
    cat("Levels <- ")
    dput(Levels)
  }
  return(Levels)
}


#' @title Calculate Average Metadata for Seurat Object
#'
#' @description Computes specified metrics (e.g., median, mean) for given metadata features across each category
#' defined by an identity column in a Seurat object's metadata. This function allows for flexible
#' metric calculation on specified features, providing insights into the data distribution.
#'
#' @param obj A Seurat object containing metadata to be analyzed. Defaults to `combined.obj`.
#' @param meta.features A character vector specifying which metadata features to calculate metrics for.
#'   Defaults to c("nFeature_RNA", "percent.ribo", "percent.mito").
#' @param ident The name of the identity column used to group the data before calculating metrics.
#'   The default is the second entry from `GetNamedClusteringRuns()`.
#' @param metrics A list of named metrics to calculate for the metadata features, where names are
#'   the metric names (e.g., 'median', 'mean') and values are the corresponding functions.
#'   Defaults to list('median' = median, 'mean' = mean).
#' @param verbose Logical flag indicating whether to print detailed information about the metrics
#'   calculation process. Defaults to TRUE.
#' @param max.categ max number of groups in ident.
#'
#' @return A list containing data frames with calculated metrics for each specified metadata feature,
#'   grouped by the identity categories. Each data frame corresponds to one of the specified metrics.
#'
#' @examples
#' # Assuming `combined.obj` is a Seurat object with relevant metadata columns:
#' results <- calculateAverageMetaData(
#'   obj = combined.obj,
#'   meta.features = c("nFeature_RNA", "nCount_RNA"),
#'   metrics = list("median" = median, "mean" = mean),
#'   verbose = TRUE
#' )
#' # This will return a list with data frames containing the median and mean
#' # of "nFeature_RNA" and "percent.ribo" for each category in "ident_column_name".
#'
#' @export
calculateAverageMetaData <- function(
    obj = combined.obj,
    meta.features = c("nFeature_RNA", "percent.ribo", "percent.mito"),
    ident = GetClusteringRuns()[1],
    metrics = list("median" = median, "mean" = mean),
    verbose = TRUE, max.categ = 30) {
  stopifnot(
    is(obj, "Seurat"),
    "ident not found in object" = ident %in% colnames(obj@meta.data),
    "Not all meta.features found in object" = all(meta.features %in% colnames(obj@meta.data)),
    length(unique(obj@meta.data[, ident])) < max.categ
  )

  # browser()
  # Initialize list to store results
  results <- list()

  # Calculate metrics for each meta.feature within each ident category
  for (m in names(metrics)) {
    results[[m]] <- obj@meta.data %>%
      group_by(!!sym(ident)) %>%
      summarise(across(all_of(meta.features), metrics[[m]], na.rm = TRUE), .groups = "drop")
  }

  # Verbose output
  if (verbose) {
    cat("Calculated metrics:", paste(names(metrics), collapse = ", "), "\n")
    cat("For features:", paste(meta.features, collapse = ", "), "\n")
    cat("Based on identifier:", ident, "\n")
  }
  return(results)
}


# _________________________________________________________________________________________________
#' @title Calculate the Percentage of Matches per Category
#'
#' @description This function calculates the percentage of matches for specified metadata features
#' against provided match values within each category of an identifier in a Seurat object.
#'
#' @param obj A Seurat object containing the data to be analyzed. Default: combined.obj.
#' @param ident A string specifying the column in the metadata that identifies the categories.
#' Default: first element of `GetClusteringRuns()`.
#' @param meta.features A vector of strings specifying which metadata features to analyze.
#' Default: c("AAV.detected.min2", "AAV.detected").
#' @param match.values A named vector where names correspond to `meta.features` and values are
#' the strings to match against. Default: c("AAV.detected.min2" = "AAV", "AAV.detected" = "AAV").
#' @param verbose A logical value indicating whether to print detailed output. Default: TRUE.
#' @param max.categ The maximum number of categories allowed before stopping. Default: 30.
#'
#' @return A data frame with the category as the first column and the subsequent columns showing
#' the percentage of matches for each metadata feature.
#' @export
#'
#' @examples
#' calculatePercentageMatch(obj = combined.obj, ident = "Simple_Celltypes")
calculatePercentageMatch <- function(
    obj,
    ident = GetClusteringRuns()[1],
    meta.features = c("AAV.detected.min2", "AAV.detected"),
    match.values = c("AAV.detected.min2" = "AAV", "AAV.detected" = "AAV"), # Named vector for matches
    verbose = TRUE,
    max.categ = 100) {
  # Check for preconditions
  stopifnot(
    is(obj, "Seurat"),
    "ident not found in object" = ident %in% colnames(obj@meta.data),
    "Not all meta.features found in object" = all(meta.features %in% colnames(obj@meta.data)),
    "Too many categories" = length(unique(obj@meta.data[, ident])) < max.categ,
    length(match.values) == length(meta.features), # Check if match.values has the same length as meta.features
    all(names(match.values) == meta.features) # Ensure match.values has names corresponding to meta.features
  )

  # Initialize a data frame to store results
  results <- data.frame(Category = unique(obj@meta.data[[ident]]))

  # Calculate the percentage of matches for each meta.feature within each ident category
  for (feature in meta.features) {
    results[[paste0("pct_match_", feature)]] <- sapply(results$Category, function(cat) {
      idx <- obj@meta.data[[ident]] == cat
      subset_data <- obj@meta.data[idx, feature, drop = TRUE]
      pct_match <- mean(subset_data == match.values[feature], na.rm = TRUE)
      return(pct_match)
    })
  }

  # Verbose output
  if (verbose) {
    cat("Calculated percentage of matches for values:", paste(match.values, collapse = ", "), "\n")
    cat("Corresponding to features:", paste(meta.features, collapse = ", "), "\n")
    cat("Based on identifier:", ident, "\n")
  }

  return(results)
}



# _________________________________________________________________________________________________
#' @title getMedianMetric.lsObj
#'
#' @description Get the median values of different columns in meta.data, can iterate over a list of Seurat objects.
#' @param ls.obj List of Seurat objects, Default: ls.Seurat
#' @param n.datasets lenght of list (n objects), Default: length(ls.Seurat)
#' @param mColname Metadata column name to calculate on. Default: 'percent.mito'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   ls.Seurat <- getMedianMetric.lsObj(
#'     ls.obj = ls.Seurat, n.datasets = length(ls.Seurat),
#'     mColname = "percent.mito"
#'   )
#' }
#' }
#' @export
getMedianMetric.lsObj <- function(ls.obj = ls.Seurat, n.datasets = length(ls.Seurat), mColname = "percent.mito") {
  medMetric <- vec.fromNames(names(ls.obj))
  for (i in 1:n.datasets) {
    medMetric[i] <- median(ls.obj[[i]]@meta.data[, mColname])
  }
  return(medMetric)
}



# _________________________________________________________________________________________________
#' @title getCellIDs.from.meta
#'
#' @description Retrieves cell IDs from a specified metadata column of a Seurat object, where the cell ID matches a provided list of values. The matching operation uses the `%in%` operator.
#' @param ident A string specifying the name of the metadata column from which to retrieve cell IDs. Default: 'res.0.6'.
#' @param ident_values A vector of values to match in the metadata column. Default: NA.
#' @param obj The Seurat object from which to retrieve the cell IDs. Default: combined.obj.
#' @param inverse A boolean value indicating whether to inverse the match, i.e., retrieve cell IDs that do not match the provided list of ident_values. Default: FALSE.
#' @return A vector of cell IDs that match (or don't match, if `inverse = TRUE`) the provided list of values.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Example usage:
#'   getCellIDs.from.meta()
#' }
#' }
#' @export
getCellIDs.from.meta <- function(ident = GetClusteringRuns()[1],
                                 ident_values = NA, obj = combined.obj,
                                 inverse = FALSE) {
  # browser()
  mdat <- obj@meta.data[, ident]
  cells.pass <- mdat %in% ident_values
  if (inverse) cells.pass <- !cells.pass

  iprint(sum(cells.pass), "cells found.")
  return(rownames(obj@meta.data)[which(cells.pass)])
}


# _________________________________________________________________________________________________
# Add new metadata  ______________________________ ----
# _________________________________________________________________________________________________


#' @title Add Metadata to a Seurat object, safely with Checks
#'
#' @description Wrapper function for `AddMetaData` that includes additional checks and assertions.
#'
#' @param obj Seurat object to which metadata will be added.
#' @param metadata The metadata to be added.
#' @param col.name The name of the new metadata column.
#' @param overwrite Logical; if TRUE, overwrites the existing column.
#'
#' @return Modified Seurat object with additional metadata.
#' @importFrom Seurat AddMetaData
#' @export
addMetaDataSafe <- function(obj, metadata, col.name, overwrite = FALSE) {

  message("Running addMetaDataSafe...")
  # browser()
  stopifnot(
    is(obj, "Seurat"),
    is.vector(metadata),
    is.character(col.name),
    is.logical(overwrite),
    "Column already exists" = ((!col.name %in% colnames(obj@meta.data)) | overwrite),
    "Check length" = (length(metadata) == ncol(obj))
  )

  if (!is.null(names(metadata))) {
    print(head(names(metadata)))
    print(head(colnames(obj)))
    stopifnot(names(metadata) == colnames(obj))
  } else {
    message("No CBCs associated with new metadata. Assuming exact match.")
    names(metadata) <- colnames(obj)
  }


  # Perform the operation
  obj <- Seurat::AddMetaData(object = obj, metadata = metadata, col.name = col.name)

  # Check for NA or NaN values
  if (all(is.na(metadata) | is.nan(metadata))) {
    warning("New metadata column contains only NA or NaN values.")
  } else if (any(is.na(metadata) | is.nan(metadata))) {
    message("New metadata column contains NA or NaN values.")
  }

  return(obj)
}


# _________________________________________________________________________________________________
#' @title seu.add.meta.from.vector
#'
#' @description Adds a new metadata column to a Seurat object.
#' @param obj A Seurat object to which the new metadata column will be added. Default: combined.obj.
#' @param metaD.colname A string specifying the name of the new metadata column. Default: metaD.colname.labeled.
#' @param Label.per.cell A vector of labels for each cell, to be added as new metadata. Default: Cl.Label.per.cell.
#' @return A Seurat object with the new metadata column added.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Example usage:
#'   combined.obj <- seu.add.meta.from.vector(
#'     obj = combined.obj,
#'     metaD.colname = metaD.colname.labeled,
#'     Label.per.cell = Cl.Label.per.cell
#'   )
#' }
#' }
#' @export
seu.add.meta.from.vector <- function(obj = combined.obj, metaD.colname = metaD.colname.labeled, Label.per.cell = Cl.Label.per.cell) { # Add a new metadata column to a Seurat  object
  obj@meta.data[, metaD.colname] <- Label.per.cell
  iprint(metaD.colname, "contains the named identitites. Use Idents(combined.obj) = '...'. The names are:", unique(Label.per.cell))
  return(obj)
}



# _________________________________________________________________________________________________
#' @title Create a Metadata Vector
#'
#' @description This function creates a metadata vector from an input vector and a Seurat object.
#' The resulting vector contains values from 'vec' for the intersecting cell names between 'vec' and 'obj'.
#' It also checks if the intersection between the cell names in 'vec' and 'obj' is more than a minimum intersection size.
#' @param vec A named vector where the names represent cell IDs. This vector should have partial overlap with the cells in a Seurat object. Default is 'All.UVI'.
#' @param obj A Seurat object that contains cell IDs which partially overlap with 'vec'. Default is 'combined.obj'.
#' @param min.intersect The minimum number of cells to find in both 'vec' and 'obj'. The function will stop if the intersection is less than this number. Default is 100.
#' @return A named vector of length equal to the number of cells in 'obj', with names from 'obj' and values from 'vec' for intersecting cell names.
#' @examples
#' \dontrun{
#' create.metadata.vector(vec = my_vector, obj = my_seurat_object, min.intersect = 50)
#' }
#' @export
create.metadata.vector <- function(vec = All.UVI, obj = combined.obj, min.intersect = 100) {
  cells.vec <- names(vec)
  cells.obj <- colnames(obj)
  cells.in.both <- intersect(cells.vec, cells.obj)

  # iprint("intersect:", length(cells.in.both), head(cells.in.both))
  iprint(
    length(cells.in.both), "cells in both;",
    length(cells.vec), "cells in vec;",
    length(cells.obj), "cells in obj",
    "intersect, e.g.:", head(cells.in.both, 5)
  )
  stopifnot(length(cells.in.both) > min.intersect)

  new_assignment <- vec.fromNames(cells.obj)
  new_assignment[cells.in.both] <- vec[cells.in.both]
  return(new_assignment)
}


# _________________________________________________________________________________________________
#' @title addMetaFraction
#'
#' @description Add a new metadata column to a Seurat object, representing the fraction of a gene set in the transcriptome (expressed as a percentage).
#' @param col.name Name of the new metadata column to be added. Default: 'percent.mito'
#' @param gene.symbol.pattern Regular expression pattern to match gene symbols. Default: c("^MT\\.|^MT-", FALSE)[1]
#' @param gene.set A set of gene symbols. If specified, it will be used instead of gene.symbol.pattern. Default: FALSE
#' @param obj Seurat object to which the new metadata column will be added. Default: ls.Seurat[[1]]
#' @param verbose Logical indicating whether to display detailed messages (TRUE) or not (FALSE). Default: TRUE
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   ls.Seurat[[1]] <- addMetaFraction(col.name = "percent.mito", gene.symbol.pattern = "^MT\\.|^MT-")
#'   ls.Seurat[[1]] <- addMetaFraction(col.name = "percent.ribo", gene.symbol.pattern = "^RPL|^RPS")
#'   ls.Seurat[[1]] <- addMetaFraction(col.name = "percent.AC.GenBank", gene.symbol.pattern = "^AC[0-9]{6}\\.")
#'   ls.Seurat[[1]] <- addMetaFraction(col.name = "percent.AL.EMBL", gene.symbol.pattern = "^AL[0-9]{6}\\.")
#'   ls.Seurat[[1]] <- addMetaFraction(col.name = "percent.LINC", gene.symbol.pattern = "^LINC0")
#'   ls.Seurat[[1]] <- addMetaFraction(col.name = "percent.MALAT1", gene.symbol.pattern = "^MALAT1")
#'   colnames(ls.Seurat[[1]]@meta.data)
#'   HGA_MarkerGenes <- c(
#'     "ENO1", "IGFBP2", "WSB1", "DDIT4", "PGK1", "BNIP3", "FAM162A", "TPI1",
#'     "VEGFA", "PDK1", "PGAM1", "IER2", "FOS", "BTG1", "EPB41L4A-AS1", "NPAS4", "HK2", "BNIP3L",
#'     "JUN", "ENO2", "GAPDH", "ANKRD37", "ALDOA", "GADD45G", "TXNIP"
#'   )
#'   sobj <- addMetaFraction(col.name = "percent.HGA", gene.set = HGA_MarkerGenes, obj = sobj)
#' }
#' }
#' @seealso
#'  \code{\link[Matrix]{colSums}}
#' @export
#' @importFrom Matrix colSums
#' @importFrom CodeAndRoll2 grepv
addMetaFraction <- function(
    col.name = "percent.mito", gene.symbol.pattern = c("^MT\\.|^MT-", FALSE)[1],
    gene.set = FALSE, obj = ls.Seurat[[1]],
    verbose = TRUE) {
  message("Should rather use the default `Seurat::PercentageFeatureSet`")

  stopif(condition = isFALSE(gene.set) && isFALSE(gene.symbol.pattern), "Either gene.set OR gene.symbol.pattern has to be defined (!= FALSE).")
  if (!isFALSE(gene.set) && !isFALSE(gene.symbol.pattern) && verbose) print("Both gene.set AND gene.symbol.pattern are defined. Only using gene.set.")

  if (!isFALSE(gene.set)) geneset <- check.genes(list.of.genes = gene.set, obj = obj)
  total_expr <- Matrix::colSums(GetAssayData(object = obj))
  genes.matching <- if (!isFALSE(gene.set)) intersect(gene.set, rownames(obj)) else CodeAndRoll2::grepv(pattern = gene.symbol.pattern, x = rownames(obj))

  genes.expr <- GetAssayData(object = obj)[genes.matching, ]
  target_expr <- if (length(genes.matching) > 1) Matrix::colSums(genes.expr) else genes.expr

  iprint(length(genes.matching), "genes found, :", head(genes.matching))

  obj <- AddMetaData(object = obj, metadata = target_expr / total_expr, col.name = col.name)
  colnames(obj@meta.data)
  return(obj)
}


# _________________________________________________________________________________________________
#' @title add.meta.tags
#'
#' @description Add metadata tags to a Seurat object dataset.
#' @param list.of.tags A list of tags to be added as metadata. Default: tags
#' @param obj A Seurat object to which the metadata tags are to be added. Default: ls.Seurat[[1]]
#' @param n The index specifying the dataset for which the tags should be applied. Default: 1
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   ls.Seurat[[1]] <- add.meta.tags(list.of.tags = tags, obj = ls.Seurat[[1]], n = 1)
#' }
#' }
#' @export
add.meta.tags <- function(list.of.tags = tags, obj = ls.Seurat[[1]], n = 1) { # N is the for which dataset
  stopifnot(length(names(tags)) == length(tags))
  nCells <- nrow(obj@meta.data)
  for (i in 1:length(list.of.tags)) {
    tagX <- list.of.tags[[i]]
    new.meta.tag.per.cell <- rep(tagX[n], nCells)
    obj <- AddMetaData(object = obj, metadata = new.meta.tag.per.cell, col.name = names(tags)[i])
  }
  return(obj)
}


# _________________________________________________________________________________________________
#' @title seu.add.meta.from.table
#'
#' @description Add multiple new metadata columns to a Seurat object from a table. #
#' @param obj Seurat object, Default: seu.ORC
#' @param meta Metadata data frame.
#' @param suffix A suffix added to the filename, Default: '.fromMeta'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   combined.obj <- seu.add.meta.from.table()
#' }
#' }
#' @export
seu.add.meta.from.table <- function(obj = combined.obj, meta, suffix = ".fromMeta") { # Add multiple new metadata columns to a Seurat object from a table.
  NotFound <- setdiff(colnames(obj), rownames(meta))
  Found <- intersect(colnames(obj), rownames(meta))
  if (length(NotFound)) iprint(length(NotFound), "cells were not found in meta, e.g.: ", trail(NotFound, N = 10))

  mCols.new <- colnames(meta)
  mCols.old <- colnames(obj@meta.data)
  overlap <- intersect(mCols.new, mCols.old)
  if (length(overlap)) {
    iprint(length(overlap), "metadata columns already exist in the seurat object: ", overlap, ". These are tagged as: *", suffix)
    colnames(meta)[overlap] <- paste0(overlap, suffix)
  }
  mCols.add <- colnames(meta)
  obj@meta.data[Found, mCols.add] <- meta[Found, ]

  return(obj)
}

# _________________________________________________________________________________________________
#' @title seu.map.and.add.new.ident.to.meta
#'
#' @description Adds a new metadata column to a Seurat object based on an identity mapping table.
#' @param obj The Seurat object to which the new metadata column will be added. Default: combined.obj.
#' @param ident.table A data frame or matrix with identity mapping data. This parameter is used to map the old identities to the new ones. Default: clusterIDs.GO.process.
#' @param orig.ident The original identities of the Seurat object. Default: Idents(obj).
#' @param metaD.colname A string specifying the name of the new metadata column. The default value is the name of the provided ident.table.
#' @return A Seurat object with the new metadata column added.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Example usage:
#'   combined.obj <- seu.map.and.add.new.ident.to.meta(
#'     obj = combined.obj,
#'     ident.table = clusterIDs.GO.process
#'   )
#' }
#' }
#' @export
seu.map.and.add.new.ident.to.meta <- function(
    obj = combined.obj, ident.table = clusterIDs.GO.process,
    orig.ident = Idents(obj),
    metaD.colname = substitute(ident.table)) {
  # identities should match
  {
    Idents(obj) <- orig.ident
    ident.vec <- as.named.vector.df(ident.table)
    ident.X <- names(ident.vec)
    ident.Y <- as.character(ident.vec)
    ident.Seu <- gtools::mixedsort(levels(Idents(obj)))
    iprint("ident.Seu: ", ident.Seu)

    OnlyInIdentVec <- setdiff(ident.X, ident.Seu)
    OnlyInSeuratIdents <- setdiff(ident.Seu, ident.X)

    msg.IdentVec <- kollapse("Rownames of 'ident.table' have entries not found in 'Idents(obj)':",
      OnlyInIdentVec, " not found in ", ident.Seu,
      collapseby = " "
    )

    msg.Seu <- kollapse("Rownames of 'Idents(obj)' have entries not found in 'ident.table':",
      OnlyInSeuratIdents, " not found in ", ident.X,
      collapseby = " "
    )

    stopif(length(OnlyInIdentVec), message = msg.IdentVec)
    stopif(length(OnlyInSeuratIdents), message = msg.Seu)
  }
  # identity mapping
  {
    new.ident <- translate(vec = as.character(Idents(obj)), oldvalues = ident.X, newvalues = ident.Y)
    obj@meta.data[[metaD.colname]] <- new.ident
    iprint(metaD.colname, "contains the named identitites. Use Idents(combined.obj) = '...'. The names are:")
    cat(paste0("\t", ident.Y, "\n"))
  }
}




# _________________________________________________________________________________________________
# Replace / overwrite / remove metadata  ______________________________ ----
# _________________________________________________________________________________________________


# _________________________________________________________________________________________________
#' @title fix.orig.ident
#'
#' @description Remove the string "filtered_feature_bc_matrix." from "orig.ident". Helper function.
#' @param obj Seurat object, Default: merged.obj
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   merged.obj$orig.ident <- fix.orig.ident(obj = merged.obj)
#'   table(merged.obj$orig.ident)
#' }
#' }
#' @export
fix.orig.ident <- function(obj = merged.obj) {
  fixed <- sub(obj$"orig.ident", pattern = "filtered_feature_bc_matrix.", replacement = "")
  return(fixed)
}



# _________________________________________________________________________________________________
#' @title seu.RemoveMetadata
#'
#' @description Remove specified metadata columns from a Seurat object.
#' @param obj A Seurat object from which metadata columns will be removed. Default: combined.obj
#' @param cols_remove A character vector specifying metadata column names to remove. By default, it will remove all columns that do not start with "integr" or "cl.names".
#' @return A Seurat object with specified metadata columns removed.
#' @export
#' @examples
#' \dontrun{
#' combined.obj <- seu.RemoveMetadata(obj = combined.obj, cols_remove = c("column1", "column2"))
#' }
seu.RemoveMetadata <- function(
    obj = combined.obj,
    cols_remove = grepv(colnames(obj@meta.data), pattern = "^integr|^cl.names", perl = TRUE)) {
  CNN <- colnames(obj@meta.data)
  iprint("cols_remove:", cols_remove)
  print("")
  (cols_keep <- setdiff(CNN, cols_remove))
  obj@meta.data <- obj@meta.data[, cols_keep]
  iprint("meta.data colnames kept:", colnames(obj@meta.data))

  return(obj)
}


# _________________________________________________________________________________________________
# Export or Transfer metadata  ______________________________ ----
# _________________________________________________________________________________________________


#' @title Save Metadata from a List of Seurat Objects
#'
#' @description This function takes a list of Seurat objects, extracts their metadata, and saves it to a file with a specified suffix.
#'
#' @param ls.obj A list of Seurat objects.
#' @param suffix A character string to append to the filename when saving metadata.
#' @return Invisible list of metadata frames
#' @export
saveLsSeuratMetadata <- function(ls.obj, suffix) {
  stopifnot(is.list(ls.obj)) # Check if input is a list
  message(length(ls.obj), " objects")
  ls.meta <- setNames(lapply(ls.obj, function(x) x@meta.data), names(ls.obj))

  ncolz <- unique(unlapply(ls.meta, ncol))
  message(ncolz, " columns in meta.data")
  if (length(ncolz) > 1) warning("Different column counts across meta.data!", immediate. = TRUE)
  xsave(ls.meta, suffix = suffix)
  invisible(ls.meta)
}


# _________________________________________________________________________________________________
#' @title Transfer Multiple Metadata Columns Between Two Seurat Objects
#'
#' @description Transfers specified metadata columns from one Seurat object to another,
#' with options for verbose output and overwriting existing columns. Checks for cell overlap and
#' reports percentages of matching and unique cells.
#'
#' @param from The source Seurat object from which metadata will be transferred.
#' @param to The destination Seurat object to which metadata will be added.
#' @param colname_from Vector of names for the columns in the source object's metadata to transfer.
#' @param colname_to Vector of names for the columns in the destination object's metadata.
#' Defaults to the same names as `colname_from`. Must be the same length as `colname_from` unless
#' it is the same as `colname_from`.
#' @param verbose Logical, indicating whether to print details about the transfer, including the
#' number and percentage of matching cells between objects, and unique cells in each object.
#' @param overwrite Logical, indicating whether to overwrite the column in the destination object
#' if it already exists. Defaults to FALSE.
#'
#' @return Returns the destination Seurat object (`to`) with the new metadata columns added.
#'
#' @examples
#' # Assuming `object1` and `object2` are Seurat objects, and you want to transfer
#' # metadata columns named 'patientID' and 'treatmentGroup' from `object1` to `object2`:
#' object2 <- transferMetadata(
#'   from = object1, to = object2,
#'   colname_from = c("patientID", "treatmentGroup")
#' )
#'
#' @details This function is useful for merging related data from separate Seurat objects,
#' ensuring that relevant metadata is consistent across datasets. The function checks for
#' the existence of the specified columns in the source object and whether the columns
#' can be overwritten in the destination object. It also reports on cell overlap between
#' the two objects, which can be useful for understanding the relationship between datasets.
#'
#' @export
transferMetadata <- function(from, to, colname_from, colname_to = colname_from, verbose = TRUE, overwrite = FALSE) {
  stopifnot(
    is(from, "Seurat"), is(to, "Seurat"),
    is.character(colname_from), is.character(colname_to),
    all(colname_from %in% colnames(from@meta.data)),
    "Length of 'colname_from' and 'colname_to' must be equal" =
      length(colname_from) == length(colname_to)
  )

  # Check cell overlaps
  cells_in_both <- intersect(colnames(from), colnames(to))
  cells_only_in_from <- setdiff(colnames(from), colnames(to))
  cells_only_in_to <- setdiff(colnames(to), colnames(from))
  nr.cells.both <- length(cells_in_both)
  nr.cells.from <- length(cells_only_in_from)
  nr.cells.to <- length(cells_only_in_to)

  if (verbose) {
    if (verbose) {
      cat(
        "Cells matching between objects:", nr.cells.both,
        "(", sprintf("%.2f%%", nr.cells.both / length(colnames(from)) * 100), "of from and",
        sprintf("%.2f%%", nr.cells.both / length(colnames(to)) * 100), "of to)\n"
      )
      cat(
        "Cells only in obj1 (from):", length(cells_only_in_from),
        "(", sprintf("%.2f%%", nr.cells.from / length(colnames(from)) * 100), ")\n"
      )
      cat(
        "Cells only in obj2 (to):", nr.cells.to,
        "(", sprintf("%.2f%%", nr.cells.to / length(colnames(to)) * 100), ")\n"
      )
    }
  }

  for (i in seq_along(colname_from)) {
    if (!(colname_to[i] %in% colnames(to@meta.data)) || overwrite) {
      if (colname_from[i] %in% colnames(from@meta.data)) {
        # Transfer the metadata column
        to[[colname_to[i]]] <- from[[colname_from[i]]]
        message(sprintf("Transferred '%s' to '%s'.", colname_from[i], colname_to[i]))
      } else {
        warning(sprintf("Column '%s' not found in source object.", colname_from[i]), immediate. = TRUE)
      }
    } else {
      warning(sprintf(
        "Column '%s' already exists in destination object. Set 'overwrite = TRUE' to overwrite.",
        colname_to[i]
      ), immediate. = TRUE)
    }
  }
  return(to)
}



# _________________________________________________________________________________________________
# Subset metadata  ______________________________ ----
# _________________________________________________________________________________________________

# _________________________________________________________________________________________________
#' @title Sample N % of a dataframe (obj@metadata), and return rownames (cell IDs).
#'
#' @description This function samples a specified percentage of a dataframe (specifically a subset
#'   of the metadata of a Seurat object) and returns the corresponding cell IDs.
#' @param metaDF A dataframe representing a subset of the metadata of a Seurat object. Default:
#'   Subset of 'MetaData' for which 'Pass' is TRUE.
#' @param pc The percentage of the dataframe to sample, expressed as a decimal. Default: 0.1.
#' @return A vector of sampled cell IDs.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Example usage:
#'   # Suppose 'MetaData' is a dataframe and 'Pass' is a boolean vector with the same length.
#'   # The following example will sample 10% of the rows of 'MetaData' where 'Pass' is TRUE.
#'   sampleNpc(metaDF = MetaData[which(Pass), ], pc = 0.1)
#' }
#' }
#' @export
sampleNpc <- function(metaDF = MetaData[which(Pass), ], pc = 0.1) {
  cellIDs <- rownames(metaDF)
  nr_cells <- floor(length(cellIDs) * pc)
  cellIDs.keep <- sample(cellIDs, size = nr_cells, replace = FALSE)
  return(cellIDs.keep)
}



# _________________________________________________________________________________________________
# Combine metadata ______________________________ ----
# _________________________________________________________________________________________________


#' @title Combine Metadata from a list of Seurat objects and Write to TSV
#'
#' @description
#' Formerly `writeMetadataToTsv`. `writeCombinedMetadataToTsvFromLsObj` takes a list of ls.Obj, extracts their `@meta.data` slots,
#' removes specified columns, checks for column consistency, creates a barplot showing the number
#' of rows per object, and finally merges these into one large data frame.
#'
#' @param ls.Obj A list of objects, each containing a `@meta.data` slot.
#' @param cols.remove A character vector of column names to be removed from each metadata data frame.
#'        Default is an empty character vector, meaning no columns will be removed.
#'
#' @details
#' The function starts by validating the input to ensure it's a list. It then extracts the `@meta.data`
#' from each object, removing the specified columns. It checks if all data frames have the same columns
#' and issues a warning if not. A barplot is created to visualize the number of rows (cells) per object.
#' Finally, it merges all the metadata into one large data frame and prints its dimensions.
#'
#' @return A large data frame that is the row-wise merge of all `@meta.data` data frames.
#'
#' @examples
#' # Assuming a list of Seurat objects with meta.data
#' mergedMetaData <- writeMetadataToTsv(seuratObjectsList, cols.remove = c("column1", "column2"))
#'
#' @note
#' This function is intended for use with S4 objects that have a `@meta.data` slot.
#' The function currently contains a `browser()` call for debugging purposes, which should be removed in production.
#'
#' @export
writeCombinedMetadataToTsvFromLsObj <- function(ls.Obj, cols.remove = character(),
                                                save_as_qs = TRUE, save_as_tsv = TRUE, ...) {
  warning("writeMetadataToTsv is EXPERIMENTAL. It writes out subset of columns", immediate. = TRUE)
  stopifnot(is.list(ls.Obj)) # Validate that input is a list

  # Extract metadata from each object and remove specified columns
  metadataList <- lapply(ls.Obj, function(obj) {
    stopifnot("meta.data" %in% slotNames(obj)) # Check for meta.data slot
    metaData <- obj@meta.data
    metaData[, !(names(metaData) %in% cols.remove)]
  })

  # Find common columns and subset
  commonCols <- CodeAndRoll2::intersect.ls(lapply(metadataList, names))
  metadataList <- lapply(metadataList, function(df) df[, commonCols, drop = FALSE])


  # Check if qbarplot is available and create a barplot showing the number of rows per object
  metadata.cells.per.obj <- sapply(metadataList, nrow)
  print(metadata.cells.per.obj)
  pobj <- ggExpress::qbarplot(metadata.cells.per.obj,
    label = metadata.cells.per.obj, ylab = "cells",
    save = FALSE
  )
  print(pobj)

  # Merge metadata into one big data frame
  mergedMetaData <- do.call(rbind, metadataList)

  # Print dimensions of the merged data frame
  print(dim(mergedMetaData))

  if (save_as_qs) xsave(mergedMetaData)
  if (save_as_tsv) ReadWriter::write.simple.tsv(mergedMetaData, ...)

  # Return the merged data frame
  invisible(mergedMetaData)
}



# _________________________________________________________________________________________________
# Plot metadata ______________________________ ----
# _________________________________________________________________________________________________
#' @title Plot Metadata Correlation Heatmap
#'
#' @description This function plots a heatmap of metadata correlation values. It accepts a Seurat object
#' and a set of metadata columns to correlate. The correlations are calculated using either Pearson
#' or Spearman methods, and the resulting heatmap can include the principal component (PCA) values
#' and be saved with a specific suffix.
#'
#' @param columns A vector of metadata column names to calculate correlations.
#' Default: c("nCount_RNA", "nFeature_RNA", "percent.mito", "percent.ribo").
#' @param obj The main Seurat object used for calculations. No default value.
#' @param cormethod The method to calculate correlations. Can be either "pearson" or "spearman". Default: "pearson".
#' @param main The main title for the plot. Default: "Metadata correlations" followed by the correlation method.
#' @param show_numbers Logical, determines if correlation values should be displayed on the plot. Default: FALSE.
#' @param digits The number of decimal places for displayed correlation values. Default: 1.
#' @param suffix A suffix added to the output filename. Default: NULL.
#' @param add_PCA Logical, determines if PCA values should be included in the correlation calculation. Default: TRUE.
#' @param n_PCs The number of PCA components to be included if 'add_PCA' is TRUE. Default: 8.
#' @param w The width of the plot in inches. Default: ceiling((length(columns)+n_PCs)/2).
#' @param h The height of the plot in inches. Default: the value of w.
#' @param use_ggcorrplot Logical, determines if the ggcorrplot package should be used for plotting. Default: FALSE.
#' @param n_cutree The number of clusters to be used in hierarchical clustering. Default: the number of PCs.
#' @param ... Additional parameters passed to the internally called ggcorrplot function.
#'
#' @seealso
#'  \code{\link[ggcorrplot]{ggcorrplot}}
#' @importFrom ggcorrplot ggcorrplot
#' @importFrom pheatmap pheatmap
#' @export
plotMetadataCorHeatmap <- function(
    columns = c("nCount_RNA", "nFeature_RNA", "percent.mito", "percent.ribo"),
    obj,
    cormethod = c("pearson", "spearman")[1],
    main = paste("Metadata", cormethod, "correlations"),
    show_numbers = FALSE,
    digits = 1,
    suffix = NULL,
    add_PCA = TRUE,
    n_PCs = 8,
    w = ceiling((length(columns) + n_PCs) / 2), h = w,
    use_ggcorrplot = FALSE,
    n_cutree = (n_PCs),
    ...) {
  meta.data <- obj@meta.data
  columns.found <- intersect(colnames(meta.data), columns)
  columns.not.found <- setdiff(columns, colnames(meta.data))
  if (length(columns.not.found)) iprint("columns.not.found:", columns.not.found)

  meta.data <- meta.data[, columns.found]

  if (add_PCA) {
    stopif(is.null(obj@reductions$"pca"), "PCA not found in @reductions.")
    main <- paste("Metadata and PC", cormethod, "correlations")
    suffix <- FixPlotName(suffix, "w.PCA")

    PCs <- obj@reductions$pca@cell.embeddings
    stopifnot(nrow(meta.data) == nrow(PCs))
    meta.data <- cbind(PCs[, 1:n_PCs], meta.data)
  }

  corX <- cor(meta.data, method = cormethod)
  if (use_ggcorrplot) {
    pl <- ggcorrplot::ggcorrplot(corX,
      title = main,
      hc.order = TRUE,
      digits = digits,
      lab = show_numbers,
      type = "full",
      ...
    )
    ggExpress::qqSave(pl, fname = FixPlotName(make.names(main), suffix, "pdf"), w = w, h = h)
  } else {
    pl <- pheatmap::pheatmap(corX,
      main = main, treeheight_row = 2, treeheight_col = 2,
      cutree_rows = n_cutree, cutree_cols = n_cutree
    )
    wplot_save_pheatmap(
      x = pl, width = w, ,
      plotname = FixPlotName(make.names(main), suffix, "pdf")
    )
  }
  pl
}


# _________________________________________________________________________________________________
#' @title Calculate and plot heatmap of cluster medians
#'
#' @description This function calculates the median of specified variables in a dataframe,
#' grouped by a column ('ident'). The function also provides an option to scale the medians,
#' subset the ident levels, and either return a matrix of median values or plot a heatmap.
#'
#' @param meta A dataframe containing metadata from a Seurat object.
#' @param ident A character string representing the column name to be used for grouping the data.
#' @param subset_ident_levels An optional vector of ident levels to subset. Default is FALSE.
#' @param variables A character vector containing the names of columns for which to calculate the median.
#' @param scale A logical indicating whether to scale the medians. Default is TRUE.
#' @param suffix A character string added to the plot file name if not returning a matrix. Default is NULL.
#' @param return_matrix A logical indicating whether to return a matrix of medians, or plot a heatmap. Default is FALSE.
#' @param plotname A character string representing the main title for the plot. Default is "Median metadata values by cluster".
#' @param n_cut_row The number of row rows to cut the tree into on the `pheatmap`. Default is NA (none).
#' @param n_cut_col The number of column clusters to cut the tree into on the `pheatmap`. Default is NA (none)
#' @param w The width of the plot (if not returning a matrix). Default is half the number of variables.
#' @param ... Additional parameters passed to the pheatmap function.
#'
#' @return If 'return_matrix' is TRUE, a matrix where rows correspond to the unique values of 'ident',
#' and columns correspond to 'variables'. Each element of the matrix represents the median of a specific
#' variable for a specific group. If 'return_matrix' is FALSE, it saves a heatmap plot and returns the plot object.
#'
#' @importFrom dplyr group_by summarize_at
#' @importFrom ReadWriter FirstCol2RowNames
#' @importFrom pheatmap pheatmap
#' @import tidyverse
#' @export
heatmap_calc_clust_median <- function(
    meta, ident, subset_ident_levels = FALSE,
    variables, scale = TRUE,
    suffix = NULL,
    return_matrix = FALSE,
    plotname = "Median metadata values by cluster",
    n_cut_row = NA,
    n_cut_col = NA,
    w = ceiling(length(variables) / 2),
    ...) {
  # Ensure that 'meta' is a dataframe, 'ident' is a column in 'meta', and 'variables' are columns in 'meta'
  stopifnot(is.data.frame(meta))
  stopifnot(ident %in% colnames(meta))
  stopifnot(all(variables %in% colnames(meta)))

  # Group by 'ident' and calculate median for each variable
  df_cluster_medians <- meta %>%
    group_by(meta[[ident]]) %>%
    summarize_at(vars(variables), median, na.rm = TRUE)
  df_cluster_medians <- ReadWriter::FirstCol2RowNames(df_cluster_medians)

  if (!isFALSE(subset_ident_levels)) {
    stopifnot(all(subset_ident_levels %in% rownames(df_cluster_medians)))
    suffix <- FixPlotName(suffix, "subset")
    df_cluster_medians <- df_cluster_medians[subset_ident_levels, ]
  }

  df_cluster_medians <- scale(df_cluster_medians)

  if (return_matrix) {
    return(df_cluster_medians)
  } else {
    plot_name <- FixPlotName(plotname, suffix)
    pl <- pheatmap::pheatmap(df_cluster_medians,
      main = plot_name,
      cutree_rows = n_cut_row,
      cutree_cols = n_cut_col,
      ...
    )
    wplot_save_pheatmap(
      x = pl, width = w, ,
      plotname = FixPlotName(make.names(plot_name), suffix, "pdf")
    )
  }
}



# _________________________________________________________________________________________________
#' @title plotMetadataMedianFractionBarplot
#'
#' @description Generates a barplot of metadata median values.
#' @param columns A vector of column names to consider for the barplot. Default: c("percent.mito", "percent.ribo").
#' @param suffix A suffix added to the output filename. Default: NULL.
#' @param group.by The variable to group by for calculations. Default: Second result of GetClusteringRuns(obj).
#' @param method Method used for calculations, either "median" or "mean". Default: "median".
#' @param min.thr Minimum threshold percentage for a cluster. Default: 2.5.
#' @param return.matrix Logical; if TRUE, returns a matrix. Default: FALSE.
#' @param main Main title for the plot. Default: "read fractions per transcript class and cluster" followed by the method and suffix.
#' @param ylab Label for the y-axis. Default: "Fraction of transcriptome (%)".
#' @param percentify Logical; if TRUE, multiplies the fraction by 100. Default: TRUE.
#' @param subt Subtitle for the plot. Default: NULL.
#' @param position Position adjustment for geoms. Default: position_stack().
#' @param w The width of the plot. Default: 10.
#' @param h The height of the plot. Default: 6.
#' @param obj The main Seurat object used for calculations. Default: combined.obj.
#' @param ... Additional parameters passed to the internally called functions.
#' @seealso
#'  \code{\link[dplyr]{summarise_all}}
#'  \code{\link[reshape2]{melt}}
#' @export plotMetadataMedianFractionBarplot
#' @importFrom dplyr summarize_all
#' @importFrom reshape2 melt

plotMetadataMedianFractionBarplot <- function(
    columns = c("percent.mito", "percent.ribo"),
    suffix = NULL,
    group.by = GetClusteringRuns(obj = obj)[2],
    method = c("median", "mean")[1],
    min.thr = 2.5 # At least this many percent in at least 1 cluster
    , return.matrix = FALSE,
    main = paste(method, "read fractions per transcript class and cluster", suffix),
    ylab = "Fraction of transcriptome (%)",
    percentify = TRUE,
    subt = NULL,
    position = position_stack(),
    w = 10, h = 6,
    obj = combined.obj,
    ...) {
  meta.data <- obj@meta.data
  stopifnot(group.by %in% colnames(meta.data))
  columns.found <- intersect(colnames(meta.data), c(group.by, columns))

  (mat.cluster.medians1 <- meta.data[, columns.found] %>%
    group_by_at(group.by) %>%
    dplyr::summarize_all(median)
  )
  if (min.thr > 0) {
    pass.cols <- colMax(mat.cluster.medians1[, -1]) > (min.thr / 100)
    cols.OK <- which_names(pass.cols)
    cols.FAIL <- which_names(!pass.cols)
    subt <- paste(length(cols.FAIL), "classed do not reach", min.thr, "% :", kpps(cols.FAIL))
    iprint(subt)
    mat.cluster.medians1 <- mat.cluster.medians1[, c(group.by, cols.OK)]
  }


  mat.cluster.medians <- mat.cluster.medians1 %>%
    reshape2::melt(id.vars = c(group.by), value.name = "Fraction")


  if (percentify) mat.cluster.medians$"Fraction" <- 100 * mat.cluster.medians$"Fraction"

  pl <- ggbarplot(mat.cluster.medians,
    x = group.by, y = "Fraction", fill = "variable",
    position = position,
    title = main, subtitle = subt, ylab = ylab
  )
  ggExpress::qqSave(pl, fname = ppp(make.names(main), "pdf"), w = w, h = h)
  pl
  if (return.matrix) mat.cluster.medians1 else pl
}



# _________________________________________________________________________________________________
#' @title Plot Metadata Category Pie Chart
#'
#' @description Generates a pie chart visualizing the distribution of categories within a specified
#' metadata column of a Seurat object.
#'
#' @param metacol The metadata column to visualize.
#' @param plot_name Name of the plot to generate.
#' @param obj Seurat object containing the metadata. Default: `combined.obj`.
#' @param max.categs The maximum number of categories to display in the pie chart.
#' If the number of categories exceeds this value, an error is thrown.
#' @param both_pc_and_value If `TRUE`, labels on the pie chart will show both the percentage
#' and the count of each category. If `FALSE`, only the percentage is shown.
#' @param subtitle Optional subtitle for the pie chart.
#' @param ... Additional arguments to pass to the pie chart plotting function.
#'
#' @examples
#' \dontrun{
#' plotMetadataCategPie(
#'   metacol = "Singlet.status",
#'   plot_name = "Singlet Status Distribution",
#'   obj = combined.obj,
#'   max.categs = 20,
#'   both_pc_and_value = TRUE
#' )
#' }
#'
#' @return A pie chart visualizing the distribution of categories within the specified metadata column.
#' @export
plotMetadataCategPie <- function(
    metacol = "Singlet.status",
    plot_name = paste(metacol, "distribution"),
    obj = combined.obj, max.categs = 20, both_pc_and_value = TRUE,
    subtitle = NULL, ...) {
  categ_pivot <- table(obj[[metacol]])
  stopifnot(length(categ_pivot) < max.categs)
  qpie(categ_pivot,
    plotname = FixPlotName(make.names(plot_name)),
    both_pc_and_value = both_pc_and_value,
    LegendSide = FALSE, labels = NULL, LegendTitle = "", subtitle = subtitle, ...
  )
}



# _________________________________________________________________________________________________
# Label / identity transfer across objects ______________________________ ----
# _________________________________________________________________________________________________




#' @title Rename Azimuth Columns in Seurat Object
#'
#' @description Dynamically renames specified metadata columns in a Seurat object, particularly those
#' prefixed with "predicted." and the "mapping.score" column, by applying a new prefix
#' that combines a user-defined prefix and a reference name.
#'
#' @param obj A Seurat object containing metadata in `meta.data` that needs column names to be renamed.
#' @param ref A character string specifying the reference; defaults to "humancortexref".
#'            The "ref" part of the string will be removed in the new column names.
#' @param prefix A character string to prefix to the column names, defaulting to "azi".
#'               This prefix is combined with the modified `ref` to form the new column names.
#' @param azim_cols Azimuth columns
#' @return Returns the Seurat object with renamed metadata columns.
#'
#' @examples
#' # Assuming `obj` is a Seurat object with metadata columns following the "predicted." pattern:
#' obj <- renameAzimuthColumns(obj, ref = "humancortexref", prefix = "azi")
#' # This will rename columns like "predicted.class" to "azi.humancortex.class"
#' # and include "mapping.score" as "azi.humancortex.mapping.score"
#'
#' @export
renameAzimuthColumns <- function(obj, ref = c("humancortexref", "fetusref")[1],
                                 prefix = "azi",
                                 azim_cols = CodeAndRoll2::grepv(
                                   x = tail(colnames(obj@meta.data), 10),
                                   pattern = "predicted."
                                 )) {
  stopifnot(
    "obj must be a Seurat object" = is(obj, "Seurat"),
    "azim_cols must be non-empty" = length(azim_cols) > 0
  )

  ref <- sub(pattern = "ref", replacement = "", x = ref)
  iprint(length(azim_cols), "azim_cols:", azim_cols)

  # Extract the column names of meta.data
  meta_col_names <- colnames(obj@meta.data)

  # Loop through the azim_cols and replace the prefix if they exist in meta.data
  for (azim_col in azim_cols) {
    if (azim_col %in% meta_col_names) {
      # Create the new column name by replacing "predicted." with the new prefix
      new_col_name <- sub(pattern = "^predicted\\.", replacement = kpp(prefix, ref, ""), x = azim_col)
      names(obj@meta.data)[names(obj@meta.data) == azim_col] <- new_col_name
    }
  }

  if ("mapping.score" %in% colnames(obj@meta.data)) {
    names(obj@meta.data)[names(obj@meta.data) == "mapping.score"] <- kpp(prefix, ref, "mapping.score")
  }

  print(tail(colnames(obj@meta.data), 10))
  return(obj)
}


# _________________________________________________________________________________________________
#' @title Rename Small Categories in Seurat Object Metadata
#'
#' @description This function renames categories within a specified identity column of a
#' Seurat object's metadata that have fewer cells than a specified minimum threshold.
#' Categories below this threshold are renamed to a common name, typically "unclear",
#' to clean up small, potentially noisy categories.
#' @param obj A Seurat object containing the metadata with categories to be cleaned.
#' @param idents A character vector specifying the names of the identity columns within
#'   `obj@meta.data` where categories are to be renamed.
#' @param min.cells An integer specifying the minimum number of cells a category must have
#'   to retain its original name. Categories with fewer cells than this threshold will be
#'   renamed. Defaults to the greater of the total number of columns divided by 2000 or 10.
#' @param new.name A character string specifying the new name to assign to small categories.
#'   Defaults to "unclear".
#'
#' @return Returns the Seurat object with renamed categories in the specified metadata columns.
#'
#' @examples
#' # Assuming obj is a Seurat object with identity columns "ident1" and "ident2":
#' idents <- c("ident1", "ident2")
#' obj <- renameSmallCategories(obj, idents = idents)
#'
#' @export
renameSmallCategories <- function(
    obj,
    idents = c("predicted.class", "predicted.cluster", "predicted.subclass"),
    min.cells = max(round((ncol(obj)) / 2000), 10),
    new.name = "unclear") {
  stopifnot("obj must be a Seurat object" = is(obj, "Seurat"))

  for (ident in idents) {
    if (ident %in% colnames(obj@meta.data)) {
      # Count the number of cells per category in the specified identity column
      category_counts <- table(obj@meta.data[[ident]])

      # Identify categories with fewer cells than min.cells
      small_categories <- names(category_counts[category_counts < min.cells])

      # Initial number of categories
      initial_categories <- length(unique(obj@meta.data[[ident]]))

      # Rename the categories in the ident column that have fewer cells than min.cells to new.name
      obj@meta.data[[ident]] <- ifelse(obj@meta.data[[ident]] %in% small_categories, new.name, obj@meta.data[[ident]])

      # Report to console
      cells_renamed <- sum(obj@meta.data[[ident]] == new.name)
      categories_removed <- length(small_categories)
      remaining_categories <- length(unique(obj@meta.data[[ident]]))

      message(
        "For ident '", ident, "':\n",
        cells_renamed, " cells were renamed.\n",
        remaining_categories, " of initial ", initial_categories, " categories remained.\n\n",
        "Removed categories: ", paste(head(small_categories, 10), collapse = ", "), "\n"
      )
    } else {
      message("Ident column '", ident, "' does not exist in obj@meta.data.")
    }
  }

  return(obj)
}



# _________________________________________________________________________________________________
#' @title Transfer labels from a reference Seurat object to a query Seurat object
#'
#' @description Function to transfer labels from a reference Seurat object to a query Seurat object
#' using anchoring and transfer data methods from the Seurat package. It then visualizes the
#' reference and the combined objects using Uniform Manifold Approximation and Projection (UMAP).
#'
#' @param query_obj A Seurat object for which the labels are to be transferred.
#' @param reference_path A character string indicating the file path to the reference Seurat object. The path must exist.
#' @param reference_obj Alternative to `reference_path`. If provided, the path is not used to load the reference data.
#' @param named_ident A character string specifying the name of the identity class to be used from the reference Seurat object. Default is 'RNA_snn_res.0.3.ordered.ManualNames'.
#' @param new_ident A character string specifying the name of the new identity class to be created in the query Seurat object. Default is obtained by replacing 'ordered' with 'transferred' in named_ident.
#' @param predictions_col A character string specifying the column in the metadata of the transferred Seurat object containing the transferred labels. Default is 'predicted.id'.
#' @param save_anchors save anchors as RDS file.
#' @param suffix A character string to be used as a suffix in the visualization. Default is 'NEW'.
#' @param plot_suffix A string to added to the UMAP with the new identity.
#' @param h Height for the saved image. Default: 12
#' @param w Width for the saved image. Default: 9
#' @param ... Additional arguments passed to the Seurat.utils::clUMAP function.
#'
#' @return The modified query Seurat object with the transferred labels as a new identity class.
#'
#' @examples
#' # combined.objX <- transferLabelsSeurat(named_ident = 'RNA_snn_res.0.3.ordered.ManualNames',
#' #                                     reference_obj = reference_obj,
#' #                                     query_obj = combined.obj)
#'
#' @importFrom readr read_rds
#' @importFrom Seurat FindTransferAnchors TransferData AddMetaData
#'
#' @export
transferLabelsSeurat <- function(
    query_obj,
    reference_obj,
    reference_path = NULL,
    reference_ident,
    anchors = NULL,
    new_ident = gsub(
      pattern = "ordered",
      replacement = "transferred",
      x = reference_ident
    ),
    predictions_col = "predicted.id",
    predictions_score = sppp(new_ident, "score"),
    save_anchors = TRUE,
    reference_suffix = "reference",
    plot_suffix = NULL,
    plot_reference = TRUE,
    w = 12, h = 9,
    ...) {
  # Assertions
  if (is.null(reference_obj)) {
    iprint("Loading reference object:", basename(reference_path))
    stopifnot(file.exists(reference_path))
    reference_obj <- readr::read_rds(reference_path)
  } else {
    stopifnot(inherits(reference_obj, "Seurat") & min(dim(reference_obj)) > 10)
  }

  # Report
  nr.cl.ref <- CodeAndRoll2::nr.unique(reference_obj[[reference_ident]])
  message("reference_ident ", reference_ident, " has ", nr.cl.ref, " categories")

  # Visualize reference object
  if (plot_reference) {
    clUMAP(
      obj = reference_obj, ident = reference_ident,
      suffix = reference_suffix, sub = reference_suffix,
      w = w, h = h, ...
    )
  }

  # browser()
  if (is.null(anchors)) {
    message("Calculating anchors. Provide anchors in 'anchors' to speed up.")
    anchors <- Seurat::FindTransferAnchors(reference = reference_obj, query = query_obj)
    if (save_anchors) xsave(obj = anchors)
  } else {
    message("Anchors provided")
  }

  message("Transferring labels")
  transferred_clIDs <- Seurat::TransferData(
    anchorset = anchors,
    refdata = reference_obj@meta.data[, reference_ident],
  )

  # browser()
  # Add New Labels to query object
  query_obj <- Seurat::AddMetaData(
    object = query_obj, metadata = transferred_clIDs[, predictions_col],
    col.name = new_ident
  )

  # Add Labels assignment scores to query object
  query_obj <- Seurat::AddMetaData(
    object = query_obj, metadata = transferred_clIDs[, "prediction.score.max"],
    col.name = predictions_score
  )

  # Visualize combined object
  clUMAP(
    ident = new_ident, obj = query_obj, suffix = plot_suffix,
    w = w, h = h, ...
  )
  qUMAP(
    feature = predictions_score, obj = query_obj, suffix = plot_suffix,
    w = w, h = h, ...
  )

  return(query_obj)
}


# _________________________________________________________________________________________________
#' @title Extract meta.data Column Names Matching a Pattern
#'
#' @param obj A dataframe from which to extract column names.
#' @param pattern A regular expression pattern to match column names against.
#'
#' @return A character vector of column names matching the pattern.
#'
#' @examples
#' # Assuming 'df' is a dataframe with column names "azi.one", "azi.two", "other"
#' extract_matching_columns(df, "^azi\\.")
.metaColnames <- function(obj = combined.obj, pattern, perl = TRUE, ...) {
  colz <- grep(pattern, colnames(obj@meta.data), value = TRUE, perl = perl, ...)
  dput(colz)
  return(colz)
}


# _________________________________________________________________________________________________
#' @title Match and Translate Best Identity
#'
#' @description This function matches the best identity from `ident_to_rename` to `reference_ident` in an object,
#' in other words, it replaces original categories with the most frequent ones from the reference,
#' hence helps to filter out less important categories.
#'
#' @param obj The object to update. This object must have a `meta.data` attribute which is a data frame
#'   containing columns named as `ident_to_rename` and `reference_ident`.
#' @param ident_to_rename A string. The name of the column in `obj@meta.data` that is used as the source of identities.
#'   There is no default value for this parameter.
#' @param reference_ident A string. The name of the column in `obj@meta.data` that is used as the target of identities.
#'   There is no default value for this parameter.
#' @param prefix A string to add to the new identity column name. Default is prefix = Reference.
#' @param suffix ...
#' @param new_ident_name A string. The name for the newly created identity column in `obj@meta.data`.
#'   Default is a concatenation: kpp(prefix, ident_to_rename, "match.to", reference_ident) .
#' @param plot_suffix A string. The suffix to add to the final UMAP.
#' @param h Height for the saved image. Default: 12
#' @param w Width for the saved image. Default: 9
#' @param ... Additional parameters to be passed to `.replace_by_most_frequent_categories` function.
#'
#' @return An updated version of `obj` with an additional column in `obj@meta.data` named as `new_ident_name`
#'   representing the new identity. The function also generates a UMAP plot based on this new identity.
#'
#' @seealso \code{\link[clUMAP]{clUMAP}}, \code{\link[kpp]{kpp}}, \code{\link[FixPlotName]{FixPlotName}},
#'   \code{\link[.replace_by_most_frequent_categories]{.replace_by_most_frequent_categories}}
#'
#' @examples
#' \dontrun{
#' updated_obj <- matchBestIdentity(my_obj, "origin_identity", "target_identity")
#' }
#' @export
matchBestIdentity <- function(
    obj, ident_to_rename,
    reference_ident = GetOrderedClusteringRuns(obj)[1],
    prefix = Reference,
    suffix = gsub(prefix, "", x = reference_ident),
    # to_suffix = "matched",
    # to_suffix = FixPlotName(gsub(pattern = "[a-zA-Z_]", replacement = "", x = ident_to_rename)),
    new_ident_name = kpp(prefix, ident_to_rename, "match.to", suffix),
    plot_suffix = prefix,
    w = 12, h = 9,
    ...) {
  stopifnot("colname prefix undefined" = !is.null(prefix))

  dictionary <- obj@meta.data[, c(ident_to_rename, reference_ident)]

  translation <- .replace_by_most_frequent_categories(
    df = dictionary, show_plot = TRUE, suffix_barplot = ident_to_rename, ...
  )

  obj@meta.data[, new_ident_name] <- translation[, 1]

  imessage("new ident name:", new_ident_name)
  px <- clUMAP(ident = new_ident_name, obj = obj, suffix = plot_suffix, w = w, h = h, ...)
  print(px)
  return(obj)
}



# _________________________________________________________________________________________________
#' @title Find Best Match: Replace Categories by the Most Frequent Match
#'
#' @description Used for mapping identity columns across objects. This function replaces each
#'   category in a query column of a data frame with the most frequently corresponding category in a
#'   reference column. It calculates the assignment quality, reports it, and optionally plots it.
#' @param df A data frame containing the data.
#' @param query_col The name of the column in 'df' whose categories are to be replaced. By default,
#'   the first column of 'df' is used.
#' @param ref_col The name of the column in 'df' used as reference for replacement. By default, the
#'   second column of 'df' is used.
#' @param show_plot Logical, whether to plot assignment quality. Defaults to TRUE.
#' @param suffix_barplot Suffix for barplot.
#' @param ... Additional parameters passed to the qbarplot function.
#' @return A data frame with categories in 'query_col' replaced by the most frequent match from
#'   'ref_col'.
#'
#' @importFrom dplyr group_by summarise arrange filter
#' @examples
#' \dontrun{
#' .replace_by_most_frequent_categories(df = my_data)
#' (MXX <- as.tibble(structure(
#'   c(
#'     "Adjut", "Adjut", "Yearn", "Adjut", "Dwarf", "Adjut",
#'     "Dwarf", "Adjut", "Dwarf", "Yearn", "Dwarf", "Dwarf", "Dwarf",
#'     "Yearn", "Dwarf", "Dwarf", "Dwarf", "Zebra", "Yucca", "Plyer",
#'     "Blaze", "Blaze", "Dazed", "Blaze", "Swept", "Bold", "Vixen",
#'     "Bold", "Swept", "Dazed", "Mirth", "Witch", "Vixen", "Dazed",
#'     "Swept", "Mirth", "Swept", "Vexed", "Query", "Yolk"
#'   ),
#'   .Dim = c(20L, 2L), .Dimnames =
#'     list(NULL, c("RNA_snn_res.0.1.ordered", "RNA_snn_res.0.3.ordered"))
#' )))
#'
#' z <- .replace_by_most_frequent_categories(df = MXX)
#' head(cbind(MXX[, 1], z[, 1]))
#' }
.replace_by_most_frequent_categories <- function(
    df, query_col = colnames(df)[1],
    ref_col = colnames(df)[2],
    show_plot = TRUE,
    suffix_barplot = NULL,
    ext = "png",
    ...) {
  # Convert to data frame if it is not
  if (!is.data.frame(df)) {
    df <- as.data.frame(df)
  }

  cat_query <- unique(df[[query_col]])
  imessage(length(cat_query), "categories to rename in", query_col, ":", head(cat_query), "...")

  cat_ref <- unique(df[[ref_col]])
  imessage(length(cat_ref), "reference categories in", ref_col, ":", head(cat_ref), "...")

  # Create a table of the most frequent reference values for each query category
  replacement_table <- df %>%
    dplyr::group_by(!!sym(query_col), !!sym(ref_col)) %>%
    dplyr::summarise(n = n(), .groups = "drop") %>%
    dplyr::arrange(!!sym(query_col), desc(n)) %>%
    dplyr::filter(!duplicated(!!sym(query_col)))

  replacement_table[[ref_col]] <- make.unique(replacement_table[[ref_col]])

  # Calculate assignment quality
  total_counts <- table(df[[query_col]])
  quality <- replacement_table$n / total_counts[replacement_table[[query_col]]]
  names(quality) <- paste0(names(quality), "->", replacement_table[[ref_col]])

  # Report assignment quality
  message("Assignment quality (proportion of total matches):")
  # print(setNames(quality, replacement_table[[query_col]]))

  # Plot assignment quality
  if (show_plot) {
    px <- ggExpress::qbarplot(quality,
      label = percentage_formatter(quality, digitz = 1),
      ext = ext,
      suffix = suffix_barplot,
      plotname = "Assignment Quality",
      filename = make.names(kpp("Assignment Quality", suffix_barplot, ext)),
      subtitle = paste(
        "From", colnames(df)[1], "->", colnames(df)[2], "| median",
        percentage_formatter(median(quality)), "\n",
        sum(quality > 0.5), "clusters above 50% match"
      ),
      hline = 0.5, filtercol = -1,
      xlab = paste("Best query match to reference"),
      ylab = "Proportion of Total Matches",
      ...
    )
    print(px)
  }

  # Replace the query values with the most frequent reference values
  df[[query_col]] <- replacement_table[[ref_col]][match(df[[query_col]], replacement_table[[query_col]])]

  return(df)
}
