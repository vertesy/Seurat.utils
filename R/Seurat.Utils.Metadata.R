# ____________________________________________________________________
# Seurat.utils ----
# ____________________________________________________________________
# source("/Users/abel.vertesy/GitHub/Packages/Seurat.utils/R/Seurat.Utils.Metadata.R")




# _________________________________________________________________________________________________
# metadata.manipulation.R ______________________________ ----
# _________________________________________________________________________________________________
# source('~/GitHub/Packages/Seurat.utils/Functions/Seurat.object.manipulations.etc.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Metadata.manipulation.R"))
# Source: self + web

# _________________________________________________________________________________________________
#' @title Check if a Column Exists in the Metadata of an S4 Object
#'
#' @description This function checks whether a given column exists in the meta.data of a Seurat object.
#' @param obj A Seurat object.
#' @param col_name A character string specifying the name of the column.
#'
#' @return A logical value indicating whether the column exists (TRUE) or not (FALSE).
#' @export
meta_col_exists <- function(col_name, obj) {
  col_name %in% colnames(obj@meta.data)
}


# _________________________________________________________________________________________________
#' @title getMedianMetric
#'
#' @description Get the median values of different columns in meta.data, can iterate over a list of Seurat objects.
#' @param ls.obj List of Seurat objects, Default: ls.Seurat
#' @param n.datasets lenght of list (n objects), Default: length(ls.Seurat)
#' @param mColname Metadata column name to calculate on. Default: 'percent.mito'
#' @examples
#' \dontrun{
#' if(interactive()){
#'  ls.Seurat <- getMedianMetric(ls.obj = ls.Seurat, n.datasets = length(ls.Seurat),
#'  mColname = "percent.mito")
#'  }
#' }
#' @export
getMedianMetric <- function(ls.obj = ls.Seurat, n.datasets = length(ls.Seurat), mColname = "percent.mito") {
  medMetric <- vec.fromNames(names(ls.obj))
  for(i in 1:n.datasets ) {
    medMetric[i] <- median(ls.obj[[i]]@meta.data[,mColname])
  }
  return(medMetric)
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
#' if(interactive()){
#'  ls.Seurat[[1]] <- add.meta.tags(list.of.tags = tags, obj = ls.Seurat[[1]], n = 1)
#'  }
#' }
#' @export
add.meta.tags <- function(list.of.tags = tags, obj = ls.Seurat[[1]], n = 1) {  # N is the for which dataset
  stopifnot( length(names(tags)) == length(tags) )
  nCells = nrow(obj@meta.data)
  for (i in 1:length(list.of.tags)) {
    tagX <- list.of.tags[[i]]
    new.meta.tag.per.cell <- rep(tagX[n], nCells)
    obj <- AddMetaData(object = obj, metadata = new.meta.tag.per.cell, col.name = names(tags)[i])
  }
  return(obj)
}



# _________________________________________________________________________________________________
#' @title add.meta.fraction
#'
#' @description Add a new metadata column to a Seurat object, representing the fraction of a gene set in the transcriptome (expressed as a percentage).
#' @param col.name Name of the new metadata column to be added. Default: 'percent.mito'
#' @param gene.symbol.pattern Regular expression pattern to match gene symbols. Default: c("^MT\\.|^MT-", F)[1]
#' @param gene.set A set of gene symbols. If specified, it will be used instead of gene.symbol.pattern. Default: F
#' @param obj Seurat object to which the new metadata column will be added. Default: ls.Seurat[[1]]
#' @param verbose Logical indicating whether to display detailed messages (TRUE) or not (FALSE). Default: T
#' @examples
#' \dontrun{
#' if(interactive()){
#'  ls.Seurat[[1]] <- add.meta.fraction(col.name = "percent.mito", gene.symbol.pattern = "^MT\\.|^MT-")
#'  ls.Seurat[[1]] <- add.meta.fraction(col.name = "percent.ribo", gene.symbol.pattern = "^RPL|^RPS")
#'  ls.Seurat[[1]] <- add.meta.fraction(col.name = "percent.AC.GenBank", gene.symbol.pattern = "^AC[0-9]{6}\\.")
#'  ls.Seurat[[1]] <- add.meta.fraction(col.name = "percent.AL.EMBL", gene.symbol.pattern = "^AL[0-9]{6}\\.")
#'  ls.Seurat[[1]] <- add.meta.fraction(col.name = "percent.LINC", gene.symbol.pattern = "^LINC0")
#'  ls.Seurat[[1]] <- add.meta.fraction(col.name = "percent.MALAT1", gene.symbol.pattern = "^MALAT1")
#'  colnames(ls.Seurat[[1]]@meta.data)
#'  HGA_MarkerGenes <- c("ENO1", "IGFBP2", "WSB1", "DDIT4", "PGK1", "BNIP3", "FAM162A", "TPI1",
#'  "VEGFA", "PDK1", "PGAM1", "IER2", "FOS", "BTG1", "EPB41L4A-AS1","NPAS4", "HK2", "BNIP3L",
#'  "JUN", "ENO2", "GAPDH", "ANKRD37", "ALDOA", "GADD45G", "TXNIP")
#'  sobj <- add.meta.fraction(col.name = "percent.HGA", gene.set = HGA_MarkerGenes, obj =  sobj)
#'  }
#' }
#' @seealso
#'  \code{\link[Matrix]{colSums}}
#' @export
#' @importFrom Matrix colSums
add.meta.fraction <- function(col.name = "percent.mito", gene.symbol.pattern = c("^MT\\.|^MT-", F)[1]
                              , gene.set = F, obj = ls.Seurat[[1]], verbose = T) {
  stopif2(condition = isFALSE(gene.set) && isFALSE(gene.symbol.pattern), "Either gene.set OR gene.symbol.pattern has to be defined (!= FALSE).")
  if (!isFALSE(gene.set) && !isFALSE(gene.symbol.pattern) && verbose) print("Both gene.set AND gene.symbol.pattern are defined. Only using gene.set.")

  if (!isFALSE(gene.set)) geneset <- check.genes(list.of.genes = gene.set, obj = obj)
  total_expr <- Matrix::colSums(GetAssayData(object = obj))
  genes.matching <- if (!isFALSE(gene.set)) intersect(gene.set, rownames(obj)) else CodeAndRoll2::grepv(pattern = gene.symbol.pattern, x = rownames(obj))

  genes.expr = GetAssayData(object = obj)[genes.matching, ]
  target_expr <- if (length(genes.matching) >1) Matrix::colSums(genes.expr) else genes.expr

  iprint(length(genes.matching), "genes found, :", head(genes.matching))

  obj <- AddMetaData(object = obj, metadata = target_expr / total_expr, col.name = col.name)
  colnames(obj@meta.data)
  return(obj)
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
#'   combined.obj <- seu.RemoveMetadata(obj = combined.obj, cols_remove = c("column1", "column2"))
#' }
seu.RemoveMetadata <- function(obj = combined.obj
                               , cols_remove = grepv(colnames(obj@meta.data), pattern = "^integr|^cl.names", perl = T)
) {

  CNN <- colnames(obj@meta.data)
  iprint('cols_remove:', cols_remove); print('')
  (cols_keep <- setdiff(CNN, cols_remove))
  obj@meta.data <- obj@meta.data[, cols_keep]
  iprint('meta.data colnames kept:', colnames(obj@meta.data))

  return(obj)
}



# _________________________________________________________________________________________________
#' @title getMetadataColumn
#'
#' @description Retrieves a specified metadata column from a Seurat object and returns it as a named vector.
#' @param ColName.metadata A string specifying the name of the metadata column to be retrieved. Default: 'batch'.
#' @param obj A Seurat object from which the metadata column will be retrieved. Default: combined.obj.
#' @param as_numeric A logical flag indicating whether the returned values should be converted to numeric format. Default: F (FALSE).
#' @return A named vector containing the values from the specified metadata column. If 'as_numeric' is TRUE, the values are converted to numeric format.
#' @examples
#' \dontrun{
#' if(interactive()){
#'   # Example usage:
#'   batch_metadata <- getMetadataColumn(ColName.metadata = 'batch', obj = combined.obj, as_numeric = T)
#'  }
#' }
#' @export
getMetadataColumn <- mmeta <- function(ColName.metadata = 'batch', obj = combined.obj, as_numeric =F) { # Get a metadata column from a Seurat object as a named vector
  stopifnot(ColName.metadata %in% colnames(obj@meta.data))

  x = as.named.vector.df(obj@meta.data[ ,ColName.metadata, drop = F])
  if (as_numeric) {
    as.numeric.wNames(x)+1
  } else {x}
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
#' if(interactive()){
#'   # Example usage:
#'   combined.obj <- seu.add.meta.from.vector(obj = combined.obj,
#'                                            metaD.colname = metaD.colname.labeled,
#'                                            Label.per.cell = Cl.Label.per.cell)
#'  }
#' }
#' @export
seu.add.meta.from.vector <- function(obj = combined.obj, metaD.colname = metaD.colname.labeled, Label.per.cell = Cl.Label.per.cell ) { # Add a new metadata column to a Seurat  object
  obj@meta.data[, metaD.colname ] = Label.per.cell
  iprint(metaD.colname, "contains the named identitites. Use Idents(combined.obj) = '...'. The names are:", unique(Label.per.cell))
  return(obj)
}



# _________________________________________________________________________________________________
#' create.metadata.vector
#'
#' @param vec cell_ID vector with partial overlap to cells in a Seurat object.
#' @param obj Seurat object
#' @param min.intersect Min number of cells to find in both.
#'
#' @export

create.metadata.vector <- function(vec = All.UVI, obj = combined.obj, min.intersect = 100) {
  cells.vec <- names(vec)
  cells.obj <- colnames(obj)
  cells.in.both <- intersect(cells.vec, cells.obj)

  # iprint("intersect:", l(cells.in.both), head(cells.in.both))
  iprint(l(cells.in.both), 'cells in both;'
         , l(cells.vec), 'cells in vec;'
         , l(cells.obj), 'cells in obj'
         , "intersect, e.g.:", head(cells.in.both, 5))
  stopifnot(length(cells.in.both) > min.intersect )

  new_assignment <- vec.fromNames(cells.obj)
  new_assignment[cells.in.both] <- vec[cells.in.both]
  return(new_assignment)
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
#' if(interactive()){
#'   # Example usage:
#'   combined.obj <- seu.map.and.add.new.ident.to.meta(obj = combined.obj,
#'                                                     ident.table = clusterIDs.GO.process)
#'  }
#' }
#' @export
seu.map.and.add.new.ident.to.meta <- function(obj = combined.obj, ident.table = clusterIDs.GO.process, orig.ident = Idents(obj)
                                              , metaD.colname = substitute(ident.table) ) { # Add a new metadata column to a Seurat  object
  # identities should match
  {
    Idents(obj) <- orig.ident
    ident.vec <- as.named.vector.df(ident.table)
    ident.X <- names(ident.vec)
    ident.Y <- as.character(ident.vec)
    ident.Seu <- sort.natural(levels(Idents(obj)))
    iprint("ident.Seu: ", ident.Seu)

    OnlyInIdentVec      <- setdiff(ident.X, ident.Seu)
    OnlyInSeuratIdents  <- setdiff(ident.Seu, ident.X)

    msg.IdentVec <- kollapse("Rownames of 'ident.table' have entries not found in 'Idents(obj)':"
                             , OnlyInIdentVec, " not found in ", ident.Seu, collapseby = " ")

    msg.Seu <- kollapse("Rownames of 'Idents(obj)' have entries not found in 'ident.table':"
                        , OnlyInSeuratIdents, " not found in ", ident.X, collapseby = " ")

    stopif (length(OnlyInIdentVec), message = msg.IdentVec)
    stopif (length(OnlyInSeuratIdents), message = msg.Seu)
  }
  # identity mapping
  {
    new.ident <- translate(vec = as.character(Idents(obj)), oldvalues = ident.X, newvalues = ident.Y)
    obj@meta.data[[metaD.colname]] = new.ident
    iprint(metaD.colname, "contains the named identitites. Use Idents(combined.obj) = '...'. The names are:"); cat(paste0("\t", ident.Y, "\n"))
  }
}




# _________________________________________________________________________________________________
#' @title getCellIDs.from.meta
#'
#' @description Retrieves cell IDs from a specified metadata column of a Seurat object, where the cell ID matches a provided list of values. The matching operation uses the `%in%` operator.
#' @param ColName.meta A string specifying the name of the metadata column from which to retrieve cell IDs. Default: 'res.0.6'.
#' @param values A vector of values to match in the metadata column. Default: NA.
#' @param obj The Seurat object from which to retrieve the cell IDs. Default: combined.obj.
#' @param inverse A boolean value indicating whether to inverse the match, i.e., retrieve cell IDs that do not match the provided list of values. Default: F.
#' @return A vector of cell IDs that match (or don't match, if `inverse = TRUE`) the provided list of values.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  # Example usage:
#'  getCellIDs.from.meta()
#'  }
#' }
#' @export
getCellIDs.from.meta <- function(ColName.meta = 'res.0.6', values = NA, obj = combined.obj, inverse = F ) { # Get cellIDs from a metadata column, matching a list of values (using %in%).
  mdat <- obj@meta.data[ , ColName.meta]
  cells <- if (inverse) {mdat %!in% values} else {mdat %in% values}
  idx.matching.cells = which(cells)
  iprint(length(idx.matching.cells), 'cells found.')
  return(rownames(obj@meta.data)[idx.matching.cells])
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
#' if(interactive()){
#'  combined.obj <- seu.add.meta.from.table()
#'  }
#' }
#' @export
seu.add.meta.from.table <- function(obj = combined.obj, meta, suffix = ".fromMeta") { # Add multiple new metadata columns to a Seurat object from a table.
  NotFound  = setdiff(colnames(obj), rownames(meta))
  Found     = intersect(colnames(obj), rownames(meta))
  if (length(NotFound)) iprint(length(NotFound), 'cells were not found in meta, e.g.: ', trail(NotFound, N = 10))

  mCols.new = colnames(meta)
  mCols.old = colnames(obj@meta.data)
  overlap = intersect(mCols.new, mCols.old)
  if (length(overlap)) {
    iprint(length(overlap), 'metadata columns already exist in the seurat object: ', overlap, '. These are tagged as: *', suffix)
    colnames(meta)[overlap] = paste0(overlap, suffix)
  }
  mCols.add = colnames(meta)
  obj@meta.data[Found, mCols.add] = meta[ Found,]

  return(obj)
}


# _________________________________________________________________________________________________
#' @title sampleNpc
#'
#' @description This function samples a specified percentage of a dataframe (specifically a subset of the metadata of a Seurat object) and returns the corresponding cell IDs.
#' @param metaDF A dataframe representing a subset of the metadata of a Seurat object. Default: Subset of 'MetaData' for which 'Pass' is TRUE.
#' @param pc The percentage of the dataframe to sample, expressed as a decimal. Default: 0.1.
#' @return A vector of sampled cell IDs.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  # Example usage:
#'  # Suppose 'MetaData' is a dataframe and 'Pass' is a boolean vector with the same length.
#'  # The following example will sample 10% of the rows of 'MetaData' where 'Pass' is TRUE.
#'  sampleNpc(metaDF = MetaData[which(Pass),], pc = 0.1)
#'  }
#' }
#' @export
sampleNpc <- function(metaDF = MetaData[which(Pass),], pc = 0.1) { # Sample N % of a dataframe (obj@metadata), and return the cell IDs.
  cellIDs = rownames(metaDF)
  nr_cells = floor(length(cellIDs) * pc)
  cellIDs.keep = sample(cellIDs, size = nr_cells, replace = FALSE)
  return(cellIDs.keep)
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
#' @param set.all.genes Create the "all.genes" variable in the global env?, Default: TRUE
#' @param show Show plot? Default: TRUE
#' @examples
#' \dontrun{
#' if(interactive()){
#'  combined.obj <- calc.q99.Expression.and.set.all.genes(obj = combined.obj, quantileX = 0.9,
#'  max.cells =  25000)
#'  head(sort(as.numeric.wNames(obj@misc$expr.q90), decreasing = T))
#'  combined.obj <- calc.q99.Expression.and.set.all.genes(obj = combined.obj, quantileX = 0.95,
#'  max.cells =  25000, set.all.genes = FALSE)
#'  }
#' }
#' @seealso
#'  \code{\link[sparseMatrixStats]{character(0)}}
#' @export
#' @importFrom tictoc tic toc
#' @importFrom sparseMatrixStats rowQuantiles
calc.q99.Expression.and.set.all.genes <- function(obj = combined.obj # Calculate the gene expression of the e.g.: 90th quantile (expression in the top 10% cells).
                                                  , quantileX = 0.99, max.cells =  1e5
                                                  , slot = "data", assay = c("RNA", "integrated")[1]
                                                  , set.all.genes = TRUE, show = TRUE) {
  tictoc::tic()
  x = GetAssayData(object = obj, assay = assay, slot = slot) #, assay = 'RNA'
  if (ncol(x) > max.cells) {
    dsampled = sample(x = 1:ncol(x), size = max.cells)
    x = x[ , dsampled]
  }
  qname = paste0("q", quantileX * 100)
  slot_name = kpp("expr", qname)

  # expr.q99 = iround(apply(x, 1, quantile, probs = quantileX) )
  print("Calculating Gene Quantiles")
  expr.q99.df = sparseMatrixStats::rowQuantiles(x, probs = quantileX)
  expr.q99 = iround(expr.q99.df)

  log2.gene.expr.of.the.90th.quantile <- as.numeric(log2(expr.q99 + 1)) # strip names
  n.cells <- floor(ncol(obj) * (1-quantileX) )
  qnameP <- p0(100*quantileX,'th quantile')
  try(
    ggExpress::qhistogram(log2.gene.expr.of.the.90th.quantile, ext = "pdf", breaks = 30
                          , plotname = paste("Gene expression in the", qnameP )
                          , subtitle = kollapse(pc_TRUE(expr.q99 > 0, NumberAndPC = T), " genes have ", qname ," expr. > 0.")
                          , caption = paste(n.cells, 'cells in', qnameP)
                          , xlab = paste0("log2(expr. in the ", qnameP,"quantile+1) [UMI]")
                          , ylab = "Nr. of genes"
                          , plot = show, save = TRUE
                          , vline  = .15
                          , filtercol = T
                          , palette_use = 'npg'
                          )
    , silent = TRUE)

  {
    all.genes = percent_rank(expr.q99);
    names(all.genes) = names(expr.q99);
    all.genes <- sort(all.genes, decreasing = TRUE)
    if (set.all.genes) obj@misc$'all.genes' = all.genes = as.list(all.genes)
    assign('all.genes', all.genes, envir = as.environment(1))
  }

  obj@misc[[slot_name]] <-  expr.q99

  iprint('Quantile', quantileX ,'is now stored under obj@misc$all.genes and $', slot_name, ' Please execute all.genes <- obj@misc$all.genes.')
  return(obj)
}




# _________________________________________________________________________________________________
#' @title fix.orig.ident
#'
#' @description Remove the string "filtered_feature_bc_matrix." from "orig.ident". Helper function.
#' @param obj Seurat object, Default: merged.obj
#' @examples
#' \dontrun{
#' if(interactive()){
#'  merged.obj$orig.ident <- fix.orig.ident(obj = merged.obj); table(merged.obj$orig.ident)
#'  }
#' }
#' @export
fix.orig.ident <- function(obj = merged.obj) {
  fixed <- sub(obj$'orig.ident', pattern = 'filtered_feature_bc_matrix.', replacement = '')
  return(fixed)
}


# _________________________________________________________________________________________________
#' @title set.all.genes
#'
#' @description It is just a reminder to use calc.q99.Expression.and.set.all.genes to create the all.genes variable
#' @param obj Seurat object, Default: combined.obj
#' @examples
#' \dontrun{
#' if(interactive()){
#'  set.all.genes(); all.genes
#'  }
#' }
#' @export
set.all.genes <- function(obj = combined.obj) iprint("Use calc.q99.Expression.and.set.all.genes()")



# _________________________________________________________________________________________________
#' Create.MiscSlot
#'
#' @param obj Seurat object
#' @param NewSlotName Name of the new element inside obj@misc.
#' @export

Create.MiscSlot <- function(obj, NewSlotName = "UVI.tables", SubSlotName = NULL ) {
  if (is.null(obj@misc[[NewSlotName]])) obj@misc[[NewSlotName]] <- list() else iprint(NewSlotName, "already exists in @misc.")
  if (is.null(obj@misc[[NewSlotName]][[SubSlotName]])) obj@misc[[NewSlotName]][[SubSlotName]] <- list() else iprint(SubSlotName, "subslot already exists in @misc$NewSlot.")
  return(obj)
}



# _________________________________________________________________________________________________
#' Transfer labels from a reference Seurat object to a query Seurat object
#'
#' This function takes a query Seurat object and a path to a reference Seurat object,
#' finds transfer anchors, transfers labels, and visualizes the combined object.
#'
#' @title Transfer Labels in Seurat
#'
#' @description Function to transfer labels from a reference Seurat object to a query Seurat object using anchoring and transfer data methods from the Seurat package.
#' It then visualizes the reference and the combined objects using Uniform Manifold Approximation and Projection (UMAP).
#'
#' @param query_obj A Seurat object for which the labels are to be transferred.
#' @param reference_path A character string indicating the file path to the reference Seurat object. The path must exist.
#' @param reference_obj Alternative to `reference_path`. If provided, the path is not used to load the reference data.
#' @param named_ident A character string specifying the name of the identity class to be used from the reference Seurat object. Default is 'RNA_snn_res.0.3.ordered.ManualNames'.
#' @param new_ident A character string specifying the name of the new identity class to be created in the query Seurat object. Default is obtained by replacing 'ordered' with 'transferred' in named_ident.
#' @param predictions_col A character string specifying the column in the metadata of the transferred Seurat object containing the transferred labels. Default is 'predicted.id'.
#' @param save_anchors save anchors as RDS file.
#' @param suffix A character string to be used as a suffix in the visualization. Default is 'NEW'.
#' @param ... Additional arguments passed to the Seurat.utils::clUMAP function.
#'
#' @return The modified query Seurat object with the transferred labels as a new identity class.
#'
#' @export
#'
#' @examples
#' # combined.objX <- transfer_labels_seurat(named_ident = 'RNA_snn_res.0.3.ordered.ManualNames',
#' #                                     reference_path = '~/Dropbox (VBC)/Abel.IMBA/Metadata.D/CON.meta/label.transfer/sc6/reference.obj.sc6.DIET.2023.07.19_13.24.Rds.gz',
#' #                                     query_obj = combined.obj)
transfer_labels_seurat <- function(query_obj, reference_path
                                   , reference_obj = NULL
                                   , anchors = NULL
                                   , named_ident = 'RNA_snn_res.0.3.ordered.ManualNames'
                                   , new_ident = gsub(pattern = 'ordered'
                                                      , replacement = 'transferred'
                                                      , x = named_ident)
                                   , predictions_col = 'predicted.id'
                                   , save_anchors = TRUE
                                   , suffix = "NEW"
                                   , plot_reference = TRUE
                                   , ...) {

  print(named_ident)
  if (is.null(reference_obj)) {
    iprint("Loading reference object:", basename(reference_path))
    stopifnot(file.exists(reference_path))
    reference_obj <- readr::read_rds(reference_path)
  } else {
    stopifnot(inherits(reference_obj, "Seurat") & min(dim(reference_obj)) > 10)
  }

  # Visualize reference object
  # Seurat.utils::
  if (plot_reference) clUMAP(obj = reference_obj, ident = named_ident, suffix = 'REFERENCE', sub = 'REFERENCE', ...)


  if (is.null(anchors)) {
    print("Find anchors")
    anchors <- Seurat::FindTransferAnchors(reference = reference_obj, query = query_obj)
    if (save_anchors) isave.RDS(obj = anchors, inOutDir = TRUE)
  } else { print("Anchors provided") }


  print("Transfer labels")
  transferred_clIDs <- Seurat::TransferData(anchorset = anchors
                                            , refdata = reference_obj@meta.data[, named_ident])

  # Add metadata to combined object
  query_obj <- Seurat::AddMetaData(object = query_obj, metadata = transferred_clIDs[, predictions_col]
                                   , col.name = new_ident)

  # Visualize combined object
  # Seurat.utils::
    clUMAP(ident = new_ident, obj = query_obj, suffix = , ...)

  return(query_obj)
}


# _________________________________________________________________________________________________
#' Match and Translate Best Identity
#'
#' @title Match and Translate Best Identity
#'
#' @description This function matches the best identity from `ident_from` to `ident_to` in an object,
#' updates the metadata of the object with this new identity and returns the updated object. Additionally,
#' it generates a UMAP plot based on the new identity. The function replaces original categories with
#' the most frequent ones, hence helps to filter out less important categories.
#'
#' @param obj The object to update. This object must have a `meta.data` attribute which is a data frame
#'   containing columns named as `ident_from` and `ident_to`.
#' @param ident_from A string. The name of the column in `obj@meta.data` that is used as the source of identities.
#'   There is no default value for this parameter.
#' @param ident_to A string. The name of the column in `obj@meta.data` that is used as the target of identities.
#'   There is no default value for this parameter.
#' @param to_suffix A string. The suffix to add to the new identity name. Default is the output of the `FixPlotName`
#'   function applied to the `ident_from` string, with all alphabetical and underscore characters removed.
#' @param new_ident_name A string. The name for the newly created identity column in `obj@meta.data`.
#'   Default is a concatenation of `ident_from`, "best.match", and `to_suffix` using `kpp` function.
#' @param ... Additional parameters to be passed to `replace_by_most_frequent_categories` function.
#'
#' @return An updated version of `obj` with an additional column in `obj@meta.data` named as `new_ident_name`
#'   representing the new identity. The function also generates a UMAP plot based on this new identity.
#'
#' @seealso \code{\link[clUMAP]{clUMAP}}, \code{\link[kpp]{kpp}}, \code{\link[FixPlotName]{FixPlotName}},
#'   \code{\link[replace_by_most_frequent_categories]{replace_by_most_frequent_categories}}
#'
#' @examples
#' \dontrun{
#' updated_obj <- match_best_identity(my_obj, "origin_identity", "target_identity")
#' }
#' @export
match_best_identity <- function(obj, ident_from
                                , ident_to = gsub(pattern = 'ordered', replacement = 'transferred', x = ident_from)
                                , to_suffix = FixPlotName(gsub(pattern = '[a-zA-Z_]', replacement = "", x = ident_from))
                                , new_ident_name = kpp(ident_from, "best.match", to_suffix)
                                , ...){
  dictionary <- obj@meta.data[, c(ident_from, ident_to)]


  translation <- replace_by_most_frequent_categories(
    df = dictionary, show_plot = TRUE, suffix_barplot = ident_from, ...)

  obj@meta.data[, new_ident_name] <- translation[,1]

  px <- clUMAP(ident = new_ident_name, obj = obj)
  print(px)
  return(obj)
}


# _________________________________________________________________________________________________
#' Replace Categories by the Most Frequent Match
#'
#' This function replaces each category in a query column of a data frame with the most
#' frequently corresponding category in a reference column. It calculates the assignment quality,
#' reports it, and optionally plots it.
#'
#' @title replace_by_most_frequent_categories
#' @description Replace each category in the 'query_col' column of a data frame
#'              with the most frequently corresponding category in the 'ref_col' column.
#'
#' @param df A data frame containing the data.
#' @param query_col The name of the column in 'df' whose categories are to be replaced.
#'                  By default, the first column of 'df' is used.
#' @param ref_col The name of the column in 'df' used as reference for replacement.
#'                By default, the second column of 'df' is used.
#' @param show_plot Logical, whether to plot assignment quality. Defaults to TRUE.
#' @param suffix_barplot Suffix for barplot.
#' @param ... Additional parameters passed to the qbarplot function.
#'
#' @return A data frame with categories in 'query_col' replaced by the most frequent match
#'         from 'ref_col'.
#'
#' @examples
#' \dontrun{
#' replace_by_most_frequent_categories(df = my_data)
#' (MXX <- as.tibble(structure(c("Adjut", "Adjut", "Yearn", "Adjut", "Dwarf", "Adjut",
#' "Dwarf", "Adjut", "Dwarf", "Yearn", "Dwarf", "Dwarf", "Dwarf",
#' "Yearn", "Dwarf", "Dwarf", "Dwarf", "Zebra", "Yucca", "Plyer",
#' "Blaze", "Blaze", "Dazed", "Blaze", "Swept", "Bold", "Vixen",
#' "Bold", "Swept", "Dazed", "Mirth", "Witch", "Vixen", "Dazed",
#' "Swept", "Mirth", "Swept", "Vexed", "Query", "Yolk")
#'  , .Dim = c(20L, 2L), .Dimnames =
#'  list(NULL, c("RNA_snn_res.0.1.ordered", "RNA_snn_res.0.3.ordered" )))))
#'
#' z <- replace_by_most_frequent_categories(df = MXX)
#' head(cbind(MXX[,1],z[,1]))
#' }
#'
replace_by_most_frequent_categories <- function(df, query_col = colnames(df)[1]
                                                , ref_col = colnames(df)[2]
                                                , show_plot = TRUE
                                                , suffix_barplot = NULL
                                                , ...) {
  # Convert to data frame if it is not
  if(!is.data.frame(df)) {
    df <- as.data.frame(df)
  }

  cat_query <- unique(df[[query_col]])
  iprint(l(cat_query), "query categories:", head(cat_query))
  cat_ref <- unique(df[[ref_col]])
  iprint(l(cat_ref), "query categories:", head(cat_ref))

  # Create a table of the most frequent reference values for each query category
  replacement_table <- df %>%
    dplyr::group_by(!!sym(query_col), !!sym(ref_col)) %>%
    dplyr::summarise(n = n(), .groups = 'drop') %>%
    dplyr::arrange(!!sym(query_col), desc(n)) %>%
    dplyr::filter(!duplicated(!!sym(query_col)))

  replacement_table[[ref_col]] <- make.unique(replacement_table[[ref_col]])

  # Calculate assignment quality
  total_counts <- table(df[[query_col]])
  quality <- replacement_table$n / total_counts[replacement_table[[query_col]]]
  names(quality) <- paste0(names(quality), '->' , replacement_table[[ref_col]])

  # Report assignment quality
  message("Assignment quality (proportion of total matches):")
  # print(setNames(quality, replacement_table[[query_col]]))

  # Plot assignment quality
  if (show_plot) {
    # barplot(quality, main = "Assignment Quality", xlab = query_col, ylab = "Proportion of Total Matches")
    px <- qbarplot(quality, label = percentage_formatter(quality)
                   , suffix = suffix_barplot
                   , plotname = "Assignment Quality"
                   , filename = make.names(kpp("Assignment Quality", suffix_barplot, "pdf"))
                   , subtitle = paste("From", colnames(df)[1], "->", colnames(df)[2])
                   , xlab = paste("Best query match to reference")
                   , ylab = "Proportion of Total Matches"
                   , ...)
    print(px)

  }

  # Replace the query values with the most frequent reference values
  df[[query_col]] <- replacement_table[[ref_col]][match(df[[query_col]], replacement_table[[query_col]])]

  return(df)
}



# _________________________________________________________________________________________________
# Plot metadata ______________________________ ----
# _________________________________________________________________________________________________

#' @title plot.Metadata.Cor.Heatmap
#'
#' @description Plots a heatmap of metadata correlation values.
#' @param columns A vector of column names for which to calculate correlations. Default: c("nCount_RNA", "nFeature_RNA", "percent.mito", "percent.ribo").
#' @param cormethod The method to calculate correlations, either "pearson" or "spearman". Default: "pearson".
#' @param main The main title for the plot. Default: "Metadata correlations" followed by the correlation method.
#' @param obj The main Seurat object used for calculations. Default: combined.obj.
#' @param w The width of the plot. Default: 10.
#' @param h The height of the plot. Default: width of the plot (w).
#' @param ... Additional parameters passed to the internally called functions.
#' @seealso
#'  \code{\link[ggcorrplot]{ggcorrplot}}
#' @export plot.Metadata.Cor.Heatmap
#' @importFrom ggcorrplot ggcorrplot

plot.Metadata.Cor.Heatmap <- function(
    columns = c( "nCount_RNA", "nFeature_RNA", "percent.mito", "percent.ribo")
    , cormethod = c('pearson', 'spearman')[1]
    , main =paste( "Metadata", cormethod,"correlations")
    , obj = combined.obj
    , w = 10, h = w
    , ...){

  meta.data <- obj@meta.data
  columns.found <- intersect(colnames(meta.data), columns)

  corX <- cor(meta.data[ , columns.found], method = cormethod)
  pl <- ggcorrplot::ggcorrplot(corX, hc.order = TRUE, title = main
                               , type = "full", lab = T)
  ggExpress::qqSave(pl, fname = ppp(make.names(main),'pdf'), w = w, h = h)
  pl
}




# _________________________________________________________________________________________________
#' @title plot.Metadata.median.fraction.barplot
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
#' @export plot.Metadata.median.fraction.barplot
#' @importFrom dplyr summarize_all
#' @importFrom reshape2 melt

plot.Metadata.median.fraction.barplot <- function(
    columns = c(  "percent.mito", "percent.ribo")
    , suffix =  NULL
    , group.by = GetClusteringRuns(obj = obj)[2]
    , method = c('median', 'mean' )[1]
    , min.thr = 2.5 # At least this many percent in at least 1 cluster
    , return.matrix = F
    , main = paste( method, "read fractions per transcript class and cluster", suffix)
    , ylab = "Fraction of transcriptome (%)"
    , percentify = T
    , subt = NULL
    , position = position_stack()
    , w = 10, h = 6
    , obj = combined.obj
    , ...){

  meta.data <- obj@meta.data
  stopifnot(group.by %in% colnames(meta.data))
  columns.found <- intersect(colnames(meta.data), c(group.by, columns) )

  (mat.cluster.medians1 <- meta.data[ , columns.found] %>%
      group_by_at(group.by) %>%
      dplyr::summarize_all(median)
  )
  if (min.thr>0) {
    pass.cols <- colMax(mat.cluster.medians1[,-1]) > (min.thr/100)
    cols.OK <- which_names(pass.cols)
    cols.FAIL <- which_names(!pass.cols)
    subt = paste(length(cols.FAIL), "classed do not reach", min.thr, "% :", kpps(cols.FAIL))
    iprint(subt)
    mat.cluster.medians1 <- mat.cluster.medians1[ , c( group.by, cols.OK) ]
  }


  mat.cluster.medians <- mat.cluster.medians1 %>%
    reshape2::melt(id.vars = c(group.by), value.name = "Fraction")


  if (percentify)  mat.cluster.medians$'Fraction' = 100*mat.cluster.medians$'Fraction'

  pl <- ggbarplot(mat.cluster.medians, x = group.by, y = 'Fraction', fill = 'variable'
                  , position = position
                  , title = main, subtitle = subt ,ylab = ylab)
  ggExpress::qqSave(pl, fname = ppp(make.names(main),'pdf'), w = w, h = h)
  pl
  if (return.matrix) mat.cluster.medians1 else pl
}


# _________________________________________________________________________________________________
#' plot.Metadata.categ.pie
#'
#' @param metacol Which column in metadata should be used?
#' @param plot_name Plot name
#' @param obj Seurat object, Default: combined.obj
#' @param max.categs Maximum number of categories allowed for this pie chart. Error otherwise.
#' @param both_pc_and_value Report both percentage AND number.
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#'
#' @export

plot.Metadata.categ.pie <- function(metacol = 'Singlet.status'
                                    , plot_name = paste(metacol, "distribution")
                                    , obj = combined.obj, max.categs = 20, both_pc_and_value = T, ...) {
  categ_pivot <- table(obj[[metacol]])
  stopifnot(length(categ_pivot) < max.categs)
  qpie(categ_pivot, plotname = plot_name
       , both_pc_and_value = both_pc_and_value
       , LegendSide = F, labels = NULL, LegendTitle = '', ...)
}




