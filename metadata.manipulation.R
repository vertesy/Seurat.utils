######################################################################
# metadata.manipulation.R
######################################################################
# source ('~/GitHub/SeuratUtil/Seurat.object.manipulations.etc.R')
# Source: self + web

# Requirements ------------------------
# May also require
# try (source ('~/GitHub/TheCorvinas/R/CodeAndRoll.R'),silent= F) # generic utilities funtions
# require('MarkdownReportsDev') # require("devtools") # plotting related utilities functions # devtools::install_github(repo = "vertesy/MarkdownReportsDev")


# get Cells from metadata  ------------------------------------------------
getMetadataColumn <- mmeta <- function(ColName.metadata = 'batch', obj = combined.obj, as_numeric =F) { # Get a metadata column from a Seurat object as a named vector
  stopifnot(ColName.metadata %in% colnames(obj@meta.data))

  x = as.named.vector(obj@meta.data[ ,ColName.metadata, drop=F])
  if (as_numeric) {
    as.numeric.wNames(x)+1
  } else {x}
}

# GetCellIDs from metadata ---------------
getCellIDs.from.meta <- function(obj=org, ColName.meta = 'res.0.6', values = NA) { # Get cellIDs from a metadata column, matching a list of values (using %in%).
  idx.matching.cells = which(obj@meta.data[ , ColName.meta] %in% values)
  iprint(l(idx.matching.cells), 'cells found.')
  return(rownames(obj@meta.data)[idx.matching.cells])
}
# getCellIDs.from.meta()

# seu.add.meta.from.vector------------------------------------------------------------------------
seu.add.meta.from.vector <- function(obj = combined.obj, metaD.colname = metaD.colname.labeled, Label.per.cell=Cl.Label.per.cell ) { # Add a new metadata column to a Seurat  object
  obj@meta.data[, metaD.colname ] = Label.per.cell
  iprint(metaD.colname, "contains the named identitites. Use Idents(combined.obj) = '...'. The names are:", unique(Label.per.cell))
  return(obj)
}
# combined.obj <- add.Cl.Label.2.Metadata(obj = combined.obj, metaD.colname = metaD.colname.labeled, Label.per.cell=Cl.Label.per.cell )
# formerly add.Cl.Label.2.Metadata


# Add to obj@metadata from an external table ------------------------------------------------------------------------
seu.add.meta.from.table <- function(obj = seu.ORC, meta = MetaData.ORC, suffix = ".fromMeta") { # Add multiple new metadata columns to a Seurat object from a table.
  NotFound  = setdiff(colnames(obj), rownames(meta))
  Found     = intersect(colnames(obj), rownames(meta))
  if (l(NotFound)) iprint(l(NotFound), 'cells were not found in meta, e.g.: ', trail(NotFound, N=10))

  mCols.new = colnames(meta)
  mCols.old = colnames(obj@meta.data)
  overlap = intersect(mCols.new, mCols.old)
  if (l(overlap)) {
    iprint(l(overlap), 'metadata columns already exist in the seurat object: ', overlap, '. These are tagged as: *', suffix)
    colnames(meta)[overlap] = p0(overlap, suffix)
  }
  mCols.add = colnames(meta)
  obj@meta.data[Found, mCols.add] = meta[ Found,]

  return(obj)
} # x=seu.add.meta.from.table()
