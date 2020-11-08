######################################################################
# metadata.manipulation.R
######################################################################
# source('~/GitHub/Packages/Seurat.utils/Seurat.object.manipulations.etc.R')
# Source: self + web

# Requirements ------------------------
# May also require
# try (source('/GitHub/Packages/CodeAndRoll/CodeAndRoll.R'),silent= F) # generic utilities funtions
# require('MarkdownReportsDev') # require("devtools") # plotting related utilities functions # devtools::install_github(repo = "vertesy/MarkdownReportsDev")

# ------------------------------------------------------------------------------------
GetClusteringRuns <- function(obj = combined.obj, res = F, pat = "*snn_res.*[0,1]\\.[0-9]$") { # Get Clustering Runs: metadata column names
  if (res) pat = gsub(x = pat, pattern = '\\[.*\\]', replacement = res)
  clustering.results <- grepv(x = colnames(obj@meta.data), pattern = pat)
  if ( identical(clustering.results, character(0)) ) warning("No matching column found!")
  return(clustering.results)
}
# GetClusteringRuns()


# ------------------------------------------------------------------------------------
GetNamedClusteringRuns <- function(obj = combined.obj  # Get Clustering Runs: metadata column names
                                   , res = c(F, 0.5)[1], topgene = F, pat = "^cl.names.Known.*[0,1]\\.[0-9]$") {
  if (res) pat = gsub(x = pat, pattern = '\\[.*\\]', replacement = res)
  if (topgene) pat = gsub(x = pat, pattern = 'Known', replacement = 'top')
  clustering.results <- grepv(x = colnames(obj@meta.data), pattern = pat)
  if ( identical(clustering.results, character(0)) ) warning("No matching column found!")
  return(clustering.results)
}
# GetNamedClusteringRuns()


# ------------------------------------------------------------------------------------
GetOrderedClusteringRuns <- function(obj = combined.obj, res = F, pat = "*snn_res.*[0,1]\\.[0-9].*ordered$") { # Get Clustering Runs: metadata column names
  if (res) pat = gsub(x = pat, pattern = '\\[.*\\]', replacement = res)
  clustering.results <- grepv(x = colnames(obj@meta.data), pattern = pat)
  if ( identical(clustering.results, character(0)) ) warning("No matching column found!")
  return(clustering.results)
}
# GetOrderedClusteringRuns(); GetOrderedClusteringRuns(res = 0.7)


# ------------------------------------------------------------------------------------
GetNumberOfClusters <- function(obj = combined.obj) { # Get Number Of Clusters
  clustering.results <- GetClusteringRuns(obj)
  print("## Number of clusters: ---------")
  for (cc in clustering.results) {
    NrCl <- length(unique(obj@meta.data[[cc]]))
    iprint( cc, "   ", NrCl)
  }
}
# GetNumberOfClusters()

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
  iprint(length(idx.matching.cells), 'cells found.')
  return(rownames(obj@meta.data)[idx.matching.cells])
}
# getCellIDs.from.meta()

# seu.add.meta.from.vector ------------------------------------------------------------------------
seu.add.meta.from.vector <- function(obj = combined.obj, metaD.colname = metaD.colname.labeled, Label.per.cell=Cl.Label.per.cell ) { # Add a new metadata column to a Seurat  object
  obj@meta.data[, metaD.colname ] = Label.per.cell
  iprint(metaD.colname, "contains the named identitites. Use Idents(combined.obj) = '...'. The names are:", unique(Label.per.cell))
  return(obj)
}
# combined.obj <- add.Cl.Label.2.Metadata(obj = combined.obj, metaD.colname = metaD.colname.labeled, Label.per.cell=Cl.Label.per.cell )
# formerly add.Cl.Label.2.Metadata


# seu.map.and.add.new.ident.to.meta ------------------------------------------------------------------------

seu.map.and.add.new.ident.to.meta <- function(obj = combined.obj, ident.table = clusterIDs.GO.process
                                              , metaD.colname = substitute(ident.table) ) { # Add a new metadata column to a Seurat  object
  ident.vec <- as.named.vector(ident.table)

  # identities should match ----------------
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

  stopif(l(OnlyInIdentVec), message = msg.IdentVec)
  stopif(l(OnlyInSeuratIdents), message = msg.Seu)

  # identity mapping ----------------
  new.ident <- translate(vec = as.character(Idents(obj)), oldvalues = ident.X, newvalues = ident.Y)
  obj@meta.data[[metaD.colname]] = new.ident
  iprint(metaD.colname, "contains the named identitites. Use Idents(combined.obj) = '...'. The names are:"); cat(paste0("\t", ident.Y, "\n"))
}
# combined.obj <- seu.map.and.add.new.ident.to.meta(obj = combined.obj, ident.table = clusterIDs.GO.process)


# calc.cluster.averages ------------------------------------------------
calc.cluster.averages <- function(obj =  combined.obj, simplify=T, plotit = T
                                  , col_name = "Score.GO.0006096"
                                  , split_by = GetClusteringRuns()[l(GetClusteringRuns())]
                                  , quantile.thr = 0.9
                                  , ylab.text = "Glycolytic Process"
                                  , title = paste("Cluster Average", col_name, ylab.text)
                                  , subtitle = paste("Clusters above the",percentage_formatter(0.9),"quantile")
                                  , ylb = paste(ylab.text, col_name)
                                  , xlb = paste ( "Clusters", split_by)
                                  , fname = ppp(col_name,split_by,"cluster.average.barplot.pdf")
) { # calc.cluster.averages of a m

  df.summary <-
    obj@meta.data %>%
    select_at(c(col_name, split_by)) %>%
    group_by_at(split_by) %>%
    summarize('nr.cells' = n()
              , 'median' =  median(!!sym(col_name), na.rm = TRUE)
    )



  if (simplify) {
    av.score <- as.named.vector(df_col = df.summary[,"median"])
    names(av.score) <- (as.numeric(names(av.score))-1)

    if (plotit) {
      wbarplot(av.score, hline = quantile(av.score, quantile.thr)
               , plotname = fname
               , xlab = xlb
               , ylab = ylb
               , main = title
               , sub = subtitle)
    }
    av.score

  } else {
    df.summary
  }
}

# calc.cluster.averages(col_name = "Score.GO.0006096", split_by = grepv(pattern = "0.6", GetClusteringRuns())                        )
# av.score <- as.named.vector(df_col = calc.cluster.averages()[,"median"])
# names(av.score) <- as.numeric(names(av.score))
# wbarplot(av.score, hline = quantile(av.score, 0.9)
#          , ylab = "Glycolytic Process Score (GO.0006096)"
#          , main = "Two clusters fall above the 90% quantile")



# Add to obj@metadata from an external table ------------------------------------------------------------------------
seu.add.meta.from.table <- function(obj = seu.ORC, meta = MetaData.ORC, suffix = ".fromMeta") { # Add multiple new metadata columns to a Seurat object from a table.
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
} # x=seu.add.meta.from.table()


# sampleNpc ------------------------------------------------------------------------
sampleNpc <- function(metaDF = MetaData[which(Pass),], pc=0.1) { # Sample N % of a dataframe (obj@metadata), and return the cell IDs.
  cellIDs = rownames(metaDF)
  nr_cells = floor(length(cellIDs) * pc)
  cellIDs.keep = sample(cellIDs, size = nr_cells, replace = FALSE)
  return(cellIDs.keep)
}


# Calcq90Expression ------------------------------------------------------------------------
Calcq90Expression <- function(obj = combined.obj # Calculate the gene expression of the e.g.: 90th quantile (expression in the top 10% cells).
                              , quantileX=0.9, assay = c("RNA", "integrated")[1]
                              , slot = "data", max.cells =  100000) {
  tic()
  x = GetAssayData(object = obj, assay = assay, slot = slot) #, assay = 'RNA'
  if (ncol(x) > max.cells) {
    dsampled = sample(x = 1:ncol(x), size = max.cells)
    x = x[ , dsampled]
  }
  expr.q90 = iround(apply(x, 1, quantile, probs = quantileX) )
  toc();

  log2.gene.expr.of.the.90th.quantile <- log2(expr.q90 + 1)
  suppressWarnings(
    whist(log2.gene.expr.of.the.90th.quantile, breaks = 30
          , xlab = "log2(expr.q90+1) [UMI]", ylab = "Cells", vline  = .2, filtercol = T)
  )

  all.genes = percent_rank(expr.q90); names(all.genes) = names(expr.q90); all.genes <- sort.decreasing(all.genes)

  obj@misc$'all.genes' = all.genes = as.list(all.genes)
  obj@misc$'expr.q90' = expr.q90
  assign('all.genes', all.genes, envir = as.environment(1))

  iprint('Quantile', quantileX ,'is now stored under obj@misc$all.genes and $expr.q90. Please execute all.genes <- obj@misc$all.genes.')
  return(obj)
}
# combined.obj <- Calcq90Expression(obj = combined.obj)
# head(sort(as.numeric.wNames(obj@misc$expr.q90), decreasing = T))

# PlotTopGenes ------------------------------------------------------------------------
PlotTopGenes <- function(obj = combined.obj, n=32 ){ # Plot the highest expressed genes on umaps, in a subfolder. Requires calling Calcq90Expression before.
  Highest.Expressed.Genes = names(head(sort(obj@misc$expr.q90, decreasing = T), n = n))
  multiFeaturePlot.A4(list.of.genes = Highest.Expressed.Genes, foldername = "Highest.Expressed.Genes" )
}
# PlotTopGenes()


# fix.orig.ident ------------------------------------------------------------------------
fix.orig.ident <- function(obj = merged.obj) {
  fixed <- sub(obj$'orig.ident', pattern = 'filtered_feature_bc_matrix.', replacement = '')
  return(fixed)
}
# merged.obj$orig.ident <- fix.orig.ident(obj = merged.obj); table(merged.obj$orig.ident)

# set.all.genes ------------------------------------------------------------------------
set.all.genes <- function(obj = combined.obj) iprint("Use Calcq90Expression()")
# set.all.genes(); all.genes

# recall.all.genes ------------------------------------------------------------------------
recall.all.genes <- function(obj = combined.obj) {
  if(!exists('all.genes')) {
    all.genes <- obj@misc$all.genes
    print(head(unlist(all.genes)))
    ww.assign_to_global(name = "all.genes", value = all.genes)
  } else {print("variable 'all.genes' exits in the global namespace")}
}
# recall.all.genes(); all.genes


# recall.parameters ------------------------------------------------------------------------
recall.parameters <- function(obj = combined.obj) {
  if(exists('p')) iprint("variable 'p' exits in the global namespace:"); print(p); print("Now it will be overwritten.")
  p <- obj@misc$p
  print(head(unlist(p)))
  ww.assign_to_global(name = "p", value = p)
  # if(!exists('p')) {
  #   all.genes <- obj@misc$all.genes
  #   ww.assign_to_global(name = "all.genes", value = all.genes)
  # } else {
  #   print("variable 'p' exits in the global namespace")
  # }
}
# recall.parameters(); p


# plot.expression.rank.q90 ------------------------------------------------------------------------
plot.expression.rank.q90 <- function(obj = combined.obj, gene="ACTB", filterZero=T) {
  expr.GOI <- obj@misc$expr.q90[gene]
  expr.all <- unlist(obj@misc$expr.q90)
  gene.found <- gene %in% names(expr.all)
  stopifnot(gene.found)

  if (expr.GOI==0) iprint(gene, "is not expressed. q90-av.exp:", expr.GOI) else
    if (expr.GOI<0.05) iprint(gene, "is lowly expressed. q90-av.exp:", expr.GOI)
    if (filterZero) {
      iprint("Zero 'q90 expression' genes (", pc_TRUE(expr.all == 0), ") are removed.")
      expr.all <- expr.all[expr.all > 0]
    }
  counts <- sum(obj@assays$RNA@counts[gene,])
  if (expr.GOI == 0) {
    quantile.GOI <- 0
    title <- paste(gene, "is too lowly expressed: q90-av.exp is zero. \n There are", counts,"counts." )
  } else {
    pos.GOI <- which(names(expr.all) == gene)
    quantile.GOI <- ecdf(expr.all)(expr.all)[pos.GOI]
    title <- paste(gene, "is in the", percentage_formatter(quantile.GOI), "quantile of 'q90-av' expression. \n There are", counts,"counts" )
  }
  suppressWarnings(
    whist(expr.all, vline = expr.GOI, breaks = 100, main = title, plotname =   make.names(title)
          , ylab = "Genes"
          , xlab = "Av. mRNA in the 10% top expressing cells (q90 av.exp.)")
  )
}
# plot.expression.rank.q90(gene = "SATB2")




# FlipReductionCoordinates ------------------------------------------------------------------------
FlipReductionCoordinates <- function(obj = combined.obj, dim=2, reduction="umap", flip=c('x', 'y', 'xy', NULL)[1]) { # Set active UMAP to `obj@reductions$umap` from `obj@misc$reductions.backup`.
  coordinates <- Embeddings(obj, reduction = reduction)
  stopifnot(ncol(coordinates) == dim )

  if (flip %in% c('x', 'xy')) coordinates[,1] = coordinates[,1] * -1
  if (flip %in% c('y', 'xy')) coordinates[,2] = coordinates[,2] * -1
  obj@reductions[[reduction]]@cell.embeddings <- coordinates
  return(obj)
}


# SeuratColorVector ------------------------------------------------------------------------
SeuratColorVector <- function(obj = combined.obj) {
  colorlevels <- hue_pal()(length(levels(obj@active.ident)))
  translate(vec = as.character(obj@active.ident)
            , oldvalues = levels(obj@active.ident)
            , newvalues = colorlevels)
}


#  ------------------------------------------------------------------------
#  ------------------------------------------------------------------------
#  ------------------------------------------------------------------------
#  ------------------------------------------------------------------------
