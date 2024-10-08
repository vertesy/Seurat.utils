######################################################################
# metadata.manipulation.R
######################################################################
# source('~/GitHub/Packages/Seurat.utils/Functions/Seurat.object.manipulations.etc.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Metadata.manipulation.R"))
# Source: self + web

# - getMedianMetric.lsObj
# - add.meta.tags
# - addMetaFraction
# - GetClusteringRuns
# - GetNamedClusteringRuns
# - GetOrderedClusteringRuns
# - GetNumberOfClusters
# - getMetadataColumn <- mmeta
# - getCellIDs.from.meta
# - seu.add.meta.from.vector
# - seu.map.and.add.new.ident.to.meta
# - calc.cluster.averages
# - seu.add.meta.from.table
# - sampleNpc
# - calc.q90.Expression.and.set.all.genes
# - PlotTopGenes
# - fix.orig.ident
# - set.all.genes
# - recallAllGenes
# - recall.meta.tags.n.datasets
# - recall.parameters
# - recall.genes.l
# - save.parameters
# - plot.expression.rank.q90
# - FlipReductionCoordinates
# - SeuratColorVector
# - getClusterColors
# - get.clustercomposition

# getMedianMetric.lsObj ------------------------------------------------------------------------------------------------
getMedianMetric.lsObj <- function(ls.Obj = ls.Seurat, n.datasets = length(ls.Seurat), mColname = "percent.mito") {
  medMetric <- vec.fromNames(names(ls.Seurat))
  for(i in 1:n.datasets ) {
    medMetric[i] <- median(ls.Seurat[[i]]@meta.data[,mColname])
  }
  return(medMetric)
}
# ls.Seurat <- getMedianMetric.lsObj(ls.Obj = ls.Seurat, n.datasets = length(ls.Seurat), mColname = "percent.mito")


# add.meta.tags ------------------------------------------------------------------------------------------------
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
# ls.Seurat[[1]] <- add.meta.tags(list.of.tags = tags, obj = ls.Seurat[[1]], n = 1)


# addMetaFraction ------------------------------------------------------------------------------------------------
addMetaFraction <- function(col.name = "percent.mito", gene.symbol.pattern = c("^MT\\.|^MT-", F)[1]
                              , gene.set = F, obj = ls.Seurat[[1]], verbose = T) {
  stopif(condition = isFALSE(gene.set) && isFALSE(gene.symbol.pattern), "Either gene.set OR gene.symbol.pattern has to be defined (!= FALSE).")
  if (!isFALSE(gene.set) && !isFALSE(gene.symbol.pattern) && verbose) print("Both gene.set AND gene.symbol.pattern are defined. Only using gene.set.")

  geneset <- check.genes(list.of.genes = gene.set, obj = obj)
  total_expr <- Matrix::colSums(GetAssayData(object = obj))
  genes.matching <- if (!isFALSE(gene.set)) intersect(gene.set, rownames(obj)) else grepv(pattern = gene.symbol.pattern, x = rownames(obj))

  genes.expr = GetAssayData(object = obj)[genes.matching, ]
  target_expr <- if (l(genes.matching) >1) Matrix::colSums(genes.expr) else genes.expr
  obj <- AddMetaData(object = obj, metadata = target_expr / total_expr, col.name = col.name)
  colnames(obj@meta.data)
  return(obj)
}


# ls.Seurat[[1]] <- addMetaFraction(col.name = "percent.mito", gene.symbol.pattern = "^MT\\.|^MT-")
# ls.Seurat[[1]] <- addMetaFraction(col.name = "percent.ribo", gene.symbol.pattern = "^RPL|^RPS")
# ls.Seurat[[1]] <- addMetaFraction(col.name = "percent.AC.GenBank", gene.symbol.pattern = "^AC[0-9]{6}\\.")
# ls.Seurat[[1]] <- addMetaFraction(col.name = "percent.AL.EMBL", gene.symbol.pattern = "^AL[0-9]{6}\\.")
# ls.Seurat[[1]] <- addMetaFraction(col.name = "percent.LINC", gene.symbol.pattern = "^LINC0")
# ls.Seurat[[1]] <- addMetaFraction(col.name = "percent.MALAT1", gene.symbol.pattern = "^MALAT1")
# colnames(ls.Seurat[[1]]@meta.data)
# HGA_MarkerGenes <- c("ENO1", "IGFBP2", "WSB1", "DDIT4", "PGK1", "BNIP3", "FAM162A", "TPI1", "VEGFA", "PDK1", "PGAM1", "IER2", "FOS", "BTG1", "EPB41L4A-AS1","NPAS4", "HK2", "BNIP3L", "JUN", "ENO2", "GAPDH", "ANKRD37", "ALDOA", "GADD45G", "TXNIP")
# sobj <- addMetaFraction(col.name = "percent.HGA", gene.set = HGA_MarkerGenes, obj =  sobj)


# ------------------------------------------------------------------------------------
GetClusteringRuns <- function(obj = combined.obj, res = F, pat = "*snn_res.*[0-9]$") { # Get Clustering Runs: metadata column names
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
  if ( identical(clustering.results, character(0)) ) {
    print("Warning: NO matching column found! Trying GetClusteringRuns(..., pat = '*_res.*[0,1]\\.[0-9]$)")
    clustering.results <- GetClusteringRuns(obj = obj, res = F, pat = "*_res.*[0,1]\\.[0-9]$")
  }
  return(clustering.results)
}
# GetNamedClusteringRuns()


# ------------------------------------------------------------------------------------
GetOrderedClusteringRuns <- function(obj = combined.obj, res = F, pat = "*snn_res.*[0,1]\\.[0-9]\\.ordered$") { # Get Clustering Runs: metadata column names
  if (res) pat = gsub(x = pat, pattern = '\\[.*\\]', replacement = res)
  clustering.results <- grepv(x = colnames(obj@meta.data), pattern = pat)
  if ( identical(clustering.results, character(0)) ) warning("No matching column found!")
  return(clustering.results)
}
# GetOrderedClusteringRuns(); GetOrderedClusteringRuns(res = 0.5)


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
getCellIDs.from.meta <- function(ident = 'res.0.6', values = NA, obj=combined.obj, inverse = F ) { # Get cellIDs from a metadata column, matching a list of values (using %in%).
  mdat <- obj@meta.data[ , ident]
  cells <- if (inverse) {mdat %!in% values} else {mdat %in% values}
  idx.matching.cells = which(cells)
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
  ident.Seu <- gtools::mixedsort(levels(Idents(obj)))
  iprint("ident.Seu: ", ident.Seu)

  OnlyInIdentVec      <- setdiff(ident.X, ident.Seu)
  OnlyInSeuratIdents  <- setdiff(ident.Seu, ident.X)

  msg.IdentVec <- kollapse("Rownames of 'ident.table' have entries not found in 'Idents(obj)':"
                           , OnlyInIdentVec, " not found in ", ident.Seu, collapseby = " ")

  msg.Seu <- kollapse("Rownames of 'Idents(obj)' have entries not found in 'ident.table':"
                      , OnlyInSeuratIdents, " not found in ", ident.X, collapseby = " ")

  stopif (l(OnlyInIdentVec), message = msg.IdentVec)
  stopif (l(OnlyInSeuratIdents), message = msg.Seu)

  # identity mapping ----------------
  new.ident <- translate(vec = as.character(Idents(obj)), old = ident.X, new = ident.Y)
  obj@meta.data[[metaD.colname]] = new.ident
  iprint(metaD.colname, "contains the named identitites. Use Idents(combined.obj) = '...'. The names are:"); cat(paste0("\t", ident.Y, "\n"))
}
# combined.obj <- seu.map.and.add.new.ident.to.meta(obj = combined.obj, ident.table = clusterIDs.GO.process)




# calc.cluster.averages ------------------------------------------------
calc.cluster.averages <- function(col_name = "Score.GO.0006096"
                                  , plot.UMAP.too = TRUE
                                  , return.plot = F
                                  , obj =  combined.obj
                                  , split_by = GetNamedClusteringRuns()[1]
                                  , scale.zscore = FALSE
                                  , simplify=T, plotit = T
                                  , histogram = FALSE, nbins = 50
                                  , suffix = NULL
                                  , stat = c("mean", "median")[2]
                                  , quantile.thr = 0.9
                                  , absolute.thr = FALSE
                                  , filter = c(FALSE, 'above', 'below')[1]
                                  , ylab.text = paste("Cluster", stat, "score")
                                  , title = paste("Cluster", stat, col_name)
                                  , subtitle = NULL
                                  , width = 8, height =6
                                  , ...
                                  # , ylb = paste(ylab.text, col_name)
                                  # , xlb = paste("Clusters >",percentage_formatter(quantile.thr),"quantile are highlighted. |", split_by)
                                  , xlb = if (absolute.thr) paste("Threshold at", absolute.thr) else paste(
                                    "Black lines: " , kppd(percentage_formatter(c(1-quantile.thr, quantile.thr))) ,"quantiles |"
                                    , "Cl. >",percentage_formatter(quantile.thr),"are highlighted. |", split_by
                                  )

                                  , fname = ppp(col_name,split_by,"cluster.average.barplot.pdf", ...)
) { # calc.cluster.averages of a m
  iprint(substitute(obj), "split by", split_by)
  if(absolute.thr) iprint('In case of the absolute threshold, only the returned values are correct, the plot annotations are not!')

  if (plot.UMAP.too) qUMAP(obj = obj, feature = col_name)

  df.summary <-
    obj@meta.data %>%
    select_at(c(col_name, split_by)) %>%
    group_by_at(split_by) %>%
    summarize('nr.cells' = n()
              , 'mean' = mean(!!sym(col_name), na.rm = TRUE)
              , 'SEM' = sem(!!sym(col_name), na.rm = TRUE)
              , 'median' = median(!!sym(col_name), na.rm = TRUE)
              , 'SE.median' = 1.2533 * sem(!!sym(col_name), na.rm = TRUE)
    )

  if (simplify) {
    av.score <- df.summary[[stat]]
    names(av.score) <- ppp("cl",df.summary[[1]])
    av.score <- sortbyitsnames(av.score)
    if (scale.zscore) av.score <- (scale(av.score)[,1])

    cutoff <- if(absolute.thr) absolute.thr else quantile(av.score, quantile.thr)
    cutoff.low <- if(absolute.thr) NULL else  quantile(av.score, (1-quantile.thr) )

    if (plotit) {
      if (histogram) {
        p <- qhistogram(vec = as.numeric(av.score), save = F
                        , vline = cutoff
                        , plotname = ppp(title, quantile.thr)
                        , bins = nbins
                        , subtitle = paste(subtitle, "| median in blue/dashed")
                        , ylab = ylab.text
                        , xlab = xlb # Abused
                        , xlab.angle = 45
                        # , ylim=c(-1,1)
                        , ...
                        # , ext = "png", w = 7, h = 5
        ) + geom_vline(xintercept = cutoff.low, lty=2)
        print(p)
        title_ <- ppp(title, suffix, flag.nameiftrue(scale.zscore))
        qqSave(ggobj = p, title = title_, ext = "png", w = width, h = height)
      } else {
        p <- qbarplot(vec = av.score, save = F
                      , hline = cutoff
                      , title = title
                      , suffix = quantile.thr
                      , subtitle = subtitle
                      , ylab = ylab.text
                      , xlab = xlb # Abused
                      , xlab.angle = 45
                      # , ylim=c(-1,1)
                      , ...
                      # , ext = "png", w = 7, h = 5
        ) + geom_hline(yintercept = cutoff.low , lty = 2)

        print(p)
        title_ <- ppp(title, suffix, flag.nameiftrue(scale.zscore))
        qqSave(ggobj = p, title = title_, fname = ppp(title_, split_by, "png"),  w = width, h = height)
      }
    }
    print(quantile.thr)

    if (filter == 'below') {
      return(filter_LP(av.score, threshold = cutoff, plot.hist = F))
    } else if (filter == 'above') {
      return(filter_HP(av.score, threshold = cutoff, plot.hist = F))
    } else {
      return(av.score)
    }
  } else if (return.plot) { # if /not/ simplify
    return(p)
  } else {
    return(df.summary)
  }
}

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


# calc.q90.Expression.and.set.all.genes ------------------------------------------------------------------------
calc.q90.Expression.and.set.all.genes <- function(obj = combined.obj # Calculate the gene expression of the e.g.: 90th quantile (expression in the top 10% cells).
                              , quantileX=0.9, max.cells =  100000
                              , slot = "data", assay = c("RNA", "integrated")[1]
                              , set.all.genes = TRUE, show = TRUE) {
  tic()
  x = GetAssayData(object = obj, assay = assay, slot = slot) #, assay = 'RNA'
  if (ncol(x) > max.cells) {
    dsampled = sample(x = 1:ncol(x), size = max.cells)
    x = x[ , dsampled]
  }
  qname = p0("q", quantileX * 100)
  slot_name = kpp("expr", qname)

  # expr.q90 = iround(apply(x, 1, quantile, probs = quantileX) )
  expr.q90.df = sparseMatrixStats::rowQuantiles(x, probs = quantileX)
  expr.q90 = iround(as.named.vector(expr.q90.df))
  toc();

  log2.gene.expr.of.the.90th.quantile <- as.numeric(log2(expr.q90 + 1)) # strip names
  try(
    qhistogram(log2.gene.expr.of.the.90th.quantile, ext = "pdf", breaks = 30
               , plotname = kpp("log2.gene.expr.of.the ", qname," quantile")
               , subtitle = kollapse(pc_TRUE(expr.q90 > 0, NumberAndPC = T), " genes have ", qname ," expr. > 0.")
               , xlab = p0("log2(expr.",qname,"+1) [UMI]"), ylab = "Genes"
               , plot = show, save = TRUE, vline  = .2)
    , silent = TRUE)

  {
    all.genes = percent_rank(expr.q90); names(all.genes) = names(expr.q90); all.genes <- sort.decreasing(all.genes)
    if (set.all.genes) obj@misc$'all.genes' = all.genes = as.list(all.genes)
    assign('all.genes', all.genes, envir = as.environment(1))
  }

  obj@misc[[slot_name]] <-  expr.q90

  iprint('Quantile', quantileX ,'is now stored under obj@misc$all.genes and $', slot_name, ' Please execute all.genes <- obj@misc$all.genes.')
  return(obj)
}

# combined.obj <- calc.q90.Expression.and.set.all.genes(obj = combined.obj, quantileX=0.9, max.cells =  25000)
# head(sort(as.numeric.wNames(obj@misc$expr.q90), decreasing = T))
# combined.obj <- calc.q90.Expression.and.set.all.genes(obj = combined.obj, quantileX=0.95, max.cells =  25000, set.all.genes = FALSE)

# PlotTopGenes ------------------------------------------------------------------------
PlotTopGenes <- function(obj = combined.obj, n=32 ){ # Plot the highest expressed genes on umaps, in a subfolder. Requires calling calc.q90.Expression.and.set.all.genes before.
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
set.all.genes <- function(obj = combined.obj) iprint("Use calc.q90.Expression.and.set.all.genes()")
# set.all.genes(); all.genes

# recallAllGenes ------------------------------------------------------------------------
recallAllGenes <- function(obj = combined.obj) { # all.genes set by calc.q90.Expression.and.set.all.genes()
  if (!exists('all.genes')) {
    all.genes <- obj@misc$all.genes
    print(head(unlist(all.genes)))
    ww.assign_to_global(name = "all.genes", value = all.genes)
  } else {print("variable 'all.genes' exits in the global namespace")}
}
# recallAllGenes(); all.genes

# recall.meta.tags.n.datasets ------------------------------------------------------------------------
recall.meta.tags.n.datasets <- function(obj = combined.obj) {
  if (!exists('n.datasets')) {
    n.datasets <- obj@misc$n.datasets
    print(head(unlist(n.datasets)))
    ww.assign_to_global(name = "n.datasets", value = n.datasets)
  } else {print("variable 'n.datasets' exits in the global namespace")}

  if (!exists('meta.tags')) {
    meta.tags <- obj@misc$meta.tags
    print(head(unlist(meta.tags)))
    ww.assign_to_global(name = "meta.tags", value = meta.tags)
  } else {print("variable 'meta.tags' exits in the global namespace")}

}
# recall.n.datasets(); n.datasets

# recall.parameters ------------------------------------------------------------------------
recall.parameters <- function(obj = combined.obj, overwrite = FALSE) {
  if (is.null(obj@misc$'p')) {
    print("obj does not have: obj@misc$p")
  } else {
    p <- obj@misc$'p'
    print(head(p))
    if (exists('p')) iprint("variable 'p' exits in the global namespace:");

    if (!exists('p') | (exists('p') & overwrite == TRUE) ) {
      ww.assign_to_global(name = "p", value = p); print("Overwritten.")
    } else {
      print("Not overwritten.")
    }
  } # else if obj@misc$'p'
}
# recall.parameters(); p


# recall.genes.ls ------------------------------------------------------------------------
recall.genes.ls<- function(obj = combined.obj) { # genes.ls
  if (!exists('genes.ls')) {
    genes.ls <- obj@misc$genes.ls
    print(head(unlist(genes.ls)))
    ww.assign_to_global(name = "genes.ls", value = genes.ls)
  } else {print("variable 'genes.ls' exits in the global namespace")}
}
# recall.genes.ls(); genes.ls


# save.parameters ------------------------------------------------------------------------
save.parameters <- function(obj = combined.obj, params = p) {
  if (!is.null(obj@misc$'p')) print("Overwriting already existing obj@misc$p. Old version:") ; print(head(unlist(obj@misc$'p')))
  obj@misc$p <- params
}
# save.parameters(obj = combined.obj, params = p);


# plot.expression.rank.q90 ------------------------------------------------------------------------
plot.expression.rank.q90 <- function(obj = combined.obj, gene="ACTB", filterZero=T) {
  expr.GOI <- obj@misc$expr.q90[gene]
  expr.all <- unlist(obj@misc$expr.q90)
  gene.found <- gene %in% names(expr.all)
  stopifnot(gene.found)

  if (expr.GOI == 0) iprint(gene, "is not expressed. q90-av.exp:", expr.GOI) else
    if (expr.GOI < 0.05) iprint(gene, "is lowly expressed. q90-av.exp:", expr.GOI)
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
FlipReductionCoordinates <- function(obj = combined.obj, dim=2, reduction="umap"
                                     , flip=c('x', 'y', 'xy', NULL)[1], FlipReductionBackupToo = TRUE) { # Set active UMAP to `obj@reductions$umap` from `obj@misc$reductions.backup`.
  coordinates <- Embeddings(obj, reduction = reduction)
  stopifnot(ncol(coordinates) == dim )

  if (flip %in% c('x', 'xy')) coordinates[,1] = coordinates[,1] * -1
  if (flip %in% c('y', 'xy')) coordinates[,2] = coordinates[,2] * -1
  obj@reductions[[reduction]]@cell.embeddings <- coordinates

  if (FlipReductionBackupToo) {
    bac.slot <- p0(reduction,dim,"d")
    if (length(obj@misc$reductions.backup[[bac.slot]])) {
      obj@misc$reductions.backup[[bac.slot]]@cell.embeddings <- coordinates
      iprint(dim, "dimensional",reduction,"backup flipped too.")
    }
  }
  return(obj)
}
# clUMAP(); combined.obj <- FlipReductionCoordinates(combined.obj); clUMAP()



# SeuratColorVector ------------------------------------------------------------------------
SeuratColorVector <- function(ident = NULL, obj = combined.obj, plot.colors = F) {
  if (!is.null(ident)) {
    print(ident)
    ident.vec <- obj[[ident]][,1]
  } else {
    ident.vec <- obj@active.ident
  }
  ident.vec <- as.factor(ident.vec)
  print(table(ident.vec))
  colorlevels <- scales::hue_pal()(length(levels(ident.vec)))
  if (plot.colors) color_check(colorlevels)
  translate(vec = as.character(ident.vec)
            , old = levels(ident.vec)
            , new = colorlevels)
}
# SeuratColorVector()
# SeuratColorVector(ident = GetNamedClusteringRuns()[2], plot.colors = T)

# getClusterColors ------------------------------------------------------------------------
getClusterColors <- function(obj = combined.obj
                             , ident = GetClusteringRuns()[1]
                             , show = T) {
  (identities <- levels(as.factor(obj[[ident]][,1])))
  color_palette <- scales::hue_pal()(length(identities))
  # color_check(color_palette)
  # names(color_palette) <- sort(as.factor(identities))
  names(color_palette) <- ((identities))
  identvec <- obj[[ident]][,1]
  colz <- color_palette[identvec]
  names(colz) <- identvec
  if (show) color_check(unique(colz)) # MarkdownReports
  colz
}
# getClusterColors(obj = combined.obj, ident = GetClusteringRuns()[2] )


#  get.clustercomposition ------------------------------------------------------------------------
get.clustercomposition <- function(obj = combined.obj, x = 'integrated_snn_res.0.3', y = 'project', color = y, ...) {
  setwd(OutDir)
  clUMAP(obj = obj, ident = x, save.plot = T, suffix = "as.in.barplot")
  categ.per.cluster <- ggbarplot(obj@meta.data
                                 , x = x
                                 , y = y
                                 , color = y
                                 , ...
  )
  qqSave(categ.per.cluster)
}

# get.clustercomposition()
# get.clustercomposition(, ylim=c(0,10000))


#  ------------------------------------------------------------------------
#  ------------------------------------------------------------------------
