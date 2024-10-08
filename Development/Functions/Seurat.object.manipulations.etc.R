######################################################################
# Seurat.object.manipulations.etc.R
######################################################################
# source('~/GitHub/Packages/Seurat.utils/Functions/Seurat.object.manipulations.etc.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Seurat.object.manipulations.etc.R"))

# ------------------------------------------------------------------------
clip10Xcellname <- function(cellnames) str_split_fixed(cellnames, "_", n = 2)[,1] # Clip all suffices after underscore (10X adds it per chip-lane, Seurat adds in during integration).

# ------------------------------------------------------------------------
make10Xcellname <- function(cellnames, suffix="_1") paste0(cellnames, suffix) # Add a suffix to cell names, so that it mimics the lane-suffix, e.g.: "_1".


# ------------------------------------------------------------------------
seu.Make.Cl.Label.per.cell <- function(TopGenes, clID.per.cell) { # Take a named vector (of e.g. values ="gene names", names = clusterID), and a vector of cell-IDs and make a vector of "GeneName.ClusterID".
  Cl.names_class= TopGenes[ clID.per.cell ]
  Cl.names_wNr = paste0(Cl.names_class,' (',names(Cl.names_class),')')
  return(Cl.names_wNr)
}
# seu.Make.Cl.Label.per.cell(TopGenes = TopGenes.Classic,
#                            clID.per.cell = getMetadataColumn(ColName.metadata = metaD.CL.colname)  )

# FeaturePlot with different defaults ------------------------------------------------------------------
GetMostVarGenes <- function(obj=org, nGenes = p$nVarGenes) { # Get the most variable rGenes
  head(rownames(slot(object = obj, name = "hvg.info")), n = nGenes)
}

# gene.name.check for read .mtx /write .rds script ---------------------------------------
gene.name.check <- function(Seu.obj = ls.Seurat[[1]] ) { # Check gene names in a seurat object, for naming conventions (e.g.: mitochondrial reads have - or .). Use for reading .mtx & writing .rds files.
  rn = rownames(GetAssayData(object = Seu.obj, slot = "counts"))
  llprint("### Gene name pattern")

  llogit('`rn = rownames(GetAssayData(object = ls.Seurat[[1]], slot = "counts"))`')
  llogit('`head(grepv(rn, pattern = "-"), 10)`')
  print('pattern = -')
  llprint(head(grepv(rn, pattern = "-"), 10))

  llogit('`head(grepv(rn, pattern = "_"), 10)`')
  print('pattern = _')
  llprint(head(grepv(rn, pattern = "_"), 10))

  llogit('`head(grepv(rn, pattern = "\\."), 10)`')
  print('pattern = \\.')
  llprint(head(grepv(rn, pattern = "\\."), 10))

  llogit('`head(grepv(rn, pattern = "\\.AS[1-9]"), 10)`')
  print('pattern = \\.AS[1-9]')
  llprint(head(grepv(rn, pattern = "\\.AS[1-9]"), 10))
}


# check.genes ---------------------------------------
check.genes <- function(list.of.genes = ClassicMarkers, makeuppercase = FALSE, verbose  = TRUE, HGNC.lookup = FALSE
                        , obj = combined.obj
                        , assay.slot=c('RNA', 'integrated')[1]
                        , dataslot = c("counts", "data")[2]) { # Check if genes exist in your dataset.
  if (makeuppercase) list.of.genes <- toupper(list.of.genes)
  all_genes = rownames(GetAssayData(object = obj, assay = assay.slot, slot = dataslot)); length(all_genes)
  missingGenes = setdiff(list.of.genes, all_genes)
  if (length(missingGenes) > 0) {
    if (verbose) { iprint(l(missingGenes), "or", percentage_formatter(l(missingGenes) / l(list.of.genes)), "genes not found in the data, e.g:", head(missingGenes, n = 10))  }
    if (HGNC.lookup) {
      if (exists('qHGNC', mode='function')) { try(qHGNC(missingGenes)) } else { print("load qHGNC() function, see database.linker")}
    }
  }
  intersect(list.of.genes, all_genes)
}
# check.genes("top2a", makeuppercase = TRUE)
# check.genes("VGLUT2", verbose = F, HGNC.lookup = T)


# replace zero indexed clusternames ------------------------------------------------
fixZeroIndexing.seurat <- function(ColName.metadata = 'res.0.6', obj=org) { # Fix zero indexing in seurat clustering, to 1-based indexing
  obj@meta.data[ ,ColName.metadata] =  as.numeric(obj@meta.data[ ,ColName.metadata])+1
  print(obj@meta.data[ ,ColName.metadata])
  return(obj)
}


# CalculateFractionInTrome ------------------------------------------------
CalculateFractionInTrome <- function(geneset = c("MALAT1") # Calculate the fraction of a set of genes within the full Transcriptome of each cell.
                                     , obj = combined.obj
                                     , dataslot = c("counts", "data")[2]
) {
  print("    >>>> Use addMetaFraction() <<<<")
  geneset <- check.genes(geneset)
  stopifnot(length(geneset)>0)

  mat <- as.matrix(slot(obj@assays$RNA, name=dataslot))
  mat.sub <- mat[geneset,,drop=F]
  RC.per.cell.geneset <- colSums(mat.sub)

  RC.per.cell <- colSums(mat)
  gene.fraction.per.cell <- 100*RC.per.cell.geneset / RC.per.cell
  return(gene.fraction.per.cell)
}

# ------------------------------------------------------------------------
AddNewAnnotation <- function(obj = obj # Create a new metadata column based on an exisiting metadata column and a list of mappings (name <- IDs).
                             , source = "RNA_snn_res.0.5", named.list.of.identities = ls.Subset.ClusterLists) {
  NewID <- as.named.vector(obj[[source]])

  for (i in 1:length(named.list.of.identities)) {
    lx <- as.character(named.list.of.identities[[i]])
    name.lx <- names(named.list.of.identities)[i]
    NewID <- translate(vec = NewID, old = lx, new = name.lx)
  }
  print(table(NewID))
  return(NewID)
}
# ls.Subset.ClusterLists = list( "hESC.h9" = c("4", "10", "14"), "hESC.176" = c("0", "1", "2")); AddNewAnnotation()

# whitelist.subset.ls.Seurat ------------------------------------------------------------------------
whitelist.subset.ls.Seurat <- function(ls.obj = ls.Seurat
                                , metadir = p$'cellWhiteList' #  '~/Dropbox/Abel.IMBA/MetadataD/POL.meta/cell.lists/'
                                , whitelist.file = "NonStressedCellIDs.2020.10.21_18h.tsv"
) {
  cells.before <- unlapply(ls.obj, ncol)
  # Find file
  df.cell.whitelist <- read.simple.tsv(metadir, whitelist.file)
  dsets <- table(df.cell.whitelist[,1])

  ls.orig.idents <- lapply(lapply(ls.Seurat, getMetadataColumn, ColName.metadata = "orig.ident"), unique)
  stopif(any(unlapply(ls.orig.idents, l) == l(ls.Seurat)), message = "Some ls.Seurat objects have 1+ orig identity.")

  dsets.in.lsSeu <- unlist(ls.orig.idents)
  isMathced <- all(dsets.in.lsSeu == names(dsets)) # Stop if either ls.Seurat OR the metadata has identities not found in the other, in the same order.
  stopif(!isMathced, message = paste("either ls.Seurat OR the metadata has identities not found in the other, or they are not in same order."
                                     , kpps(dsets.in.lsSeu),"vs.", kpps(names(dsets) ) )
  )

  # identX <- ls.orig.idents[[1]]
  for (i in 1:l(ls.orig.idents)) {
    identX <- ls.orig.idents[[i]]; print(identX)

    # Extract and process cellIDs ----
    idx.match <- which(df.cell.whitelist[,1] == identX)
    cell.whitelist <- rownames(df.cell.whitelist)[idx.match]
    cell.whitelist <- substr(x = cell.whitelist
                             , start = 1 ,stop = nchar(cell.whitelist)-2)

    # Extract and process cellIDs ----
    ls.obj[[i]] <- subset(x = ls.obj[[i]], cells = cell.whitelist)
  }
  cells.after <- unlapply(ls.obj, ncol)
  iprint("cells.before",cells.before,"cells.after",cells.after)
  return(ls.obj)
}

# FindCorrelatedGenes ------------------------------------------------------------------------
FindCorrelatedGenes <- function(gene ="N.RabV.N2c", obj = combined.obj, assay = "RNA", slot = "data"
                                , HEonly =F , minExpr = 1, minCells = 1000
                                , trailingNgenes = 1000) {
  tic()
  AssayData <- GetAssayData(object = obj, assay = assay, slot = slot )
  matrix_mod <- iround(as.matrix(AssayData))
  if (HEonly) {
    idx.pass <- (matrixStats::rowSums2(matrix_mod > minExpr) > minCells)
    pc_TRUE(idx.pass)
    matrix_mod <- matrix_mod[ which(idx.pass), ]
  }
  geneExpr <- as.numeric(matrix_mod[gene, ])
  correlations <- apply(matrix_mod, 1, cor, geneExpr)
  topGenes <- trail(sort(correlations, decreasing = T), N = trailingNgenes)
  toc()
  wbarplot(head(topGenes, n =25))
  topGenes
}
# FindCorrelatedGenes(gene ="N.RabV.N2c", obj = combined.obj)
# write_clip(names(head(topGenes[-(1:6)], n=50)))


# ------------------------------------------------------------------------
Calc.Cor.Seurat <- function(assay.use = "RNA", slot.use = "data", geneset = FALSE
                            , quantileX = 0.95, max.cells =  10000, seed = p$"seed"
                            , digits = 2, obj = combined.obj) {
  expr.mat <- GetAssayData(slot = slot.use, assay = assay.use, object = obj)
  if (ncol(expr.mat) > max.cells) {
    set.seed(seed = seed)
    cells.use <- sample(x = colnames(expr.mat), size = max.cells)
  } else {
    cells.use <- colnames(obj)
  }

  qname = p0("q", quantileX * 100)
  quantile_name = kpp("expr", qname)
  if (is.null(obj@misc[[quantile_name]])) { iprint("Quantile data missing! Call: combined.obj <- calc.q90.Expression.and.set.all.genes(combined.obj, quantileX =",quantileX,") first!"); stop()}

  genes.HE  <- if (isFALSE(geneset)) {  which_names(obj@misc[[quantile_name]] > 0) } else {
                                        check.genes(geneset)  }
  iprint("Pearson correlation is calculated for", l(genes.HE), "HE genes with expr."
         , qname,": > 0 on a sample of", max.cells, " cells.")
  tic(); ls.cor <- sparse.cor(smat = t(expr.mat[genes.HE, cells.use])); toc()

  ls.cor <- lapply(ls.cor, round, digits = 2)

  slot__name <- kpp(slot.use, assay.use, quantile_name)
  obj@misc[[kpp('cor', slot__name)]] <- ls.cor$'cor'
  obj@misc[[kpp('cov', slot__name)]] <- ls.cor$'cov'
  iprint("Stored under obj@misc$", kpp('cor', slot.use, assay.use), "or cov... ." )
  return(obj)
}



# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
