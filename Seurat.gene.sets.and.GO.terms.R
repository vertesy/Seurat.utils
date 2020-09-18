######################################################################
# Seurat.gene.sets.and.GO.terms.R
######################################################################
# source('~/GitHub/Seurat.utils/Seurat.gene.sets.and.GO.terms.R')

# require(MarkdownReports)
# source ('~/GitHub/CodeAndRoll/CodeAndRoll.R')

# Setup ------------------------------------------------------------
library(biomaRt)
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl") #uses human ensembl annotations


# ------------------------------------------------------------------------
PlotGoTermScores <- function(obj = combined.obj # Automate retrieving, processing and plotting GO term based gene scores.
                             , openBrowser = F
                             , GO = "GO:0061621", desc = "canonical.glycolysis") {
  GO.wDot<- make.names(GO)

  obj <- GetGOTerms(obj = obj, GO = GO, web.open = openBrowser);
  iprint(desc, obj@misc$GO[[ GO.wDot ]])
  obj <- AddGOScore(obj = obj, GO = GO);
  FeaturePlotSave(obj = obj, GO = paste0("Score.", GO.wDot), name_desc = desc)
  return(obj)
}
# "GO:0061621"  "canonical.glycolysis"
# PlotGoTermScores(GO = "GO:0061621", desc = "canonical.glycolysis")


# ------------------------------------------------------------------------
IntersectWithExpressed <- function(genes, obj=combined.obj) { # Intersect a set of genes with genes in the Seurat object.
  # print(head(genes, n=15))
  diff = setdiff(genes, rownames(obj))
  iprint(length(diff),"genes (of",length(genes), ") are MISSING from the Seurat object:",diff)
  return(intersect(rownames(obj), genes))
}
# GO.0010941.regulation.of.cell.death <- IntersectWithExpressed(GO.0010941.regulation.of.cell.death)

# ------------------------------------------------------------------------
GetGOTerms <- function(obj = combined.obj, GO = 'GO:0034976', web.open = T) { # Get GO terms via Biomart package
  genes <- getBM(attributes=c('hgnc_symbol'), #  'ensembl_transcript_id', 'go_id'
                 filters = "go_parent_term",  uniqueRows = TRUE,
                 values = GO, mart = ensembl)[,1]

  (GO.wDot<- make.names(GO))
  iprint(length(genes), "Gene symbols downloaded:", head(genes, n = 25))
  genes <- IntersectWithExpressed(obj = obj, genes = genes)

  if (is.null(obj@misc$GO)) obj@misc$GO <- list()
  obj@misc$GO[[ GO.wDot ]] <- genes
  iprint("Genes in", GO, "are saved under obj@misc$GO$", GO.wDot)
  if (web.open) system(paste0("open https://www.ebi.ac.uk/QuickGO/search/", GO))
  return(obj)
}
# combined.obj <- GetGOTerms(obj = combined.obj, GO = 'GO:0034976'); combined.obj@misc$GO$GO.0034976

# ------------------------------------------------------------------------
AddGOGeneList.manual <- function(obj = combined.obj, GO = 'GO:0034976', web.open=F  # Add GO terms via Biomart package.
                                 , genes =  c("A0A140VKG3", "ARX", "CNTN2", "DRD1", "DRD2", "FEZF2", "LHX6")) {
  print(head(genes, n =15))
  genes <- IntersectWithExpressed(obj = obj, genes = genes)

  if (is.null(obj@misc$GO)) obj@misc$GO <- list()
  obj@misc$GO[[make.names(GO)]] <- genes
  iprint("Genes in", GO, "are saved under obj@misc$GO$", make.names(GO))
  if (web.open) system(paste0("open https://www.ebi.ac.uk/QuickGO/search/", GO))
  return(obj)
}
# combined.obj <- AddGOGeneList.manual(obj = combined.obj, GO = 'GO:1904936'
#       , genes =  c("A0A140VKG3", "ARX", "CNTN2", "DRD1", "DRD2", "FEZF2", "LHX6")); combined.obj@misc$GO$GO.0034976


# ------------------------------------------------------------------------
AddGOScore <- function(obj = combined.obj, GO = "GO:0034976", FixName = TRUE ) { # Call after GetGOTerms. Calculates Score for gene set. Fixes name.
  GO.wDot<- make.names(GO)
  (genes.GO = list(obj@misc$GO[[GO.wDot]]))
  # print(genes.GO)
  (ScoreName = paste0("Score.", make.names(GO)))
  if (!is.list(genes.GO)) genes.GO<- list(genes.GO) # idk why this structure is not consistent...
  obj <- AddModuleScore(object = obj, features = genes.GO, name = ScoreName)

  if (FixName) {
    colnames(obj@meta.data) <-
      gsub(x = colnames(obj@meta.data)
           , pattern = paste0(ScoreName,1)
           , replacement = ScoreName
      )
    iprint("Trailing '1' in metadata column name is removed. Column name:", ScoreName)
  }
  return(obj)
}
# combined.obj <- AddGOScore(obj = combined.obj, GO = "GO:0034976", FixName = TRUE)
# combined.obj$Score.GO.0034976


# ------------------------------------------------------------------------
# name_desc="esponse to endoplasmic reticulum stress"
FeaturePlotSave <- function(obj = combined.obj, GO.score = "Score.GO.0034976", name_desc=NULL, h=7, PNG =T) { # Plot and save a FeaturePlot, e.g. showing gene set scores.
  proper.GO <- paste(sstrsplit(GO.score, pattern = "\\.", n = 3)[2:3], collapse = ":")
  (genes.GO = obj@misc$GO[[make.names(proper.GO)]])

  ggplot.obj <-
    FeaturePlot(obj, features = GO.score, min.cutoff = "q05", max.cutoff = "q95", reduction = 'umap') +
    labs(title = paste(GO.score, name_desc), caption = paste("Score calc. from",length(genes.GO), "expr. genes from BioMart.", paste0("https://www.ebi.ac.uk/QuickGO/search/", proper.GO)))
  pname = paste0("FeaturePlot.",(GO.score))
  fname = ww.FnP_parser(kpp(pname,name_desc), if (PNG) "png" else "pdf")
  save_plot(filename =fname, plot = ggplot.obj, base_height=h)
  ggplot.obj
}
# FeaturePlotSave()

# ------------------------------------------------------------------------
PasteUniqueGeneList <- function() {
  dput(sort(unique(clipr::read_clip())))
}


# ------------------------------------------------------------------------
CalcTranscriptomePercentage <- function(obj = combined.obj, genes = genes.GO.0061621.can.glyc) {
  total_expr = Matrix::colSums(GetAssayData(object = obj))
  Matrix::colSums(obj[ genes, ]) / total_expr
}

# genes.GO.0061621.canonical.glycolysis <- c("ENO1", "PKLR", "HK2", "HK3", "FOXK1", "PGAM2", "GCK", "BPGM",
#   "PGK1", "ALDOB", "PGM2L1", "PFKP", "HK1", "PGAM1", "GAPDH", "TPI1",
#   "ENO2", "PFKM", "PKM", "ADPGK", "ALDOA", "ENO3", "ALDOC", "FOXK2",
#   "GPI", "GAPDHS", "PFKL")
# Percentile.GO.0061621 <- CalcTranscriptomePercentage(genes = genes.GO.0061621.canonical.glycolysis)
# vioplot::vioplot(Percentile.GO.0061621)

# ------------------------------------------------------------------------
CalcTranscriptomePercentageGO <- function(obj = combined.obj, GO.score = "GO.0061621") {
  total_expr = Matrix::colSums(GetAssayData(object = obj))
  Matrix::colSums(obj[ obj@misc$GO[[GO.score]], ]) / total_expr
}

# Percentile.GO.0061621 <- CalcTranscriptomePercentageGO(GO.score = "GO.0061621")
# vioplot::vioplot(Percentile.GO.0061621)



# ------------------------------------------------------------------------
# xxGetGenesGo <- function(obj= combined.obj, GO = 'GO:0034976' ) obj@misc$GO[[  make.names(GO) ]]
