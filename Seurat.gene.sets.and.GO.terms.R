######################################################################
# Seurat.gene.sets.and.GO.terms.R
######################################################################
# source ('~/GitHub/Seurat.utils/Seurat.gene.sets.and.GO.terms.R')

# require(MarkdownReports)
# source ('~/GitHub/CodeAndRoll/CodeAndRoll.R')

# Setup ------------------------------------------------------------
library(biomaRt)
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl") #uses human ensembl annotations

# ------------------------------------------------------------------------

IntersectWithExpressed <- function(genes, obj=combined.obj) { # Intersect a set of genes with genes in the Seurat object.
  # print(head(genes, n=15))
  diff = setdiff(genes, rownames(obj))
  iprint(length(diff),"genes (of",l(genes), ") are MISSING from the Seurat object:",diff)
  return(intersect(rownames(obj), genes))
}
# GO.0010941.regulation.of.cell.death <- IntersectWithExpressed(GO.0010941.regulation.of.cell.death)

# ------------------------------------------------------------------------
GetGOTerms <- function(obj = combined.obj, GO = 'GO:0034976', web.open = T) { # Get GO terms via Biomart package
  genes <- getBM(attributes=c('hgnc_symbol'), #  'ensembl_transcript_id', 'go_id'
                 filters = "go",  uniqueRows = TRUE,
                 values = GO, mart = ensembl)[,1]
  iprint("Gene symbols downloaded:", genes)
  genes <- IntersectWithExpressed(obj = obj, genes = genes)

  if (is.null(obj@misc$GO)) obj@misc$GO <- list()
  obj@misc$GO[[make.names(GO)]] <- genes
  iprint("Genes in", GO, "are saved under obj@misc$GO$", make.names(GO))
  if (web.open) system(paste0("open https://www.ebi.ac.uk/QuickGO/search/", GO))
  return(obj)
}
# combined.obj <- GetGOTerms(obj = combined.obj, GO = 'GO:0034976'); combined.obj@misc$GO$GO.0034976

# ------------------------------------------------------------------------
AddGOGeneList.manual <- function(obj = combined.obj, GO = 'GO:0034976', web.open=F  # Add GO terms via Biomart package
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
AddGOScore <- function(obj = combined.obj, GO = "GO:0034976", FixName = TRUE ) { # Call after GetGOTerms. Calculates Score for gene set. Fixes name
  (genes.GO = list(obj@misc$GO[make.names(GO)]))
  (ScoreName = paste0("Score.",make.names(GO)))
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
FeaturePlotSave <- function(obj = combined.obj, GO = "Score.GO.0034976", name_desc=NULL, h=7, PNG =F) { # Plot and save a FeaturePlot, e.g. showing gene set scores.
  ggplot.obj <-
    FeaturePlot(obj, features = make.names(GO)
                , min.cutoff = "q05", max.cutoff = "q95", reduction = 'umap') +
    labs(title = paste(GO, name_desc))
  pname = paste0("FeaturePlot.",make.names(GO))
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

# ------------------------------------------------------------------------


