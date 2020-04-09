######################################################################
# Seurat.gene.sets.and.GO.terms.R
######################################################################
# source ('~/GitHub/Seurat.utils/Seurat.gene.sets.and.GO.terms.R')



# ------------------------------------------------------------------------

IntersectWithExpressed <- function(genes, obj=combined.obj) { # Intersect a set of genes with genes in the Seurat object.
  diff = setdiff(genes, rownames(obj))
  iprint(length(diff),"genes (of",l(genes), ") are not found in the Seurat object:",diff)
  return(intersect(rownames(obj), genes))
}
# GO.0010941.regulation.of.cell.death <- IntersectWithExpressed(GO.0010941.regulation.of.cell.death)

# ------------------------------------------------------------------------

GetGOTerms <- function(obj = combined.obj, GO = 'GO:0034976', web.open = T) { # Get GO terms via Biomart package
  genes <- getBM(attributes=c('hgnc_symbol'), #  'ensembl_transcript_id', 'go_id'
                 filters = "go",  uniqueRows = TRUE,
                 values = GO, mart = ensembl)[,1]
  genes <- IntersectWithExpressed(obj = obj, genes = genes)

  if (is.null(obj@misc$GO)) obj@misc$GO <- list()
  obj@misc$GO[[make.names(GO)]] <- genes
  iprint("Genes in", GO, "are saved under obj@misc$GO$", make.names(GO))
  if (web.open) system(paste0("open https://www.ebi.ac.uk/QuickGO/search/", GO))
  return(obj)
}
# combined.obj <- GetGOTerms(obj = combined.obj, GO = 'GO:0034976'); combined.obj@misc$GO$GO.0034976


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
    print("Trailing '1' in metadata column name is removed")
  }
  return(obj)
}
# combined.obj <- AddGOScore(obj = combined.obj, GO = "GO:0034976", FixName = TRUE)
# combined.obj$Score.GO.0034976


# ------------------------------------------------------------------------
FeaturePlotSave <- function(obj = combined.obj, GO = "Score.GO.0034976", h=7, PNG =F) {
  ggplot.obj <-
    FeaturePlot(obj, features = make.names(GO)
                , min.cutoff = "q10", max.cutoff = "q90", reduction = 'umap')
  pname = paste0("FeaturePlot.",make.names(GO))
  fname = ww.FnP_parser(pname, if (PNG) "png" else "pdf")
  save_plot(filename =fname, plot = ggplot.obj, base_height=h)
  ggplot.obj
}
# FeaturePlotSave()

# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------
