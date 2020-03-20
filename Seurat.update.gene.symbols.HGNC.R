######################################################################
# Seurat.update.gene.symbols.HGNC.R
######################################################################
# source ('~/GitHub/SeuratUtil/Seurat.update.gene.symbols.HGNC.R')


# updateHGNC helper ------------------------------------------------------------------------------------
RenameGenesSeurat <- function(SeuObj = ls.Seurat[[i]], newnames = HGNC.updated[[i]]) {
  print("Run this before integration. It only changes SeuObj@assays$RNA@counts, @data and @scale.data")
  RNA <- SeuObj@assays$RNA

  if (nrow(RNA) == nrow(newnames)) {
    if(l(RNA@counts)) RNA@counts@Dimnames[[1]] <-         newnames$Suggested.Symbol
    if(l(RNA@data)) RNA@data@Dimnames[[1]] <-             newnames$Suggested.Symbol
    if(l(RNA@scale.data)) RNA@scale.data@Dimnames[[1]] <- newnames$Suggested.Symbol
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  SeuObj@assays$RNA <- RNA
  return(SeuObj)
}

# updateHGNC ------------------------------------------------------------------------------------
UpdateGenesSeurat <- function(seu, species_="human") {
  HGNC.updated <- checkGeneSymbols(rownames(seu), unmapped.as.na = FALSE, map = NULL, species = species_)
  seu <- RenameGenesSeurat(seu, newnames = HGNC.updated)
}


# updateHGNC plot ------------------------------------------------------------------------------------
plot.UpdateStats <- function(genes = HGNC.updated[[i]]) {
  (MarkedAsUpdated <- genes[genes$Approved == FALSE, ])
  (AcutallyUpdated <- sum(MarkedAsUpdated[,1] != MarkedAsUpdated[,3]))
  (UpdateStats = c((AcutallyUpdated / nrow(genes)), AcutallyUpdated, nrow(genes)))
  return(UpdateStats)
}

