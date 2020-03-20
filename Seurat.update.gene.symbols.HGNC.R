######################################################################
# Seurat.update.gene.symbols.HGNC.R
######################################################################
# source ('~/GitHub/Seurat.utils/Seurat.update.gene.symbols.HGNC.R')


# updateHGNC helper ------------------------------------------------------------------------------------
RenameGenesSeurat <- function(SeuObj = ls.Seurat[[i]], newnames = HGNC.updated[[i]]) { # Replace gene names in different slots of a Seurat object. Run this before integration. It only changes SeuObj@assays$RNA@counts, @data and @scale.data.
  print("Run this before integration. It only changes SeuObj@assays$RNA@counts, @data and @scale.data")
  RNA <- SeuObj@assays$RNA

  if (nrow(RNA) == nrow(newnames)) {
    if(length(RNA@counts)) RNA@counts@Dimnames[[1]] <-         newnames$Suggested.Symbol
    if(length(RNA@data)) RNA@data@Dimnames[[1]] <-             newnames$Suggested.Symbol
    if(length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]] <- newnames$Suggested.Symbol
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  SeuObj@assays$RNA <- RNA
  return(SeuObj)
}

# updateHGNC ------------------------------------------------------------------------------------
UpdateGenesSeurat <- function(seu, species_="human") { # Update genes symbols that are stored in a Seurat object. It returns a data frame. The last column are the updated gene names.
  HGNC.updated <- checkGeneSymbols(rownames(seu), unmapped.as.na = FALSE, map = NULL, species = species_)
  seu <- RenameGenesSeurat(seu, newnames = HGNC.updated)
}


# updateHGNC plot ------------------------------------------------------------------------------------
plot.UpdateStats <- function(genes = HGNC.updated[[i]]) { # Plot the Update Statistcs. Works on the data frame returned by `UpdateGenesSeurat()`.
  (MarkedAsUpdated <- genes[genes$Approved == FALSE, ])
  (AcutallyUpdated <- sum(MarkedAsUpdated[,1] != MarkedAsUpdated[,3]))
  (UpdateStats = c((AcutallyUpdated / nrow(genes)), AcutallyUpdated, nrow(genes)))
  return(UpdateStats)
}

