######################################################################
# Seurat.update.gene.symbols.HGNC.R
######################################################################
# source ('~/GitHub/Seurat.utils/Seurat.update.gene.symbols.HGNC.R')
require(HGNChelper)


# updateHGNC ------------------------------------------------------------------------------------
UpdateGenesSeurat <- function(seu, species_="human", EnforceUnique = T, ShowStats=F ) { # Update genes symbols that are stored in a Seurat object. It returns a data frame. The last column are the updated gene names.
  HGNC.updated <- HGNChelper::checkGeneSymbols(rownames(seu), unmapped.as.na = FALSE, map = NULL, species = species_)
  if (EnforceUnique) HGNC.updated <- HGNC.EnforceUnique(HGNC.updated)
  if (ShowStats) print(GetUpdateStats(HGNC.updated))
  seu <- RenameGenesSeurat(seu, newnames = HGNC.updated)
  return(seu)
}

# HELPER updateHGNC  ------------------------------------------------------------------------------------
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

# HELPER Enforce Unique names ------------------------------------------------------------------------------------
HGNC.EnforceUnique <- function(updatedSymbols) { # Enforce Unique names after HGNC symbol update. updatedSymbols is the output of HGNChelper::checkGeneSymbols.
    NGL <- updatedSymbols[,3]
    if (any.duplicated(NGL)) {
      updatedSymbols[,3] <- make.unique(NGL); "Unique names are enforced by suffixing .1, .2, etc."
      }
    return(updatedSymbols)
}
# x <- HGNC.EnforceUnique(updatedSymbols = SymUpd)
# While "make.unique" is not the ideal solution, because it generates mismatched, in my integration example it does reduce the mismatching genes from ~800 to 4


# update stats HGNC  ------------------------------------------------------------------------------------
GetUpdateStats <- function(genes = HGNC.updated[[i]]) { # Plot the Symbol-update statistics. Works on the data frame returned by `UpdateGenesSeurat()`.
  (MarkedAsUpdated <- genes[genes$Approved == FALSE, ])
  (AcutallyUpdated <- sum(MarkedAsUpdated[,1] != MarkedAsUpdated[,3]))
  (UpdateStats = c("Updated (%)"=(AcutallyUpdated / nrow(genes)), "Updated Genes"=AcutallyUpdated, "Total Genes"=nrow(genes)))
  return(UpdateStats)
}

# update stats HGNC plot ------------------------------------------------------------------------------------
PlotUpdateStats <- function(mat = UpdateStatMat) { # Scatter plot of update stats.
  HGNC.UpdateStatistics <- mat[, c("Updated (%)",  "Updated (Nr.)") ]
  HGNC.UpdateStatistics[, "Updated (%)"] <- 100*HGNC.UpdateStatistics[, "Updated (%)"]
  colnames(HGNC.UpdateStatistics) <-  c("Gene Symbols updated (% of Total Genes)",  "Number of Gene Symbols updated")
  lll <- wcolorize(vector = rownames(HGNC.UpdateStatistics))
  wplot(HGNC.UpdateStatistics, col = lll
        , xlim = c(0,max(HGNC.UpdateStatistics[,1]))
        , ylim = c(0,max(HGNC.UpdateStatistics[,2])) )
  wlegend(NamedColorVec = lll, poz = 1)
}
