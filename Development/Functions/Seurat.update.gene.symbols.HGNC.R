######################################################################
# Seurat.update.gene.symbols.HGNC.R
######################################################################
# source('~/GitHub/Packages/Seurat.utils/Functions/Seurat.update.gene.symbols.HGNC.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Seurat.update.gene.symbols.HGNC.R"))
require(HGNChelper)

# updateHGNC ------------------------------------------------------------------------------------
UpdateGenesSeurat <- function(obj = ls.Seurat[[i]], species_="human", EnforceUnique = T, ShowStats=F ) { # Update genes symbols that are stored in a Seurat object. It returns a data frame. The last column are the updated gene names.
  HGNC.updated <- HGNChelper::checkGeneSymbols(rownames(obj), unmapped.as.na = FALSE, map = NULL, species = species_)
  if (EnforceUnique) HGNC.updated <- HGNC.EnforceUnique(HGNC.updated)
  if (ShowStats) print(GetUpdateStats(HGNC.updated))
  obj <- RenameGenesSeurat(obj, newnames = HGNC.updated$Suggested.Symbol)
  return(obj)
}
# UpdateGenesSeurat()

# HELPER RenameGenesSeurat  ------------------------------------------------------------------------------------
RenameGenesSeurat <- function(obj = ls.Seurat[[i]], newnames = HGNC.updated[[i]]$Suggested.Symbol) { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
  RNA <- obj@assays$RNA

  if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
    if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- newnames
    # if (length(obj@meta.data)) rownames(obj@meta.data)          <- newnames
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}
# RenameGenesSeurat(obj = SeuratObj, newnames = HGNC.updated.genes$Suggested.Symbol)


# RemoveGenesSeurat ------------------------------------------------------------------------------------
RemoveGenesSeurat <- function(obj = ls.Seurat[[i]], symbols2remove = c("TOP2A")) { # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes metadata; obj@assays$RNA@counts, @data and @scale.data.
  print("Run this as the first thing after creating the Seurat object. It only removes genes from: metadata; obj@assays$RNA@counts, @data and @scale.data.")
  RNA <- obj@assays$RNA

  if (length(RNA@counts)) {
    NotFound <- setdiff(symbols2remove, RNA@counts@Dimnames[[1]])
    if (length(NotFound) == 0)  {
      RNA@counts@Dimnames[[1]] <- symbols2remove
      print("Genes removed from RNA@counts")
    } else {print("Not All Genes Found in RNA@counts. Missing:"); print(NotFound)}
  }
  if (length(RNA@data)) {
    if (length(setdiff(symbols2remove, RNA@data@Dimnames[[1]]) ) == 0)  {
      RNA@data@Dimnames[[1]] <- symbols2remove
      print("Genes removed from RNA@data.")
    } else {print("Not All Genes Found in RNA@data")}
  }
  if (length(RNA@scale.data)) {
    if (length(setdiff(symbols2remove, RNA@scale.data@Dimnames[[1]]) ) == 0)  {
      RNA@scale.data@Dimnames[[1]] <- symbols2remove
      print("Genes removed from RNA@scale.data.")
    } else {print("Not All Genes Found in RNA@scale.data")}
  }
  if (length(obj@meta.data)) {
    if (length(setdiff(symbols2remove, rownames(obj@meta.data)) ) == 0)  {
      rownames(obj@meta.data) <- symbols2remove
      print("Genes removed from @meta.data.")
    } else {print("Not All Genes Found in @metadata")}
  }
  obj@assays$RNA <- RNA
  return(obj)
}
# RemoveGenesSeurat(obj = SeuratObj, symbols2remove = "TOP2A")


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
  (UpdateStats = c("Updated (%)"=percentage_formatter(AcutallyUpdated / nrow(genes)), "Updated Genes"=floor(AcutallyUpdated), "Total Genes"=floor(nrow(genes))))
  return(UpdateStats)
}
# GetUpdateStats(genes = HGNC.updated.genes)

# update stats HGNC plot ------------------------------------------------------------------------------------
PlotUpdateStats <- function(mat = UpdateStatMat, column.names = c("Updated (%)",  "Updated (Nr.)")) { # Scatter plot of update stats.
  stopifnot(column.names %in% colnames(UpdateStatMat))
  HGNC.UpdateStatistics <- mat[, column.names]
  HGNC.UpdateStatistics[, "Updated (%)"] <- 100*HGNC.UpdateStatistics[, "Updated (%)"]
  colnames(HGNC.UpdateStatistics) <-  c("Gene Symbols updated (% of Total Genes)",  "Number of Gene Symbols updated")
  lll <- wcolorize(vector = rownames(HGNC.UpdateStatistics))
  wplot(HGNC.UpdateStatistics, col = lll
        , xlim = c(0,max(HGNC.UpdateStatistics[,1]))
        , ylim = c(0,max(HGNC.UpdateStatistics[,2])) )
  wlegend(NamedColorVec = lll, poz = 1)
}
# PlotUpdateStats(mat = result.of.GetUpdateStats)
