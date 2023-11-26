# ____________________________________________________________________
# Seurat.utils.less.used.R ----
# ____________________________________________________________________
# source("~/GitHub/Packages/Seurat.utils/R/Seurat.utils.less.used.R")




# _________________________________________________________________________________________________
#' @title Convert10Xfolders.old
#'
#' @description This function takes a parent directory with a number of subfolders, each containing the standard output of 10X Cell Ranger. It (1) loads the filtered data matrices, (2) converts them to Seurat objects, and (3) saves them as .RDS files.
#' @param InputDir A character string specifying the input directory.
#' @param folderPattern A character vector specifying the pattern of folder names to be searched. Default is 'filtered'.
#' @param min.cells An integer value specifying the minimum number of cells. Default is 10.
#' @param min.features An integer value specifying the minimum number of features. Default is 200.
#' @param updateHGNC A logical value indicating whether to update the HGNC. Default is TRUE.
#' @param ShowStats A logical value indicating whether to show statistics. Default is TRUE.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  Convert10Xfolders.old(InputDir = InputDir)
#'  }
#' }
#' @export
Convert10Xfolders.old <- function(InputDir
                                  , folderPattern = c("filtered", "SoupX_decont")[1]
                                  , min.cells = 10, min.features = 200
                                  , updateHGNC = TRUE, ShowStats = TRUE) {
  # ... function body ...
}

Convert10Xfolders.old <- function(InputDir # Take a parent directory with a number of subfolders, each containing the standard output of 10X Cell Ranger. (1.) It loads the filtered data matrices; (2.) converts them to Seurat objects, and (3.) saves them as *.RDS files.
                                  , folderPattern = c("filtered", "SoupX_decont")[1]
                                  , min.cells = 10, min.features = 200, updateHGNC = T, ShowStats = T) {
  fin <- list.dirs(InputDir, recursive = F)
  fin <- CodeAndRoll2::grepv(x = fin, pattern = folderPattern, perl = F)

  for (i in 1:length(fin)) {
    pathIN = fin[i]; print(pathIN)
    fnameIN = basename(fin[i])
    fnameOUT = ppp(paste0(InputDir, '/', fnameIN), 'min.cells', min.cells, 'min.features', min.features,"Rds")
    count_matrix <- Read10X(pathIN)

    if ( !is.list(count_matrix) | length(count_matrix) == 1) {
      seu <- CreateSeuratObject(counts = count_matrix, project = fnameIN,
                                min.cells = min.cells, min.features = min.features)
    } else if (is.list(count_matrix) & length(count_matrix) == 2)  {
      seu <- CreateSeuratObject(counts = count_matrix[[1]], project = fnameIN,
                                min.cells = min.cells, min.features = min.features)

      # LSB, Lipid Sample barcode (Multi-seq) --- --- --- --- --- ---
      LSB <- CreateSeuratObject(counts = count_matrix[[2]], project = fnameIN)
      LSBnameOUT = ppp(paste0(InputDir, '/LSB.', fnameIN),"Rds")
      saveRDS(LSB, file = LSBnameOUT)
    } else {
      print('More than 2 elements in the list of matrices')
    }
    # update --- --- --- ---
    if (updateHGNC) seu <- UpdateGenesSeurat(seu, EnforceUnique = T, ShowStats = T)
    saveRDS(seu, file = fnameOUT)
  }
}

