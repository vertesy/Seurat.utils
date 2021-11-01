require(Seurat)
InputDir = "/Volumes/SEO/Geschwind_Fetal/dbGaP-27527.dropseq"
require(tictoc)

# ConvertDropSeqfolders ------------------------------------------------------------------------
ConvertDropSeqfolders <- function(InputDir # Take a parent directory with a number of subfolders, each containing the standard output of 10X Cell Ranger. (1.) It loads the filtered data matrices; (2.) converts them to Seurat objects, and (3.) saves them as *.RDS files.
                                  , folderPattern = "SRR*", filePattern = "expression.tsv.gz"
                                  , min.cells=10, min.features=200, updateHGNC=T, ShowStats=T) {
  fin <- list.dirs(InputDir, recursive = F)
  fin <- grepv(x = fin, pattern = folderPattern, perl = F)

  for (i in 1:length(fin)) {
    pathIN <- fin[i]; print(pathIN)
    fnameIN <- basename(fin[i])
    subdir <- kpps(InputDir, fnameIN)
    fnameOUT = ppp(subdir, 'min.cells', min.cells, 'min.features', min.features,"Rds"); print(fnameOUT)

    CountTable <- list.files(subdir, pattern = filePattern,recursive = F)
    stopifnot(length(CountTable) == 1 )

    fileX <- kpps(subdir, CountTable)
    file.exists(fileX)
    Sys.setenv("VROOM_CONNECTION_SIZE" = 2^20)
    tic()
    count_matrix <- vroom::vroom(file = fileX)
    toc()

    tic()
    count_matrix2 <- readr::read_tsv(file=fileX )
    toc()

    count_matrix <- FirstCol2RowNames(count_matrix)[,-1] # remove 1st "Cell column" # https://github.com/vertesy/SEO/issues/63
    seu <- CreateSeuratObject(counts = count_matrix, project = fnameIN,
                              min.cells = min.cells, min.features = min.features)
    # update----
    if (updateHGNC) seu <- UpdateGenesSeurat(seu, EnforceUnique = T, ShowStats = T)
    saveRDS(seu, file = fnameOUT)
  }
}

