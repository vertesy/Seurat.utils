######################################################################
# Read.Write.Save.Load.functions.R
######################################################################
# source('~/GitHub/Packages/Seurat.utils/Read.Write.Save.Load.functions.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Read.Write.Save.Load.functions.R"))

"Multicore read / write (I/O) functions are https://github.com/vertesy/Seurat.multicore"
"Single core read / write (I/O) functions are in https://github.com/vertesy/Seurat.utils/"

# Convert10Xfolders ------------------------------------------------------------------------
Convert10Xfolders <- function(InputDir # Take a parent directory with a number of subfolders, each containing the standard output of 10X Cell Ranger. (1.) It loads the filtered data matrices; (2.) converts them to Seurat objects, and (3.) saves them as *.RDS files.
                              , folderPattern = c("filtered", "SoupX_decont")[1]
                              , min.cells=10, min.features=200, updateHGNC=T, ShowStats=T) {
  fin <- list.dirs(InputDir, recursive = F)
  fin <- grepv(x = fin, pattern = folderPattern, perl = F)

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

      # LSB, Lipid Sample barcode (Multi-seq) --------------------
      LSB <- CreateSeuratObject(counts = count_matrix[[2]], project = fnameIN)
      LSBnameOUT = ppp(paste0(InputDir, '/LSB.', fnameIN),"Rds")
      saveRDS(LSB, file = LSBnameOUT)
    } else {
      print('More than 2 elements in the list of matrices')
    }
    # update----
    if (updateHGNC) seu <- UpdateGenesSeurat(seu, EnforceUnique = T, ShowStats = T)
    saveRDS(seu, file = fnameOUT)
  }
}
# Convert10Xfolders(InputDir = InputDir)

# ConvertDropSeqfolders ------------------------------------------------------------------------
ConvertDropSeqfolders <- function(InputDir # Take a parent directory with a number of subfolders, each containing the standard output of 10X Cell Ranger. (1.) It loads the filtered data matrices; (2.) converts them to Seurat objects, and (3.) saves them as *.RDS files.
                                  , folderPattern = "SRR*", filePattern = "expression.tsv.gz"
                                  , useVroom = T, col_types.vroom = list("GENE" = "c", .default = "d")
                                  , min.cells=10, min.features=200, updateHGNC=T, ShowStats=T, minDimension = 10, overwrite = FALSE) {
  InputDir <- FixPath(InputDir)
  fin <- list.dirs(InputDir, recursive = F)
  fin <- grepv(x = fin, pattern = folderPattern, perl = F)

  for (i in 1:length(fin)) { print(i)
    pathIN <- FixPath(fin[i]); print(pathIN)
    fnameIN <- basename(fin[i])
    subdir <- p0(InputDir, fnameIN)
    fnameOUT <- ppp(subdir, 'min.cells', min.cells, 'min.features', min.features,"Rds"); print(fnameOUT)
    if (!overwrite) {
      OutFile <- list.files(InputDir, pattern = basename(fnameOUT), recursive = T)
      if (length(OutFile) > 0) {
        if (grepl(pattern = ".Rds$", OutFile, perl = T)) {
          iprint("      RDS OBJECT ALREADY EXISTS.");
          next
        }
      } # if length
    }
    CountTable <- list.files(subdir, pattern = filePattern,recursive = F)
    stopifnot(length(CountTable) == 1 )
    count_matrix <- if (useVroom) {
      vroom::vroom(file = kpps(subdir, CountTable), col_types = col_types.vroom)
    } else {
      readr::read_tsv(file = kpps(subdir, CountTable))
    }

    if (nrow(count_matrix) < minDimension | ncol(count_matrix) < minDimension ) {
      iprint(""); iprint("      EXPRESSION MATRIX TOO SMALL.", nrow(count_matrix), "x", ncol(count_matrix),". Not processed.");
    } else {
      count_matrix <- FirstCol2RowNames(count_matrix)[,-1] # remove 1st "Cell column" # https://github.com/vertesy/SEO/issues/63
      seu <- CreateSeuratObject(counts = count_matrix, project = fnameIN,
                                min.cells = min.cells, min.features = min.features)
      if (ncol(seu) < 1000) print("Only", ncol(seu), "cells survived filtering in the Seurat obj!")
      if (nrow(seu) < 1000) print("Only", nrow(seu), "genes found in the Seurat obj!")

      # update HGNC ----
      Sys.setenv('R_MAX_VSIZE' = 32000000000)
      if (updateHGNC) seu <- UpdateGenesSeurat(seu, EnforceUnique = T, ShowStats = T)
      saveRDS(seu, file = fnameOUT)
    }
  }
}
# ConvertDropSeqfolders(InputDir)

# LoadAllSeurats ------------------------------------------------------------------------
LoadAllSeurats <- function(InputDir # Load all Seurat objects found in a directory. Also works with symbolic links (but not with aliases).
                           , file.pattern = "^filtered.+Rds$"
                           , string.remove1 = c(F, "filtered_feature_bc_matrix.", "raw_feature_bc_matrix." )[2]
                           , string.remove2 = c(F, ".min.cells.10.min.features.200.Rds")[2]) {
  tic()
  InputDir <- AddTrailingSlash(InputDir) # add '/' if necessary

  fin.orig <- list.files(InputDir, include.dirs = F, pattern = file.pattern)
  print(fin.orig)
  fin <- if (!isFALSE(string.remove1)) sapply(fin.orig, gsub, pattern = string.remove1, replacement = "") else fin.orig
  fin <- if (!isFALSE(string.remove2)) sapply(fin, gsub, pattern = string.remove2, replacement = "") else fin

  ls.Seu <- list.fromNames(fin)
  for (i in 1:length(fin)) {print(fin[i]); ls.Seu[[i]] <- readRDS(paste0(InputDir, fin.orig[i]))}
  print(toc())
  return(ls.Seu)
}
# ls.Seurat <- LoadAllSeurats(InputDir)


# ------------------------------------------------------------------------------------------------
read10x <- function(dir) { # read10x from gzipped matrix.mtx, features.tsv and barcodes.tsv
  tictoc::tic()
  names <- c("barcodes.tsv", "features.tsv", "matrix.mtx")
  for (i in 1:length(names)) {
    R.utils::gunzip(paste0(dir, "/", names[i], ".gz"))
  }
  file.copy(paste0(dir, "/features.tsv"), paste0(dir, "/genes.tsv"))
  mat <- Seurat::Read10X(dir)
  file.remove(paste0(dir, "/genes.tsv"))
  for (i in 1:length(names)) {
    R.utils::gzip(paste0(dir, "/", names[i]))
  }
  tictoc::toc()
  mat
}

#### Functions in Saving.and.loading.R

# Save an object -----------------------------------------------
isave.RDS <- function(object, prefix =NULL, suffix=NULL, showMemObject=T, saveParams =T, inOutDir = F){ # Faster saving of workspace, and compression outside R, when it can run in the background. Seemingly quite CPU hungry and not very efficient compression.
  path_rdata = if (inOutDir) { OutDir } else { paste0( "~/Documents/RDS.files/", basename(OutDir) ) }
  dir.create(path_rdata)

  if (showMemObject) { memory.biggest.objects() }
  if ( "seurat" %in% is(object) & saveParams) {
    try(object@misc$p <- p, silent = T)
    try(object@misc$all.genes  <- all.genes, silent = T)
  }
  fnameBase = kppu(prefix, substitute(object), suffix, idate(Format = "%Y.%m.%d_%H.%M"))
  fnameBase = trimws(fnameBase, whitespace = '_')
  fname = MarkdownReportsDev::kollapse(path_rdata, "/",fnameBase , ".Rds")
  tictoc::tic()
  saveRDS(object, file = fname, compress = F)
  tictoc::toc()
  MarkdownReportsDev::iprint("Saved, being compressed", fname)
  say()
  system(paste("gzip", fname),  wait = FALSE) # execute in the background
}

# Save workspace -----------------------------------------------
# requires MarkdownReportsDev (github) and defining OutDir
# requires github/vertesy/CodeAndRoll.r

isave.image <- function(..., showMemObject=T, options=c("--force", NULL)[1]){ # Faster saving of workspace, and compression outside R, when it can run in the background. Seemingly quite CPU hungry and not very efficient compression.
  path_rdata = paste0("~/Documents/Rdata.files/", basename(OutDir))
  dir.create(path_rdata)

  if (showMemObject) { memory.biggest.objects() }
  fname = MarkdownReportsDev::kollapse(path_rdata, "/",idate(),...,".Rdata")
  print(fname)
  if (nchar(fname) > 2000) stop()
  save.image(file = fname, compress = F)
  MarkdownReportsDev::iprint("Saved, being compressed", fname)
  system(paste("gzip", options, fname),  wait = FALSE) # execute in the background
}


# Save workspace -----------------------------------------------
# requires MarkdownReportsDev (github) and defining OutDir
# requires github/vertesy/CodeAndRoll.r

qsave.image <- function(..., showMemObject=T, options=c("--force", NULL)[1]){ # Faster saving of workspace, and compression outside R, when it can run in the background. Seemingly quite CPU hungry and not very efficient compression.
  fname = MarkdownReportsDev::kollapse(getwd(), "/",basename(OutDir),idate(),...,".Rdata")
  print(fname)
  if (nchar(fname) > 2000) stop()
  tic()
  save.image(file = fname, compress = F)
  MarkdownReportsDev::iprint("Saved, being compressed", fname)
  system(paste("gzip", options, fname),  wait = FALSE) # execute in the background
  cat(toc)
}

# subsetSeuObj -----------------------------------------------------------------------
subsetSeuObj <- function(obj=ls.Seurat[[i]], fraction_ = 0.25, nCells = F, seed_ = 1989 ) { # Subset a compressed Seurat Obj and save it in wd.
  set.seed(seed_)
  if (isFALSE(nCells)) {
    cellIDs.keep = sampleNpc(metaDF = obj@meta.data, pc = fraction_)
    iprint(length(cellIDs.keep), "or",percentage_formatter(fraction_),"of the cells are kept. Seed:", seed_)
  } else if (nCells > 1) {
    nKeep = min(ncol(obj), nCells)
    # print(nKeep)
    cellIDs.keep = sample(colnames(obj), size = nKeep, replace = F)
    if (nKeep < nCells) iprint("Only",nCells,"cells were found in the object, so downsampling is not possible.")
  }
  obj <- subset(x = obj, cells = cellIDs.keep) # downsample
  return(obj)
}

# subsetSeuObj.and.Save ------------------------------------------------------------------------
subsetSeuObj.and.Save <- function(obj=ORC, fraction = 0.25, seed = 1989 ) { # Subset a compressed Seurat Obj and save it in wd.
  obj_Xpc <- subsetSeuObj(obj = obj, fraction_ =  fraction, seed_ = seed)
  saveRDS(obj_Xpc, compress = TRUE,
          file = ppp(paste0(InputDir, 'seu.ORC'), length(cellIDs.keep), 'cells.with.min.features', p$min.features,"Rds" ) )
  say()
}


# Downsample.Seurat.Objects ------------------------------------------------------------------------
Downsample.Seurat.Objects <- function(ls.obj = ls.Seurat, NrCells = p$"dSample.Organoids") {
  names.ls = names(ls.obj)
  n.datasets = length(ls.obj)
  iprint(NrCells, "cells")
  tic()
  if (getDoParRegistered() ) {
    ls.obj.downsampled <- foreach(i = 1:n.datasets ) %dopar% {
      iprint(names(ls.obj)[i], percentage_formatter(i/n.datasets, digitz = 2))
      subsetSeuObj(obj = ls.obj[[i]], nCells = NrCells)
    }; names(ls.obj.downsampled)  <- names.ls
  } else {
    ls.obj.downsampled <- list.fromNames(names.ls)
    for (i in 1:n.datasets ) {
      iprint(names(ls.obj)[i], percentage_formatter(i/n.datasets, digitz = 2))
      ls.obj.downsampled[[i]] <- subsetSeuObj(obj = ls.obj[[i]], nCells = NrCells)
    };
  } # else
  toc();

  print(head(unlapply(ls.obj, ncol)))
  print(head(unlapply(ls.obj.downsampled, ncol)))

  isave.RDS(object = ls.obj.downsampled, suffix = ppp(NrCells, "cells"), inOutDir = T)

}
# Downsample.Seurat.Objects(NrCells = 2000)
# Downsample.Seurat.Objects(NrCells = 1000)
# Downsample.Seurat.Objects(NrCells = 500)
# Downsample.Seurat.Objects(NrCells = 200)


