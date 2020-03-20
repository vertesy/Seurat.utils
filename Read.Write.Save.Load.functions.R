######################################################################
# Read.Write.Save.Load.functions.R
######################################################################
# source ('~/GitHub/Seurat.utils/Read.Write.Save.Load.functions.R')

"Multicore read / write (I/O) functions are https://github.com/vertesy/Seurat.multicore"
"Single core read / write (I/O) functions are in https://github.com/vertesy/Seurat.utils/"


# Convert10Xfolders ------------------------------------------------------------------------
Convert10Xfolders <- function(InputDir # Take a parent directory with a number of subfolders, each containing the standard output of 10X Cell Ranger. (1.) It loads the filtered data matrices; (2.) converts them to Seurat objects, and (3.) saves them as *.RDS files.
  , min.cells=10, min.features=200, updateHGNC=T) {
  fin <- list.dirs(InputDir)[-1]
  for (i in 1:length(fin)) { print(fin[i])
    pathIN = fin[i]
    fnameIN = basename(fin[i])
    fnameOUT = ppp(p0(InputDir, 'filtered.', fnameIN), 'min.cells', min.cells, 'min.features', min.features,"Rds")
    x <- Read10X(pathIN)
    seu <- CreateSeuratObject(counts = x, project = fnameIN,
                              min.cells = min.cells, min.features = min.features)
    # update----
    if (updateHGNC) seu <- UpdateGenesSeurat(seu)
    saveRDS(seu, file = fnameOUT)
  }
}
# Convert10Xfolders(InputDir = InputDir)

# LoadAllSeurats ------------------------------------------------------------------------
LoadAllSeurats <- function(InputDir) { # Load a Seurat objects found in a directory. Also works with symbolic links (but not with aliases).
  fin <- list.files(InputDir, include.dirs = F, pattern = "*.Rds")
  ls.Seu <- list.fromNames(fin)
  for (i in 1:length(fin)) {print(fin[i]); ls.Seu[[i]] <- readRDS(p0(InputDir, fin[i]))}
  return(ls.Seu)
}
# ls.Seu <- LoadAllSeurats(InputDir = InputDir)

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
isave.RDS <- function(object, prefix =NULL, suffix=NULL, showMemObject=T, saveParams =T){ # Faster saving of workspace, and compression outside R, when it can run in the background. Seemingly quite CPU hungry and not veryefficient compression.
  path_rdata = paste0("~/Documents/RDS.files/", basename(OutDir))
  dir.create(path_rdata)

  if (showMemObject) { memory.biggest.objects() }
  if ( "seurat" %in% is(object) & saveParams) {
    try(object@misc$p <- p, silent = T)
    try(object@misc$all.genes  <- all.genes, silent = T)
  }
  fnameBase = kppu(prefix, substitute(object), suffix, idate())
  fname = MarkdownReportsDev::kollapse(path_rdata, "/",fnameBase , ".Rds")
  tictoc::tic()
  saveRDS(object, file = fname, compress=F)
  tictoc::toc()
  MarkdownReportsDev::iprint("Saved, being compressed", fname)
  say()
  system(paste("gzip", fname),  wait = FALSE) # execute in the background
}

# Save workspace -----------------------------------------------
# requires MarkdownReportsDev (github) and defining OutDir
# requires github/vertesy/CodeAndRoll.r

isave.image <- function(..., showMemObject=T, options=c("--force", NULL)[1]){ # Faster saving of workspace, and compression outside R, when it can run in the background. Seemingly quite CPU hungry and not veryefficient compression.
  path_rdata = paste0("~/Documents/Rdata.files/", basename(OutDir))
  dir.create(path_rdata)

  if (showMemObject) { memory.biggest.objects() }
  fname = MarkdownReportsDev::kollapse(path_rdata, "/",idate(),...,".Rdata")
  print(fname)
  if (nchar(fname) > 2000) stop()
  save.image( file = fname, compress=F)
  MarkdownReportsDev::iprint("Saved, being compressed", fname)
  system(paste("gzip", options, fname),  wait = FALSE) # execute in the background
}

# ------------------------------------------------------------------------

subsetSeuObj.and.Save <- function(obj=ORC, fraction = 0.25 ) { # Subset a compressed Seurat Obj and save it in wd.
  cellIDs.keep = sampleNpc(metaDF = obj@meta.data, pc = fraction)

  obj_Xpc <- subset(obj, cells = cellIDs.keep) # downsample
  saveRDS(obj_Xpc, compress = TRUE,
          file = ppp(p0(InputDir, 'seu.ORC'), length(cellIDs.keep), 'cells.with.min.features', p$min.features,"Rds" ) )
  say()
}

