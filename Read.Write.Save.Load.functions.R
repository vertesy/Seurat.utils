######################################################################
# Read.Write.Save.Load.functions.R
######################################################################
# source ('~/GitHub/SeuratUtil/Read.Write.Save.Load.functions.R')

# Convert10Xfolders ------------------------------------------------------------------------
Convert10Xfolders <- function(InputDir # Take a parent directory with a number of subfolders, each containing the standard output of 10X Cell Ranger. (1.) It loads the filtered data matrices; (2.) converts them to Seurat objects, and (3.) saves them as *.RDS files.
  , min.cells=10, min.features=200, updateHGNC=T) {
  fin <- list.dirs(InputDir)[-1]
  for (i in 1:l(fin)) { print(fin[i])
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
  for (i in 1:l(fin)) {print(fin[i]); ls.Seu[[i]] <- readRDS(p0(InputDir, fin[i]))}
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


# Save an object -----------------------------------------------
isave.RDS.pigz <- function(object, prefix =NULL, suffix=NULL, showMemObject=T, saveParams =T){ # Faster saving of workspace, and compression outside R, when it can run in the background. Seemingly quite CPU hungry and not veryefficient compression.
  path_rdata = paste0("~/Documents/RDS.files/", basename(OutDir))
  dir.create(path_rdata)

  if ( "seurat" %in% is(object) & saveParams) {
    try(object@misc$p <- p, silent = T)
    try(object@misc$all.genes  <- all.genes, silent = T)
  }
  if (showMemObject) { memory.biggest.objects() }
  fnameBase = kppu(prefix, substitute(object), suffix, idate())
  fname = MarkdownReportsDev::kollapse(path_rdata, "/",fnameBase , ".Rds")
  tictoc::tic()
  ssaveRDS(object = object, file = fname)
  tictoc::toc()
  say()
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
          file = ppp(p0(InputDir, 'seu.ORC'), l(cellIDs.keep), 'cells.with.min.features', p$min.features,"Rds" ) )
  say()
}


# ------------------------------------------------------------------------

seuSaveRds <- function(object = ls.Seurat, tags = setupFlags, use_Original_OutDir = F) { # Save a compressed Seurat Object, with parallel gzip by pgzip
  if (use_Original_OutDir) create_set_Original_OutDir()
  fname.comb.rds = ppp(substitute(object), tags, idate(), ".Rds")
  iprint( fname.comb.rds)
  ssaveRDS(object = object, filename = fname.comb.rds)
  say()
}

# ------------------------------------------------------------------------

sampleNpc <- function(metaDF = MetaData[which(Pass),], pc=0.1) { # Sample N % of a dataframe (obj@metadata), and return the cell IDs.
  cellIDs = rownames(metaDF)
  nr_cells = floor(l(cellIDs) * pc)
  cellIDs.keep = sample(cellIDs, size = nr_cells, replace = FALSE)
  return(cellIDs.keep)
}

# Wrapper layer 2 (top) -----------------------------------------------------------------------------------

# read / load multiple objects
rrRDS <- function(list_of_objectnames = c("ls.Seurat", "ls2", "org"), ...) { # Load a list of RDS files with parallel ungzip by pgzip.
  tictoc::tic()
  path_rdata = paste0("~/Documents/RDS.files/", basename(OutDir))
  iprint("Looking for files under:" , path_rdata)
  iprint("(Base on OutDir)")
  dir.exists(path_rdata)
  # lsX <- foreach(obj = list_of_objectnames) %dopar% {
  for (obj in list_of_objectnames) {
    fname = MarkdownReportsDev::kollapse(path_rdata , "/", obj)
    tmp.obj = rreadRDS(filename = fname, ..., cores=2)
    # pigzx multi-core decompression is not so useful, with 2 cores ist actually faster:

    assigned.name = strsplit(obj, split = '\\.201*')[[1]][1]
    print(paste(" --> Loaded as object: ",assigned.name))
    assign(x = assigned.name,value = tmp.obj, envir = .GlobalEnv)
  }
  tictoc::toc()
}


# require(clipr)
# clip2clip.vector()
# x = c("org.2019.03.12_10h.Rds", "ls2.2019.03.12_10h.Rds", "ls.Seurat.2019.03.12_10h.Rds")
# rrRDS(x)


# Save multiple objects using pigz by default ---------------------------------------------
sssRDS <- function(list_of_objectnames = c("ls.Seurat", "ls2", "org.ALL", "org"), name.suffix =NULL, ...) { # Save multiple objects into a list of RDS files using parallel gzip by pgzip (optional).
  tictoc::tic()
  base_name <- character()
  path_rdata = paste0("~/Documents/RDS.files/", basename(OutDir))
  iprint("Files saved under:" , path_rdata)

  try(dir.create(path_rdata))
  for (obj in list_of_objectnames) {
    print(obj)
    if ( "seurat" %in% is(obj)) {
      obj@misc$p = p
      obj@misc$all.genes = all.genes
    }
    base_name[i] = paste0(obj, '.', name.suffix, '.', idate(),".Rds", collapse = "")
    fname = paste0(path_rdata , "/", base_name[i], collapse = "")
    ssaveRDS( object = get(obj), filename = fname, ...)
  }
  tictoc::toc()
  dput(base_name)
}

### Wrapper layer 1 -----------------------------------------------------------------------------------
ssaveRDS <- function(object, filename, con_func = list(pigz_pipe, snappy_pipe)[[1]], func_type = c("pipe", "builtin")[1], ...) { # Save an object with parallel gzip by pgzip.
  tictoc::tic()
  if(func_type == "builtin") {
    con <- con_func(filename)
  } else {
    con <- con_func(filename, mode="write")
  }
  on.exit(close(con))
  saveRDS(object, con)
  tictoc::toc()
}

### rreadRDS -----------------------------------------------------------------------------------
rreadRDS <- function(filename, con_func = list(pigz_pipe, snappy_pipe)[[1]], func_type = c("pipe", "builtin")[1], ...) {  # Read an object with parallel ungzip by pgzip.
  tictoc::tic()
  if(func_type == "builtin") {
    con <- con_func(filename)
  } else {
    con <- con_func(filename, mode="read", ...)
  }
  on.exit(close(con))
  obj=readRDS(con)
  tictoc::toc()
  return(obj)
}

# snappy_pipe -----------------------------------------------------------------------------------
snappy_pipe <- function(filename, mode="read") { # Alternative, fast compression. Low compression rate, lightning fast.
  if(mode == "read") {
    con <- pipe(paste0("cat ", filename, " | snzip -dc"), "rb")
  } else {
    con <- pipe(paste0("snzip > ", filename), "wb")
  }
  con
}

# pigz_pipe -----------------------------------------------------------------------------------
pigz_pipe <- function(filename, mode="read", cores=6) { # Alternative: normal gzip output (& compression rate), ~*cores faster in zipping.
  if(mode == "read") {
    con <- pipe(paste0("cat ", filename, " | pigz -dcp ", cores), "rb")
  } else {
    con <- pipe(paste0("pigz -p ", cores, " > ", filename), "wb")
  }
  con
}
