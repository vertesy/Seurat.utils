######################################################################################################
# Create_the_Seurat.utils_Package.v0.1.R
# 31 10 2021
######################################################################################################
# source("/Users/abel.vertesy/GitHub/Packages/Seurat.utils/Development/Create_the_Seurat.utils_Package.v0.1.R")
rm(list = ls(all.names = TRUE));
try(dev.off(), silent = TRUE)

# Functions ------------------------
# install_version("devtools", version = "2.0.2", repos = "http://cran.at.r-project.org") # install.packages("devtools")
require("devtools")
require("roxygen2")
require("stringr")

# devtools::install_github(repo = "vertesy/CodeAndRoll2")
require('CodeAndRoll2')
require('Stringendo')



# Setup ------------------------
PackageName = 	"Seurat.utils"
setwd("~/GitHub/Packages/")

RepositoryDir = kollapse("~/GitHub/Packages/", PackageName, "/")
fname = 	kollapse(PackageName, ".R")
Package_FnP = 	kollapse(RepositoryDir, "R/", fname)

BackupDir = "~/GitHub/Packages/Seurat.utils/Development/"
dir.create(BackupDir)

DESCRIPTION <- list("Title" = "Seurat.utils - utility functions for Seurat"
    , "Author" = person(given = "Abel", family = "Vertesy", email = "abel.vertesy@imba.oeaw.ac.at", role =  c("aut", "cre") )
    , "Authors@R" = 'person(given = "Abel", family = "Vertesy", email = "a.vertesy@imba.oeaw.ac.at", role =  c("aut", "cre") )'
    , "Description" = "Seurat.utils Is a collection of utility functions for Seurat v3.
    Functions allow the automation / multiplexing of plotting, 3D plotting, visualisation of statistics &
    QC, interaction with the Seurat object, etc. Some functionalities require functions from CodeAndRoll and MarkdownReports libraries."
    , "License" = "GPL-3 + file LICENSE"
    , "Version" = "1.1.0"
    , "Packaged" =  Sys.time()
    , "Repository" =  "CRAN"
    , "Depends" =  "Stringendo, CodeAndRoll2, ggExpress, ggplot2"
    , "Imports" = "base, cowplot, dplyr, ggcorrplot, ggpubr, ggrepel, graphics, grDevices, HGNChelper, htmlwidgets, MarkdownHelpers, MarkdownReports, Matrix, matrixStats, methods, princurve, ReadWriter, R.utils, readr, reshape2, scales, Seurat, sparseMatrixStats, stats, stringr, tibble, tictoc, utils, vroom"
    , "Suggests" = "SoupX"
    , "BugReports"= "https://github.com/vertesy/Seurat.utils/issues"
)


setwd(RepositoryDir)
if ( !dir.exists(RepositoryDir) ) { create(path = RepositoryDir, description = DESCRIPTION, rstudio = TRUE)
} else {
    getwd()
    try(file.remove(c("DESCRIPTION","NAMESPACE", "Seurat.utils.Rproj")))
    create_package(path = RepositoryDir, fields = DESCRIPTION, open = F)
}


# go and write fun's ------------------------------------------------------------------------
# file.edit(Package_FnP)

# Create Roxygen Skeletons ------------------------
# RoxygenReady(Package_FnP)

# replace output files ------------------------------------------------
BackupOldFile = 	kollapse(BackupDir, "Development", ".bac", print = FALSE)
AnnotatedFile = 	kollapse(BackupDir, "Development", ".annot.R", print = FALSE)
file.copy(from = Package_FnP, to = BackupOldFile, overwrite = TRUE)
# file.copy(from = AnnotatedFile, to = Package_FnP, overwrite = TRUE)

# Manual editing of descriptors ------------------------------------------------
# file.edit(Package_FnP)

# Compile a package ------------------------------------------------
setwd(RepositoryDir)
getwd()
document()


# Install your package ------------------------------------------------
# # setwd(RepositoryDir)
install(RepositoryDir, upgrade = F)
# require("Seurat.utils")
# # remove.packages("Seurat.utils")
# # Test your package ------------------------------------------------
# help("wplot")
# cat("\014")
# devtools::run_examples()


# Test if you can install from github ------------------------------------------------
# devtools::install_github(repo = "vertesy/Seurat.utils")

# require("Seurat.utils")

# Clean up if not needed anymore ------------------------------------------------
# View(installed.packages())
# remove.packages("Seurat.utils")

check(RepositoryDir, cran = TRUE)
# as.package(RepositoryDir)
#
#
# # source("https://install-github.me/r-lib/desc")
# # library(desc)
# # desc$set("Seurat.utils", "foo")
# # desc$get(Seurat.utils)
#
#
# system("cd ~/GitHub/Seurat.utils/; ls -a; open .Rbuildignore")

# Check package dependencies ------------------------------------------------
depFile = paste0(RepositoryDir, 'Development/Dependencies.R')

(f.deps <- NCmisc::list.functions.in.file(filename = Package_FnP))
# clipr::write_clip(f.deps)

sink(file = depFile); print(f.deps); sink()
p.deps <- gsub(x = names(f.deps), pattern = 'package:', replacement = '')
write(x = p.deps, file = depFile, append = T)



# NCmisc::list.functions.in.file(Package_FnP)
