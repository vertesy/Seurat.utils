######################################################################################################
# Create_the_Seurat.Utils_Package.v0.1.R
# 31 10 2021
######################################################################################################
# source("/Users/abel.vertesy/GitHub/Packages/Seurat.Utils/Development/Create_the_Seurat.Utils_Package.v0.1.R")
rm(list = ls(all.names = TRUE));
try(dev.off(), silent = TRUE)
# install.packages("devtools")
# Functions ------------------------
try (source('~/GitHub/Packages/CodeAndRoll/CodeAndRoll.R'),silent= FALSE)

# irequire("devtools")
# install_version("devtools", version = "2.0.2", repos = "http://cran.at.r-project.org")
irequire("devtools")
irequire("roxygen2")
irequire("stringr")

kollapse <-function(..., print = TRUE) {
if (print == TRUE) {
    print(paste0(c(...), collapse = ""))
  }
  paste0(c(...), collapse = "")
}

# Setup ------------------------
PackageName = 	"Seurat.Utils"
setwd("~/GitHub/Packages/")

RepositoryDir = kollapse("~/GitHub/Packages/", PackageName, "/")
fname = 	kollapse(PackageName, ".R")
Package_FnP = 	kollapse(RepositoryDir, "R/", fname)

BackupDir = "~/GitHub/Packages/Seurat.Utils/Development/"
dir.create(BackupDir)

DESCRIPTION <- list("Title" = "Seurat.Utils - utility functions for Seurat"
    , "Author" = person(given = "Abel", family = "Vertesy", email = "abel.vertesy@imba.oeaw.ac.at", role =  c("aut", "cre") )
    , "Authors@R" = 'person(given = "Abel", family = "Vertesy", email = "a.vertesy@imba.oeaw.ac.at", role =  c("aut", "cre") )'
    , "Description" = "Seurat.utils Is a collection of utility functions for Seurat v3.
    Functions allow the automation / multiplexing of plotting, 3D plotting, visualisation of statistics &
    QC, interaction with the Seurat object, etc. Some functionalities require functions from CodeAndRoll and MarkdownReports libraries."
    , "License" = "GPL-3 + file LICENSE"
    , "Version" = "0.9.1"
    # , "Version" = "4.0.0"
    , "Packaged" =  Sys.time()
    , "Repository" =  "CRAN"
    , "Imports" = "readr, gtools, openxlsx"
    # , "Suggests" = ""
    , "BugReports"= "https://github.com/vertesy/Seurat.Utils/issues"
)


setwd(RepositoryDir)
if ( !dir.exists(RepositoryDir) ) { create(path = RepositoryDir, description = DESCRIPTION, rstudio = TRUE)
} else {
    getwd()
    try(file.remove(c("DESCRIPTION","NAMESPACE", "Seurat.Utils.Rproj")))
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
install(RepositoryDir)
# require("Seurat.Utils")
# # remove.packages("Seurat.Utils")
# # Test your package ------------------------------------------------
# help("wplot")
# cat("\014")
# devtools::run_examples()


# Test if you can install from github ------------------------------------------------
# devtools::install_github(repo = "vertesy/Seurat.Utils")

# require("Seurat.Utils")

# Clean up if not needed anymore ------------------------------------------------
# View(installed.packages())
# remove.packages("Seurat.Utils")

check(RepositoryDir, cran = TRUE)
# as.package(RepositoryDir)
#
#
# # source("https://install-github.me/r-lib/desc")
# # library(desc)
# # desc$set("Seurat.Utils", "foo")
# # desc$get(Seurat.Utils)
#
#
# system("cd ~/GitHub/Seurat.Utils/; ls -a; open .Rbuildignore")
#
