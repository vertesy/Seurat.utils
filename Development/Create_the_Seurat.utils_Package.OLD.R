######################################################################################################
# Create_the_Seurat.utils_Package.R
######################################################################################################
# source("~/GitHub/Packages/Seurat.utils/Development/Create_the_Seurat.utils_Package.R")
# rm(list = ls(all.names = TRUE));
try(dev.off(), silent = TRUE)

# Functions ------------------------
# require("devtools")


# Setup ------------------------
package.name <- "Seurat.utils"
package.version <- "2.2.0"
setwd("~/GitHub/Packages/")

RepositoryDir <- paste0("~/GitHub/Packages/", package.name, "/")
fname <- paste0(package.name, ".R")
package.FnP <- paste0(RepositoryDir, "R/", fname)

BackupDir <- "~/GitHub/Packages/Seurat.utils/Development/"
dir.create(BackupDir)

DESCRIPTION <- list(
  "Title" = "Seurat.utils - utility functions for Seurat",
  "Author" = person(given = "Abel", family = "Vertesy", email = "av@imba.oeaw.ac.at", role = c("aut", "cre")),
  "Authors@R" = 'person(given = "Abel", family = "Vertesy", email = "av@imba.oeaw.ac.at", role =  c("aut", "cre") )',
    "Description" = "Seurat.utils is a collection of utility functions for Seurat v3.
      Functions allow the automation / multiplexing of plotting, 3D plotting, visualisation of statistics &
      QC, interaction with the Seurat object, etc. Some functionalities require functions from CodeAndRoll2 and MarkdownReports libraries.",
  "License" = "GPL-3 + file LICENSE",
  "Version" = package.version,
  "Packaged" = Sys.time()
  # , "Repository" =  "CRAN"
  , "Depends" = "ggplot2, Seurat, Stringendo, CodeAndRoll2, ggExpress, magrittr",
  "Imports" = "cowplot, Seurat, ReadWriter, dplyr, ggcorrplot, ggpubr, ggrepel, grDevices, HGNChelper,
                          htmlwidgets, MarkdownHelpers, MarkdownReports, Matrix, matrixStats, princurve, pheatmap,
                          R.utils, readr, reshape2, scales, SoupX, sparseMatrixStats, stringr, tibble, tictoc,
                          vroom, job, qs, foreach, tidyverse",
  "Suggests" = "SoupX, princurve, EnhancedVolcano",
  "BugReports" = "https://github.com/vertesy/Seurat.utils/issues"
)



setwd(RepositoryDir)
if (!dir.exists(RepositoryDir)) {
  create(path = RepositoryDir, description = DESCRIPTION, rstudio = TRUE)
} else {
  getwd()
  try(file.remove(c("DESCRIPTION", "NAMESPACE"))) # , "Seurat.utils.Rproj"
  usethis::create_package(path = RepositoryDir, fields = DESCRIPTION, open = F)
}



# go and write fun's ------------------------------------------------------------------------
# file.edit(package.FnP)

# Create Roxygen Skeletons ------------------------
# RoxygenReady(package.FnP)

# replace output files ------------------------------------------------
BackupOldFile <- paste0(BackupDir, "Development", ".bac")
AnnotatedFile <- paste0(BackupDir, "Development", ".annot.R")
file.copy(from = package.FnP, to = BackupOldFile, overwrite = TRUE)
# file.copy(from = AnnotatedFile, to = package.FnP, overwrite = TRUE)

# Manual editing of descriptors ------------------------------------------------
# file.edit(package.FnP)

# Compile a package ------------------------------------------------
setwd(RepositoryDir)
getwd()
devtools::document()
warnings()

{
  "update cff version"
  citpath <- paste0(RepositoryDir, "CITATION.cff")
  xfun::gsub_file(
    file = citpath, perl = T,
    "^version: v.+", paste0("version: v", package.version)
  )
}


# Install your package ------------------------------------------------
devtools::install(RepositoryDir, upgrade = F)

# Test if you can install from github ------------------------------------------------
pak::pkg_install("vertesy/Seurat.utils")
# unload("Seurat.utils")
# require("Seurat.utils")
# # remove.packages("Seurat.utils")


# Check CRAN ------------------------------------------------
checkres <- devtools::check(RepositoryDir, cran = TRUE)
head(checkres)
# as.package(RepositoryDir)
# # source("https://install-github.me/r-lib/desc")
# # library(desc)
# # desc$set("Seurat.utils", "foo")
# # desc$get(Seurat.utils)
# system("cd ~/GitHub/Seurat.utils/; ls -a; open .Rbuildignore")



# Check package dependencies ------------------------------------------------
{
  depFile <- paste0(RepositoryDir, "Development/Dependencies.R")

  (f.deps <- NCmisc::list.functions.in.file(filename = package.FnP))
  clipr::write_clip(f.deps$`character(0)`)

  sink(file = depFile)
  print(f.deps)
  sink()
  p.deps <- gsub(x = names(f.deps), pattern = "package:", replacement = "")
  write(x = p.deps, file = depFile, append = T)
}



# Package styling, and visualization ------------------------------------------------
{
  styler::style_pkg(RepositoryDir)
  # styler::style_file("~/GitHub/Packages/Seurat.utils/Development/Create_the_Seurat.utils_Package.R")

  {
    # Exploring the Structure and Dependencies of my R Package:
    "works on an installed package!"
    pkgnet_result <- pkgnet::CreatePackageReport(package.name)
    fun_graph <- pkgnet_result$FunctionReporter$pkg_graph$"igraph"

    # devtools::load_all('~/GitHub/Packages/PackageTools/R/DependencyTools.R')
    convert_igraph_to_mermaid(graph = fun_graph, openMermaid = T, copy_to_clipboard = T)
  }

  if (F) {
    # Add @importFrom statements
    (FNP <- package.FnP)
    PackageTools::add_importFrom_statements(FNP, exclude_packages = "")
    add_importFrom_statements(FNP, exclude_packages = "")
  }
}


if (F) {
  findFunctionOrigin <- function(function_names) {
    # Get a list of all installed packages
    installed_packages <- rownames(installed.packages())


    # Initialize an empty list to store the results
    function_origins <- list()

    # Iterate over each function name
    for (func_name in function_names) {
      print(func_name)
      # Initialize a vector to store the packages where the function is found
      found_in_packages <- c()

      # Check each package
      for (pkg in installed_packages) {
        cat(".", append = T)
        # Check if the package contains the function
        if (func_name %in% ls(getNamespace(pkg), all.names = TRUE)) {
          found_in_packages <- c(found_in_packages, pkg)
        }
      }

      # Add the vector of packages to the result list, named by the function
      function_origins[[func_name]] <- found_in_packages
    }

    # Return the list
    return(function_origins)
  }

  FunctionOrigins <- findFunctionOrigin(f.deps$`character(0)`)
}

# NCmisc::list.functions.in.file(package.FnP)
