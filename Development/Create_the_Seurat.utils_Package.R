######################################################################################################
# Create_the_Seurat.utils_Package.R
######################################################################################################
# source("~/GitHub/Packages/Seurat.utils/Development/Create_the_Seurat.utils_Package.R")
# rm(list = ls(all.names = TRUE));
try(dev.off(), silent = TRUE)

# Functions ------------------------
require(PackageTools)
devtools::load_all("~/GitHub/Packages/PackageTools")

# Setup ------------------------
repository.dir <- "~/GitHub/Packages/Seurat.utils/"
config.path <- file.path(repository.dir, "Development/config.R")

"TAKE A LOOK AT"
file.edit(config.path)
source(config.path)

PackageTools::document_and_create_package(repository.dir, config_file = 'config.R')
'git add commit push to remote'

# Install your package ------------------------------------------------
"disable rprofile by"
rprofile()
devtools::install_local(repository.dir, upgrade = F)


# Test if you can install from github ------------------------------------------------
remote.path <- file.path(DESCRIPTION$'github.user', DESCRIPTION$'package.name')
pak::pkg_install(remote.path)
# unload(DESCRIPTION$'package.name')
# require(DESCRIPTION$'package.name')
# # remove.packages(DESCRIPTION$'package.name')


# CMD CHECK ------------------------------------------------
checkres <- devtools::check(repository.dir, cran = FALSE)



# Automated Codebase linting to tidyverse style ------------------------------------------------
styler::style_pkg(repository.dir)


# Extract package dependencies ------------------------------------------------
PackageTools::extract_package_dependencies(repository.dir)


# Visualize function dependencies within the package------------------------------------------------
{
  warning("works only on the installed version of the package!")
  pkgnet_result <- pkgnet::CreatePackageReport(DESCRIPTION$'package.name')
  fun_graph     <- pkgnet_result$FunctionReporter$pkg_graph$'igraph'
  PackageTools::convert_igraph_to_mermaid(graph = fun_graph, openMermaid = T, copy_to_clipboard = T)
}


# Try to find and add missing @importFrom statements------------------------------------------------
devtools::load_all("~/GitHub/Packages/PackageTools/")
if (F) {
  (excluded.packages <- unlist(strsplit(DESCRIPTION$'depends', split = ", ")))
  (ls.scripts.full.path <- list.files(file.path(repository.dir, "R"), full.names = T))
  for (scriptX in ls.scripts.full.path) {
    PackageTools::add_importFrom_statements(scriptX, exclude_packages = excluded.packages)
  }
}


# Generate the list of functions ------------------------------------------------
for (scriptX in ls.scripts.full.path) {
  PackageTools::list_of_funs_to_markdown(scriptX)
}

PackageTools::copy_github_badge("active") # Add badge to readme via clipboard


# Replaces T with TRUE and F with FALSE ------------------------------------------------
for (scriptX in ls.scripts.full.path) {
  PackageTools::replace_tf_with_true_false(scriptX)
}


