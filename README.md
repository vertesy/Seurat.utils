# Seurat.utils ![status: active](https://raw.githubusercontent.com/vertesy/TheCorvinas/master/GitHub/Badges/active.svg) [![DOI](https://zenodo.org/badge/248721133.svg)](https://zenodo.org/badge/latestdoi/248721133) 

`Seurat.utils` Is a collection of utility functions for Seurat. Functions allow the automation / multiplexing of plotting, 3D plotting, visualisation of statistics & QC, interaction with the Seurat object, etc.  Some functionalities require functions from [CodeAndRoll2](https://github.com/vertesy/CodeAndRoll2), [ReadWriter](https://github.com/vertesy/ReadWriter), [Stringendo](https://github.com/vertesy/Stringendo), [ggExpressDev](https://github.com/vertesy/ggExpressDev), [MarkdownReports](https://github.com/vertesy/MarkdownReports), and the [Rocinante](https://github.com/vertesy/Rocinante) (See installation).

[TOC]



# Installation

Seurat.utils relies on:

- [Stringendo](https://github.com/vertesy/Stringendo)
- [ReadWriter](https://github.com/vertesy/ReadWriter)
- [CodeAndRoll2](https://github.com/vertesy/CodeAndRoll2)
- [MarkdownHelpers](https://github.com/vertesy/MarkdownHelpers)
- [Markdownreports](https://github.com/vertesy/Markdownreports)
- [ggExpress](https://github.com/vertesy/ggExpress)

... and provides functions for

- [Seurat.pipeline](https://github.com/vertesy/Seurat.pipeline)



You can install all of them directly from **GitHub** via **devtools** with one R command:

```R
# install.packages("devtools"); # If you don't have it.
BiocManager::install("sparseMatrixStats")
require("devtools")

# Install dependencies
devtools::install_github(repo = "vertesy/Stringendo", upgrade = F)
devtools::install_github(repo = "vertesy/CodeAndRoll2", upgrade = F)
devtools::install_github(repo = "vertesy/ReadWriter", upgrade = F)
devtools::install_github(repo = "vertesy/MarkdownHelpers", upgrade = F)
devtools::install_github(repo = "vertesy/Markdownreports", upgrade = F)
devtools::install_github(repo = "vertesy/ggExpress", upgrade = F)

# Recommended
devtools::install_github(repo = "vertesy/DatabaseLinke.R", upgrade = F)

# Install Seurat.utils
devtools::install_github(repo = "vertesy/Seurat.utils", upgrade = F)

# You may need to install
# install.packages("pheatmap")
# install.packages("checkmate")

```

...then simply load the package:

```R
require("Seurat.utils")
```

**NOTE: If you type 'all' when R asks to update dependencies, you may get into installation errors / infinite loops. If updating fails, type 'no' when prompted. **

Alternatively, you simply source it from the web. 
*This way function help will not work, and you will have no local copy of the code on your hard drive.*

```R
source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/R/Seurat.Utils.R")
```

<br>

## Troubleshooting

*If you encounter a **bug**, something doesn't work or unclear, please let me know by raising an issue on [Seurat.utils](https://github.com/vertesy/Seurat.utils/issues) – Please check if it has been asked.*



### Error on install

- Dependencies should be installed in the specific order.
- `Error: Failed to install 'unknown package' from GitHub:  HTTP error 403.  API rate limit exceeded`   -> Try connecting later or via another wifi
- Try updating all packages `update.packages()`



### Error during usage

#### #1 Check and reinstall dependencies

If you have an older installation, it is quite possible that some of the above packages are out of date. Please reinstall all of them, in order (see: `Installation`). 

Some functionalities are coming from dependencies. If the error points to a function outside or `Seurat.utils`, check in the packages below:

- [CodeAndRoll2](https://github.com/vertesy/CodeAndRoll2)
- [ReadWriter](https://github.com/vertesy/ReadWriter)
- [Stringendo](https://github.com/vertesy/Stringendo)
- [MarkdownReports](https://github.com/vertesy/MarkdownReports)

##### #2 Maybe you need [Rocinante](https://github.com/vertesy/Rocinante).

It needs to be sourced, cannot be installed:

```R
source("https://raw.githubusercontent.com/vertesy/Rocinante/master/R/Rocinante.R")
```

# List of Functions 

## List of Functions in Seurat.Utils.R (94) 

Updated: 2024/03/22 15:41

- #### 1 `parallel.computing.by.future()`

  parallel.computing.by.future. Run gc(), load multi-session computing and extend memory limits.

- #### 2 `IntersectGeneLsWithObject()`

  Intersect Genes with Seurat Object. Intersects a set of gene names with those found in a Seurat object.

- #### 3 `SelectHighlyExpressedGenesq99()`

  Intersect Genes with the List of Noticeably Expressed Genes. Intersects a vector of gene names with a Seurat object to find genes that are both  in the input list and have expression levels in the top quantiles as defined by the object's  q99 expression data. It aims to filter genes based on their expression levels being above a  specified threshold. Additionally, it offers an option to sort the genes by their expression  levels in decreasing order. 

- #### 4 `SmallestNonAboveX()`

  SmallestNonAboveX. replace small values with the next smallest value found, which is >X.

- #### 5 `Create.MiscSlot()`

  AreTheseCellNamesTheSame. Assert and compare two character vectors (e.g.: cell IDs) how much they overlap and  plot a Venn Diagram. The function aborts with an error if overlap is too small.

- #### 6 `addToMiscOrToolsSlot()`

  Add to Misc or Tools Slot. This function creates and adds a sub-slot to either the 'misc' or 'tools' slot of a  Seurat object. If the sub-slot already exists, it can either be overwritten or a warning will be issued. 

- #### 7 `showToolsSlots()`

  Display Slots in the @tools of an Seurat Object.   `showToolsSlots` prints the names of slots in the `@tools` of a given object.  It specifically targets list elements, skipping over data frames and other non-list objects. 

- #### 8 `showMiscSlots()`

  Display Slots in the @misc of an Seurat Object. See `showToolsSlots` for details. Prints the names of slots in the `@misc` of a given object.  It specifically targets list elements, skipping over data frames and other non-list objects. 

- #### 9 `calc.q99.Expression.and.set.all.genes()`

  calc.q99.Expression.and.set.all.genes. Calculate the gene expression of the e.g.: 99th quantile (expression in the top 1% cells).

- #### 10 `RenameClustering()`

  RenameClustering. Rename clustering in a Seurat object.

- #### 11 `getClusterNames()`

  Shorten Clustering Names. This function takes in a string representing a clustering name,  and shortens it according to specific rules. It replaces "snn_res." with "",  "best.matching.names" with "bmatch", "ordered" with "ord",  "ManualNames" with "mNames", and ".long" at the end of the string with ".L". 

- #### 12 `GetClusteringRuns()`

  GetClusteringRuns. Get Clustering Runs: metadata column names.

- #### 13 `GetNamedClusteringRuns()`

  GetNamedClusteringRuns. Get Clustering Runs: metadata column names

- #### 14 `GetOrderedClusteringRuns()`

  GetOrderedClusteringRuns. Get Clustering Runs: metadata column names.

- #### 15 `GetNumberOfClusters()`

  GetNumberOfClusters. Get Number Of Clusters #

- #### 16 `calc.cluster.averages()`

  calc.cluster.averages. Calculates the average of a metadata column (numeric) per cluster.

- #### 17 `plot.expression.rank.q90()`

  plot.expression.rank.q90. Plot gene expression based on the expression at the 90th quantile  (so you will not lose genes expressed in few cells).

- #### 18 `BackupReduction()`

  BackupReduction. Backup UMAP to `obj@misc$reductions.backup` from `obj@reductions$umap`. #

- #### 19 `SetupReductionsNtoKdimensions()`

  SetupReductionsNtoKdimensions. Function to calculate N-to-K dimensional umaps (default = 2:3); and back them up to  slots `obj@misc$reductions.backup` from @reductions$umap

- #### 20 `RecallReduction()`

  RecallReduction. Set active UMAP to `obj@reductions$umap` from `obj@misc$reductions.backup`.

- #### 21 `set.mm()`

  set.mm. Helps to find metadata columns. It creates a list with the names of of 'obj@meta.data'.

- #### 22 `ww.get.1st.Seur.element()`

  Get the First Seurat Object from a List of Seurat Objects. #' If provided with a list of Seurat objects, this function returns the first  Seurat object in the list. If the input is a single Seurat object, it returns  the object itself. It is assumed that all elements of the list are Seurat  objects if the input is a list. 

- #### 23 `recallAllGenes()`

  recallAllGenes. all.genes set by calc.q99.Expression.and.set.all.genes() #

- #### 24 `recall.meta.tags.n.datasets()`

  recall.meta.tags.n.datasets. Recall  meta.tags from obj@misc to "meta.tags" in the global environment.

- #### 25 `recall.parameters()`

  recall.parameters. Recall parameters from obj@misc to "p" in the global environment.

- #### 26 `recall.genes.ls()`

  recall.genes.ls. Recall genes.ls from obj@misc to "genes.ls" in the global environment.

- #### 27 `save.parameters()`

  Save Parameters to Seurat Object. Stores a list of parameters within the `@misc$p` slot of a Seurat object,  allowing for easy reference and tracking of analysis parameters used. 

- #### 28 `create_scCombinedMeta()`

  Create Single-Cell Metadata Object for a collection of Seurat Objects. This function creates a metadata object to correspond to a list of    single-cell experiments, for storing parent level information.    It initializes the object with the experiment and project name, and the    creation date. The created object is of class 'scMetadata_class'.

- #### 29 `copyMiscElements()`

  Copy Specified Elements from One Seurat Object's @misc to Another's. Copies specified elements from the `@misc` slot of one Seurat object to the `@misc` slot  of another. It warns if some specified elements are missing in the source object or if elements are  overwritten in the destination object, depending on the `overwrite` argument. 

- #### 30 `copyCompleteToolsSlots()`

  Copy Tools Slots from Multiple Seurat Objects. This function copies the `@tools` slots from a list of Seurat objects into a new slot  of a target Seurat object. This allows for the aggregation of tools information from multiple  experiments or datasets into a single, consolidated Seurat object. 

- #### 31 `subsetSeuObjByIdent()`

  Subset a Seurat Object by Identity. Subsets a Seurat object based on a specified identity column and values. It allows    for an optional inversion of the selection. 

- #### 32 `downsampleSeuObj()`

  downsampleSeuObj. Subset a compressed Seurat object and save it in the working directory.

- #### 33 `downsampleSeuObj.and.Save()`

  downsampleSeuObj.and.Save. Subset a compressed Seurat Obj and save it in wd. #

- #### 34 `downsampleSeuObjByIdentAndMaxcells()`

  Sample Cells From Identifiers in Seurat Object. This function samples a specified maximum number of cells from each identity class  in a Seurat object, in the meta.data. It ensures that the sampling does not exceed the total  number of cells available per identity. 

- #### 35 `removeResidualSmallClusters()`

  Remove Residual Small Clusters from Seurat Object. Removes clusters containing fewer cells than specified by `max.cells`  from a Seurat object. This function is particularly useful after subsetting a dataset,  where small, possibly unrepresentative clusters may remain. 

- #### 36 `dropLevelsSeurat()`

  dropLevelsSeurat. Drop unused levels from factor variables in a Seurat object.

- #### 37 `removeClustersAndDropLevels()`

  Remove Clusters and Drop Levels. This function removes residual small clusters from specified Seurat objects and  drops levels in factor-like metadata.

- #### 38 `removeCellsByUmap()`

  Remove Cells by Dimension Reduction. This function applies a cutoff in the specified dimension of a given  dimension reduction (UMAP, PCA, or t-SNE) to remove cells.

- #### 39 `downsampleListSeuObjsNCells()`

  Downsample a List of Seurat Objects to a Specific Number of Cells. Downsampling each Seurat object in a list to a specified number of cells. This function is  particularly useful for creating smaller, more manageable subsets of large single-cell datasets for  preliminary analyses or testing. 

- #### 40 `downsampleListSeuObjsPercent()`

  Downsample a List of Seurat Objects to a Fraction. Downsampling a list of Seurat objects to a specified fraction of their original size.  This is useful for reducing dataset size for quicker processing or testing workflows. 

- #### 41 `Add.DE.combined.score()`

  Add.DE.combined.score. Add a combined score to differential expression (DE) results. The score is  calculated as log-fold change (LFC) times negative logarithm of scaled  p-value (LFC * -log10( p_cutoff / pval_scaling )).

- #### 42 `StoreTop25Markers()`

  Save Top 25 Markers per Cluster. Stores the top 25 markers for each cluster identified in a Seurat object, based on  the `avg_log2FC` from the output table of `FindAllMarkers()`. The result is saved under `@misc$df.markers$res...`,  rounding insignificant digits to three decimal places. 

- #### 43 `StoreAllMarkers()`

  Store All Differential Expression Markers. Saves the complete output table from `FindAllMarkers()` to a Seurat object, facilitating  easy access to differential expression analysis results. This function rounds numerical values to a  specified number of digits to maintain readability and manage file sizes. 

- #### 44 `GetTopMarkersDF()`

  Get Top Differential Expression Genes Data Frame. Retrieves a data frame of the top N differentially expressed genes from  differential gene expression analysis results, offering an option to exclude certain genes  based on patterns. 

- #### 45 `GetTopMarkers()`

  Get Top Differential Expression Markers from DGEA Results. Retrieves the top N differentially expressed genes from the results of a differential  gene expression analysis, such as that provided by `FindAllMarkers()`. 

- #### 46 `AutoLabelTop.logFC()`

  AutoLabelTop.logFC. Create a new "named identity" column in the metadata of a Seurat object,  with `Ident` set to a clustering output matching the `res` parameter of the function.  It requires the output table of `FindAllMarkers()`.  If you used `StoreAllMarkers()`  is stored under `@misc$df.markers$res...`, which location is assumed by default.

- #### 47 `AutoLabel.KnownMarkers()`

  AutoLabel.KnownMarkers. Creates a new "named identity" column in the metadata of a Seurat object,   setting 'Ident' to a clustering output matching the 'res' parameter.   This function requires the output table of `FindAllMarkers()`.  If you used `StoreAllMarkers()`, the output is stored under `@misc$df.markers$res...`,  which is the default location.

- #### 48 `# sparse.cor4()`

  Calculate Sparse Correlation Matrix. Computes a sparse correlation matrix from a given sparse matrix input. This function is  useful for efficiently handling large datasets where most values are zero, facilitating the calculation  of both covariance and correlation matrices without converting to a dense format. 

- #### 49 `Calc.Cor.Seurat()`

  Calc.Cor.Seurat. Calculate gene correlation on a Seurat object.

- #### 50 `plot.Gene.Cor.Heatmap()`

  Plot Gene Correlation Heatmap. Generates a heatmap visualization of gene correlations based on expression data.  Useful for identifying groups of genes that exhibit similar expression patterns across different conditions  or cell types in a Seurat object. 

- #### 51 `prefix_cells_seurat()`

  Add Prefixes to Cell Names in Seurat Objects. Adds prefixes derived from a vector of identifiers to cell names in a list of Seurat objects.  This is useful for ensuring unique cell names across multiple samples or conditions when combining or comparing datasets. 

- #### 52 `find_prefix_in_cell_IDs()`

  Check Prefix in Seurat Object Cell IDs. This function checks if a prefix has been added to the standard  cell-IDs (16 characters of A,TRUE,C,G) in a Seurat object. If so, it prints the number of unique prefixes found,  issues a warning if more than one unique prefix is found, and returns the identified prefix(es). 

- #### 53 `seu.Make.Cl.Label.per.cell()`

  Create Cluster Labels for Each Cell. Generates labels for each cell by combining gene names and cluster IDs. This function  takes a named vector, typically representing top genes for clusters (values) and their corresponding  cluster IDs (names), along with a vector of cell IDs. It then creates a new vector where each cell  is labeled with its top gene and cluster ID in the format "GeneName.ClusterID". 

- #### 54 `GetMostVarGenes()`

  Retrieve the Top Variable Genes from a Seurat Object. Retrieves the names of the most variable genes from a Seurat object,  typically used to focus subsequent analyses on genes with the greatest variation across cells. 

- #### 55 `gene.name.check()`

  Check Gene Names in Seurat Object. Examines gene names in a Seurat object for specific naming conventions,  such as the presence of hyphens (-) or dots (.) often found in mitochondrial gene names.  This function is useful for ensuring gene names conform to expected patterns,  especially when preparing data for compatibility with other tools or databases. 

- #### 56 `check.genes()`

  Check if Gene Names exist in Seurat Object or HGNC Database. Verifies the presence of specified gene names within a Seurat object or  queries them against the HGNC database. This function is useful for ensuring gene names are  correctly formatted and exist within the dataset or are recognized gene symbols. 

- #### 57 `fixZeroIndexing.seurat()`

  Fix Zero Indexing in Seurat Clustering. Adjusts Seurat object metadata to fix zero-based cluster indexing, converting it to one-based indexing.  This function modifies a specified metadata column in the Seurat object to replace zero-indexed cluster names with one-based indexing. 

- #### 58 `CalculateFractionInTrome()`

  Calculate Fraction of Genes in Transcriptome. Calculates the fraction of specified genes within the entire transcriptome of  each cell in a Seurat object.  This function is useful for assessing the relative abundance of a set of genes across cells,  such as identifying cells with high expression of marker genes. 

- #### 59 `AddNewAnnotation()`

  AddNewAnnotation. This function creates a new metadata column based on an existing metadata column  and a list of mappings (name <- IDs).

- #### 60 `whitelist.subset.ls.Seurat()`

  whitelist.subset.ls.Seurat. Subsets cells in a list of Seurat objects based on an externally provided list of cell IDs.

- #### 61 `FindCorrelatedGenes()`

  FindCorrelatedGenes. Find correlated genes in a Seurat object

- #### 62 `UpdateGenesSeurat()`

  Update Gene Symbols in a Seurat Object. This function updates gene symbols in a Seurat object based on current gene  nomenclature guidelines, using HGNChelper(). It checks and updates gene symbols to their  latest approved versions,ensuring that gene annotations are current and consistent.  The function optionally enforces unique gene symbols and provides statistics on the update process. 

- #### 63 `RenameGenesSeurat()`

  Rename Gene Symbols in a Seurat Object. This function replaces gene names across various slots within a specified assay  of a Seurat object. It is designed to be run prior to any data integration or downstream analysis  processes. The function targets the `@counts`, `@data`, and `@meta.features` slots within  the specified assay, ensuring consistency in gene nomenclature across the object. 

- #### 64 `.check_and_rename()`

  Check and Rename Gene Names in Seurat Assay Object. This function renames rows (genes) in a specified slot of a Seurat assay object.  It supports slots storing data as either a dense or a sparse matrix (dgCMatrix) or data.frame. 

- #### 65 `RemoveGenesSeurat()`

  Remove Specific Genes from a Seurat Object. Removes specified genes from the metadata, counts, data, and scale.data slots of a Seurat object.  This operation is typically performed prior to data integration to ensure that gene sets are consistent  across multiple datasets. The function modifies the Seurat object in place. 

- #### 66 `HGNC.EnforceUnique()`

  Enforce Unique HGNC Gene Symbols. Ensures that gene symbols are unique after being updated with HGNC symbols. This function  applies a suffix to duplicate gene symbols to enforce uniqueness. While using `make.unique` might not  be the ideal solution due to potential mismatches, it significantly reduces the number of mismatching  genes in certain scenarios, making it a practical approach for data integration tasks. 

- #### 67 `GetUpdateStats()`

  Gene Symbol Update Statistics. Generates statistics on the gene symbol updates performed by `UpdateGenesSeurat()`.  This function analyzes the data frame of gene symbols before and after the update process,  providing insights into the proportion and total number of genes that were updated. 

- #### 68 `PlotUpdateStats()`

  PlotUpdateStats. Creates a scatter plot of update statistics.

- #### 69 `Convert10Xfolders()`

  Convert10Xfolders. This function takes a parent directory with a number of subfolders, each  containing the standard output of 10X Cell Ranger. It (1) loads the filtered data matrices,  (2) converts them to Seurat objects, and (3) saves them as .RDS files.

- #### 70 `ConvertDropSeqfolders()`

  ConvertDropSeqfolders. This function takes a parent directory with a number of subfolders, each  containing the standard output of 10X Cell Ranger. It (1) loads the filtered data matrices,  (2) converts them to Seurat objects, and (3) saves them as .RDS files.

- #### 71 `LoadAllSeurats()`

  LoadAllSeurats. This function loads all Seurat objects found in a directory. It also works with  symbolic links (but not with aliases).

- #### 72 `read10x()`

  Load 10X Genomics Data as Seurat Object. Reads 10X Genomics dataset files (gzipped) including matrix, features, and barcodes,  to a single expression matrix. This function handles the unzipping of these files, reads the data,  and re-compresses the files back to their original gzipped format. 

- #### 73 `.saveRDS.compress.in.BG()`

  .saveRDS.compress.in.BG. Save and RDS object and compress resulting file in the background using system(gzip). OS X or unix.

- #### 74 `isave.RDS()`

  isave.RDS. Save an RDS object, using a faster and efficient compression method that runs in the background.

- #### 75 `xsave()`

  Save an R Object Using 'qs' Package for Fast Compressed Saving. This function saves an R object to a file in a quick and efficient format using the 'qs' package.  It constructs the file name based on various inputs and stores additional metadata if the object is a Seurat object.  The saving path can be adjusted by the presence of 'OutDir' in the global environment or defaults to the working directory. 

- #### 76 `xread()`

  Read an R Object Using 'qs' Package for Fast Decompression. This function reads an R object from a file saved in a format specific to the 'qs' package,  which is designed for quick and efficient compression and decompression of R objects.  It also times the read operation, providing feedback on the duration of the operation. 

- #### 77 `.getCPUsCBE()`

  Get the number of CPUs to use for CBE processing. This function checks for the presence of a global `CBE.params` list and,  if found and contains a `cpus` entry, returns the number of CPUs specified by `cpus` minus one.  Otherwise, it returns a default number of CPUs. 

- #### 78 `isave.image()`

  isave.image. Save an image of the current workspace using a faster and efficient compression  method that runs in the background.

- #### 79 `qsave.image()`

  Save workspace - qsave.image. Faster saving of workspace, and compression outside R, when it can run in the background.  Seemingly quite CPU hungry and not very efficient compression. #

- #### 80 `# findBamFilesInSubdirs()`

  Find 'Outs' Subdirectories in Specified Subdirectories. This function searches through specified subdirectories within a root directory  to find all subdirectories named 'outs' and returns a character vector with their full paths. 

- #### 81 `clip10Xcellname()`

  Clip Suffixes from 10X Cell Names. Removes suffixes from cell names that are added by 10X technology and Seurat during data processing. 

- #### 82 `make10Xcellname()`

  Add Suffix to Cell Names (e.g. lane suffix: _1). Appends a specified suffix to cell names to mimic lane suffixes used in 10X datasets. 

- #### 83 `plotTheSoup()`

  plotTheSoup. Plot stats about the ambient RNA content in a 10X experiment.

- #### 84 `jJaccardIndexVec()`

  jJaccardIndexVec. Calculate jaccard similarity for 2 vecotrs. Helper to jPairwiseJaccardIndexList.

- #### 85 `jPairwiseJaccardIndexList()`

  jPairwiseJaccardIndexList. Create a pairwise jaccard similarity matrix across all combinations of columns in  binary.presence.matrix. Modified from:  https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/

- #### 86 `jPresenceMatrix()`

  jPresenceMatrix. Make a binary presence matrix from a list. Source:  https://stackoverflow.com/questions/56155707/r-how-to-create-a-binary-relation-matrix-from-a-list-of-strings #

- #### 87 `jJaccardIndexBinary()`

  jJaccardIndexBinary. Calculate Jaccard Index. Modified from:  https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/ #

- #### 88 `jPairwiseJaccardIndex()`

  jPairwiseJaccardIndex. Create a pairwise jaccard similarity matrix across all combinations of columns in  binary.presence.matrix. Modified from:  https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/

- #### 89 `compareVarFeaturesAndRanks()`

  Compare variable features and their ranks in two Seurat objects.. This function compares variable features (genes) between two Seurat objects,    reporting the number of genes in each, the percentage of common genes, the percentage    of unique genes in each object, and the similarity in the ranking of overlapping genes    using Spearman's rank correlation coefficient. Optionally, it can also generate a scatterplot    of the ranks of common genes using ggpubr's ggscatter. The function returns the common genes    and the Spearman's rank correlation coefficient. 

- #### 90 `processSeuratObject()`

  Process Seurat Objects in Parallel. Applies a series of Seurat processing steps to each Seurat object in a list.               The operations include scaling data, running PCA, UMAP, finding neighbors, and finding clusters.               This is done in parallel using multiple cores. 

- #### 91 `regress_out_and_recalculate_seurat()`

  Regress Out and Recalculate Seurat. The function performs a series of calculations and manipulations on a Seurat object,  including identifying variable features, scaling data, running PCA, setting up reductions, finding neighbors,  and finding clusters. It optionally performs t-SNE and saves the object. 

- #### 92 `.parseRegressionVariablesForScaleData()`

  Check List Elements. Tests if list elements are defined and reports the value or issues a warning. 

- #### 93 `.parseKeyParams()`

  Parse key parameters from an object and format as a string. This function extracts the number of scaled features, the number of principal components,  and formats additional information including regression variables.

- #### 94 `removeScaleData()`

  Parse basic obj stats. Parse cell and feature count to a string.



## List of Functions in Seurat.Utils.Visualization.R (57) 

Updated: 2024/03/22 15:41

- #### 1 `PlotFilters()`

  PlotFilters. Plot filtering threshold and distributions, using four panels to highlight the relation between Gene- and UMI-count, ribosomal- and mitochondrial-content.

- #### 2 `scCalcPCAVarExplained()`

  Calculate the percent of variation explained by individual PC's. This function calculates the percentage of variation each principal component (PC)  accounts for in a Seurat object. It's specifically tailored for Seurat objects and provides a  convenient way to understand the variance distribution across PCs. For similar calculations on  standard PCA objects, refer to github.com/vertesy/Rocinante `PCA.percent.var.explained()`. 

- #### 3 `scPlotPCAvarExplained()`

  Plot the percent of variation explained by individual PC's. This function plots the percentage of variation explained by each principal  component (PC) in a Seurat object. It allows for a visual assessment of how much variance is  captured by each PC, which is crucial for dimensionality reduction analyses. Users can choose  between two plotting methods: one using `MarkdownReports` and the other using `ggExpress`. 

- #### 4 `Percent.in.Trome()`

  Gene Expression as Fraction of Total UMI Counts. This function computes and visualizes gene expression levels as a fraction of total  UMI (Unique Molecular Identifier) counts across all genes in a Seurat object. It aims to highlight  the relative contribution of the most highly expressed genes to the overall transcriptome. 

- #### 5 `geneExpressionLevelPlots()`

  Histogram of Gene Expression Levels. This function generates a histogram to visualize the expression level distribution  of a specified gene across all cells in a Seurat object. It highlights the position of the gene  of interest within the overall distribution. 

- #### 6 `PrctCellExpringGene()`

  Proportion of Cells Expressing Given Genes. Calculates the proportion of cells expressing one or more specified genes. 

- #### 7 `ww.calc_helper()`

  Helper to calculate Cell Expression Proportion for Gene. Computes the proportion of cells expressing a specific gene within a Seurat object. 

- #### 8 `get.clustercomposition()`

  Cluster Composition Analysis. Analyzes and visualizes the composition of clusters in a Seurat object, indicating  the contribution of different datasets to each cluster. 

- #### 9 `scBarplot.CellFractions()`

  Generate Barplot of Cell Fractions. This function generates a bar plot of cell fractions per cluster from a Seurat object.  It offers the option to downsample data, equalizing the number of cells in each group  to the number in the smallest group. The plot's bars are grouped by one variable and filled by another.  The function supports custom color palettes, drawing numerical values on bars, and saving the plot. 

- #### 10 `scBarplot.CellsPerCluster()`

  Barplot of Fraction of Cells per Cluster. Visualizes the fraction of cells within each cluster through a barplot. 

- #### 11 `scBarplot.CellsPerObject()`

  Barplot of Cells Per Seurat Object. Visualizes the number of cells in each Seurat object within a list, showing the  distribution of cell counts across different datasets or experimental conditions. 

- #### 12 `plotClustSizeDistr()`

  Cluster Size Distribution Plot (Barplot or Histogram). Generates a bar plot or histogram to visualize the size distribution of clusters  within a Seurat object, based on the specified clustering identity. 

- #### 13 `scBarplot.FractionAboveThr()`

  Barplot the Fraction of Cells Above Threshold per Cluster. Generates a bar plot depicting the percentage of cells within each cluster that  exceed a specified threshold, based on a selected metadata column. 

- #### 14 `scBarplot.FractionBelowThr()`

  Fraction of Cells Below Threshold per Cluster. Generates a bar plot to visualize the percentage of cells within each cluster that  fall below a specified threshold, according to a metadata column value. 

- #### 15 `scBarplotStackedMetaCateg_List()`

  Stacked Barplot of Metadata Categories for List of Seurat Objects. Creates and saves a stacked barplot for a specified metadata category  from a list of Seurat objects. 

- #### 16 `gg_color_hue()`

  Reproduce the ggplot2 default color palette. Generates a vector of colors that emulates the default color palette used by ggplot2.  This function is useful for creating color sets for custom plotting functions or for applications  outside of ggplot2 where a similar aesthetic is desired. 

- #### 17 `# getDiscretePalette()`

  Safely generate a discrete color palette (NA).. Safe wrapper around Seurat's DiscretePalette(), which returns NA's if too many  categories are requested

- #### 18 `getDiscretePaletteObj()`

  Generate a Discrete Color Palette for Seurat Clusters. Generates a discrete color palette for visualizing clusters in a Seurat object,  using a specified identity column to determine the number of unique clusters. 

- #### 19 `DiscretePaletteSafe()`

  Safely generate a Discrete color palette.. Generates a discrete color palette, ensuring no NA values are included, suitable  for visualizations where a specific number of distinct, reproducible colors is needed. 

- #### 20 `getClusterColors()`

  Regenerate Cluster Colors from a Seurat Object. Regenerate and optionally displays the color scheme associated with the clusters  in a Seurat object as defined by a specified identity column. 

- #### 21 `SeuratColorVector()`

  Regenerate Color Scheme for Clusters in Seurat Object as a vector. Extracts and optionally displays the color scheme assigned to cluster identities  within a Seurat object, facilitating consistent color usage across visualizations. You can  check results in a barplot with `MarkdownHelpers::color_check()`. 

- #### 22 `plotAndSaveHeatmaps()`

  Plot and Save Heatmaps from Metadata Calculation Results. Generates and saves heatmap visualizations for each metric in the results obtained  from metadata calculations, such as mean or median values of specified features  across different categories. 

- #### 23 `qFeatureScatter()`

  Scatter Plot of Two Features in Seurat Object. Generates a scatter plot comparing two features (genes or metrics) from a Seurat  object and optionally saves it. The function wraps around Seurat's `FeatureScatter` for  enhanced usability, including optional logarithmic transformations and saving capabilities. 

- #### 24 `qSeuViolin()`

  Create a Violin Plot for a Seurat Object Feature and save the file.. Generates a violin plot for a specified feature in a Seurat object,  allowing for the data to be split by a specified grouping variable.  The function supports customization options such as logarithmic scaling, custom titles, and more. 

- #### 25 `plotGeneExpHist()`

  Histogram of Gene Expression in Seurat Object. Creates and optionally saves a histogram showing expression levels of specified genes  within a Seurat object. Provides options for aggregate gene expression, expression threshold filtering,  and quantile clipping for count data. 

- #### 26 `qUMAP()`

  Quick UMAP Visualization of Gene Expression and automatically save the plot. Generates a UMAP visualization for a specific feature from a Seurat object, and  automatically saves it. Offers options for custom titles, subtitles, saving, and more. Assumes  default options for custom titles, subtitles, saving, and more. 

- #### 27 `clUMAP()`

  Quick Visualization of Clustering Results with UMAP and automatically save the plot. Generates a UMAP visualization based on clustering results from a Seurat object,  and automatically saves it. Offers options for custom titles, subtitles, saving, and more. Assumes  default options for custom titles, subtitles, saving, and more. 

- #### 28 `umapHiLightSel()`

  Highlight Selected Clusters on UMAP. Generates a UMAP plot from a Seurat object with specified clusters highlighted.  It saves the resulting UMAP plot directly to the current working directory. 

- #### 29 `DimPlot.ClusterNames()`

  DimPlot.ClusterNames. Plot UMAP with Cluster names.

- #### 30 `multiFeaturePlot.A4()`

  multiFeaturePlot.A4. Save multiple FeaturePlots, as jpeg, on A4 for each gene, which are stored as a list of gene names.

- #### 31 `multiSingleClusterHighlightPlots.A4()`

  Generate Cluster Highlight UMAPs compiled into A4 pages. This function generates and saves cluster highlight plots for both single and multiple  clusters using UMAP or other dimensionality reduction techniques. It supports saving plots in various  formats and allows customization of plot appearance and layout. 

- #### 32 `multiFeatureHeatmap.A4()`

  multiFeatureHeatmap.A4. Save multiple FeatureHeatmaps from a list of genes on A4 jpeg.

- #### 33 `multi_clUMAP.A4()`

  Plot multiple categorical variables in combined UMAPs. Generates and saves multiple UMAP plots for clustering results, adjusting the  layout and plot dimensions. Supports the generation of plots in different  formats and customization of the visual appearance. 

- #### 34 `qClusteringUMAPS()`

  Quick Clustering UMAPs on A4 Page. Generates and arranges UMAP plots for up to four specified clustering resolutions  from a Seurat object onto an A4 page, facilitating comparative visualization. 

- #### 35 `plotQUMAPsInAFolder()`

  Plot qUMAPs for Genes in a Folder. This function plots qUMAPs for a specified set of genes, storing the results in a  specified folder. If no folder name is provided, it defaults to using the gene set name. 

- #### 36 `PlotTopGenesPerCluster()`

  Plot Top N Differentially Expressed Genes Per Cluster. Visualizes the top N differentially expressed (DE) genes for each cluster within a  specified clustering resolution of a Seurat object, facilitating the exploration of gene  expression patterns across clusters. 

- #### 37 `qQC.plots.BrainOrg()`

  Quickly Plot Key QC Markers in Brain Organoids. Generates and arranges UMAP plots for specified QC features  from a Seurat object on an A4 page, facilitating a quick quality control (QC) overview. 

- #### 38 `qMarkerCheck.BrainOrg()`

  Quickly Plot Key Markers in Brain Organoids. Generates plots for a predefined or custom set of gene markers within brain organoids,  aiding in the quick assessment of their expression across different cells or clusters. 

- #### 39 `PlotTopGenes()`

  Plot Top Genes. This function plots the highest expressed genes on UMAPs, saving the plots in a  subfolder. It requires the prior execution of `calc.q99.Expression.and.set.all.genes`. 

- #### 40 `FlipReductionCoordinates()`

  Flip Reduction Coordinates. Flips dimensionality reduction coordinates (such as UMAP or tSNE) vertically or  horizontally to change the visualization perspective. 

- #### 41 `AutoNumber.by.UMAP()`

  Relabel Cluster Numbers Along a UMAP (or tSNE) Axis. Automatically renumbers clusters based on their position along a specified dimension  in a UMAP (or tSNE or PCA) plot, potentially enhancing interpretability by ordering clusters. 

- #### 42 `.adjustLayout()`

  Adjust Layout Parameters for multi* plotting fucntions. Adjusts layout dimensions and properties based on the specified layout type.               Updates the provided environment with new dimensions and layout configuration. 

- #### 43 `save2plots.A4()`

  Save Two Plots on One A4 Page. Arranges and saves two UMAP plots (or any plots) side-by-side or one above  the other on a single A4 page. 

- #### 44 `save4plots.A4()`

  Save Four Plots on One A4 Page. Arranges and saves four plots (e.g. UMAPs) onto a single A4 page, allowing for a  compact comparison of different visualizations or clustering results. 

- #### 45 `qqSaveGridA4()`

  qqSaveGridA4. Saves a grid of 2 or 4 ggplot objects onto an A4 page.

- #### 46 `ww.check.if.3D.reduction.exist()`

  ww.check.if.3D.reduction.exist. ww.check.if.3D.reduction.exist in backup slot #

- #### 47 `()`

  ww.check.quantile.cutoff.and.clip.outliers. Function to check a specified quantile cutoff and clip outliers from a given expression vector.

- #### 48 `plot3D.umap.gene()`

  plot3D.umap.gene. Plot a 3D umap with gene expression. Uses plotly. Based on github.com/Dragonmasterx87.

- #### 49 `plot3D.umap()`

  plot3D.umap. Plot a 3D umap based on one of the metadata columns. Uses plotly. Based on github.com/Dragonmasterx87.

- #### 50 `SavePlotlyAsHtml()`

  SavePlotlyAsHtml. Save a Plotly 3D scatterplot as an HTML file.

- #### 51 `BackupReduction()`

  Backup Dimensionality Reduction Data. Stores a backup of specified dimensionality reduction data (e.g., UMAP, tSNE, PCA)  within the Seurat object, from `obj@reductions$umap` to the `@misc$reductions.backup` slot. 

- #### 52 `SetupReductionsNtoKdimensions()`

  Compute and Backup Dimensionality Reductions. Executes specified dimensionality reduction (UMAP, tSNE, or PCA) over a range of dimensions  and backs up the results within a Seurat object. This function allows for exploration of data structure  at varying levels of granularity and ensures that reduction results are preserved for future reference. 

- #### 53 `RecallReduction()`

  Recall Dimensionality Reduction from backup slot. Restores dimensionality reduction data (e.g., UMAP, tSNE, PCA) from a backup  stored within `obj@misc$reductions.backup` to the active `obj@reductions` slot. 

- #### 54 `Annotate4Plotly3D()`

  Annotate4Plotly3D. Create annotation labels for 3D plots. Source https://plot.ly/r/text-and-annotations/#3d-annotations.

- #### 55 `Plot3D.ListOfGenes()`

  Plot3D.ListOfGenes. Plot and save list of 3D UMAP or tSNE plots using plotly.

- #### 56 `Plot3D.ListOfCategories()`

  Plot3D.ListOfCategories. This function plots and saves a list of 3D UMAP or tSNE plots using plotly.

- #### 57 `# panelCorPearson()`

  Display Correlation Values in Pairs Plot. This function displays the correlation coefficient and significance level within  a scatterplot generated by the `pairs()` function. The default correlation method is Pearson,  but Kendall or Spearman methods can also be selected. 



## List of Functions in Seurat.Utils.Metadata.R (27) 

Updated: 2024/03/22 15:41

- #### 1 `getMetaColnames()`

  Get Metadata Column Names Matching Pattern. Retrieves column names from an object's metadata that match a specified pattern. 

- #### 2 `metaColnameExists()`

  Check if a Column Exists in the Metadata of an S4 Object. This function checks whether a given column exists in the meta.data of a Seurat object.

- #### 3 `getMetadataColumn()`

  getMetadataColumn. Retrieves a specified metadata column from a Seurat object and returns it as a named vector.

- #### 4 `get_levels_seu()`

  Get Unique Levels of a Seurat Object Ident Slot. This function extracts the unique levels present in the 'ident' slot of a Seurat object.  The function throws an error if the number of levels exceeds 'max_levels'.  The function optionally prints the R code to recreate the 'Levels' vector using 'dput'. 

- #### 5 `calculateAverageMetaData()`

  Calculate Average Metadata for Seurat Object. Computes specified metrics (e.g., median, mean) for given metadata features across each category  defined by an identity column in a Seurat object's metadata. This function allows for flexible  metric calculation on specified features, providing insights into the data distribution. 

- #### 6 `getMedianMetric.lsObj()`

  getMedianMetric.lsObj. Get the median values of different columns in meta.data, can iterate over a list of Seurat objects.

- #### 7 `getCellIDs.from.meta()`

  getCellIDs.from.meta. Retrieves cell IDs from a specified metadata column of a Seurat object, where the cell ID matches a provided list of values. The matching operation uses the `%in%` operator.

- #### 8 `addMetaDataSafe()`

  Add Metadata to a Seurat object, safely with Checks. Wrapper function for `AddMetaData` that includes additional checks and assertions. 

- #### 9 `seu.add.meta.from.vector()`

  seu.add.meta.from.vector. Adds a new metadata column to a Seurat object.

- #### 10 `create.metadata.vector()`

  Create a Metadata Vector. This function creates a metadata vector from an input vector and a Seurat object.  The resulting vector contains values from 'vec' for the intersecting cell names between 'vec' and 'obj'.  It also checks if the intersection between the cell names in 'vec' and 'obj' is more than a minimum intersection size.

- #### 11 `addMetaFraction()`

  addMetaFraction. Add a new metadata column to a Seurat object, representing the fraction of a gene set in the transcriptome (expressed as a percentage).

- #### 12 `add.meta.tags()`

  add.meta.tags. Add metadata tags to a Seurat object dataset.

- #### 13 `seu.add.meta.from.table()`

  seu.add.meta.from.table. Add multiple new metadata columns to a Seurat object from a table. #

- #### 14 `seu.map.and.add.new.ident.to.meta()`

  seu.map.and.add.new.ident.to.meta. Adds a new metadata column to a Seurat object based on an identity mapping table.

- #### 15 `fix.orig.ident()`

  fix.orig.ident. Remove the string "filtered_feature_bc_matrix." from "orig.ident". Helper function.

- #### 16 `seu.RemoveMetadata()`

  seu.RemoveMetadata. Remove specified metadata columns from a Seurat object.

- #### 17 `saveLsSeuratMetadata()`

  Save Metadata from a List of Seurat Objects. This function takes a list of Seurat objects, extracts their metadata, and saves it to a file with a specified suffix. 

- #### 18 `# transferMetadataV1()`

  Transfer Multiple Metadata Columns Between Two Seurat Objects. Transfers specified metadata columns from one Seurat object to another,  with options for verbose output and overwriting existing columns. Checks for cell overlap and  reports percentages of matching and unique cells. 

- #### 19 `sampleNpc()`

  Sample N % of a dataframe (obj@metadata), and return rownames (cell IDs).. This function samples a specified percentage of a dataframe (specifically a subset    of the metadata of a Seurat object) and returns the corresponding cell IDs.

- #### 20 `writeCombinedMetadataToTsvFromLsObj()`

  Combine Metadata from a list of Seurat objects and Write to TSV.   Formerly `writeMetadataToTsv`. `writeCombinedMetadataToTsvFromLsObj` takes a list of ls.Obj, extracts their `@meta.data` slots,  removes specified columns, checks for column consistency, creates a barplot showing the number  of rows per object, and finally merges these into one large data frame. 

- #### 21 `plotMetadataCorHeatmap()`

  Plot Metadata Correlation Heatmap. This function plots a heatmap of metadata correlation values. It accepts a Seurat object  and a set of metadata columns to correlate. The correlations are calculated using either Pearson  or Spearman methods, and the resulting heatmap can include the principal component (PCA) values  and be saved with a specific suffix. 

- #### 22 `heatmap_calc_clust_median()`

  Calculate and plot heatmap of cluster medians. This function calculates the median of specified variables in a dataframe,  grouped by a column ('ident'). The function also provides an option to scale the medians,  subset the ident levels, and either return a matrix of median values or plot a heatmap. 

- #### 23 `plotMetadataCategPie()`

  plotMetadataMedianFractionBarplot. Generates a barplot of metadata median values.

- #### 24 `renameAzimuthColumns()`

  Rename Azimuth Columns in Seurat Object. Dynamically renames specified metadata columns in a Seurat object, particularly those  prefixed with "predicted." and the "mapping.score" column, by applying a new prefix  that combines a user-defined prefix and a reference name. 

- #### 25 `renameSmallCategories()`

  Rename Small Categories in Seurat Object Metadata. This function renames categories within a specified identity column of a  Seurat object's metadata that have fewer cells than a specified minimum threshold.  Categories below this threshold are renamed to a common name, typically "unclear",  to clean up small, potentially noisy categories.

- #### 26 `.metaColnames()`

  Transfer labels from a reference Seurat object to a query Seurat object. Function to transfer labels from a reference Seurat object to a query Seurat object  using anchoring and transfer data methods from the Seurat package. It then visualizes the  reference and the combined objects using Uniform Manifold Approximation and Projection (UMAP). 

- #### 27 `matchBestIdentity()`

  Match and Translate Best Identity. This function matches the best identity from `ident_to_rename` to `reference_ident` in an object,  in other words, it replaces original categories with the most frequent ones from the reference,  hence helps to filter out less important categories. 



## List of Functions in Seurat.utils.less.used.R (5) 

Updated: 2024/03/22 15:41

- #### 1 `plot.UMAP.tSNE.sidebyside()`

  plot.UMAP.tSNE.sidebyside. Plot a UMAP and tSNE side by side.

- #### 2 `umapNamedClusters()`

  Plot and Save UMAP without legend. Generates a UMAP plot colored by a specified metadata column and saves the plot to a file. 

- #### 3 `AutoNumber.by.PrinCurve()`

  AutoNumber.by.PrinCurve. Relabel cluster numbers along the principal curve of 2 UMAP (or tSNE) dimensions. #

- #### 4 `load10Xv3()`

  Load 10X Genomics Version 3 Data. Loads 10X Genomics data from a specified directory containing output folders for raw and filtered data.  This function is designed to handle data from 10X Genomics Chromium Single Cell technologies (version 3). 

- #### 5 `Convert10Xfolders.old()`

  Convert10Xfolders.old. This function takes a parent directory with a number of subfolders, each containing the standard output of 10X Cell Ranger. It (1) loads the filtered data matrices, (2) converts them to Seurat objects, and (3) saves them as .RDS files.




# Usage

You can use most functions at relevant steps of a standard Seurat analysis.

We are preparing a vignette.




-----------
[Get Seurat.utils](https://github.com/vertesy/Seurat.utils). Vertesy, 2024. [![DOI](https://zenodo.org/badge/248721133.svg)](https://zenodo.org/badge/latestdoi/248721133) 
*If you use these functions, please star the repo, or cite via `DOI`. Thanks!*

<br>

