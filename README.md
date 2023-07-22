[![DOI](https://zenodo.org/badge/248721133.svg)](https://zenodo.org/badge/latestdoi/248721133) *If you use these functions, please star the repo, or cite via `DOI`. Thanks!*

# Seurat.utils

`Seurat.utils` Is a collection of utility functions for Seurat. Functions allow the automation / multiplexing of plotting, 3D plotting, visualisation of statistics & QC, interaction with the Seurat object, etc.  Some functionalities require functions from [CodeAndRoll2](https://github.com/vertesy/CodeAndRoll2), [ReadWriter](https://github.com/vertesy/ReadWriter), [Stringendo](https://github.com/vertesy/Stringendo), [ggExpressDev](https://github.com/vertesy/ggExpressDev), [MarkdownReports](https://github.com/vertesy/MarkdownReports), and the [Rocinante](https://github.com/vertesy/Rocinante) (See installation).



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

## Seurat.Utils.R
Updated: 2023/07/22 12:18

- #### 1 `parallel.computing.by.future()`

  Run gc(), load multi-session computing and extend memory limits.

- #### 2 `IntersectWithExpressed()`


- #### 3 `SmallestNonAboveX()`

  replace small values with the next smallest value found, which is >X. #

- #### 4 `AreTheseCellNamesTheSame()`

  Compare two character vectors (e.g.: cell IDs) how much they overlap and plot a Venn Diagramm.

- #### 5 `getProject()`

  Try to get the project name you are wokring on in Rstudio.

- #### 6 `PlotFilters()`

  Plot filtering threshold and distributions, using four panels to highlight the relation between Gene- and UMI-count, ribosomal- and mitochondrial-content.

- #### 7 `seu.PC.var.explained()`

  Determine percent of variation associated with each PC. For normal prcomp objects, see: PCA.percent.var.explained().

- #### 8 `seu.plot.PC.var.explained()`

  Plot the percent of variation associated with each PC.

- #### 9 `Percent.in.Trome()`

  Gene expression as fraction of all UMI's

- #### 10 `gene.expression.level.plots()`

  Histogram of gene expression levels.

- #### 11 `PrctCellExpringGene()`

  Function to calculate the proportion of cells expressing a given set of genes.

- #### 12 `ww.calc_helper()`

  Helper function for PrctCellExpringGene() to calculate the proportion of cells in a Seurat object that express a given gene.

- #### 13 `scBarplot.FractionAboveThr()`

  Create a bar plot showing the fraction of cells, within each cluster, that exceed a certain threshold based on a metadata column.

- #### 14 `scBarplot.FractionBelowThr()`

  Create a bar plot showing the fraction of cells, within each cluster, that are below a certain threshold based on a metadata column.

- #### 15 `getClusterNames()`

  Rename clustering in a Seurat object.

- #### 16 `GetClusteringRuns()`

  Get Clustering Runs: metadata column names #

- #### 17 `GetNamedClusteringRuns()`

  Get Clustering Runs: metadata column names #

- #### 18 `GetOrderedClusteringRuns()`

  Get Clustering Runs: metadata column names #

- #### 19 `GetNumberOfClusters()`

  Get Number Of Clusters #

- #### 20 `calc.cluster.averages()`

  Calculates the average of a metadata column (numeric) per cluster.

- #### 21 `plot.expression.rank.q90()`

  Plot gene expression based on the expression at the 90th quantile (so you will not lose genes expressed in few cells).

- #### 22 `set.mm()`

  Helps to find metadata columns. It creates a list with the names of of 'obj@meta.data'.

- #### 23 `recall.all.genes()`

  all.genes set by calc.q99.Expression.and.set.all.genes() #

- #### 24 `recall.meta.tags.n.datasets()`

  Recall  meta.tags from obj@misc to "meta.tags" in the global environment.

- #### 25 `recall.parameters()`

  Recall parameters from obj@misc to "p" in the global environment.

- #### 26 `recall.genes.ls()`

  Recall genes.ls from obj@misc to "genes.ls" in the global environment.

- #### 27 `save.parameters()`

  Save parameters to obj@misc$p

- #### 28 `subsetSeuObj()`

  Subset a compressed Seurat object and save it in the working directory.

- #### 29 `subsetSeuObj.and.Save()`

  Subset a compressed Seurat Obj and save it in wd. #

- #### 30 `subsetSeuObj.ident.class()`

  Subset a Seurat Obj to a given column

- #### 31 `Downsample.Seurat.Objects()`

  Downsample a list of Seurat objects

- #### 32 `Downsample.Seurat.Objects.PC()`

  Downsample a list of Seurat objects, by fraction

- #### 33 `remove.residual.small.clusters()`

  E.g.: after subsetting often some residual cells remain in clusters originally defined in the full dataset.

- #### 34 `drop.levels.Seurat()`

  Drop unused levels from factor variables in a Seurat object.

- #### 35 `remove_clusters_and_drop_levels()`

  This function removes residual small clusters from specified Seurat objects and drops levels in factor-like metadata.

- #### 36 `remove.cells.by.UMAP()`

  This function applies a cutoff in the specified dimension of a given dimension reduction (UMAP, PCA, or t-SNE) to remove cells.

- #### 37 `FlipReductionCoordinates()`

  Flip reduction coordinates (like UMAP upside down).

- #### 38 `AutoNumber.by.UMAP()`

  Relabel cluster numbers along a UMAP (or tSNE) axis #

- #### 39 `AutoNumber.by.PrinCurve()`

  Relabel cluster numbers along the principal curve of 2 UMAP (or tSNE) dimensions. #

- #### 40 `Add.DE.combined.score()`

  Add a combined score to differential expression (DE) results. The score is calculated as log-fold change (LFC) times negative logarithm of scaled p-value (LFC * -log10( p_cutoff / pval_scaling )).

- #### 41 `StoreTop25Markers()`

  Save the top 25 makers based on `avg_log2FC` output table of `FindAllMarkers()` (df_markers) under `@misc$df.markers$res...`. By default, it rounds up insignificant digits up to 3. #

- #### 42 `StoreAllMarkers()`

  Save the output table of `FindAllMarkers()` (df_markers) under `@misc$df.markers$res...`. By default, it rounds up insignificant digits up to 3. #

- #### 43 `GetTopMarkersDF()`

  Get the vector of N most diff. exp. genes. #

- #### 44 `GetTopMarkers()`

  Get the vector of N most diff. exp. genes. #

- #### 45 `AutoLabelTop.logFC()`

  Create a new "named identity" column in the metadata of a Seurat object, with `Ident` set to a clustering output matching the `res` parameter of the function. It requires the output table of `FindAllMarkers()`. If you used `StoreAllMarkers()` is stored under `@misc$df.markers$res...`, which location is assumed by default. #

- #### 46 `scEnhancedVolcano()`

  Creates a new "named identity" column in the metadata of a Seurat object, setting 'Ident' to a clustering output matching the 'res' parameter. This function requires the output table of `FindAllMarkers()`. If you used `StoreAllMarkers()`, the output is stored under `@misc$df.markers$res...`, which is the default location.

- #### 47 `BulkGEScatterPlot()`

  Plots scatterplots of bulk gene expression to identify differentially expressed genes across conditions.

- #### 48 `get.clustercomposition()`

  Get cluster composition: which datasets contribute to each cluster?

- #### 49 `scBarplot.CellFractions()`

  Generates a bar plot of cell fractions per cluster.

- #### 50 `  scBarplot.CellsPerCluster()`

  Barplot the Fraction of cells per cluster. (dupl?)

- #### 51 `scBarplot.CellsPerObject()`

  Creates a bar plot for the number of cells per object from a list of Seurat objects.

- #### 52 `plot.clust.size.distr()`

  Creates a bar plot or histogram of the cluster size distribution from a given Seurat object.

- #### 53 `gg_color_hue()`

  Emulates the default color palette of ggplot2. Source:  https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette

- #### 54 `getDiscretePalette()`

  Generate a discrete color palette.

- #### 55 `getClusterColors()`

  get Seurat's cluster colors.

- #### 56 `SeuratColorVector()`

  Recall a Seurat color vector.

- #### 57 `plot.GeneExpHist()`

  Generates and optionally saves a scatter plot of two features from a Seurat object.

- #### 58 `qUMAP()`

  The quickest way to draw a gene expression UMAP.

- #### 59 `clUMAP()`

  The quickest way to draw a clustering result UMAP.

- #### 60 `umapNamedClusters()`

  Plot and save umap based on a metadata column. #

- #### 61 `umapHiLightSel()`

  Generates a UMAP plot from a Seurat object with a subset of cells highlighted.

- #### 62 `multiFeaturePlot.A4()`

  Save multiple FeaturePlots, as jpeg, on A4 for each gene, which are stored as a list of gene names.

- #### 63 `multiFeatureHeatmap.A4()`

  Save multiple FeatureHeatmaps from a list of genes on A4 jpeg.

- #### 64 `plot.UMAP.tSNE.sidebyside()`

  Plot a UMAP and tSNE side by side.

- #### 65 `PlotTopGenesPerCluster()`

  Plot the top N diff. exp. genes in each cluster.

- #### 66 `qQC.plots.BrainOrg()`

  Quickly plot key QC markers in brain organoids

- #### 67 `qMarkerCheck.BrainOrg()`

  Quickly plot key markers in brain organoids

- #### 68 `PlotTopGenes()`

  Plot the highest expressed genes on umaps, in a subfolder. Requires calling calc.q99.Expression.and.set.all.genes before. #

- #### 69 `DimPlot.ClusterNames()`

  Plot UMAP with Cluster names. #

- #### 70 `save2umaps.A4()`

  Save 2 umaps on 1 A4

- #### 71 `save4umaps.A4()`

  Save 4 umaps on 1 A4

- #### 72 `qqSaveGridA4()`

  Saves a grid of 2 or 4 ggplot objects onto an A4 page.

- #### 73 `ww.check.if.3D.reduction.exist()`

  ww.check.if.3D.reduction.exist in backup slot #

- #### 74 `ww.check.quantile.cutoff.and.clip.outliers()`

  Function to check a specified quantile cutoff and clip outliers from a given expression vector.

- #### 75 `plot3D.umap.gene()`

  Plot a 3D umap with gene expression. Uses plotly. Based on github.com/Dragonmasterx87.

- #### 76 `plot3D.umap()`

  Plot a 3D umap based on one of the metadata columns. Uses plotly. Based on github.com/Dragonmasterx87.

- #### 77 `SavePlotlyAsHtml()`

  Save a Plotly 3D scatterplot as an HTML file.

- #### 78 `BackupReduction()`

  Backup UMAP to `obj@misc$reductions.backup` from `obj@reductions$umap`. #

- #### 79 `SetupReductionsNtoKdimensions()`

  Function to compute dimensionality reductions for a given Seurat object and backup the computed reductions.

- #### 80 `RecallReduction()`

  Set active UMAP to `obj@reductions$umap` from `obj@misc$reductions.backup`. #

- #### 81 `Annotate4Plotly3D()`

  Create annotation labels for 3D plots. Source https://plot.ly/r/text-and-annotations/#3d-annotations.

- #### 82 `Plot3D.ListOfGenes()`

  Plot and save list of 3D UMAP or tSNE plots using plotly.

- #### 83 `Plot3D.ListOfCategories()`

  This function plots and saves a list of 3D UMAP or tSNE plots using plotly.

- #### 84 `# sparse.cor4()`

  Calculate a sparse correlation matrix.

- #### 85 `Calc.Cor.Seurat()`

  Calculate gene correlation on a Seurat object.

- #### 86 `plot.Gene.Cor.Heatmap()`

  Plot a gene correlation heatmap.

- #### 87 `prefix_cells_seurat()`

  This function adds prefixes from 'obj_IDs' to cell names in Seurat S4 objects from 'ls_obj' 

- #### 88 `find_prefix_in_cell_IDs()`

  This function checks if a prefix has been added to the standard cell-IDs (16 characters of A,T,C,G)  in a Seurat object. If so, it prints the number of unique prefixes found,  issues a warning if more than one unique prefix is found, and returns the identified prefix(es).

- #### 89 `seu.Make.Cl.Label.per.cell()`

  Take a named vector (of e.g. values ="gene names", names = clusterID), and a vector of cell-IDs and make a vector of "GeneName.ClusterID".

- #### 90 `GetMostVarGenes()`

  Get the N most variable Genes

- #### 91 `gene.name.check()`

  Check gene names in a seurat object, for naming conventions (e.g.: mitochondrial reads have - or .). Use for reading .mtx & writing .rds files. #

- #### 92 `check.genes()`

  Check if a gene name exists in a Seurat object, or in HGNC?

- #### 93 `fixZeroIndexing.seurat()`

  Fix zero indexing in seurat clustering, to 1-based indexing. replace zero indexed clusternames.

- #### 94 `CalculateFractionInTrome()`

  This function calculates the fraction of a set of genes within the full transcriptome of each cell.

- #### 95 `AddNewAnnotation()`

  This function creates a new metadata column based on an existing metadata column and a list of mappings (name <- IDs).

- #### 96 `whitelist.subset.ls.Seurat()`

  Subsets cells in a list of Seurat objects based on an externally provided list of cell IDs.

- #### 97 `FindCorrelatedGenes()`

  Find correlated genes in a Seurat object

- #### 98 `UpdateGenesSeurat()`

  Update genes symbols that are stored in a Seurat object. It returns a data frame. The last column are the updated gene names.

- #### 99 `  check_and_rename()`

  Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data. #

- #### 100 `RemoveGenesSeurat()`

  Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes metadata; obj@assays$RNA@counts, @data and @scale.data. #

- #### 101 `HGNC.EnforceUnique()`

  Enforce Unique names after HGNC symbol update.

- #### 102 `GetUpdateStats()`

  Plot the Symbol-update statistics. Works on the data frame returned by `UpdateGenesSeurat()`. #

- #### 103 `PlotUpdateStats()`

  Creates a scatter plot of update statistics.

- #### 104 `calculate.observable.multiplet.rate.10X.LT()`

  Calculate the observable multiplet rate for 10X standard lane.

- #### 105 `SNP.demux.fix.GT.table()`

  This function cleans and standardizes a Genotype assignment table obtained from the SoupOrCell tool.

- #### 106 `Convert10Xfolders()`

  This function takes a parent directory with a number of subfolders, each containing the standard output of 10X Cell Ranger. It (1) loads the filtered data matrices, (2) converts them to Seurat objects, and (3) saves them as .RDS files.

- #### 107 `ConvertDropSeqfolders()`

  This function takes a parent directory with a number of subfolders, each containing the standard output of 10X Cell Ranger. It (1) loads the filtered data matrices, (2) converts them to Seurat objects, and (3) saves them as .RDS files.

- #### 108 `LoadAllSeurats()`

  This function loads all Seurat objects found in a directory. It also works with symbolic links (but not with aliases).

- #### 109 `read10x()`

  This function reads a 10X dataset from gzipped matrix.mtx, features.tsv and barcodes.tsv files.

- #### 110 `load10Xv3()`

  Load 10X output folders.

- #### 111 `saveRDS.compress.in.BG()`

  Save and RDS object and compress resulting file in the background using system(gzip). OS X or unix.

- #### 112 `isave.RDS()`

  Save an RDS object, using a faster and efficient compression method that runs in the background.

- #### 113 `isave.image()`

  Save an image of the current workspace using a faster and efficient compression method that runs in the background.

- #### 114 `qsave.image()`

  Faster saving of workspace, and compression outside R, when it can run in the background. Seemingly quite CPU hungry and not very efficient compression. #

- #### 115 `clip10Xcellname()`

  Clip all suffices after underscore (10X adds it per chip-lane, Seurat adds in during integration). #

- #### 116 `make10Xcellname()`

  Add a suffix to cell names, so that it mimics the lane-suffix, e.g.: "_1". #

- #### 117 `plotTheSoup()`

  Plot stats about the ambient RNA content in a 10X experiment.

- #### 118 `jJaccardIndexVec()`

  Calculate jaccard similarity for 2 vecotrs. Helper to jPairwiseJaccardIndexList.

- #### 119 `jPairwiseJaccardIndexList()`

  Create a pairwise jaccard similarity matrix across all combinations of columns in binary.presence.matrix. Modified from: https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/ #

- #### 120 `jPresenceMatrix()`

  Make a binary presence matrix from a list. Source: https://stackoverflow.com/questions/56155707/r-how-to-create-a-binary-relation-matrix-from-a-list-of-strings #

- #### 121 `jJaccardIndexBinary()`

  Calculate Jaccard Index. Modified from: https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/ #

------------------------------------------------------------------------------------------------------------
## Seurat.Utils.Metadata.R
Updated: 2023/07/22 12:18

- #### 1 `meta_col_exists()`

  This function checks whether a given column exists in the meta.data of a Seurat object.

- #### 2 `getMedianMetric()`

  Get the median values of different columns in meta.data, can iterate over a list of Seurat objects.

- #### 3 `add.meta.tags()`

  Add metadata tags to a Seurat object dataset.

- #### 4 `add.meta.fraction()`

  Add a new metadata column to a Seurat object, representing the fraction of a gene set in the transcriptome (expressed as a percentage).

- #### 5 `seu.RemoveMetadata()`

  Remove specified metadata columns from a Seurat object.

- #### 6 `getMetadataColumn()`

  Retrieves a specified metadata column from a Seurat object and returns it as a named vector.

- #### 7 `create.metadata.vector()`

  Adds a new metadata column to a Seurat object.

- #### 8 `seu.map.and.add.new.ident.to.meta()`

  Adds a new metadata column to a Seurat object based on an identity mapping table.

- #### 9 `getCellIDs.from.meta()`

  Retrieves cell IDs from a specified metadata column of a Seurat object, where the cell ID matches a provided list of values. The matching operation uses the `%in%` operator.

- #### 10 `seu.add.meta.from.table()`

  Add multiple new metadata columns to a Seurat object from a table. #

- #### 11 `sampleNpc()`

  This function samples a specified percentage of a dataframe (specifically a subset of the metadata of a Seurat object) and returns the corresponding cell IDs.

- #### 12 `calc.q99.Expression.and.set.all.genes()`

  Calculate the gene expression of the e.g.: 90th quantile (expression in the top 10% cells). #

- #### 13 `fix.orig.ident()`

  Remove the string "filtered_feature_bc_matrix." from "orig.ident". Helper function.

- #### 14 `Create.MiscSlot()`

  It is just a reminder to use calc.q99.Expression.and.set.all.genes to create the all.genes variable

- #### 15 `transfer_labels_seurat()`

  Function to transfer labels from a reference Seurat object to a query Seurat object using anchoring and transfer data methods from the Seurat package.  It then visualizes the reference and the combined objects using Uniform Manifold Approximation and Projection (UMAP). 

- #### 16 `match_best_identity()`

  This function matches the best identity from `ident_from` to `ident_to` in an object,  updates the metadata of the object with this new identity and returns the updated object. Additionally,  it generates a UMAP plot based on the new identity. The function replaces original categories with  the most frequent ones, hence helps to filter out less important categories. 

- #### 17 `replace_by_most_frequent_categories()`

  This function replaces each category in a query column of a data frame with the most    frequently corresponding category in a reference column. It calculates the assignment quality,    reports it, and optionally plots it. 

- #### 18 `plot.Metadata.Cor.Heatmap()`

  Plots a heatmap of metadata correlation values.




# Usage

You can use most functions at relevant steps of a standard Seurat analysis.


-----------
[Get Seurat.utils](https://github.com/vertesy/Seurat.utils). Vertesy, 2021. [![DOI](https://zenodo.org/badge/248721133.svg)](https://zenodo.org/badge/latestdoi/248721133) 
*If you use these functions, please star the repo, or cite via `DOI`. Thanks!*

<br>

