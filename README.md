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

## Seurat.Utils.R (127) 
Updated: 2023/11/24 16:40

- #### 1 `parallel.computing.by.future()`

  parallel.computing.by.future. Run gc(), load multi-session computing and extend memory limits.

- #### 2 `IntersectWithExpressed()`

  IntersectWithExpressed. Intersect a set of genes with genes found in the Seurat object.

- #### 3 `SmallestNonAboveX()`

  SmallestNonAboveX. replace small values with the next smallest value found, which is >X.

- #### 4 `AreTheseCellNamesTheSame()`

  AreTheseCellNamesTheSame. Compare two character vectors (e.g.: cell IDs) how much they overlap and plot a Venn Diagram.

- #### 5 `getProject()`

  getProject. Try to get the project name you are wokring on in Rstudio.

- #### 6 `Create.MiscSlot()`

  Shorten Clustering Names. This function takes in a string representing a clustering name,  and shortens it according to specific rules. It replaces "snn_res." with "",  "best.matching.names" with "bmatch", "ordered" with "ord",  "ManualNames" with "mNames", and ".long" at the end of the string with ".L". 

- #### 7 `calc.q99.Expression.and.set.all.genes()`

  calc.q99.Expression.and.set.all.genes. Calculate the gene expression of the e.g.: 90th quantile (expression in the top 10% cells). #

- #### 8 `PlotFilters()`

  PlotFilters. Plot filtering threshold and distributions, using four panels to highlight the relation between Gene- and UMI-count, ribosomal- and mitochondrial-content.

- #### 9 `seu.PC.var.explained()`

  PCA percent of variation associated with each PC. Determine percent of variation associated with each PC. For normal prcomp objects, see: PCA.percent.var.explained().

- #### 10 `seu.plot.PC.var.explained()`

  seu.plot.PC.var.explained. Plot the percent of variation associated with each PC.

- #### 11 `Percent.in.Trome()`

  Percent.in.Trome. Gene expression as fraction of all UMI's

- #### 12 `gene.expression.level.plots()`

  gene.expression.level.plots. Histogram of gene expression levels.

- #### 13 `PrctCellExpringGene()`

  PrctCellExpringGene. Function to calculate the proportion of cells expressing a given set of genes.

- #### 14 `ww.calc_helper()`

  ww.calc_helper. Helper function for PrctCellExpringGene() to calculate the proportion of cells in a Seurat object that express a given gene.

- #### 15 `scBarplot.FractionAboveThr()`

  scBarplot.FractionAboveThr. Create a bar plot showing the fraction of cells, within each cluster, that exceed a certain threshold based on a metadata column.

- #### 16 `scBarplot.FractionBelowThr()`

  scBarplot.FractionBelowThr. Create a bar plot showing the fraction of cells, within each cluster, that are below a certain threshold based on a metadata column.

- #### 17 `getClusterNames()`

  RenameClustering. Rename clustering in a Seurat object.

- #### 18 `GetClusteringRuns()`

  GetClusteringRuns. Get Clustering Runs: metadata column names #

- #### 19 `GetNamedClusteringRuns()`

  GetNamedClusteringRuns. Get Clustering Runs: metadata column names #

- #### 20 `GetOrderedClusteringRuns()`

  GetOrderedClusteringRuns. Get Clustering Runs: metadata column names #

- #### 21 `GetNumberOfClusters()`

  GetNumberOfClusters. Get Number Of Clusters #

- #### 22 `calc.cluster.averages()`

  calc.cluster.averages. Calculates the average of a metadata column (numeric) per cluster.

- #### 23 `plot.expression.rank.q90()`

  plot.expression.rank.q90. Plot gene expression based on the expression at the 90th quantile (so you will not lose genes expressed in few cells).

- #### 24 `set.mm()`

  set.mm. Helps to find metadata columns. It creates a list with the names of of 'obj@meta.data'.

- #### 25 `ww.get.1st.Seur.element()`

  Get the First Seurat Object from a List of Seurat Objects. #' If provided with a list of Seurat objects, this function returns the first  Seurat object in the list. If the input is a single Seurat object, it returns  the object itself. It is assumed that all elements of the list are Seurat  objects if the input is a list. 

- #### 26 `recall.all.genes()`

  recall.all.genes. all.genes set by calc.q99.Expression.and.set.all.genes() #

- #### 27 `recall.meta.tags.n.datasets()`

  recall.meta.tags.n.datasets. Recall  meta.tags from obj@misc to "meta.tags" in the global environment.

- #### 28 `recall.parameters()`

  recall.parameters. Recall parameters from obj@misc to "p" in the global environment.

- #### 29 `recall.genes.ls()`

  recall.genes.ls. Recall genes.ls from obj@misc to "genes.ls" in the global environment.

- #### 30 `save.parameters()`

  save.parameters. Save parameters to obj@misc$p

- #### 31 `create_scCombinedMeta()`

  Create Single-Cell Metadata Object for a collection of Seurat Objects. This function creates a metadata object to correspond to a list of    single-cell experiments, for storing parent level information.    It initializes the object with the experiment and project name, and the    creation date. The created object is of class 'scMetadata_class'.

- #### 32 `subsetSeuObj()`

  subsetSeuObj. Subset a compressed Seurat object and save it in the working directory.

- #### 33 `subsetSeuObj.and.Save()`

  subsetSeuObj.and.Save. Subset a compressed Seurat Obj and save it in wd. #

- #### 34 `subsetSeuObj.ident.class()`

  subsetSeuObj.ident.class. Subset a Seurat Obj to a given column

- #### 35 `Downsample.Seurat.Objects()`

  Downsample.Seurat.Objects. Downsample a list of Seurat objects

- #### 36 `Downsample.Seurat.Objects.PC()`

  Downsample.Seurat.Objects.PC. Downsample a list of Seurat objects, by fraction

- #### 37 `remove.residual.small.clusters()`

  remove.residual.small.clusters. E.g.: after subsetting often some residual cells remain in clusters originally defined in the full dataset.

- #### 38 `dropLevelsSeurat()`

  dropLevelsSeurat. Drop unused levels from factor variables in a Seurat object.

- #### 39 `remove_clusters_and_drop_levels()`

  Remove Clusters and Drop Levels. This function removes residual small clusters from specified Seurat objects and drops levels in factor-like metadata.

- #### 40 `remove.cells.by.UMAP()`

  Remove Cells by Dimension Reduction. This function applies a cutoff in the specified dimension of a given dimension reduction (UMAP, PCA, or t-SNE) to remove cells.

- #### 41 `FlipReductionCoordinates()`

  FlipReductionCoordinates. Flip reduction coordinates (like UMAP upside down).

- #### 42 `AutoNumber.by.UMAP()`

  AutoNumber.by.UMAP. Relabel cluster numbers along a UMAP (or tSNE) axis #

- #### 43 `AutoNumber.by.PrinCurve()`

  AutoNumber.by.PrinCurve. Relabel cluster numbers along the principal curve of 2 UMAP (or tSNE) dimensions. #

- #### 44 `Add.DE.combined.score()`

  Add.DE.combined.score. Add a combined score to differential expression (DE) results. The score is calculated as log-fold change (LFC) times negative logarithm of scaled p-value (LFC * -log10( p_cutoff / pval_scaling )).

- #### 45 `StoreTop25Markers()`

  StoreTop25Markers. Save the top 25 makers based on `avg_log2FC` output table of `FindAllMarkers()` (df_markers) under `@misc$df.markers$res...`. By default, it rounds up insignificant digits up to 3. #

- #### 46 `StoreAllMarkers()`

  StoreAllMarkers. Save the output table of `FindAllMarkers()` (df_markers) under  `@misc$df.markers$res...`. By default, it rounds up insignificant digits up to 3.

- #### 47 `GetTopMarkersDF()`

  GetTopMarkersDF. Get the vector of N most diff. exp. genes. #

- #### 48 `GetTopMarkers()`

  GetTopMarkers. Get the vector of N most diff. exp. genes. #

- #### 49 `AutoLabelTop.logFC()`

  AutoLabelTop.logFC. Create a new "named identity" column in the metadata of a Seurat object,  with `Ident` set to a clustering output matching the `res` parameter of the function.  t requires the output table of `FindAllMarkers()`.   If you used `StoreAllMarkers()` is stored under `@misc$df.markers$res...`,   which location is assumed by default.

- #### 50 `AutoLabel.KnownMarkers()`

  AutoLabel.KnownMarkers. Creates a new "named identity" column in the metadata of a Seurat object,   setting 'Ident' to a clustering output matching the 'res' parameter.   This function requires the output table of `FindAllMarkers()`.  If you used `StoreAllMarkers()`, the output is stored under `@misc$df.markers$res...`,  which is the default location.

- #### 51 `scEnhancedVolcano()`

  scEnhancedVolcano. This function creates an enhanced volcano plot. 

- #### 52 `BulkGEScatterPlot()`

  BulkGEScatterPlot. Plots scatterplots of bulk gene expression to identify differentially expressed genes across conditions.

- #### 53 `get.clustercomposition()`

  get.clustercomposition. Get cluster composition: which datasets contribute to each cluster?

- #### 54 `scBarplot.CellFractions()`

  Generate Barplot of Cell Fractions. This function generates a bar plot of cell fractions per cluster  from a Seurat object. It offers the option to downsample data, which equalizes  the number of cells in each group to the number in the smallest group.  The plot's bars are grouped by one variable and filled by another. 

- #### 55 `  scBarplot.CellsPerCluster()`

  scBarplot.CellsPerCluster. Barplot the Fraction of cells per cluster. (dupl?)

- #### 56 `scBarplot.CellsPerObject()`

  scBarplot.CellsPerObject. Creates a bar plot for the number of cells per object from a list of Seurat objects.

- #### 57 `plot.clust.size.distr()`

  plot.clust.size.distr. Creates a bar plot or histogram of the cluster size distribution from a given Seurat object.

- #### 58 `gg_color_hue()`

  gg_color_hue. Emulates the default color palette of ggplot2. Source:  https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette

- #### 59 `getDiscretePalette()`

  getDiscretePalette. Generate a discrete color palette.

- #### 60 `getClusterColors()`

  getClusterColors. get Seurat's cluster colors.

- #### 61 `SeuratColorVector()`

  SeuratColorVector. Recall a Seurat color vector.

- #### 62 `plotGeneExpHist()`

  qFeatureScatter. Generates and optionally saves a scatter plot of two features from a Seurat object.

- #### 63 `qUMAP()`

  qUMAP. The quickest way to draw a gene expression UMAP.

- #### 64 `clUMAP()`

  clUMAP. The quickest way to draw a clustering result UMAP.

- #### 65 `umapNamedClusters()`

  umapNamedClusters. Plot and save umap based on a metadata column. #

- #### 66 `umapHiLightSel()`

  umapHiLightSel. Generates a UMAP plot from a Seurat object with a subset of cells highlighted.

- #### 67 `multiFeaturePlot.A4()`

  multiFeaturePlot.A4. Save multiple FeaturePlots, as jpeg, on A4 for each gene, which are stored as a list of gene names.

- #### 68 `multiFeatureHeatmap.A4()`

  multiFeatureHeatmap.A4. Save multiple FeatureHeatmaps from a list of genes on A4 jpeg.

- #### 69 `plot.UMAP.tSNE.sidebyside()`

  plot.UMAP.tSNE.sidebyside. Plot a UMAP and tSNE side by side.

- #### 70 `PlotTopGenesPerCluster()`

  PlotTopGenesPerCluster. Plot the top N diff. exp. genes in each cluster.

- #### 71 `qQC.plots.BrainOrg()`

  qQC.plots.BrainOrg. Quickly plot key QC markers in brain organoids

- #### 72 `qMarkerCheck.BrainOrg()`

  qMarkerCheck.BrainOrg. Quickly plot key markers in brain organoids

- #### 73 `PlotTopGenes()`

  PlotTopGenes. Plot the highest expressed genes on umaps, in a subfolder. Requires calling calc.q99.Expression.and.set.all.genes before. #

- #### 74 `DimPlot.ClusterNames()`

  DimPlot.ClusterNames. Plot UMAP with Cluster names. #

- #### 75 `save2umaps.A4()`

  save2umaps.A4. Save 2 umaps on 1 A4

- #### 76 `save4umaps.A4()`

  save4umaps.A4. Save 4 umaps on 1 A4

- #### 77 `qqSaveGridA4()`

  qqSaveGridA4. Saves a grid of 2 or 4 ggplot objects onto an A4 page.

- #### 78 `ww.check.if.3D.reduction.exist()`

  ww.check.if.3D.reduction.exist. ww.check.if.3D.reduction.exist in backup slot #

- #### 79 `ww.check.quantile.cutoff.and.clip.outliers()`

  ww.check.quantile.cutoff.and.clip.outliers. Function to check a specified quantile cutoff and clip outliers from a given expression vector.

- #### 80 `plot3D.umap.gene()`

  plot3D.umap.gene. Plot a 3D umap with gene expression. Uses plotly. Based on github.com/Dragonmasterx87.

- #### 81 `plot3D.umap()`

  plot3D.umap. Plot a 3D umap based on one of the metadata columns. Uses plotly. Based on github.com/Dragonmasterx87.

- #### 82 `SavePlotlyAsHtml()`

  SavePlotlyAsHtml. Save a Plotly 3D scatterplot as an HTML file.

- #### 83 `BackupReduction()`

  BackupReduction. Backup UMAP to `obj@misc$reductions.backup` from `obj@reductions$umap`. #

- #### 84 `SetupReductionsNtoKdimensions()`

  SetupReductionsNtoKdimensions. Function to compute dimensionality reductions for a given Seurat object and backup the computed reductions.

- #### 85 `RecallReduction()`

  RecallReduction. Set active UMAP to `obj@reductions$umap` from `obj@misc$reductions.backup`. #

- #### 86 `Annotate4Plotly3D()`

  Annotate4Plotly3D. Create annotation labels for 3D plots. Source https://plot.ly/r/text-and-annotations/#3d-annotations.

- #### 87 `Plot3D.ListOfGenes()`

  Plot3D.ListOfGenes. Plot and save list of 3D UMAP or tSNE plots using plotly.

- #### 88 `Plot3D.ListOfCategories()`

  Plot3D.ListOfCategories. This function plots and saves a list of 3D UMAP or tSNE plots using plotly.

- #### 89 `# sparse.cor4()`

  sparse.cor. Calculate a sparse correlation matrix.

- #### 90 `Calc.Cor.Seurat()`

  Calc.Cor.Seurat. Calculate gene correlation on a Seurat object.

- #### 91 `plot.Gene.Cor.Heatmap()`

  plot.Gene.Cor.Heatmap. Plot a gene correlation heatmap.

- #### 92 `prefix_cells_seurat()`

  prefix_cells_seurat. This function adds prefixes from 'obj_IDs' to cell names in Seurat S4 objects from 'ls_obj' 

- #### 93 `find_prefix_in_cell_IDs()`

  Check Prefix in Seurat Object Cell IDs. This function checks if a prefix has been added to the standard cell-IDs (16 characters of A,T,C,G)  in a Seurat object. If so, it prints the number of unique prefixes found,  issues a warning if more than one unique prefix is found, and returns the identified prefix(es).

- #### 94 `seu.Make.Cl.Label.per.cell()`

  seu.Make.Cl.Label.per.cell. Take a named vector (of e.g. values ="gene names", names = clusterID), and a vector of cell-IDs and make a vector of "GeneName.ClusterID".

- #### 95 `GetMostVarGenes()`

  GetMostVarGenes. Get the N most variable Genes

- #### 96 `gene.name.check()`

  gene.name.check. Check gene names in a seurat object, for naming conventions (e.g.: mitochondrial reads have - or .). Use for reading .mtx & writing .rds files. #

- #### 97 `check.genes()`

  check.genes. Check if a gene name exists in a Seurat object, or in HGNC?

- #### 98 `fixZeroIndexing.seurat()`

  fixZeroIndexing.seurat. Fix zero indexing in seurat clustering, to 1-based indexing. replace zero indexed clusternames.

- #### 99 `CalculateFractionInTrome()`

  CalculateFractionInTranscriptome. This function calculates the fraction of a set of genes within the full transcriptome of each cell.

- #### 100 `AddNewAnnotation()`

  AddNewAnnotation. This function creates a new metadata column based on an existing metadata column and a list of mappings (name <- IDs).

- #### 101 `whitelist.subset.ls.Seurat()`

  whitelist.subset.ls.Seurat. Subsets cells in a list of Seurat objects based on an externally provided list of cell IDs.

- #### 102 `FindCorrelatedGenes()`

  FindCorrelatedGenes. Find correlated genes in a Seurat object

- #### 103 `UpdateGenesSeurat()`

  UpdateGenesSeurat. Update genes symbols that are stored in a Seurat object. It returns a data frame. The last column are the updated gene names.

- #### 104 `  check_and_rename()`

  RenameGenesSeurat. Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data. #

- #### 105 `RemoveGenesSeurat()`

  RemoveGenesSeurat. Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes metadata; obj@assays$RNA@counts, @data and @scale.data. #

- #### 106 `HGNC.EnforceUnique()`

  HGNC.EnforceUnique. Enforce Unique names after HGNC symbol update.

- #### 107 `GetUpdateStats()`

  GetUpdateStats. Plot the Symbol-update statistics. Works on the data frame returned by `UpdateGenesSeurat()`. #

- #### 108 `PlotUpdateStats()`

  PlotUpdateStats. Creates a scatter plot of update statistics.

- #### 109 `calculate.observable.multiplet.rate.10X.LT()`

  calculate.observable.multiplet.rate.10X.LT. Calculate the observable multiplet rate for 10X standard lane.

- #### 110 `SNP.demux.fix.GT.table()`

  SNP.demux.fix.GT.table. This function cleans and standardizes a Genotype assignment table obtained from the SoupOrCell tool.

- #### 111 `Convert10Xfolders()`

  Convert10Xfolders. This function takes a parent directory with a number of subfolders, each containing the standard output of 10X Cell Ranger. It (1) loads the filtered data matrices, (2) converts them to Seurat objects, and (3) saves them as .RDS files.

- #### 112 `ConvertDropSeqfolders()`

  ConvertDropSeqfolders. This function takes a parent directory with a number of subfolders, each containing the standard output of 10X Cell Ranger. It (1) loads the filtered data matrices, (2) converts them to Seurat objects, and (3) saves them as .RDS files.

- #### 113 `LoadAllSeurats()`

  LoadAllSeurats. This function loads all Seurat objects found in a directory. It also works with symbolic links (but not with aliases).

- #### 114 `read10x()`

  read10x. This function reads a 10X dataset from gzipped matrix.mtx, features.tsv and barcodes.tsv files.

- #### 115 `load10Xv3()`

  load10Xv3. Load 10X output folders.

- #### 116 `.saveRDS.compress.in.BG()`

  .saveRDS.compress.in.BG. Save and RDS object and compress resulting file in the background using system(gzip). OS X or unix.

- #### 117 `xread()`

  isave.RDS. Save an RDS object, using a faster and efficient compression method that runs in the background.

- #### 118 `isave.image()`

  isave.image. Save an image of the current workspace using a faster and efficient compression method that runs in the background.

- #### 119 `qsave.image()`

  Save workspace - qsave.image. Faster saving of workspace, and compression outside R, when it can run in the background. Seemingly quite CPU hungry and not very efficient compression. #

- #### 120 `clip10Xcellname()`

  clip10Xcellname. Clip all suffices after underscore (10X adds it per chip-lane, Seurat adds in during integration). #

- #### 121 `make10Xcellname()`

  make10Xcellname. Add a suffix to cell names, so that it mimics the lane-suffix, e.g.: "_1". #

- #### 122 `plotTheSoup()`

  plotTheSoup. Plot stats about the ambient RNA content in a 10X experiment.

- #### 123 `jJaccardIndexVec()`

  jJaccardIndexVec. Calculate jaccard similarity for 2 vecotrs. Helper to jPairwiseJaccardIndexList.

- #### 124 `jPairwiseJaccardIndexList()`

  jPairwiseJaccardIndexList. Create a pairwise jaccard similarity matrix across all combinations of columns in binary.presence.matrix. Modified from: https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/ #

- #### 125 `jPresenceMatrix()`

  jPresenceMatrix. Make a binary presence matrix from a list. Source: https://stackoverflow.com/questions/56155707/r-how-to-create-a-binary-relation-matrix-from-a-list-of-strings #

- #### 126 `jJaccardIndexBinary()`

  jJaccardIndexBinary. Calculate Jaccard Index. Modified from: https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/ #

- #### 127 `jPairwiseJaccardIndex()`

  jPairwiseJaccardIndex. Create a pairwise jaccard similarity matrix across all combinations of columns in binary.presence.matrix. Modified from: https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/ #



------------------------------------------------------------------------------------------------------------
## Seurat.Utils.Metadata.R


- #### 1 `meta_col_exists()`

  Check if a Column Exists in the Metadata of an S4 Object. This function checks whether a given column exists in the meta.data of a Seurat object.

- #### 2 `getMetadataColumn()`

  getMetadataColumn. Retrieves a specified metadata column from a Seurat object and returns it as a named vector.

- #### 3 `get_levels_seu()`

  Get Unique Levels of a Seurat Object Ident Slot. This function extracts the unique levels present in the 'ident' slot of a Seurat object.  The function throws an error if the number of levels exceeds 'max_levels'.  The function optionally prints the R code to recreate the 'Levels' vector using 'dput'. 

- #### 4 `getMedianMetric()`

  getMedianMetric. Get the median values of different columns in meta.data, can iterate over a list of Seurat objects.

- #### 5 `getCellIDs.from.meta()`

  getCellIDs.from.meta. Retrieves cell IDs from a specified metadata column of a Seurat object, where the cell ID matches a provided list of values. The matching operation uses the `%in%` operator.

- #### 6 `seu.add.meta.from.vector()`

  seu.add.meta.from.vector. Adds a new metadata column to a Seurat object.

- #### 7 `create.metadata.vector()`

  Create a Metadata Vector. This function creates a metadata vector from an input vector and a Seurat object.  The resulting vector contains values from 'vec' for the intersecting cell names between 'vec' and 'obj'.  It also checks if the intersection between the cell names in 'vec' and 'obj' is more than a minimum intersection size.

- #### 8 `add.meta.fraction()`

  add.meta.fraction. Add a new metadata column to a Seurat object, representing the fraction of a gene set in the transcriptome (expressed as a percentage).

- #### 9 `add.meta.tags()`

  add.meta.tags. Add metadata tags to a Seurat object dataset.

- #### 10 `seu.add.meta.from.table()`

  seu.add.meta.from.table. Add multiple new metadata columns to a Seurat object from a table. #

- #### 11 `seu.map.and.add.new.ident.to.meta()`

  seu.map.and.add.new.ident.to.meta. Adds a new metadata column to a Seurat object based on an identity mapping table.

- #### 12 `fix.orig.ident()`

  fix.orig.ident. Remove the string "filtered_feature_bc_matrix." from "orig.ident". Helper function.

- #### 13 `seu.RemoveMetadata()`

  seu.RemoveMetadata. Remove specified metadata columns from a Seurat object.

- #### 14 `sampleNpc()`

  Sample N % of a dataframe (obj@metadata), and return rownames (cell IDs).. This function samples a specified percentage of a dataframe (specifically a subset    of the metadata of a Seurat object) and returns the corresponding cell IDs.

- #### 15 `set.all.genes()`

  set.all.genes. It is just a reminder to use calc.q99.Expression.and.set.all.genes to create the all.genes variable

- #### 16 `plotMetadataCorHeatmap()`

  Plot Metadata Correlation Heatmap. This function plots a heatmap of metadata correlation values. It accepts a Seurat object  and a set of metadata columns to correlate. The correlations are calculated using either Pearson  or Spearman methods, and the resulting heatmap can include the principal component (PCA) values  and be saved with a specific suffix. 

- #### 17 `heatmap_calc_clust_median()`

  Calculate and plot heatmap of cluster medians. This function calculates the median of specified variables in a dataframe,  grouped by a column ('ident'). The function also provides an option to scale the medians,  subset the ident levels, and either return a matrix of median values or plot a heatmap. 

- #### 18 `plotMetadataCategPie()`

  plotMetadataMedianFractionBarplot. Generates a barplot of metadata median values.

- #### 19 `transfer_labels_seurat()`

  Transfer Labels in Seurat. Function to transfer labels from a reference Seurat object to a query Seurat object using anchoring and transfer data methods from the Seurat package.  It then visualizes the reference and the combined objects using Uniform Manifold Approximation and Projection (UMAP). 

- #### 20 `match_best_identity()`

  Match and Translate Best Identity. This function matches the best identity from `ident_from` to `ident_to` in an object,  updates the metadata of the object with this new identity and returns the updated object. Additionally,  it generates a UMAP plot based on the new identity. The function replaces original categories with  the most frequent ones, hence helps to filter out less important categories. 




# Usage

You can use most functions at relevant steps of a standard Seurat analysis.


-----------
[Get Seurat.utils](https://github.com/vertesy/Seurat.utils). Vertesy, 2021. [![DOI](https://zenodo.org/badge/248721133.svg)](https://zenodo.org/badge/latestdoi/248721133) 
*If you use these functions, please star the repo, or cite via `DOI`. Thanks!*

<br>

