[![DOI](https://zenodo.org/badge/248721133.svg)](https://zenodo.org/badge/latestdoi/248721133) *If you use these functions, please star the repo, or cite via `DOI`. Thanks!*

# Seurat.utils

`Seurat.utils` Is a collection of utility functions for Seurat v3. Functions allow the automation / multiplexing of plotting, 3D plotting, visualisation of statistics & QC, interaction with the Seurat object, etc.  Some functionalities require functions from [CodeAndRoll2](https://github.com/vertesy/CodeAndRoll2), [ReadWriter](https://github.com/vertesy/ReadWriter), [Stringendo](https://github.com/vertesy/Stringendo), [ggExpressDev](https://github.com/vertesy/ggExpressDev), [MarkdownReports](https://github.com/vertesy/MarkdownReports), and the [Rocinante](https://github.com/vertesy/Rocinante) (See installation).



# Installation

Seurat.utils relies on:

- [Stringendo](https://github.com/vertesy/Stringendo)
- [ReadWriter](https://github.com/vertesy/ReadWriter)
- [CodeAndRoll2](https://github.com/vertesy/CodeAndRoll2)
- [MarkdownHelpers](https://github.com/vertesy/MarkdownHelpers)
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




# Usage

You can use most functions at relevant steps of a standard Seurat analysis.

# Function descriptions

Updated: `01 11 2021`

- #### SmallestNonAboveX
replace small values with the next smallest value found, which is >X. #


- #### Add.DE.combined.score
Add combined score to DE results. (LFC * -log10( p_cutoff / pval_scaling ) )


- #### StoreTop25Markers
Save the top 25 makers based on `avg_log2FC` output table of `FindAllMarkers()` (df_markers) under `@misc$df.markers$res...`. By default, it rounds up insignificant digits up to 3. #


- #### StoreAllMarkers
Save the output table of `FindAllMarkers()` (df_markers) under `@misc$df.markers$res...`. By default, it rounds up insignificant digits up to 3. #


- #### GetTopMarkersDF
Get the vector of N most diff. exp. genes. #


- #### GetTopMarkers
Get the vector of N most diff. exp. genes. #


- #### AutoLabelTop.logFC
Create a new "named identity" column in the metadata of a Seurat object, with `Ident` set to a clustering output matching the `res` parameter of the function. It requires the output table of `FindAllMarkers()`. If you used `StoreAllMarkers()` is stored under `@misc$df.markers$res...`, which location is assumed by default. #


- #### AutoLabel.KnownMarkers
Create a new "named identity" column in the metadata of a Seurat object, with `Ident` set to a clustering output matching the `res` parameter of the function. It requires the output table of `FindAllMarkers()`. If you used `StoreAllMarkers()` is stored under `@misc$df.markers$res...`, which location is assumed by default. #


- #### DimPlot.ClusterNames
Plot UMAP with Cluster names. #


- #### AutoNumber.by.UMAP
Relabel cluster numbers along a UMAP (or tSNE) axis #


- #### AutoNumber.by.PrinCurve
Relabel cluster numbers along the principal curve of 2 UMAP (or tSNE) dimensions. #


- #### ggplotColours
Generate ggplot colours for slingshot


- #### points_on_curve
Helper to visualize points_on_curve in slingshot.


- #### points_on_curve.principal_curve
Helper to visualize points_on_curve in slingshot.


- #### points_on_curve.SlingshotDataSet
Helper to visualize points_on_curve in slingshot.


- #### gg_plot
Adjusted gg_plot for slingshot


- #### jJaccardIndexVec
Calculate jaccard similarity for 2 vecotrs. Helper to jPairwiseJaccardIndexList.


- #### jPairwiseJaccardIndexList
Create a pairwise jaccard similarity matrix across all combinations of columns in binary.presence.matrix. Modified from: https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/ #


- #### jPresenceMatrix
Make a binary presence matrix from a list. Source: https://stackoverflow.com/questions/56155707/r-how-to-create-a-binary-relation-matrix-from-a-list-of-strings #


- #### jJaccardIndexBinary
Calculate Jaccard Index. Modified from: https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/ #


- #### jPairwiseJaccardIndex
Create a pairwise jaccard similarity matrix across all combinations of columns in binary.presence.matrix. Modified from: https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/ #


- #### getMedianMetric
Get the median values of different columns in meta.data, can iterate over a list of Seurat objects.


- #### add.meta.tags
N is the for which dataset #


- #### add.meta.fraction
Add a new metadata column, with the fraction of gene set in the transcripome (percentage).


- #### GetClusteringRuns
Get Clustering Runs: metadata column names #


- #### GetNamedClusteringRuns
Get Clustering Runs: metadata column names #


- #### GetOrderedClusteringRuns
Get Clustering Runs: metadata column names #


- #### GetNumberOfClusters
Get Number Of Clusters #


- #### getMetadataColumn <- mmeta
Get a metadata column from a Seurat object as a named vector #


- #### getCellIDs.from.meta
Get cellIDs from a metadata column, matching a list of values (using %in%). #


- #### seu.add.meta.from.vector
Add a new metadata column to a Seurat  object


- #### seu.map.and.add.new.ident.to.meta
Add a new metadata column to a Seurat  object


- #### calc.cluster.averages
Calculate the average of a metadata column (numeric) per cluster.


- #### seu.add.meta.from.table
Add multiple new metadata columns to a Seurat object from a table. #


- #### sampleNpc
Sample N % of a dataframe (obj@metadata), and return the cell IDs. #


- #### calc.q90.Expression.and.set.all.genes
Calculate the gene expression of the e.g.: 90th quantile (expression in the top 10% cells). #


- #### PlotTopGenes
Plot the highest expressed genes on umaps, in a subfolder. Requires calling calc.q90.Expression.and.set.all.genes before. #


- #### fix.orig.ident
Remove the string "filtered_feature_bc_matrix." from "orig.ident". Helper function.


- #### set.all.genes
It is just a reminder to use calc.q90.Expression.and.set.all.genes to create the all.genes variable


- #### recall.all.genes
all.genes set by calc.q90.Expression.and.set.all.genes() #


- #### recall.meta.tags.n.datasets
Recall  meta.tags from obj@misc to "meta.tags" in the global environment.


- #### recall.parameters
Recall parameters from obj@misc to "p" in the global environment.


- #### recall.genes.ls
Recall genes.ls from obj@misc to "genes.ls" in the global environment.


- #### save.parameters
Save parameters to obj@misc$p


- #### plot.expression.rank.q90
Plot gene expression based on the expression at the 90th quantile (so you will not lose genes expressed in few cells).


- #### FlipReductionCoordinates
Flip reduction coordinates (like UMAP upside down).


- #### SeuratColorVector
Recall a Seurat color vector.


- #### getClusterColors
get Seurat's cluster colors.


- #### get.clustercomposition
Get cluster composition: which datasets contribute to each cluster?


- #### mplotGene
Plot genes in Monocle.


- #### mplotManyGenes
Plot many genes in Monocle.


- #### m3DplotGene
Plot a gene in 3D in Monocle.


- #### m3DplotKeyGenes
Plot many genes in 3D in Monocle.


- #### subsetMonocleObject
Subset a compressed Seurat Obj and save it in wd.


- #### m3.get.umap
Fetch the umap coordinates from obj@int_colData@listData$reducedDims[[slot]]


- #### m3.backup.umap
Backup umap coordinates to obj@int_colData@listData$reducedDims[[new.slot]]


- #### m3.recall.umap
Fetch UMAP coordinates.


- #### m3.export.umap.2.Seurat
Export umap coordinates.


- #### BarTableSweepList
BarTableSweepList


- #### mSeq.map.all96.BCs
mSeq.map.all96.BCs


- #### aux.plotAllMseqBCs
aux.plotAllMseqBCs


- #### qUMAP
The quickest way to draw a gene expression UMAP #


- #### clUMAP
The quickest way to draw a clustering result  UMAP #


- #### gg_color_hue
reproduce the ggplot2 default color palette #


- #### save2umaps.A4
Save 2 umaps on 1 A4


- #### save4umaps.A4
Save 4 umaps on 1 A4


- #### umapNamedClusters
Plot and save umap based on a metadata column. #


- #### qqSaveGridA4
Save 2 or 4 ggplot objects using plot_grid() on an A4 page #


- #### umapHiLightSel
Highlight a set of cells based on clusterIDs provided. #


- #### multiFeaturePlot.A4
Save multiple FeaturePlots, as jpeg, on A4 for each gene, which are stored as a list of gene names. #


- #### multiFeatureHeatmap.A4
Save multiple FeatureHeatmaps from a list of genes on A4 jpeg #


- #### plot.UMAP.tSNE.sidebyside
Plot a UMAP and tSNE sidebyside #


- #### PlotTopGenesPerCluster
Plot the top N diff. exp. genes in each cluster


- #### qFeatureScatter
Quickly plot and save a FeatureScatter plot.


- #### qMarkerCheck.BrainOrg
Quickly plot key markers in brain organoids


- #### getDiscretePalette
Generate a Discrete color Palette.


- #### ww.check.if.3D.reduction.exist
ww.check.if.3D.reduction.exist in backup slot #


- #### ww.check.quantile.cutoff.and.clip.outliers
Helper function.


- #### plot3D.umap.gene
Plot a 3D umap with gene expression. Uses plotly. Based on github.com/Dragonmasterx87. #


- #### plot3D.umap
Plot a 3D umap based on one of the metadata columns. Uses plotly. Based on github.com/Dragonmasterx87. #


- #### SavePlotlyAsHtml
Save Plotly 3D scatterplot as an html file. #


- #### BackupReduction
Backup UMAP to `obj@misc$reductions.backup` from `obj@reductions$umap`. #


- #### SetupReductionsNtoKdimensions
Calculate N-to-K dimensional umaps (default = 2:3); and back them up UMAP to `obj@misc$reductions.backup` from @reductions$umap #


- #### RecallReduction
Set active UMAP to `obj@reductions$umap` from `obj@misc$reductions.backup`. #


- #### Annotate4Plotly3D
Create annotation labels for 3D plots. Source https://plot.ly/r/text-and-annotations/#3d-annotations #


- #### Plot3D.ListOfGenes
Plot and save list of 3D UMAP ot tSNE plots using plotly. #


- #### Plot3D.ListOfCategories
Plot and save list of 3D UMAP ot tSNE plots using plotly. #


- #### PlotFilters
Plot filtering threshold and distributions, using four panels to highlight the relation between Gene- and UMI-count, ribosomal- and mitochondrial-content. #


- #### seu.PC.var.explained
Determine percent of variation associated with each PC. #


- #### seu.plot.PC.var.explained
Plot the percent of variation associated with each PC. #


- #### BarplotCellsPerObject
Take a List of Seurat objects and draw a barplot for the number of cells per object. #


- #### CellFractionsBarplot2
Barplot the Fraction of cells per cluster.


- #### barplot.cells.per.cluster
Barplot the Fraction of cells per cluster. (dupl?)


- #### BulkGEScatterPlot
Plot bulk scatterplots to identify differential expressed genes across conditions #


- #### sparse.cor
Sparse, fast correlation.


- #### Calc.Cor.Seurat
Calculate gene correlation on a Seurat object.


- #### plot.Metadata.Cor.Heatmap
Plot a heatmap with Metadata correlation values.


- #### plot.Metadata.median.fraction.barplot
Barplot Metadata median values


- #### plot.Gene.Cor.Heatmap
Plot a gene correlation heatmap.


- #### plot.clust.size.distr
Barplot of Histogram of cluster size distribution


- #### gene.expression.level.plots
Histogram of gene expression levels.


- #### PrctCellExpringGene
From Github/Ryan-Zhu https://github.com/satijalab/seurat/issues/371 #


- #### ww.calc_helper
From Github/Ryan-Zhu https://github.com/satijalab/seurat/issues/371 #


- #### scBarplotFractionAboveThr
Barplot the fraction of cell above a threshold value (based on a meta.data column), in each cluster.


- #### scBarplotFractionBelowThr
Barplot the fraction of cell below a threshold value (based on a meta.data column), in each cluster.


- #### Convert10Xfolders
Take a parent directory with a number of subfolders, each containing the standard output of 10X Cell Ranger. (1.) It loads the filtered data matrices; (2.) converts them to Seurat objects, and (3.) saves them as *.RDS files. #


- #### Convert10Xfolders.old
Take a parent directory with a number of subfolders, each containing the standard output of 10X Cell Ranger. (1.) It loads the filtered data matrices; (2.) converts them to Seurat objects, and (3.) saves them as *.RDS files. #


- #### ConvertDropSeqfolders
Take a parent directory with a number of subfolders, each containing the standard output of 10X Cell Ranger. (1.) It loads the filtered data matrices; (2.) converts them to Seurat objects, and (3.) saves them as *.RDS files. #


- #### LoadAllSeurats
Load all Seurat objects found in a directory. Also works with symbolic links (but not with aliases). #


- #### read10x
read10x from gzipped matrix.mtx, features.tsv and barcodes.tsv #


- #### saveRDS.compress.in.BG
Save and RDS object and compress resulting file in the background using system(gzip). OS X or unix.


- #### isave.RDS
Save and RDS object.


- #### isave.image
Save and RData image.


- #### qsave.image
Faster saving of workspace, and compression outside R, when it can run in the background. Seemingly quite CPU hungry and not very efficient compression. #


- #### subsetSeuObj
Subset a compressed Seurat Obj and save it in wd. #


- #### subsetSeuObj.and.Save
Subset a compressed Seurat Obj and save it in wd. #


- #### Downsample.Seurat.Objects
Downsample a list of Seurat objects


- #### clip10Xcellname
Clip all suffices after underscore (10X adds it per chip-lane, Seurat adds in during integration). #


- #### make10Xcellname
Add a suffix to cell names, so that it mimics the lane-suffix, e.g.: "_1". #


- #### seu.Make.Cl.Label.per.cell
Take a named vector (of e.g. values ="gene names", names = clusterID), and a vector of cell-IDs and make a vector of "GeneName.ClusterID". #


- #### GetMostVarGenes
Get the most variable rGenes #


- #### gene.name.check
Check gene names in a seurat object, for naming conventions (e.g.: mitochondrial reads have - or .). Use for reading .mtx & writing .rds files. #


- #### check.genes
Check if a gene name exists in a Seurat object, or in HGNC?


- #### fixZeroIndexing.seurat
Fix zero indexing in seurat clustering, to 1-based indexing #


- #### CalculateFractionInTrome
Calculate the fraction of a set of genes within the full Transcriptome of each cell. #


- #### AddNewAnnotation
Create a new metadata column based on an exisiting metadata column and a list of mappings (name <- IDs). #


- #### whitelist.subset.ls.Seurat
Subset cells in a (list of) Seurat objects, based on an externally provided list of cell IDs.


- #### FindCorrelatedGenes
Find correlated genes in a Seurat object


- #### UpdateGenesSeurat
Update genes symbols that are stored in a Seurat object. It returns a data frame. The last column are the updated gene names. #


- #### RenameGenesSeurat
Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data. #


- #### RemoveGenesSeurat
Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes metadata; obj@assays$RNA@counts, @data and @scale.data. #


- #### HGNC.EnforceUnique
Enforce Unique names after HGNC symbol update. updatedSymbols is the output of HGNChelper::checkGeneSymbols. #


- #### GetUpdateStats
Plot the Symbol-update statistics. Works on the data frame returned by `UpdateGenesSeurat()`. #


- #### PlotUpdateStats
Scatter plot of update stats. #


- #### plotTheSoup
Plot stats about the ambient RNA content in a 10X experiment.


- #### load10Xv3
Load 10X output folders.



CURRENTLY NOT THE CASE:

### The content of this repo is organized into files per rough functionalities:

1. `Seurat.Utils.Load.R`								→ *Top level wrapper to source each file below.*
2. `metadata.manipulation.R`						→ *Metadata manipulation.*
3. `plotting.dim.reduction.2D.R`				→ *Plotting dimensionality reduction in 2D (UMAP, tSNE, PCA).*
4. `plotting.dim.reduction.3D.R`				→ *Plotting dimensionality reduction in 3D (UMAP, tSNE, PCA).*
5. `plotting.filtering.R`								→ *Plotting filtering.*
6. `plotting.statistics.and.QC.R`				→ *Plotting statistics and QC.*
7. `Read.Write.Save.Load.functions.R`		→ *Read Write Save Load functions.*
8. `Seurat.object.manipulations.etc.R`	→ *Seurat object manipulations etc.*
9. `Seurat.update.gene.symbols.HGNC.R`	→ *Seurat update gene symbols HGNC.*



-----------
[Get Seurat.utils](https://github.com/vertesy/Seurat.utils). Vertesy, 2021. [![DOI](https://zenodo.org/badge/248721133.svg)](https://zenodo.org/badge/latestdoi/248721133) 
*If you use these functions, please star the repo, or cite via `DOI`. Thanks!*

<br>

