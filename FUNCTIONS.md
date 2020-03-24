## Read.Write.Save.Load.functions.R
Updated: `Tue Mar 24 10:01 2020`


- #### `Convert10Xfolders()`
  Take a parent directory with a number of subfolders, each containing the standard output of 10X Cell Ranger. (1.) It loads the filtered data matrices; (2.) converts them to Seurat objects, and (3.) saves them as *.RDS files.

- #### `LoadAllSeurats()`
  Load all Seurat objects found in a directory. Also works with symbolic links (but not with aliases).

- #### `read10x()`
  read10x from gzipped matrix.mtx, features.tsv and barcodes.tsv

- #### `isave.RDS()`
  Faster saving of workspace, and compression outside R, when it can run in the background. Seemingly quite CPU hungry and not very efficient compression.

- #### `isave.RDS.pigz()`
  Faster saving of workspace, and compression outside R, when it can run in the background. Seemingly quite CPU hungry and not very efficient compression.

- #### `isave.image()`
  Faster saving of workspace, and compression outside R, when it can run in the background. Seemingly quite CPU hungry and not very efficient compression.

- #### `subsetSeuObj.and.Save()`
Subset a compressed Seurat Obj and save it in wd.

#### *The functions below are now part of the [Seurat.multicore](https://github.com/vertesy/Seurat.multicore) library:*


- #### `seuSaveRds()`
  Save a compressed Seurat Object, with parallel gzip by pgzip

- #### `rrRDS()`
  Load a list of RDS files with parallel ungzip by pgzip.

- #### `sssRDS()`
  Save multiple objects into a list of RDS files using parallel gzip by pgzip (optional).

- #### `ssaveRDS()`
  Save an object with parallel gzip by pgzip.

- #### `rreadRDS()`
  Read an object with parallel ungzip by pgzip.

- #### `snappy_pipe()`
  Alternative, fast compression. Low compression rate, lightning fast.

- #### `pigz_pipe()`
Alternative
normal gzip output (& compression rate), ~*cores faster in zipping.

## Seurat.object.manipulations.etc.R

- #### `clip10Xcellname()`
Clip all suffices after underscore (10X adds it per chip-lane, Seurat adds in during integration).
- #### `make10Xcellname()`
Add a suffix to cell names, so that it mimics the lane-suffix, e.g.
"_1".
- #### `seu.Make.Cl.Label.per.cell()`
Take a named vector (of e.g. values ="gene names", names = clusterID), and a vector of cell-IDs and make a vector of "GeneName.ClusterID".
- #### `GetMostVarGenes()`
Get the most variable rGenes
- #### `gene.name.check()`
Check gene names in a seurat object, for naming conventions (e.g.
mitochondrial reads have - or .). Use for reading .mtx & writing .rds files.
- #### `check.genes()`
Check if genes exist in your dataset.
- #### `fixZeroIndexing.seurat()`
Fix zero indexing in seurat clustering, to 1-based indexing

## Seurat.update.gene.symbols.HGNC.R

- #### `RenameGenesSeurat()`
Replace gene names in different slots of a Seurat object. Run this before integration. It only changes SeuObj@assays$RNA@counts, @data and @scale.data.
- #### `UpdateGenesSeurat()`
Update genes symbols that are stored in a Seurat object. It returns a data frame. The last column are the updated gene names.
- #### `plot.UpdateStats()`
Plot the Symbol-update statistics. Works on the data frame returned by `UpdateGenesSeurat()`.

## metadata.manipulation.R

- #### `getMetadataColumn <- mmeta()`
  Get a metadata column from a Seurat object as a named vector

- #### `getCellIDs.from.meta()`
  Get cellIDs from a metadata column, matching a list of values (using %in%).

- #### `seu.add.meta.from.vector()`
  Add a new metadata column to a Seurat  object

- #### `seu.add.meta.from.table()`
Add multiple new metadata columns to a Seurat object from a table.

- #### `sampleNpc()`

  Sample N % of a dataframe (obj@metadata), and return the cell IDs.

## plotting.dim.reduction.2D.R

- #### `qUMAP()`
  The quickest way to a draw a UMAP

- #### `gg_color_hue()`
  reproduce the ggplot2 default color palette

- #### `save2umaps.A4()`
  Save 2 umaps on an A4 page.

- #### `save4umaps.A4()`
  Save 4 umaps on an A4 page.

- #### `umapNamedClusters()`
  Plot and save umap based on a metadata column.

- #### `umapHiLightSel()`
  Highlight a set of cells based on clusterIDs provided.

- #### `multiFeaturePlot.A4()`
  Save multiple FeaturePlots, as jpeg, on A4 for each gene, which are stored as a list of gene names.

- #### `multiFeatureHeatmap.A4()`
  Save multiple FeatureHeatmaps from a list of genes on A4 jpeg

- #### `plot.UMAP.tSNE.sidebyside()`
Plot a UMAP and tSNE sidebyside



## plotting.dim.reduction.3D.R

 - #### `plot3D.umap.gene` 
 Plot a 3D umap with gene expression. Uses plotly. Based on [Dragonmasterx87](https://github.com/Dragonmasterx87/Interactive-3D-Plotting-in-Seurat-3.0.0).
 
 - #### `plot3D.umap` 
    Plot a 3D umap based on one of the metadata columns. Uses plotly. Based on [Dragonmasterx87](https://github.com/Dragonmasterx87/Interactive-3D-Plotting-in-Seurat-3.0.0).
    
 - #### `SavePlotlyAsHtml` 
 Save Plotly 3D scatterplot as an html file.
 
 - #### `BackupReduction` 
 Backup UMAP to `obj@misc$reductions.backup` from `obj@reductions$umap`.
 
 - #### `SetupReductionsNtoKdimensions` 
 Calculate N-to-K dimensional umaps (default 
 
 - #### `RecallReduction` 
 Set active UMAP to `obj@reductions$umap` from `obj@misc$reductions.backup`.
 
 - #### `Annotate4Plotly3D` 
 Create annotation labels for 3D plots. Source https://plot.ly/r/text-and-annotations/#3d-annotations
 
 

## plotting.filtering.R

- #### `PlotFilters()`
Plot filtering threshold and distributions, using four panels to highlight the relation between Gene- and UMI-count, ribosomal- and mitochondrial-content.

## plotting.statistics.and.QC.R

- #### `seu.PC.var.explained()`
Determine percent of variation associated with each PC.
- #### `seu.plot.PC.var.explained()`
Plot the percent of variation associated with each PC.
- #### `BarplotCellsPerObject()`
Take a List of Seurat objects and draw a barplot for the number of cells per object.
- #### `sgCellFractionsBarplot.Mseq()`
Cell fractions Barplot for MULTI-seq. sg stands for "seurat ggplot".
- #### `ssgCellFractionsBarplot.CORE()`
Cell Fractions Barplots, basic. sg stands for "seurat ggplot".
- #### `sgCellFractionsBarplot()`
Cell Fractions Barplots. sg stands for "seurat ggplot".
- #### `plotTheSoup()`
Plot the ambient RNA content of droplets without a cell (background droplets).
