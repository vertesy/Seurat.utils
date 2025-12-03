## List of Functions in Seurat.Utils.Visualization.R (61) 
Updated: 2025/12/03 10:19
- #### 1 `PlotFilters()`
Plot filtering thresholds and distributions. This function plots the filtering thresholds and distributions for Seurat objects,  using four panels to highlight the relationship between gene- and UMI-counts, and the  ribosomal- and mitochondrial-content.  !! Default arguments assume that `p` is a list of  parameters, present in the global environment, with elements `thr.lp.mito`, `thr.hp.mito`,  `thr.lp.ribo`, `thr.hp.ribo`, `thr.lp.nFeature_RNA`, and `thr.hp.nFeature_RNA`. 

- #### 2 `scCalcPCAVarExplained()`
Calculate the percent of variation explained by individual PC's. This function calculates the percentage of variation each principal component (PC)  accounts for in a Seurat object. It's specifically tailored for Seurat objects and provides a  convenient way to understand the variance distribution across PCs. For similar calculations on  standard PCA objects, refer to github.com/vertesy/Rocinante `PCA.percent.var.explained()`. 

- #### 3 `scPlotPCAvarExplained()`
Plot the percent of variation explained by individual PC's. This function plots the percentage of variation explained by each principal  component (PC) in a Seurat object. It allows for a visual assessment of how much variance is  captured by each PC, which is crucial for dimensionality reduction analyses. Users can choose  between two plotting methods: one using `MarkdownReports` and the other using `ggExpress`. 

- #### 4 `PercentInTranscriptome()`
Gene Expression as Fraction of Total UMI Counts. This function computes and visualizes gene expression levels as a fraction of total  UMI (Unique Molecular Identifier) counts across all genes in a Seurat object. It aims to highlight  the relative contribution of the most highly expressed genes to the overall transcriptome. 

- #### 5 `plotGeneExpressionInBackgroundHist()`
Histogram All Genes' Expression Level and a Highlighted Gene. Shows a comparison of the expression level of the chosen gene to all genes.  Very useful to see if the gene has a meaningful expression level. This function generates a  histogram to visualize the expression level distribution of a specified gene across all cells in  a Seurat object. It highlights the position of the gene of interest within the overall distribution. 

- #### 6 `plotGeneExprHistAcrossCells()`
Histogram of Gene / Geneset Aggregate Expression Across Cells. Creates and optionally saves a histogram showing expression levels of specified genes  within a Seurat object. Provides options for aggregate gene expression, expression threshold filtering,  and quantile clipping for count data. 

- #### 7 `PctCellsAboveX()`
Compute % of Cells Above a Threshold for a Metadata or Gene Feature.   Computes the fraction of cells with values above a specified threshold for  either a metadata column or an assay feature. Supports optional subsetting  and optional regrouping into boxplot-style categories. 

- #### 8 `  calc_proportions()`
PctCellsExpressingGenes. Calculates the proportion of cells expressing one or more specified genes using a Seurat  object as input. 

- #### 9 `scBarplot.CellFractions()`
Generate Barplot of Cell Fractions. This function generates a bar plot of cell fractions per cluster from a Seurat object.  It offers the option to downsample data, equalizing the number of cells in each group  to the number in the smallest group. The plot's bars are grouped by one variable and filled by another.  The function supports custom color palettes, drawing numerical values on bars, and saving the plot. 

- #### 10 `scBarplot.CellsPerCluster()`
Barplot of Fraction of Cells per Cluster. Visualizes the fraction of cells within each cluster through a barplot. 

- #### 11 `scBarplot.FractionAboveThr()`
Barplot the Fraction of Cells Above Threshold per Cluster. Generates a bar plot depicting the percentage of cells within each cluster that  exceed a specified threshold, based on a selected metadata column. 

- #### 12 `scBarplot.FractionBelowThr()`
Fraction of Cells Below Threshold per Cluster. Generates a bar plot to visualize the percentage of cells within each cluster that  fall below a specified threshold, according to a metadata column value.  Inherits all parameters from `scBarplot.FractionAboveThr` with the exception that `above` is set to FALSE. 

- #### 13 `scPieClusterDistribution()`
scPieClusterDistribution. This function generates a pie chart of cluster distributions for a given clustering  identity in a single-cell RNA-seq object. 

- #### 14 `scBarplot.CellsPerObject()`
Barplot of Cells Per Seurat Object. Visualizes the number of cells in each Seurat object within a list, showing the  distribution of cell counts across different datasets or experimental conditions. 

- #### 15 `scBarplotStackedMetaCateg_List()`
Stacked Barplot of Metadata Categories for List of Seurat Objects. Creates and saves a stacked barplot for a specified metadata category  from a list of Seurat objects. 

- #### 16 `gg_color_hue()`
Reproduce the ggplot2 default color palette. Generates a vector of colors that emulates the default color palette used by ggplot2.  This function is useful for creating color sets for custom plotting functions or for applications  outside of ggplot2 where a similar aesthetic is desired. 

- #### 17 `getDiscretePalette()`
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
Plot and Save Heatmaps from Metadata Calculation Results. Generates and saves heatmap visualizations for each metric in the results obtained  from metadata calculations, such as  `calculateAverageMetaData() - mean or median values of  specified features across different categories. 

- #### 23 `qFeatureScatter()`
Scatter Plot of Two Features in Seurat Object. Generates a scatter plot comparing two features (genes or metrics) from a Seurat  object and optionally saves it. The function wraps around Seurat's `FeatureScatter` for  enhanced usability, including optional logarithmic transformations and saving capabilities. 

- #### 24 `qSeuViolin()`
Create a Violin Plot for a Seurat Object Feature and save the file.. Generates a violin plot for a specified feature in a Seurat object,  allowing for the data to be split by a specified grouping variable.  The function supports customization options such as logarithmic scaling, custom titles, and more. 

- #### 25 `qUMAP()`
Quick UMAP Visualization of Gene Expression and automatically save the plot. Generates a UMAP visualization for a specific feature from a Seurat object, and  automatically saves it. Offers options for custom titles, subtitles, saving, and more. Assumes  default options for custom titles, subtitles, saving, and more. 

- #### 26 `clUMAP()`
clUMAP - Quick Visualization of Clustering Results with UMAP and automatically save the plot. Generates a UMAP visualization based on clustering results from a Seurat object,  and automatically saves it. Offers options for custom titles, subtitles, saving, and more. Assumes  default options for custom titles, subtitles, saving, and more. 

- #### 27 `umapHiLightSel()`
Highlight Selected Clusters on UMAP. Generates a UMAP plot from a Seurat object with specified clusters highlighted.  It saves the resulting UMAP plot directly to the current working directory. 

- #### 28 `DimPlot.ClusterNames()`
DimPlot.ClusterNames. Plot UMAP with Cluster names.

- #### 29 `multiFeaturePlot.A4()`
multiFeaturePlot.A4. Save multiple FeaturePlots, as jpeg, on A4 for each gene, which are stored as a list of gene names.

- #### 30 `multiSingleClusterHighlightPlots.A4()`
Generate Cluster Highlight UMAPs compiled into A4 pages. This function generates and saves cluster highlight plots for both single and multiple  clusters using UMAP or other dimensionality reduction techniques. It supports saving plots in various  formats and allows customization of plot appearance and layout. 

- #### 31 `qClusteringUMAPS()`
Quick Clustering UMAPs on A4 Page. Generates and arranges UMAP plots for up to four specified clustering resolutions  from a Seurat object onto an A4 page, facilitating comparative visualization. 

- #### 32 `qGeneExpressionUMAPS()`
Quickly Draw 4 Gene Expression UMAPs on an A4 Page. Generates and arranges UMAP plots for up to four specified gene expressions  from a Seurat object onto an A4 page, facilitating comparative visualization. 

- #### 33 `plotQUMAPsInAFolder()`
Plot qUMAPs for Genes in a Folder. This function plots qUMAPs for a specified set of genes, storing the results in a  specified folder. If no folder name is provided, it defaults to using the gene set name. 

- #### 34 `PlotTopGenesPerCluster()`
Plot Top N Differentially Expressed Genes Per Cluster. Visualizes the top N differentially expressed (DE) genes for each cluster within a  specified clustering resolution of a Seurat object, facilitating the exploration of gene  expression patterns across clusters. 

- #### 35 `qQC.plots.BrainOrg()`
Quickly Plot Key QC Markers in Brain Organoids. Generates and arranges UMAP plots for specified QC features  from a Seurat object on an A4 page, facilitating a quick quality control (QC) overview. 

- #### 36 `qMarkerCheck.BrainOrg()`
Quickly Plot Key Markers in Brain Organoids. Generates plots for a predefined or custom set of gene markers within brain organoids,  aiding in the quick assessment of their expression across different cells or clusters. 

- #### 37 `PlotTopGenes()`
Plot Top Genes. This function plots the highest expressed genes on UMAPs, saving the plots in a  subfolder. It requires the prior execution of `calc.q99.Expression.and.set.all.genes`. 

- #### 38 `FlipReductionCoordinates()`
Flip Reduction Coordinates. Flips dimensionality reduction coordinates (such as UMAP or tSNE) vertically or  horizontally to change the visualization perspective. 

- #### 39 `AutoNumber.by.UMAP()`
Relabel Cluster Numbers Along a UMAP (or tSNE) Axis. Automatically renumbers clusters based on their position along a specified dimension  in a UMAP (or tSNE or PCA) plot, potentially enhancing interpretability by ordering clusters. 

- #### 40 `scEnhancedVolcano()`
scEnhancedVolcano. This function creates an enhanced volcano plot. 

- #### 41 `.estMinimumFC()`
Estimate Minimum Log2-Based Fold Change. This function estimates the minimum log2-based fold change from a data frame column. 

- #### 42 `countRelevantEnrichments()`
Count Relevant Enrichments. This function counts the number of relevantly expressed genes from a differential  gene expression table. It considers genes to be relevant if they fall under a maximum p-value  cutoff and are above a minimum log2 fold change cutoff. The function reports the number of  enriched and depleted genes. 

- #### 43 `scGOEnrichment()`
Perform GO Enrichment Analysis. This function performs Gene Ontology (GO) enrichment analysis using the  `clusterProfiler::enrichGO` function. It takes the gene list, universe, organism database,  gene identifier type, and ontology type as inputs and returns the enrichment results. 

- #### 44 `scBarplotEnrichr()`
Barplot GO Enrichment Results by enrichplot. This function creates a bar plot of GO enrichment analysis results using the  `enrichplot::barplot.enrichResult` function. It also allows saving the plot to a file. 

- #### 45 `filterGoEnrichment()`
Filter GO Enrichment Results. This function filters GO enrichment results based on adjusted p-value and q-value  cutoffs, and retrieves the descriptions of the filtered results. 

- #### 46 `countEnrichedDepletedGenes()`
Count Enriched and Depleted Genes. This function counts the number of significantly enriched and depleted genes  based on the provided criteria. It filters the genes based on adjusted p-value and  logarithm of fold change. 

- #### 47 `.adjustLayout()`
Adjust Layout Parameters for multi* plotting functions. Adjusts layout dimensions and properties based on the specified layout type.               Updates the provided environment with new dimensions and layout configuration. 

- #### 48 `save2plots.A4()`
Save Two Plots on One A4 Page. Arranges and saves two UMAP plots (or any plots) side-by-side or one above  the other on a single A4 page. 

- #### 49 `save4plots.A4()`
Save Four Plots on One A4 Page. Arranges and saves four plots (e.g. UMAPs) onto a single A4 page, allowing for a  compact comparison of different visualizations or clustering results. 

- #### 50 `qqSaveGridA4()`
qqSaveGridA4. Saves a grid of 2 or 4 ggplot objects onto an A4 page.

- #### 51 `ww.check.quantile.cutoff.and.clip.outliers()`
Check Quantile Cutoff and Clip Outliers. Checks a specified quantile cutoff and clips outliers from an expression vector,  ensuring that a minimum number of cells expressing a gene remain. 

- #### 52 `plot3D.umap.gene()`
plot3D.umap.gene. Plot a 3D umap with gene expression. Uses plotly. Based on github.com/Dragonmasterx87.

- #### 53 `plot3D.umap()`
plot3D.umap. Plot a 3D umap based on one of the metadata columns. Uses plotly. Based on github.com/Dragonmasterx87.

- #### 54 `SavePlotlyAsHtml()`
SavePlotlyAsHtml. Save a Plotly 3D scatterplot as an HTML file.

- #### 55 `BackupReduction()`
Backup Dimensionality Reduction Data. This function is mostly used internally.It stores a backup of specified  dimensionality reduction data (e.g., UMAP, tSNE, PCA)  within the Seurat object, from `obj@reductions$umap` to the `@misc$reductions.backup` slot. This  allows to store 2D and 3D UMAP visualizations in parallel and easily switch between them via  the `RecallReduction` function. 

- #### 56 `SetupReductionsNtoKdimensions()`
SetupReductionsNtoKdimensions. Function to calculate N-to-K dimensional umaps (default = 2:3); and back them up to  slots `obj@misc$reductions.backup` from @reductions$umap

- #### 57 `RecallReduction()`
Recall Dimensionality Reduction from backup slot. Restores dimensionality reduction data (e.g., UMAP, tSNE, PCA) from a backup  stored within `obj@misc$reductions.backup` to the active `obj@reductions` slot. 

- #### 58 `.Annotate4Plotly3D()`
.Annotate4Plotly3D. Internal helper function. Create annotation labels for 3D plots.  Source https://plot.ly/r/text-and-annotations/#3d-annotations.

- #### 59 `Plot3D.ListOfGenes()`
Plot3D.ListOfGenes. Plot and save list of 3D UMAP or tSNE plots using plotly.

- #### 60 `Plot3D.ListOfCategories()`
Plot3D.ListOfCategories. This function plots and saves a list of 3D UMAP or tSNE plots using plotly.

- #### 61 `panelCorPearson()`
Display Correlation Values in Pairs Plot. This function displays the correlation coefficient and significance level within  a scatterplot generated by the `pairs()` function. The default correlation method is Pearson,  but Kendall or Spearman methods can also be selected. 

