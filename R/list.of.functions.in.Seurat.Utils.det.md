## List of Functions in Seurat.Utils.R (98) 
Updated: 2025/12/03 10:19
- #### 1 `processSeuratObject()`
Process Seurat Objects in Parallel. Applies a series of Seurat processing steps to each Seurat object in a list.               The operations include scaling data, running PCA, UMAP, finding neighbors, and finding clusters.               This is done in parallel using multiple cores. 

- #### 2 `runDGEA()`
Run Differential Gene Expression Analysis (DGEA). Runs a differential gene expression analysis based on specified parameters,  reorders clusters if needed, and optionally saves results. Supports output and plotting configurations. 

- #### 3 `UpdateSeuratObjectProperly()`
Update Seurat Object Properly, including Assays and DimReducs. This function is an extension on `SeuratObject::UpdateSeuratObject()`. It  first calls `UpdateSeuratObject()`, to updates the class definitions of of a (v3) Seurat object,  then it updates its assays to the 'Assay5' class, and updates the UMAP DimReduc to keys. 

- #### 4 `parallel.computing.by.future()`
parallel.computing.by.future. Run gc(), load multi-session computing and extend memory limits.

- #### 5 `IntersectGeneLsWithObject()`
Intersect Genes with Seurat Object. Intersects a set of gene names with those found in a Seurat object.

- #### 6 `SelectHighlyExpressedGenesq99()`
Intersect Genes with the List of Noticeably Expressed Genes. Intersects a vector of gene names with a Seurat object to find genes that are both  in the input list and have expression levels in the top quantiles as defined by the object's  q99 expression data. It aims to filter genes based on their expression levels being above a  specified threshold. Additionally, it offers an option to sort the genes by their expression  levels in decreasing order. 

- #### 7 `SmallestNonAboveX()`
SmallestNonAboveX. replace small values with the next smallest value found, which is >X.

- #### 8 `AreTheseCellNamesTheSame()`
AreTheseCellNamesTheSame. Assert and compare two character vectors (e.g.: cell IDs) how much they overlap and  plot a Venn diagram. The function aborts with an error if overlap is too small.

- #### 9 `addToMiscOrToolsSlot()`
Add to Misc or Tools Slot. This function creates and adds a sub-slot to either the 'misc' or 'tools' slot of a  Seurat object. If the sub-slot already exists, it can either be overwritten or a warning will be issued. 

- #### 10 `showToolsSlots()`
Display Slots in the @tools of an Seurat Object. `showToolsSlots` prints the names of slots in the `@tools` of a given object.  It specifically targets list elements, skipping over data frames and other non-list objects. 

- #### 11 `showMiscSlots()`
Display Slots in the @misc of an Seurat Object. See `showToolsSlots` for details. Prints the names of slots in the `@misc` of a given object.  It specifically targets list elements, skipping over data frames and other non-list objects. 

- #### 12 `calc.q99.Expression.and.set.all.genes()`
calc.q99.Expression.and.set.all.genes. Calculate the gene expression of the e.g.: 99th quantile (expression in the top 1% cells).

- #### 13 `filterCodingGenes()`
Filter Coding Gene Symbols (or any matching input Patterns). This function filters out gene names that match specified patterns. It reports  the original and final number of gene symbols and the percentage remaining after filtering.  It filters out non-coding gene symbols by default. 

- #### 14 `filterExpressedGenes()`
Filter and Sort Gene Expression List Based on Specified Genes and Expression Threshold. This function takes a named list of gene expression values and a character vector of gene  symbols. It identifies the intersection of gene symbols with names in the list, filters genes based on a  specified expression threshold, and returns a character vector of genes that meet the criteria, sorted  by expression in descending order. 

- #### 15 `RenameClustering()`
RenameClustering. Rename clustering in a Seurat object.

- #### 16 `shorten_clustering_names()`
Shorten Clustering Names. This function takes in a string representing a clustering name,  and shortens it according to specific rules. It replaces "snn_res." with "",  "best.matching.names" with "bmatch", "ordered" with "ord",  "ManualNames" with "mNames", and ".long" at the end of the string with ".L". 

- #### 17 `getClusterNames()`
Retrieve Cluster Names. Extracts cluster names based on a specified identity class from a Seurat object. 

- #### 18 `GetClusteringRuns()`
GetClusteringRuns. The `GetClusteringRuns` function retrieves metadata column names associated with   clustering runs, based on a pattern to match, `"*snn_res.[0-9].[0-9]$"`, by default.

- #### 19 `GetNamedClusteringRuns()`
GetNamedClusteringRuns. The `GetNamedClusteringRuns` function retrieves metadata column names associated with   non-numeric ("named") clustering runs, based on a pattern to match, `"Name|name"`, by default.

- #### 20 `GetOrderedClusteringRuns()`
GetOrderedClusteringRuns. Get Clustering Runs: metadata column names.

- #### 21 `GetNumberOfClusters()`
GetNumberOfClusters. Print the number of clusters for each stored clustering run.

- #### 22 `calc.cluster.averages()`
calc.cluster.averages. Calculates the average of a metadata column (numeric) per cluster.

- #### 23 `plot.expression.rank.q90()`
plot.expression.rank.q90. Plot gene expression based on the expression at the 90th quantile  (so you will not lose genes expressed in few cells).

- #### 24 `set.mm()`
set.mm. Helps to find metadata columns. It creates a list with the names of of 'obj@meta.data'.

- #### 25 `ww.get.1st.Seur.element()`
Get the First Seurat Object from a List of Seurat Objects. Return the first Seurat object from a list, or the object itself    if a single Seurat object is supplied. 

- #### 26 `recall.all.genes()`
Recall all.genes global variable from a Seurat object. Recall \code{all.genes} from a Seurat object's \code{misc} slot,    which is stored by \code{calc.q99.Expression.and.set.all.genes()}, and optionally reset the global variable.

- #### 27 `recall.meta.tags.n.datasets()`
recall.meta.tags.n.datasets. Recall  meta.tags from obj@misc to "meta.tags" in the global environment.

- #### 28 `recall.parameters()`
recall.parameters. Recall parameters from obj@misc to "p" in the global environment.

- #### 29 `recall.genes.ls()`
recall.genes.ls. Recall genes.ls from obj@misc to "genes.ls" in the global environment.

- #### 30 `save.parameters()`
Save Parameters to Seurat Object. Stores a list of parameters within the `@misc$p` slot of a Seurat object,  allowing for easy reference and tracking of analysis parameters used. 

- #### 31 `create_scCombinedMeta()`
Create Single-Cell Metadata Object for a collection of Seurat Objects. This function creates a metadata object to correspond to a list of    single-cell experiments, for storing parent level information.    It initializes the object with the experiment and project name, and the    creation date. The created object is of class 'scMetadata_class'.

- #### 32 `copyMiscElements()`
Copy Specified Elements from One Seurat Object's @misc to Another's. Copies specified elements from the `@misc` slot of one Seurat object to the `@misc` slot  of another. It warns if some specified elements are missing in the source object or if elements are  overwritten in the destination object, depending on the `overwrite` argument. 

- #### 33 `copyCompleteToolsSlots()`
Copy Tools Slots from Multiple Seurat Objects. This function copies the `@tools` slots from a list of Seurat objects into a new slot  of a target Seurat object. This allows for the aggregation of tools information from multiple  experiments or datasets into a single, consolidated Seurat object. 

- #### 34 `subsetSeuObjByIdent()`
Subset a Seurat Object by Identity. Subsets a Seurat object based on a specified identity column and values. It allows    for an optional inversion of the selection. 

- #### 35 `downsampleSeuObj()`
downsampleSeuObj. Subset a compressed Seurat object and save it in the working directory.

- #### 36 `downsampleSeuObj.and.Save()`
downsampleSeuObj.and.Save. Downsample a Seurat object to a target fraction and save it.

- #### 37 `downsampleSeuObjByIdentAndMaxcells()`
Sample max number of Cells From each identity in a Seurat Object. This function samples a specified maximum number of cells from each identity class  in a Seurat object, in the meta.data. It ensures that the sampling does not exceed the total  number of cells available per identity. 

- #### 38 `RelabelSmallCategories()`
Relabel Small Categories / Clusters.   Relabels small categories in a specified metadata column of a Seurat object. Categories with  cell counts less than a minimum count are relabeled to a specified label. The function adds  a new metadata column with the updated labels. 

- #### 39 `removeResidualSmallClusters()`
Remove Residual Small Clusters from a Seurat Object. Removes clusters containing fewer cells than specified by `max.cells`  from a Seurat object. This function is particularly useful after subsetting a dataset,  where small, possibly unrepresentative clusters may remain. 

- #### 40 `dropLevelsSeurat()`
dropLevelsSeurat. Drop unused levels from `factor` variables in a Seurat object's meta.data.

- #### 41 `removeClustersAndDropLevels()`
Remove Clusters and Drop Levels from a List of Seurat Objects. This function removes residual small clusters from specified Seurat objects and  drops levels in factor-like metadata.

- #### 42 `removeCellsByUmap()`
Remove Cells by Dimension Reduction. This function applies a cutoff in the specified dimension of a given  dimension reduction (UMAP, PCA, or t-SNE) to remove cells.

- #### 43 `downsampleListSeuObjsNCells()`
Downsample a List of Seurat Objects to a Specific Number of Cells. Downsampling each Seurat object in a list to a specified number of cells. This function is  particularly useful for creating smaller, more manageable subsets of large single-cell datasets for  preliminary analyses or testing. 

- #### 44 `downsampleListSeuObjsPercent()`
Downsample a List of Seurat Objects to a Fraction. Downsampling a list of Seurat objects to a specified fraction of their original size.  This is useful for reducing dataset size for quicker processing or testing workflows. 

- #### 45 `Add.DE.combined.score()`
Add.DE.combined.score. Add a combined score to differential expression (DE) results. The score is  calculated as log-fold change (LFC) times negative logarithm of scaled  p-value (LFC * -log10( p_cutoff / pval_scaling )).

- #### 46 `StoreTop25Markers()`
Save Top 25 Markers per Cluster. Stores the top 25 markers for each cluster identified in a Seurat object, based on  the `avg_log2FC` from the output table of `FindAllMarkers()`. The result is saved under `@misc$df.markers$res...`,  rounding insignificant digits to three decimal places. 

- #### 47 `StoreAllMarkers()`
Store All Differential Expression Markers. Saves the complete output table from `FindAllMarkers()` to a Seurat object, facilitating  easy access to differential expression analysis results. This function rounds numerical values to a  specified number of digits to maintain readability and manage file sizes. 

- #### 48 `GetTopMarkersDF()`
Get Top Differential Expression Genes Data Frame. Retrieves a data frame of the top N differentially expressed genes from  differential gene expression analysis results, offering an option to exclude certain genes  based on patterns. 

- #### 49 `GetTopMarkers()`
Get Top Differential Expression Markers from DGEA Results. Retrieves the top N differentially expressed genes from the results of a differential  gene expression analysis, such as that provided by `FindAllMarkers()`. 

- #### 50 `AutoLabelTop.logFC()`
AutoLabelTop.logFC. Create a new "named identity" column in the metadata of a Seurat object,  with `Ident` set to a clustering output matching the `res` parameter of the function.  It requires the output table of `FindAllMarkers()`.  If you used `StoreAllMarkers()`  is stored under `@misc$df.markers$res...`, which location is assumed by default.

- #### 51 `AutoLabel.KnownMarkers()`
AutoLabel.KnownMarkers. Creates a new "named identity" column in the metadata of a Seurat object,   setting 'Ident' to a clustering output matching the 'res' parameter.   This function requires the output table of `FindAllMarkers()`.  If you used `StoreAllMarkers()`, the output is stored under `@misc$df.markers$res...`,  which is the default location.

- #### 52 `sparse.cor()`
Calculate Sparse Correlation Matrix. Computes a sparse correlation matrix from a given sparse matrix input. This function is  useful for efficiently handling large datasets where most values are zero, facilitating the calculation  of both covariance and correlation matrices without converting to a dense format. 

- #### 53 `Calc.Cor.Seurat()`
Calc.Cor.Seurat. Calculate gene correlation on a Seurat object.

- #### 54 `plot.Gene.Cor.Heatmap()`
Plot Gene Correlation Heatmap. Generates a heatmap visualization of gene correlations based on expression data.  Useful for identifying groups of genes that exhibit similar expression patterns across different conditions  or cell types in a Seurat object. 

- #### 55 `prefix_cells_seurat()`
Add Prefixes to Cell Names in Seurat Objects. Adds prefixes derived from a vector of identifiers to cell names in a list of Seurat objects.  This is useful for ensuring unique cell names across multiple samples or conditions when combining or comparing datasets. 

- #### 56 `find_prefix_in_cell_IDs()`
Check Prefix in Seurat Object Cell IDs. This function checks if a prefix has been added to the standard  cell-IDs (16 characters of A,TRUE,C,G) in a Seurat object. If so, it prints the number of unique prefixes found,  issues a warning if more than one unique prefix is found, and returns the identified prefix(es). 

- #### 57 `seu.Make.Cl.Label.per.cell()`
Create Cluster Labels for Each Cell. Generates labels for each cell by combining gene names and cluster IDs. This function  takes a named vector, typically representing top genes for clusters (values) and their corresponding  cluster IDs (names), along with a vector of cell IDs. It then creates a new vector where each cell  is labeled with its top gene and cluster ID in the format "GeneName.ClusterID". 

- #### 58 `GetMostVarGenes()`
Retrieve the Top Variable Genes from a Seurat Object. Retrieves the names of the most variable genes from a Seurat object,  typically used to focus subsequent analyses on genes with the greatest variation across cells. 

- #### 59 `gene.name.check()`
Check Gene Names in Seurat Object. Examines gene names in a Seurat object for specific naming conventions,  such as the presence of hyphens (-) or dots (.) often found in mitochondrial gene names.  This function is useful for ensuring gene names conform to expected patterns,  especially when preparing data for compatibility with other tools or databases. 

- #### 60 `check.genes()`
Check if Gene Names exist in Seurat Object or HGNC Database. Verifies the presence of specified gene names within a Seurat object or  queries them against the HGNC database. This function is useful for ensuring gene names are  correctly formatted and exist within the dataset or are recognized gene symbols. 

- #### 61 `fixZeroIndexing.seurat()`
Fix Zero Indexing in Seurat Clustering. Adjusts Seurat object metadata to fix zero-based cluster indexing, converting it to one-based indexing.  This function modifies a specified metadata column in the Seurat object to replace zero-indexed cluster names with one-based indexing. 

- #### 62 `CalculateFractionInTrome()`
Calculate Fraction of Genes in Transcriptome. Calculates the fraction of specified genes within the entire transcriptome of  each cell in a Seurat object.  This function is useful for assessing the relative abundance of a set of genes across cells,  such as identifying cells with high expression of marker genes. 

- #### 63 `AddNewAnnotation()`
AddNewAnnotation. This function creates a new metadata column based on an existing metadata column  and a list of mappings (name <- IDs).

- #### 64 `whitelist.subset.ls.Seurat()`
whitelist.subset.ls.Seurat. Subsets cells in a list of Seurat objects based on an externally provided list of cell IDs.

- #### 65 `FindCorrelatedGenes()`
FindCorrelatedGenes. Find correlated genes in a Seurat object

- #### 66 `UpdateGenesSeurat()`
Update Gene Symbols in a Seurat Object. This function updates gene symbols in a Seurat object based on current gene  nomenclature guidelines, using HGNChelper(). It checks and updates gene symbols to their  latest approved versions,ensuring that gene annotations are current and consistent.  The function optionally enforces unique gene symbols and provides statistics on the update process. 

- #### 67 `RenameGenesSeurat()`
Rename Gene Symbols in a Seurat Object. This function replaces gene names across various slots within a specified assay  of a Seurat object. It is designed to be run prior to any data integration or downstream analysis  processes. The function targets the `@counts`, `@data`, and `@meta.features` slots within  the specified assay, ensuring consistency in gene nomenclature across the object. 

- #### 68 `.check_and_rename()`
Check and Rename Gene Names in Seurat Assay Object. This function renames rows (genes) in a specified slot of a Seurat assay object.  It supports slots storing data as either a dense or a sparse matrix (dgCMatrix) or data.frame. 

- #### 69 `RemoveGenesSeurat()`
Remove Specific Genes from a Seurat Object. Removes specified genes from the metadata, counts, data, and scale.data slots of a Seurat object.  This operation is typically performed prior to data integration to ensure that gene sets are consistent  across multiple datasets. The function modifies the Seurat object in place. 

- #### 70 `HGNC.EnforceUnique()`
Enforce Unique HGNC Gene Symbols. Ensures that gene symbols are unique after being updated with HGNC symbols. This function  applies a suffix to duplicate gene symbols to enforce uniqueness. While using `make.unique` might not  be the ideal solution due to potential mismatches, it significantly reduces the number of mismatching  genes in certain scenarios, making it a practical approach for data integration tasks. 

- #### 71 `GetUpdateStats()`
Gene Symbol Update Statistics. Generates statistics on the gene symbol updates performed by `UpdateGenesSeurat()`.  This function analyzes the data frame of gene symbols before and after the update process,  providing insights into the proportion and total number of genes that were updated. 

- #### 72 `PlotUpdateStats()`
PlotUpdateStats. Creates a scatter plot of update statistics.

- #### 73 `Convert10Xfolders()`
Convert10Xfolders. This function takes a parent directory with a number of subfolders, each  containing the standard output of 10X Cell Ranger. It (1) loads the (filtered) data matrices,  (2) converts them to Seurat objects, and (3) saves them as .qs files 

- #### 74 `ConvertDropSeqfolders()`
ConvertDropSeqfolders. This function takes a parent directory with a number of subfolders, each  containing the standard output of 10X Cell Ranger. It (1) loads the filtered data matrices,  (2) converts them to Seurat objects, and (3) saves them as .RDS files.

- #### 75 `LoadAllSeurats()`
LoadAllSeurats. This function loads all Seurat objects found in a directory. It also works with  symbolic links (but not with aliases).

- #### 76 `read10x()`
Load 10X Genomics Data as Seurat Object. Reads 10X Genomics dataset files (gzipped) including matrix, features, and barcodes,  to a single expression matrix. This function handles the unzipping of these files, reads the data,  and re-compresses the files back to their original gzipped format. 

- #### 77 `.saveRDS.compress.in.BG()`
.saveRDS.compress.in.BG. Save and RDS object and compress resulting file in the background using system(gzip). OS X or unix.

- #### 78 `isave.RDS()`
isave.RDS. Save an RDS object, using a faster and efficient compression method that runs in the background.

- #### 79 `xsave()`
Save an R Object Using 'qs' Package for Fast Compressed Saving. This function saves an R object to a file in a quick and efficient format using the 'qs' package.  It constructs the file name based on various inputs and stores additional metadata if the object is a Seurat object.  The saving path can be adjusted by the presence of 'OutDir' in the global environment or defaults to the working directory. 

- #### 80 `xread()`
Read an R Object Using 'qs' Package for Fast Decompression. This function reads an R object from a file saved in a format specific to the 'qs' package,  which is designed for quick and efficient compression and decompression of R objects.  It also times the read operation, providing feedback on the duration of the operation. 

- #### 81 `  get_slurm_limit()`
Load a .qs object with optional SLURM-based safe-memory check (CBE only).   Loads a `.qs` serialized object (e.g., Seurat object or list) with an optional memory-safety  check that prevents loading objects larger than the SLURM job's allocated memory. The memory  check runs **only** when:  1) `safe_load = TRUE`  2) the global variable `onCBE` exists and is `TRUE`  3) the session is running inside a SLURM job with a defined memory limit   If no SLURM memory limit is detected, the safety check is skipped and a warning is shown. 

- #### 82 `isave.image()`
isave.image. Save an image of the current workspace using a faster and efficient compression  method that runs in the background.

- #### 83 `qsave.image()`
Save workspace - qsave.image. Save the workspace with external gzip compression (CPU intensive).

- #### 84 `find10XoutputFolders()`
Find 'Outs' Subdirectories in Specified Subdirectories. This function searches through specified subdirectories within a root directory  to find all subdirectories named 'outs' and returns a character vector with their full paths. 

- #### 85 `clip10Xcellname()`
Clip Suffixes from 10X Cell Names. Removes suffixes from cell names that are added by 10X technology and Seurat during data processing. 

- #### 86 `make10Xcellname()`
Add Suffix to Cell Names (e.g. lane suffix: _1). Appends a specified suffix to cell names to mimic lane suffixes used in 10X datasets. 

- #### 87 `plotTheSoup()`
plotTheSoup. Plot stats about the ambient RNA content in a 10X experiment. 

- #### 88 `jJaccardIndexVec()`
jJaccardIndexVec. Calculate jaccard similarity for 2 vectors. Helper to jPairwiseJaccardIndexList.

- #### 89 `jPairwiseJaccardIndexList()`
jPairwiseJaccardIndexList. Create a pairwise jaccard similarity matrix across all combinations of columns in  binary.presence.matrix. Modified from:  https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/

- #### 90 `jPresenceMatrix()`
jPresenceMatrix. Make a binary presence matrix from a list of vectors. Source:    https://stackoverflow.com/questions/56155707/r-how-to-create-a-binary-relation-matrix-from-a-list-of-strings

- #### 91 `jJaccardIndexBinary()`
jJaccardIndexBinary. Calculate the Jaccard Index for two binary vectors. Modified from:    https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/

- #### 92 `jPairwiseJaccardIndex()`
jPairwiseJaccardIndex. Create a pairwise jaccard similarity matrix across all combinations of columns in  binary.presence.matrix. Modified from:  https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/

- #### 93 `compareVarFeaturesAndRanks()`
Compare variable features and their ranks in two Seurat objects.. This function compares variable features (genes) between two Seurat objects,    reporting the number of genes in each, the percentage of common genes, the percentage    of unique genes in each object, and the similarity in the ranking of overlapping genes    using Spearman's rank correlation coefficient. Optionally, it can also generate a scatterplot    of the ranks of common genes using ggpubr's ggscatter. The function returns the common genes    and the Spearman's rank correlation coefficient. 

- #### 94 `.getNrCores()`
Get the number of CPUs to use for CBE processing. This function checks for the presence of a global `CBE.params` list and,  if found and contains a `cpus` entry, returns the number of CPUs specified by `cpus` minus one.  Otherwise, it returns a default number of CPUs. 

- #### 95 `.getNrPCs()`
Check List Elements. Tests if list elements are defined and reports the value or issues a warning. 

- #### 96 `.getRegressionVariablesForScaleData()`
Parse regression variables for name tagging. This function extracts the regression variables from the `@commands` slot of a Seurat object.  If no regression variables are found, a message is printed.

- #### 97 `.parseKeyParams()`
Parse key parameters from an object and format as a string. This function extracts the number of scaled features, the number of principal components,  and formats additional information including regression variables.

- #### 98 `.FindCommandInObject()`
Find Command in Seurat Object by Partial Match.   This function searches for commands in a list within a Seurat object using a partial match  (e.g., pattern matching) on the command names. It returns the content of the first match if only  one match is found. If multiple matches are found, it outputs the number of hits and their names. 

