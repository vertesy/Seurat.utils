## List of Functions in Seurat.Utils.Metadata.R (31) 
Updated: 2025/12/03 10:19
- #### 1 `addTranslatedMetadata()`
Add Translated Metadata to a Seurat Object.   This function translates a specified metadata vector in a Seurat object using a named vector of old  and new values and adds it to the Seurat object. A suffix can optionally be appended to the name of  the newly created metadata column. The function also generates UMAP plots for the new metadata. 

- #### 2 `getMetaColnames()`
Get Metadata Column Names Matching Pattern. Retrieves column names from an object's metadata that match a specified pattern. 

- #### 3 `metaColnameExists()`
Check if a Column Exists in the Metadata of an S4 Object. This function checks whether a given column exists in the meta.data of a Seurat object.

- #### 4 `getMetadataColumn()`
getMetadataColumn. Retrieves a specified metadata column from a Seurat object and returns it as a named vector.

- #### 5 `get_levels_seu()`
Get Unique Levels of a Seurat Object Ident Slot. This function extracts the unique levels present in the 'ident' slot of a Seurat object.  The function throws an error if the number of levels exceeds 'max_levels'.  The function optionally prints the R code to recreate the 'Levels' vector using 'dput'. 

- #### 6 `calculateAverageMetaData()`
Calculate Average Metadata for Seurat Object. Computes specified metrics (e.g., median, mean) for given metadata features across each category  defined by an identity column in a Seurat object's metadata. This function allows for flexible  metric calculation on specified features, providing insights into the data distribution. 

- #### 7 `calculatePercentageMatch()`
Calculate the Percentage of Matches per Category. This function calculates the percentage of matches for specified metadata features  against provided match values within each category of an identifier in a Seurat object. 

- #### 8 `getMedianMetric.lsObj()`
getMedianMetric.lsObj. Get the median values of different columns in meta.data, can iterate over a list of Seurat objects.

- #### 9 `getCellIDs.from.meta()`
getCellIDs.from.meta. Retrieves cell IDs from a specified metadata column of a Seurat object, where the cell ID matches a provided list of values.    The matching operation uses the `%in%` operator and can handle `NA`/`NaN` values.

- #### 10 `addMetaDataSafe()`
Add Metadata to a Seurat object, safely with Checks. Wrapper function for `AddMetaData` that includes additional checks and assertions. 

- #### 11 `create.metadata.vector()`
Create a Metadata Vector. This function creates a metadata vector from an input vector and a Seurat object.  The resulting vector contains values from 'vec' for the intersecting cell names between 'vec' and 'obj'.  It also checks if the intersection between the cell names in 'vec' and 'obj' is more than a  minimum intersection size.

- #### 12 `addMetaFraction()`
addMetaFraction. Add a new metadata column to a Seurat object, representing the fraction of a gene set in the transcriptome (expressed as a percentage).

- #### 13 `addGeneClassFractions()`
Add Metadata for Gene-Class Fractions.   This function adds metadata for various gene-class fractions such as percent.mito, percent.ribo,  percent.AC.GenBank, percent.AL.EMBL, percent.LINC, percent.MALAT1, and percent.HGA to a Seurat object.  If the metadata already exists, a message will be displayed. 

- #### 14 `add.meta.tags()`
add.meta.tags. Add metadata tags to a Seurat object dataset.

- #### 15 `seu.add.meta.from.table()`
seu.add.meta.from.table. Add metadata columns from a table keyed by cell names.

- #### 16 `seu.map.and.add.new.ident.to.meta()`
seu.map.and.add.new.ident.to.meta. Adds a new metadata column to a Seurat object based on an identity mapping table.

- #### 17 `fix.orig.ident()`
fix.orig.ident. Remove the string "filtered_feature_bc_matrix." from "orig.ident". Helper function.

- #### 18 `seu.RemoveMetadata()`
seu.RemoveMetadata. Remove specified metadata columns from a Seurat object.

- #### 19 `saveLsSeuratMetadata()`
Save Metadata from a List of Seurat Objects. This function takes a list of Seurat objects, extracts their metadata, and saves it to a file with a specified suffix. 

- #### 20 `transferMetadata()`
Transfer Multiple Metadata Columns Between Two Seurat Objects. Transfers specified metadata columns from one Seurat object to another,  with options for verbose output and overwriting existing columns. Checks for cell overlap and  reports percentages of matching and unique cells. 

- #### 21 `sampleNpc()`
Sample N % of a dataframe (obj@metadata), and return rownames (cell IDs).. This function samples a specified percentage of a dataframe (specifically a subset    of the metadata of a Seurat object) and returns the corresponding cell IDs.

- #### 22 `merge_seurat_metadata()`
Merge Seurat Metadata. Merges the `@metadata` from a list of Seurat objects, binds them by row, and applies optional inclusion/exclusion of columns.

- #### 23 `writeCombinedMetadataToTsvFromLsObj()`
Combine Metadata from a list of Seurat objects and Write to TSV.   Formerly `writeMetadataToTsv`. `writeCombinedMetadataToTsvFromLsObj` takes a list of ls.Obj, extracts their `@meta.data` slots,  removes specified columns, checks for column consistency, creates a barplot showing the number  of rows per object, and finally merges these into one large data frame. 

- #### 24 `plotMetadataCorHeatmap()`
Plot Metadata Correlation Heatmap. This function plots a heatmap of metadata correlation values. It accepts a Seurat object  and a set of metadata columns to correlate. The correlations are calculated using either Pearson  or Spearman methods, and the resulting heatmap can include the principal component (PCA) values  and be saved with a specific suffix. 

- #### 25 `heatmap_calc_clust_median()`
Calculate and plot heatmap of cluster medians. This function calculates the median of specified variables in a dataframe,  grouped by a column ('ident'). The function also provides an option to scale the medians,  subset the ident levels, and either return a matrix of median values or plot a heatmap. 

- #### 26 `plotMetadataMedianFractionBarplot()`
plotMetadataMedianFractionBarplot. Generates a barplot of metadata median values.

- #### 27 `plotMetadataCategPie()`
Plot Metadata Category Pie Chart. Generates a pie chart visualizing the distribution of categories within a specified  metadata column of a Seurat object. 

- #### 28 `renameAzimuthColumns()`
Rename Azimuth Columns in Seurat Object. Dynamically renames specified metadata columns in a Seurat object, particularly those  prefixed with "predicted." and the "mapping.score" column, by applying a new prefix  that combines a user-defined prefix and a reference name. 

- #### 29 `renameSmallCategories()`
Rename Small Categories in Seurat Object Metadata. This function renames categories within a specified identity column of a  Seurat object's metadata that have fewer cells than a specified minimum threshold.  Categories below this threshold are renamed to a common name, typically "unclear",  to clean up small, potentially noisy categories.

- #### 30 `.metaColnames()`
Transfer labels from a reference Seurat object to a query Seurat object. Function to transfer labels from a reference Seurat object to a query Seurat object  using anchoring and transfer data methods from the Seurat package. It then visualizes the  reference and the combined objects using Uniform Manifold Approximation and Projection (UMAP). 

- #### 31 `matchBestIdentity()`
Match and Translate Best Identity. This function matches the best identity from `ident_to_rename` to `reference_ident` in an object,  in other words, it replaces original categories with the most frequent ones from the reference,  hence helps to filter out less important categories. 

