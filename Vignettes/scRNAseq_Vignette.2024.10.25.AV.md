
---
title: "Visualization Techniques for scRNAseq Analysis with Seurat.utils"
subtitle: "!WORK IN PROGRESS! - Demonstration of the Visualization section of the Seurat.utils package."
author: "Richard and Abel"
date: "2024-10-29"
knit_root_dir: ~/Documents/Analysis_R/Vignette_SU/
output:
  html_document:
    keep_md: true    # Optional, keeps .md output for debugging
  pdf_document:
    includes:
      in_header: preamble.tex
---

<!-- file.edit('~/GitHub/Zacc/Vignette.seurat.utils/scRNAseq_Vignette.2024.10.25.AV.Rmd') -->
<!-- knitr::purl("~/GitHub/Zacc/Vignette.seurat.utils/scRNAseq_Vignette.2024.10.25.AV.Rmd", documentation = 2, -->
<!--             output = "~/GitHub/Zacc/Vignette.seurat.utils/scRNAseq_Vignette.2024.10.25.AV.R") -->
<!-- file.edit("~/GitHub/Zacc/Vignette.seurat.utils/scRNAseq_Vignette.2024.10.25.AV.R") -->

<!-- "Replace:" -->
<!-- FROM knitr::include_graphics\("([^"]*)"  -->
<!-- TO knitr::include_graphics(\n # "\1"\npaste0(OutDir, ""  -->


<!-- git remote set-url origin git@github.com:vertesy/Vignette.seurat.utils.git -->




# Vignette for Visualtization with the `Seurat.utils` package

## Introduction

`Seurat.utils` enhances scRNAseq analysis by building on the `Seurat` package, and introducing advanced and convenient tools for typical analytic needs.
This vignette focuses on the `Seurat.Utils.Visualization.R` component, demonstrating how it streamlines and enriches the exploration of single-cell RNA sequencing data.

### Known caveats

I wanted to demonstrate `Seurat.utils` on a real, published Seurat object. Therefore, I chose a an integrated object from our previous publication, [Gruffi](https://www.embopress.org/doi/full/10.15252/embj.2022111118).

This object created the time before Seurat v5 existed. I used `SeuratObject::UpdateSeuratObject()` to update the object to v5.During the creation of this vignettes, we realized that this function actually creates a "hybrid" object between v5 and v3, where only certain assays, etc. are updated. I contacted the authors, but they declared this is a feature not bug.

As I did not want to write `UpdateSeuratObjectThoroughly()`, so I was stuck with this mutant v3-v5 object, and therefore I had to hack some of the
functions which would otherwise work perfectly on a clean v5 object. You will see signs of this below.

Apologies to the readers.


## Differences Between Seurat and Seurat.utils

While `Seurat` offers comprehensive tools for scRNAseq data analysis, `Seurat.utils` extends these functions with specialized or more convenient tools.

Examples:
- **Custom UMAP Plots**: Automatic file saving (png/pdf/jpeg), automatic annotation, better defaults, and more.
- **PCA Variance Explained**: Tools for a more detailed examination of variance explained by principal components, aiding in the interpretation of dimensionality reduction results.

## Loading Essential Packages

If `Seurat.utils` is not installed, installation instructions are available at [github.com/vertesy/Seurat.utils](https://github.com/vertesy/Seurat.utils).

Before diving into data analysis, ensure that `Seurat.utils` is installed and loaded alongside `MarkdownReports` for comprehensive reporting capabilities.


``` r
rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T); gc()
```

```
## null device 
##           1
```

```
##           used (Mb) gc trigger (Mb) limit (Mb) max used (Mb)
## Ncells  901636 48.2    1758464   94         NA  1280591 68.4
## Vcells 1513305 11.6    8388608   64     102400  2523139 19.3
```

``` r
library(dplyr)
# library(png)
require("Stringendo")
require("CodeAndRoll2")
require("MarkdownHelpers")
require("MarkdownReports")
require("ggExpress")
# require("Seurat.utils")
source('~/GitHub/Packages/Rocinante/R/Rocinante.R')
```

```
## [1] "Loading Rocinante custom function library."
## [1] "Depends on CodeAndRoll2, MarkdownReports, gtools, readr, gdata, clipr. Some functions depend on other libraries."
```

``` r
r$Seurat.utils()
OutDir <- "~/Documents/Analysis_R/Vignette_SU/"
# OutDir <- "~/GitHub/Zacc/Vignette.seurat.utils"
setwd(OutDir)
```

## Reading scRNAseq Data


``` r
# combined.obj <- xread("C:\\gruffi\\obj.Fig.4C.clean_8200.cells_Vignette.seurat.utils_2024.03.27_10.15.qs")
# combined.obj <- read("~/Downloads/obj_8200.cells_Fig.4C.clean_Vignette.seurat.utils_2024.06.03_21.00.qs")
# obj.small <- downsampleSeuObj(combined.obj, nCells = 2500)
# xsave(obj.small, v = T, showMemObject = F)
combined.obj <- xread('~/Documents/Analysis_R/Vignette_SU/obj.small_2500.cells_Vignette.seurat.utils_2024.10.25_13.46.qs')
```

```
##    MALAT1    TMSB4X    MT-CO1    TUBA1A    TMSB10     RPS19 
## 1.0000000 0.9999625 0.9999251 0.9998876 0.9998501 0.9998127 
## [1] "Seurat with 2500 cells & 41 meta colums."
## xread: 1.869 sec elapsed
```


``` r
# cc.genes <- list(
#   s.genes = c("MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", 
#   "RRM1", "UNG", "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", 
#   "MLF1IP", "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", 
#   "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", 
#   "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", 
#   "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8"), 
#   g2m.genes = c("HMGB2", 
#   "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", 
#   "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "FAM64A", 
#   "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", 
#   "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "HN1", 
#   "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", 
#   "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", "PSRC1", 
#   "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", 
#   "CBX5", "CENPA")
#   )
# combined.obj <- CellCycleScoring(combined.obj, s.features = cc.genes$'s.genes', g2m.features = cc.genes$'g2m.genes')
```





``` r
(identX <- GetClusteringRuns(combined.obj)[1])
```

```
## c("integrated_snn_res.0.1", "integrated_snn_res.0.2", "integrated_snn_res.0.3", 
## "integrated_snn_res.0.4", "integrated_snn_res.0.5")
```

```
## [1] "integrated_snn_res.0.1"
```

``` r
(ident2 <- GetNamedClusteringRuns(combined.obj)[2])
```

```
## c("cl.names.top.gene.res.0.2", "cl.names.KnownMarkers.0.2", "cl.names.top.gene.res.0.5", 
## "cl.names.KnownMarkers.0.5")
```

```
## [1] "cl.names.KnownMarkers.0.2"
```

``` r
raster <- if (ncol(combined.obj) > 1e+05) TRUE else FALSE
nr.Col <- 2
nr.Row <- 4
wA4 <- 8.27
hA4 <- 11.69
list.of.genes <- c("MALAT1", "TMSB4X", "MT-CO1")

ls.Seu <- list("Exp1" = combined.obj,
               "Exp2" = downsampleSeuObj(combined.obj, fractionCells = .2),
               "Exp3" = downsampleSeuObj(combined.obj, fractionCells = .3))
```

```
## [1] "500 or 20% of the cells are kept. Seed: 1989"
## [1] "750 or 30% of the cells are kept. Seed: 1989"
```

A quick summary to understand dataset dimensions and composition:


``` r
stopifnot(exists("ident2"))
combined.obj
```

```
## An object of class Seurat 
## 28690 features across 2500 samples within 2 assays 
## Active assay: integrated (2000 features, 2000 variable features)
##  2 layers present: data, scale.data
##  1 other assay present: RNA
##  2 dimensional reductions calculated: pca, umap
```



<!-- # ---------------------------------------------------------------------------------------------------- -->
<!-- # TEST AREA BELOW ------------------------------------------------------------------------------------ -->
<!-- # ---------------------------------------------------------------------------------------------------- -->



``` r
scBarplot.CellFractions(obj = combined.obj, group.by = ident2, fill.by = "Phase")
```

```
##           
##             G1 G2M   S
##   0.MKI67    0 137  89
##   1.ID2    224  53  75
##   10.OPCML 162   9  37
##   2.HES6   111  53  74
##   3.MAF     15   4   3
##   4.BNIP3  163  52  97
##   5.MEIS2  251  34 104
##   6.POLR2A 103  34  72
##   7.ERBB4  141  48 114
##   8.MEF2C  139  23  23
##   9.ZFHX3   40   1  15
## [1] "/Users/abel.vertesy/Documents/Analysis_R/Vignette_SU/Cell.proportions.of.Phase.by.cl.names.KnownMarkers.0.2.fr.barplot.png"
```

```
## 1.035 sec elapsed
```

![](scRNAseq_Vignette.2024.10.25.AV_files/figure-html/unnamed-chunk-1-1.png)<!-- -->

``` r
# knitr::include_graphics("Cell.proportions.of.Phase.by.cl.names.top.gene.res.0.2.downsampled.fr.barplot.png")
```

<!-- # ---------------------------------------------------------------------------------------------------- -->
<!-- # TEST AREA ABOVE ------------------------------------------------------------------------------------ -->
<!-- # ---------------------------------------------------------------------------------------------------- -->

## UMAPs

### `SetupReductionsNtoKdimensions()`

SetupReductionsNtoKdimensions is a function that calculates N-to-K dimensional UMAPs, executing specified dimensionality reduction (UMAP, tSNE, or PCA) over a range of dimensions and backing up the results within a Seurat object.


``` r
# combined.obj <- SetupReductionsNtoKdimensions(
#   obj = combined.obj,
#   nPCs = 30,
#   dimensions = 3:2,
#   reduction_input =  "pca",
#   reduction_output =  "umap"
# )
```


### `BackupReduction()`

`BackupReduction` stores a backup of specified dimensionality reduction data (e.g., UMAP, tSNE, PCA) within the *Seurat* object.


``` r
combined.obj <- BackupReduction(obj = combined.obj, dim = 2, reduction = "umap")
```

### `RecallReduction()`

`RecallReduction` restores dimensionality reduction data (e.g., UMAP, tSNE, PCA) from a backup stored within **obj@misc$reductions.backup** to the active **obj@reductions slot**.


``` r
combined.obj <- RecallReduction(obj = combined.obj, dim = 2, reduction = "umap")
```

```
## [1] "2 dimensional umap from obj@misc$reductions.backup is set active. "
```



### `qUMAP()`

`qUMAP` is a wrapper function for `Seurat::FeaturePlot` that allows for quick visualization of UMAPs colored by a **numeric** features (metadata columns) or *genes*.


``` r
qUMAP(feature = "nFeature_RNA")
```

```
## [1] "/Users/abel.vertesy/Documents/Analysis_R/Vignette_SU/UMAP.nFeature_RNA.RNA.2500c.png"
```

![](scRNAseq_Vignette.2024.10.25.AV_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

``` r
qUMAP(feature = "TOP2A", PNG = FALSE, save.plot = TRUE)
```

```
## [1] "/Users/abel.vertesy/Documents/Analysis_R/Vignette_SU/UMAP.TOP2A.RNA.2500c.pdf"
```

![](scRNAseq_Vignette.2024.10.25.AV_files/figure-html/unnamed-chunk-5-2.png)<!-- -->


### `clUMAP()`

`clUMAP` is a wrapper function for `Seurat::DimPlot` that allows for quick visualization of UMAPs colored by categorical features from metadata columns, e.g clustering results.


``` r
clUMAP(ident = "integrated_snn_res.0.1", cols = RColorBrewer::brewer.pal(7, "Dark2"))
```

```
## [1] "/Users/abel.vertesy/Documents/Analysis_R/Vignette_SU/UMAP.integrated_snn_res.0.1.2500c.png"
```

```
## 0.918 sec elapsed
```

![](scRNAseq_Vignette.2024.10.25.AV_files/figure-html/unnamed-chunk-6-1.png)<!-- -->


### `FlipReductionCoordinates()`

`FlipReductionCoordinates` reverses the dimensionality reduction coordinates (such as UMAP or tSNE) vertically or horizontally to alter the visualization perspective.


``` r
clUMAP(obj = combined.obj)
```

```
## [1] "/Users/abel.vertesy/Documents/Analysis_R/Vignette_SU/UMAP.cl.names.top.gene.res.0.2.2500c.png"
```

```
## 0.964 sec elapsed
```

![](scRNAseq_Vignette.2024.10.25.AV_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

``` r
combined.obj <- FlipReductionCoordinates(
  obj = combined.obj,
  dim = 2,
  reduction = "umap",
  flip = c("x", "y", "xy", NULL)[1],
  FlipReductionBackupToo = FALSE # If you want to keep backup of the original coordinates, set to FALSE
)
clUMAP(obj = combined.obj, sub = "flipped x axis in the UMAP coordinates")
```

```
## [1] "/Users/abel.vertesy/Documents/Analysis_R/Vignette_SU/UMAP.cl.names.top.gene.res.0.2.2500c.flipped.x.axis.in.the.UMAP.coordinates.png"
```

```
## 0.941 sec elapsed
```

![](scRNAseq_Vignette.2024.10.25.AV_files/figure-html/unnamed-chunk-7-2.png)<!-- -->

### `AutoNumber.by.UMAP()`

`AutoNumber.by.UMAP` automatically renumbers clusters based on their position along a specified dimension in a UMAP (or tSNE or PCA) plot, potentially enhancing interpretability by ordering clusters.


``` r
combined.obj <- AutoNumber.by.UMAP(
  obj = combined.obj, reduction = "umap", 
  dim = 1, swap = TRUE,  # Swap the order along dimension 1, to get a pseudotime-like ordering from progenitor to differentiated cells
  ident = ident2,
  obj.version = 3,  # You will not needed unless updating via `SeuratObject::UpdateSeuratObject()`.
  plot = TRUE
)
```

```
## [1] "NewMetaCol: cl.names.KnownMarkers.0.2.ordered"
## [1] "/Users/abel.vertesy/Documents/Analysis_R/Vignette_SU/UMAP.cl.names.KnownMarkers.0.2.ordered.2500c.png"
```

```
## 0.911 sec elapsed
```

![](scRNAseq_Vignette.2024.10.25.AV_files/figure-html/unnamed-chunk-8-1.png)<!-- -->


### `multiSingleClusterHighlightPlots.A4()`

`multiSingleClusterHighlightPlots.A4` generates and saves cluster highlight plots for both single and multiple clusters using UMAP or other dimensionality reduction techniques, supporting various plot formats and customization options, and ensuring an A4 paper size for the output.


``` r
# scBarplot.CellsPerCluster(ident = identX, obj = combined.obj)
r$Seurat.utils()
multiSingleClusterHighlightPlots.A4(ident = identX, obj = combined.obj)
```

```
## [1] "integrated_snn_res.0.1-umap/"
## [1] "All files will be saved under 'NewOutDir':  /Users/abel.vertesy/Documents/Analysis_R/Vignette_SU/integrated_snn_res.0.1-umap/"
## [1] "ParentDir will be:"
## [1] "ParentDir defined as:"
## [1] "Call *create_set_Original_OutDir()* when chaning back to the main dir."
## [1] "OutDir defined as:"
## [1] "b.Subdirname defined as:"
## [1] "page: 1 | clusters 0, 1, 2, 3"
## 0.05 sec elapsed
## 0.052 sec elapsed
## 0.048 sec elapsed
## 0.047 sec elapsed
```

```
## [1] "page: 2 | clusters 4"
## 0.049 sec elapsed
```

```
## [1] "All files will be saved under 'OutDir':"
## [1] "OutDir defined as:"
## 2.124 sec elapsed
```

<img src="RNA_snn_res.0.1-umap/umap.1.clusters.1.2.3.4.jpg" width="100%" /><img src="RNA_snn_res.0.1-umap/umap.2.clusters.5.6.7.8.jpg" width="100%" />
 

<img src="../../../Documents/Analysis_R/Vignette_SU/umap.2.clusters.5.jpg" width="100%" />


### `qClusteringUMAPS()`

`qClusteringUMAPS` generates and arranges UMAP plots for up to four specified clustering resolutions from a *Seurat* object onto an A4 page, facilitating comparative visualization.


``` r
ident2 <- na.omit.strip(GetClusteringRuns(combined.obj)[1:2])
```

```
## c("integrated_snn_res.0.1", "integrated_snn_res.0.2", "integrated_snn_res.0.3", 
## "integrated_snn_res.0.4", "integrated_snn_res.0.5")
```

``` r
qClusteringUMAPS( obj = combined.obj, idents = ident2)
```

```
## 0.387 sec elapsed
## 0.094 sec elapsed
## c("integrated_snn_res.0.1", "integrated_snn_res.0.2", "integrated_snn_res.0.3", 
## "integrated_snn_res.0.4", "integrated_snn_res.0.5")
## [1] "Identity not found. Plotting integrated_snn_res.0.1 \n"
## 0.085 sec elapsed
## c("integrated_snn_res.0.1", "integrated_snn_res.0.2", "integrated_snn_res.0.3", 
## "integrated_snn_res.0.4", "integrated_snn_res.0.5")
## [1] "Identity not found. Plotting integrated_snn_res.0.1 \n"
## 0.089 sec elapsed
```

```
## [1] "/Users/abel.vertesy/Documents/Analysis_R/Vignette_SU/Clustering.UMAP.Res_0.1_0.2.png"
```

<img src="../../../Documents/Analysis_R/Vignette_SU/Clustering.UMAP.Res_0.1_0.2.png" width="100%" />

### `plotQUMAPsInAFolder()`

`plotQUMAPsInAFolder` plots qUMAPs for a specified set of genes, storing the results in a specified folder, with the option to default to using the gene set name if no folder name is provided.


``` r
plotQUMAPsInAFolder(
  genes = c("MT-CO1", "TUBA1A"), obj = combined.obj,
  foldername = "MyGenePlots", intersectionAssay = "RNA",
  plot.reduction = "umap"
)
```

```
## [1] "MyGenePlots-umap/"
## [1] "All files will be saved under 'NewOutDir':  /Users/abel.vertesy/Documents/Analysis_R/Vignette_SU/MyGenePlots-umap/"
## [1] "ParentDir was defined as:"
## [1] "ParentDir will be:"
## [1] "ParentDir defined as:"
## [1] "Call *create_set_Original_OutDir()* when chaning back to the main dir."
## [1] "OutDir defined as:"
## [1] "b.Subdirname defined as:"
## check.genes: 0.007 sec elapsed
## [1] "/Users/abel.vertesy/Documents/Analysis_R/Vignette_SU/MyGenePlots-umap/UMAP.MT-CO1.RNA.2500c.png"
```

```
## [1] "/Users/abel.vertesy/Documents/Analysis_R/Vignette_SU/MyGenePlots-umap/UMAP.TUBA1A.RNA.2500c.png"
```

```
## [1] "All files will be saved under 'OutDir':"
## [1] "OutDir defined as:"
```

``` r
"Creates individual plots in */MyGenePlots/*.png"
```

```
## [1] "Creates individual plots in */MyGenePlots/*.png"
```

<img src="../../../Documents/Analysis_R/Vignette_SU/MyGenePlots-umap/UMAP.MT-CO1.RNA.2500c.png" width="100%" />

<img src="../../../Documents/Analysis_R/Vignette_SU/TUBA1A.png" width="100%" />

### `umapHiLightSel()`

The `umapHiLightSel` function generates a UMAP plot from a **Seurat** object, highlighting specified clusters, and saves the resulting plot directly to the current working directory.


``` r
umapHiLightSel(
  obj = combined.obj,
  COI = c("0", "2", "4"),
  ident = GetClusteringRuns()[1]
  )
```

```
## c("integrated_snn_res.0.1", "integrated_snn_res.0.2", "integrated_snn_res.0.3", 
## "integrated_snn_res.0.4", "integrated_snn_res.0.5")
## [1] "1514 cells found."
```

```
## [1] "cells.0.2.4"
```

![](scRNAseq_Vignette.2024.10.25.AV_files/figure-html/unnamed-chunk-17-1.png)<!-- -->

# ----------------------------------------------------------------------------------------------------

## Barplots and Histograms

### `scBarplot.FractionAboveThr()`

`scBarplot.FractionAboveThr` draws a barplot of the fraction of cells above a certain threshold for a given feature.


``` r
scBarplot.FractionAboveThr(id.col = identX, value.col = "percent.ribo", thrX = 0.1)
```

```
## # A tibble: 5 × 3
##   value names colour
##   <dbl> <chr> <lgl> 
## 1  47.2 0     TRUE  
## 2  37.5 1     FALSE 
## 3  64.6 2     TRUE  
## 4  32.5 3     FALSE 
## 5  59.1 4     TRUE  
## [1] "/Users/abel.vertesy/Documents/Analysis_R/Vignette_SU/Pc.cells.above.percent.ribo.of.0.1.integrated_snn_res.0.1.png"
```

```
## 0.371 sec elapsed
```

![](scRNAseq_Vignette.2024.10.25.AV_files/figure-html/unnamed-chunk-18-1.png)<!-- -->![](scRNAseq_Vignette.2024.10.25.AV_files/figure-html/unnamed-chunk-18-2.png)<!-- -->

### `scBarplot.FractionBelowThr()`

`scBarplot.FractionBelowThr` generates a bar plot to visualize the percentage of cells within each cluster that fall below a specified threshold, according to a metadata column value.


``` r
scBarplot.FractionBelowThr(id.col = identX, value.col = "percent.ribo", thrX = 0.1)
```

```
## # A tibble: 5 × 3
##   value names colour
##   <dbl> <chr> <lgl> 
## 1  52.8 0     FALSE 
## 2  62.5 1     TRUE  
## 3  35.4 2     FALSE 
## 4  67.5 3     TRUE  
## 5  40.9 4     FALSE 
## [1] "/Users/abel.vertesy/Documents/Analysis_R/Vignette_SU/Pc.cells.below.percent.ribo.of.0.1.integrated_snn_res.0.1.png"
```

```
## 0.322 sec elapsed
```

![](scRNAseq_Vignette.2024.10.25.AV_files/figure-html/unnamed-chunk-19-1.png)<!-- -->![](scRNAseq_Vignette.2024.10.25.AV_files/figure-html/unnamed-chunk-19-2.png)<!-- -->
### `PercentInTranscriptome()`


`PercentInTranscriptome` computes and visualizes gene expression levels relative to total UMI counts, emphasizing highly expressed genes' contribution to the transcriptome. The first argument of the `readPNG()` function is the image path. Additionally, you can provide only the file name if the working directory is set to the folder containing the image.


``` r
PercentInTranscriptome(combined.obj, n.genes.barplot = 25)
```

```
##    FTH1  MT-CO1  MT-CO3   GAPDH     FTL  MT-CO2  MT-ND4   FABP7  MT-CYB C1orf61 
##   0.758   0.746   0.605   0.580   0.577   0.570   0.564   0.476   0.470   0.452 
## MT-ATP6  TUBA1B    TPI1   STMN2     CKB   HMGB1     JUN   H2AFZ  MT-ND2   SOX11 
##   0.451   0.443   0.422   0.400   0.391   0.387   0.372   0.363   0.357   0.355 
##   SERF2     PKM    ENO1   STMN4   DDIT4 
##   0.354   0.352   0.351   0.336   0.336 
## [1] "/Users/abel.vertesy/Documents/Analysis_R/Vignette_SU/Gene.expression.as.fraction.of.all.transcripts.UMI.s.logY.hist.png"
```

```
## 0.368 sec elapsed
## # A tibble: 25 × 3
##    value names   colour
##    <dbl> <chr>   <chr> 
##  1 0.758 FTH1    1     
##  2 0.746 MT-CO1  1     
##  3 0.605 MT-CO3  1     
##  4 0.58  GAPDH   1     
##  5 0.577 FTL     1     
##  6 0.57  MT-CO2  1     
##  7 0.564 MT-ND4  1     
##  8 0.476 FABP7   1     
##  9 0.47  MT-CYB  1     
## 10 0.452 C1orf61 1     
## # ℹ 15 more rows
## [1] "/Users/abel.vertesy/Documents/Analysis_R/Vignette_SU/Percentage.of.highest.expressed.genes.bar.png"
```

```
## 0.413 sec elapsed
```

```
## An object of class Seurat 
## 28690 features across 2500 samples within 2 assays 
## Active assay: integrated (2000 features, 2000 variable features)
##  2 layers present: data, scale.data
##  1 other assay present: RNA
##  2 dimensional reductions calculated: pca, umap
```

``` r
# showing the generated plot
```

<img src="Percentage.of.highest.expressed.genes.bar.png" width="100%" />

<img src="Gene.expression.as.fraction.of.all.UMI.s.logY.hist.png" width="100%" />

### `plotGeneExpressionInBackgroundHist()`



``` r
plotGeneExpressionInBackgroundHist(gene = "HMGB2", obj = combined.obj)
```

```
## [1] "/Users/abel.vertesy/Documents/Analysis_R/Vignette_SU/HMGB2.and.the.normalised.logtransformed.transcript.count.distribution.logX.hist.png"
```

```
## 0.375 sec elapsed
```

![](scRNAseq_Vignette.2024.10.25.AV_files/figure-html/unnamed-chunk-23-1.png)<!-- -->

### `scBarplot.FractionAboveThr()`

scBarplot.FractionAboveThr generates a bar plot depicting the percentage of cells within each cluster that exceed a specified threshold, based on a selected metadata column.


``` r
scBarplot.FractionAboveThr(id.col = "cl.names.top.gene.res.0.2", 
                           , value.col = "percent.ribo", thrX = 0.1)
```

```
## # A tibble: 11 × 3
##    value names        colour
##    <dbl> <chr>        <lgl> 
##  1  62.8 0.HMGB2      TRUE  
##  2  61.1 1.CLU        TRUE  
##  3  28.8 10.OPCML     FALSE 
##  4  71.0 2.HES6       TRUE  
##  5  68.2 3.TTR        TRUE  
##  6  61.2 4.BNIP3      TRUE  
##  7  31.9 5.CNTNAP2    FALSE 
##  8  36.8 6.AL118516.1 FALSE 
##  9  40.9 7.DLX6-AS1   FALSE 
## 10  31.4 8.MEF2C      FALSE 
## 11  10.7 9.GRIA4      FALSE 
## [1] "/Users/abel.vertesy/Documents/Analysis_R/Vignette_SU/Pc.cells.above.percent.ribo.of.0.1.cl.names.top.gene.res.0.2.png"
```

```
## 0.335 sec elapsed
```

![](scRNAseq_Vignette.2024.10.25.AV_files/figure-html/unnamed-chunk-24-1.png)<!-- -->![](scRNAseq_Vignette.2024.10.25.AV_files/figure-html/unnamed-chunk-24-2.png)<!-- -->

``` r
# knitr::include_graphics("Pc.cells.above.percent.ribo.of.0.1.cl.names.top.gene.res.0.2.png")
```

### `scBarplot.FractionBelowThr()`

`scBarplot.FractionBelowThr` generates a bar plot to visualize the percentage of cells within each cluster that fall below a specified threshold, based on a selected metadata column value.


``` r
scBarplot.FractionBelowThr(id.col = "cl.names.top.gene.res.0.2", value.col = "percent.ribo", thrX = 0.1, palette_use ="npg" )
```

```
## # A tibble: 11 × 3
##    value names        colour
##    <dbl> <chr>        <lgl> 
##  1  37.2 0.HMGB2      FALSE 
##  2  38.9 1.CLU        FALSE 
##  3  71.2 10.OPCML     TRUE  
##  4  29.0 2.HES6       FALSE 
##  5  31.8 3.TTR        FALSE 
##  6  38.8 4.BNIP3      FALSE 
##  7  68.1 5.CNTNAP2    TRUE  
##  8  63.2 6.AL118516.1 TRUE  
##  9  59.1 7.DLX6-AS1   TRUE  
## 10  68.6 8.MEF2C      TRUE  
## 11  89.3 9.GRIA4      TRUE  
## [1] "/Users/abel.vertesy/Documents/Analysis_R/Vignette_SU/Pc.cells.below.percent.ribo.of.0.1.cl.names.top.gene.res.0.2.png"
```

```
## 0.364 sec elapsed
```

![](scRNAseq_Vignette.2024.10.25.AV_files/figure-html/unnamed-chunk-25-1.png)<!-- -->![](scRNAseq_Vignette.2024.10.25.AV_files/figure-html/unnamed-chunk-25-2.png)<!-- -->

``` r
# knitr::include_graphics("Pc.cells.below.percent.ribo.of.0.1.cl.names.top.gene.res.0.2.png")
```

### `scBarplot.CellFractions()`

`scBarplot.CellFractions` generates a bar plot of cell fractions per cluster from a *Seurat* object, allowing downsampling, grouping by one variable, filling by another, custom color palettes, numerical value display on bars, and saving the plot.


``` r
# scBarplot.CellFractions(obj = combined.obj, group.by = ident2, fill.by = "Phase")
# knitr::include_graphics("Cell.proportions.of.Phase.by.cl.names.top.gene.res.0.2.downsampled.fr.barplot.png")
```

### `scBarplot.CellsPerCluster()`

`scBarplot.CellsPerCluster` generates a bar plot visualizing the fraction of cells within each cluster.


``` r
# ident2<- GetNamedClusteringRuns(combined.obj)[2]
scBarplot.CellsPerCluster(ident = ident2, sort = TRUE)
```

```
## 0# A tibble: 5 × 3
##   value names colour
##   <int> <chr> <chr> 
## 1   203 4     1     
## 2   228 3     2     
## 3   483 2     3     
## 4   758 1     4     
## 5   828 0     5     
## [1] "/Users/abel.vertesy/Documents/Analysis_R/Vignette_SU/Cells.per.Identity.Group.integrated_snn_res.0.1.integrated_snn_res.0.2.2500.c.bar.png"
```

```
## 0.305 sec elapsed
```

![](scRNAseq_Vignette.2024.10.25.AV_files/figure-html/unnamed-chunk-27-1.png)<!-- -->

### `plotClustSizeDistr()`

`plotClustSizeDistr` generates a bar plot or histogram to visualize the size distribution of clusters within a Seurat object, based on the specified clustering identity.


``` r
plotClustSizeDistr(obj = combined.obj, ident=identX, plot = TRUE, thr.hist = 30)
```

```
## 
##   0   1   2   3   4 
## 828 758 483 228 203 
## # A tibble: 5 × 3
##   value names colour
##   <int> <chr> <chr> 
## 1   828 0     1     
## 2   758 1     1     
## 3   483 2     1     
## 4   228 3     1     
## 5   203 4     1     
## [1] "/Users/abel.vertesy/Documents/Analysis_R/Vignette_SU/Cluster.sizes.at.integrated_snn_res.0.1.bar.png"
```

```
## 0.317 sec elapsed
```

![](scRNAseq_Vignette.2024.10.25.AV_files/figure-html/unnamed-chunk-28-1.png)<!-- -->
### `plotGeneExpHist()`

`plotGeneExpHist` creates and optionally saves a histogram displaying expression levels of specified genes within a Seurat object, with features for aggregate gene expression, expression threshold filtering, and quantile clipping for count data.

<!-- ```{r}
plotGeneExpHist(obj = combined.obj, genes = c("RPL15", "FTH1"))
``` -->

### `scBarplot.CellsPerObject()`

The `scBarplot.CellsPerObject` function visualizes the number of cells in each **Seurat** object within a list, displaying the distribution of cell counts across different datasets or experimental conditions.


``` r
scBarplot.CellsPerObject(ls.Seu)
```

```
## # A tibble: 3 × 3
##   value names colour
##   <int> <chr> <chr> 
## 1  2500 Exp1  1     
## 2   500 Exp2  1     
## 3   750 Exp3  1     
## [1] "/Users/abel.vertesy/Documents/Analysis_R/Vignette_SU/Nr.Cells.After.Filtering.bar.png"
```

```
## 0.287 sec elapsed
```

![](scRNAseq_Vignette.2024.10.25.AV_files/figure-html/unnamed-chunk-29-1.png)<!-- -->

----------------------------------------------------------------------------------------------------

## PCA

### `scPlotPCAvarExplained()`
Visualize the variance explained by PCA components, integrating seamlessly with `MarkdownReports` for documentation.


``` r
scPlotPCAvarExplained(combined.obj, plotname = "Variance Explained by Principal Components")
```

```
## # A tibble: 50 × 3
##    value names colour
##    <dbl> <chr> <lgl> 
##  1  9.33 1     TRUE  
##  2  6.51 2     TRUE  
##  3  5.87 3     TRUE  
##  4  3.91 4     TRUE  
##  5  3.62 5     TRUE  
##  6  3.22 6     TRUE  
##  7  3.01 7     TRUE  
##  8  2.68 8     TRUE  
##  9  2.59 9     TRUE  
## 10  2.29 10    TRUE  
## # ℹ 40 more rows
## [1] "/Users/abel.vertesy/Documents/Analysis_R/Vignette_SU/Variance.Explained.by.Principal.Components.bar.png"
```

```
## 0.448 sec elapsed
```

![](scRNAseq_Vignette.2024.10.25.AV_files/figure-html/unnamed-chunk-30-1.png)<!-- -->


### `scCalcPCAVarExplained()`

Calculate the variance explained by principal components, providing insights into the data structure and dimensionality reduction.
Used by `scPlotPCAvarExplained()`.


``` r
var_explained <- scCalcPCAVarExplained(combined.obj)
print(var_explained)
```

```
##        1        2        3        4        5        6        7        8 
## 9.332360 6.512810 5.872723 3.912707 3.615792 3.218808 3.007419 2.678719 
##        9       10       11       12       13       14       15       16 
## 2.585753 2.289343 2.249712 1.970911 1.947799 1.919924 1.816933 1.697605 
##       17       18       19       20       21       22       23       24 
## 1.614767 1.569853 1.537237 1.492872 1.492427 1.462241 1.410362 1.397870 
##       25       26       27       28       29       30       31       32 
## 1.394652 1.360058 1.339489 1.328720 1.308753 1.306202 1.302914 1.291201 
##       33       34       35       36       37       38       39       40 
## 1.285678 1.284520 1.279822 1.275696 1.274921 1.273328 1.270788 1.265798 
##       41       42       43       44       45       46       47       48 
## 1.264938 1.261772 1.259973 1.257276 1.255964 1.253753 1.252349 1.249852 
##       49       50 
## 1.247522 1.247114
```

----------------------------------------------------------------------------------------------------

## Scatter plots and Pair Plots

### `qFeatureScatter()`

`qFeatureScatter` generates a scatter plot comparing two features (genes or metrics) from a *Seurat* object with optional logarithmic transformations and saving capabilities, wrapping around *Seurat*'s FeatureScatter for enhanced usability.


``` r
qFeatureScatter(feature1 = "TOP2A", feature2 = "ID2", obj = combined.obj)
```

```
## [1] "/Users/abel.vertesy/Documents/Analysis_R/Vignette_SU/TOP2A.VS.ID2.png"
```

```
## 0.541 sec elapsed
```

![](scRNAseq_Vignette.2024.10.25.AV_files/figure-html/unnamed-chunk-32-1.png)<!-- -->

### `suPlotVariableFeatures()`

`suPlotVariableFeatures` generates a Variable Feature Plot for a specified *Seurat* object, labeling points with the top 20 variable genes, and saves the plot to a PDF file.


<!-- ```{r}
suPlotVariableFeatures(combined.obj)
``` -->

### `multiFeaturePlot.A4()`

The `multiFeaturePlot.A4` function saves multiple FeaturePlots, each representing a gene from a list of gene names, as JPEG images on A4 size paper.


``` r
combined.obj@version <- package_version("3.0.0") # Because of `SeuratObject::UpdateSeuratObject()`...
multiFeaturePlot.A4(
  list.of.genes,
  obj = combined.obj,
  foldername = substitute(list.of.genes),
  plot.reduction = "umap",
  intersectionAssay = c("RNA", "integrated")[1],
  layout = c("tall", "wide", FALSE)[2],
  colors = c("grey", "red"),
  nr.Col = nr.Col,
  nr.Row = nr.Row,
  raster = raster,
  cex = round(0.1/(nr.Col * nr.Row), digits = 2),
  cex.min = if (raster) TRUE else FALSE,
  gene.min.exp = "q01",
  gene.max.exp = "q99",
  subdir = TRUE,
  prefix = NULL,
  suffix = NULL,
  background_col = "white",
  aspect.ratio = c(FALSE, 0.6)[2],
  saveGeneList = FALSE,
  w = wA4,
  h = hA4,
  scaling = 1,
  format = c("jpg", "pdf", "png")[1]
)
```

```
## [1] "list.of.genes-umap/"
## [1] "ParentDir defined as:"
## [1] "OutDir defined as:"
## [1] "b.Subdirname defined as:"
## check.genes: 0.007 sec elapsed
## [1] "layout active, nr.Col ignored."
## [1] "1 MALAT1 TMSB4X MT-CO1"
```

```
## [1] "OutDir defined as:"
## 1.333 sec elapsed
```

``` r
combined.obj@version <- package_version("5.0.1")
```

<!-- ```{r echo=FALSE, out.width="100%"} -->
<!-- knitr::include_graphics( -->
<!--  # "C:\\images\\cells.0.2.4.png" -->
<!-- paste0(OutDir, "cells.0.2.4.png"), -->
<!--                         error=FALSE) -->
<!-- ``` -->

### `PlotTopGenesPerCluster()`

The `PlotTopGenesPerCluster` function visualizes the top N differentially expressed (DE) genes for each cluster within a specified clustering resolution of a **Seurat** object, aiding in the exploration of gene expression patterns across clusters.


``` r
PlotTopGenesPerCluster(obj = combined.obj, cl_res = 0.5, nrGenes = 2)
```

```
## [1] "TopGenes.umaps-umap/"
## [1] "ParentDir defined as:"
## [1] "OutDir defined as:"
## [1] "b.Subdirname defined as:"
## check.genes: 0.006 sec elapsed
## [1] "Old names contained duplicated elements 0"
## [1] "layout active, nr.Col ignored."
## [1] "1 HIST1H4C TUBA1B"
```

```
## [1] "OutDir defined as:"
## 0.939 sec elapsed
## [1] "TopGenes.umaps-umap/"
## [1] "ParentDir defined as:"
## [1] "OutDir defined as:"
## [1] "b.Subdirname defined as:"
## check.genes: 0.006 sec elapsed
## [1] "Old names contained duplicated elements 1"
## [1] "layout active, nr.Col ignored."
## [1] "1 UBE2C PTTG1"
```

```
## [1] "OutDir defined as:"
## 0.789 sec elapsed
## [1] "TopGenes.umaps-umap/"
## [1] "ParentDir defined as:"
## [1] "OutDir defined as:"
## [1] "b.Subdirname defined as:"
## check.genes: 0.006 sec elapsed
## [1] "Old names contained duplicated elements 10"
## [1] "layout active, nr.Col ignored."
## [1] "1 HS6ST3 RBFOX1"
```

```
## [1] "OutDir defined as:"
## 0.758 sec elapsed
## [1] "TopGenes.umaps-umap/"
## [1] "ParentDir defined as:"
## [1] "OutDir defined as:"
## [1] "b.Subdirname defined as:"
## check.genes: 0.006 sec elapsed
## [1] "Old names contained duplicated elements 11"
## [1] "layout active, nr.Col ignored."
## [1] "1 AL118516.1 POLR2A"
```

```
## [1] "OutDir defined as:"
## 0.842 sec elapsed
## [1] "TopGenes.umaps-umap/"
## [1] "ParentDir defined as:"
## [1] "OutDir defined as:"
## [1] "b.Subdirname defined as:"
## check.genes: 0.007 sec elapsed
## [1] "Old names contained duplicated elements 12"
## [1] "layout active, nr.Col ignored."
## [1] "1 DLX6-AS1 NRXN3"
```

```
## [1] "OutDir defined as:"
## 0.763 sec elapsed
## [1] "TopGenes.umaps-umap/"
## [1] "ParentDir defined as:"
## [1] "OutDir defined as:"
## [1] "b.Subdirname defined as:"
## check.genes: 0.007 sec elapsed
## [1] "Old names contained duplicated elements 13"
## [1] "layout active, nr.Col ignored."
## [1] "1 MEF2C NKAIN2"
```

```
## [1] "OutDir defined as:"
## 0.781 sec elapsed
## [1] "TopGenes.umaps-umap/"
## [1] "ParentDir defined as:"
## [1] "OutDir defined as:"
## [1] "b.Subdirname defined as:"
## check.genes: 0.006 sec elapsed
## [1] "Old names contained duplicated elements 14"
## [1] "layout active, nr.Col ignored."
## [1] "1 ERBB4 SCGN"
```

```
## [1] "OutDir defined as:"
## 0.773 sec elapsed
## [1] "TopGenes.umaps-umap/"
## [1] "ParentDir defined as:"
## [1] "OutDir defined as:"
## [1] "b.Subdirname defined as:"
## check.genes: 0.006 sec elapsed
## [1] "Old names contained duplicated elements 15"
## [1] "layout active, nr.Col ignored."
## [1] "1 GRIA4 CRABP1"
```

```
## [1] "OutDir defined as:"
## 0.846 sec elapsed
## [1] "TopGenes.umaps-umap/"
## [1] "ParentDir defined as:"
## [1] "OutDir defined as:"
## [1] "b.Subdirname defined as:"
## check.genes: 0.006 sec elapsed
## [1] "Old names contained duplicated elements 16"
## [1] "layout active, nr.Col ignored."
## [1] "1 OPCML KAZN"
```

```
## [1] "OutDir defined as:"
## 0.787 sec elapsed
## [1] "TopGenes.umaps-umap/"
## [1] "ParentDir defined as:"
## [1] "OutDir defined as:"
## [1] "b.Subdirname defined as:"
## check.genes: 0.007 sec elapsed
## [1] "Old names contained duplicated elements 2"
## [1] "layout active, nr.Col ignored."
## [1] "1 CRYAB CST3"
```

```
## [1] "OutDir defined as:"
## 0.751 sec elapsed
## [1] "TopGenes.umaps-umap/"
## [1] "ParentDir defined as:"
## [1] "OutDir defined as:"
## [1] "b.Subdirname defined as:"
## check.genes: 0.007 sec elapsed
## [1] "Old names contained duplicated elements 3"
## [1] "layout active, nr.Col ignored."
## [1] "1 MYC HES6"
```

```
## [1] "OutDir defined as:"
## 0.764 sec elapsed
## [1] "TopGenes.umaps-umap/"
## [1] "ParentDir defined as:"
## [1] "OutDir defined as:"
## [1] "b.Subdirname defined as:"
## check.genes: 0.007 sec elapsed
## [1] "Old names contained duplicated elements 4"
## [1] "layout active, nr.Col ignored."
## [1] "1 AC092957.1 APOE"
```

```
## [1] "OutDir defined as:"
## 0.815 sec elapsed
## [1] "TopGenes.umaps-umap/"
## [1] "ParentDir defined as:"
## [1] "OutDir defined as:"
## [1] "b.Subdirname defined as:"
## check.genes: 0.006 sec elapsed
## [1] "Old names contained duplicated elements 5"
## [1] "layout active, nr.Col ignored."
## [1] "1 NNAT AL589740.1"
```

```
## [1] "OutDir defined as:"
## 0.817 sec elapsed
## [1] "TopGenes.umaps-umap/"
## [1] "ParentDir defined as:"
## [1] "OutDir defined as:"
## [1] "b.Subdirname defined as:"
## check.genes: 0.007 sec elapsed
## [1] "Old names contained duplicated elements 6"
## [1] "layout active, nr.Col ignored."
## [1] "1 TTR COL3A1"
```

```
## [1] "OutDir defined as:"
## 0.906 sec elapsed
## [1] "TopGenes.umaps-umap/"
## [1] "ParentDir defined as:"
## [1] "OutDir defined as:"
## [1] "b.Subdirname defined as:"
## check.genes: 0.007 sec elapsed
## [1] "Old names contained duplicated elements 7"
## [1] "layout active, nr.Col ignored."
## [1] "1 MT1X BNIP3"
```

```
## [1] "OutDir defined as:"
## 0.97 sec elapsed
## [1] "TopGenes.umaps-umap/"
## [1] "ParentDir defined as:"
## [1] "OutDir defined as:"
## [1] "b.Subdirname defined as:"
## check.genes: 0.01 sec elapsed
## [1] "Old names contained duplicated elements 8"
## [1] "layout active, nr.Col ignored."
## [1] "1 CNTNAP2 KCNQ3"
```

```
## [1] "OutDir defined as:"
## 1.187 sec elapsed
## [1] "TopGenes.umaps-umap/"
## [1] "ParentDir defined as:"
## [1] "OutDir defined as:"
## [1] "b.Subdirname defined as:"
## check.genes: 0.01 sec elapsed
## [1] "Old names contained duplicated elements 9"
## [1] "layout active, nr.Col ignored."
## [1] "1 IGFBP2 PGK1"
```

```
## [1] "OutDir defined as:"
## 1.102 sec elapsed
```

<img src="../../../Documents/Analysis_R/Vignette_SU/PlotTopGenesPerCluster/DEG.markers.res.0.5.cluster.0.umap.1.HIST1H4C.TUBA1B.jpg" width="100%" />

<img src="../../../Documents/Analysis_R/Vignette_SU/PlotTopGenesPerCluster/DEG.markers.res.0.5.cluster.1.umap.1.UBE2C.PTTG1.jpg" width="100%" />

<!-- ```{r echo=FALSE, out.width="100%"} -->
<!-- knitr::include_graphics( -->
<!--  # "C:\\images\\PlotTopGenesPerCluster\\DEG.markers.res.0.5.cluster.2.umap.1.CRYAB.CST3.jpg" -->
<!-- paste0(OutDir, "PlotTopGenesPerCluster/"), -->
<!--                         error=FALSE) -->
<!-- ``` -->

<!-- ```{r echo=FALSE, out.width="100%"} -->
<!-- knitr::include_graphics( -->
<!--  # "C:\\images\\PlotTopGenesPerCluster\\DEG.markers.res.0.5.cluster.3.umap.1.MYC.HES6.jpg" -->
<!-- paste0(OutDir, "PlotTopGenesPerCluster/"), -->
<!--                         error=FALSE) -->
<!-- ``` -->

### `qQC.plots.BrainOrg()`

The `qQC.plots.BrainOrg` function generates and arranges UMAP plots for specified quality control (QC) features from a **Seurat** object on an A4 page, providing a quick overview of quality control metrics specific to brain organization data.


``` r
qQC.plots.BrainOrg(
  obj = combined.obj,
  QC.Features = c("nFeature_RNA", "percent.ribo", "percent.mito", "nCount_RNA"),
  nrow = 2,
  ncol = 2
)
```

```
## [1] "/Users/abel.vertesy/Documents/Analysis_R/Vignette_SU/QC.markers.4.UMAP_nFeature_RNA_percent.ribo_percent.mito_nCount_RNA.png"
```

<img src="../../../Documents/Analysis_R/Vignette_SU/QC.markers.4.UMAP_nFeature_RNA_percent.ribo_percent.mito_nCount_RNA.png" width="100%" />

### `qMarkerCheck.BrainOrg()`

The `qMarkerCheck.BrainOrg` function generates plots for a predefined or custom set of gene markers within brain organoids, facilitating the quick assessment of their expression across different cells or clusters.


``` r
qMarkerCheck.BrainOrg(combined.obj)
```

```
##                   dl-EN                   ul-EN        Immature neurons 
##                  "KAZN"                 "SATB2"                   "SLA" 
##            Interneurons            Interneurons            Interneurons 
##              "DLX6-AS1"                 "ERBB4"                  "SCGN" 
## Intermediate progenitor                 S-phase               G2M-phase 
##                 "EOMES"                 "TOP2A"                  "H4C3" 
##                     oRG               Astrocyte          Hypoxia/Stress 
##                  "HOPX"                 "S100B"                 "DDIT4" 
##          Choroid.Plexus             Low-Quality              Mesenchyme 
##                   "TTR"                "POLR2A"                   "DCN" 
##              Glycolytic 
##                  "PDK1" 
##                         [,1]      
## dl-EN                   "KAZN"    
## ul-EN                   "SATB2"   
## Immature neurons        "SLA"     
## Interneurons...4        "DLX6-AS1"
## Interneurons...5        "ERBB4"   
## Interneurons...6        "SCGN"    
## Intermediate progenitor "EOMES"   
## S-phase                 "TOP2A"   
## G2M-phase               "H4C3"    
## oRG                     "HOPX"    
## Astrocyte               "S100B"   
## Hypoxia/Stress          "DDIT4"   
## Choroid.Plexus          "TTR"     
## Low-Quality             "POLR2A"  
## Mesenchyme              "DCN"     
## Glycolytic              "PDK1"    
## [1] "Signature.Genes.Top16-umap/"
## [1] "ParentDir defined as:"
## [1] "OutDir defined as:"
## [1] "b.Subdirname defined as:"
## check.genes: 0.006 sec elapsed
## [1] "Old names contained duplicated elements Interneurons Interneurons"
## [1] "layout active, nr.Col ignored."
## [1] "1 KAZN SATB2 SLA DLX6-AS1 ERBB4 SCGN EOMES TOP2A"
```

```
## [1] "2 HOPX S100B DDIT4 TTR POLR2A DCN PDK1"
```

```
## [1] "OutDir defined as:"
## 4.049 sec elapsed
```

<img src="../../../Documents/Analysis_R/Vignette_SU/umap.1.KAZN.SATB2.SLA.DLX6-AS1.ERBB4.SCGN.EOMES.TOP2A.jpg" width="100%" />

<img src="../../../Documents/Analysis_R/Vignette_SU/umap.2.HOPX.S100B.DDIT4.TTR.POLR2A.DCN.PDK1.jpg" width="100%" />

### `PlotTopGenes()`

The `PlotTopGenes` function generates UMAP plots showcasing the highest expressed genes, saving the plots in a subfolder. Prior execution of `calc.q99.Expression.and.set.all.genes` is required for this function.


``` r
combined.obj <- calc.q99.Expression.and.set.all.genes(obj = combined.obj,
                                                      obj.version = 3  # You will not needed unless updating via `SeuratObject::UpdateSeuratObject()`.
                                                      )
```

```
## [1] "Calculating Gene Quantiles"
## [1] "59.6% or 15915 of 26690 genes have q99 expr. > 0 (in 25 cells)."
## [1] "/Users/abel.vertesy/Documents/Analysis_R/Vignette_SU/Gene.expression.in.the.99th.quantile.in.combined.obj.combined.obj.hist.pdf"
```

```
## calc.q99.Expression.and.set.all.genes: 0.21 sec elapsed
## 1.002 sec elapsed
```

![](scRNAseq_Vignette.2024.10.25.AV_files/figure-html/unnamed-chunk-42-1.png)<!-- -->

```
## [1] "Quantile 0.99 is now stored under obj@misc$all.genes and $ expr.q99  Please execute all.genes <- obj@misc$all.genes."
```

``` r
PlotTopGenes(obj = combined.obj, n = 4)
```

```
## [1] "Highest.Expressed.Genes-umap/"
## [1] "ParentDir defined as:"
## [1] "OutDir defined as:"
## [1] "b.Subdirname defined as:"
## check.genes: 0.007 sec elapsed
## [1] "layout active, nr.Col ignored."
## [1] "1 MALAT1 TMSB4X MT-CO1 FTH1"
```

```
## [1] "OutDir defined as:"
## 1.235 sec elapsed
```

<img src="../../../Documents/Analysis_R/Vignette_SU/umap.7.RPL7.RPL15.ACTG1.FTL.jpg" width="100%" />

<img src="../../../Documents/Analysis_R/Vignette_SU/umap.4.RPS18.ACTB.MT-CYB.NRXN3.jpg" width="100%" />

<!-- ```{r echo=FALSE, out.width="100%"} -->
<!-- knitr::include_graphics( -->
<!--  # "C:\\images\\PlotTopGenes\\umap.5.MT-ATP6.RPL13A.MT-ND4.RPL10.jpg" -->
<!-- paste0(OutDir, "umap.5.MT-ATP6.RPL13A.MT-ND4.RPL10.jpg"), -->
<!--                         error=FALSE) -->
<!-- ``` -->

<!-- ```{r echo=FALSE, out.width="100%"} -->
<!-- knitr::include_graphics( -->
<!--  # "C:\\images\\PlotTopGenes\\umap.6.RPS14.OOEP.RPLP1.RPL13.jpg" -->
<!-- paste0(OutDir, "umap.6.RPS14.OOEP.RPLP1.RPL13.jpg"), -->
<!--                         error=FALSE) -->
<!-- ``` -->

### `save2plots.A4()`

The `save2plots.A4` function arranges and saves two plots, such as UMAP plots or any other types of plots, side-by-side or one above the other on a single A4 page.


``` r
p1 <- ggplot(iris, aes(Sepal.Length, Sepal.Width, color = Species)) +
  geom_point()
p2 <- ggplot(iris, aes(Petal.Length, Petal.Width, color = Species)) +
  geom_point()
st.malo <- list(p1, p2)
save2plots.A4(plot_list = st.malo)
```

```
## [1] "Saved as: st.malo"
```

<img src="../../../Documents/Analysis_R/Vignette_SU/st.malo.png" width="100%" />

### `save4plots.A4()`

The `save4plots.A4` function arranges and saves four plots, such as UMAPs or any other visualizations, onto a single A4 page, facilitating a compact comparison of different visualizations or clustering results.


``` r
p1 <- ggplot(iris, aes(Sepal.Length, Sepal.Width, color = Species)) +
  geom_point()
p2 <- ggplot(mtcars, aes(mpg, disp, color = as.factor(cyl))) +
  geom_point()
p3 <- ggplot(mpg, aes(displ, hwy, color = class)) +
  geom_point()
p4 <- ggplot(diamonds, aes(carat, price, color = cut)) +
  geom_point()
nantes <- list(p1, p2, p3, p4)
save4plots.A4(plot_list = nantes)
```

```
## [1] "Saved as: nantes"
```

<img src="../../../Documents/Analysis_R/Vignette_SU/nantes.png" width="100%" />

### `qqSaveGridA4()`

The `qqSaveGridA4` function saves a grid of 2 or 4 **ggplot** objects onto an A4 page, enabling efficient visualization arrangements for analysis or presentation purposes.


``` r
p1 <- ggplot(iris, aes(Sepal.Length, Sepal.Width, color = Species)) +
  geom_point()
p2 <- ggplot(iris, aes(Petal.Length, Petal.Width, color = Species)) +
  geom_point()
pl <- list(p1,p2)
qqSaveGridA4(plotlist = pl, plots = 1:2, fname = "Fractions.per.Cl.png")
```

```
## [1] "2 plots found, 1 2 are saved."
```

```
## [1] "/Users/abel.vertesy/Documents/Analysis_R/Vignette_SU/Fractions.per.Cl.png"
```

<img src="../../../Documents/Analysis_R/Vignette_SU/Fractions.per.Cl.png" width="100%" />

----------------------------------------------------------------------------------------------------
## Violin Plot

`qSeuViolin` generates a violin plot for a specified feature in a Seurat object, allowing for the data to be split by a specified grouping variable, with support for customization options such as logarithmic scaling, custom titles, and more.


``` r
qSeuViolin(obj = combined.obj, feature = "nFeature_RNA", caption = "Test")
```

```
## [1] "/Users/abel.vertesy/Documents/Analysis_R/Vignette_SU/nFeature_RNA.by.cl.names.top.gene.res.0.2.logY.violin.png"
```

```
## 0.938 sec elapsed
```

![](scRNAseq_Vignette.2024.10.25.AV_files/figure-html/unnamed-chunk-51-1.png)<!-- -->


----------------------------------------------------------------------------------------------------
## Colors


### `getClusterColors()`

The `getClusterColors` function retrieves and, if desired, displays the color scheme linked to clusters in a Seurat object based on a specified identity column.


``` r
colors <- head(getClusterColors(ident = GetClusteringRuns(combined.obj)[1]), 10)
```

```
## c("integrated_snn_res.0.1", "integrated_snn_res.0.2", "integrated_snn_res.0.3", 
## "integrated_snn_res.0.4", "integrated_snn_res.0.5")
```

![](scRNAseq_Vignette.2024.10.25.AV_files/figure-html/unnamed-chunk-52-1.png)<!-- -->

### `SeuratColorVector()`

The `SeuratColorVector` function extracts and, if specified, displays the color scheme associated with cluster identities within a Seurat object, ensuring consistent color representation in visualizations.


``` r
head(SeuratColorVector(ident = "integrated_snn_res.0.2", plot.colors = TRUE), 4)
```

```
## [1] "integrated_snn_res.0.2"
## ident.vec
##   0   1   2   3   4   5   6   7 
## 619 390 300 298 254 229 202 208
```

![](scRNAseq_Vignette.2024.10.25.AV_files/figure-html/unnamed-chunk-53-1.png)<!-- -->

```
## [1] "#CD9600" "#00BE67" "#7CAE00" "#00A9FF"
```

### `getDiscretePaletteObj()`

The `getDiscretePaletteObj` function creates a discrete color palette for visualizing clusters in a Seurat object, adjusting the palette size based on a specified identity column to accommodate the number of unique clusters.


``` r
colors <- getDiscretePaletteObj(ident.used = "integrated_snn_res.0.1", obj = combined.obj)
print(colors)
```

```
##         0         2         3         4         1 
## "#AA0DFE" "#3283FE" "#85660D" "#782AB6" "#565656"
```

### `gg_color_hue()`

The `gg_color_hue` function produces a vector of colors mimicking the default color palette of ggplot2, facilitating the creation of color sets for custom plotting functions or other applications requiring a similar aesthetic.


``` r
print(gg_color_hue(5))
```

```
## [1] "#F8766D" "#A3A500" "#00BF7D" "#00B0F6" "#E76BF3"
```

### `DiscretePaletteSafe()`

The `DiscretePaletteSafe` function generates a discrete color palette excluding any NA values, making it suitable for visualizations requiring a fixed number of distinct and reproducible colors.


``` r
colors <- DiscretePaletteSafe(n = 10)
print(colors)
```

```
##  [1] "#AA0DFE" "#3283FE" "#85660D" "#782AB6" "#565656" "#1C8356" "#16FF32"
##  [8] "#F7E1A0" "#E2E2E2" "#1CBE4F"
```



----------------------------------------------------------------------------------------------------
## 3D plots

### `plot3D.umap()`

The `plot3D.umap` function plots a 3D UMAP (Uniform Manifold Approximation and Projection) based on one of the metadata columns of a **Seurat** object. It utilizes Plotly for interactive visualization.


``` r
# plot3D.umap(combined.obj, category = "Phase")
```



### `plot3D.umap.gene()`

`plot3D.umap.gene` plots a three-dimensional UMAP with gene expression using the *Plotly* library.


``` r
# plot3D.umap.gene(obj = combined.obj, gene = "TOP2A")
```


----------------------------------------------------------------------------------------------------
## Miscellaneous

### `GetClusteringRuns()`

The `GetClusteringRuns` function retrieves metadata column names associated with clustering runs.


``` r
head(getClusterColors(obj = combined.obj, ident = GetClusteringRuns(combined.obj)[1]), 4)
```

```
## c("integrated_snn_res.0.1", "integrated_snn_res.0.2", "integrated_snn_res.0.3", 
## "integrated_snn_res.0.4", "integrated_snn_res.0.5")
```

![](scRNAseq_Vignette.2024.10.25.AV_files/figure-html/unnamed-chunk-59-1.png)<!-- -->

```
##         0         0         2         3 
## "#0000FF" "#0000FF" "#00FF00" "#000033"
```

