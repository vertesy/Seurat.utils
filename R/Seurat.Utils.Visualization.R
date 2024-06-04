# ____________________________________________________________________
# Seurat.Utils.Visualization.R ----
# ____________________________________________________________________
# source("~/GitHub/Packages/Seurat.utils/R/Seurat.Utils.Visualization.R")
# devtools::load_all("~/GitHub/Packages/Seurat.utils")
# devtools::document("~/GitHub/Packages/Seurat.utils"); devtools::load_all("~/GitHub/Packages/Seurat.utils")


# _________________________________________________________________________________________________
#' @title Plot filtering thresholds and distributions
#'
#' @description This function plots the filtering thresholds and distributions for Seurat objects,
#' using four panels to highlight the relationship between gene- and UMI-counts, and the
#' ribosomal- and mitochondrial-content.
#'
#' @param ls.obj A list of Seurat objects to be analyzed. Default: `ls.Seurat`.
#' @param parentdir The parent directory where the plots will be stored. Default: `OutDirOrig`.
#' @param suffices Suffixes used in the output plot file names. Default: Names of the Seurat objects in `ls.obj`.
#' @param filetype The file type of the output plot images. Default: `'.png'`.
#' @param below.mito Lower bound of mitochondrial content. Default: `p$thr.lp.mito`.
#' @param above.mito Upper bound of mitochondrial content. Default: `p$thr.hp.mito`.
#' @param below.ribo Lower bound of ribosomal content. Default: `p$thr.lp.ribo`.
#' @param above.ribo Upper bound of ribosomal content. Default: `p$thr.hp.ribo`.
#' @param below.nFeature_RNA Lower bound of RNA features. Default: `p$thr.lp.nFeature_RNA`.
#' @param above.nFeature_RNA Upper bound of RNA features. Default: `p$thr.hp.nFeature_RNA`.
#' @param subdir Subdirectory within `parentdir` where plots will be stored. Default is generated using a call to `kpp()`.
#' @param transparency Point transparency on scatter plots. Default: 0.25.
#' @param cex Size of points on scatter plots. Default: 0.75.
#' @param theme.used A `ggplot2` theme for all plots. Default: `theme_bw(base_size = 18)`.
#' @param LabelDistFromTop Distance from top for label placement. Default: 200.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   PlotFilters(ls.Seurat)
#' }
#' }
#'
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link[ggplot2]{labs}},
#' \code{\link[ggplot2]{geom_point}}
#'
#' @importFrom ggplot2 ggplot ggtitle geom_point
#' @importFrom Stringendo percentage_formatter
#' @importFrom MarkdownHelpers llprint
#' @importFrom cowplot plot_grid save_plot
#'
#' @export
PlotFilters <- function(
    ls.obj = ls.Seurat,
    parentdir = OutDirOrig,
    suffices = names(ls.obj),
    filetype = ".png",
    below.mito = p$"thr.lp.mito",
    above.mito = p$"thr.hp.mito",
    below.ribo = p$"thr.lp.ribo",
    above.ribo = p$"thr.hp.ribo",
    below.nFeature_RNA = if ("quantile.thr.lp.nFeature_RNA" %in% names(p)) p$"quantile.thr.lp.nFeature_RNA" else p$"thr.lp.nFeature_RNA",
    above.nFeature_RNA = p$"thr.hp.nFeature_RNA",
    subdir = Stringendo::FixPlotName(
      "Filtering.plots",
      "mito", p$"thr.hp.mito", p$"thr.lp.mito",
      "ribo", p$"thr.hp.ribo", p$"thr.lp.ribo",
      "nFeature", p$"thr.hp.nFeature_RNA", below.nFeature_RNA
    ),
    transparency = 0.25,
    cex = 0.75,
    theme.used = theme_bw(base_size = 18),
    LabelDistFromTop = 200 # for barplot_label
    ) {
  stopif(is.null(below.nFeature_RNA))
  MarkdownHelpers::llprint(
    "We filtered for high quality cells based on the number of genes detected [", above.nFeature_RNA, ";", below.nFeature_RNA,
    "] and the fraction of mitochondrial [", Stringendo::percentage_formatter(above.mito), ";", Stringendo::percentage_formatter(below.mito),
    "] and ribosomal [", Stringendo::percentage_formatter(above.ribo), ";", Stringendo::percentage_formatter(below.ribo), "] reads."
  )

  theme_set(theme.used)
  OutDir <- Stringendo::FixPath(parentdir, subdir)
  # print(OutDir)
  print(subdir)
  # stop()
  MarkdownReports::create_set_OutDir(OutDir)
  stopifnot(length(suffices) == length(ls.obj))

  Calculate_nFeature_LowPass <- if (below.nFeature_RNA < 1) below.nFeature_RNA else FALSE
  for (i in 1:length(ls.obj)) {
    print(suffices[i])
    mm <- ls.obj[[i]]@meta.data

    if (Calculate_nFeature_LowPass < 1) {
      below.nFeature_RNA <- floor(quantile(ls.obj[[i]]$"nFeature_RNA", probs = Calculate_nFeature_LowPass))
      iprint("below.nFeature_RNA at", percentage_formatter(Calculate_nFeature_LowPass), "percentile:", below.nFeature_RNA)
    }

    AllMetaColumnsPresent <- all(c("nFeature_RNA", "percent.mito", "percent.ribo") %in% colnames(mm))
    if (!AllMetaColumnsPresent) {
      print(c("nFeature_RNA", "percent.mito", "percent.ribo"))
      print(c("nFeature_RNA", "percent.mito", "percent.ribo") %in% colnames(mm))
      print("Try to run:")
      print('objX <- addMetaFraction(obj = objX, col.name = "percent.mito", gene.symbol.pattern =  "^MT\\.|^MT-")')
      print('objX <- addMetaFraction(obj = objX, col.name = "percent.ribo", gene.symbol.pattern =  "^RPL|^RPS")')
      stop()
    }

    filt.nFeature_RNA <- (mm$"nFeature_RNA" < below.nFeature_RNA & mm$"nFeature_RNA" > above.nFeature_RNA)
    filt.below.mito <- (mm$"percent.mito" < below.mito & mm$"percent.mito" > above.mito)
    filt.below.ribo <- (mm$"percent.ribo" < below.ribo & mm$"percent.ribo" > above.ribo)

    mm <- cbind(mm, filt.nFeature_RNA, filt.below.mito, filt.below.ribo)

    mm$colour.thr.nFeature <- cut(mm$"nFeature_RNA",
      breaks = c(-Inf, above.nFeature_RNA, below.nFeature_RNA, Inf),
      labels = c(
        paste0("LQ (<", above.nFeature_RNA, ")"),
        paste0("HQ (", above.nFeature_RNA, "< X <", below.nFeature_RNA, ")"),
        paste0("Dbl/Outlier (>", below.nFeature_RNA, ")")
      )
    )

    A <- ggplot(data = mm, aes(x = nFeature_RNA, fill = colour.thr.nFeature)) +
      geom_histogram(binwidth = 100) +
      ggtitle(paste("Cells between", above.nFeature_RNA, "and", below.nFeature_RNA, " UMIs are selected (", pc_TRUE(filt.nFeature_RNA), ")")) +
      geom_vline(xintercept = below.nFeature_RNA) +
      geom_vline(xintercept = above.nFeature_RNA)
    # A

    B <- ggplot2::ggplot(mm, aes(x = nFeature_RNA, y = percent.mito)) +
      ggplot2::ggtitle(paste(
        "Cells below", Stringendo::percentage_formatter(below.mito),
        "mito reads are selected (with A:", pc_TRUE(filt.nFeature_RNA & filt.below.mito), ")"
      )) +
      ggplot2::geom_point(
        alpha = transparency, size = cex, show.legend = FALSE,
        aes(color = filt.nFeature_RNA & filt.below.mito)
      ) +
      scale_x_log10() + # scale_y_log10() +
      # annotation_logticks() +
      geom_hline(yintercept = below.mito) +
      geom_hline(yintercept = above.mito) +
      geom_vline(xintercept = below.nFeature_RNA) +
      geom_vline(xintercept = above.nFeature_RNA)
    # B

    C <- ggplot(mm, aes(x = nFeature_RNA, y = percent.ribo)) +
      ggtitle(paste(
        "Cells below", Stringendo::percentage_formatter(below.ribo),
        "ribo reads are selected (with A:",
        pc_TRUE(filt.nFeature_RNA & filt.below.ribo), ")"
      )) +
      geom_point(
        alpha = transparency, size = cex, show.legend = FALSE,
        aes(color = filt.nFeature_RNA & filt.below.ribo)
      ) +
      scale_x_log10() + # scale_y_log10() +
      annotation_logticks() +
      geom_hline(yintercept = below.ribo) +
      geom_hline(yintercept = above.ribo) +
      geom_vline(xintercept = below.nFeature_RNA) +
      geom_vline(xintercept = above.nFeature_RNA)
    # C

    D <- ggplot(mm, aes(x = percent.ribo, y = percent.mito)) +
      ggtitle(paste(
        "Cells w/o extremes selected (with A,B,C:",
        pc_TRUE(filt.nFeature_RNA & filt.below.mito & filt.below.ribo), ")"
      )) +
      geom_point(
        alpha = transparency, size = cex, show.legend = FALSE,
        aes(color = filt.nFeature_RNA & filt.below.mito & filt.below.ribo)
      ) +
      scale_x_log10() +
      scale_y_log10() +
      annotation_logticks() +
      geom_hline(yintercept = below.mito) +
      geom_hline(yintercept = above.mito) +
      geom_vline(xintercept = below.ribo) +
      geom_vline(xintercept = above.ribo)
    # D

    plot_list <- list(A, B, C, D)
    px <- cowplot::plot_grid(plotlist = plot_list, nrow = 2, ncol = 2, labels = LETTERS[1:4])
    fname <- kpps(OutDir, FixPlotName("Filtering.thresholds", suffices[i], filetype))
    # print(fname)
    cowplot::save_plot(filename = fname, plot = px, base_height = 12, ncol = 1, nrow = 1) # Figure 2
    stopifnot(file.exists(fname))
  } # for
  # _________________________________________________________________________________________________
  create_set_Original_OutDir()
}


# _________________________________________________________________________________________________
# plotting.statistics.and.QC.R ______________________________ ----
# _________________________________________________________________________________________________
# source('~/GitHub/Packages/Seurat.utils/Functions/Plotting.statistics.and.QC.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Plotting.statistics.and.QC.R"))

# Source: self + web

# Requirements __________________________________________
# require(Seurat)
# require(ggplot2)
# tools for tools::toTitleCase

# May also require
# try (source('/GitHub/Packages/CodeAndRoll/CodeAndRoll.R'),silent= FALSE) # generic utilities funtions
# require('MarkdownReports') # require("devtools")

# _________________________________________________________________________________________________
#' @title Calculate the percent of variation explained by individual PC's
#'
#' @description This function calculates the percentage of variation each principal component (PC)
#' accounts for in a Seurat object. It's specifically tailored for Seurat objects and provides a
#' convenient way to understand the variance distribution across PCs. For similar calculations on
#' standard PCA objects, refer to github.com/vertesy/Rocinante `PCA.percent.var.explained()`.
#'
#' @param obj A Seurat object from which to calculate the percentage of variation explained by each
#' PC. Default: `combined.obj`.
#'
#' @return A named vector with the percentage of variation explained by each principal component.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   data("combined.obj") # Example Seurat object
#'   var_explained <- scCalcPCAVarExplained(combined.obj)
#'   print(var_explained)
#' }
#' }
#'
#' @export
scCalcPCAVarExplained <- function(obj = combined.obj) { # Determine percent of variation associated with each PC.
  pct <- obj@reductions$pca@stdev / sum(obj@reductions$pca@stdev) * 100
  names(pct) <- 1:length(obj@reductions$pca@stdev)
  return(pct)
}

# _________________________________________________________________________________________________
#' @title Plot the percent of variation explained by individual PC's
#'
#' @description This function plots the percentage of variation explained by each principal
#' component (PC) in a Seurat object. It allows for a visual assessment of how much variance is
#' captured by each PC, which is crucial for dimensionality reduction analyses. Users can choose
#' between two plotting methods: one using `MarkdownReports` and the other using `ggExpress`.
#'
#' @param obj A Seurat object from which to plot the percentage of variation explained by each PC.
#' Default: `combined.obj`.
#' @param plotname The title of the plot to be generated. Default: "Variance Explained by Principal
#' Components".
#' @param sub Subtitle for the plot, typically including information about the number of cells and
#' features analyzed. Default: A string generated from `obj` stating the number of cells and
#' features.
#' @param use.MarkdownReports Boolean indicating whether to use `MarkdownReports` for plotting.
#' If `FALSE`, `ggExpress` is used. Default: `FALSE`.
#'
#' @return Generates a plot showing the percent of variation each PC accounts for. This function
#' does not return a value but instead generates a plot directly.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   data("combined.obj") # Example Seurat object
#'   scPlotPCAvarExplained(combined.obj, use.MarkdownReports = TRUE)
#' }
#' }
#'
#' @importFrom MarkdownReports wbarplot
#' @importFrom ggExpress qbarplot
#'
#' @export
scPlotPCAvarExplained <- function(obj = combined.obj,
                                  plotname = "Variance Explained by Principal Components",
                                  sub = paste(ncol(obj), "cells, ", nrow(obj), "features."),
                                  use.MarkdownReports = FALSE,
                                  # caption = .parseKeyParams(obj, suffix = "| hline at 1%"),
                                  caption = "hline at 1%",
                                  ...) {

  pct <- scCalcPCAVarExplained(obj)
  if (use.MarkdownReports) {
    MarkdownReports::wbarplot(pct, xlab = "Principal Components", ylab = "% of variation explained")
    barplot_label(round(pct, digits = 2), barplotted_variable = pct, cex = .5)
  } else {
    ggExpress::qbarplot(
      vec = pct, plotname = plotname, subtitle = sub,
      xlab = "Principal Components", ylab = "% of variation explained",
      w = 10, h = 5, hline = 1, caption = caption
    )
  }
}



# _________________________________________________________________________________________________
# Gene Expression Plots ______________________________ ----
# _________________________________________________________________________________________________


# _________________________________________________________________________________________________
#' @title Gene Expression as Fraction of Total UMI Counts
#'
#' @description This function computes and visualizes gene expression levels as a fraction of total
#' UMI (Unique Molecular Identifier) counts across all genes in a Seurat object. It aims to highlight
#' the relative contribution of the most highly expressed genes to the overall transcriptome.
#'
#' @param obj A Seurat object containing gene expression data.
#' @param n.genes.barplot The number of top genes to be displayed in the final barplot, showing
#' their expression as a percentage of the total UMIs. Default is 25.
#' @param width.barplot The width of the barplot that visualizes the highest expressed genes.
#' Default is a quarter of `n.genes.barplot`.
#'
#' @return The same Seurat object passed as input, but with an additional list in the `@misc` slot
#' named `'TotalReadFraction'` that contains the relative total expression of each gene as a
#' fraction of total UMIs.
#'
#' @examples
#' \dontrun{
#' combined.obj <- PercentInTranscriptome(combined.obj)
#' }
#'
#' @export
PercentInTranscriptome <- function(
    obj = combined.obj,
    assay = DefaultAssay(obj),
    n.genes.barplot = 25,
    width.barplot = round(n.genes.barplot / 4),
    ...) {
  #
  message("Obj. version: ", obj@version)
  message("assay: ", assay)

  m.expr <- if (obj@version < "5") {
    obj@assays$RNA@counts
  } else {
    SeuratObject::GetAssayData(object = obj, assay = assay)
  }


  total.Expr <- sort(rowSums(m.expr), decreasing = TRUE)
  relative.total.Expr <- total.Expr / sum(total.Expr)
  print(head(iround(100 * relative.total.Expr), n = n.genes.barplot))

  Relative.of.Total.Gene.Expression <- relative.total.Expr * 100

  qhistogram(Relative.of.Total.Gene.Expression,
    logX = FALSE, logY = TRUE,
    plotname = "Gene expression as fraction of all UMI's",
    subtitle = "Percentage in RNA-counts",
    xlab = "Percent in Transcriptome (total per gene)",
    ylab = "Number of genes",
    xlab.angle = 45,
    w = 7, h = 5,
    ...)

  Highest.Expressed.Genes <- head(iround(100 * relative.total.Expr), n = n.genes.barplot)
  qbarplot(Highest.Expressed.Genes,
    plotname = "Percentage of highest expressed genes",
    subtitle = "Total, in RNA-counts",
    xlab = "",
    ylab = "Gene expression as percent of all UMI's",
    xlab.angle = 45,
    w = width.barplot,  h = 5,
    ...)

  message("!!! \nTotalReadFraction is now stored under combined.obj@misc$'TotalReadFraction'.")

  obj@misc$"TotalReadFraction" <- relative.total.Expr
  return(obj)
}



# _________________________________________________________________________________________________
#' @title Histogram All Genes' Expression Level and a Highlighted Gene
#'
#' @description Shows a comparison of the expression level of the chose gene to all genes.
#' Very useful to see if the gene has a meaningful expression level. This function generates a
#' histogram to visualize the expression level distribution of a specified gene across all cells in
#' a Seurat object. It highlights the position of the gene of interest within the overall distribution.
#'
#' @param gene The gene of interest for which the expression level distribution is to be plotted.
#' Default: 'TOP2A'.
#' @param obj A Seurat object containing the expression data. Default: The first Seurat object in `ls.Seurat`.
#' @param assay The assay from which to retrieve the expression data. Default: "RNA".
#' @param slot The slot in the Seurat object from which to retrieve the expression data. Options
#' include "counts" for raw counts and "data" for normalized (and possibly log-transformed) data.
#' Default: "data".
#' @param w The width of the plot. Default: 7.
#' @param h The height of the plot. Default: 4.
#' @param ... Any other parameter that can be passed to the internally called functions.
#'
#' @export
plotGeneExpressionInBackgroundHist <- function(
    gene = "TOP2A",
    obj = combined.obj,
    assay = "RNA",
    slot = c("counts", "data")[2],
    w = 7, h = 4,
    ...) {
  message("gene: ", gene)
  stopifnot(gene %in% rownames(obj))


  GEX.Counts <- GetAssayData(object = obj, assay = assay, slot = slot)

  GEX.Counts.total <- rowSums(GEX.Counts)
  genes.expression <- GEX.Counts.total[gene]
  mean.expr <- iround(mean(GEX.Counts[gene, ]))

  suffx <- if (slot == "counts") "raw" else "normalised, logtransformed"
  (pname <- paste(gene, "and the", suffx, "transcript count distribution"))

  ggExpress::qhistogram(GEX.Counts.total,
    vline = genes.expression, logX = TRUE, w = w, h = h,
    subtitle = paste("It belong to the top", pc_TRUE(GEX.Counts.total > genes.expression), "of genes (black line). Mean expr:", mean.expr),
    plotname = pname, xlab = "Total Transcripts in Dataset", ylab = "Number of Genes",
    ...
  )
}




# _________________________________________________________________________________________________
#' @title Histogram of Gene / Geneset Aggregate Expression Across Cells
#'
#' @description Creates and optionally saves a histogram showing expression levels of specified genes
#' within a Seurat object. Provides options for aggregate gene expression, expression threshold filtering,
#' and quantile clipping for count data.
#'
#' @param obj Seurat object to analyze; Default: `combined.obj`.
#' @param genes Vector of gene names to include in the analysis; Default: c("MALAT1", "MT-CO1").
#' @param assay Assay to use from the Seurat object; Default: "RNA".
#' @param slot_ Data slot to use ('data' or 'counts'); Default: "data".
#' @param thr_expr Expression threshold for highlighting in the plot; Default: 10.
#' @param suffix Additional text to append to the plot title; Default: NULL.
#' @param xlab Label for the x-axis; Default: "log10(Summed UMI count @data)".
#' @param return_cells_passing If TRUE, returns count of cells exceeding the expression threshold; Default: TRUE.
#' @param quantile_thr Quantile threshold for clipping count data; Default: 0.95.
#' @param return_quantile If TRUE, returns cell count exceeding the quantile threshold; Default: FALSE.
#' @param w Width of the plot in inches; Default: 9.
#' @param h Height of the plot in inches; Default: 5.
#' @param show_plot If TRUE, displays the generated plot; Default: TRUE.
#' @param ... Additional arguments for customization.
#'
#' @return Depending on the parameters, can return a ggplot object, the number of cells passing
#' the expression threshold, or the number of cells exceeding the quantile threshold.
#'
#' @examples
#' \dontrun{
#' plotGeneExprHistAcrossCells(obj = yourSeuratObject, genes = c("GeneA", "GeneB"))
#' }
#'
#' @return Depending on the parameters, can return a ggplot object, the number of cells passing
#' the expression threshold, or the number of cells exceeding the quantile threshold.
#'
#' @export
#' @importFrom scales hue_pal
#' @importFrom Seurat GetAssayData
#' @importFrom ggplot2 geom_vline labs
#' @importFrom ggExpress qhistogram
plotGeneExprHistAcrossCells <- function(
    genes = c("MALAT1", "MT-CO1", "MT-CO2", "MT-CYB", "TMSB4X", "KAZN"),
    obj = combined.obj,
    assay = "RNA", slot_ = "data",
    thr_expr = 10,
    suffix = NULL,
    xlab = paste0("Expression -log10(Summed UMI count @", slot_, ")"),
    return_cells_passing = TRUE,
    quantile_thr = 0.95,
    return_quantile,
    w = 9, h = 5,
    show_plot = TRUE,
    ...) {
  stopifnot(
    length(genes) > 0,
    slot_ %in% c("data", "counts")
  )

  # browser()
  # Aggregate genes if necessary
  aggregate <- length(genes) > 1
  SummedExpressionPerCell <- colSums(GetAssayData(object = obj, assay = assay,
                                                  slot = slot_)[genes, , drop = F])
  head(SummedExpressionPerCell)

  # Clip counts if necessary
  if (slot_ == "counts") {
    SummedExpressionPerCell <- CodeAndRoll2::clip.at.fixed.value(
      distribution = SummedExpressionPerCell,
      thr = quantile(SummedExpressionPerCell, probs = .95)
    )
  }

  # Create annotation
  CPT <- paste("slot:", slot_, "| assay:", assay, "| cutoff at", iround(thr_expr))

  # Add a subtitle with the number of genes and the expression threshold
  SUBT <- filter_HP(SummedExpressionPerCell, threshold = thr_expr, return_conclusion = TRUE, plot.hist = FALSE)
  if (aggregate) {
    SUBT <- paste(SUBT, "\n", length(genes), "genes summed up, \n e.g:", kppc(head(genes)))
    TTL <- paste("Summed Gene-set Expression -", suffix)
  } else {
    TTL <- trimws(paste("Gene Expression -", paste(genes), suffix))
  }

  # Create the plot
  pobj <- ggExpress::qhistogram(SummedExpressionPerCell,
    plotname = TTL,
    subtitle = SUBT,
    caption = CPT,
    suffix = suffix,
    vline = thr_expr[1], filtercol = -1,
    xlab = xlab,
    ylab = "# of cells",
    w = w, h = h,
    ...
  )

  # draw additional vlines if needed
  if (length(thr_expr) > 1) {
    pobj <- pobj +
      ggplot2::geom_vline(xintercept = thr_expr[-1], col = 2, lty = 2, lwd = 1) +
      ggplot2::labs(caption = "Red line marks original estimate")
    ggExpress::qqSave(ggobj = pobj, title = sppp(TTL, "w.orig")) # , ext = '.png'
  }


  # Print the plot
  if (show_plot) print(pobj)

  # Return the number of cells passing the filter
  if (return_cells_passing) {
    return(MarkdownHelpers::filter_HP(SummedExpressionPerCell, threshold = thr_expr, plot.hist = FALSE))
  }
}


# _________________________________________________________________________________________________
#' @title Proportion of Cells Expressing Given Genes
#'
#' @description Calculates the proportion of cells expressing one or more specified genes.
#'
#' @param genes Character vector of gene names of interest.
#' @param group.by Optional grouping variable for analysis (e.g., cell type). Default: 'all'.
#' @param obj Seurat object to analyze. Default: `combined.obj`.
#'
#' @return Data frame with genes and their cell expression proportion, optionally grouped.
#'
#' @examples
#' \dontrun{
#' PrctCellExpringGene(genes = c("LTB", "GNLY"), obj = combined.obj)
#' }
#'
#' @source Adapted from Ryan-Zhu on GitHub.
#'
#' @export
PrctCellExpringGene <- function(genes, group.by = "all", obj = combined.obj) {
  stopifnot("Some genes not foun!." = all(genes %in% row.names(obj)))

  if (group.by == "all") {
    prct <- 1:length(genes)
    for (i in seq(prct)) prct[i] <- ww.calc_helper(genes = genes[1], obj = obj)
    result <- data.frame("Markers" = genes, "Cell_proportion" = prct)
    return(result)
  } else {
    list <- Seurat::SplitObject(object = obj, split.by = group.by)
    factors <- names(list)
    results <- lapply(list, PrctCellExpringGene, genes = genes)
    for (i in 1:length(factors)) {
      results[[i]]$Feature <- factors[i]
    }
    combined <- do.call("rbind", results)
    return(combined)
  }
}


# _________________________________________________________________________________________________
#' @title Helper to calculate Cell Expression Proportion for Gene
#'
#' @description Computes the proportion of cells expressing a specific gene within a Seurat object.
#'
#' @param obj Seurat object with cell data.
#' @param genes Single gene name as a character string.
#' @param slot Slot to use for the analysis. Default: 'RNA'.
#'
#' @return Proportion of cells expressing the gene. Returns `NA` if the gene is not found.
#'
#' @examples
#' \dontrun{
#' ww.calc_helper(obj = seurat_object, genes = "Gene1")
#' }
#'
#' @source Adapted from Ryan-Zhu on GitHub.
#'
#' @export
ww.calc_helper <- function(obj, genes, slot = "RNA") {
  # stopifnot("Some genes not found!." = all(genes %in% row.names(obj)))
  counts <- obj[[slot]]@counts
  ncells <- ncol(counts)
  if (genes %in% row.names(counts)) {
    sum(counts[genes, ] > 0) / ncells
  } else {
    return(NA)
  }
}



# _________________________________________________________________________________________________
# Barplots / Compositional analysis ______________________________ ----
# _________________________________________________________________________________________________


# _________________________________________________________________________________________________

#' @title Generate Barplot of Cell Fractions
#'
#' @description This function generates a bar plot of cell fractions per cluster from a Seurat object.
#' It offers the option to downsample data, equalizing the number of cells in each group
#' to the number in the smallest group. The plot's bars are grouped by one variable and filled by another.
#' The function supports custom color palettes, drawing numerical values on bars, and saving the plot.
#'
#' @param group.by The variable to group by for the bar plot.
#' @param fill.by The variable to fill by for the bar plot.
#' @param downsample Logical indicating whether to downsample data to equalize group sizes.
#' @param min.nr.sampled.cells The minimal number of cells to sample from each identity class. Defaults to 200 cells.
#' @param dsample.to.repl.thr Logical indicating if sampling should be done with replacement. Defaults to FALSE.
#' @param plotname The title of the plot.
#' @param suffix Optional suffix for the plot title.
#' @param sub_title Optional subtitle for the plot.
#' @param hlines Numeric vector specifying y-intercepts of horizontal lines to add to the plot.
#' @param return_table Logical; if TRUE, returns a contingency table instead of plotting.
#' @param save_plot Logical; if TRUE, saves the generated plot.
#' @param seedNr Seed for random number generation to ensure reproducibility.
#' @param w Width of the plot in inches.
#' @param h Height of the plot in inches.
#' @param draw_plot Logical; if FALSE, suppresses plotting (useful if only the table is desired).
#' @param show_numbers Logical; if TRUE, adds count numbers on top of each bar in the plot.
#' @param min_frequency Minimum fraction to display individually in the plot; smaller fractions
#' are aggregated into an "Other" category.
#' @param custom_col_palette Specifies whether to use a standard or custom color palette.
#' @param color_scale Defines the color scale to use for the plot if a custom palette is selected.
#' @param also.pdf Save plot in both png and pdf formats.
#' @param min.pct Show % Labels above this threshold. Default = 0.05, or above 5 pct.
#' @param cex.pct Font size of pct labels.
#' @param show.total.cells Show total cells
#' @param cex.total Label size for total cells
#' @param ... Additional parameters passed to internally called functions.
#'
#' @return Depending on the value of `return_table`, either returns a ggplot object or a list
#' containing values and percentages tables.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   scBarplot.CellFractions(obj = combined.obj, group.by = "integrated_snn_res.0.1", fill.by = "Phase", downsample = TRUE)
#'   scBarplot.CellFractions(obj = combined.obj, group.by = "integrated_snn_res.0.1", fill.by = "Phase", downsample = FALSE)
#' }
#' }
#' @seealso \code{\link[tools]{toTitleCase}}, \code{\link[ggplot2]{ggplot}}, \code{\link[dplyr]{group_by}}, \code{\link[dplyr]{summarise}}
#' @importFrom tools toTitleCase
#' @importFrom dplyr group_by summarise sample_n
#' @importFrom ggplot2 ggplot geom_bar geom_hline labs theme_classic theme element_text scale_fill_manual geom_text
#'
#' @export
scBarplot.CellFractions <- function(
    fill.by,
    group.by = GetNamedClusteringRuns()[1],
    obj = combined.obj,
    downsample = FALSE,
    min.nr.sampled.cells = 200,
    plotname = kpp(toTitleCase(fill.by), "proportions.by", group.by),
    suffix = NULL,
    prefix = NULL,
    sub_title = suffix,
    hlines = c(.25, .5, .75),
    return_table = FALSE,
    save_plot = TRUE,
    also.pdf = FALSE,
    seedNr = 1989,
    w = 10, h = ceiling(0.7 * w),
    draw_plot = TRUE,
    show_numbers = FALSE,
    min.pct = 0.05,
    cex.pct = 2.5,
    min_frequency = 0, # 0.025,
    custom_col_palette = FALSE,
    color_scale = getDiscretePaletteObj(ident.used = fill.by, obj = obj, palette.used = "glasbey"),
    rnd_colors = FALSE,
    show.total.cells = TRUE,
    cex.total = 2,
    xlab.angle = 45,
    ...) {
  # Input assertions
  stopifnot(
    inherits(obj, "Seurat"), # obj must be a Seurat object
    is.numeric(min_frequency) && length(min_frequency) == 1 && min_frequency >= 0 && min_frequency < 1, # min_frequency must be between 0 and 1
    group.by %in% colnames(obj@meta.data), # group.by must be a valid column in the meta.data slot of the Seurat object
    fill.by %in% colnames(obj@meta.data) # fill.by must be a valid column in the meta.data slot of the Seurat object
  )

  META <- obj@meta.data

  set.seed(seedNr)
  pname.suffix <- capt.suffix <- NULL

  if (downsample) {
    tbl_X <- table(META[[fill.by]])
    n_smallest_group <- min(tbl_X)
    largest_grp <- max(tbl_X)

    message("The size of the smallest group is: ", n_smallest_group, " cells.")

    dsample.to.repl.thr <- (n_smallest_group < min.nr.sampled.cells) # if less than 200 cells are sampled, sample with replacement
    if (dsample.to.repl.thr) {
      message(paste(
        "If smallest category is <", min.nr.sampled.cells,
        "of total cells, than down- or up-sampling, with replacement to that minimum."
      ))
    }

    # Update plot name and caption to reflect downsampling
    plotname <- kpp(plotname, "downsampled")
    pname.suffix <- "(downsampled)"

    capt.suffix <- paste0(
      "\nDownsampled all groups in ", fill.by, " (Y) to ", min.nr.sampled.cells,
      " cells before splitting by X. \nThis number is max(smallest group, 5% of total cells). Largest groups previosly was: ", largest_grp
    )
  }


  # Construct the caption based on downsampling and minimum frequency
  PFX <- if (show_numbers) "Numbers denote # cells." else percentage_formatter(min.pct, prefix = "Labeled above")
  caption_ <- paste("Top: Total cells per bar. |", PFX, capt.suffix)

  if (min_frequency > 0) caption_ <- paste(caption_, "\nCategories <", percentage_formatter(min_frequency), "are shown together as 'Other'")
  pname_ <- paste(plotname, pname.suffix)


  # Create a contingency table of the data
  contingency.table <- table(META[, group.by], META[, fill.by])
  print(contingency.table)


  if (show.total.cells) {
    # First, calculate the total counts per group
    totals <- META %>%
      group_by(!!sym(group.by)) %>%
      summarise(Total = n()) %>%
      ungroup()

    # Merge totals back with the original data for labeling
    group_by_column <- group.by
    META <- META %>%
      left_join(totals, by = setNames(nm = group_by_column, group_by_column))
  }


  if (draw_plot) {
    # calculate the proportions and add up small fractions
    prop_table <- META %>%
      group_by(!!as.name(fill.by)) %>%
      summarise(proportion = n() / nrow(META)) %>%
      mutate("category" = ifelse(proportion < min_frequency, "Other", as.character(!!as.name(fill.by))))

    categories <- unique(prop_table$"category")
    n.categories <- length(categories)
    message(n.categories, " Y-categories present: ", kppc(sort(categories)))

    # join the proportions back to the original data
    META <- left_join(META, prop_table, by = fill.by)

    subtt <- kppws(group.by, "|", ncol(obj), "cells", sub_title)


    if (downsample) {
      # Downsample the data
      META <-
        META %>%
        group_by(!!sym(fill.by)) %>%
        sample_n(
          size = max(n_smallest_group, min.nr.sampled.cells),
          replace = dsample.to.repl.thr
        ) %>%
        ungroup()

      contingency.table <- table(META[[group.by]], META[[fill.by]])
      contingency.table <- addmargins(contingency.table)
      print(contingency.table)
    }

    # Plot the data
    pl <- META %>%
      group_by(!!sym(group.by)) %>%
      ggplot(aes(fill = category, x = !!sym(group.by))) +
      geom_hline(yintercept = hlines, lwd = 1.5) +
      geom_bar(position = "fill") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(
        title = pname_, subtitle = subtt,
        x = group.by, y = "Fraction of Cells",
        fill = fill.by, caption = caption_
      ) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
      scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.1))) +
      theme(axis.text.x = element_text(angle = xlab.angle, hjust = 1))


    # Apply custom color palette if specified
    if (custom_col_palette) {
      palette_x <- color_scale[seq(categories)]
      message("palette: ", kppc(palette_x))
      pl <- pl + scale_fill_manual(values = palette_x)
    } else if (rnd_colors) {
      colzz <- sample(rainbow(n.categories))
      pl <- pl + scale_fill_manual(values = colzz)
    }

    if (show_numbers) {
      pl <- pl + geom_text(aes(label = ..count..),
        stat = "count", position = position_fill(vjust = 0.5)
      )
    } else {
      pl <- pl + geom_text(
        aes(label = ifelse((..count.. / tapply(..count.., ..x.., sum)[..x..]) >= min.pct,
          scales::percent(..count.. / tapply(..count.., ..x.., sum)[..x..], accuracy = 1),
          ""
        )),
        stat = "count", position = position_fill(vjust = 0.5),
        size = cex.pct
      )
    }

    if (show.total.cells) {
      pl <- pl + geom_text(
        data = totals, aes(x = !!sym(group.by), y = 1, label = Total),
        vjust = -0.5, size = cex.total, inherit.aes = FALSE
      )
    }

    if (save_plot) {
      sfx <- shorten_clustering_names(group.by)
      if (!is.null(suffix)) sfx <- sppp(sfx, suffix)
      if (min_frequency) sfx <- sppp(sfx, min_frequency)
      qqSave(
        ggobj = pl, title = plotname, also.pdf = also.pdf, w = w, h = h,
        suffix = sppp(sfx, "fr.barplot"), ...
      )
    } # save_plot
  } # draw_plot

  # Return contingency table or plot based on return_table flag
  if (return_table) {
    ls.tables <- list(
      "values" = contingency.table,
      "percentages" = CodeAndRoll2::rowDivide(mat = contingency.table, vec = rowSums(contingency.table))
    )
    return(ls.tables)
  } else {
    # if(show_plot)
    return(pl)
  } # else barplot
}





# _________________________________________________________________________________________________
#' @title Barplot of Fraction of Cells per Cluster
#'
#' @description Visualizes the fraction of cells within each cluster through a barplot.
#'
#' @param obj Seurat object for analysis. Default: `combined.obj`.
#' @param ident Cluster identity. Used to specify which clustering results to visualize.
#' Default: First entry from ordered clustering runs.
#' @param sort If TRUE, sorts clusters by size. Default: FALSE.
#' @param title Title for the plot. Default: "Cells per Identity Group".
#' @param sub Subtitle for the plot. Default: "identity".
#' @param label If TRUE, shows cell count or percentage based on the label vector. Default: TRUE.
#' @param palette Color palette for the barplot. Default: 'glasbey'.
#' @param return_table If TRUE, returns the data used for plotting instead of the plot itself. Default: FALSE.
#' @param min.cells Minimum cell count threshold for categories. Adjusted by data size.
#' @param suffix Optional suffix for file naming. Used in conjunction with `kpp`.
#' @param ylab_adj Adjustment factor for y-axis label positioning. Default: 1.1.
#' @param ... Additional parameters for internal function calls.
#'
#' @examples
#' \dontrun{
#' scBarplot.CellsPerCluster()
#' scBarplot.CellsPerCluster(sort = TRUE)
#' }
#' @export scBarplot.CellsPerCluster
#'
#' @importFrom ggExpress qbarplot

scBarplot.CellsPerCluster <- function(
    obj = combined.obj,
    ident = GetOrderedClusteringRuns(obj = obj)[1],
    sort = FALSE,
    title = "Cells per Identity Group",
    sub = ident,
    label = list(TRUE, "percent")[[1]],
    suffix = if (label == "percent") "percent" else NULL,
    palette = c("alphabet", "alphabet2", "glasbey", "polychrome", "stepped")[3],
    return_table = FALSE,
    ylab_adj = 1.1,
    min.cells = round(ncol(obj) / 100),
    ...) {
  stopifnot(ident %in% colnames(obj@meta.data))

  cell.per.cl <- obj[[ident]][, 1]
  cell.per.cluster <- (table(cell.per.cl))
  if (sort) cell.per.cluster <- sort(cell.per.cluster)
  lbl <- if (isFALSE(label)) {
    NULL
  } else if (label == "percent") {
    percentage_formatter(cell.per.cluster / sum(cell.per.cluster))
  } else if (isTRUE(label)) {
    cell.per.cluster
  } else {
    label
  }

  imessage("min cell thr:", min.cells)
  n.clusters <- length(cell.per.cluster)
  nr.cells.per.cl <- table(obj[[ident]][, 1])
  SBT <- pc_TRUE(nr.cells.per.cl < min.cells,
    NumberAndPC = TRUE,
    suffix = paste("of identites are < 1% (below min.cells:", min.cells, ")")
  )

  pl <- ggExpress::qbarplot(cell.per.cluster,
    plotname = title,
    subtitle = paste0(sub, "\n", SBT),
    suffix = kpp(ident, suffix),
    col = 1:n.clusters,
    xlab.angle = 45,
    ylim = c(0, ylab_adj * max(cell.per.cluster)),
    label = lbl,
    ylab = "Cells",
    palette_use = DiscretePaletteSafe(n = n.clusters, palette.used = palette),
    ...
  )

  if (return_table) {
    return(cell.per.cluster)
  } else {
    return(pl)
  }
}

# _________________________________________________________________________________________________
#' @title Cluster Size Distribution Plot (Barplot or Histogram)
#'
#' @description Generates a bar plot or histogram to visualize the size distribution of clusters
#' within a Seurat object, based on the specified clustering identity.
#'
#' @param obj Seurat object for analysis. Default: `combined.obj`.
#' @param ident Clustering identity to base the plot on.
#' Default: The second entry from `GetClusteringRuns()`.
#' @param plot Whether to display the plot (TRUE) or return cluster sizes (FALSE). Default: TRUE.
#' @param thr.hist Threshold for switching from a bar plot to a histogram based on the number of
#' clusters. Default: 30.
#' @param ... Extra parameters for the plot.
#'
#' @examples
#' \dontrun{
#' plotClustSizeDistr()
#' }
#'
#' @importFrom ggExpress qbarplot qhistogram
#'
#' @export
plotClustSizeDistr <- function(
    obj = combined.obj, ident,
    plot = TRUE, thr.hist = 30, ...) {
  stopifnot(ident %in% colnames(obj@meta.data))

  clust.size.distr <- table(obj@meta.data[, ident])
  print(clust.size.distr)
  resX <- gsub(pattern = ".*res\\.", replacement = "", x = ident)
  ptitle <- paste("Cluster sizes at ", ident)
  psubtitle <- paste(
    "Nr.clusters:", length(clust.size.distr),
    "| median size:", median(clust.size.distr),
    "| CV:", Stringendo::percentage_formatter(cv(clust.size.distr))
  )
  xlb <- "Cluster size (cells)"
  ylb <- "Nr of Clusters"
  xlim <- c(0, max(clust.size.distr))

  if (plot) {
    if (length(clust.size.distr) < thr.hist) {
      ggExpress::qbarplot(clust.size.distr,
        plotname = ptitle, subtitle = psubtitle,
        label = clust.size.distr, xlab = "Clusters", ylab = xlb, ...
      )
    } else {
      ggExpress::qhistogram(
        vec = clust.size.distr, plotname = ptitle, subtitle = psubtitle,
        xlab = xlb, ylab = ylb, xlim = xlim, ...
      )
    }
  } else {
    "return vector"
    clust.size.distr
  }
}


# _________________________________________________________________________________________________
#' @title Barplot the Fraction of Cells Above Threshold per Cluster
#'
#' @description Generates a bar plot depicting the percentage of cells within each cluster that
#' exceed a specified threshold, based on a selected metadata column.
#'
#' @param value.col Column in metadata with values to assess against `thrX`. Default: 'percent.ribo'.
#' @param thrX Threshold for calculating the fraction of cells. Default: 0.3.
#' @param obj Seurat object with single-cell data. Default: `combined.obj`.
#' @param id.col Cluster identity column in metadata. Default: 'cl.names.top.gene.res.0.3'.
#' @param return.df Whether to return the underlying data frame instead of the plot. Default: FALSE.
#' @param label Whether to add labels to the bar plot. Default: NULL.
#' @param subtitle Optional subtitle for the plot.
#' @param suffix Suffix for the output file name.
#' @param above Whether to calculate the fraction of cells above or below the threshold. Default: TRUE.
#' @param ... Additional parameters for plotting functions.
#'
#' @examples
#' \dontrun{
#' scBarplot.FractionAboveThr(id.col = "cl.names.top.gene.res.0.3", value.col = "percent.ribo", thrX = 0.2)
#' }
#'
#' @seealso \code{\link[dplyr]{select}}, \code{\link[dplyr]{group_by}}
#'
#' @importFrom dplyr select group_by summarize
#'
#' @export
scBarplot.FractionAboveThr <- function(
    value.col = "percent.ribo",
    thrX = 0.1,
    obj = combined.obj,
    id.col = GetClusteringRuns(obj)[1],
    subtitle = id.col,
    return.df = FALSE,
    label = NULL,
    suffix = NULL,
    above = TRUE,
    ylim = c(0, 100), # set to null for relative y axis
    ...) {
  stopifnot(value.col %in% colnames(obj@meta.data))

  meta <- obj@meta.data
  metacol <- meta %>%
    dplyr::select(c(id.col, value.col))

  (df_cells_above <- metacol %>%
    dplyr::group_by(!!sym(id.col)) %>%
    summarize(
      n_cells = n(),
      n_cells_above = sum(!!sym(value.col) > thrX),
      fr_n_cells_above = n_cells_above / n_cells,
      fr_n_cells_below = 1 - fr_n_cells_above
    )
  )


  pass <-
    if (above) {
      metacol[, value.col] > thrX
    } else {
      metacol[, value.col] < thrX
    }
  total_average <- iround(100 * mean(pass))

  df_2vec <-
    if (above) {
      df_cells_above[, c(1, 4)]
    } else {
      df_cells_above[, c(1, 5)]
    }


  (v.fr_n_cells_above <- 100 * deframe(df_2vec))

  tag <- if (above) "above" else "below"
  if (is.null(label)) label <- percentage_formatter(deframe(df_2vec), digitz = 2)

  pname <- paste("Pc. cells", tag, value.col, "of", thrX)
  ggobj <- ggExpress::qbarplot(v.fr_n_cells_above,
    label = label,
    plotname = pname,
    filename = FixPlotName(kpp(pname, id.col, ".pdf")),
    suffix = suffix,
    subtitle = subtitle,
    caption = paste(
      "Overall average (black line):", iround(total_average), "% |",
      substitute(obj)
    ),
    xlab.angle = 45,
    xlab = "Clusters",
    ylab = paste("% Cells", tag, "thr. (", value.col, ")"),
    ylim = ylim,
    hline = total_average,
    ...
  )
  print(ggobj)
  if (return.df) {
    return(df_cells_above)
  } else {
    ggobj
  }
}


# _________________________________________________________________________________________________
#' @title Fraction of Cells Below Threshold per Cluster
#'
#' @description Generates a bar plot to visualize the percentage of cells within each cluster that
#' fall below a specified threshold, according to a metadata column value.
#' Inherits all parameters from `scBarplot.FractionAboveThr` with the exception that `above` is set to FALSE.
#'
#' @param thrX Threshold value for assessing cell counts. Default: 0.01.
#' @param value.col Metadata column with values for threshold comparison. Default: 'percent.ribo'.
#' @param id.col Cluster identifier in metadata. Default: 'cl.names.top.gene.res.0.3'.
#' @param obj Seurat object with cell data. Default: `combined.obj`.
#' @param return.df If TRUE, returns the data frame instead of the plot. Default: FALSE.
#'
#' @examples
#' \dontrun{
#' scBarplot.FractionBelowThr(id.col = "cl.names.top.gene.res.0.3", value.col = "percent.ribo", thrX = 0.01)
#' }
#'
#' @seealso `scBarplot.FractionAboveThr`
#' @seealso \code{\link[dplyr]{select}}, \code{\link[dplyr]{group_by}}
#'
#' @importFrom dplyr select group_by summarize
#'
#' @export
scBarplot.FractionBelowThr <- function(
    thrX = 0.2,
    value.col = "percent.ribo",
    id.col = "cl.names.top.gene.res.0.3",
    obj = combined.obj,
    return.df = FALSE,
    subtitle = id.col,
    suffix = NULL,
    ...) {
  scBarplot.FractionAboveThr(
    thrX = thrX,
    value.col = value.col,
    id.col = id.col,
    obj = obj,
    return.df = return.df,
    subtitle = subtitle,
    suffix = suffix,
    above = FALSE # Set `above` argument to FALSE to get fraction below threshold
  )
}

# _________________________________________________________________________________________________
# Pie Charts / Compositional analysis ______________________________ ----
# _________________________________________________________________________________________________
#' @title scPieClusterDistribution
#'
#' @description This function generates a pie chart of cluster distributions for a given clustering
#' identity in a single-cell RNA-seq object.
#'
#' @param obj A single-cell RNA-seq object. Default: `combined.obj`.
#' @param ident A character string specifying the clustering identity to be used. Default: the first
#' clustering run in the object.
#' @param ... Additional arguments passed to other methods.
#'
#' @return A pie chart displaying the cluster distribution.
#' @importFrom ggplot2 ggplot aes geom_bar coord_polar theme_minimal
#' @importFrom scales percent
#'
#' @examples
#' \dontrun{
#' scPieClusterDistribution(obj = combined.obj, ident = 'cluster_identity')
#' }
scPieClusterDistribution <- function(obj = combined.obj, ident = GetClusteringRuns(obj)[1],
                                     ...) {

  # Input assertions
  stopifnot(
    is(obj, "Seurat"), is.character(ident), length(ident) == 1,
    ident %in% colnames(obj@meta.data)
  )

  # Compute cluster sizes
  cluster_sizes <- table(obj[[ident]])
  print(cluster_sizes)

  # Create pie chart
  qpie(cluster_sizes, caption = .parseBasicObjStats(obj), subtitle = ident)
}




# _________________________________________________________________________________________________
# List of Seurat Objects ______________________________ ----
# _________________________________________________________________________________________________






# _________________________________________________________________________________________________
#' @title Barplot of Cells Per Seurat Object
#'
#' @description Visualizes the number of cells in each Seurat object within a list, showing the
#' distribution of cell counts across different datasets or experimental conditions.
#'
#' @param ls.obj List of Seurat objects to analyze. Default: `ls.Seurat`.
#' @param plotname Title for the plot. Default: 'Nr.Cells.After.Filtering'.
#' @param xlab.angle Angle for x-axis labels, enhancing readability. Default: 45.
#' @param names Optionally provide custom names for x-axis labels. If FALSE, uses object names
#' from `ls.obj`. Default: FALSE.
#' @param ... Extra parameters passed to `qbarplot`.
#'
#' @examples
#' \dontrun{
#' ls.obj <- list(obj, obj)
#' scBarplot.CellsPerObject(ls.obj)
#' }
#'
#' @export scBarplot.CellsPerObject

scBarplot.CellsPerObject <- function(
    ls.obj = ls.Seurat,
    plotname = "Nr.Cells.After.Filtering",
    xlab.angle = 45,
    names = FALSE, ...) {
  stopifnot(
    "Should be run on a list of Seu. objects" = length(ls.obj) > 1,
    "Should be run on a list of Seu. objects" = all(sapply(ls.obj, is, "Seurat"))
  )

  cellCounts <- unlapply(ls.obj, ncol)
  names(cellCounts) <- if (length(names) == length(ls.obj)) names else names(ls.obj)
  ggExpress::qbarplot(cellCounts,
    plotname = plotname,
    subtitle = paste(sum(cellCounts), "cells in total"),
    label = cellCounts,
    xlab.angle = xlab.angle,
    ylab = "Cells",
    ...
  )
}



# _________________________________________________________________________________________________
#' @title Stacked Barplot of Metadata Categories for List of Seurat Objects
#'
#' @description Creates and saves a stacked barplot for a specified metadata category
#' from a list of Seurat objects.
#'
#' @param ls.obj List of Seurat objects.
#' @param meta.col The metadata column name to be used for the barplot.
#' @return A ggplot object representing the stacked barplot.
#'
#' @examples
#' \dontrun{
#' ls.obj <- list(obj, obj)
#' scBarplotStackedMetaCateg_List(ls.obj, meta.col = orig.ident)
#' }
#'
#' @importFrom ggExpress qbarplot.df
#' @importFrom dplyr group_by summarise select
#'
#' @export
scBarplotStackedMetaCateg_List <- function(
    ls.obj,
    meta.col, # e.g. orig.ident
    ...) {
  stopifnot(
    "Should be run on a list of Seu. objects" = length(ls.obj) > 1,
    "Should be run on a list of Seu. objects" = all(sapply(ls.obj, is, "Seurat")),
    length(meta.col) == 1,
    "meta.col not found in 1st object" = meta.col %in% colnames(ls.obj[[1]]@meta.data)
  )

  # Creating a data frame for the stacked bar plot
  df <- do.call(rbind, lapply(seq_along(ls.obj), function(x) {
    data.frame(
      Sample = names(ls.obj)[x],
      Category = ls.obj[[x]]@meta.data[[meta.col]],
      stringsAsFactors = FALSE
    )
  }))

  # Summarizing to get counts of cells per category for each sample
  df <- df %>%
    dplyr::group_by(Sample, Category) %>%
    dplyr::summarise(Cells = n(), .groups = "drop") %>%
    dplyr::select(Sample, Cells, Category)

  TTL <- paste(meta.col, "per object")
  p <- ggExpress::qbarplot.df(df,
    plotname = TTL,
    scale = TRUE, hide.legend = FALSE,
    ...
  )
  print(p)
  return(df)
}


# _________________________________________________________________________________________________
# Colors ______________________________ ----
# _________________________________________________________________________________________________
#' @title Reproduce the ggplot2 default color palette
#'
#' @description Generates a vector of colors that emulates the default color palette used by ggplot2.
#' This function is useful for creating color sets for custom plotting functions or for applications
#' outside of ggplot2 where a similar aesthetic is desired.
#'
#' @param n Integer specifying the number of distinct colors to generate.
#'
#' @return A character vector of color values in hexadecimal format, replicating the default hue
#' color scale of ggplot2.
#'
#' @examples
#' \dontrun{
#' # Generate a palette of 5 colors
#' print(gg_color_hue(5))
#' }
#'
#' @export
gg_color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}



# _________________________________________________________________________________________________
#' @title Safely generate a discrete color palette (NA).
#'
#' @description Safe wrapper around Seurat's DiscretePalette(), which returns NA's if too many
#' categories are requested
#' @param ident.used The identity column used for determining the number of clusters, Default: GetClusteringRuns()[1]
#' @param obj Seurat object, Default: combined.obj
#' @param palette.used The name of the palette to use, Default: c("alphabet", "alphabet2",
#' "glasbey", "polychrome", "stepped")[1]
#' @param show.colors Whether to display the colors in the palette, Default: FALSE
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   getDiscretePalette()
#' }
#' }
#' @importFrom MarkdownHelpers color_check
#' @importFrom gplots rich.colors
#'
#' @export
getDiscretePalette <- function() .Deprecated("DiscretePaletteSafe and DiscretePaletteObj")


# _________________________________________________________________________________________________
#' @title Generate a Discrete Color Palette for Seurat Clusters
#'
#' @description Generates a discrete color palette for visualizing clusters in a Seurat object,
#' using a specified identity column to determine the number of unique clusters.
#'
#' @param ident.used Identity column in the Seurat object to base the color palette on.
#' @param obj Seurat object containing clustering information.
#' @param palette.used The palette name to use for color generation. Options include "alphabet",
#' "alphabet2", "glasbey", "polychrome", and "stepped". Default: "alphabet2".
#' @param show.colors If TRUE, displays the generated colors. Default: FALSE.
#' @param seed Seed for random color generation, ensuring reproducibility. Default: 1989.
#'
#' @return A character vector of color values corresponding to the number of clusters.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   ident.used <- "resolution_1"
#'   obj <- YourSeuratObject
#'   colors <- getDiscretePaletteObj(ident.used = ident.used, obj = obj)
#'   print(colors)
#' }
#' }
#'
#' @export
getDiscretePaletteObj <- function(ident.used,
                                  obj,
                                  palette.used = c("alphabet", "alphabet2", "glasbey", "polychrome", "stepped")[2],
                                  show.colors = FALSE,
                                  seed = 1989) {
  stopifnot(
    is.character(ident.used), is(obj, "Seurat"),
    is.character(palette.used), is.logical(show.colors), is.numeric(seed)
  )

  n.clusters <- CodeAndRoll2::nr.unique(obj[[ident.used]])
  # browser()
  colorz <- DiscretePaletteSafe(
    n = n.clusters,
    palette.used = palette.used,
    show.colors = show.colors,
    seed = seed
  )

  return(colorz)
}


# _________________________________________________________________________________________________
#' @title Safely generate a Discrete color palette.
#'
#' @description Generates a discrete color palette, ensuring no NA values are included, suitable
#' for visualizations where a specific number of distinct, reproducible colors is needed.
#'
#' @param n Number of colors to generate.
#' @param palette.used Palette name to use for generating colors. Options include "alphabet",
#' "alphabet2", "glasbey", "polychrome", "stepped". Default: "alphabet2".
#' @param show.colors If TRUE, displays the generated color palette. Default: FALSE.
#' @param seed Seed value for reproducibility, especially when random color generation is involved.
#' Default: 1989.
#'
#' @return Character vector of HEX color values.
#'
#' @examples
#' \dontrun{
#' colors <- DiscretePaletteSafe(n = 10)
#' print(colors)
#' }
#'
#' @importFrom gplots rich.colors
#' @importFrom Seurat DiscretePalette
#'
#' @export
DiscretePaletteSafe <- function(n,
                                palette.used = c("alphabet", "alphabet2", "glasbey", "polychrome", "stepped")[2],
                                show.colors = FALSE,
                                seed = 1989) {
  stopifnot(
    is.numeric(n), n > 0, is.character(palette.used),
    is.logical(show.colors), is.numeric(seed)
  )

  colorz <- Seurat::DiscretePalette(n = n, palette = palette.used)

  if (anyNA(colorz)) {
    colorsOK <- colorz[!is.na(colorz)]
    n.colz <- length(colorsOK)

    msg <- paste(
      "More categories then present in the palette", n, "vs.",
      n.colz, "in", palette.used, "-> recycling."
    )
    warning(msg, immediate. = TRUE)

    set.seed(seed)
    if (n > 10 * n.colz) {
      colorz <- sample(gplots::rich.colors(n))
    } else {
      colorz <- sample(x = colorsOK, size = n, replace = TRUE)
    }

    stopifnot(!anyNA(colorz))
  }

  if (show.colors) MarkdownHelpers::color_check(colorz)
  return(colorz)
}


# _________________________________________________________________________________________________
#' @title Regenerate Cluster Colors from a Seurat Object
#'
#' @description Regenerate and optionally displays the color scheme associated with the clusters
#' in a Seurat object as defined by a specified identity column.
#'
#' @param obj Seurat object containing clustering information.
#' @param use_new_palettes Logical indicating whether to use custom palettes defined in
#' `DiscretePalette` function. Default: TRUE.
#' @param palette Name of the color palette to use if `use_new_palettes` is TRUE.
#' Options: "alphabet", "alphabet2", "glasbey", "polychrome", "stepped". Default: "glasbey".
#' @param ident Clustering identity to use for coloring. Retrieved from the first entry
#' of `GetClusteringRuns()` by default.
#' @param show If TRUE, displays a plot showing the color mapping for each cluster. Default: TRUE.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   getClusterColors(obj = combined.obj, ident = GetClusteringRuns(combined.obj)[1])
#' }
#' }
#' @export
#'
#' @importFrom scales hue_pal
getClusterColors <- function(
    obj = combined.obj,
    use_new_palettes = TRUE,
    palette = c("alphabet", "alphabet2", "glasbey", "polychrome", "stepped")[3],
    ident = GetClusteringRuns()[1],
    show = TRUE) {
  (identities <- levels(as.factor(obj[[ident]][, 1])))
  n.clusters <- length(unique(obj[[ident]][, 1]))
  color_palette <- if (use_new_palettes) {
    DiscretePalette(n = n.clusters, palette = palette)
  } else {
    scales::hue_pal()(length(identities))
  }

  names(color_palette) <- (identities)
  identvec <- obj[[ident]][, 1]
  colz <- color_palette[identvec]
  names(colz) <- identvec
  if (show) MarkdownHelpers::color_check(unique(colz))
  colz
}



# _________________________________________________________________________________________________
#' @title Regenerate Color Scheme for Clusters in Seurat Object as a vector
#'
#' @description Extracts and optionally displays the color scheme assigned to cluster identities
#' within a Seurat object, facilitating consistent color usage across visualizations. You can
#' check results in a barplot with `MarkdownHelpers::color_check()`.
#'
#' @param ident Specific clustering identity to use for color extraction.
#' If NULL, the active identity in `obj` is used. Default: NULL.
#' @param obj Seurat object from which to extract cluster colors.
#' Default: `combined.obj`.
#' @param plot.colors If TRUE, visually displays the color scheme.
#' Default: FALSE.
#' @param simple If TRUE, returns only the unique set of colors used.
#' If FALSE, returns a named vector mapping cluster identities to colors.
#' Default: FALSE.
#'
#' @examples
#' \dontrun{
#' # Display colors for the active identity
#' SeuratColorVector()
#' # Retrieve and plot colors for a specified clustering identity
#' SeuratColorVector(ident = "RNA_snn_res.1", plot.colors = TRUE)
#' }
#'
#' @seealso \code{\link[scales]{hue_pal}}
#'
#' @export
#' @importFrom scales hue_pal
SeuratColorVector <- function(ident = NULL, obj = combined.obj, plot.colors = FALSE, simple = FALSE) {
  if (!is.null(ident)) {
    print(ident)
    ident.vec <- obj[[ident]][, 1]
  } else {
    print(obj@active.ident)
    ident.vec <- obj@active.ident
  }
  ident.vec <- as.factor(ident.vec)
  print(table(ident.vec))
  colorlevels <- scales::hue_pal()(length(levels(ident.vec)))
  if (plot.colors) color_check(colorlevels)
  if (simple) {
    colorlevels
  } else {
    translate(
      vec = as.character(ident.vec),
      oldvalues = levels(ident.vec),
      newvalues = colorlevels
    )
  }
}


# _________________________________________________________________________________________________
# Metadata Heatmaps ______________________________ ----
# _________________________________________________________________________________________________


#' @title Plot and Save Heatmaps from Metadata Calculation Results
#'
#' @description Generates and saves heatmap visualizations for each metric in the results obtained
#' from metadata calculations, such as  `calculateAverageMetaData() - mean or median values of
#' specified features across different categories.
#'
#' @param results A list containing data frames with calculated metrics for each specified
#'   metadata feature, grouped by categories. Typically, this is the output from a
#'   function like `calculateAverageMetaData()`.
#' @param path The directory path where the heatmap images will be saved.
#'   Defaults to the current working directory (`getwd()`).
#' @param file.prefix A prefix for the filenames of the saved heatmap images.
#'   Defaults to "heatmap_".
#' @param scale Character indicating if the values should be scaled in the row direction,
#'   column direction, both ('row', 'column', 'none'). Defaults to "column".
#' @param cluster_rows Logical indicating whether to cluster rows. Defaults to FALSE.
#' @param show_rownames Logical indicating whether to show row names. Defaults to TRUE.
#' @param show_colnames Logical indicating whether to show column names. Defaults to TRUE.
#' @param ... Additional arguments passed to `pheatmap::pheatmap`.
#'
#' @details This function loops through each metric in the `results`, creates a heatmap
#' for it using `pheatmap`, and saves the heatmap as a PNG file in the specified path.
#' The file names will start with the provided `file.prefix`, followed by the metric name.
#'
#' @examples
#' # Assuming `results` is the output from `calculateAverageMetaData`:
#' plotAndSaveHeatmaps(results, path = "path/to/save/heatmaps", file.prefix = "myData_")
#'
#' @return Invisible. The function primarily generates and saves files without returning data.
#'
#' @export
plotAndSaveHeatmaps <- function(results, path = getwd(),
                                file.prefix = "heatmap_",
                                scale = "column",
                                cluster_rows = FALSE,
                                display_numbers = TRUE,
                                show_rownames = TRUE,
                                show_colnames = TRUE,
                                ...) {
  stopifnot(is.list(results), is.character(file.prefix), is.character(path))

  for (mt in names(results)) {
    # Generate heatmap plot
    pobj <- pheatmap::pheatmap(ReadWriter::column.2.row.names(results[[mt]]),
      main = paste("Heatmap of", mt, "values"),
      scale = "column",
      cluster_rows = cluster_rows,
      display_numbers = display_numbers,
      show_rownames = show_rownames,
      show_colnames = show_colnames
    )

    # Construct file name
    file_name <- paste0(file.prefix, mt, ".png")
    file_path <- file.path(path, file_name)

    # Save plot
    MarkdownReports::wplot_save_pheatmap(x = pobj, plotname = file_name, png = TRUE, pdf = FALSE, ...)
    cat("Saved heatmap for", mt, "to", file_path, "\n")
  } # for
}


# _________________________________________________________________________________________________
# plotting generic, misc ______________________________ ----
# _________________________________________________________________________________________________


# _________________________________________________________________________________________________
#' @title Scatter Plot of Two Features in Seurat Object
#'
#' @description Generates a scatter plot comparing two features (genes or metrics) from a Seurat
#' object and optionally saves it. The function wraps around Seurat's `FeatureScatter` for
#' enhanced usability, including optional logarithmic transformations and saving capabilities.
#'
#' @param feature1 The first feature for the scatter plot's x-axis.
#' @param feature2 The second feature for the scatter plot's y-axis.
#' @param obj Seurat object containing the data for features.
#' @param ext File extension for saving the plot, if enabled.
#' @param plot Flag to display the plot within the R session.
#' @param logX Apply logarithmic transformation to x-axis values.
#' @param logY Apply logarithmic transformation to y-axis values.
#' @param ... Additional parameters passed to Seurat's `FeatureScatter`.
#'
#' @return A `ggplot` object of the feature scatter plot if `plot` is TRUE.
#'
#' @examples
#' \dontrun{
#' # Generate and display a scatter plot for features TOP2A and ID2
#' qFeatureScatter(feature1 = "TOP2A", feature2 = "ID2", obj = yourSeuratObject)
#' }
#'
#' @seealso \code{\link[Seurat]{FeatureScatter}}, \code{\link[ggplot2]{ggplot}}
#'
#' @export
#' @importFrom ggExpress qqSave
#' @importFrom Seurat FeatureScatter
#' @importFrom ggplot2 ggtitle theme_linedraw scale_x_log10 scale_y_log10
qFeatureScatter <- function(
    feature1 = "TOP2A", feature2 = "ID2", obj = combined.obj,
    ext = "png", plot = TRUE,
    logX = FALSE, logY = FALSE,
    ...) {
  plotname <- kpp(feature1, "VS", feature2)
  p <- FeatureScatter(object = obj, feature1 = feature1, feature2 = feature2, ...) +
    ggtitle(paste("Correlation", plotname)) +
    theme_linedraw()

  if (logX) p <- p + scale_x_log10()
  if (logY) p <- p + scale_y_log10()

  fname <- kpp("FeatureScatter", plotname)
  ggExpress::qqSave(ggobj = p, title = plotname, ext = ext, w = 8, h = 5)
  if (plot) p
}


# _________________________________________________________________________________________________
#' @title Create a Violin Plot for a Seurat Object Feature and save the file.
#'
#' @description Generates a violin plot for a specified feature in a Seurat object,
#' allowing for the data to be split by a specified grouping variable.
#' The function supports customization options such as logarithmic scaling, custom titles, and more.
#'
#' @param obj A Seurat object to be plotted.
#' @param features A character string specifying the name of the feature to plot.
#' @param idents A character vector specifying the identities to be used in the plot.
#' @param split.by A character string specifying the grouping variable for splitting the plot.
#' @param replace.na A logical indicating whether NA values should be replaced.
#' @param suffix An optional string to append to the title of the plot.
#' @param suffix.2.title A logical indicating whether to append the suffix to the plot title.
#' @param logY A logical indicating whether to use a logarithmic scale for the y-axis.
#' @param hline A numeric or logical value; if numeric, the value where a horizontal line should be drawn.
#' @param caption A character string or logical for the plot caption. If FALSE, no caption is displayed.
#' @param show_plot A logical indicating whether to display the plot.
#' @param w Width of the plot.
#' @param h Height of the plot.
#' @param ... Additional arguments passed to `VlnPlot`.
#'
#' @return A ggplot object representing the violin plot.
#'
#' @examples
#' # Assuming `seurat_obj` is a valid Seurat object
#' qSeuViolin(obj = seurat_obj, features = "nFeature_RNA")
#'
#' @export
qSeuViolin <- function(
    obj,
    features = "nFeature_RNA",
    idents = GetNamedClusteringRuns(obj)[1],
    split.by = NULL,
    replace.na = FALSE,
    suffix = NULL,
    suffix.2.title = FALSE,
    caption = .parseKeyParams(obj),
    logY = TRUE,
    hline = FALSE,
    ylab = "Expression",
    show_plot = TRUE,
    ylimit = NULL,
    w = 9, h = 5,
    ...) {
  #
  stopifnot(
    "Seurat" %in% class(obj), # object must be a Seurat object
    is.logical(logY), # logY must be logical (TRUE or FALSE)
    is.logical(hline) || is.numeric(hline), # hline must be logical or numeric
    is.logical(caption) || is.character(caption), # caption must be logical or character
    is.logical(suffix.2.title), # suffix.2.title must be logical
    is.character(split.by) | is.null(split.by), # split.by must be a character or NULL
    split.by %in% colnames(obj@meta.data),
    is.character(idents),
    idents %in% colnames(obj@meta.data),
    is.character(features),
    features %in% colnames(obj@meta.data) || features %in% rownames(obj)
  )

  ttl <- if (suffix.2.title) {
    paste(features, "|", suffix)
  } else {
    as.character(features)
  }
  subt <- paste(features, "- by -", idents)

  if (replace.na) {
    warning("NA's are not, but zeros are displayed on the plot. Avoid replace.na when possible", immediate. = TRUE)
    obj@meta.data[[features]] <- na.replace(x = obj@meta.data[[features]], replace = 0)
  }

  p <- VlnPlot(object = obj, features = features, split.by = split.by, group.by = idents, ...) +
    theme(axis.title.x = element_blank()) +
    labs(y = ylab) +
    ggtitle(label = ttl, subtitle = subt )

  # Add additional customization, if needed..
  if (!is.null(ylimit)) p <- p + ylim(ylimit[1], ylimit[2])
  if (logY) p <- p + ggplot2::scale_y_log10()
  if (hline) p <- p + ggplot2::geom_hline(yintercept = hline)
  if (!isFALSE(caption)) p <- p + ggplot2::labs(caption = caption)

  # Save the plot.
  TTL <- ppp(as.character(features), suffix, flag.nameiftrue(logY))
  qqSave(p, title = TTL, w = w, h = h)
  if (show_plot) p
}




# _________________________________________________________________________________________________
# Plotting 2D UMAPs, etc. ______________________________ ----
# _________________________________________________________________________________________________

# _________________________________________________________________________________________________
#' @title Quick UMAP Visualization of Gene Expression and automatically save the plot
#'
#' @description Generates a UMAP visualization for a specific feature from a Seurat object, and
#' automatically saves it. Offers options for custom titles, subtitles, saving, and more. Assumes
#' default options for custom titles, subtitles, saving, and more.
#'
#' @param feature Feature to visualize on the UMAP; Default: 'TOP2A'.
#' @param obj Seurat object containing single-cell RNA-seq data; Default: `combined.obj`.
#' @param title Title of the plot; Default: `feature`.
#' @param sub Subtitle of the plot; Default: NULL.
#' @param reduction Dimension reduction technique to be used ('umap', 'tsne', or 'pca'); Default: 'umap'.
#' @param splitby Column in the metadata to split the cells by; Default: NULL.
#' @param prefix Prefix added before the filename; Default: NULL.
#' @param suffix Suffix added to the end of the filename; Default: `sub`.
#' @param save.plot If TRUE, the plot is saved into a file; Default: TRUE.
#' @param PNG If TRUE, the file is saved as a .png; Default: TRUE.
#' @param h Height of the plot in inches; Default: 7.
#' @param w Width of the plot in inches; Default: NULL.
#' @param nr.cols Number of columns to combine multiple feature plots, ignored if `split.by` is not NULL; Default: NULL.
#' @param assay Which assay to use ('RNA' or 'integrated'); Default: 'RNA'.
#' @param axes If TRUE, axes are shown on the plot; Default: FALSE.
#' @param aspect.ratio Ratio of height to width. If TRUE, the ratio is fixed at 0.6; Default: FALSE.
#' @param HGNC.lookup If TRUE, HGNC gene symbol lookup is performed; Default: TRUE.
#' @param make.uppercase If TRUE, feature names are converted to uppercase; Default: FALSE.
#' @param qlow Lower quantile for the color scale; Default: 'q10'.
#' @param qhigh Upper quantile for the color scale; Default: 'q90'.
#' @param check_for_2D If TRUE, checks if UMAP is 2 dimensional; Default: TRUE.
#' @param caption Adds a caption to the ggplot object; Default: dynamically generated from `obj`.
#' @param ... Additional parameters to pass to the internally called functions.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   qUMAP(feature = "nFeature_RNA", obj = yourSeuratObject)
#'   qUMAP(feature = "TOP2A", obj = yourSeuratObject, PNG = FALSE, save.plot = TRUE)
#' }
#' }
#'
#' @export
#' @importFrom Seurat FeaturePlot NoLegend NoAxes
#' @importFrom ggplot2 ggtitle coord_fixed labs
#'
qUMAP <- function(
    feature = "TOP2A", obj = combined.obj,
    title = feature, sub = NULL,
    reduction = "umap", splitby = NULL,
    prefix = NULL,
    suffix = make.names(sub),
    save.plot = MarkdownHelpers::TRUE.unless("b.save.wplots", v = FALSE),
    PNG = TRUE,
    h = 7, w = NULL, nr.cols = NULL,
    assay = c("RNA", "integrated")[1],
    axes = FALSE,
    aspect.ratio = c(FALSE, 0.6)[1],
    HGNC.lookup = TRUE,
    make.uppercase = FALSE,
    check_for_2D = TRUE,
    qlow = "q10", qhigh = "q90",
    caption = .parseBasicObjStats(obj, simple = TRUE),
    ...) {
  # Checks
  if (check_for_2D) {
    umap_dims <- ncol(obj@reductions[[reduction]]@cell.embeddings)
    if (umap_dims != 2) warning(">>> UMAP is not 2 dimensional! \n Check obj@reductions[[reduction]]@cell.embeddings")
  }

  if (feature %in% colnames(obj@meta.data)) {
    message("feature found in meta.data")
    stopifnot(is.numeric(obj@meta.data[, feature]))
  }

  if (!(feature %in% colnames(obj@meta.data) | feature %in% rownames(obj))) {
    feature <- check.genes(
      list.of.genes = feature, obj = obj, verbose = FALSE,
      HGNC.lookup = HGNC.lookup, makeuppercase = make.uppercase
    )
  }

  DefaultAssay(obj) <- assay
  gg.obj <- Seurat::FeaturePlot(obj,
    features = feature,
    reduction = reduction,
    min.cutoff = qlow, max.cutoff = qhigh,
    ncol = nr.cols,
    split.by = splitby,
    ...
  ) +
    ggtitle(label = title, subtitle = sub) +
    if (!axes) NoAxes() else NULL

  if (aspect.ratio) gg.obj <- gg.obj + ggplot2::coord_fixed(ratio = aspect.ratio)
  if (!isFALSE(caption)) gg.obj <- gg.obj + ggplot2::labs(caption = caption)

  if (save.plot) {
    fname <- ww.FnP_parser(Stringendo::sppp(prefix, toupper(reduction), feature, assay, paste0(ncol(obj),"c"), suffix), if (PNG) "png" else "pdf")
    try(save_plot(filename = fname, plot = gg.obj, base_height = h, base_width = w)) # , ncol = 1, nrow = 1
  }
  return(gg.obj)
}



# _________________________________________________________________________________________________
#' @title Quick Visualization of Clustering Results with UMAP and automatically save the plot
#'
#' @description Generates a UMAP visualization based on clustering results from a Seurat object,
#' and automatically saves it. Offers options for custom titles, subtitles, saving, and more. Assumes
#' default options for custom titles, subtitles, saving, and more.
#'
#' @param ident Cluster identity for visualization; Default: 'integrated_snn_res.0.5'.
#' @param obj Seurat object containing single-cell data; Default: `combined.obj`.
#' @param reduction Dimension reduction method ('umap', 'tsne', 'pca'); Default: 'umap'.
#' @param splitby Metadata column to split cells by; optional; Default: NULL.
#' @param title Main title of the plot; Default: `ident`.
#' @param sub Subtitle of the plot; optional; Default: NULL.
#' @param prefix Prefix for saved filename; optional; Default: NULL.
#' @param suffix Suffix for saved filename; defaults to plot subtitle; Default: NULL.
#' @param label.cex Size of cluster labels; Default: 7.
#' @param h Height of plot in inches; Default: 7.
#' @param w Width of plot in inches; optional; Default: NULL.
#' @param nr.cols Number of columns for facet wrap if `splitby` is not NULL; Default: NULL.
#' @param plotname Custom plot name for saving; Default: dynamically generated from `reduction` and `ident`.
#' @param cols Custom color vector for clusters; optional; Default: NULL.
#' @param palette Color palette for generating cluster colors; Default: 'glasbey'.
#' @param highlight.clusters Specific clusters to be highlighted; optional; Default: NULL.
#' @param cells.highlight Specific cells to be highlighted; optional; Default: NULL.
#' @param label Show cluster labels; Default: TRUE.
#' @param repel Repel labels to avoid overlap; Default: TRUE.
#' @param legend Show legend; Default: opposite of `label`.
#' @param legend.pos Position of legend; Default: 'NULL'.
#' @param axes Show axes; Default: FALSE.
#' @param aspect.ratio Fixed aspect ratio for the plot; Default: TRUE.
#' @param MaxCategThrHP Maximum number of categories before simplification; Default: 200.
#' @param save.plot Save plot to file; Default: TRUE.
#' @param PNG Save as PNG (TRUE) or PDF (FALSE); Default: TRUE.
#' @param check_for_2D Ensure UMAP is 2D; Default: TRUE.
#' @param caption Plot caption; optional; Default: dynamically generated from `obj`.
#' @param ... Additional parameters for `DimPlot`.
#'
#' @examples
#' \dontrun{
#' clUMAP(ident = "integrated_snn_res.0.5", obj = yourSeuratObj)
#' clUMAP(ident = "integrated_snn_res.0.5", obj = yourSeuratObj, cols = RColorBrewer::brewer.pal(8, "Dark2"))
#' }
#'
#' @importFrom ggplot2 ggtitle labs coord_fixed ggsave
#' @importFrom Seurat DimPlot NoLegend NoAxes
#' @importFrom RColorBrewer brewer.pal
#'
#' @export
clUMAP <- function(
    ident = GetNamedClusteringRuns()[1],
    obj = combined.obj,
    reduction = "umap", splitby = NULL,
    title = ident,
    sub = NULL,
    prefix = NULL,
    suffix = make.names(sub),
    label.cex = 7,
    h = 7, w = NULL, nr.cols = NULL,
    plotname = ppp(toupper(reduction), ident),
    cols = NULL,
    palette = c("alphabet", "alphabet2", "glasbey", "polychrome", "stepped")[3],
    highlight.clusters = NULL, cells.highlight = NULL,
    label = TRUE, repel = TRUE,
    legend = !label,
    legend.pos = NULL, # c("top", "bottom", "left", "right", "none")[2],
    MaxCategThrHP = 200,
    axes = NULL,
    aspect.ratio = c(FALSE, 0.6)[2],
    save.plot = MarkdownHelpers::TRUE.unless("b.save.wplots", v = FALSE),
    PNG = TRUE,
    check_for_2D = TRUE,
    caption = .parseKeyParams(obj),
    # caption = NULL,
    ...) {
  #
  if (check_for_2D) {
    umap_dims <- ncol(obj@reductions[[reduction]]@cell.embeddings)
    if (umap_dims != 2) warning(">>> UMAP is not 2 dimensional! \n Check obj@reductions[[reduction]]@cell.embeddings")
  }

  IdentFound <- (ident %in% colnames(obj@meta.data))
  if (!IdentFound) {
    ident <- GetClusteringRuns(obj = obj, pat = "_res.*[0,1]\\.[0-9]$")[1]
    iprint("Identity not found. Plotting", ident)
  }
  identity <- obj[[ident]]
  NtCategs <- length(unique(identity[, 1]))
  if (NtCategs > 1000) warning("More than 1000 levels! qUMAP?", immediate. = TRUE)


  if (!missing(highlight.clusters)) {
    if (!(all(highlight.clusters %in% identity[, 1]))) {
      MSG <- paste(
        "Some clusters not found in the object! Missing:",
        kppc(setdiff(highlight.clusters, unique(identity[, 1]))), "\nFrom:\n",
        kppc(sort(unique(identity[, 1])))
      )
      warning(MSG, immediate. = TRUE)
    }

    idx.ok <- identity[, 1] %in% highlight.clusters
    stopifnot("minimum 10 cells are needed" = sum(idx.ok) > 10)

    highlight.these <- rownames(identity)[idx.ok]
    PCT <- percentage_formatter(length(highlight.these) / ncol(obj), suffix = "or")
    if (is.null(sub)) sub <- paste(PCT, length(highlight.these), "cells in ", ident)

    title <- kppc(highlight.clusters)
  } else {
    highlight.these <- NULL
  }
  if (!missing(cells.highlight)) {
    highlight.these <- cells.highlight
  } # overwrite, if directly defined

  if (is.null(cols)) {
    # browser()
    cols <- if (NtCategs > 7) {
      getDiscretePaletteObj(
        ident.used = ident, palette.used = palette,
        obj = obj, show.colors = FALSE
      )
    }
  }
  if (!is.null(highlight.these)) {
    cols <- "lightgrey"
  }

  if (NtCategs > MaxCategThrHP) {
    iprint("Too many categories (", NtCategs, ") in ", ident, "- use qUMAP for continous variables.")
  } else {
    if (length(unique(identity)) < MaxCategThrHP) {
      gg.obj <-
        Seurat::DimPlot(
          object = obj, group.by = ident,
          cols = cols,
          reduction = reduction, split.by = splitby,
          ncol = nr.cols, cells.highlight = highlight.these,
          label = label, repel = repel, label.size = label.cex, ...
        ) +
        ggtitle(label = title, subtitle = sub) +
        if (!legend) NoLegend() else NULL
    }

    if (is.null(axes)) gg.obj <- gg.obj + NoAxes()
    if (!is.null(caption)) gg.obj <- gg.obj + labs(caption = caption)
    if (!is.null(legend.pos)) gg.obj <- gg.obj + theme(legend.position = legend.pos)
    if (aspect.ratio) gg.obj <- gg.obj + ggplot2::coord_fixed(ratio = aspect.ratio)
    if (legend) suffix <- paste0(suffix, ".lgnd")

    if (save.plot) {
      pname <- Stringendo::sppp(prefix, plotname, paste0(ncol(obj),"c"), suffix, sppp(highlight.clusters))
      fname <- ww.FnP_parser(pname, if (PNG) "png" else "pdf")
      try(save_plot(filename = fname, plot = gg.obj, base_height = h, base_width = w)) # , ncol = 1, nrow = 1
    }
    return(gg.obj)
  } # if not too many categories
}






# _________________________________________________________________________________________________
#' @title Highlight Selected Clusters on UMAP
#'
#' @description Generates a UMAP plot from a Seurat object with specified clusters highlighted.
#' It saves the resulting UMAP plot directly to the current working directory.
#'
#' @param obj Seurat object to be visualized; Default: `combined.obj`.
#' @param COI Vector of cluster IDs to highlight on the UMAP plot;
#' Default: `c("0", "2", "4")`.
#' @param ident Name of the metadata column containing cluster IDs;
#' Default: 'integrated_snn_res.0.3'.
#'
#' @return Saves a UMAP plot highlighting specified clusters to the current working directory.
#' The function itself does not return an object within R.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # GetClusteringRuns()[1] "integrated_snn_res.0.1"
#'   umapHiLightSel(obj = combined.obj, COI = c("0", "1"), ident = GetClusteringRuns()[1])
#' }
#' }
#'
#' @seealso \code{\link[Seurat]{DimPlot}}
#'
#' @export
#' @importFrom Seurat DimPlot
#' @importFrom ggplot2 ggsave
umapHiLightSel <- function(obj = combined.obj,
                           COI = c("0", "2", "4"),
                           ident = GetClusteringRuns()[1],
                           h = 7, w = 5,
                           show_plot = T,
                           ...) {
  stopifnot(is(obj, "Seurat"),
            "Ident no found the object!" =ident %in% colnames(obj@meta.data),
            "Not all clusters in COI are found the object!" = all(COI %in% unique(obj@meta.data[[ident]]))
            )

  cellsSel <- getCellIDs.from.meta(ident = ident, ident_values = COI, obj = obj)
  pl <- Seurat::DimPlot(obj,
    reduction = "umap",
    group.by = ident,
    label = TRUE,
    cells.highlight = cellsSel,
    ...
  )
  if(show_plot) print(pl)

  ggplot2::ggsave(filename = extPNG(kollapse("cells", COI, collapseby = ".")),
                  height = h, width = w)
}



# _________________________________________________________________________________________________
#' @title DimPlot.ClusterNames
#'
#' @description Plot UMAP with Cluster names.
#' @param obj Seurat object, Default: combined.obj
#' @param ident identity used, Default: 'cl.names.top.gene.res.0.5'
#' @param reduction UMAP, tSNE, or PCA (Dim. reduction to use), Default: 'umap'
#' @param title Title of the plot, Default: ident
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   DimPlot.ClusterNames()
#' }
#' }
#' @export
DimPlot.ClusterNames <- function(
    obj = combined.obj,
    ident = "cl.names.top.gene.res.0.5",
    reduction = "umap", title = ident, ...) {
  Seurat::DimPlot(
    object = obj, reduction = reduction, group.by = ident,
    label = TRUE, repel = TRUE, ...
  ) + NoLegend() + ggtitle(title)
}




# _________________________________________________________________________________________________
# Multiplex 2D UMAPs, etc. ______________________________ ----


# _________________________________________________________________________________________________
#' @title multiFeaturePlot.A4
#'
#' @description Save multiple FeaturePlots, as jpeg, on A4 for each gene, which are stored as a list of gene names.
#' @param list.of.genes List of gene names for which the plots are to be generated. No default.
#' @param obj Seurat object, Default: combined.obj
#' @param foldername Folder name to save the generated plots. Default: The name of the list of genes.
#' @param plot.reduction Dimension reduction technique to use for plots. Default: 'umap'
#' @param intersectionAssay The assay to intersect with, either 'RNA' or 'integrated'. Default: 'RNA'
#' @param layout Layout orientation of the plot. Default: 'wide'
#' @param colors Vector of colors to be used in the plot. Default: c("grey", "red")
#' @param nr.Col Number of columns in the plot grid. Default: 2
#' @param nr.Row Number of rows in the plot grid. Default: 4
#' @param cex Point size in the plot. Default: round(0.1/(nr.Col * nr.Row), digits = 2)
#' @param gene.min.exp Minimum gene expression level for plotting. Default: 'q01'
#' @param gene.max.exp Maximum gene expression level for plotting. Default: 'q99'
#' @param subdir Should plots be saved in a sub-directory? Default: TRUE
#' @param prefix Prefix for the plot filenames. Default: NULL
#' @param suffix Suffix for the plot filenames. Default: NULL
#' @param background_col Background color of the plots. Default: "white"
#' @param saveGeneList Should the list of genes be saved? Default: FALSE
#' @param w Width of the plot. Default: 8.27
#' @param h Height of the plot. Default: 11.69
#' @param scaling Scaling factor for plot size. Default: 1
#' @param aspect.ratio Should the aspect ratio be fixed? Default: Yes, at 0.6
#' @param format Format to save the plot file. Default: 'jpg'
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @seealso
#'  \code{\link[tictoc]{tic}}
#'  \code{\link[cowplot]{plot_grid}}
#' @importFrom tictoc tic toc
#' @importFrom cowplot plot_grid
#' @importFrom MarkdownReports create_set_OutDir
#'
#' @export
multiFeaturePlot.A4 <- function(
    list.of.genes,
    obj = combined.obj,
    foldername = substitute(list.of.genes), plot.reduction = "umap",
    intersectionAssay = c("RNA", "integrated")[1],
    layout = c("tall", "wide", FALSE)[2],
    colors = c("grey", "red"),
    nr.Col = 2, nr.Row = 4,
    raster = if (ncol(obj) > 1e5) TRUE else FALSE,
    cex = round(0.1 / (nr.Col * nr.Row), digits = 2),
    cex.min = if (raster) TRUE else FALSE,
    gene.min.exp = "q01", gene.max.exp = "q99", subdir = TRUE,
    prefix = NULL, suffix = NULL,
    background_col = "white",
    aspect.ratio = c(FALSE, 0.6)[2],
    saveGeneList = FALSE,
    w = 8.27, h = 11.69, scaling = 1,
    format = c("jpg", "pdf", "png")[1],
    ...) {
  tictoc::tic()
  ParentDir <- OutDir
  if (is.null(foldername)) foldername <- "genes"
  final.foldername <- FixPlotName(paste0(foldername, "-", plot.reduction, suffix))
  if (subdir) create_set_SubDir(final.foldername, "/")

  list.of.genes.found <- check.genes(
    list.of.genes = list.of.genes, obj = obj,
    assay.slot = intersectionAssay, makeuppercase = FALSE
  )
  DefaultAssay(obj) <- intersectionAssay

  if (!is.null(cex.min)) cex <- max(cex.min, cex)

  if (layout == "tall") {
    w <- 8.27 * scaling
    h <- 11.69 * scaling
    nr.Col <- 2
    nr.Row <- 4
    print("layout active, nr.Col ignored.")
  }
  if (layout == "wide") {
    w <- 11.69 * scaling
    h <- 8.27 * scaling
    nr.Col <- 2
    nr.Row <- 2
    print("layout active, nr.Col ignored.")
  }

  lsG <- CodeAndRoll2::split_vec_to_list_by_N(1:length(list.of.genes.found), by = nr.Row * nr.Col)
  for (i in 1:length(lsG)) {
    genes <- list.of.genes.found[lsG[[i]]]
    iprint(i, genes)
    plotname <- kpp(c(prefix, plot.reduction, i, genes, suffix, format))

    plot.list <- Seurat::FeaturePlot(
      object = obj, features = genes, reduction = plot.reduction, combine = FALSE,
      ncol = nr.Col, cols = colors, raster = raster,
      min.cutoff = gene.min.exp, max.cutoff = gene.max.exp,
      pt.size = cex, ...
    )

    for (i in 1:length(plot.list)) {
      plot.list[[i]] <- plot.list[[i]] + NoLegend() + NoAxes()
      if (aspect.ratio) plot.list[[i]] <- plot.list[[i]] + ggplot2::coord_fixed(ratio = aspect.ratio)
    }

    pltGrid <- cowplot::plot_grid(plotlist = plot.list, ncol = nr.Col, nrow = nr.Row)
    # cowplot::ggsave2(filename = plotname, width = w, height = h, bg = background_col, plot = pltGrid)
    cowplot::save_plot(
      plot = pltGrid, filename = plotname,
      base_width = w, base_height = h,
      bg = background_col
    )
  }

  if (subdir) MarkdownReports::create_set_OutDir(ParentDir)
  if (saveGeneList) {
    if (is.null(obj@misc$gene.lists)) obj@misc$gene.lists <- list()
    obj@misc$gene.lists[[substitute(list.of.genes)]] <- list.of.genes.found
    print("Genes saved under: obj@misc$gene.lists")
    return(obj)
  }
  tictoc::toc()
}



# ____________________________________________________________________________________
#' @title Generate Cluster Highlight UMAPs compiled into A4 pages
#'
#' @description This function generates and saves cluster highlight plots for both single and multiple
#' clusters using UMAP or other dimensionality reduction techniques. It supports saving plots in various
#' formats and allows customization of plot appearance and layout.
#'
#' @param ident The name of the metadata column in the Seurat object `obj` to use for identifying clusters.
#' @param obj A Seurat object combining multiple datasets. Default: `combined.obj`.
#' @param foldername Name of the folder to save the plots in. Default: Value of `ident`.
#' @param plot.reduction The dimensionality reduction technique to use for the plots. Default: `"umap"`.
#' @param intersectionAssay The assay to use when calculating intersections. Default: `"RNA"`.
#' @param layout Plot layout, can be `"tall"`, `"wide"`, or `FALSE` for no specific layout. Default: `"wide"`.
#' @param colors A vector of colors to use for non-highlighted and highlighted clusters. Default: `c("grey", "red")`.
#' @param nr.Col Number of columns in the plot grid. Default: 2.
#' @param nr.Row Number of rows in the plot grid. Default: 4.
#' @param cex Size of the text in the plot, calculated based on the number of rows and columns. Default: Calculated value.
#' @param subdir Logical flag indicating whether to create a subdirectory for the plots. Default: `TRUE`.
#' @param prefix Optional prefix for the plot file names. Default: `NULL`.
#' @param suffix Optional suffix for the plot file names. Default: `NULL`.
#' @param background_col Background color of the plots. Default: `"white"`.
#' @param aspect.ratio Aspect ratio of the plots, can be `FALSE` for default ratio or a numeric value. Default: 0.6.
#' @param saveGeneList Logical flag indicating whether to save the list of genes used in the plots. Default: `FALSE`.
#' @param w Width of the plots, in inches. Default: `8.27`.
#' @param h Height of the plots, in inches. Default: `11.69`.
#' @param scaling Scaling factor for adjusting the size of the plots. Default: 1.
#' @param format Format to save the plots in, can be `"jpg"`, `"pdf"`, or `"png"`. Default: `"jpg"`.
#' @param ... Additional arguments passed to lower-level plotting functions.
#'
#' @return Invisible. This function primarily saves plots to files.
#' @examples
#' multiSingleClusterHighlightPlots.A4(ident = "cluster_id", obj = yourSeuratObject)
#'
#' @importFrom ggplot2 ggplot geom_point
#' @importFrom cowplot plot_grid ggsave2
#' @importFrom tictoc tic toc
#' @importFrom MarkdownReports create_set_OutDir
#'
#' @export
multiSingleClusterHighlightPlots.A4 <- function(
    obj = combined.obj,
    ident = GetClusteringRuns(obj)[1],
    foldername = ident,
    plot.reduction = "umap",
    intersectionAssay = DefaultAssay(combined.obj), # c("RNA", "integrated")[1],
    layout = c("tall", "wide", FALSE)[2],
    colors = c("grey", "red"),
    nr.Col = 2, nr.Row = 4,
    cex = round(0.1 / (nr.Col * nr.Row), digits = 2),
    subdir = TRUE,
    prefix = NULL, suffix = NULL,
    background_col = "white",
    aspect.ratio = c(FALSE, 0.6)[2],
    saveGeneList = FALSE,
    w = 8.27, h = 11.69, scaling = 1,
    format = c("jpg", "pdf", "png")[1],
    ...) {
  tictoc::tic()
  ParentDir <- OutDir
  if (is.null(foldername)) foldername <- "clusters"
  if (subdir) create_set_SubDir(paste0(foldername, "-", plot.reduction), "/")

  clusters <- unique(obj@meta.data[[ident]])

  DefaultAssay(obj) <- intersectionAssay

  # Adjust plot dimensions and grid layout based on specified layout
  if (layout == "tall") {
    w <- 8.27 * scaling
    h <- 11.69 * scaling
    nr.Col <- 2
    nr.Row <- 4
    message("tall layout active, nr.Col ignored.")
  }
  if (layout == "wide") {
    w <- 11.69 * scaling
    h <- 8.27 * scaling
    nr.Col <- 2
    nr.Row <- 2
    message("wide layout active, nr.Col ignored.")
  }


  # Split clusters into lists for plotting
  ls.Clust <- CodeAndRoll2::split_vec_to_list_by_N(1:length(clusters), by = nr.Row * nr.Col)
  for (i in 1:length(ls.Clust)) {
    clusters_on_this_page <- clusters[ls.Clust[[i]]]
    iprint("page:", i, "| clusters", kppc(clusters_on_this_page))
    (plotname <- kpp(c(prefix, plot.reduction, i, "clusters", ls.Clust[[i]], suffix, format)))

    plot.list <- list()
    for (i in seq(clusters_on_this_page)) {
      cl <- clusters_on_this_page[i]
      message(cl)
      plot.list[[i]] <- clUMAP(
        ident = ident, obj = obj,
        highlight.clusters = cl, label = FALSE, legend = FALSE, save.plot = FALSE,
        plotname = plotname, cols = colors, h = h, w = w, ...
      )
    }

    # Customize plot appearance
    for (i in 1:length(plot.list)) {
      plot.list[[i]] <- plot.list[[i]] + NoLegend() + NoAxes()
      if (aspect.ratio) {
        plot.list[[i]] <- plot.list[[i]] +
          ggplot2::coord_fixed(ratio = aspect.ratio)
      }
    }

    # Save plots
    pltGrid <- cowplot::plot_grid(plotlist = plot.list, ncol = nr.Col, nrow = nr.Row)
    cowplot::ggsave2(filename = plotname, width = w, height = h, bg = background_col, plot = pltGrid)
  } # for ls.Clust

  if (subdir) MarkdownReports::create_set_OutDir(ParentDir)
  tictoc::toc()
}




# _________________________________________________________________________________________________
#' @title Quick Clustering UMAPs on A4 Page
#'
#' @description Generates and arranges UMAP plots for up to four specified clustering resolutions
#' from a Seurat object onto an A4 page, facilitating comparative visualization.
#'
#' @param obj Seurat object to visualize; Default: `combined.obj`.
#' @param idents Vector of clustering resolution identifiers to plot;
#' dynamically defaults to the first 4 found by `GetClusteringRuns`.
#' @param prefix Prefix for plot titles; Default: "Clustering.UMAP.Res".
#' @param suffix Suffix for plot titles; Default: "".
#' @param title Custom title for the composite plot; dynamically generated from `prefix`, `idents`, and `suffix`.
#' @param nrow Number of rows in the plot grid; Default: 2.
#' @param ncol Number of columns in the plot grid; Default: 2.
#' @param w Width of the plot; Default: 11.69.
#' @param h Height of the plot; Default: 8.27.
#' @param ... Additional parameters for individual UMAP plots.
#'
#' @examples
#' \dontrun{
#' qClusteringUMAPS()
#' }
#'
#' @export
#' @importFrom Seurat NoAxes
#' @importFrom ggExpress qA4_grid_plot
qClusteringUMAPS <- function(
    obj = combined.obj,
    idents = na.omit.strip(GetClusteringRuns(obj)[1:4]),
    prefix = "Clustering.UMAP.Res",
    suffix = "",
    nrow = 2, ncol = 2,
    w = 11.69, h = 8.27,
    title = sppu(
      prefix,
      as.numeric(stringr::str_extract(idents, "\\d+\\.\\d+$")),
      suffix
    ),
    ...) {
  message("Plotting qClusteringUMAPS")

  # Check that the QC markers are in the object
  n.found <- intersect(idents, colnames(obj@meta.data))
  stopifnot("None of the idents found" = length(n.found) > 1,
            "Only 4 res's allowed" = length(n.found) <5)
  message(kppws(length(n.found), " found of ", idents))

  px <- list(
    "A" = clUMAP(ident = idents[1], save.plot = FALSE, obj = obj, caption = NULL, ...) + NoAxes(),
    "B" = clUMAP(ident = idents[2], save.plot = FALSE, obj = obj, caption = NULL, ...) + NoAxes(),
    "C" = clUMAP(ident = idents[3], save.plot = FALSE, obj = obj, caption = NULL, ...) + NoAxes(),
    "D" = clUMAP(ident = idents[4], save.plot = FALSE, obj = obj, ...) + NoAxes()
  )

  ggExpress::qA4_grid_plot(
    plot_list = px,
    plotname = title,
    w = w, h = h,
    nrow = nrow, ncol = ncol
  )
}

# _________________________________________________________________________________________________
#' @title Quickly Draw 4 Gene Expression UMAPs on an A4 Page
#'
#' @description Generates and arranges UMAP plots for up to four specified gene expressions
#' from a Seurat object onto an A4 page, facilitating comparative visualization.
#'
#' @param obj Seurat object to visualize; Default: `combined.obj`.
#' @param features Vector of gene identifiers to plot;
#' dynamically defaults to the first 4 found by `rownames(obj)`.
#' @param prefix Prefix for plot titles; Default: "Expression.UMAP.Gene".
#' @param suffix Suffix for plot titles; Default: "".
#' @param title Custom title for the composite plot; dynamically generated from `prefix`, `genes`, and `suffix`.
#' @param nrow Number of rows in the plot grid; Default: 2.
#' @param ncol Number of columns in the plot grid; Default: 2.
#' @param w Width of the plot; Default: 11.69.
#' @param h Height of the plot; Default: 8.27.
#' @param ... Additional parameters for individual UMAP plots.
#'
#' @examples
#' \dontrun{
#' qGeneExpressionUMAPS()
#' }
#'
#' @export
#' @importFrom Seurat NoAxes
#' @importFrom ggExpress qA4_grid_plot
qGeneExpressionUMAPS <- function(
    obj = combined.obj,
    features = rownames(obj)[1:4],
    prefix = "Expression.UMAP.Gene",
    suffix = "",
    nrow = 2, ncol = 2,
    w = 11.69, h = 8.27,
    title = paste0(prefix, " ", paste(features, collapse = ", "), " ", suffix),
    ...) {
  message("Plotting qGeneExpressionUMAPS")

  # Check that the features are in the object
  n.found <- intersect(features, c(colnames(obj@meta.data), rownames(obj)))
  stopifnot("None of the features found" = length(n.found) > 1,
            "Only 4 features are allowed" = length(n.found) <5)

  message(kppws(length(n.found), "found of", length(features), "features:", features))

  px <- list(
    "A" = qUMAP(feature = features[1], save.plot = FALSE, obj = obj, caption = NULL, ...) + NoAxes(),
    "B" = qUMAP(feature = features[2], save.plot = FALSE, obj = obj, caption = NULL, ...) + NoAxes(),
    "C" = qUMAP(feature = features[3], save.plot = FALSE, obj = obj, caption = NULL, ...) + NoAxes(),
    "D" = qUMAP(feature = features[4], save.plot = FALSE, obj = obj, ...) + NoAxes()
  )

  ggExpress::qA4_grid_plot(
    plot_list = px,
    plotname = title,
    w = w, h = h,
    nrow = nrow, ncol = ncol
  )
}


# _________________________________________________________________________________________________
#' @title Plot qUMAPs for Genes in a Folder
#'
#' @description This function plots qUMAPs for a specified set of genes, storing the results in a
#' specified folder. If no folder name is provided, it defaults to using the gene set name.
#'
#' @param genes A vector of gene names to be plotted.
#' @param obj An object containing the UMAP and gene data. Default: combined.obj.
#' @param foldername The name of the folder where the plots will be saved. If NULL, the gene set
#' name is used. Default: NULL.
#' @param intersectionAssay The assay slot to use for intersection. Default: 'RNA'.
#' @param plot.reduction The type of reduction to plot. Default: 'umap'.
#' @param ... Additional arguments passed to plotting and directory creation functions.
#'
#' @return Invisible. The function generates plots and saves them in the specified folder.
#'
#' @examples
#' plotQUMAPsInAFolder(
#'   genes = c("Gene1", "Gene2"), obj = combined.obj,
#'   foldername = "MyGenePlots", intersectionAssay = "RNA",
#'   plot.reduction = "umap"
#' )
#'
#' @importFrom MarkdownReports create_set_SubDir create_set_OutDir
#' @export

plotQUMAPsInAFolder <- function(genes, obj = combined.obj,
                                foldername = NULL,
                                intersectionAssay = DefaultAssay(obj),
                                plot.reduction = "umap",
                                ...) {
  # Input checks
  stopifnot(is.character(genes),
            is.null(foldername) || is.character(foldername),
            is.character(plot.reduction))

  ParentDir <- OutDir
  if (is.null(foldername)) foldername <- deparse(substitute(genes))

  MarkdownReports::create_set_SubDir(paste0(foldername, "-", plot.reduction), "/")

  list.of.genes.found <- check.genes(
    list.of.genes = genes, obj = obj,
    assay.slot = intersectionAssay, makeuppercase = FALSE
  )

  for (g in list.of.genes.found) {
    message(g)
    qUMAP(g, reduction = plot.reduction, obj = obj, ...)
  }

  MarkdownReports::create_set_OutDir(ParentDir)

  invisible()
}


# _________________________________________________________________________________________________
#' @title Plot Top N Differentially Expressed Genes Per Cluster
#'
#' @description Visualizes the top N differentially expressed (DE) genes for each cluster within a
#' specified clustering resolution of a Seurat object, facilitating the exploration of gene
#' expression patterns across clusters.
#'
#' @param obj Seurat object containing single-cell RNA-seq data and clustering information;
#' Default: `combined.obj`.
#' @param cl_res Cluster resolution used to identify distinct clusters for analysis; Default: `res`.
#' @param nrGenes Number of top DE genes to display for each cluster;
#' Default: GetClusteringRuns()[1].
#' @param order.by Criteria for ranking DE genes within clusters; Default: `"combined.score"`.
#' @param df_markers Data frame or list of DE genes across clusters. If not provided,
#' attempts to retrieve from `obj@misc$df.markers[[paste0("res.", cl_res)]]`;
#' Default: calculated based on `cl_res`.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   PlotTopGenesPerCluster(obj = combined.obj, cl_res = 0.5, nrGenes = 10)
#' }
#' }
#'
#' @export
PlotTopGenesPerCluster <- function(
    obj = combined.obj,
    cl_res = GetClusteringRuns()[1],
    nrGenes = p$"n.markers",
    order.by = c("combined.score", "avg_log2FC", "p_val_adj")[1],
    df_markers = obj@misc$"df.markers"[[paste0("res.", cl_res)]],
    ...) {
  message("Running PlotTopGenesPerCluster...")

  topX.markers <- GetTopMarkers(
    df = df_markers, n = nrGenes,
    order.by = order.by
  )

  ls.topMarkers <- splitbyitsnames(topX.markers)
  for (i in 1:length(ls.topMarkers)) {
    multiFeaturePlot.A4(
      list.of.genes = ls.topMarkers[[i]], obj = obj, subdir = FALSE,
      prefix = ppp("DEG.markers.res", cl_res, "cluster", names(ls.topMarkers)[i])
    )
  }
}

# _________________________________________________________________________________________________
#' @title Quickly Plot Key QC Markers in Brain Organoids
#'
#' @description Generates and arranges UMAP plots for specified QC features
#' from a Seurat object on an A4 page, facilitating a quick quality control (QC) overview.
#'
#' @param obj Seurat object to visualize; Default: `combined.obj`.
#' @param QC.Features Vector of QC feature names to plot; Default:
#' `c("nFeature_RNA", "percent.ribo", "percent.mito", "nuclear.fraction")`.
#' @param prefix Prefix for plot titles; Default: "QC.markers.4.UMAP".
#' @param suffix Suffix for plot titles; Default: "".
#' @param title Custom title for the composite plot; dynamically generated from `prefix`,
#' `QC.Features`, and `suffix`.
#' @param nrow Number of rows in the plot grid; Default: 2.
#' @param ncol Number of columns in the plot grid; Default: 2.
#' @param ... Additional parameters for individual UMAP plots.
#'
#' @examples
#' \dontrun{
#' qQC.plots.BrainOrg()
#' }
#'
#' @export
#' @importFrom ggExpress qA4_grid_plot
qQC.plots.BrainOrg <- function(
    obj = combined.obj,
    QC.Features = c("nFeature_RNA", "percent.ribo", "percent.mito", "nuclear.fraction"),
    prefix = "QC.markers.4.UMAP",
    suffix = "",
    title = sppu(prefix, QC.Features, suffix),
    nrow = 2, ncol = 2,
    ...) {
  message("> > > > > Plotting qQC.plots.BrainOrg")

  # Check that the QC markers are in the object
  n.found <- intersect(QC.Features, colnames(obj@meta.data))
  message(kppws(length(n.found), " found of ", QC.Features))
  stopifnot(length(n.found) > 1)

  # Count the number of NAs in specified columns
  na_counts <- sapply(X = obj@meta.data[, QC.Features], function(x) sum(is.na(x)))

  # Raise a warning if there are any NAs
  if (sum(na_counts) > 0) {
    warning(sprintf("There are %d NA values found\n", na_counts),
      immediate. = TRUE
    )
  }

  px <- list(
    "A" = qUMAP(QC.Features[1], save.plot = FALSE, obj = obj, ...) + NoAxes(),
    "B" = qUMAP(QC.Features[2], save.plot = FALSE, obj = obj, ...) + NoAxes(),
    "C" = qUMAP(QC.Features[3], save.plot = FALSE, obj = obj, ...) + NoAxes(),
    "D" = qUMAP(QC.Features[4], save.plot = FALSE, obj = obj, ...) + NoAxes()
  )

  ggExpress::qA4_grid_plot(
    plot_list = px,
    plotname = title,
    w = 11.69, h = 8.27,
    nrow = nrow, ncol = ncol
  )
}


# _________________________________________________________________________________________________
#' @title Quickly Plot Key Markers in Brain Organoids
#'
#' @description Generates plots for a predefined or custom set of gene markers within brain organoids,
#' aiding in the quick assessment of their expression across different cells or clusters.
#'
#' @param obj Seurat object for visualization; Default: `combined.obj`.
#' @param custom.genes Logical indicating whether to use a custom set of genes.
#' If FALSE, a predefined list of key brain organoid markers is used; Default: FALSE.
#' @param suffix Suffix for the folder name where the plots are saved; Default: "".
#'
#' @examples
#' \dontrun{
#' qMarkerCheck.BrainOrg(combined.obj)
#' qMarkerCheck.BrainOrg(combined.obj, custom.genes = c("Gene1", "Gene2"))
#' }
#'
#' @export
#' @importFrom CodeAndRoll2 as_tibble_from_namedVec

qMarkerCheck.BrainOrg <- function(obj = combined.obj, custom.genes = FALSE,
                                  suffix = "") {
  message("> > > > > Plotting qMarkerCheck.BrainOrg")

  Signature.Genes.Top16 <- if (!isFALSE(custom.genes)) {
    custom.genes
  } else {
    Signature.Genes.Top16 <- c(
      `dl-EN` = "KAZN", `ul-EN` = "SATB2" # dl-EN = deep layer excitatory neuron
      , `Immature neurons` = "SLA", Interneurons = "DLX6-AS1",
      Interneurons = "ERBB4", Interneurons = "SCGN",
      `Intermediate progenitor` = "EOMES" # ,  `Intermediate progenitor1` = "TAC3"
      , `S-phase` = "TOP2A", `G2M-phase` = "H4C3" # formerly: HIST1H4C
      , `oRG` = "HOPX" # , `oRG` = "ID4" # oRG outer radial glia
      # , Astroglia = "GFAP"
      , Astrocyte = "S100B", `Hypoxia/Stress` = "DDIT4",
      `Choroid.Plexus` = "TTR", `Low-Quality` = "POLR2A",
      `Mesenchyme` = "DCN", Glycolytic = "PDK1"
      # , `Choroid.Plexus` = "OTX2", `Mesenchyme` = "DCN"
    )
    print(Signature.Genes.Top16)
  }

  stopifnot()

  print(CodeAndRoll2::as_tibble_from_namedVec(Signature.Genes.Top16))
  multiFeaturePlot.A4(
    obj = obj, list.of.genes = Signature.Genes.Top16, layout = "tall",
    foldername = sppp("Signature.Genes.Top16", suffix)
  )
}





# _________________________________________________________________________________________________
#' @title Plot Top Genes
#'
#' @description This function plots the highest expressed genes on UMAPs, saving the plots in a
#' subfolder. It requires the prior execution of `calc.q99.Expression.and.set.all.genes`.
#'
#' @param obj A Seurat object containing the data for plotting. Default: combined.obj.
#' @param n The number of top genes to plot. Default: 32.
#' @param exp.slot The slot in the Seurat object where the expression data is stored.
#' Default: "expr.q99".
#'
#' @details This function identifies the top `n` genes based on expression levels stored in
#' `exp.slot` of the provided Seurat object. It then plots these genes using UMAPs and saves
#' the results in a subfolder named "Highest.Expressed.Genes".
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   PlotTopGenes()
#' }
#' }
#'
#' @export

PlotTopGenes <- function(obj = combined.obj, n = 32, exp.slot = "expr.q99") {
  message("Using obj@misc$", exp.slot)
  stopifnot(inherits(obj, "Seurat"),
    "Requires calling calc.q99.Expression.and.set.all.genes before. " =
      exp.slot %in% names(obj@misc)
  )

  Highest.Expressed.Genes <- names(head(sort(obj@misc[[exp.slot]], decreasing = TRUE), n = n))
  multiFeaturePlot.A4(list.of.genes = Highest.Expressed.Genes, foldername = "Highest.Expressed.Genes")
}






# _________________________________________________________________________________________________
# Manipulating UMAP and PCA  ______________________________ ----
# _________________________________________________________________________________________________

#' @title Flip Reduction Coordinates
#'
#' @description Flips dimensionality reduction coordinates (such as UMAP or tSNE) vertically or
#' horizontally to change the visualization perspective.
#'
#' @param obj Seurat object to modify; Default: `combined.obj`.
#' @param dim Number of dimensions in the reduction to consider; Default: 2.
#' @param reduction Dimension reduction technique to modify ('umap', 'tsne', or 'pca'); Default: 'umap'.
#' @param flip Axis (or axes) to flip; can be 'x', 'y', or 'xy' to flip both; Default: "x".
#' @param FlipReductionBackupToo Boolean indicating whether to also flip coordinates in the backup slot; Default: TRUE.
#'
#' @examples
#' \dontrun{
#' # Before flipping UMAP coordinates
#' clUMAP()
#' # Flip UMAP coordinates and visualize again
#' combined.obj <- FlipReductionCoordinates(combined.obj)
#' clUMAP()
#' }
#'
#' @export
#' @importFrom Seurat Embeddings
FlipReductionCoordinates <- function(
    obj = combined.obj, dim = 2, reduction = "umap",
    flip = c("x", "y", "xy", NULL)[1], FlipReductionBackupToo = TRUE) {
  coordinates <- Embeddings(obj, reduction = reduction)
  stopifnot(ncol(coordinates) == dim)

  if (flip %in% c("x", "xy")) coordinates[, 1] <- coordinates[, 1] * -1
  if (flip %in% c("y", "xy")) coordinates[, 2] <- coordinates[, 2] * -1
  obj@reductions[[reduction]]@cell.embeddings <- coordinates

  if (FlipReductionBackupToo) {
    bac.slot <- paste0(reduction, dim, "d")
    if (length(obj@misc$reductions.backup[[bac.slot]])) {
      obj@misc$reductions.backup[[bac.slot]]@cell.embeddings <- coordinates
      iprint(dim, "dimensional", reduction, "backup flipped too.")
    }
  }
  return(obj)
}



# _________________________________________________________________________________________________
#' @title Relabel Cluster Numbers Along a UMAP (or tSNE) Axis
#'
#' @description Automatically renumbers clusters based on their position along a specified dimension
#' in a UMAP (or tSNE or PCA) plot, potentially enhancing interpretability by ordering clusters.
#'
#' @param obj Seurat object containing clustering and UMAP (or other dimensional reduction) data;
#' Default: `combined.obj`.
#' @param dim Dimension along which to order clusters (1 for the first dimension, typically horizontal);
#' Default: 1.
#' @param swap If TRUE, reverses the ordering direction; Default: FALSE.
#' @param reduction Dimension reduction technique used for cluster positioning ('umap', 'tsne', or 'pca');
#' Default: 'umap'.
#' @param ident Clustering resolution identifier used to fetch cluster labels from `obj` metadata;
#' Default: 'integrated_snn_res.0.5'.
#' @param plot If TRUE, plots the UMAP with new cluster names; Default: TRUE.
#'
#' @examples
#' \dontrun{
#' combined.obj <- AutoNumber.by.UMAP(
#'   obj = combined.obj, dim = 1, reduction = "umap",
#'   ident = "integrated_snn_res.0.5"
#' )
#' DimPlot.ClusterNames(combined.obj, ident = GetClusteringRuns(combined.obj)[1])
#' }
#'
#' @export
#' @importFrom CodeAndRoll2 as.named.vector.df unlapply translate
#' @importFrom Stringendo kpp kppu iprint
#' @importFrom Seurat FetchData
AutoNumber.by.UMAP <- function(obj = combined.obj,
                               reduction = "umap",
                               dim = 1, swap = FALSE,
                               ident = GetClusteringRuns(obj = obj)[1],
                               plot = TRUE) {
  dim_name <- kppu(reduction, dim)
  if (obj@version < 5) dim_name <- toupper(dim_name)
  message("Obj. version: ", obj@version, " \ndimension name: ", dim_name)
  message("Resolution: ", ident)

  stopifnot("Identity not found." = ident %in% colnames(obj@meta.data))

  coord.umap <- obj@reductions$umap@cell.embeddings[, dim_name]

  # coord.umap <- round(coord.umap,digits = 2)
  identX <- as.character(obj@meta.data[[ident]])

  ls.perCl <- split(coord.umap, f = identX)
  MedianClusterCoordinate <- sapply(ls.perCl, median)
  # sort(MedianClusterCoordinate)

  OldLabel <- names(sort(MedianClusterCoordinate, decreasing = swap))
  NewLabel <- as.character(0:(length(MedianClusterCoordinate) - 1))
  NewMeta <- translate(vec = identX, oldvalues = OldLabel, newvalues = NewLabel)
  NewMetaCol <- kpp(ident, "ordered")
  iprint("NewMetaCol:", NewMetaCol)

  obj[[NewMetaCol]] <- NewMeta
  if (plot) {
    clUMAP(obj, ident = NewMetaCol)
  }
  return(obj)
}



# _________________________________________________________________________________________________
# Helpers ______________________________ ----
# _________________________________________________________________________________________________

#' @title Adjust Layout Parameters for multi* plotting fucntions
#'
#' @description Adjusts layout dimensions and properties based on the specified layout type.
#'              Updates the provided environment with new dimensions and layout configuration.
#'
#' @param layout A string specifying the layout type. Can be either "tall" or "wide". Default: NULL.
#' @param scaling A numeric scaling factor to adjust the dimensions. Default: 1.
#' @param wA4 The width of the A4 paper in inches. Default: 8.27.
#' @param hA4 The height of the A4 paper in inches. Default: 11.69.
#' @param env The environment where the layout dimensions and properties should be assigned.
#'            Default: parent.frame().
#'
#' @return Invisible NULL. The function operates by side effects, updating the `env` environment.
#' @examples
#' env <- new.env()
#' .adjustLayout("tall", 1, 8.27, 11.69, env)
#' print(env$w) # Should print the width based on "tall" layout scaling.
#'
.adjustLayout <- function(layout, scaling, wA4, hA4, env) {
  # Input checks
  stopifnot(
    is.character(layout), is.numeric(scaling), is.numeric(wA4),
    is.numeric(hA4), is.environment(env),
    layout %in% c("tall", "wide")
  )

  if (layout == "tall") {
    assign("w", wA4 * scaling, envir = env)
    assign("h", hA4 * scaling, envir = env)
    assign("nr.Col", 2, envir = env)
    assign("nr.Row", 4, envir = env)
    message("tall layout active, nr.Col ignored.")
  } else if (layout == "wide") {
    assign("w", hA4 * scaling, envir = env)
    assign("h", wA4 * scaling, envir = env)
    assign("nr.Col", 2, envir = env) # Adjusted for consistency with wide layout explanation
    assign("nr.Row", 2, envir = env)
    message("wide layout active, nr.Col ignored.")
  } else {
    message("No specific layout selected, defaulting to input parameters.")
  }
}


# _________________________________________________________________________________________________
# Saving plots ______________________________ ----
# _________________________________________________________________________________________________



# _________________________________________________________________________________________________
#' @title Save Two Plots on One A4 Page
#'
#' @description Arranges and saves two UMAP plots (or any plots) side-by-side or one above
#' the other on a single A4 page.
#'
#' @param plot_list A list containing ggplot objects to be arranged and saved.
#' @param pname Boolean indicating if the plot name should be automatically generated;
#' if FALSE, the name is based on `plot_list` and `suffix`; Default: FALSE.
#' @param suffix Suffix to be added to the generated filename if `pname` is FALSE; Default: NULL.
#' @param scale Scaling factor for adjusting the plot size; Default: 1.
#' @param nrow Number of rows in the plot arrangement; Default: 2.
#' @param ncol Number of columns in the plot arrangement; Default: 1.
#' @param h Height of the plot, calculated as A4 height times `scale`;
#' calculated dynamically based on `scale`.
#' @param w Width of the plot, calculated as A4 width times `scale`;
#' calculated dynamically based on `scale`.
#' @param ... Additional parameters passed to `plot_grid`.
#'
#' @examples
#' \dontrun{
#' p1 <- ggplot(iris, aes(Sepal.Length, Sepal.Width, color = Species)) +
#'   geom_point()
#' p2 <- ggplot(iris, aes(Petal.Length, Petal.Width, color = Species)) +
#'   geom_point()
#' save2plots.A4(plot_list = list(p1, p2))
#' }
#'
#' @export
#' @importFrom cowplot plot_grid save_plot ggdraw
#' @importFrom ggplot2 theme
save2plots.A4 <- function(
    plot_list, pname = FALSE, suffix = NULL, scale = 1,
    nrow = 2, ncol = 1,
    h = 11.69 * scale, w = 8.27 * scale, ...) {
  if (pname == FALSE) pname <- Stringendo::sppp(substitute(plot_list), suffix)
  p1 <- cowplot::plot_grid(
    plotlist = plot_list, nrow = nrow, ncol = ncol,
    labels = LETTERS[1:length(plot_list)], ...
  )
  p1 <- cowplot::ggdraw(p1) +
    theme(plot.background = element_rect(fill = "white", color = NA))

  iprint("Saved as:", pname)

  save_plot(plot = p1, filename = extPNG(pname), base_height = h, base_width = w)
}

# _________________________________________________________________________________________________
#' @title Save Four Plots on One A4 Page
#'
#' @description Arranges and saves four plots (e.g. UMAPs) onto a single A4 page, allowing for a
#' compact comparison of different visualizations or clustering results.
#'
#' @param plot_list A list containing ggplot objects to be arranged and saved; each object represents one panel.
#' @param pname Plot name; if FALSE, a name is generated automatically based on `plot_list` and `suffix`; Default: FALSE.
#' @param suffix Suffix to be added to the filename; Default: NULL.
#' @param scale Scaling factor for adjusting the size of the overall plot canvas; Default: 1.
#' @param nrow Number of rows to arrange the plots in; Default: 2.
#' @param ncol Number of columns to arrange the plots in; Default: 2.
#' @param h Height of the plot canvas, calculated as the height of an A4 page times `scale`; Default: `8.27 * scale`.
#' @param w Width of the plot canvas, calculated as the width of an A4 page times `scale`; Default: `11.69 * scale`.
#' @param ... Additional parameters passed to `plot_grid`.
#'
#' @examples
#' \dontrun{
#' p1 <- ggplot(iris, aes(Sepal.Length, Sepal.Width, color = Species)) +
#'   geom_point()
#' p2 <- ggplot(mtcars, aes(mpg, disp, color = as.factor(cyl))) +
#'   geom_point()
#' p3 <- ggplot(mpg, aes(displ, hwy, color = class)) +
#'   geom_point()
#' p4 <- ggplot(diamonds, aes(carat, price, color = cut)) +
#'   geom_point()
#' save4plots.A4(plot_list = list(p1, p2, p3, p4))
#' }
#'
#' @export
#' @importFrom cowplot plot_grid save_plot ggdraw
#' @importFrom ggplot2 theme
save4plots.A4 <- function(
    plot_list, pname = FALSE, suffix = NULL, scale = 1,
    nrow = 2, ncol = 2,
    h = 8.27 * scale, w = 11.69 * scale,
    ...) {
  if (pname == FALSE) pname <- Stringendo::sppp(substitute(plot_list), suffix)
  p1 <- cowplot::plot_grid(
    plotlist = plot_list, nrow = nrow, ncol = ncol,
    labels = LETTERS[1:length(plot_list)], ...
  )
  # https://stackoverflow.com/questions/13691415/change-the-background-color-of-grid-arrange-output
  p1 <- cowplot::ggdraw(p1) +
    theme(plot.background = element_rect(fill = "white", color = NA))

  iprint("Saved as:", pname)
  # fname <- MarkdownHelpers::ww.FnP_parser(extPNG(pname) )
  save_plot(plot = p1, filename = extPNG(pname), base_height = h, base_width = w)
}



# _________________________________________________________________________________________________
#' @title qqSaveGridA4
#'
#' @description Saves a grid of 2 or 4 ggplot objects onto an A4 page.
#' @param plotlist A list of ggplot objects. Default: pl.
#' @param plots A numeric vector indicating the indices of the plots to save from the 'plotlist'. Default: 1:2.
#' @param NrPlots Number of plots to save. Default: length(plots).
#' @param height Height for the saved image. Default: 11.69.
#' @param width Width for the saved image. Default: 8.27.
#' @param fname File name for the saved image. Default: "Fractions.Organoid-to-organoid variation.png".
#' @param ... Additional arguments passed to the plot_grid function.
#' @return This function does not return a value. It saves a grid plot of ggplot objects to the specified file.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   qqSaveGridA4(plotlist = pl, plots = 1:2, fname = "Fractions.per.Cl.png")
#'   qqSaveGridA4(plotlist = pl, plots = 1:4, fname = "Fractions.per.Cl.4.png")
#' }
#' }
#' @seealso
#' \code{\link[cowplot]{plot_grid}}
#' @importFrom cowplot plot_grid
#'
#' @export
qqSaveGridA4 <- function(
    plotlist = pl,
    plots = 1:2, NrPlots = length(plots), height = 11.69, width = 8.27,
    fname = "Fractions.Organoid-to-organoid variation.png", ...) {
  stopifnot(NrPlots %in% c(2, 4))
  iprint(NrPlots, "plots found,", plots, "are saved.")
  pg.cf <- cowplot::plot_grid(plotlist = plotlist[plots], nrow = 2, ncol = NrPlots / 2, labels = LETTERS[1:NrPlots], ...)
  if (NrPlots == 4) list2env(list(height = width, width = height), envir = as.environment(environment()))
  save_plot(
    filename = fname,
    plot = pg.cf, base_height = height, base_width = width
  )
  MarkdownHelpers::ww.FnP_parser(fname)
}



# _________________________________________________________________________________________________
# plotting.dim.reduction.3D.R ______________________________ ----
# _________________________________________________________________________________________________
# source('~/GitHub/Packages/Seurat.utils/Functions/plotting.dim.reduction.3D.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Plotting.dim.reduction.3D.R"))
# Source: self + https://github.com/Dragonmasterx87/Interactive-3D-Plotting-in-Seurat-3.0.0

# Requirements __________________________________________
# try(library(plotly), silent = TRUE)
# try(library(MarkdownReports), silent = TRUE)
# try(library(htmlwidgets), silent = TRUE)

# May also require
# try (source('~/GitHub/Packages/CodeAndRoll/CodeAndRoll.R'),silent= TRUE) # generic utilities funtions
# require('MarkdownReports') # require("devtools")



# _________________________________________________________________________________________________
#' @title Check Quantile Cutoff and Clip Outliers
#'
#' @description Checks a specified quantile cutoff and clips outliers from an expression vector,
#' ensuring that a minimum number of cells expressing a gene remain.
#'
#' @param expr.vec A numeric vector representing gene expression data.
#' @param quantileCutoffX The quantile cutoff for clipping outliers.
#' @param min.cells.expressing The minimum number of cells that should remain expressing after clipping.
#' @return The expression vector with outliers clipped, ensuring the minimum number of cells expressing.
#' @examples
#' \dontrun{
#' expr.vec <- c(...)
#' quantileCutoff <- 0.99
#' min.cells.expressing <- 10
#' ww.check.quantile.cutoff.and.clip.outliers(expr.vec, quantileCutoff, min.cells.expressing)
#' }
#' @export
#' @importFrom CodeAndRoll2 clip.outliers.at.percentile
#'
ww.check.quantile.cutoff.and.clip.outliers <- function(expr.vec = plotting.data[, gene],
                                                       quantileCutoffX = quantileCutoff,
                                                       min.cells.expressing = 10) {
  expr.vec.clipped <-
    CodeAndRoll2::clip.outliers.at.percentile(expr.vec, probs = c(1 - quantileCutoffX, quantileCutoffX))
  if (sum(expr.vec.clipped > 0) > min.cells.expressing) {
    expr.vec <- expr.vec.clipped
  } else {
    iprint("WARNING: quantile.cutoff too stringent, would leave <", min.cells.expressing, "cells. It is NOT applied.")
  }
  return(expr.vec)
}


# _________________________________________________________________________________________________
#' @title plot3D.umap.gene
#'
#' @description Plot a 3D umap with gene expression. Uses plotly. Based on github.com/Dragonmasterx87.
#' @param gene The gene of interest. Default: 'TOP2A'
#' @param obj The Seurat object for which the 3D umap plot will be generated. Default: combined.obj
#' @param quantileCutoff Cutoff value for the quantile for clipping outliers in the gene expression data. Default: 0.99
#' @param def.assay The default assay to be used. Choose between "integrated" and "RNA". Default: "RNA"
#' @param suffix A suffix added to the filename. Default: NULL
#' @param annotate.by The cluster or grouping to be used for automatic annotation. Default: First returned result from GetNamedClusteringRuns(obj) function.
#' @param alpha Opacity of the points in the plot. Default: 0.5
#' @param dotsize The size of the dots in the plot. Default: 1.25
#' @param ... Pass any other parameter to the internally called `plotly::plot_ly`.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   plot3D.umap.gene(obj = combined.obj, gene = "DDIT4", quantileCutoff = .95)
#'   plot3D.umap.gene(obj = combined.obj, gene = "percent.mito", quantileCutoff = .95) # for continous meta variables
#'   plot3D.umap.gene(obj = combined.obj, gene = "nFeature_RNA", quantileCutoff = .95) # for continous meta variables
#' }
#' }
#' @importFrom plotly plot_ly layout
#' @importFrom Seurat FetchData
#'
#' @export

plot3D.umap.gene <- function(
    gene = "TOP2A",
    obj = combined.obj,
    annotate.by = GetNamedClusteringRuns(obj)[1],
    quantileCutoff = .99,
    def.assay = c("integrated", "RNA")[2],
    suffix = NULL,
    alpha = .5,
    dotsize = 1.25,
    col.names = c("umap_1", "umap_2", "umap_3"),
    ...) {
  # Input assertions ____________________________________

  stopifnot(
    is(obj, "Seurat"),
    is.character(gene),
    "gene or feature not found in obj" = (gene %in% Features(obj) | gene %in% colnames(obj@meta.data)),
    "annotate.by not found in @meta" = (annotate.by %in% colnames(obj@meta.data) | annotate.by == FALSE),
    "reductions.backup is missing from @misc" = is.list(obj@misc$"reductions.backup"),
    "umap3d is missing from @misc$reductions.backup" = is(obj@misc$reductions.backup$"umap3d", class2 = "DimReduc"),
    "reductionn has 3 columns" = (ncol(obj@misc$reductions.backup$"umap3d") == 3),
    "3D reduction has >/< cells than object" = (ncol(obj) == nrow(obj@misc$reductions.backup$"umap3d"@cell.embeddings))
  )

  if (obj@version < 5) col.names <- toupper(col.names)
  message("Obj. version: ", obj@version, " \ndim names: ", kppc(col.names))

  DefaultAssay(object = obj) <- def.assay
  iprint(DefaultAssay(object = obj), "assay")

  # Get and format 3D plotting data ____________________________________
  plotting.data <- obj@misc$reductions.backup$"umap3d"@cell.embeddings
  colnames(plotting.data) <- toupper(col.names)

  # browser()
  Expression <- Seurat::FetchData(object = obj, vars = gene)
  plotting.data <- cbind(plotting.data, Expression)

  plotting.data$"Expression" <- ww.check.quantile.cutoff.and.clip.outliers(
    expr.vec = plotting.data[, gene],
    quantileCutoffX = quantileCutoff, min.cells.expressing = 10
  )
  # browser()
  # CodeAndRoll2::clip.outliers.at.percentile(plotting.data[, gene], probs = c(1 - quantileCutoff, quantileCutoff))
  plotting.data$"label" <- paste(rownames(plotting.data), " - ", plotting.data[, gene], sep = "")

  ls.ann.auto <- if (annotate.by != FALSE) {
    .Annotate4Plotly3D(obj = obj, plotting.data. = plotting.data, annotation.category = annotate.by)
  } else {
    NULL
  }

  plt <- plotly::plot_ly(
    data = plotting.data,
    x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,
    type = "scatter3d",
    mode = "markers",
    marker = list(size = dotsize),
    text = ~label,
    color = ~Expression,
    opacity = alpha
    # , colors = c('darkgrey', 'red')
    , colorscale = "Viridis"
    # , hoverinfo="text"
    , ...
  ) %>%
    plotly::layout(title = gene, scene = list(annotations = ls.ann.auto))

  SavePlotlyAsHtml(plt, category. = gene, suffix. = suffix)
  return(plt)
}




# _________________________________________________________________________________________________
#' @title plot3D.umap
#'
#' @description Plot a 3D umap based on one of the metadata columns. Uses plotly. Based on github.com/Dragonmasterx87.
#' @param category The metadata column based on which the 3D UMAP will be plotted.
#' Default: First returned result from GetNamedClusteringRuns(obj) function.
#' @param obj The Seurat object for which the 3D umap plot will be generated. Default: combined.obj
#' @param suffix A suffix added to the filename. Default: NULL
#' @param annotate.by The cluster or grouping to be used for automatic annotation.
#' Default: First returned result from GetNamedClusteringRuns(obj) function.
#' @param dotsize The size of the dots in the plot. Default: 1.25
#' @param ... Pass any other parameter to the internally called `plotly::plot_ly`.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   plot3D.umap(category = "integrated_snn_res.0.1", obj = combined.obj)
#' }
#' }
#' @importFrom plotly plot_ly layout
#' @importFrom Seurat FetchData
#'
#' @export

plot3D.umap <- function(
    category,
    obj = combined.obj,
    annotate.by = category,
    suffix = NULL,
    dotsize = 1.25,
    col.names = c("umap_1", "umap_2", "umap_3"),
    ...) {

  message("category: ", category)
  message("annotate.by: ", annotate.by)

  # browser()

  # Input assertions ____________________________________
  stopifnot(
    is(obj, "Seurat"),
    category %in% colnames(obj@meta.data),
    annotate.by %in% colnames(obj@meta.data),
    "reductions.backup is missing from @misc" = is.list(obj@misc$"reductions.backup"),
    "umap3d is missing from @misc$reductions.backup" = is(obj@misc$reductions.backup$"umap3d", class2 = "DimReduc"),
    "reductionn has 3 columns" = (ncol(obj@misc$reductions.backup$"umap3d") == 3),
    "3D reduction has >/< cells than object" = (ncol(obj) == nrow(obj@misc$reductions.backup$"umap3d"@cell.embeddings))
  )

  if (obj@version < 5) col.names <- toupper(col.names)
  message("Obj. version: ", obj@version, " \ndim names: ", kppc(col.names))

  # Get and format 3D plotting data ____________________________________
  plotting.data <- obj@misc$reductions.backup$"umap3d"@cell.embeddings # plotting.data <- Seurat::FetchData(object = obj, vars = c(col.names, category))
  colnames(plotting.data) <- toupper(col.names)

  plotting.data <- cbind(plotting.data, obj[[category]])
  colnames(plotting.data)[4] <- "category"
  plotting.data$label <- paste(rownames(plotting.data)) # Make a column of row name identities (these will be your cell/barcode names)

  ls.ann.auto <- if (annotate.by != FALSE) {
    .Annotate4Plotly3D(obj = obj, plotting.data. = plotting.data, annotation.category = annotate.by)
  } else {
    NULL
  }

  plt <- plotly::plot_ly(
    data = plotting.data,
    x = ~UMAP_1, y = ~UMAP_2, z = ~UMAP_3,
    type = "scatter3d",
    mode = "markers",
    marker = list(size = dotsize),
    text = ~label,
    color = ~category,
    colors = gg_color_hue(length(unique(plotting.data$"category")))
    # , hoverinfo="text"
    , ...
  ) %>%
    plotly::layout(title = category, scene = list(annotations = ls.ann.auto))

  SavePlotlyAsHtml(plt, category. = category, suffix. = suffix)
  return(plt)
}


# _________________________________________________________________________________________________
#' @title SavePlotlyAsHtml
#'
#' @description Save a Plotly 3D scatterplot as an HTML file.
#' @param plotly_obj The Plotly object to save.
#' @param category The category of the plot.
#' @param suffix A suffix to add to the filename.
#' @param OutputDir The output directory.
#' @seealso
#'  \code{\link[htmlwidgets]{saveWidget}}
#' @examples \dontrun{
#' plt <- plotly::plot_ly("some stuff")
#' SavePlotlyAsHtml(plt, category. = "label.categ", suffix. = "test")
#' }
#'
#' @export
#' @importFrom htmlwidgets saveWidget
SavePlotlyAsHtml <- function(plotly_obj, category. = category, suffix. = NULL) { # Save Plotly 3D scatterplot as an html file.
  OutputDir <- if (exists("OutDir")) OutDir else getwd()
  name.trunk <- kpp("umap.3D", category., suffix., idate(), "html")
  fname <- kpps(OutputDir, name.trunk)
  iprint("Plot saved as:", fname)
  htmlwidgets::saveWidget(plotly_obj, file = fname, selfcontained = TRUE, title = category.)
}


# _________________________________________________________________________________________________
#' @title Backup Dimensionality Reduction Data
#'
#' @description This function is mostly used internally.It stores a backup of specified
#' dimensionality reduction data (e.g., UMAP, tSNE, PCA)
#' within the Seurat object, from `obj@reductions$umap` to the `@misc$reductions.backup` slot. This
#' allows to store 2D and 3D UMAP visualizations in parallel and easily switch between them via
#' the `RecallReduction` function.
#'
#' @param obj Seurat object containing dimensionality reduction data; Default: `combined.obj`.
#' @param dim Number of dimensions to include in the backup; Default: 2.
#' @param reduction Type of dimensionality reduction to backup ('umap', 'tsne', 'pca'); Default: 'umap'.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   obj <- BackupReduction(obj = obj, dim = 2, reduction = "umap")
#' }
#' }
#'
#' @export
BackupReduction <- function(obj = combined.obj, dim = 2, reduction = "umap") { # Backup UMAP to `obj@misc$reductions.backup` from `obj@reductions$umap`.
  if (is.null(obj@misc$"reductions.backup")) obj@misc$"reductions.backup" <- list()
  dslot <- paste0(reduction, dim, "d")
  obj@misc$reductions.backup[[dslot]] <- obj@reductions[[reduction]]
  return(obj)
}


# _________________________________________________________________________________________________
#' @title Compute and Backup Dimensionality Reductions
#'
#' @description Executes specified dimensionality reduction (UMAP, tSNE, or PCA) over a range of dimensions
#' and backs up the results within a Seurat object. This function allows for exploration of data structure
#' at varying levels of granularity and ensures that reduction results are preserved for future reference.
#'
#' @param obj Seurat object to process; Default: `combined.obj`.
#' @param nPCs Number of principal components to consider in the reduction; Default: `p$n.PC`.
#' @param dimensions Numeric vector specifying target dimensions for the reductions; Default: `3:2`.
#' @param reduction Type of dimensionality reduction to apply ('umap', 'tsne', or 'pca'); Default: 'umap'.
#' @param ... Additional parameters to pass to the dimensionality reduction functions.
#'
#' @return Modified Seurat object with added dimensionality reduction data and backups.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   combined.obj <- SetupReductionsNtoKdimensions(
#'     obj = combined.obj, nPCs = 10,
#'     dimensions = 2:3, reduction = "umap"
#'   )
#' }
#' }
#'
#' @export
#' @importFrom Seurat RunUMAP RunTSNE RunPCA
SetupReductionsNtoKdimensions <- function(obj = combined.obj, nPCs = p$"n.PC", dimensions = 3:2,
                                          reduction = "umap", ...) {
  red <- reduction
  for (d in dimensions) {
    iprint(d, "dimensional", red, "is calculated")
    obj <- if (reduction == "umap") {
      RunUMAP(obj, dims = 1:nPCs, n.components = d, ...)
    } else if (reduction == "tsne") {
      RunTSNE(obj, dims = 1:nPCs, n.components = d, ...)
    } else if (reduction == "pca") {
      RunPCA(obj, dims = 1:nPCs, n.components = d, ...)
    }
    obj <- BackupReduction(obj = obj, dim = d, reduction = red)
  }
  return(obj)
}


# _________________________________________________________________________________________________
#' @title Recall Dimensionality Reduction from backup slot
#'
#' @description Restores dimensionality reduction data (e.g., UMAP, tSNE, PCA) from a backup
#' stored within `obj@misc$reductions.backup` to the active `obj@reductions` slot.
#'
#' @param obj Seurat object from which the backup will be restored; Default: `combined.obj`.
#' @param dim Number of dimensions of the reduction data to restore; Default: 2.
#' @param reduction Type of dimensionality reduction to be restored ('umap', 'tsne', 'pca'); Default: 'umap'.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   combined.obj <- RecallReduction(obj = combined.obj, dim = 2, reduction = "umap")
#'   qUMAP()
#'   combined.obj <- RecallReduction(obj = combined.obj, dim = 3, reduction = "umap")
#'   qUMAP()
#' }
#' }
#'
#' @export
RecallReduction <- function(obj = combined.obj, dim = 2, reduction = "umap") {
  dslot <- paste0(reduction, dim, "d")
  reduction.backup <- obj@misc$reductions.backup[[dslot]]
  msg <- paste(dim, "dimensional", reduction, "from obj@misc$reductions.backup")
  stopif(is.null(reduction.backup), message = paste0(msg, " is NOT FOUND"))
  iprint(msg, "is set active. ")
  stopifnot(dim == ncol(reduction.backup))
  obj@reductions[[reduction]] <- reduction.backup
  return(obj)
}



# _________________________________________________________________________________________________
#' @title .Annotate4Plotly3D
#'
#' @description Internal helper function. Create annotation labels for 3D plots.
#' Source https://plot.ly/r/text-and-annotations/#3d-annotations.
#' @param obj The Seurat object for which the 3D plot annotations will be generated. Default: combined.obj
#' @param plotting.data. The data frame containing plotting data.
#' @param annotation.category The category for which the annotation is generated.
#' @export
#' @importFrom dplyr group_by summarise
#' @importFrom Seurat FetchData

.Annotate4Plotly3D <- function(
    obj = combined.obj,
    plotting.data.,
    annotation.category) {
  stopifnot(
    "annotation.category is missing" = !is.null(annotation.category),
    "plotting.data. is missing" = !is.null(plotting.data.),
    "annotation.category is not in meta.data" = annotation.category %in% colnames(obj@meta.data)
  )

  plotting.data.$"annot" <- Seurat::FetchData(object = obj, vars = c(annotation.category))[, 1]
  auto_annot <-
    plotting.data. %>%
    group_by(annot) %>%
    summarise(
      showarrow = FALSE,
      xanchor = "left",
      xshift = 10,
      opacity = 0.7,
      "x" = mean(UMAP_1),
      "y" = mean(UMAP_2),
      "z" = mean(UMAP_3)
    )
  names(auto_annot)[1] <- "text"
  ls.ann.auto <- apply(auto_annot, 1, as.list)
  return(ls.ann.auto)
}

# _________________________________________________________________________________________________
#' @title Plot3D.ListOfGenes
#'
#' @description Plot and save list of 3D UMAP or tSNE plots using plotly.
#' @param obj Seurat object to be used for the plot. Default: combined.obj
#' @param annotate.by Variable to annotate the clusters by. Default: 'integrated_snn_res.0.7'
#' @param opacity Opacity for the plot points. Default: 0.5
#' @param cex Point size for the plot. Default: 1.25
#' @param default.assay Default assay to be used from the Seurat object. Default: second entry from c("integrated", "RNA")
#' @param ListOfGenes List of genes to be plotted. Default: c("BCL11B", "FEZF2", "EOMES", "DLX6-AS1", "HOPX", "DDIT4")
#' @param SubFolderName Name of the subfolder where the plots will be saved. Default: a subfolder named 'plot3D' concatenated with the list of genes.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   CellTypeMarkers <- c("PGK1", "CTIP2" = "BCL11B", "FEZF2", "EOMES", "DLX6-AS1", "HOPX", "DDIT4", "TOP2A", "PTGDS", "EDNRB", "EGFR", "SCGN", "NR2F2", "EMX2", "GAD2", "DLX2", "SATB2")
#'   Plot3D.ListOfGenes(obj = combined.obj, ListOfGenes = CellTypeMarkers)
#' }
#' }
#' @export
Plot3D.ListOfGenes <- function(
    obj = combined.obj # Plot and save list of 3D UMAP ot tSNE plots using plotly.
    , annotate.by = "integrated_snn_res.0.7", opacity = 0.5, cex = 1.25, default.assay = c("integrated", "RNA")[2],
    ListOfGenes = c("BCL11B", "FEZF2", "EOMES", "DLX6-AS1", "HOPX", "DDIT4"),
    SubFolderName = ppp("plot3D", substitute(ListOfGenes))) {
  try(create_set_SubDir(SubFolderName))
  obj. <- obj
  rm("obj")
  stopifnot(annotate.by %in% c(colnames(obj.@meta.data), FALSE))

  DefaultAssay(object = obj.) <- default.assay
  MissingGenes <- setdiff(ListOfGenes, rownames(obj.))
  if (length(MissingGenes)) iprint("These genes are not found, and omitted:", MissingGenes, ". Try to change default assay.")
  ListOfGenes <- intersect(ListOfGenes, rownames(obj.))

  for (i in 1:length(ListOfGenes)) {
    g <- ListOfGenes[i]
    print(g)
    plot3D.umap.gene(obj = obj., gene = g, annotate.by = annotate.by, alpha = opacity, def.assay = default.assay, dotsize = cex)
  }
  try(oo())
  try(create_set_Original_OutDir(NewOutDir = ParentDir))
}


# _________________________________________________________________________________________________
#' @title Plot3D.ListOfCategories
#'
#' @description This function plots and saves a list of 3D UMAP or tSNE plots using plotly.
#' @param obj A Seurat object for which the plot is to be created. Default is 'combined.obj'.
#' @param annotate.by Character vector specifying the metadata column to be used for annotating the plot. Default is 'integrated_snn_res.0.7'.
#' @param cex Numeric value specifying the point size on the plot. Default is 1.25.
#' @param default.assay Character vector specifying the assay to be used. Default is 'RNA' (second element in the vector c("integrated", "RNA")).
#' @param ListOfCategories Character vector specifying the categories to be included in the plot. Default categories are "v.project", "experiment", "Phase", "integrated_snn_res.0.7".
#' @param SubFolderName String specifying the name of the subfolder where the plots will be saved. By default, it's created using the function ppp("plot3D", substitute(ListOfCategories)).
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   categ3Dplots <- c("v.project", "experiment", "Phase", "integrated_snn_res.0.7", "Area", "Individual", "Type")
#'   Plot3D.ListOfCategories(obj = combined.obj, ListOfCategories = categ3Dplots)
#' }
#' }
#' @export
Plot3D.ListOfCategories <- function(
    obj = combined.obj # Plot and save list of 3D UMAP ot tSNE plots using plotly.
    , annotate.by = "integrated_snn_res.0.7", cex = 1.25, default.assay = c("integrated", "RNA")[2],
    ListOfCategories = c("v.project", "experiment", "Phase", "integrated_snn_res.0.7"),
    SubFolderName = ppp("plot3D", substitute(ListOfCategories))) {
  try(create_set_SubDir(SubFolderName))
  obj. <- obj
  rm("obj")
  stopifnot(annotate.by %in% colnames(obj.@meta.data))
  DefaultAssay(object = obj.) <- default.assay

  MissingCateg <- setdiff(ListOfCategories, colnames(obj.@meta.data))
  if (length(MissingCateg)) iprint("These metadata categories are not found, and omitted:", MissingCateg, ". See colnames(obj@meta.data).")
  ListOfCategories <- intersect(ListOfCategories, colnames(obj.@meta.data))

  for (i in 1:length(ListOfCategories)) {
    categ <- ListOfCategories[i]
    print(categ)
    plot3D.umap(obj = obj., category = categ, annotate.by = annotate.by, dotsize = cex)
  }
  try(oo())
  try(create_set_Original_OutDir(NewOutDir = ParentDir))
}


# _________________________________________________________________________________________________
# TEMPORARY ______________________________ ----
# _________________________________________________________________________________________________

# _________________________________________________________________________________________________
#' @title Display Correlation Values in Pairs Plot
#'
#' @description This function displays the correlation coefficient and significance level within
#' a scatterplot generated by the `pairs()` function. The default correlation method is Pearson,
#' but Kendall or Spearman methods can also be selected.
#'
#' @param x Numeric vector or the first half of the data pair.
#' @param y Numeric vector or the second half of the data pair.
#' @param digits Number of significant digits to display in the correlation coefficient.
#' Default: 2.
#' @param prefix A string prefix added before the correlation coefficient. Default: "".
#' @param cex.cor The character expansion factor for the correlation coefficient text.
#' This argument directly influences the text size. Default: 2.
#' @param method The method of correlation coefficient calculation. It can be "pearson" (default),
#' "kendall", or "spearman".
#'
#' @return This function does not return a value but modifies the current plot by adding the
#' correlation coefficient and its significance level.
#'
#' @examples
#' \dontrun{
#' pairs(mtcars[, 1:4], panel = panelCorPearson)
#' }
#' @importFrom graphics text par
#' @importFrom stats cor cor.test
#' @export

panelCorPearson <- function(x, y, digits = 2, prefix = "", cex.cor = 2, method = "pearson") {
  # Input validation
  stopifnot(
    is.numeric(x), is.numeric(y),
    is.numeric(digits) && digits > 0,
    is.character(prefix),
    is.numeric(cex.cor) && cex.cor > 0,
    method %in% c("pearson", "kendall", "spearman")
  )

  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, method = method, use = "complete.obs"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste(prefix, txt, sep = "")
  if (missing(cex.cor)) cex <- 0.8 / strwidth(txt)

  test <- cor.test(x, y, method = method)
  Signif <- symnum(test$p.value,
    corr = FALSE, na = FALSE,
    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
    symbols = c("***", "**", "*", ".", " ")
  )

  cex <- ifelse(missing(cex.cor), 0.8 / strwidth(txt), cex.cor)
  text(0.5, 0.5, txt, cex = cex * r)
  text(.8, .8, Signif, cex = cex, col = 2)
}



# _________________________________________________________________________________________________
#' @title suPlotVariableFeatures for Single Seurat Object
#'
#' @description Generates a Variable Feature Plot for a specified Seurat object, labels points with
#' the top 20 variable genes, and saves the plot to a PDF file.
#'
#' @param obj A single Seurat object.
#' @param NrVarGenes A vector containing the top 20 variable genes for the Seurat object.
#' @param sampleName A string specifying the sample name, used to generate the filename for saving
#' the plot.
#' @param ppp A function for constructing the path and filename for saving the plot. It takes three
#' arguments: a prefix for the filename, a sample name, and the file extension ('pdf').
#' @param repel A logical value indicating whether to repel the labels to avoid overlap. Default: TRUE.
#' @param plotWidth Numeric value specifying the width of the plot when saved. Default: 7.
#' @param plotHeight Numeric value specifying the height of the plot when saved. Default: 5.
#'
#' @examples
#' \dontrun{
#' suPlotVariableFeatures(combined.obj)
#' }
#' @export
suPlotVariableFeatures <- function(obj = combined.obj, NrVarGenes = 15,
                                   repel = TRUE, plotWidth = 7, plotHeight = 5, save = TRUE,
                                   # suffix = kpp("nVF", .getNrScaledFeatures(obj)),
                                   suffix = NULL,
                                   ...) {
  # Input validation
  stopifnot(
    is(obj, "Seurat"), is.function(ppp), is.logical(repel),
    is.numeric(plotWidth), is.numeric(plotHeight)
  )

  obj.name <- deparse(substitute(obj))

  plot1 <- Seurat::VariableFeaturePlot(obj) +
    theme(panel.background = element_rect(fill = "white")) +
    ggtitle(label = "Variable Genes", subtitle = kppws(obj.name, suffix))


  # Assuming LabelPoints is defined elsewhere and available for use.
  TopVarGenes <- VariableFeatures(obj)[1:NrVarGenes]
  labeledPlot <- LabelPoints(
    plot = plot1, points = TopVarGenes, repel = repel,
    xnudge = 0, ynudge = 0, max.overlaps = 15
  )
  print(labeledPlot)
  filename <- ppp("Var.genes", obj.name, suffix, idate(), "png")

  # if (save) ggplot2::ggsave(plot = labeledPlot, filename = filename, width = plotWidth, height = plotHeight)
  if (save) {
    qqSave(
      ggobj = labeledPlot,
      # title = plotname,
      fname = filename, ext = ext,
      w = plotWidth, h = plotHeight, also.pdf = FALSE
    )
  }
}





