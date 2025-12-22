# ____________________________________________________________________
# Seurat.Utils.Visualization.R ----
# ____________________________________________________________________
# file.edit("~/GitHub/Packages/Seurat.utils/R/Seurat.Utils.Visualization.R")
# file.edit("~/GitHub/Packages/Seurat.utils/R/Seurat.utils.less.used.R")
# devtools::load_all("~/GitHub/Packages/Seurat.utils")
# devtools::document("~/GitHub/Packages/Seurat.utils"); devtools::load_all("~/GitHub/Packages/Seurat.utils")


# _________________________________________________________________________________________________
#' @title Plot filtering thresholds and distributions
#'
#' @description This function plots the filtering thresholds and distributions for Seurat objects,
#' using four panels to highlight the relationship between gene- and UMI-counts, and the
#' ribosomal- and mitochondrial-content.  !! Default arguments assume that `p` is a list of
#' parameters, present in the global environment, with elements `thr.lp.mito`, `thr.hp.mito`,
#' `thr.lp.ribo`, `thr.hp.ribo`, `thr.lp.nFeature_RNA`, and `thr.hp.nFeature_RNA`.
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
#' @param subdir Subdirectory within `parentdir` where plots will be stored. Default: generated using a call to `kpp()`.
#' @param transparency Point transparency on scatter plots. Default: 0.25.
#' @param cex Size of points on scatter plots. Default: 0.75.
#' @param theme.used A `ggplot2` theme for all plots. Default: `theme_bw(base_size = 18)`.
# #' @param LabelDistFromTop Distance from top for label placement. Default: 200.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # !! Default arguments assume that `p` is a list of parameters, present in the global
#'   # environment, with elements `thr.lp.mito`, `thr.hp.mito`, `thr.lp.ribo`, `thr.hp.ribo`,
#'   # `thr.lp.nFeature_RNA`, and `thr.hp.nFeature_RNA`.
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
    par.ls = p,
    parentdir = OutDirOrig,
    suffices = names(ls.obj),
    filetype = ".jpg",
    below.mito = par.ls$"thr.lp.mito",
    above.mito = par.ls$"thr.hp.mito",
    below.ribo = par.ls$"thr.lp.ribo",
    above.ribo = par.ls$"thr.hp.ribo",
    below.nFeature_RNA = if ("quantile.thr.lp.nFeature_RNA" %in% names(par.ls)) par.ls$"quantile.thr.lp.nFeature_RNA" else par.ls$"thr.lp.nFeature_RNA",
    above.nFeature_RNA = par.ls$"thr.hp.nFeature_RNA",
    subdir = FixPlotName(
      "Filtering.plots",
      "mito", par.ls$"thr.hp.mito", par.ls$"thr.lp.mito",
      "ribo", par.ls$"thr.hp.ribo", par.ls$"thr.lp.ribo",
      "nFeature", below.nFeature_RNA, above.nFeature_RNA
    ),
    transparency = 0.25,
    cex = 0.75,
    theme.used = theme_bw(base_size = 18)
    # LabelDistFromTop = 200 # for barplot_label
    ) {
  message("Expects a list of Seurat objects, `ls.obj` with names, and a list of parameters, `par.ls` with a defined structure.")

  # Create names based on the Seurat objects, catenating "dataset" and numbers 1:n
  if (is.null(suffices)) {
    suffices <- paste0("obj_", 1:length(ls.obj))
    message("Provide suffixes unique to each dataset, ideally as names of the list of Seu objects!")
  }

  stopifnot(
    is.list(ls.obj), is.list(par.ls) | is.null(par.ls),
    is.numeric(above.nFeature_RNA), is.numeric(below.nFeature_RNA),
    (below.nFeature_RNA > above.nFeature_RNA) | below.nFeature_RNA < 1, # either an absolute feature count or a quantile
    is.numeric(above.mito), is.numeric(below.mito), below.mito > above.mito,
    is.numeric(above.ribo), is.numeric(below.ribo), below.ribo > above.ribo,
    is.character(parentdir), is.character(subdir), is.character(filetype), is.numeric(transparency), is.numeric(cex),
    is.character(suffices), length(suffices) == length(ls.obj)
  )

  MarkdownHelpers::llprint(
    "We filtered for high quality cells based on the number of genes detected [", above.nFeature_RNA, ";", below.nFeature_RNA,
    "] and the fraction of mitochondrial [", percentage_formatter(above.mito), ";", percentage_formatter(below.mito),
    "] and ribosomal [", percentage_formatter(above.ribo), ";", percentage_formatter(below.ribo), "] reads."
  )

  theme_set(theme.used)
  OutDir <- FixPath(parentdir, subdir)

  print(subdir)

  MarkdownReports::create_set_OutDir(OutDir)
  stopifnot(length(suffices) == length(ls.obj))

  Calculate_nFeature_LowPass <- below.nFeature_RNA > 0 && below.nFeature_RNA < 1 # Use quantile low pass threshold
  if (Calculate_nFeature_LowPass) qval <- below.nFeature_RNA

  for (i in 1:length(ls.obj)) {
    print(suffices[i])
    metadata_df <- ls.obj[[i]]@meta.data

    if (Calculate_nFeature_LowPass) {
      nFtr <- ls.obj[[i]]$"nFeature_RNA"
      nFtr_valid <- nFtr[nFtr > above.nFeature_RNA]                                # restrict to high-quality cells, otherwise the quantile may be Influenced by the amount of junk in the library.

      below.nFeature_RNA <- floor(quantile(nFtr_valid, probs = qval))
      message(pc_TRUE(nFtr_valid < below.nFeature_RNA, NumberAndPC = T, suffix = paste("cells below thr.", below.nFeature_RNA, "at quantile:", qval)))
      stopifnot(below.nFeature_RNA > above.nFeature_RNA)
    }

    AllMetaColumnsPresent <- all(c("nFeature_RNA", "percent.mito", "percent.ribo") %in% colnames(metadata_df))
    if (!AllMetaColumnsPresent) {
      print(c("nFeature_RNA", "percent.mito", "percent.ribo"))
      print(c("nFeature_RNA", "percent.mito", "percent.ribo") %in% colnames(metadata_df))
      print("Try to run:")
      print('objX <- addMetaFraction(obj = objX, col.name = "percent.mito", gene.symbol.pattern =  "^MT\\.|^MT-")')
      print('objX <- addMetaFraction(obj = objX, col.name = "percent.ribo", gene.symbol.pattern =  "^RPL|^RPS")')
      stop()
    }

    filt.nFeature_RNA <- (metadata_df$"nFeature_RNA" < below.nFeature_RNA & metadata_df$"nFeature_RNA" > above.nFeature_RNA)
    filt.below.mito <- (metadata_df$"percent.mito" < below.mito & metadata_df$"percent.mito" > above.mito)
    filt.below.ribo <- (metadata_df$"percent.ribo" < below.ribo & metadata_df$"percent.ribo" > above.ribo)

    metadata_df <- cbind(metadata_df, filt.nFeature_RNA, filt.below.mito, filt.below.ribo)

    # Define colour thresholds for nFeature_RNA
    metadata_df$colour.thr.nFeature <- cut(metadata_df$"nFeature_RNA",
      breaks = c(-Inf, above.nFeature_RNA, below.nFeature_RNA, Inf),
      labels = c(
        paste0("LQ (<", above.nFeature_RNA, ")"),
        paste0("HQ (", above.nFeature_RNA, "< X <", below.nFeature_RNA, ")"),
        paste0("Dbl/Outlier (>", below.nFeature_RNA, ")")
      )
    )

    boolean_LC_cells <- metadata_df$"nFeature_RNA" <= above.nFeature_RNA
    LQ <- pc_TRUE(boolean_LC_cells)
    Doublets <- pc_TRUE(metadata_df$"nFeature_RNA"[!boolean_LC_cells] >= below.nFeature_RNA)

    A <- ggplot(data = metadata_df, aes(x = nFeature_RNA, fill = colour.thr.nFeature)) +
      geom_histogram(binwidth = 100) +
      ggtitle(paste(
        "Cells between", above.nFeature_RNA, "and", below.nFeature_RNA,
        "UMIs are selected \n(", pc_TRUE(filt.nFeature_RNA), "), with",
        LQ, "low-quality and", Doublets, "doublet cells excluded."
      )) +
      scale_y_log10() + annotation_logticks() +
      geom_vline(xintercept = below.nFeature_RNA) + geom_vline(xintercept = above.nFeature_RNA) +
      theme(legend.position = "none") # "top"
    # A

    B <- ggplot2::ggplot(metadata_df, aes(x = nFeature_RNA, y = percent.mito)) +
      ggplot2::ggtitle(paste(
        "Cells below", percentage_formatter(below.mito),
        "mito reads are selected \n(with A:", pc_TRUE(filt.nFeature_RNA & filt.below.mito), ")"
      )) +
      ggplot2::geom_point(
        alpha = transparency, size = cex, show.legend = FALSE,
        aes(color = filt.nFeature_RNA & filt.below.mito)
      ) +
      scale_x_log10() + annotation_logticks() +
      geom_hline(yintercept = below.mito) + geom_hline(yintercept = above.mito) +
      geom_vline(xintercept = below.nFeature_RNA) + geom_vline(xintercept = above.nFeature_RNA)
    # B

    C <- ggplot(metadata_df, aes(x = nFeature_RNA, y = percent.ribo)) +
      ggtitle(paste(
        "Cells below", percentage_formatter(below.ribo),
        "ribo reads are selected \n(with A:",
        pc_TRUE(filt.nFeature_RNA & filt.below.ribo), ")"
      )) +
      geom_point(
        alpha = transparency, size = cex, show.legend = FALSE,
        aes(color = filt.nFeature_RNA & filt.below.ribo)
      ) +
      scale_x_log10() + annotation_logticks() +
      geom_hline(yintercept = below.ribo) + geom_hline(yintercept = above.ribo) +
      geom_vline(xintercept = below.nFeature_RNA) + geom_vline(xintercept = above.nFeature_RNA)
    # C

    D <- ggplot(metadata_df, aes(x = percent.ribo, y = percent.mito)) +
      ggtitle(paste(
        "Final: All cells w/o extreme values are selected \n(with A,B,C:",
        pc_TRUE(filt.nFeature_RNA & filt.below.mito & filt.below.ribo), ")"
      )) +
      geom_point(
        alpha = transparency, size = cex, show.legend = FALSE,
        aes(color = filt.nFeature_RNA & filt.below.mito & filt.below.ribo)
      ) +
      scale_x_log10() + scale_y_log10() + annotation_logticks() +
      geom_hline(yintercept = below.mito) + geom_hline(yintercept = above.mito) +
      geom_vline(xintercept = below.ribo) + geom_vline(xintercept = above.ribo)
    # D

    # Add title to A and caption to D
    main_title <- paste("Object", suffices[i])
    caption_text <- paste0(
      "pct.mito [", below.mito, "-", above.mito, "] | ",
      "pct.ribo [", below.ribo, "-", above.ribo, "] | ",
      "nFeature_RNA [", above.nFeature_RNA, "-", below.nFeature_RNA, "] | ",
      "Doublet % and quantile cutoff is calculated on cells above min. nFeature_RNA"
    )


    # Combine plots in 2x2 grid
    p_grid <- cowplot::plot_grid(
      plotlist = list(A, B, C, D),
      nrow = 2, ncol = 2,
      labels = LETTERS[1:4],
      label_size = 20
    )

    # Add overall title and caption cleanly, with white background
    px <- cowplot::ggdraw() +
      theme(plot.background = element_rect(fill = "white", colour = NA)) +          # white canvas
      cowplot::draw_label(main_title, x = 0.02, y = 0.98, hjust = 0, vjust = 1,
                          fontface = "bold", size = 22) +                           # slightly closer to top
      cowplot::draw_plot(p_grid, y = 0.055, height = 0.89) +                        # raise the grid, taller height
      cowplot::draw_label(caption_text, x = 0.98, y = 0.02, hjust = 1, vjust = 0,
                          fontface = "italic", size = 12)                           # closer to bottom


    # Save figure
    fname <- kpps(OutDir, FixPlotName("Filtering.thresholds", suffices[i], filetype))
    cowplot::save_plot(filename = fname, plot = px, base_height = 15)
    stopifnot(file.exists(fname))
  } # for
  # _________________________________________________________________________________________________
  create_set_OutDir(parentdir)
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
# try (source('/GitHub/Packages/CodeAndRoll/CodeAndRoll.R'),silent= FALSE) # generic utilities functions
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
  stopifnot(
    inherits(obj, "Seurat"),
    "pca" %in% names(obj@reductions)
  )
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
#' @param caption A caption for the plot. Default: "hline at 1%".
#' @param use.MarkdownReports Boolean indicating whether to use `MarkdownReports` for plotting.
#' If `FALSE`, `ggExpress` is used. Default: `FALSE`.
#' @param ... Additional arguments to be passed to `ggExpress::qbarplot` or `MarkdownReports::wbarplot`.
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
                                  caption = "hline at 1%",
                                  # caption = .parseKeyParams(obj, suffix = "| hline at 1%"),
                                  use.MarkdownReports = FALSE,
                                  ...) {
  stopifnot(
    inherits(obj, "Seurat"),
    "pca" %in% names(obj@reductions)
  )
  message(" > Running scPlotPCAvarExplained...")

  pct <- scCalcPCAVarExplained(obj)
  if (use.MarkdownReports) {
    MarkdownReports::wbarplot(pct, xlab = "Principal Components", ylab = "% of variation explained", ...)
    barplot_label(round(pct, digits = 2), barplotted_variable = pct, cex = .5)
  } else {
    ggExpress::qbarplot(
      vec = pct, plotname = plotname, subtitle = sub,
      xlab = "Principal Components", ylab = "% of variation explained",
      w = 10, h = 5, hline = 1, caption = caption, ...
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
#' their expression as a percentage of the total UMIs. Default: 25.
#' @param width.barplot The width of the barplot that visualizes the highest expressed genes.
#' Default: a quarter of `n.genes.barplot`.
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
    plotname = "Gene expression as fraction of all transcripts *UMI's)",
    subtitle = "Percentage in RNA-counts",
    xlab = "Percent in Transcriptome (total per gene)",
    ylab = "Number of genes",
    xlab.angle = 45,
    w = 7, h = 5,
    ...
  )

  Highest.Expressed.Genes <- head(iround(100 * relative.total.Expr), n = n.genes.barplot)
  qbarplot(Highest.Expressed.Genes,
    plotname = "Percentage of highest expressed genes",
    subtitle = "Total, in RNA-counts",
    xlab = "",
    ylab = "Gene expression as percent of all UMI's",
    xlab.angle = 45,
    w = width.barplot, h = 5,
    ...
  )

  message("!!! \nTotalReadFraction is now stored under combined.obj@misc$'TotalReadFraction'.")

  obj@misc$"TotalReadFraction" <- relative.total.Expr
  return(obj)
}



# _________________________________________________________________________________________________
#' @title Histogram All Genes' Expression Level and a Highlighted Gene
#'
#' @description Shows a comparison of the expression level of the chosen gene to all genes.
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
#' @param layerX Data slot to use ('data' or 'counts'); Default: "data".
#' @param thr_expr Expression threshold for highlighting in the plot; Default: 10.
#' @param suffix Additional text to append to the plot title; Default: NULL.
#' @param prefix Additional text to prepend to the plot title; Default: NULL.
#' @param xlab Label for the x-axis; Default: "log10(Summed UMI count @data)".
#' @param return_cells_passing If TRUE, returns count of cells exceeding the expression threshold; Default: `TRUE`.
#' @param clip_count_qtl_thr Quantile threshold for clipping if using count data; Default: 0.95.
#' Needed for visualization (to avoid x axis compression).
#' @param log10_counts If TRUE, log10-transforms the COUNT expression values; Default: `TRUE`.
#' @param return_quantile If TRUE, returns cell count exceeding the quantile threshold; Default: `FALSE`.
#' @param w Width of the plot in inches; Default: 9.
#' @param h Height of the plot in inches; Default: 5.
#' @param show_plot If TRUE, displays the generated plot; Default: `TRUE`.
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
    assay = "RNA", layerX = "data",
    thr_expr = 10,
    suffix = NULL,
    prefix = NULL,
    plotname = c("Summed Gene-set Expression -", "Expression of"),
    xlab = paste0("Expression -log10(Summed UMIs @", layerX, ")"),
    return_cells_passing = TRUE,
    clip_count_qtl_thr = 0.99,
    log10_counts = TRUE,
    return_quantile,
    w = 9, h = 5,
    show_plot = TRUE,
    ...) {
  #
  stopifnot(
    length(genes) > 0,
    layerX %in% c("data", "counts")
  )

  # Aggregate genes if necessary
  aggregate <- length(genes) > 1
  SummedExpressionPerCell <- colSums(LayerData(
    object = obj, assay = assay,
    layer = layerX
  )[genes, , drop = FALSE])

  # Clip counts if necessary
  if (layerX == "counts") {
    SummedExpressionPerCell <- CodeAndRoll2::clip.at.fixed.value(
      x = SummedExpressionPerCell,
      thr = quantile(SummedExpressionPerCell, probs = clip_count_qtl_thr)
    )
    if (log10_counts) SummedExpressionPerCell <- log10(SummedExpressionPerCell + 1)
  }

  # Create annotation
  CPT <- paste("layer:", layerX, "| assay:", assay, "| cutoff at", iround(thr_expr))

  # Add a subtitle with the number of genes and the expression threshold
  SUBT <- filter_HP(SummedExpressionPerCell, threshold = thr_expr, return_conclusion = TRUE, plot.hist = FALSE)

  if (aggregate) {
    SUBT <- paste0(SUBT, "\n", length(genes), " genes summed up, e.g: ", kppc(head(genes)))
    TTL <- kppd(prefix, plotname[1], suffix)
  } else {
    TTL <- trimws(paste(prefix, plotname[length(plotname)], paste(genes), suffix))
  }

  # Create the plot
  pobj <- ggExpress::qhistogram(SummedExpressionPerCell,
    plotname = TTL,
    subtitle = SUBT,
    caption = CPT,
    prefix = prefix,
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
#' @title Compute % of Cells Above a Threshold for a Metadata or Gene Feature
#'
#' @description
#' Computes the fraction of cells with values above a specified threshold for
#' either a metadata column or an assay feature. Supports optional subsetting
#' and optional regrouping into boxplot-style categories.
#'
#' @param object Seurat object.
#' @param feature Character. Metadata column or gene/feature name.
#' @param ident Character. Metadata column used for grouping.
#' @param threshold Numeric. Threshold above which cells are counted.
#' @param box Logical. If TRUE, regroup by \code{ident.box}.
#' @param ident.box Character or NULL. Secondary grouping variable for box mode.
#' @param subset_ident Character or NULL. Metadata column used for subsetting.
#' @param subset_values Vector or NULL. Values of \code{subset_ident} to keep.
#' @param omit.na Logical. Whether to remove NA values in the feature.
#' @param assay Character. Assay name passed to \code{FetchData}.
#' @param slot Character. Slot name passed to \code{FetchData}.
#' @param plot Logical. Whether to generate a plot.
#' @param caption Character or NULL. Caption for the plot.
#' @param ylab Character. Y-axis label.
#' @param ... Additional arguments passed to plotting functions.
#'
#' @return
#' A named numeric vector of percentages, or a list of such vectors if
#' \code{box = TRUE}.
#'
#' @export

PctCellsAboveX <- function(
    object,
    feature,
    ident,
    threshold = 1,
    box = FALSE,
    ident.box = NULL,
    subset_ident = NULL,
    subset_values = NULL,
    omit.na = TRUE,
    assay = "RNA",
    slot = "data",
    plot = TRUE,
    caption = NULL,
    ylab = "% cells above threshold",
    palette = "jco",
    ...
) {

  # 1. Validate input ________________________________________
  stopifnot(
    inherits(object, "Seurat"),
    is.character(feature),
    is.character(ident),
    ident %in% colnames(object@meta.data),
    is.null(subset_ident) || subset_ident %in% colnames(object@meta.data),
    box || is.null(ident.box)
  )

  # 2. Fetch metadata/gene values _____________________________
  vars_to_get <- unique(c(feature, ident, subset_ident, ident.box))
  df <- Seurat::FetchData(
    object = object,
    assay = assay,
    slot = slot,
    vars = vars_to_get
  )
  df$expr <- df[[feature]]

  # 3. Subset rows if required ________________________________
  if (!is.null(subset_ident)) {
    keep_rows <- df[[subset_ident]] %in% subset_values
    df <- df[keep_rows, , drop = FALSE]
  }

  # 4. Remove NA values _______________________________________
  if (omit.na) df <- df[!is.na(df$expr), , drop = FALSE]

  # 5. Split by ident or ident.box ____________________________
  split_col <- if (box) ident.box else ident
  ls_feat <- split(df$expr, df[[split_col]])

  # 6. Compute percentages ____________________________________
  pct_vec <- vapply(ls_feat, function(x) mean(x > threshold), numeric(1))
  # browser()


  # 7. Re-group for box mode __________________________________
  if (box) {
    from_to <- split(df[[ident.box]], df[[ident]])
    from_to <- list.2.replicated.name.vec(from_to)
    pct_vec <- pct_vec[names(from_to)]
    pct_list <- split(pct_vec, f = from_to)
  }

  # 8. Plotting (your exact original code) ____________________
  if (plot) {
    if (is.null(caption)) {
      caption <- pc_TRUE(is.na(pct_vec),
        suffix = "of idents yielded NA/NaN & excluded from plot."
      )
    }

    TTL <- paste("Percentage of Cells Above Threshold for", feature)
    STL <- paste("Cells above threshold for", feature, "above", threshold)
    SFX <- ppp(feature, "by", ident, "thr", threshold, "subset_ident", subset_ident)

    if (box) {
      pobj <- qboxplot(
        pct_list,
        add = "dotplot",
        xlab.angle = 45,
        hide.legend = TRUE,
        plotname = TTL,
        subtitle = STL,
        caption = caption,
        suffix = SFX,
        ylab = ylab,
        palette_use = palette,
        ...
      )
    } else {
      pobj <- qbarplot(
        pct_vec,
        label = percentage_formatter(pct_vec),
        plotname = TTL,
        subtitle = STL,
        caption = caption,
        suffix = SFX,
        ylab = ylab,
        palette_use = palette,
        ...
      )
    }

    print(pobj)
  }

  # 9. Return final output ____________________________________
  if (box) return(pct_list)
  return(pct_vec)
}




# _________________________________________________________________________________________________
#' @title PctCellsExpressingGenes
#'
#' @description Calculates the proportion of cells expressing one or more specified genes using a Seurat
#' object as input.
#'
#' @param genes A character vector specifying the genes of interest. Must be a non-empty character vector.
#' @param obj A Seurat object containing single-cell data.
#' @param assay The assay to use for expression data. Default: "RNA".
#' @param min.expr The minimum expression level to consider a gene as "expressed". Default: 1.
#' @param ident A categorical variable from the metadata of the Seurat object. If NULL, returns overall
#' proportions. Default: NULL.
#' @param max.idents Maximum number of unique values allowed in the `ident` variable. Default: 100.
#'
#' @return A named vector if `ident` is NULL, containing the proportion of cells co-expressing all genes
#' (AND), the proportion expressing any gene (OR), and the proportion expressing each gene individually.
#' If `ident` is provided, returns a matrix with rows representing categories and columns representing
#' expression proportions.
#'
#' @examples
#' \dontrun{
#' # Load the Seurat object (example)
#' library(Seurat)
#' `combined.obj <- readRDS("path/to/your/seurat_object.rds")`
#'
#' # Define genes of interest
#' # Define genes of interest
#' genes <- c("TOP2A", "MAP2")
#' # Call the function
#' PctCellsExpressingGenes(genes = genes, obj = combined.obj)
#' # Call the function with ident
#' #' PctCellsExpressingGenes(genes = genes, obj = combined.obj, ident = "cluster")
#' }
#'
#' @importFrom Seurat GetAssayData
#' @export

PctCellsExpressingGenes <- function(genes, obj, assay = "RNA", min.expr = 1,
                                    ident = NULL, max.idents = 100) {
  # Input assertions
  stopifnot(
    is.character(genes) && length(genes) > 0, # genes must be a non-empty character vector
    inherits(obj, "Seurat"), # obj must be a Seurat object
    is.character(assay) && length(assay) == 1, # assay must be a single character string
    is.numeric(min.expr) && length(min.expr) == 1, # min.expr must be a single numeric value
    is.null(ident) || (is.character(ident) && length(ident) == 1), # ident must be NULL or a single character string
    is.numeric(max.idents) && length(max.idents) == 1 && max.idents > 0 # max.idents must be a single positive numeric value
  )

  # Message parameters to console
  message("Parameters:")
  message("  genes: ", paste(genes, collapse = ", "))
  message("  assay: ", assay)
  message("  min.expr: ", min.expr)
  message("  ident: ", ifelse(is.null(ident), "NULL", paste(ident, length(ident), "-", head(ident))))
  message("  max.idents: ", max.idents)

  # Get the expression data
  expr.data <- Seurat::GetAssayData(obj, assay = assay, slot = "data")

  # Check if the genes are in the expression data
  genes <- intersect(genes, rownames(expr.data))
  if (length(genes) == 0) {
    stop("None of the specified genes are present in the expression data.")
  }

  # Define a function to calculate proportions
  calc_proportions <- function(expr.data, genes, min.expr) {
    # Calculate the proportion of cells expressing each gene
    expr.prop <- sapply(genes, function(gene) {
      sum(expr.data[gene, ] >= min.expr) / ncol(expr.data)
    })

    # Calculate the proportion of cells co-expressing all genes (AND)
    coexpr.prop <- sum(apply(expr.data[genes, ] >= min.expr, 2, all)) / ncol(expr.data)

    # Calculate the proportion of cells expressing any gene (OR)
    orexpr.prop <- sum(apply(expr.data[genes, ] >= min.expr, 2, any)) / ncol(expr.data)

    # Return the proportions
    return(c(coexpr.prop, orexpr.prop, expr.prop))
  }

  # Calculate proportions
  proportions <- calc_proportions(expr.data, genes, min.expr)

  # Message to console
  message(sprintf("Percentage of cells co-expressing all genes (AND): %.2f%%", proportions[1] * 100))
  message(sprintf("Percentage of cells expressing any gene (OR): %.2f%%", proportions[2] * 100))
  for (i in seq_along(genes)) {
    message("gene: ", genes[i], " ...")
    message(sprintf("Percentage of cells expressing %s: %.2f%%", genes[i], proportions[2 + i] * 100))
  }

  # If ident is NULL, return the proportions vector
  if (is.null(ident)) {
    names(proportions) <- c("Expr.ALL", "Expr.ANY", genes)
    return(proportions)
  }

  # Check ident
  ident_vals <- obj@meta.data[[ident]]
  if (is.null(ident_vals)) {
    stop(sprintf("The ident '%s' is not present in the metadata.", ident))
  }

  if (length(unique(ident_vals)) > max.idents) {
    stop(sprintf("The number of unique values in ident '%s' exceeds max.idents.", ident))
  }

  # Message ident details
  message("Ident details:")
  message("  Length of idents: ", length(unique(ident_vals)))
  message("  Head of ident values: ", paste(head(unique(ident_vals)), collapse = ", "))

  # Calculate proportions per ident
  idents <- unique(ident_vals)
  result_matrix <- matrix(NA, nrow = length(idents), ncol = length(proportions))
  rownames(result_matrix) <- idents
  colnames(result_matrix) <- c("Expr.ALL", "Expr.ANY", genes)

  for (ident_value in idents) {
    message("Cluster: ", ident_value)
    cells_in_ident <- which(ident_vals == ident_value)
    expr.data_subset <- expr.data[, cells_in_ident, drop = FALSE]
    result_matrix[ident_value, ] <- calc_proportions(expr.data_subset, genes, min.expr)
  }

  return(result_matrix)
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
#' @param fill.by The variable to fill by for the bar plot.
#' @param group.by The variable to group by for the bar plot.
#' @param obj A Seurat object.
#' @param plotname The title of the plot.
#' @param min.nr.sampled.cells The minimal number of cells to sample from each identity class. Defaults to 200 cells.
#' @param downsample Logical indicating whether to downsample data to equalize group sizes.
#' @param prefix Optional prefix for the plot title.
#' @param suffix Optional suffix for the plot title.
#' @param sub_title Optional subtitle for the plot.
#' @param hlines Numeric vector specifying y-intercepts of horizontal lines to add to the plot.
#' @param return_table Logical; if TRUE, returns a contingency table instead of plotting.
#' @param save_table Logical; if TRUE, saves the table behind the plot.
#' @param save_plot Logical; if TRUE, saves the generated plot.
#' @param also.pdf Save plot in both png and pdf formats.
#' @param seedNr Seed for random number generation to ensure reproducibility.
#' @param draw_plot Logical; if FALSE, suppresses plotting (useful if only the table is desired).
#' @param show_numbers Logical; if TRUE, adds count numbers on top of each bar in the plot.
#' @param min.pct Show % Labels above this threshold. Default = 0.05, or above 5 pct.
#' @param cex.pct Font size of pct labels.
#' @param min_frequency Minimum fraction to display individually in the plot; smaller fractions
#' are aggregated into an "Other" category.
#' @param custom_col_palette Specifies whether to use a standard or custom color palette.
#' @param color_scale Defines the color scale to use for the plot if a custom palette is selected.
#' @param show.total.cells Show total cells
#' @param cex.total Label size for total cells
#' @param xlab.angle Angle of x-axis labels.
#' @param show_plot Logical; if TRUE, shows the plot.
#' @param w Width of the plot in inches. Default: `NULL`
#' @param h Height of the plot in inches. Default: `6`
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
    plotname = kppws("Cell proportions of", fill.by, "by", group.by),
    suffix = NULL,
    prefix = NULL,
    sub_title = suffix,
    hlines = c(.25, .5, .75),
    return_table = FALSE,
    save_table = TRUE,
    save_plot = TRUE,
    also.pdf = FALSE,
    seedNr = 1989,
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
    show_plot = TRUE,
    w = NULL,
    h = 6,
    ...) {
  # Input assertions
  stopifnot(
    inherits(obj, "Seurat"), # obj must be a Seurat object
    is.numeric(min_frequency) && length(min_frequency) == 1 && min_frequency >= 0 && min_frequency < 1, # min_frequency must be between 0 and 1
    group.by %in% colnames(obj@meta.data), # group.by must be a valid column in the meta.data slot of the Seurat object
    fill.by %in% colnames(obj@meta.data), # fill.by must be a valid column in the meta.data slot of the Seurat object
    "To many categories for X axis (group.by)" = nr.unique(obj@meta.data[, group.by]) < 100
  )

  META <- obj@meta.data

  if (is.null(w)) {
    categ_X <- nr.unique(META[, group.by])
    categ_Y <- nr.unique(META[, fill.by])
    w <- ceiling(max(7, categ_Y / 4, categ_X / 1.5))
  }

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
      "\nDownsampled all groups in ", fill.by, " (Y) to ", n_smallest_group,
      " cells before splitting by X. \nIt is calculated as max(smallest group, 5% of total cells). Largest group previously contained: ", largest_grp
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
    totals <- META |>
      group_by(!!sym(group.by)) |>
      summarise(Total = n()) |>
      ungroup()

    # Merge totals back with the original data for labeling
    group_by_column <- group.by
    META <- META |>
      left_join(totals, by = setNames(nm = group_by_column, group_by_column))
  }


  if (draw_plot) {
    # calculate the proportions and add up small fractions
    prop_table <- META |>
      group_by(!!as.name(fill.by)) |>
      summarise(proportion = n() / nrow(META)) |>
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
        META |>
        group_by(!!sym(fill.by)) |>
        sample_n(
          size = max(n_smallest_group, min.nr.sampled.cells),
          replace = dsample.to.repl.thr
        ) |>
        ungroup()

      contingency.table <- table(META[[group.by]], META[[fill.by]])
      contingency.table <- addmargins(contingency.table)
      print(contingency.table)
    }

    # Plot the data
    pl <- META |>
      group_by(!!sym(group.by)) |>
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
    if (!isFALSE(custom_col_palette)) {
      stopifnot("Length(custom_col_palette) should be >= nr. categories displayed." = length(custom_col_palette) >= n.categories)

      all_categs_have_a_col <- all(categories %in% names(custom_col_palette))
      if (all_categs_have_a_col) {
        colz_manual <- custom_col_palette[categories]
      } else {
        colz_manual <- custom_col_palette[seq(categories)]
      } # end if all_categs_have_a_col
      pl <- pl + scale_fill_manual(values = colz_manual)
    } else if (rnd_colors) {
      colz_manual <- sample(rainbow(n.categories))
      pl <- pl + scale_fill_manual(values = colz_manual)
    } # end if custom_col_palette / rnd_colors



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
      # sfx <- shorten_clustering_names(group.by)
      sfx <- if (!is.null(suffix)) suffix else NULL
      if (min_frequency) sfx <- sppp(sfx, min_frequency)
      qqSave(
        ggobj = pl, title = FixPlotName(plotname), also.pdf = also.pdf, w = w, h = h,
        suffix = sppp(sfx, "fr.barplot")
        # , ...
      )
    } # save_plot
  } # draw_plot

  # Compile contingency table and its frequencies
  CT_freq_sc <- list(
    "values" = contingency.table,
    "percentages" = CodeAndRoll2::rowDivide(mat = contingency.table, vec = rowSums(contingency.table))
  )

  if (save_table) {
    ReadWriter::write.simple.xlsx(CT_freq_sc,
      filename = sppp(FixPlotName(plotname), suffix, "fr.barplot")
      # suffix = sppp(FixPlotName(plotname), "fr.barplot")
    )
  }

  # Return contingency table or plot based on return_table flag
  if (show_plot) print(pl)

  if (return_table) {
    return(CT_freq_sc)
  } else {
    invisible(pl)
  } # end if return_table
}





# _________________________________________________________________________________________________
#' @title Barplot of Fraction of Cells per Cluster
#'
#' @description Visualizes the fraction of cells within each cluster through a barplot.
#'
#' @param obj Seurat object for analysis. Default: `combined.obj`.
#' @param ident Cluster identity. Used to specify which clustering results to visualize.
#' Default: First entry from ordered clustering runs.
#' @param sort If TRUE, sorts clusters by size. Default: `FALSE`.
#' @param title Title for the plot. Default: "Cells per Identity Group".
#' @param sub Subtitle for the plot. Default: "identity".
#' @param label If TRUE, shows cell count or percentage based on the label vector. Default: `TRUE`.
#' @param palette Color palette for the barplot. Default: 'glasbey'.
#' @param return_table If TRUE, returns the data used for plotting instead of the plot itself. Default: `FALSE`.
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
    plotname = "Cells per Identity Group",
    sub = ident,
    label = list(TRUE, "percent")[[1]],
    suffix = if (label == "percent") "percent" else NULL,
    col = NULL,
    palette = c("alphabet", "alphabet2", "glasbey", "polychrome", "stepped")[3],
    return_table = FALSE,
    ylab_adj = 1.1,
    min.cells = round(ncol(obj) / 100),
    ...) {
  #
  stopifnot(
    inherits(obj, "Seurat"), is.character(ident), is.logical(sort), is.character(plotname), is.character(sub),
    is.logical(label) | label == "percent", is.character(suffix) | is.null(suffix), is.character(palette), is.logical(return_table),
    is.numeric(ylab_adj), is.numeric(min.cells), ident %in% colnames(obj@meta.data)
  )

  1
  cell.per.cl <- obj[[ident]][, 1]
  cell.per.cluster <- (table(cell.per.cl, useNA = "ifany"))
  if (sort) cell.per.cluster <- sort(cell.per.cluster)
  lbl <- if (isFALSE(label)) {
    NULL
  } else if (label == "percent") {
    percentage_formatter(cell.per.cluster / sum(cell.per.cluster), digitz = 2)
  } else if (isTRUE(label)) {
    cell.per.cluster
  } else {
    label
  }

  min.PCT.cells <- min.cells / ncol(obj)
  message("min cell thr: ", min.cells, " corresponding to min: ", percentage_formatter(min.PCT.cells))

  n.clusters <- length(cell.per.cluster)
  nr.cells.per.cl <- table(obj[[ident]][, 1], useNA = "ifany")

  SBT <- pc_TRUE(nr.cells.per.cl < min.cells,
    NumberAndPC = TRUE,
    suffix = paste("of identities are below:", min.cells, "cells, or", percentage_formatter(min.PCT.cells), "of all cells.")
  )

  color <- if (is.null(col)) 1:n.clusters else col

  # Fix NA names, if any
  names(cell.per.cluster)[is.na(names(cell.per.cluster))] <- "NA"


  pl <- ggExpress::qbarplot(cell.per.cluster,
    plotname = plotname,
    subtitle = paste0(sub, "\n", SBT),
    suffix = kpp(ident, ncol(obj), "c", suffix),
    col = color,
    caption = .parseBasicObjStats(obj = obj),
    xlab.angle = 45,
    ylim = c(0, ylab_adj * max(cell.per.cluster)),
    label = lbl,
    ylab = "Cells",
    palette_use = DiscretePaletteSafe(n = n.clusters, palette.used = palette),
    ...
  )

  if (return_table) {
    print(pl)
    return(cell.per.cluster)
  } else {
    return(pl)
  }
}

# _________________________________________________________________________________________________

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
#' @param return.df Whether to return the underlying data frame instead of the plot. Default: `FALSE`.
#' @param label Whether to add labels to the bar plot. Default: NULL.
#' @param subtitle Optional subtitle for the plot.
#' @param ext File extension for saving the plot. Default: '.png'.
#' @param suffix Suffix for the output file name.
#' @param above Whether to calculate the fraction of cells above or below the threshold. Default: `TRUE`.
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
    ext = ".png",
    return.df = FALSE,
    label = NULL,
    suffix = NULL,
    above = TRUE,
    ylim = c(0, 100), # set to null for relative y axis
    ...) {
  stopifnot(value.col %in% colnames(obj@meta.data))

  meta <- obj@meta.data
  metacol <- meta |>
    dplyr::select(c(id.col, value.col))

  (df_cells_above <- metacol |>
    dplyr::group_by(!!sym(id.col)) |>
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
    filename = FixPlotName(kpp(pname, id.col, ext)),
    suffix = suffix,
    subtitle = subtitle,
    caption = paste(
      "Overall average (black line):", iround(total_average), "% |",
      substitute_deparse(obj),
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
#' @param return.df If TRUE, returns the data frame instead of the plot. Default: `FALSE`.
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
  #
  scBarplot.FractionAboveThr(
    thrX = thrX,
    value.col = value.col,
    id.col = id.col,
    obj = obj,
    return.df = return.df,
    subtitle = subtitle,
    suffix = suffix,
    above = FALSE,  # Set `above` to FALSE to get fraction below threshold
    ...
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
#' scPieClusterDistribution(obj = combined.obj, ident = "cluster_identity")
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
#' from `ls.obj`. Default: `FALSE`.
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

  cellCounts <- sapply(ls.obj, ncol)
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
#' @param ... Additional arguments passed to `ggExpress::qbarplot.df`.

#'
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
  df <- df |>
    dplyr::group_by(Sample, Category) |>
    dplyr::summarise(Cells = n(), .groups = "drop") |>
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
#' @param show.colors Whether to display the colors in the palette, Default: `FALSE`.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   getDiscretePalette()
#' }
#' }
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
#' @param show.colors If TRUE, displays the generated colors. Default: `FALSE`.
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
                                  palette.used = c("alphabet", "alphabet2", "glasbey", "polychrome", "parade", "stepped")[2],
                                  show.colors = FALSE,
                                  seed = 1989) {
  stopifnot(
    is.character(ident.used), is(obj, "Seurat"),
    is.character(palette.used), is.logical(show.colors), is.numeric(seed)
  )

  categs <- unique(unlist(obj[[ident.used]]))
  n.clusters <- length(categs)

  colorz <- DiscretePaletteSafe(
    n = n.clusters,
    palette.used = palette.used,
    show.colors = show.colors,
    seed = seed
  )
  names(colorz) <- categs

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
#' @param show.colors If TRUE, displays the generated color palette. Default: `FALSE`.
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
                                palette.used = c("alphabet", "alphabet2", "glasbey", "polychrome", "parade", "stepped")[2],
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
#' `DiscretePalette` function. Default: `TRUE`.
#' @param palette Name of the color palette to use if `use_new_palettes` is TRUE.
#' Options: "alphabet", "alphabet2", "glasbey", "polychrome", "stepped". Default: "glasbey".
#' @param ident Clustering identity to use for coloring. Retrieved from the first entry
#' of `GetClusteringRuns()` by default.
#' @param show If TRUE, displays a plot showing the color mapping for each cluster. Default: `TRUE`.
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
#' Default: `FALSE`.
#' @param simple If TRUE, returns only the unique set of colors used.
#' If FALSE, returns a named vector mapping cluster identities to colors.
#' Default: `FALSE`.
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
  if (plot.colors) MarkdownHelpers::color_check(colorlevels)
  if (simple) {
    colorlevels
  } else {
    CodeAndRoll2::translate(
      vec = as.character(ident.vec),
      old = levels(ident.vec),
      new = colorlevels
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
#' results <- calculateAverageMetaData(obj = combined.obj)
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
                                rowname_column = 1,
                                ...) {
  stopifnot(is.list(results), is.character(file.prefix), is.character(path))

  for (mt in names(results)) {
    res <- results[[mt]]
    stopifnot(
      !anyNA(res[[rowname_column]]),
      !anyNaN(res[[rowname_column]])
    )

    # Generate heatmap plot
    x <- ReadWriter::column.2.row.names(results[[mt]], rowname_column = rowname_column)
    pobj <- pheatmap::pheatmap(
      mat = x,
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
    MarkdownReports::wplot_save_pheatmap(
      x = pobj, data = x, plotname = file_name,
      png = TRUE, pdf = FALSE, ...
    )
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

  # fname <- kpp("FeatureScatter", plotname)
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
#' @param feature A character string specifying the name of the feature to plot.
#' @param ident A character vector specifying the identities to be used in the plot.
#' @param split.by A character string specifying the grouping variable for splitting the plot.
#' @param colors A character vector specifying the colors to use for the plot.
#' @param clip.outliers A logical indicating whether to clip outliers.
#' @param replace.na A logical indicating whether NA values should be replaced.
#' @param pt.size The size of the individual datapoints in the plot. Default: 0 for simple violin plot.
#' @param sub Subtitle of the plot. Default: feature by ident.
#' @param suffix An optional string to append to the title of the plot.
#' @param suffix.2.title A logical indicating whether to append the suffix to the plot title.
#' @param logY A logical indicating whether to use a logarithmic scale for the y-axis.
#' @param hline A numeric or logical value; if numeric, the value where a horizontal line should be drawn.
#' @param caption A character string or logical for the plot caption. If FALSE, no caption is displayed.
#' @param ylab Y-axis label. Default: "Expression".
#' @param ylimit A numeric vector specifying the limits of the y-axis.
#' @param legend Show legend; Default: opposite of `label`.
#' @param legend.pos Position of legend; Default: 'NULL'.
#' @param legend.title Title of legend; Default: 'split.by'.
#' @param show_plot A logical indicating whether to display the plot.
#' @param grid A logical indicating whether to display grid lines.
#' @param w Width of the plot.
#' @param h Height of the plot.
#' @param ... Additional arguments passed to `VlnPlot`.
#'
#' @return A ggplot object representing the violin plot.
#'
#' @examples
#' # Assuming `seurat_obj` is a valid Seurat object
#' qSeuViolin(obj = seurat_obj, feature = "nFeature_RNA")
#'
#' @export
qSeuViolin <- function(
    obj,
    feature = "nFeature_RNA",
    ident = GetNamedClusteringRuns(obj = obj, v = FALSE)[1],
    assay = "RNA",
    slot = "data",
    split.by = NULL,
    colors = NULL,
    clip.outliers = TRUE,
    replace.na = FALSE,
    pt.size = 0,
    sub = NULL,
    suffix = NULL,
    suffix.2.title = FALSE,
    caption = try(.parseKeyParams(obj), silent = TRUE),
    logY = TRUE,
    hline = FALSE,
    ylab = "Expression",
    ylimit = NULL,
    legend = TRUE,
    legend.pos = NULL, # c("top", "bottom", "left", "right", "none")[2],
    legend.title = NULL,
    show_plot = TRUE,
    grid = TRUE,
    w = NULL, h = 7,
    ...) {
  message(" > Running qSeuViolin...")
  #
  stopifnot(
    "Seurat" %in% class(obj), # object must be a Seurat object
    is.logical(logY), # logY must be logical (TRUE or FALSE)
    is.logical(hline) || is.numeric(hline), # hline must be logical or numeric
    is.logical(caption) || is.character(caption), # caption must be logical or character
    is.logical(suffix.2.title), # suffix.2.title must be logical
    is.character(split.by) | is.null(split.by), # split.by must be a character or NULL
    split.by %in% colnames(obj@meta.data),
    is.character(ident),
    ident %in% colnames(obj@meta.data),
    is.character(feature),
    feature %in% colnames(obj@meta.data) || feature %in% rownames(obj) # This version of the function works on a single feature only.
  )

  is_meta_feature <- feature %in% colnames(obj@meta.data) # check if feature is in meta.data, otherwise it is a  gene name

  value_vector <-
    if(is_meta_feature){
      obj@meta.data[[feature]]
    } else {
      # Extract normalized expression values ("data" slot) from the RNA assay
      GetAssayData(obj, assay = assay, slot = slot)[feature, ]
    }


  if (exists("idents")) warning("Use arg. ident instead of idents!\n", immediate. = TRUE)
  if (exists("features")) warning("Use arg. feature instead of features!\n", immediate. = TRUE)

  split_col <- unlist(obj[[ident]])
  if (is.null(w)) w <- ceiling(length(unique(split_col)) / 6) + 6
  message("Plot width: ", w)


  ttl <- if (suffix.2.title) {
    paste(feature, "|", suffix)
  } else {
    as.character(feature)
  }
  subt <- paste(feature, "- by -", ident)
  if (!is.null(sub)) subt <- paste0(subt, "\n", sub)

  if (replace.na) {
    warning("NA's are not, but zeros are displayed on the plot. Avoid replace.na when possible", immediate. = TRUE)
    value_vector <- na.replace(x = value_vector, replace = 0)
  }


  if (clip.outliers) {
    warning("Outliers are clipped at percentiles 0.5% and 99.5%", immediate. = TRUE)
    obj@meta.data[[feature]] <- CodeAndRoll2::clip.outliers.at.percentile(
      x = value_vector, percentiles = c(.005, .995)
    )
  }

  if (!is.null(colors)) {
    stopifnot(colors %in% colnames(obj@meta.data))
    col_long <- as.factor(unlist(obj[[colors]]))
    colors <- as.factor.numeric(sapply(split(col_long, split_col), unique))
    stopifnot("colors cannot be uniquely split by ident. Set colors = NULL!" = length(colors) == length(unique(split_col)))
  }
  # browser()

  p.obj <- Seurat::VlnPlot(
    object = obj,
    features = feature, group.by = ident,
    cols = colors, split.by = split.by,
    pt.size = pt.size, ...
  ) +
    theme(axis.title.x = element_blank()) +
    labs(y = ylab) +
    ggtitle(label = ttl, subtitle = subt)

  if (!legend) p.obj <- p.obj + NoLegend()
  if (!is.null(legend.title)) p.obj <- p.obj + guides(fill = guide_legend(legend.title)) else NULL
  if (grid) p.obj <- p.obj + ggpubr::grids(axis = "y")

  # Add additional customization, if needed..
  if (!is.null(ylimit)) p.obj <- p.obj + ylim(ylimit[1], ylimit[2])
  if (logY) p.obj <- p.obj + ggplot2::scale_y_log10()
  if (hline[1]) p.obj <- p.obj + ggplot2::geom_hline(yintercept = hline)
  if (!isFALSE(caption)) p.obj <- p.obj + ggplot2::labs(caption = caption)
  if (!is.null(legend.pos)) p.obj <- p.obj + theme(legend.position = legend.pos)

  # Save the plot.
  TTL <- ppp(as.character(feature), "by", ident, suffix)
  qqSave(p.obj, title = TTL, suffix = ppp(flag.nameiftrue(logY), "violin"), w = w, h = h, limitsize = FALSE)
  if (show_plot) p.obj
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
#' @param save.plot If TRUE, the plot is saved into a file; Default: `TRUE`.
#' @param PNG If TRUE, the file is saved as a .png; Default: `TRUE`.
#' @param h Height of the plot in inches; Default: 7.
#' @param w Width of the plot in inches; Default: NULL.
#' @param nr.cols Number of columns to combine multiple feature plots, ignored if `split.by` is not NULL; Default: NULL.
#' @param assay Which assay to use ('RNA' or 'integrated'); Default: 'RNA'.
#' @param axes If TRUE, axes are shown on the plot; Default: `FALSE`.
#' @param aspect.ratio Ratio of height to width. If TRUE, the ratio is fixed at 0.6; Default: `FALSE`.
#' @param HGNC.lookup If TRUE, HGNC gene symbol lookup is performed; Default: `TRUE`.
#' @param make.uppercase If TRUE, feature names are converted to uppercase; Default: `FALSE`.
#' @param qlow Lower quantile for the color scale; Default: 'q10'.
#' @param qhigh Upper quantile for the color scale; Default: 'q90'.
#' @param check_for_2D If TRUE, checks if UMAP is 2 dimensional; Default: `TRUE`.
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
  #
  message("Feature: ", feature, " | Assay: ", assay)
  META <- obj@meta.data
  stopifnot(is(obj) == "Seurat", is.character(feature),
            "Feature not found in genes / meta.data."  = feature %in% c(Features(obj, assay = assay) , colnames(META)),
            "meta.data column is not numeric"  = if(feature %in% colnames(META)) is.numeric(META[, feature]) else TRUE,
            "UMAP is not 2 dimensional! \n Check obj@reductions[[reduction]]@cell.embeddings" =
              if (check_for_2D) ncol(obj@reductions[[reduction]]@cell.embeddings) == 2,
            reduction %in% names(obj@reductions),
            assay %in% names(obj@assays),
            "split.by column not found in meta.data / not categorical" =
              if (!is.null(splitby)) {splitby %in% colnames(META) && is.factor(META[[splitby]]) || is.character(META[[splitby]])} else TRUE
            )

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
    fname <- ww.FnP_parser(sppp(prefix, toupper(reduction), feature, assay, paste0(ncol(obj), "c"), suffix), if (PNG) "png" else "pdf")
    try(save_plot(filename = fname, plot = gg.obj, base_height = h, base_width = w)) # , ncol = 1, nrow = 1
  }
  return(gg.obj)
}



# _________________________________________________________________________________________________
#' @title clUMAP - Quick Visualization of Clustering Results with UMAP and automatically save the plot
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
#' @param caption Plot caption; optional; Default: dynamically generated from `obj`.
#' @param label.cex Size of cluster labels; Default: 7.
#' @param h Height of plot in inches; Default: 7.
#' @param w Width of plot in inches; optional; Default: NULL.
#' @param nr.cols Number of columns for facet wrap if `splitby` is not NULL; Default: NULL.
#' @param plotname Custom plot name for saving; Default: dynamically generated from `reduction` and `ident`.
#' @param cols Custom color vector for clusters; optional; Default: NULL.
#' @param palette Color palette for generating cluster colors; Default: 'glasbey'.
#' @param highlight.clusters Specific clusters to be highlighted; optional; Default: NULL.
#' @param cells.highlight Specific cells to be highlighted; optional; Default: NULL.
#' @param cols.highlight Color for highlighted cells; Default: 'red'.
#' @param sizes.highlight Size of highlighted cells; Default: 1.
#' @param label Show cluster labels; Default: `TRUE`.
#' @param repel Repel labels to avoid overlap; Default: `TRUE`.
#' @param legend Show legend; Default: opposite of `label`.
#' @param legend.pos Position of legend; Default: 'NULL'.
#' @param axes Show axes; Default: `FALSE`.
#' @param aspect.ratio Fixed aspect ratio for the plot; Default: `TRUE`.
#' @param MaxCategThrHP Maximum number of categories before simplification; Default: 200.
#' @param save.plot Save plot to file; Default: `TRUE`.
#' @param PNG Save as PNG (TRUE) or PDF (FALSE); Default: `TRUE`.
#' @param check_for_2D Ensure UMAP is 2D; Default: `TRUE`.
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
    ident = NULL,
    obj = combined.obj,
    title = ident,
    sub = NULL,
    prefix = NULL,
    suffix = make.names(sub),
    caption = .parseBasicObjStats(obj, simple = TRUE), # try(.parseKeyParams(obj = obj), silent = TRUE),
    reduction = "umap", splitby = NULL,
    label.cex = 7,
    h = 7, w = NULL,
    nr.cols = NULL,
    plotname = ppp(toupper(reduction), ident),
    cols = NULL,
    palette = c("alphabet", "alphabet2", "glasbey", "polychrome", "stepped")[4],
    max.cols.for.std.palette = 7,
    highlight.clusters = NULL, cells.highlight = NULL,
    cols.highlight = "red",
    sizes.highlight = 1,
    label = TRUE, repel = TRUE,
    legend = !label,
    legend.pos = NULL, # c("top", "bottom", "left", "right", "none")[2],
    MaxCategThrHP = 200,
    axes = NULL,
    aspect.ratio = c(FALSE, 0.6)[2],
    save.plot = MarkdownHelpers::TRUE.unless("b.save.wplots", v = FALSE),
    PNG = TRUE,
    check_for_2D = TRUE,
    ...) {
  #
  stopifnot(
    inherits(obj, "Seurat"),
    is.character(reduction) && length(reduction) == 1,
    reduction %in% names(obj@reductions),
    is.character(caption) | is.null(caption),
    is.logical(save.plot),
    is.character(suffix) | is.null(suffix),
    is.null(ident) || (is.character(ident) && length(ident) == 1 && ident %in% colnames(obj@meta.data)),
    is.null(splitby) || splitby %in% colnames(obj@meta.data)
  )
  tictoc::tic()

  if (is.null(ident)) {
    ident <- GetNamedClusteringRuns(obj, v = FALSE)[1]
    message("Identity not provided. Plotting: ", ident, "\n")
  }

  if (check_for_2D) {
    umap_dims <- ncol(obj@reductions[[reduction]]@cell.embeddings)
    if (umap_dims != 2) warning(">>> UMAP is not 2 dimensional! \n Check obj@reductions[[reduction]]@cell.embeddings")
  }

  IdentFound <- (ident %in% colnames(obj@meta.data))
  if (!IdentFound) {
    ident <- GetClusteringRuns(obj = obj, pat = "_res.*[0,1]\\.[0-9]$")[1]
    iprint("Identity not found. Plotting", ident, "\n")
  }
  identity <- obj[[ident]]
  Ident_categories <- unique(identity[, 1])
  NtCategs <- length(Ident_categories)
  if (NtCategs > 1000) warning("More than 1000 levels! qUMAP?", immediate. = TRUE)

  # Highlight specific clusters if provided _____________________________________________________
  if (!missing(highlight.clusters)) {
    if (!(all(highlight.clusters %in% identity[, 1]))) {
      MSG <- paste(
        "Some clusters not found in the object! Missing:",
        kppc(setdiff(highlight.clusters, Ident_categories)), "\nFrom:\n",
        kppc(sort(Ident_categories))
      )
      warning(MSG, immediate. = TRUE)
    }

    idx.ok <- identity[, 1] %in% highlight.clusters
    stopifnot("minimum 10 cells are needed" = sum(idx.ok) > 10)

    highlight.these <- rownames(identity)[idx.ok]
    PCT <- percentage_formatter(length(highlight.these) / ncol(obj), suffix = "or")

    # Annotation to subtitle _________________________________________________________________
    sub2 <- paste(PCT, length(highlight.these), "cells in", ident, "are highlighted")
    sub3 <- paste("Highlighted clusters:", kppc(highlight.clusters))
    sub <- if (is.null(sub)) ppnl(sub2, sub3) else ppnl(sub, sub2, sub3)

    # title <- kpipe(ident, )
  } else {
    highlight.these <- NULL
  }

  # Message if highlighting cells _____________________________________________________________
  if (!missing(cells.highlight)) {
    highlight.these <- cells.highlight
    message("Highlighting ", length(highlight.these), " cells, e.g.: ", head(highlight.these))
    message("cols.highlight: ", cols.highlight, " | sizes.highlight: ", sizes.highlight)
  }

  if (is.null(cols)) {
    cols <- if (NtCategs > max.cols.for.std.palette) {
      getDiscretePaletteObj(
        ident.used = ident, palette.used = palette,
        obj = obj, show.colors = FALSE
      )
    }
  }


  # if (FALSE) cols <- adjustcolor(cols, alpha.f = alpha)

  if (!is.null(highlight.these)) {
    cols <- "lightgrey"
  }

  # Plot _________________________________________________________________________________________
  if (NtCategs > MaxCategThrHP) {
    iprint("Too many categories (", NtCategs, ") in ", ident, "- use qUMAP for continuous variables.")
  } else {
    if (length(unique(identity)) < MaxCategThrHP) {
      gg.obj <-
        Seurat::DimPlot(
          object = obj, group.by = ident,
          cols = cols,
          reduction = reduction, split.by = splitby,
          ncol = nr.cols,
          cells.highlight = highlight.these,
          cols.highlight = cols.highlight,
          sizes.highlight = sizes.highlight,
          label = label, repel = repel, label.size = label.cex,
          ...
        ) +
        ggtitle(label = title, subtitle = sub) +
        if (!legend) NoLegend() else NULL
    }

    # Additional options ___________________________________________________
    if (is.null(axes)) gg.obj <- gg.obj + Seurat::NoAxes()
    if (!is.null(caption)) gg.obj <- gg.obj + ggplot2::labs(caption = caption)
    if (!is.null(legend.pos)) gg.obj <- gg.obj + ggplot2::theme(legend.position = legend.pos)
    if (aspect.ratio) gg.obj <- gg.obj + ggplot2::coord_fixed(ratio = aspect.ratio)
    if (legend) suffix <- paste0(suffix, ".lgnd")

    # Save plot ___________________________________________________________
    if (save.plot) {
      pname <- sppp(prefix, plotname, paste0(ncol(obj), "c"), suffix, sppp(highlight.clusters))
      fname <- ww.FnP_parser(pname, if (PNG) "png" else "pdf")
      try(save_plot(filename = fname, plot = gg.obj, base_height = h, base_width = w))
    }
    tictoc::toc()
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
#' Default: `GetClusteringRuns()[1]`.
#' @param h Height of the plot; Default: `7`.
#' @param w Width of the plot; Default: `5`.
#' @param show_plot Logical; if `TRUE`, the plot will be displayed in the RStudio viewer;
#' Default: `TRUE`.
#' @param ... Additional arguments to be passed to the `DimPlot` function.#'
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
                           show_plot = TRUE,
                           ...) {
  stopifnot(is(obj, "Seurat"),
    "Ident no found the object!" = ident %in% colnames(obj@meta.data),
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
  if (show_plot) print(pl)

  ggplot2::ggsave(
    filename = extPNG(kollapse("cells", COI, collapseby = ".")),
    height = h, width = w
  )
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
#'   DimPlot.ClusterNames(obj = combined.obj)
#' }
#' }
#' @export
DimPlot.ClusterNames <- function(
    obj = combined.obj,
    ident = GetNamedClusteringRuns(obj = obj, v = FALSE)[1],
    reduction = "umap",
    title = ident,
    ...) {
  #
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
#' @param subdir Should plots be saved in a sub-directory? Default: `TRUE`.
#' @param foldername Folder name to save the generated plots. Default: The name of the list of genes.
#' @param subtitle.from.names Should the subtitle be extracted from the names of the gene symbols,
#' eg: `c("Astrocytes" = "AQP4")` ? Default: `TRUE`.
#' @param plot.reduction Dimension reduction technique to use for plots. Default: 'umap'
#' @param intersectionAssay The assay to intersect with, either 'RNA' or 'integrated'. Default: 'RNA'
#' @param layout Layout orientation of the plot. Default: 'wide'
#' @param colors Vector of colors to be used in the plot. Default: c("grey", "red")
#' @param nr.Col Number of columns in the plot grid. Default: 2
#' @param nr.Row Number of rows in the plot grid. Default: 4
#' @param cex Point size in the plot. Default: round(0.1/(nr.Col * nr.Row), digits = 2)
#' @param gene.min.exp Minimum gene expression level for plotting. Default: 'q01'
#' @param gene.max.exp Maximum gene expression level for plotting. Default: 'q99'
#' @param prefix Prefix for the plot filenames. Default: NULL
#' @param suffix Suffix for the plot filenames. Default: NULL
#' @param background_col Background color of the plots. Default: "white"
#' @param saveGeneList Should the list of genes be saved? Default: `FALSE`.
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
    subdir = TRUE,
    foldername = substitute(list.of.genes),
    subtitle.from.names = TRUE,
    plot.reduction = "umap",
    intersectionAssay = c("RNA", "integrated")[1],
    layout = c("tall", "wide", FALSE)[2],
    colors = c("grey", "red"),
    nr.Col = 2, nr.Row = 4,
    raster = if (ncol(obj) > 1e5) TRUE else FALSE,
    cex = round(0.1 / (nr.Col * nr.Row), digits = 2),
    cex.min = if (raster) TRUE else FALSE,
    gene.min.exp = "q01", gene.max.exp = "q99",
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
  if (subdir) create_set_SubDir(final.foldername, "/", verbose = FALSE)

  if (is.null(names(list.of.genes))) subtitle.from.names <- FALSE

  list.of.genes.found <- check.genes(
    genes = list.of.genes, obj = obj,
    assay.slot = intersectionAssay, makeuppercase = FALSE
  )
  if (subtitle.from.names) names(list.of.genes.found) <- as.character(flip_value2name(list.of.genes)[list.of.genes.found])
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

    # Remove the legend and axes
    for (j in 1:length(plot.list)) {
      plot.list[[j]] <- plot.list[[j]] + NoLegend() + NoAxes()
      if (aspect.ratio) plot.list[[j]] <- plot.list[[j]] + ggplot2::coord_fixed(ratio = aspect.ratio)
      if (subtitle.from.names) plot.list[[j]] <- plot.list[[j]] + ggplot2::ggtitle(label = genes[j], subtitle = names(genes)[j])
    }
    # browser()

    pltGrid <- cowplot::plot_grid(plotlist = plot.list, ncol = nr.Col, nrow = nr.Row)
    # cowplot::ggsave2(filename = plotname, width = w, height = h, bg = background_col, plot = pltGrid)
    cowplot::save_plot(
      plot = pltGrid, filename = plotname,
      base_width = w, base_height = h,
      bg = background_col
    )
  }

  if (subdir) MarkdownReports::create_set_OutDir(ParentDir, verbose = FALSE)
  if (saveGeneList) {
    if (is.null(obj@misc$gene.lists)) obj@misc$gene.lists <- list()
    obj@misc$gene.lists[[substitute_deparse(list.of.genes)]] <- list.of.genes.found
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
#' @param sizes.highlight Size of highlighted cells; Default: 1.
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
    sizes.highlight = 1,
    subdir = TRUE,
    prefix = NULL, suffix = NULL,
    background_col = "white",
    aspect.ratio = c(FALSE, 0.6)[2],
    saveGeneList = FALSE,
    w = 8.27, h = 11.69, scaling = 1,
    format = c("jpg", "pdf", "png")[1],
    ...) {
  message(" > Running multiSingleClusterHighlightPlots.A4...")

  NrCellsPerCluster <- sort(table(obj[[ident]]), decreasing = TRUE)
  stopifnot(
    ident %in% colnames(obj@meta.data),
    "Some clusters too small (<20 cells). See: table(obj[[ident]]) | Try: removeResidualSmallClusters()" =
      all(NrCellsPerCluster > 20)
  )

  tictoc::tic()
  ParentDir <- OutDir
  if (is.null(foldername)) foldername <- "clusters"
  if (subdir) MarkdownReports::create_set_SubDir(paste0(foldername, "-", plot.reduction), "/")

  clusters <- sort(unique(obj@meta.data[[ident]]))

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
  # browser()
  for (i in 1:length(ls.Clust)) { # for each page

    clusters_on_this_page <- as.character(clusters[ls.Clust[[i]]])
    iprint("page:", i, "| clusters", kppc(clusters_on_this_page))
    (plotname <- kpp(c(prefix, plot.reduction, i, "clusters", ls.Clust[[i]], suffix, format)))

    plot.list <- list()
    for (j in seq(clusters_on_this_page)) {  # for each cluster
      cl <- clusters_on_this_page[j]
      message(cl)
      plot.list[[j]] <- clUMAP(
        ident = ident, obj = obj,
        highlight.clusters = cl, label = FALSE, legend = FALSE, save.plot = FALSE,
        plotname = plotname, cols = colors,
        sizes.highlight = sizes.highlight,
        h = h, w = w, ...
      )
    } # for j

    # Customize plot appearance
    for (j in 1:length(plot.list)) {
      plot.list[[j]] <- plot.list[[j]] + NoLegend() + NoAxes()
      if (aspect.ratio) {
        plot.list[[j]] <- plot.list[[j]] +
          ggplot2::coord_fixed(ratio = aspect.ratio)
      }
    } # for j2

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
    w = 2 * 11.69, h = 2 * 8.27,
    title = sppu(
      prefix,
      as.numeric(stringr::str_extract(idents, "\\d+\\.\\d+$")),
      suffix
    ),
    ...) {
  message(" > Running qClusteringUMAPS...")

  # Check that the QC markers are in the object
  idents.found <- intersect(idents, colnames(obj@meta.data))
  n.found <- length(idents.found)
  stopifnot("None of the idents found" = n.found > 0)
  message(kppws(n.found, " found of ", idents))

  if (n.found > 5) {
    idents.found <- idents.found[1:4]
    message("Only the first 4 idents will be plotted: ", idents.found)
  }

  px <- list(
    "A" = clUMAP(ident = idents[1], save.plot = FALSE, obj = obj, caption = NULL, ...) + NoAxes(),
    "B" = clUMAP(ident = idents[2], save.plot = FALSE, obj = obj, caption = NULL, ...) + NoAxes(),
    "C" = clUMAP(ident = idents[3], save.plot = FALSE, obj = obj, caption = NULL, ...) + NoAxes(),
    "D" = clUMAP(ident = idents[4], save.plot = FALSE, obj = obj, caption = NULL, ...) + NoAxes()
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
  features.found <- intersect(features, c(colnames(obj@meta.data), rownames(obj)))
  n.found <- length(features.found)
  stopifnot(
    "None of the features found" = n.found > 1,
    "Only 4 features are allowed" = n.found < 5
  )

  message(kppws(n.found, "found of", length(features), "features:", features))

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
                                # subtitles = NULL, # assume the names of the vector.
                                foldername = NULL,
                                intersectionAssay = DefaultAssay(obj),
                                plot.reduction = "umap",
                                ...) {
  message(" > Running plotQUMAPsInAFolder...")

  # Input checks
  stopifnot(
    is.character(genes),
    is.null(foldername) || is.character(foldername),
    is.character(plot.reduction)
  )

  ParentDir <- OutDir
  if (is.null(foldername)) foldername <- substitute_deparse(genes)

  MarkdownReports::create_set_SubDir(paste0(foldername, "-", plot.reduction), "/")

  list.of.genes.found <- check.genes(
    genes = genes, obj = obj,
    assay.slot = intersectionAssay, makeuppercase = FALSE
  )

  for (i in seq_along(list.of.genes.found)) {
    g <- list.of.genes.found[i];   message(g)
    qUMAP(feature = g, sub = names(list.of.genes.found)[i],
          reduction = plot.reduction, obj = obj, ...)
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
#' @param filt_coding_genes Logical indicating whether to filter out non-coding genes; Default: `TRUE`.
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
    filt_coding_genes = TRUE,
    ...) {
  message(" > Running PlotTopGenesPerCluster...")

  if(filt_coding_genes) {
    genes <- df_markers$gene
    df_markers <- df_markers[genes %in% filterCodingGenes(genes), ]
  }

  topX.markers <- GetTopMarkers(
    dfDE = df_markers, n = nrGenes,
    order.by = order.by
  )

  ls.topMarkers <- splitbyitsnames(topX.markers)
  for (i in 1:length(ls.topMarkers)) {
    multiFeaturePlot.A4(
      list.of.genes = ls.topMarkers[[i]], obj = obj, subdir = TRUE, foldername = ppp("TopGenes.umaps"),
      prefix = ppp("DEG.markers.res", cl_res, "cluster", names(ls.topMarkers)[i]),
      ...
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
    QC.Features = c("nFeature_RNA", "percent.ribo", "percent.mito", "nuclear.fraction", "percent.HGA"),
    prefix = "QC.markers.4.UMAP",
    suffix = "",
    title = sppu(prefix, QC.Features, suffix),
    nrow = 2, ncol = 2,
    ...) {
  message(" > Plotting qQC.plots.BrainOrg...")

  # Check that the QC markers are in the object
  QC.Features.Found <- intersect(QC.Features, union(colnames(obj@meta.data), rownames(obj)))
  QC.Features.Found.in.meta <- intersect(QC.Features, colnames(obj@meta.data))
  n.found <- length(QC.Features.Found)
  message(kppws(n.found, " found: ", QC.Features.Found))
  stopifnot(n.found > 1)


  # Count the number of NAs in specified columns
  na_counts <- sapply(X = obj@meta.data[, QC.Features.Found.in.meta], function(x) sum(is.na(x)))

  # Raise a warning if there are any NAs
  if(length(na_counts)){
    if (sum(na_counts) > 0) {
      warning(sprintf("There are %d NA values found\n", na_counts),
              immediate. = TRUE
      )
    }
  }

  px <- list(
    "A" = qUMAP(QC.Features.Found[1], save.plot = FALSE, obj = obj, ...) + NoAxes(),
    "B" = qUMAP(QC.Features.Found[2], save.plot = FALSE, obj = obj, ...) + NoAxes(),
    "C" = qUMAP(QC.Features.Found[3], save.plot = FALSE, obj = obj, ...) + NoAxes(),
    "D" = qUMAP(QC.Features.Found[4], save.plot = FALSE, obj = obj, ...) + NoAxes()
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
#' If FALSE, a predefined list of key brain organoid markers is used; Default: `FALSE`.
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
  message(" > Running qMarkerCheck.BrainOrg...")

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

PlotTopGenes <- function(obj = combined.obj, n = 8
                         , exp.slot = "expr.q99") {
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
#' @param FlipReductionBackupToo Boolean indicating whether to also flip coordinates in the backup slot; Default: `TRUE`.
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
#' @param swap If TRUE, reverses the ordering direction; Default: `FALSE`.
#' @param reduction Dimension reduction technique used for cluster positioning ('umap', 'tsne', or 'pca');
#' Default: 'umap'.
#' @param ident Clustering resolution identifier used to fetch cluster labels from `obj` metadata;
#' Default: 'integrated_snn_res.0.5'.
#' @param plot If TRUE, plots the UMAP with new cluster names; Default: `TRUE`.
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
#' @importFrom CodeAndRoll2 translate
#' @importFrom Stringendo kpp kppu iprint
#' @importFrom Seurat FetchData
AutoNumber.by.UMAP <- function(obj = combined.obj,
                               reduction = "umap",
                               dim = 1, swap = FALSE,
                               ident = GetClusteringRuns(obj = obj)[1],
                               plot = TRUE,
                               obj.version = obj@version) {
  dim_name <- kppu(reduction, dim)
  if (obj.version < "5") dim_name <- toupper(dim_name)
  message("Obj. version: ", obj.version, " \ndimension name: ", dim_name)
  message("Resolution: ", ident)

  stopifnot("Identity not found." = ident %in% colnames(obj@meta.data))

  coord.umap <- obj@reductions$umap@cell.embeddings[, dim_name]

  # Get cluster labels
  identX <- as.character(obj@meta.data[[ident]])

  # Order clusters by median coordinate
  ls.perCl <- split(coord.umap, f = identX)
  MedianClusterCoordinate <- sapply(ls.perCl, median)

  OldLabel <- names(sort(MedianClusterCoordinate, decreasing = swap))
  NewLabel <- as.character(0:(length(MedianClusterCoordinate) - 1))
  NewMeta <- CodeAndRoll2::translate(vec = identX, old = OldLabel, new = NewLabel)
  NewMetaCol <- kpp(ident, "ordered")
  iprint("NewMetaCol:", NewMetaCol)

  obj[[NewMetaCol]] <- NewMeta
  if (plot) {
    x <- clUMAP(obj, ident = NewMetaCol)
    print(x)
  }
  return(obj)
}


# _________________________________________________________________________________________________
# General ______________________________ ----
# _________________________________________________________________________________________________



# _________________________________________________________________________________________________
# DGEA ______________________________ ----
# _________________________________________________________________________________________________


#' @title scEnhancedVolcano
#'
#' @description This function creates an enhanced volcano plot.
#'
#' @param toptable A data frame with the results of differential gene expression analysis.
#' @param x The x-axis, which is typically the average log2 fold change.
#' @param y The y-axis, which is typically the adjusted p-value.
#' @param title The title of the plot.
#' @param lab A vector of gene symbols to label on the plot.
#' @param selectLab A vector of gene symbols to select for labeling.
#' @param min.p The minimum p-value, to trim high values on the Y-axis.
#' @param max.l2fc The maximum log2 fold change, to trim high values on the X-axis.
#' @param min.pct.cells The minimum percentage of cells in which a gene must be expressed to be included in the plot.
#' @param pCutoffCol The column in the toptable that contains the p-value cutoff.
#' @param pCutoff The p-value cutoff.
#' @param FCcutoff The fold change cutoff.
#'
#' @param suffix A string to append to the filename/title of the plot.
#' @param suffix A string to append to the filename/title of the plot.
#' @param caption The first line of caption of the plot.
#' @param caption2 The second line of caption of the plot.
#' @param count_stats Logical. Calculates a data frame with the count statistics.
#' @param drawConnectors Whether to draw connectors between the labels and the points.
#' @param max.overlaps The maximum number of labels that can overlap.
#' @param also.pdf Logical. Whether to save the plot as a PDF in addition to the default png format.
#' @param h The height of the plot.
#' @param w The width of the plot.
#' @param ... Pass any other parameter to `EnhancedVolcano::EnhancedVolcano()`.
#'
#' @return A ggplot object.
#'
#' @importFrom EnhancedVolcano EnhancedVolcano
#'
#' @export scEnhancedVolcano

scEnhancedVolcano <- function(
    toptable,
    x = "avg_log2FC",
    y = "p_val_adj",
    title = paste("DGEA"),
    lab = rownames(toptable),
    selectLab = trail(filterCodingGenes(lab), 10),
    min.p = 1e-50,
    max.l2fc = Inf,
    min.pct.cells = 0.1,
    pCutoffCol = "p_val_adj",
    pCutoff = 1e-3,
    FCcutoff = 1,

    suffix = NULL,
    caption = paste("Min. Fold Change in Input:", .estMinimumFC(toptable)),
    caption2 = paste("min p_adj:", min.p, "(Y-axis values clipped at)"),
    count_stats = TRUE,
    drawConnectors = TRUE,
    max.overlaps = Inf,
    also.pdf = FALSE,
    h = 9, w = h,
    ...) {
  #
  message(
    "\nMin. log2fc: ", FCcutoff, "\nMax. p-adj: ", pCutoff,
    "\nMin. p-adj (trim high y-axis): ", min.p,
    "\nMin. pct cells expressing: ", min.pct.cells
  )
  stopifnot(nrow(toptable) > 5)


  # Filter min. cells expressing.
  toptable <- toptable |> dplyr::filter(pct.1 > min.pct.cells | pct.2 > min.pct.cells)

  # calculate true min pct cells expressing (maybe input prefiltered above thr. already).
  min.pct.cells <- toptable |>
    select(pct.1, pct.2) |>
    as.matrix() |>
    rowMax() |>
    min()

  # Clip p-values.
  toptable[["p_val_adj"]] <-
    clip.at.fixed.value(x = toptable[["p_val_adj"]], thr = min.p, above = FALSE)

  # Clip log2FC.
  if (max.l2fc < Inf) {
    toptable[["avg_log2FC"]] <-
      clip.at.fixed.value(x = toptable[["avg_log2FC"]], thr = -max.l2fc, above = FALSE)
    toptable[["avg_log2FC"]] <-
      clip.at.fixed.value(x = toptable[["avg_log2FC"]], thr = max.l2fc, above = TRUE)
  }

  # Add statistical information to the subtitle.
  if (count_stats) {
    enr_stats <- unlist(countRelevantEnrichments(
      df = toptable, logfc_col = x, pval_col = y,
      logfc_cutoff = FCcutoff, pval_cutoff = pCutoff
    ))
    stat_info <- kppws("Genes", intermingle2vec(names(enr_stats), enr_stats), "(red)")
    subtitle <- paste0(
      stat_info, "\n",
      paste(
        "Cutoffs: max.p_adj: ", pCutoff, " |  min.log2FC: ", FCcutoff,
        " |  min.pct.cells: ", min.pct.cells
      )
    )
  }
  caption <- paste0(caption, "\n", caption2)

  # Create an enhanced volcano plot.
  # try.dev.off();
  pobj <- EnhancedVolcano::EnhancedVolcano(
    toptable = toptable,
    x = x, y = y,
    title = title, subtitle = subtitle,
    lab = lab, selectLab = selectLab,
    caption = caption,
    pCutoffCol = pCutoffCol,
    pCutoff = pCutoff,
    FCcutoff = FCcutoff,
    drawConnectors = drawConnectors,
    max.overlaps = max.overlaps,
    ...
  )

  print(pobj)
  # Save the plot.
  qqSave(
    ggobj = pobj, title = paste0("Volcano.", make.names(title), suffix), also.pdf = also.pdf,
    h = h, w = w
  )
  return(pobj)
}

# ________________________________________________________________________
#' @title Estimate Minimum Log2-Based Fold Change
#'
#' @description This function estimates the minimum log2-based fold change from a data frame column.
#'
#' @param df A data frame containing the fold change data. Default: `df.m.UL`.
#' @param col A character string specifying the column name containing log2 fold change values. Default: "avg_log2FC".
#'
#' @return The minimum log2-based fold change, rounded and transformed from log2 to linear scale.
#'
#' @examples
#' \dontrun{
#' df <- data.frame(avg_log2FC = c(-1, -0.5, 0.5, 1))
#' .estMinimumFC(df, "avg_log2FC")
#' # .estMinimumFC(df = df.m.UL, col = "avg_log2FC")
#' }
#' @return The minimum log2-based fold change, rounded and transformed from log2 to linear scale.

.estMinimumFC <- function(df, col = "avg_log2FC") {
  lfc <- df[[col]]
  lfc_enr <- min(lfc[lfc > 0])
  lfc_depl <- abs(max(lfc[lfc < 0]))
  estim_min_l2fc <- min(lfc_enr, lfc_depl)
  return(iround(2^estim_min_l2fc))
}


# ________________________________________________________________________
#' @title Count Relevant Enrichments
#'
#' @description This function counts the number of relevantly expressed genes from a differential
#' gene expression table. It considers genes to be relevant if they fall under a maximum p-value
#' cutoff and are above a minimum log2 fold change cutoff. The function reports the number of
#' enriched and depleted genes.
#'
#' @param df Data frame containing the gene expression data.
#' @param pval_col Character. Name of the column containing the adjusted p-values. Default: "p_val_adj".
#' @param logfc_col Character. Name of the column containing the log2 fold change values. Default: "avg_log2FC".
#' @param pval_cutoff Numeric. The maximum adjusted p-value to consider a gene relevant. Default: 1e-2.
#' @param logfc_cutoff Numeric. The minimum log2 fold change to consider a gene relevant. Default: 1.
#'
#' @return A list with the counts of enriched and depleted genes.
#' @export
#'
#' @examples
#' df <- data.frame(
#'   p_val_adj = c(0.001, 0.02, 0.03, 0.0001),
#'   avg_log2FC = c(1.5, -2, 0.5, 2)
#' )
#' countRelevantEnrichments(df)
countRelevantEnrichments <- function(df,
                                     pval_col = "p_val_adj", logfc_col = "avg_log2FC",
                                     pval_cutoff = 1e-2, logfc_cutoff = 1) {
  #
  stopifnot(
    is.data.frame(df),
    pval_col %in% colnames(df),
    logfc_col %in% colnames(df),
    is.numeric(pval_cutoff),
    is.numeric(logfc_cutoff)
  )

  relevant_genes <- df |>
    dplyr::filter(!!sym(pval_col) <= pval_cutoff)

  enriched_count <- relevant_genes |>
    dplyr::filter(!!sym(logfc_col) >= logfc_cutoff) |>
    nrow()

  depleted_count <- relevant_genes |>
    dplyr::filter(!!sym(logfc_col) <= -logfc_cutoff) |>
    nrow()

  return(list(enriched = enriched_count, depleted = depleted_count))
}



# _________________________________________________________________________________________________
# GO-term enrichment ______________________________ ----
# _________________________________________________________________________________________________



# ________________________________________________________________________
#' @title Perform GO Enrichment Analysis
#'
#' @description This function performs Gene Ontology (GO) enrichment analysis using the
#' `clusterProfiler::enrichGO` function. It takes the gene list, universe, organism database,
#' gene identifier type, and ontology type as inputs and returns the enrichment results.
#'
#' @param genes Character vector. List of genes for enrichment analysis. Default: NULL.
#' @param universe Character vector. Background gene list (universe). Default: NULL.
#' @param org_db Character. Organism-specific database to use (e.g., 'org.Hs.eg.db'). Default: 'org.Hs.eg.db'.
#' @param key_type Character. Gene identifier type (e.g., 'SYMBOL', 'ENTREZID'). Default: 'SYMBOL'.
#' @param ont Character. Ontology type to use (e.g., 'BP', 'MF', 'CC'). Default: 'BP'.
#' @param pAdjustMethod Character. Method for p-value adjustment. Default: 'BH' for Benjamini-Hochberg.
#' @param pvalueCutoff Numeric. P-value cutoff for significance. Default: 0.05.
#' @param qvalueCutoff Numeric. Q-value cutoff for significance. Default: 0.2.
#' @param save Logical. Save the results as a data frame. Default: `TRUE`.
#' @param suffix Character. Suffix to append to the output file name. Default: 'GO.Enrichments'.
#' @param check.gene.symbols Logical. Check gene symbols for validity. Default: `TRUE`.
#' @param ... Additional arguments to pass to `clusterProfiler::enrichGO`.

#'
#' @return A data frame with GO enrichment results.
#' @export
#'
#' @examples
#' \dontrun{
#' gene_list <- rownames(df.m.DL.up.in.TSC)
#' background_genes <- names(all.genes)
#' go_results <- performGOEnrichment(gene_list, background_genes, "org.Hs.eg.db", "SYMBOL", "BP")
#' print(go_results)
#' }
scGOEnrichment <- function(genes, universe = NULL,
                           org_db = "org.Hs.eg.db", key_type = "SYMBOL", ont = "BP",
                           pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                           save = TRUE,
                           suffix = NULL,
                           check.gene.symbols = TRUE,
                           ...) {
  # Load required library
  stopifnot("Package 'clusterProfiler' must be installed to use this function." = require("clusterProfiler"))

  # Input assertions
  stopifnot(
    is.character(genes) | is.null(genes),
    (is.character(universe) | is.null(universe)),
    (length(universe) > 100 | is.null(universe)),
    is.character(org_db), is.character(key_type), is.character(ont),
    is.character(ont)
  )

  if (is.null(genes) | length(genes) == 0) {
    return(NULL)
  }

  # check.gene.symbols
  if (check.gene.symbols) {
    x <- checkGeneSymbols(genes, species = "human")
    genes <- x[x[, "Approved"], 1]
  }

  message("Performing enrichGO() analysis...")
  message(length(genes), " approved genes of interest, in ", length(universe), " background genes.")

  go_results <- clusterProfiler::enrichGO(
    gene = genes,
    universe = universe,
    pAdjustMethod = pAdjustMethod,
    OrgDb = org_db,
    keyType = key_type,
    pvalueCutoff = pvalueCutoff,
    qvalueCutoff = qvalueCutoff,
    ont = ont,
    ...
  )

  nr_of_enr_terms <- length(go_results@result$"ID")

  # Output assertions
  if (nrow(go_results) < 1) warning("No enriched terms found!", immediate. = TRUE)

  if (save) {
    xsave(go_results,
      suffix = kpp("enr", nr_of_enr_terms, suffix),
      showMemObject = FALSE, saveParams = FALSE, allGenes = FALSE
    )
  }
  message("\nNr of enriched terms: ", nr_of_enr_terms)

  return(go_results)
}


# ________________________________________________________________________
#' @title Filter GO Enrichment Results
#'
#' @description This function filters GO enrichment results based on adjusted p-value and q-value
#' cutoffs, and retrieves the descriptions of the filtered results.
#'
#' @param df.enrichments An object of class `enrichResult` containing the GO enrichment results.
#' @param colname Character. The name of the column containing the GO-term names, or else.
#' @param pvalueCutoff Numeric. The p-value cutoff for filtering the results. Default: NULL, meaning
#' that the default cutoff of the input object is used. It is stored in `df.enrichments@pvalueCutoff`.
#' @param qvalueCutoff Numeric. The q-value cutoff for filtering the results. Default: NULL,
#' meaning that the default cutoff of the input object is used. It is stored in `df.enrichments@qvalueCutoff`.
#'
#' @return A character vector of descriptions of the filtered GO enrichment results.
#'
#' @examples
#' # Assuming GO.Enriched.DL.Ctrl is an object of class `enrichResult` created by clusterprofiler or equivalent
#' descriptions <- filterGoEnrichment(GO.Enriched.DL.Ctrl)
#' print(descriptions)
#'
#' @importFrom dplyr filter pull
#'
#' @export
filterGoEnrichment <- function(df.enrichments,
                               pvalueCutoff = NULL,
                               qvalueCutoff = NULL,
                               colname = "Description") {
  # Input assertions
  stopifnot(
    "enrichResult" %in% class(df.enrichments),
    !is.null(df.enrichments@result),
    !is.null(df.enrichments@pvalueCutoff),
    !is.null(df.enrichments@qvalueCutoff)
  )

  pvalueCutoff <- if (is.null(pvalueCutoff)) df.enrichments@pvalueCutoff else pvalueCutoff
  qvalueCutoff <- if (is.null(qvalueCutoff)) df.enrichments@qvalueCutoff else qvalueCutoff

  message(paste(
    "Filtering GO enrichment results with \np-value cutoff",
    pvalueCutoff, "and q-value cutoff", qvalueCutoff
  ))

  # Filter and retrieve GO
  descriptions <- df.enrichments@result |>
    dplyr::filter(p.adjust < pvalueCutoff) |>
    dplyr::filter(qvalue < qvalueCutoff) |>
    dplyr::pull(!!sym(colname))

  # Output assertions
  stopifnot(is.character(descriptions))
  message("\nNr of enriched terms: ", length(descriptions))

  return(descriptions)
}

# Example usage
# Assuming GO.Enriched.DL.Ctrl is an object of class `enrichResult`
# descriptions <- filterGoEnrichment(GO.Enriched.DL.Ctrl)
# print(descriptions)



# ________________________________________________________________________
#' @title Barplot GO Enrichment Results by enrichplot
#'
#' @description This function creates a bar plot of GO enrichment analysis results using the
#' `enrichplot::barplot.enrichResult` function. It also allows saving the plot to a file.
#'
#' @param df.enrichment Data frame. Enrichment results from GO analysis. Default: NULL.
#' @param showCategory Integer. Number of categories (GO-terms as bars) to show in the plot. Default: 20.
#' @param tag Character. Tag to be added to the title of the plot. Default: "in ...".
#' @param universe Character. Background gene list (universe). Default: `df.enrichment@universe`.
#' @param title Character. Title of the plot. Default: "GO Enrichment Analysis" followed by `tag`.
#' @param subtitle Character. Subtitle of the plot. Default: NULL.
#' @param caption Character. Caption of the plot. Default: constructed from input parameters.
#' @param save Logical. Whether to save the plot to a file. Default: `TRUE`.
#' @param also.pdf Save plot in both png and pdf formats.
#' @param save.obj Logical. Whether to save the ggplot object. Default: FALSE.
#' @param h Height of the plot canvas, calculated as the height of an A4 page times `scale`; Default: `8.27 * scale`.
#' @param w Width of the plot canvas, calculated as the width of an A4 page times `scale`; Default: `11.69 * scale`.
#' @param ... Additional arguments passed to `enrichplot::barplot.enrichResult`.
#'
#' @importFrom ggplot2 labs
#'
#' @return None. The function prints the plot and optionally saves it.
#' @export
#'
#' @examples
#' \dontrun{
#' df.enrichment <- data.frame() # Example enrichment results data frame
#' plotGOEnrichment(df.enrichment)
#' }
scBarplotEnrichr <- function(df.enrichment,
                             showCategory = 20,
                             label_format = 30,
                             tag = "...",
                             universe = df.enrichment@universe,
                             title = paste("GO Enriched Terms", tag),
                             subtitle = kppws("Input: ", substitute_deparse(df.enrichment)),
                             caption = paste0(
                               "Input genes: ", length(df.enrichment@"gene"),
                               " | Enriched terms: ", nrow(df.enrichment),
                               " | Shown: ", min(showCategory, nrow(df.enrichment)),
                               " | background genes: ", length(universe)
                             ),
                             save = TRUE,
                             also.pdf = FALSE,
                             save.obj = FALSE,
                             w = 10, h = 10,
                             ...) {

  stopifnot("Package 'enrichplot' must be installed to use this function." = require("enrichplot"))

  if (tag == "...") warning("Please provide a tag describing where are the enrichments.", immediate. = TRUE)
  nr_GOENR_input_genes <- length(df.enrichment@"gene")

  pobj <-
    if (is.null(df.enrichment) || nr_terms < 1) {
      Seurat.utils:::.emptyAnnotatedPlot(label = "No enriched terms input!")

    } else if (nr_GOENR_input_genes < 5) {
      Seurat.utils:::.emptyAnnotatedPlot(label = "Too few input genes for GO enrichment (<5).")

    } else {
  # pobj <-
  #   if (nrow(df.enrichment) < 1 | is.null(df.enrichment)) {
  #     warning("No enriched terms input!", immediate. = TRUE)
  #     ggplot() +
  #       theme_void() +
  #       annotate("text",
  #         x = 1, y = 1, label = "NO ENRICHMENT",
  #         size = 8, color = "red", hjust = 0.5, vjust = 0.5
  #       )
  #   } else if (nr_GOENR_input_genes < 5) {
  #     warning("Very few inputs for GOENR", immediate. = TRUE)
  #     ggplot() +
  #       theme_void() +
  #       annotate("text",
  #         x = 1, y = 1, label = "TOO FEW GENES (<5)",
  #         size = 8, color = "red", hjust = 0.5, vjust = 0.5
  #       )
  #   } else {
      enrichplot:::barplot.enrichResult(df.enrichment, showCategory = showCategory, label_format = label_format)
    }

  pobj <- pobj + ggplot2::labs(title = title, subtitle = subtitle, caption = caption)

  if (save) {
    qqSave(pobj, title = title, w = w, h = h, also.pdf = also.pdf, save.obj = save.obj)
  }

  return(pobj)
}




# ________________________________________________________________________
#' @title Enrichment Map (GO term network) by enrichplot
#'
#' @description
#' Wrapper around `enrichplot::emapplot()` to visualize GO enrichment results
#' as a network. Nodes are enriched GO terms, edges represent gene overlap.
#' Includes safety checks, informative fallback plots, and optional saving,
#' mirroring the behavior of `scBarplotEnrichr()`.
#'
#' @param df.enrichment enrichResult object (e.g. from clusterProfiler::enrichGO).
#' @param showCategory Integer. Number of GO terms (nodes) to show. Default: 15.
#' @param min_edge Numeric. Minimum similarity (overlap) to draw edges. Default: 0.2.
#' @param tag Character. Tag added to the plot title. Default: "...".
#' @param universe Character. Background gene list. Default: `df.enrichment@universe`.
#' @param title Character. Plot title. Default: "GO Enrichment Map" + tag.
#' @param subtitle Character. Subtitle. Default: derived from input object.
#' @param caption Character. Caption. Default: constructed from input parameters.
#' @param layout Character. igraph layout name passed to emapplot. Default: "kk".
#' @param save Logical. Whether to save the plot. Default: TRUE.
#' @param also.pdf Logical. Save both png and pdf. Default: FALSE.
#' @param save.obj Logical. Whether to save the ggplot object. Default: FALSE.
#' @param w Width in inches. Default: 10.
#' @param h Height in inches. Default: 10.
#' @param ... Additional arguments passed to `enrichplot::emapplot()`.
#'
#' @importFrom ggplot2 labs theme_void annotate
#'
#' @return A ggplot object (invisibly if saved).
#' @export
#'
#' @examples
#' \dontrun{
#' scEmapplotEnrichr(df.enrichment, tag = "Cluster 3 neurons")
#' }
scEmapplotEnrichr <- function(
    df.enrichment,
    showCategory = 15,
    min_edge = 0.2,
    tag = "...",
    universe = df.enrichment@universe,
    title = paste("GO Enrichment Map", tag),
    subtitle = kppws("Input: ", substitute_deparse(df.enrichment)),
    caption = paste0( "Input genes: ", length(df.enrichment@"gene"),
                      " | Enriched terms: ", nrow(df.enrichment),
                      " | Shown: ", min(showCategory, nrow(df.enrichment)),
                      " | background genes: ", length(universe),
                      " | min edge overlap: ", min_edge),
    label_format = NULL,
    layout = "kk",
    cex_label_category = 0.8,
    save = TRUE,
    also.pdf = FALSE,
    save.obj = FALSE,
    w = 10, h = 10,
    ...
) {

  stopifnot("Package 'enrichplot' must be installed." = requireNamespace("enrichplot", quietly = TRUE) )

  if (tag == "...") {
    warning(
      "Please provide a tag describing where the enrichments come from.",
      immediate. = TRUE
    )
  }

  nr_GOENR_input_genes <- length(df.enrichment@"gene")

  pobj <-
    if (is.null(df.enrichment) || nr_terms < 1) {
      Seurat.utils:::.emptyAnnotatedPlot(label = "No enriched terms input!")

    } else if (nr_GOENR_input_genes < 5) {
      Seurat.utils:::.emptyAnnotatedPlot(label = "Too few input genes for GO enrichment (<5).")

    } else {

      # similarity matrix is computed internally by emapplot()
      enrichplot::emapplot(
        x = df.enrichment,
        showCategory = showCategory,
        layout.params  = list(layout = layout),
        edge.params    = list(min = min_edge),
        cex.params     = list(category_label = cex_label_category),
        cluster.params = list(label_format = label_format),
        ...
      )
    }

  pobj <- pobj +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      caption = caption
    )

  if (save) {
    qqSave(pobj, title = title, w = w, h = h, also.pdf = also.pdf, save.obj = save.obj)
  }

  return(pobj)
}



# ________________________________________________________________________
#' @title Gene–Concept Network Plot (cnetplot wrapper)
#'
#' @description
#' Wrapper around `enrichplot::cnetplot()` to visualize the gene–concept
#' (e.g. GO / KEGG) network for enrichment results. The plot shows which genes
#' drive which enriched terms, optionally colored by fold change (e.g. DE).
#' Behavior mirrors `scBarplotEnrichr()` and `scEmapplotEnrichr()` with
#' safety checks, informative fallbacks, and optional saving.
#'
#' @param df.enrichment enrichResult or gseaResult object.
#' @param foldChange Named numeric vector of gene-level statistics
#'   (e.g. logFC), names must match gene IDs in enrichment.
#' @param showCategory Integer. Number of enriched terms to show. Default: 10.
#' @param tag Character. Tag added to the plot title. Default: "...".
#' @param title Character. Plot title. Default: "Gene–Concept Network" + tag.
#' @param subtitle Character. Subtitle. Default: derived from input object.
#' @param caption Character. Caption. Default: constructed from input parameters.
#' @param circular Logical. Draw network in circular layout. Default: FALSE.
#' @param colorEdge Logical. Color edges by category. Default: TRUE.
#' @param cex_label_category Numeric. Size of category labels. Default: 0.8.
#' @param cex_label_gene Numeric. Size of gene labels. Default: same as cex_label_category.
#' @param node_label Character. Which nodes to label: "all", "gene",
#'   or "category". Default: "category".
#' @param save Logical. Whether to save the plot. Default: TRUE.
#' @param also.pdf Logical. Save both png and pdf. Default: FALSE.
#' @param save.obj Logical. Whether to save the ggplot object. Default: FALSE.
#' @param w Width in inches. Default: 10.
#' @param h Height in inches. Default: 10.
#' @param ... Additional arguments passed to `enrichplot::cnetplot()`.
#'
#' @importFrom ggplot2 labs theme_void annotate
#'
#' @return A ggplot object (invisibly if saved).
#' @export
#'
#' @examples
#' \dontrun{
#' edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
#' scGeneConceptNetworkEnrichr(
#'   df.enrichment = edox,
#'   foldChange = geneList,
#'   tag = "Cluster 3 neurons"
#' )
#' }
scGeneConceptNetworkEnrichr <- function(
    df.enrichment,
    showCategory = 10,
    foldChange = NULL,
    tag = NULL,
    title = paste("Gene–Concept Network", tag),
    subtitle = kppws("Input: ", substitute_deparse(df.enrichment)),
    caption = paste0(
      "Enriched terms: ", ifelse(is.null(df.enrichment), 0, nrow(df.enrichment)),
      " | Shown: ", ifelse(is.null(df.enrichment), 0,
                           min(showCategory, nrow(df.enrichment))),
      if (!is.null(foldChange))
        paste0(" | genes w/ foldChange: ", length(foldChange))
      else ""
    ),
    circular = FALSE,
    colorEdge = TRUE,
    cex_label_category = 0.8,
    cex_label_gene   = cex_label_category,
    node_label = "category",
    save = TRUE,
    also.pdf = FALSE,
    save.obj = FALSE,
    w = 10, h = 10,
    ...
) {

  stopifnot(
    "Package 'enrichplot' must be installed." =
      requireNamespace("enrichplot", quietly = TRUE)
  )

  if(is.null(tag)) warning("Please provide a tag describing where the enrichments come from.",immediate. = TRUE)

  nr_terms <- if (is.null(df.enrichment)) 0 else nrow(df.enrichment)
  nr_GOENR_input_genes <- length(df.enrichment@"gene")

  pobj <-
    if (is.null(df.enrichment) || nr_terms < 1) {
      Seurat.utils:::.emptyAnnotatedPlot(label = "No enriched terms input!")


    } else if (!is.null(foldChange) && length(nr_GOENR_input_genes) < 5) {
      Seurat.utils:::.emptyAnnotatedPlot(label = "Too few input genes for GO enrichment (<5).")

    } else {
      enrichplot::cnetplot(
        x = df.enrichment,
        showCategory = showCategory,
        node_label = node_label,
        color.params = list(
          foldChange = foldChange,
          edge = colorEdge
        ),
        cex.params = list(
          category_label = cex_label_category,
          gene_label = cex_label_gene
        ),
        circular = circular,
        ...
      )
    }

  pobj <- pobj +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      caption = caption
    )

  if (save) {
    qqSave(
      pobj,
      title = title,
      w = w, h = h,
      also.pdf = also.pdf,
      save.obj = save.obj
    )
  }

  return(pobj)
}





# ________________________________________________________________________
#' @title Count Enriched and Depleted Genes
#'
#' @description This function counts the number of significantly enriched and depleted genes
#' based on the provided criteria. It filters the genes based on adjusted p-value and
#' logarithm of fold change.
#'
#' @param df A dataframe containing the result of the differential gene expression analysis.
#' @param min_padj A numeric value specifying the minimum adjusted p-value. Default: 0.01.
#' @param min_logFC A numeric value specifying the minimum logarithm to fold change.
#' Default: 0.5. This value should be positive and will be used as a negative value for
#' depleted genes.
#' @param colname.p A string specifying the column name for the adjusted p-value in the dataframe.
#' Default: 'p_val_adj'.
#' @param colname.lFC A string specifying the column name for the logarithm to fold change in
#' the dataframe. Default: 'avg_log2FC'.
#'
#' @return A list of two elements:
#' \item{GeneCounts}{A named numeric vector containing the numbers of enriched and depleted genes.}
#' \item{Parameters}{A named numeric vector containing the parameter names and their values.}
#'
#' @examples
#' df <- data.frame(
#'   p_val = c(5.580902e-14, 4.607790e-12, 1.643436e-11),
#'   avg_log2FC = c(0.4985875, 0.4983416, 0.4977825),
#'   pct.1 = c(0.429, 0.575, 0.387),
#'   pct.2 = c(0.251, 0.396, 0.232),
#'   p_val_adj = c(1.091513e-09, 9.011916e-08, 3.214233e-07)
#' )
#' result <- countEnrichedDepletedGenes(df)
#' print(result)
#'
#' @export

countEnrichedDepletedGenes <- function(df, min_padj = 0.01, min_logFC = 0.5,
                                       # genes = rownames(df),
                                       colname.p = "p_val_adj", colname.lFC = "avg_log2FC") {
  stopifnot(
    min_padj > 0,
    min_logFC > 0,
    colname.p %in% colnames(df),
    colname.lFC %in% colnames(df)
  )

  # Filter the dataframe for enriched genes
  idx.enr <- df[[colname.p]] <= min_padj & df[[colname.lFC]] >= min_logFC
  enriched_genes <- df[idx.enr, ]
  enriched_symbols <- rownames(enriched_genes)
  # enriched_symbols <- genes[idx.enr]

  # Filter the dataframe for depleted genes
  idx.depl <- (df[[colname.p]] <= min_padj & df[[colname.lFC]] <= -min_logFC)

  depleted_genes <- df[idx.depl, ]
  depleted_symbols <- rownames(depleted_genes)
  # depleted_symbols <- genes[idx.depl]

  # Create the named numeric vectors
  gene_counts <- c("Enriched" = nrow(enriched_genes), "Depleted" = nrow(depleted_genes))
  print("gene_counts")
  print(gene_counts)
  parameters <- c("min_padj" = min_padj, "min_logFC" = min_logFC)
  print("parameters")
  print(parameters)

  # Create the list of gene symbols
  gene_symbols <- list("Enriched" = enriched_symbols, "Depleted" = depleted_symbols)

  # Return the results as a list
  result <- list(gene_counts, parameters, gene_symbols)
  names(result) <- c("GeneCounts", "Parameters", "GeneSymbols")

  return(result)
}

# ________________________________________________________________________



# _________________________________________________________________________________________________
# Helpers ______________________________ ----
# _________________________________________________________________________________________________

#' @title Adjust Layout Parameters for multi* plotting functions
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


# ________________________________________________________________________
#' @title Empty ggplot with centered annotation and optional warning
#'
#' @description
#' Create a blank ggplot with a centered text annotation.
#' Optionally emits a warning with a custom message.
#'
#' @param label Character. Text shown in the plot.
#' @param warning_msg Character or NULL. Warning text. Default: NULL.
#' @param color Character. Text color. Default: "red".
#' @param size Numeric. Text size. Default: 8.
#'
#' @return An placeholder ggplot object with annotation.
.emptyAnnotatedPlot <- function(
    label,
    warning_msg = label,
    color = "red",
    size = 8
) {

  stopifnot(
    is.character(label), is.character(color), is.numeric(size),
    is.null(warning_msg) || (is.character(warning_msg) && length(warning_msg) == 1)
  )

  # if (!is.null(warning_msg)) warning(warning_msg, immediate. = TRUE)

  ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::annotate(
      geom  = "text",
      x  = 1, y = 1,
      label = label, color = color,size  = size,
      hjust = 0.5, vjust = 0.5
    )
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
#' if FALSE, the name is based on `plot_list` and `suffix`; Default: `FALSE`.
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
  if (pname == FALSE) pname <- sppp(substitute_deparse(plot_list), suffix)
  p1 <- cowplot::plot_grid(
    plotlist = plot_list, nrow = nrow, ncol = ncol,
    labels = LETTERS[1:length(plot_list)], ...
  )
  p1 <- cowplot::ggdraw(p1) +
    theme(plot.background = element_rect(fill = "white", color = NA))

  print("Saved as:")
  MarkdownHelpers::ww.FnP_parser(extPNG(pname))

  save_plot(plot = p1, filename = extPNG(pname), base_height = h, base_width = w)
}

# _________________________________________________________________________________________________
#' @title Save Four Plots on One A4 Page
#'
#' @description Arranges and saves four plots (e.g. UMAPs) onto a single A4 page, allowing for a
#' compact comparison of different visualizations or clustering results.
#'
#' @param plot_list A list containing ggplot objects to be arranged and saved; each object represents one panel.
#' @param pname Plot name; if FALSE, a name is generated automatically based on `plot_list` and `suffix`; Default: `FALSE`.
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
  if (pname == FALSE) pname <- sppp(substitute_deparse(plot_list), suffix)
  p1 <- cowplot::plot_grid(
    plotlist = plot_list, nrow = nrow, ncol = ncol,
    labels = LETTERS[1:length(plot_list)], ...
  )
  # https://stackoverflow.com/questions/13691415/change-the-background-color-of-grid-arrange-output
  p1 <- cowplot::ggdraw(p1) +
    theme(plot.background = element_rect(fill = "white", color = NA))

  iprint("Saved as:", pname)
  MarkdownHelpers::ww.FnP_parser(extPNG(pname))
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
  MarkdownHelpers::ww.FnP_parser(fname)
  save_plot(
    filename = fname,
    plot = pg.cf, base_height = height, base_width = width
  )
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
# try (source('~/GitHub/Packages/CodeAndRoll/CodeAndRoll.R'),silent= TRUE) # generic utilities functions
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
#'   plot3D.umap.gene(obj = combined.obj, gene = "percent.mito", quantileCutoff = .95) # for continuous meta variables
#'   plot3D.umap.gene(obj = combined.obj, gene = "nFeature_RNA", quantileCutoff = .95) # for continuous meta variables
#' }
#' }
#' @importFrom plotly plot_ly layout
#' @importFrom Seurat FetchData
#'
#' @export plot3D.umap.gene

plot3D.umap.gene <- function(
    gene = "TOP2A",
    obj = combined.obj,
    annotate.by = GetNamedClusteringRuns(obj = obj, v = FALSE)[1],
    quantileCutoff = .99,
    def.assay = c("integrated", "RNA")[2],
    suffix = NULL,
    alpha = .5,
    dotsize = 1.25,
    col.names = c("umap_1", "umap_2", "umap_3"),
    assay = "RNA",
    ...) {
  # Input assertions ____________________________________
  # browser()
  stopifnot(
    is(obj, "Seurat"),
    is.character(gene),
    "gene or feature not found in obj" = (gene %in% Features(obj, assay = assay) | gene %in% colnames(obj@meta.data)),
    "annotate.by not found in @meta" = (annotate.by %in% colnames(obj@meta.data) | annotate.by == FALSE),
    "reductions.backup is missing from @misc" = is.list(obj@misc$"reductions.backup"),
    "umap3d is missing from @misc$reductions.backup" = is(obj@misc$reductions.backup$"umap3d", class2 = "DimReduc"),
    "reductionn has 3 columns" = (ncol(obj@misc$reductions.backup$"umap3d") == 3),
    "3D reduction has >/< cells than object" = (ncol(obj) == nrow(obj@misc$reductions.backup$"umap3d"@cell.embeddings))
  )
  # browser()

  if (obj@version < "5") col.names <- toupper(col.names)
  message("Obj. version: ", obj@version, " \ndim names: ", kppc(col.names))

  DefaultAssay(object = obj) <- def.assay
  iprint(DefaultAssay(object = obj), "assay")

  # Get and format 3D plotting data ____________________________________
  plotting.data <- obj@misc$reductions.backup$"umap3d"@cell.embeddings
  colnames(plotting.data) <- toupper(col.names)


  Expression <- Seurat::FetchData(object = obj, vars = gene)
  plotting.data <- cbind(plotting.data, Expression)

  plotting.data$"Expression" <- ww.check.quantile.cutoff.and.clip.outliers(
    expr.vec = plotting.data[, gene],
    quantileCutoffX = quantileCutoff, min.cells.expressing = 10
  )

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
    marker = list(
      size = dotsize,
      color = ~Expression,   # Map Expression to color
      colorscale = "Viridis",
      opacity = alpha
    ),
    text = ~label,
    ...
  )  |>
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
#' @export plot3D.umap

plot3D.umap <- function(
    obj = combined.obj,
    category = GetNamedClusteringRuns(obj = obj, v = FALSE)[1],
    annotate.by = category,
    suffix = NULL,
    dotsize = 1.25,
    col.names = c("umap_1", "umap_2", "umap_3"),
    ...) {
  message("category: ", category)
  message("annotate.by: ", annotate.by)



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

  if (obj@version < "5") col.names <- toupper(col.names)
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
  ) |>
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
#' @title SetupReductionsNtoKdimensions
#'
#' @description Function to calculate N-to-K dimensional umaps (default = 2:3); and back them up to
#' slots `obj@misc$reductions.backup` from @reductions$umap
#' @param obj A Seurat object. Default: combined.obj
#' @param nPCs A numeric value representing the number of principal components to use. Default: p$n.PC
#' @param dimensions A numeric vector specifying the dimensions to use for the dimensionality reductions. Default: 3:2
#' @param reduction_input The type of dimensionality reduction to use as input. Can be "pca", or
#' some correction results, like harmony pca. Default: 'pca'
#' @param reduction_output The type of dimensionality reduction to perform.  Can be "umap", "tsne",
#' "pca", or some correction results, like harmony pca.  Default: 'umap'
#' @param ... Additional arguments to pass to the dimensionality reduction function.
#' @return The input Seurat object with computed dimensionality reductions and backups of these reductions.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   combined.obj <- SetupReductionsNtoKdimensions(obj = combined.obj, nPCs = 10, dimensions = 2:3, reduction = "umap")
#' }
#' }
#' @importFrom tictoc tic toc
#'
#' @export
SetupReductionsNtoKdimensions <- function(obj = combined.obj, nPCs = p$"n.PC", dimensions = 3:2,
                                          reduction_input = "pca",
                                          reduction_output = "umap", ...) {
  cat("Starting dimensional reduction setup\n")
  tictoc::tic()

  for (d in dimensions) {
    iprint("Calculating", d, "dimensional", reduction_output)

    # Assign the reduction based on the output type
    obj <- switch(reduction_output,
                  umap = RunUMAP(obj, dims = 1:nPCs, reduction = reduction_input, n.components = d, ...),
                  tsne = RunTSNE(obj, dims = 1:nPCs, reduction = reduction_input, n.components = d, ...),
                  pca = RunPCA(obj, dims = 1:nPCs, n.components = d, ...)
    )

    obj <- BackupReduction(obj = obj, dim = d, reduction = reduction_output)
  }

  tictoc::toc()
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
    plotting.data. |>
    group_by(annot) |>
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
#' @export Plot3D.ListOfGenes
Plot3D.ListOfGenes <- function(
    obj = combined.obj # Plot and save list of 3D UMAP or tSNE plots using plotly.
    , annotate.by = "integrated_snn_res.0.7", opacity = 0.5, cex = 1.25, default.assay = c("integrated", "RNA")[2],
    ListOfGenes = c("BCL11B", "FEZF2", "EOMES", "DLX6-AS1", "HOPX", "DDIT4"),
    SubFolderName = ppp("plot3D", substitute_deparse(ListOfGenes))) {
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
#' @param obj A Seurat object for which the plot is to be created. Default: 'combined.obj'.
#' @param annotate.by Character vector specifying the metadata column to be used for annotating
#' the plot. Default: 'integrated_snn_res.0.7'.
#' @param cex Numeric value specifying the point size on the plot. Default: 1.25.
#' @param default.assay Character vector specifying the assay to be used. Default: 'RNA'
#' (second element in the vector c("integrated", "RNA")).
#' @param ListOfCategories Character vector specifying the categories to be included in the plot.
#' Default categories are "v.project", "experiment", "Phase", "integrated_snn_res.0.7".
#' @param SubFolderName String specifying the name of the subfolder where the plots will be saved.
#' By default, it's created using the function ppp("plot3D", substitute_deparse(ListOfCategories)).
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   categ3Dplots <- c("v.project", "experiment", "Phase", "integrated_snn_res.0.7", "Area", "Individual", "Type")
#'   Plot3D.ListOfCategories(obj = combined.obj, ListOfCategories = categ3Dplots)
#' }
#' }
#' @export Plot3D.ListOfCategories
Plot3D.ListOfCategories <- function(
    obj = combined.obj # Plot and save list of 3D UMAP or tSNE plots using plotly.
    , annotate.by = "integrated_snn_res.0.7", cex = 1.25, default.assay = c("integrated", "RNA")[2],
    ListOfCategories = c("v.project", "experiment", "Phase", "integrated_snn_res.0.7"),
    SubFolderName = ppp("plot3D", substitute_deparse(ListOfCategories))) {
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
#' @param repel A logical value indicating whether to repel the labels to avoid overlap. Default: `TRUE`.
#' @param plotWidth Numeric value specifying the width of the plot when saved. Default: 7.
#' @param plotHeight Numeric value specifying the height of the plot when saved. Default: 5.
#' @param save A logical value indicating whether to save the plot to a PDF file. Default: `TRUE`.
#' @param suffix A string suffix to append to the plot filename. Default: NULL.
#' @param assay The assay to use for the plot. Default: DefaultAssay(obj).
#' @param ... Additional arguments to pass to the Seurat::VariableFeaturePlot function.
#'
#' @examples
#' \dontrun{
#' suPlotVariableFeatures(combined.obj)
#' }
#' @export
suPlotVariableFeatures <- function(obj = combined.obj, NrVarGenes = 15,
                                   repel = TRUE, plotWidth = 7, plotHeight = 5, save = TRUE,
                                   # suffix = kpp("nVF", .getNrScaledFeatures(obj)),
                                   assay = DefaultAssay(obj),
                                   suffix = NULL,
                                   ...) {
  message(" > Running suPlotVariableFeatures()...")
  message(length(Cells(obj)), " cells | assay: ", assay, " | NrVarGenes: ", NrVarGenes)

  stopifnot(
    is(obj, "Seurat"), is.function(ppp), is.logical(repel),
    is.numeric(plotWidth), is.numeric(plotHeight)
  )

  obj.name <- substitute_deparse(obj)

  plot1 <- Seurat::VariableFeaturePlot(obj, assay = assay, ...) +
    theme(panel.background = element_rect(fill = "white")) +
    labs(
      title = "Variable Genes",
      subtitle = kppws(obj.name, suffix),
      caption = paste("Assay:", assay, "|", idate())
    )


  # Assuming LabelPoints is defined elsewhere and available for use.
  TopVarGenes <- VariableFeatures(obj, assay = assay)[1:NrVarGenes]
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


# Notes --------------------------------------------------------------------------------------------

# plotMetadataCategPie() is in Seurat.Utils.Metadata.R
