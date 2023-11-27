# ____________________________________________________________________
# Seurat.utils ----
# ____________________________________________________________________
# source("~/GitHub/Packages/Seurat.utils/R/Seurat.Utils.R")
# devtools::load_all("~/GitHub/Packages/Seurat.utils")
# devtools::document("~/GitHub/Packages/Seurat.utils"); devtools::load_all("~/GitHub/Packages/Seurat.utils")

# _________________________________________________________________________________________________
# Cluster.Auto-naming.DE.R
# _________________________________________________________________________________________________
# source('~/GitHub/Packages/Seurat.utils/Functions/Cluster.Auto-naming.DE.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Cluster.Auto-naming.DE.R"))

# _________________________________________________________________________________________________
# require(princurve) # only for AutoNumber.by.PrinCurve




# _________________________________________________________________________________________________
# plotting.filtering.R ______________________________ ----
# ____________________________________________________________________
# source('~/GitHub/Packages/Seurat.utils/Functions/plotting.filtering.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Plotting.filtering.R"))


# _________________________________________________________________________________________________
#' @title PlotFilters
#'
#' @description Plot filtering threshold and distributions, using four panels to highlight the relation between Gene- and UMI-count, ribosomal- and mitochondrial-content.
#' @param ls.obj A list of Seurat objects to be analyzed. Default is ls.Seurat.
#' @param parentdir A string representing the parent directory where the plots will be stored. Default is OutDirOrig.
#' @param suffices A vector of strings that will be used as suffixes in the output plot file names. Default is the names of the Seurat objects in ls.obj.
#' @param filetype A string indicating the file type of the output plot images. Default is '.png'.
#' @param below.mito Numeric threshold for the lower bound of mitochondrial content. Default is p$thr.lp.mito.
#' @param above.mito Numeric threshold for the upper bound of mitochondrial content. Default is p$thr.hp.mito.
#' @param below.ribo Numeric threshold for the lower bound of ribosomal content. Default is p$thr.lp.ribo.
#' @param above.ribo Numeric threshold for the upper bound of ribosomal content. Default is p$thr.hp.ribo.
#' @param below.nFeature_RNA Numeric threshold for the lower bound of RNA features. Default is p$thr.lp.nFeature_RNA.
#' @param above.nFeature_RNA Numeric threshold for the upper bound of RNA features. Default is p$thr.hp.nFeature_RNA.
#' @param subdir A string specifying the subdirectory within the parent directory where the plots will be stored. Default is generated using a call to kpp().
#' @param transparency Numeric value controlling the transparency of points on the scatter plots. Default is 0.25.
#' @param cex Numeric value controlling the size of points on the scatter plots. Default is 0.75.
#' @param theme.used A ggplot2 theme to be used for all plots. Default is theme_bw(base_size = 18).
#' @param LabelDistFromTop Numeric value specifying the distance from the top of the plot for the label placement. Default is 200.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   PlotFilters(ls.Seurat)
#' }
#' }
#' @seealso
#'  \code{\link[ggplot2]{ggplot}}, \code{\link[ggplot2]{labs}}, \code{\link[ggplot2]{geom_point}}
#' @importFrom ggplot2 ggplot ggtitle geom_point
#' @importFrom Stringendo percentage_formatter
#' @importFrom MarkdownHelpers llprint create_set_OutDir
#' @importFrom cowplot plot_grid
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
    subdir = kpp(
      "Filtering.plots",
      "mito", p$"thr.hp.mito", p$"thr.lp.mito",
      "ribo", p$"thr.hp.ribo", p$"thr.lp.ribo",
      "nFeature", p$"thr.hp.nFeature_RNA", below.nFeature_RNA, "/"
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
  MarkdownReports::create_set_OutDir(parentdir, subdir)
  # if (length(suffices) == length(ls.obj)) print("ls.obj elements have no names (required).")
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
      print('objX <- add.meta.fraction(obj = objX, col.name = "percent.mito", gene.symbol.pattern =  "^MT\\.|^MT-")')
      print('objX <- add.meta.fraction(obj = objX, col.name = "percent.ribo", gene.symbol.pattern =  "^RPL|^RPS")')
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
    fname <- ppp("Filtering.thresholds", suffices[i], filetype)
    save_plot(filename = fname, plot = px, base_height = 12, ncol = 1, nrow = 1) # Figure 2
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
#' @title PCA percent of variation associated with each PC
#'
#' @description Determine percent of variation associated with each PC. For normal prcomp objects, see: PCA.percent.var.explained().
#' @param obj Seurat object, Default: combined.obj
#' @export
seu.PC.var.explained <- function(obj = combined.obj) { # Determine percent of variation associated with each PC.
  pct <- obj@reductions$pca@stdev / sum(obj@reductions$pca@stdev) * 100
  names(pct) <- 1:length(obj@reductions$pca@stdev)
  return(pct)
}

# _________________________________________________________________________________________________
#' @title seu.plot.PC.var.explained
#'
#' @description Plot the percent of variation associated with each PC.
#' @param obj Seurat object, Default: combined.obj
#' @param use.MarkdownReports Use MarkdownReports for plotting, Default: FALSE
#' @importFrom MarkdownReports wbarplot
#' @importFrom ggExpress qbarplot
#' @export
seu.plot.PC.var.explained <- function(obj = combined.obj, use.MarkdownReports = FALSE) { # Plot the percent of variation associated with each PC.
  pct <- seu.PC.var.explained(obj)
  if (use.MarkdownReports) {
    MarkdownReports::wbarplot(pct, xlab = "Principal Components", ylab = "% of variation explained")
    barplot_label(round(pct, digits = 2), barplotted_variable = pct, cex = .5)
  } else {
    ggExpress::qbarplot(vec = pct, xlab = "Principal Components", ylab = "% of variation explained", w = 10, h = 5, hline = 1)
  }
}



# _________________________________________________________________________________________________
#' @title Percent.in.Trome
#'
#' @description Gene expression as fraction of all UMI's
#' @param obj Seurat object
#' @param n.genes.barplot number of top genes shows
#' @param width.barplot barplot width
#' @return Seurat object
#' @examples # combined.obj <- Percent.in.Trome()

#' @export
Percent.in.Trome <- function(
    obj = combined.obj, n.genes.barplot = 25,
    width.barplot = round(n.genes.barplot / 4)) {
  m.expr <- combined.obj@assays$RNA@counts
  total.Expr <- sort(rowSums(m.expr), decreasing = TRUE)
  relative.total.Expr <- total.Expr / sum(total.Expr)
  print(head(iround(100 * relative.total.Expr), n = n.genes.barplot))

  qhistogram(relative.total.Expr * 100,
    logX = FALSE, logY = TRUE,
    plotname = "Gene expression as fraction of all UMI's",
    subtitle = "Percentage in RNA-counts",
    xlab = "Percent in Transcriptome (total per gene)",
    ylab = "Number of genes",
    xlab.angle = 45
  ) # + geom_hline(yintercept = 10)

  Highest.Expressed.Genes <- head(iround(100 * relative.total.Expr), n = n.genes.barplot)
  qbarplot(Highest.Expressed.Genes,
    w = width.barplot,
    plotname = "Percentage of highest expressed genes",
    subtitle = "Total, in RNA-counts",
    xlab = "",
    ylab = "Gene expression as percent of all UMI's",
    xlab.angle = 45
  )
  print("!!!")
  print("TotalReadFraction is stored under combined.obj@misc$'TotalReadFraction'  ")
  print("!!!")
  combined.obj@misc$"TotalReadFraction" <- relative.total.Expr
  return(combined.obj)
}



# _________________________________________________________________________________________________
#' @title gene.expression.level.plots
#'
#' @description Histogram of gene expression levels.
#' @param gene gene of interest, Default: 'TOP2A'
#' @param obj Seurat object, Default: ls.Seurat[[1]]
#' @param slot slot in the Seurat object. Default: c("counts", "data")[2]
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @export
gene.expression.level.plots <- function(
    gene = "TOP2A", obj = ls.Seurat[[1]],
    slot = c("counts", "data")[2],
    ...) {
  print(gene)
  if (gene %in% rownames(obj)) {
    GEX.Counts <- GetAssayData(object = obj, assay = "RNA", slot = slot)

    GEX.Counts.total <- rowSums(GEX.Counts)
    genes.expression <- GEX.Counts.total[gene]
    mean.expr <- iround(mean(GEX.Counts[gene, ]))

    suffx <- if (slot == "counts") "raw" else "normalised, logtransformed"
    (pname <- paste(gene, "and the", suffx, "transcript count distribution"))

    ggExpress::qhistogram(GEX.Counts.total,
      vline = genes.expression, logX = TRUE, w = 7, h = 4,
      subtitle = paste("It belong to the top", pc_TRUE(GEX.Counts.total > genes.expression), "of genes (black line). Mean expr:", mean.expr),
      plotname = pname, xlab = "Total Transcripts in Dataset", ylab = "Number of Genes",
      ...
    )
  } else {
    print("     !!! Gene not found in object!")
  }
}

# _________________________________________________________________________________________________
#' @title PrctCellExpringGene
#'
#' @description Function to calculate the proportion of cells expressing a given set of genes.
#' @param genes A character vector of genes of interest.
#' @param group.by Grouping variable, Default: 'all'.
#' @param obj A Seurat object containing cell data. Default: combined.obj.
#' @return A data frame with the proportion of cells expressing each gene, grouped by the group.by variable.
#' @examples
#' \dontrun{
#' PrctCellExpringGene(genes = c("Gene1", "Gene2"), obj = seurat_object)
#' }
#' @source Adapted from code by Ryan-Zhu on Github (https://github.com/satijalab/seurat/issues/371)
#' @export
PrctCellExpringGene <- function(genes, group.by = "all", obj = combined.obj) { # From Github/Ryan-Zhu https://github.com/satijalab/seurat/issues/371
  if (group.by == "all") {
    prct <- unlist(lapply(genes, ww.calc_helper, object = obj))
    result <- data.frame(Markers = genes, Cell_proportion = prct)
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
#' @title ww.calc_helper
#'
#' @description Helper function for PrctCellExpringGene() to calculate the proportion of cells in a Seurat object that express a given gene.
#' @param obj A Seurat object containing cell data.
#' @param genes A character vector of genes of interest.
#' @return The proportion of cells in obj that express the specified gene.
#' @examples
#' \dontrun{
#' ww.calc_helper(obj = seurat_object, genes = "Gene1")
#' }
#' @source Adapted from code by Ryan-Zhu on Github (https://github.com/satijalab/seurat/issues/371)
#' @export
ww.calc_helper <- function(obj, genes) {
  counts <- obj[["RNA"]]@counts
  ncells <- ncol(counts)
  if (genes %in% row.names(counts)) {
    sum(counts[genes, ] > 0) / ncells
  } else {
    return(NA)
  }
}

# _________________________________________________________________________________________________
#' @title scBarplot.FractionAboveThr
#'
#' @description Create a bar plot showing the fraction of cells, within each cluster, that exceed a certain threshold based on a metadata column.
#' @param thrX The threshold value to determine the fraction of cells. Default: 0.3
#' @param value.col Column name from metadata which holds the values for calculating the fraction of cells. Default: 'percent.ribo'
#' @param id.col Column name from metadata to be used for identifying clusters. Default: 'cl.names.top.gene.res.0.3'
#' @param obj A Seurat object holding single cell data. Default: combined.obj
#' @param suffix An optional suffix for the filename.
#' @param return.df A logical indicating if the function should return the data frame used to create the plot. Default: FALSE
#' @param label A logical indicating if labels should be added to the bar plot. Default: FALSE
#' @param subtitle Subtitle
#' @param ... Additional parameters to pass to internally called functions.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   scBarplot.FractionAboveThr(id.col = "cl.names.top.gene.res.0.3", value.col = "percent.ribo", thrX = 0)
#' }
#' }
#' @seealso
#'  \code{\link[dplyr]{select}}, \code{\link[dplyr]{se-deprecated}}
#' @importFrom dplyr select group_by_
#'
#' @export
scBarplot.FractionAboveThr <- function(
    thrX = 0.3, value.col = "percent.ribo",
    id.col = "cl.names.top.gene.res.0.3",
    obj = combined.obj, return.df = FALSE, label = FALSE,
    suffix = NULL, subtitle = id.col,
    ...) { # Calculat the fraction of cells per cluster above a certain threhold
  meta <- obj@meta.data
  metacol <- meta %>%
    dplyr::select(c(id.col, value.col))

  (df_cells_above <- metacol %>%
    dplyr::group_by_(id.col) %>%
    summarize(
      n_cells = n(),
      n_cells_above = sum(!!as.name(value.col) > thrX),
      fr_n_cells_above = n_cells_above / n_cells
    )
  )

  total_average <- iround(100 * mean(metacol[, value.col] > thrX))

  df_2vec <- df_cells_above[, c(1, 4)]
  (v.fr_n_cells_above <- 100 * deframe(df_2vec))
  if (label == TRUE) lab <- percentage_formatter(deframe(df_2vec), digitz = 2) else lab <- NULL

  pname <- paste("Pc. cells above", value.col, "of", thrX)
  ggobj <- ggExpress::qbarplot(v.fr_n_cells_above,
    plotname = pname,
    filename = FixPlotName(kpp(pname, id.col, ".pdf")),
    suffix = suffix,
    subtitle = subtitle,
    caption = paste(
      "Overall average:", iround(total_average), "% |",
      substitute(obj)
    ) # , '\n', id.col
    , xlab.angle = 45,
    xlab = "Clusters", ylab = paste("% Cells above thr. (", value.col, ")"),
    label = lab,
    hline = total_average,
    ...
  )
  if (return.df) {
    return(df_cells_above)
  } else {
    ggobj
  }
}


# _________________________________________________________________________________________________
#' @title scBarplot.FractionBelowThr
#'
#' @description Create a bar plot showing the fraction of cells, within each cluster, that are below a certain threshold based on a metadata column.
#' @param thrX The threshold value to determine the fraction of cells. Default: 0.01
#' @param value.col Column name from metadata which holds the values for calculating the fraction of cells. Default: 'percent.ribo'
#' @param id.col Column name from metadata to be used for identifying clusters. Default: 'cl.names.top.gene.res.0.3'
#' @param obj A Seurat object holding single cell data. Default: combined.obj
#' @param return.df A logical indicating if the function should return the data frame used to create the plot. Default: FALSE
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   scBarplot.FractionBelowThr(id.col = "cl.names.top.gene.res.0.3", value.col = "percent.ribo", thrX = 0.01, return.df = TRUE)
#' }
#' }
#' @seealso
#'  \code{\link[dplyr]{select}}, \code{\link[dplyr]{se-deprecated}}
#' @importFrom dplyr select group_by_
#'
#' @export
scBarplot.FractionBelowThr <- function(
    thrX = 0.01, value.col = "percent.ribo", id.col = "cl.names.top.gene.res.0.3",
    obj = combined.obj, return.df = FALSE) { # Calculat the fraction of cells per cluster below a certain threhold
  meta <- obj@meta.data
  (df_cells_below <- meta %>%
    dplyr::select(c(id.col, value.col)) %>%
    dplyr::group_by_(id.col) %>%
    summarize(
      n_cells = n(),
      n_cells_below = sum(!!as.name(value.col) < thrX),
      fr_n_cells_below = n_cells_below / n_cells
    ) %>%
    FirstCol2RowNames())

  (v.fr_n_cells_below <- 100 * as.named.vector.df(df_cells_below[3]))

  pname <- make.names(paste("Cells with", value.col, "<", thrX, id.col))
  ggobj <- ggExpress::qbarplot(v.fr_n_cells_below,
    xlab = "Clusters", ylab = "% Cells",
    plotname = pname,
    subtitle = id.col, xlab.angle = 45
  )
  if (return.df) {
    return(df_cells_below)
  } else {
    ggobj
  }
}




# _________________________________________________________________________________________________
# Colors ______________________________ ----
# _________________________________________________________________________________________________
#' @title gg_color_hue
#'
#' @description Emulates the default color palette of ggplot2. Source:  https://stackoverflow.com/questions/8197559/emulate-ggplot2-default-color-palette
#' @param n The number of colors to generate.
#' @return A vector of colors emulating the default color palette of ggplot2.
#' @examples
#' \dontrun{
#' gg_color_hue(5)
#' }
#' @export gg_color_hue

gg_color_hue <- function(n) { # reproduce the ggplot2 default color palette
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}



# _________________________________________________________________________________________________
#' @title getDiscretePalette
#'
#' @description Generate a discrete color palette.
#' @param ident.used The identity column used for determining the number of clusters, Default: GetClusteringRuns()[1]
#' @param obj Seurat object, Default: combined.obj
#' @param palette.used The name of the palette to use, Default: c("alphabet", "alphabet2", "glasbey", "polychrome", "stepped")[1]
#' @param show.colors Whether to display the colors in the palette, Default: FALSE
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   getDiscretePalette()
#' }
#' }
#' @importFrom MarkdownHelpers MarkdownHelpers::color_check
#'
#' @export getDiscretePalette
getDiscretePalette <- function(
    ident.used = GetClusteringRuns()[1],
    obj = combined.obj,
    palette.used = c("alphabet", "alphabet2", "glasbey", "polychrome", "stepped")[1],
    show.colors = FALSE) {
  n.clusters <- nrow(unique(obj[[ident.used]]))
  colz <- DiscretePalette(n = n.clusters, palette = palette.used)
  if (anyNA(colz)) {
    colzOK <- na.omit.strip(colz)
    repNeeded <- ceiling(length(colz) / length(colzOK))
    colzFixed <- rep(colzOK, repNeeded)[1:length(colz)]
    stopif(anyNA(colzFixed))
    colz <- colzFixed
  }
  if (show.colors) MarkdownHelpers::color_check(colz)
  return(colz)
}



# _________________________________________________________________________________________________
#' @title getClusterColors
#'
#' @description get Seurat's cluster colors.
#' @param obj Seurat object, Default: combined.obj
#' @param ident identity used, Default: GetClusteringRuns()[1]
#' @param show Show plot of colors? Default: TRUE
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   getClusterColors(obj = combined.obj, ident = GetClusteringRuns()[2])
#' }
#' }
#' @seealso
#'  \code{\link[scales]{hue_pal}}
#' @export
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
  # color_check(color_palette)
  # names(color_palette) <- sort(as.factor(identities))
  names(color_palette) <- (identities)
  identvec <- obj[[ident]][, 1]
  colz <- color_palette[identvec]
  names(colz) <- identvec
  if (show) color_check(unique(colz)) # MarkdownReports
  colz
}



# _________________________________________________________________________________________________
#' @title SeuratColorVector
#'
#' @description Recall a Seurat color vector.
#' @param ident identity used, Default: NULL
#' @param obj Seurat object, Default: combined.obj
#' @param plot.colors Show colors? Default: FALSE
#' @param simple Return simply the unique colors, in order? Default: FALSE
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   SeuratColorVector()
#'   SeuratColorVector(ident = GetNamedClusteringRuns()[2], plot.colors = TRUE)
#' }
#' }
#' @seealso
#'  \code{\link[scales]{hue_pal}}
#' @export
#' @importFrom scales hue_pal

SeuratColorVector <- function(ident = NULL, obj = combined.obj, plot.colors = FALSE, simple = FALSE) {
  if (!is.null(ident)) {
    print(ident)
    ident.vec <- obj[[ident]][, 1]
  } else {
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
# plotting generic, misc ______________________________ ----
# _________________________________________________________________________________________________


# _________________________________________________________________________________________________
#' @title qFeatureScatter
#'
#' @description Generates and optionally saves a scatter plot of two features from a Seurat object.
#' @param feature1 Name of the first feature to plot. Default: 'TOP2A'.
#' @param feature2 Name of the second feature to plot. Default: 'ID2'.
#' @param obj Seurat object from which to extract feature data. Default: combined.obj.
#' @param ext File extension for the saved plot. Default: 'png'.
#' @param logX Logical indicating whether to apply a logarithmic transformation to the x-axis. Default: FALSE.
#' @param logY Logical indicating whether to apply a logarithmic transformation to the y-axis. Default: FALSE.
#' @param plot Logical indicating whether to display the plot. Default: TRUE.
#' @param ... Additional arguments to pass to the FeatureScatter function.
#' @return If `plot` is TRUE, a ggplot object containing the feature scatter plot is returned.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   qFeatureScatter(feature1 = "Gene1", feature2 = "Gene2", obj = seuratObject)
#' }
#' }
#' @seealso
#' \code{\link[Seurat]{FeatureScatter}}, \code{\link[ggplot2]{ggplot}}
#' @export

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
#' qSeuViolin
#'
#' This function creates a violin plot of a single feature in a Seurat object, split by a grouping variable.
#' @param object A Seurat object.
#' @param suffix A string to append to the title of the plot.
#' @param features The name of the feature to plot.
#' @param split.by The grouping variable to split the plot by.
#' @param logY Whether to plot the y-axis on a log scale.
#'
#' @return A ggplot object.
#'
#' @export

qSeuViolin <- function(
    object = ls.Seurat[[1]], suffix = GEX_library,
    features = "nFeature_RNA", split.by = "orig.ident", logY = TRUE) {
  # Create a violin plot of the feature, split by the grouping variable.
  p <- VlnPlot(object = object, features = features, split.by = split.by) +
    ggtitle(label = features, subtitle = paste(suffix, "by", split.by)) +
    theme(axis.title.x = element_blank()) + labs(y = "Top UVI's depth")

  # If `logY` is TRUE, plot the y-axis on a log scale.
  if (logY) p <- p + ggplot2::scale_y_log10()

  # Save the plot.
  title_ <- ppp(as.character(features), suffix, flag.nameiftrue(logY))
  qqSave(p, title = title_, w = 7, h = 5)
  p
}



# _________________________________________________________________________________________________
#' plotGeneExpHist
#'
#' This function creates a histogram of gene expression for a given set of genes in a Seurat object.
#' @param obj A Seurat object.
#' @param genes A vector of genes to plot.
#' @param assay The name of the assay to use.
#' @param slot_ The slot to use.
#' @param thr_expr The expression threshold to use for filtering.
#' @param suffix A string to append to the title of the plot.
#' @param xlab The x-axis label.
#' @param return_cells_passing Whether to return the number of cells passing the filter.
#' @param quantile_thr The quantile to use for clipping the counts slot.
#' @param return_quantile Whether to return the number of cells passing the quantile filter.
#' @param ... Additional arguments passed to `qhistogram()`.
#'
#' @return A ggplot object.
#' @importFrom MarkdownHelpers filter_HP
#'
#' @export
plotGeneExpHist <- function(
    obj = cobj.H9.L92, genes = c("MALAT1", "MT-CO1", "MT-CO2", "MT-CYB", "TMSB4X", "KAZN"),
    assay = "RNA", slot_ = "data",
    thr_expr = 10,
    suffix = NULL,
    xlab = paste0("log10(Summed UMI count @", slot_, ")"),
    return_cells_passing = TRUE,
    quantile_thr = 0.95,
    return_quantile,
    ...) {
  # Check arguments
  stopifnot(length(genes) > 0)
  stopifnot(slot_ %in% c("data", "counts"))

  # Aggregate genes if necessary
  aggregate <- length(genes) > 1
  G_expression <- colSums(GetAssayData(object = obj, assay = assay, slot = slot_)[genes, ])

  # Add a subtitle with the number of genes and the expression threshold
  subx <- filter_HP(G_expression, threshold = thr_expr, return_conclusion = TRUE, plot.hist = FALSE)
  if (aggregate) subx <- paste0(subx, "\n", length(genes), " aggregated:", paste(head(genes), collapse = " "))

  # Clip counts if necessary
  if (slot_ == "counts") G_expression <- CodeAndRoll2::clip.at.fixed.value(distribution = G_expression, thr = quantile(G_expression, probs = .95))

  # Create the plot
  title_ <- paste("Gene expression histogram", Stringendo::flag.nameiftrue(aggregate, prefix = "- "), suffix, slot_)
  pobj <- qhistogram(G_expression,
    plotname = title_,
    suffix = suffix,
    vline = thr_expr, filtercol = -1,
    xlab = xlab,
    ylab = "# of cells",
    subtitle = subx,
    caption = paste("cutoff at", iround(thr_expr)),
    w = 7, h = 6,
    ...
  )

  # Print the plot
  print(pobj)

  # Return the number of cells passing the filter
  if (return_cells_passing) {
    return(MarkdownHelpers::filter_HP(G_expression, threshold = thr_expr, plot.hist = FALSE))
  }
}



# _________________________________________________________________________________________________
# plotting.dim.reduction.2D.R ______________________________ ----
# _________________________________________________________________________________________________
# source('~/GitHub/Packages/Seurat.utils/Functions/plotting.dim.reduction.2D.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Plotting.dim.reduction.2D.R"))
# Source: self + web

# Requirements __________________________________________
# library(plotly)
# try(source("~/GitHub/Packages/ggExpressDev/ggExpress.functions.R"), silent = TRUE)
# try(source("https://raw.githubusercontent.com/vertesy/ggExpressDev/main/ggExpress.functions.R"), silent = TRUE)

# May also require
# try (source('/GitHub/Packages/CodeAndRoll/CodeAndRoll.R'),silent= FALSE) # generic utilities funtions
# require('MarkdownReports') # require("devtools")


# _________________________________________________________________________________________________
#' @title qUMAP
#'
#' @description The quickest way to draw a gene expression UMAP.
#' @param feature Feature to be visualized on the UMAP, Default: 'TOP2A'.
#' @param obj Seurat object containing single-cell RNA seq data, Default: combined.obj.
#' @param title Title of the plot, Default: feature.
#' @param sub Subtitle of the plot, Default: NULL.
#' @param reduction Dimension reduction technique to be used. Choose from 'umap', 'tsne', or 'pca'. Default: 'umap'.
#' @param splitby Column in the metadata to split the cells by, Default: NULL.
#' @param prefix A prefix added before the filename, Default: NULL.
#' @param suffix A suffix added to the end of the filename, Default: subtitle.
#' @param save.plot If TRUE, the plot is saved into a file, Default: TRUE.
#' @param PNG If TRUE, the file is saved as a .png, Default: TRUE.
#' @param h Height of the plot in inches, Default: 7.
#' @param w Width of the plot in inches, Default: NULL.
#' @param nr.cols Number of columns to combine multiple feature plots, ignored if split.by is not NULL, Default: NULL.
#' @param assay Which assay to use, 'RNA' or 'integrated', Default: 'RNA'.
#' @param axes If TRUE, the axes are shown on the plot. Default: FALSE.
#' @param aspect.ratio Ratio of height to width. If TRUE, the ratio is fixed at 0.6. Default: FALSE.
#' @param HGNC.lookup If TRUE, the HGNC gene symbol lookup is performed. Default: TRUE.
#' @param make.uppercase If TRUE, feature names are converted to uppercase. Default: FALSE.
#' @param qlow Lower quantile for the color scale, Default: 'q10'.
#' @param qhigh Upper quantile for the color scale, Default: 'q90'.
#' @param check_for_2D If TRUE, checks if UMAP is 2 dimensional. Default: TRUE.
#' @param caption Add caption to the ggplot object (e.g. a description in bottom right).
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   qUMAP("nFeature_RNA")
#'   qUMAP("TOP2A")
#' }
#' }
#' @importFrom MarkdownHelpers TRUE.unless
#'
#' @export
qUMAP <- function(
    feature = "TOP2A", obj = combined.obj # The quickest way to draw a gene expression UMAP
    , title = feature, sub = NULL,
    reduction = "umap", splitby = NULL,
    prefix = NULL,
    suffix = make.names(sub),
    save.plot = MarkdownHelpers::TRUE.unless("b.save.wplots"),
    PNG = TRUE,
    h = 7, w = NULL, nr.cols = NULL,
    assay = c("RNA", "integrated")[1],
    axes = FALSE,
    aspect.ratio = c(FALSE, 0.6)[1],
    HGNC.lookup = TRUE,
    make.uppercase = FALSE,
    check_for_2D = TRUE,
    qlow = "q10", qhigh = "q90",
    caption = FALSE,
    ...) {
  if (check_for_2D) {
    umap_dims <- ncol(obj@reductions[[reduction]]@cell.embeddings)
    if (umap_dims != 2) warning(">>> UMAP is not 2 dimensional! \n Check obj@reductions[[reduction]]@cell.embeddings")
  }


  if (!(feature %in% colnames(obj@meta.data) | feature %in% rownames(obj))) {
    feature <- check.genes(
      list.of.genes = feature, obj = obj, verbose = FALSE,
      HGNC.lookup = HGNC.lookup, makeuppercase = make.uppercase
    )
  }

  DefaultAssay(obj) <- assay
  ggplot.obj <- Seurat::FeaturePlot(obj,
    features = feature,
    reduction = reduction,
    min.cutoff = qlow, max.cutoff = qhigh
    # , plotname = ppp(toupper(reduction), feature)
    , ncol = nr.cols,
    split.by = splitby,
    ...
  ) +
    ggtitle(label = title, subtitle = sub) +
    if (!axes) NoAxes() else NULL

  if (aspect.ratio) ggplot.obj <- ggplot.obj + ggplot2::coord_fixed(ratio = aspect.ratio)
  if (!isFALSE(caption)) ggplot.obj <- ggplot.obj + labs(caption = caption)

  if (save.plot) {
    fname <- ww.FnP_parser(Stringendo::sppp(prefix, toupper(reduction), feature, assay, suffix), if (PNG) "png" else "pdf")
    try(save_plot(filename = fname, plot = ggplot.obj, base_height = h, base_width = w)) # , ncol = 1, nrow = 1
  }
  return(ggplot.obj)
}



# _________________________________________________________________________________________________
#' @title clUMAP
#'
#' @description The quickest way to draw a clustering result UMAP.
#' @param ident Identity to be used for clustering, Default: 'integrated_snn_res.0.5'.
#' @param obj Seurat object containing single-cell RNA seq data, Default: combined.obj.
#' @param reduction Dimension reduction technique to be used. Choose from 'umap', 'tsne', or 'pca'. Default: 'umap'.
#' @param splitby Column in the metadata to split the cells by, Default: NULL.
#' @param title Title of the plot, Default: ident.
#' @param sub Subtitle of the plot, Default: NULL.
#' @param prefix A prefix added before the filename, Default: NULL.
#' @param suffix A suffix added to the end of the filename, Default: sub.
#' @param label.cex Scaling factor for label sizes, Default: 7.
#' @param h Height of the plot in inches, Default: 7.
#' @param w Width of the plot in inches, Default: NULL.
#' @param nr.cols Number of columns to combine multiple feature plots, ignored if split.by is not NULL, Default: NULL.
#' @param plotname Title of the plot, Default: ppp(toupper(reduction), ident).
#' @param cols Colors to be used for the plot, Default: NULL.
#' @param palette Color palette to be used, Default: 'glasbey'.
#' @param highlight.clusters Specific clusters to be highlighted, Default: NULL.
#' @param cells.highlight Specific cells to be highlighted, Default: NULL.
#' @param label If TRUE, clusters are labeled, Default: TRUE.
#' @param repel If TRUE, labels are repelled to avoid overlap, Default: TRUE.
#' @param legend If TRUE, a legend is added to the plot, Default: !label.
#' @param axes If TRUE, the axes are shown on the plot. Default: FALSE.
#' @param aspect.ratio Ratio of height to width. If TRUE, the ratio is fixed at 0.6. Default: TRUE.
#' @param MaxCategThrHP Maximum threshold for the number of categories, Default: 200.
#' @param save.plot If TRUE, the plot is saved into a file, Default: TRUE.
#' @param PNG If TRUE, the file is saved as a .png, Default: TRUE.
#' @param check_for_2D If TRUE, checks if UMAP is 2 dimensional. Default: TRUE.
#' @param caption Add caption to the ggplot object (e.g. a description in bottom right).
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   clUMAP("cl.names.KnownMarkers.0.5")
#'   clUMAP("cl.names.KnownMarkers.0.5", cols = NULL)
#' }
#' }
#' @importFrom MarkdownHelpers TRUE.unless
#'
#' @export
clUMAP <- function(
    ident = "integrated_snn_res.0.5", obj = combined.obj # The quickest way to draw a clustering result  UMAP
    , reduction = "umap", splitby = NULL,
    title = ident, sub = NULL,
    prefix = NULL,
    suffix = make.names(sub),
    label.cex = 7,
    h = 7, w = NULL, nr.cols = NULL,
    plotname = ppp(toupper(reduction), ident),
    cols = NULL, palette = c("alphabet", "alphabet2", "glasbey", "polychrome", "stepped")[3],
    highlight.clusters = NULL, cells.highlight = NULL,
    label = TRUE, repel = TRUE, legend = !label, MaxCategThrHP = 200,
    axes = FALSE,
    aspect.ratio = c(FALSE, 0.6)[2],
    save.plot = MarkdownHelpers::TRUE.unless("b.save.wplots"),
    PNG = TRUE,
    check_for_2D = TRUE,
    caption = FALSE
    # , save.object = FALSE
    , ...) {
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

  if (!missing(highlight.clusters)) {
    idx.ok <- identity[, 1] %in% highlight.clusters
    highlight.these <- rownames(identity)[idx.ok]
  } else {
    highlight.these <- NULL
  }
  if (!missing(cells.highlight)) {
    highlight.these <- cells.highlight
  } # overwrite, if directly defined


  if (is.null(cols)) {
    cols <- if (NtCategs > 5) getDiscretePalette(ident.used = ident, palette.used = palette, obj = obj, show.colors = FALSE)
  }
  if (!is.null(cells.highlight)) {
    cols <- "lightgrey"
  }


  if (NtCategs > MaxCategThrHP) {
    iprint("Too many categories (", NtCategs, ") in ", ident, "- use qUMAP for continous variables.")
  } else {
    if (length(unique(identity)) < MaxCategThrHP) {
      ggplot.obj <-
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

    if (!axes) ggplot.obj <- ggplot.obj + NoAxes()
    if (aspect.ratio) ggplot.obj <- ggplot.obj + ggplot2::coord_fixed(ratio = aspect.ratio)
    if (!isFALSE(caption)) ggplot.obj <- ggplot.obj + labs(caption = caption)

    if (save.plot) {
      pname <- Stringendo::sppp(prefix, plotname, suffix, sppp(highlight.clusters))
      fname <- ww.FnP_parser(pname, if (PNG) "png" else "pdf")
      try(save_plot(filename = fname, plot = ggplot.obj, base_height = h, base_width = w)) # , ncol = 1, nrow = 1
    }
    # if(save.object) saveRDS(object = ggplot.obj, file = ppp(fname, 'ggobj.RDS'))
    return(ggplot.obj)
  } # if not too many categories
}




# _________________________________________________________________________________________________
#' @title umapNamedClusters
#'
#' @description Plot and save umap based on a metadata column. #
#' @param obj Seurat object, Default: combined.obj
#' @param metaD.colname Metadata column name. Default: metaD.colname.labeled
#' @param ext File extension for saving, Default: 'png'
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   umapNamedClusters(obj = combined.obj, metaD.colname = metaD.colname.labeled)
#' }
#' }
#' @export
umapNamedClusters <- function(obj = combined.obj, metaD.colname = metaD.colname.labeled, ext = "png", ...) { # Plot and save umap based on a metadata column.
  fname <- ppp("Named.clusters", metaD.colname, ext)
  p.named <-
    Seurat::DimPlot(obj, reduction = "umap", group.by = metaD.colname, label = TRUE, ...) +
    NoLegend() +
    ggtitle(metaD.colname)
  save_plot(p.named, filename = fname)
  p.named
}




# _________________________________________________________________________________________________
#' @title umapHiLightSel
#'
#' @description Generates a UMAP plot from a Seurat object with a subset of cells highlighted.
#' @param obj A Seurat object. Default: combined.obj.
#' @param COI A vector of cluster IDs to highlight in the UMAP plot. Default: c("0", "2", "4", "5",  "11").
#' @param res.cl Name of the column in the Seurat object metadata that contains the cluster IDs. Default: 'integrated_snn_res.0.3'.
#' @return This function does not return a value. It saves a UMAP plot to the current working directory.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   umapHiLightSel(obj = seuratObject, COI = c("0", "1"), res.cl = "resolution_0.8")
#' }
#' }
#' @seealso
#' \code{\link[Seurat]{DimPlot}}
#' @export

umapHiLightSel <- function(obj = combined.obj, # Highlight a set of cells based on clusterIDs provided.
                           COI = c("0", "2", "4", "5", "11"), res.cl = "integrated_snn_res.0.3") {
  cellsSel <- getCellIDs.from.meta(obj, values = COI, ColName.meta = res.cl)
  Seurat::DimPlot(obj,
    reduction = "umap", group.by = res.cl,
    label = TRUE, cells.highlight = cellsSel
  )
  ggsave(filename = extPNG(kollapse("cells", COI, collapseby = ".")))
}



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
#' @param w Width of the plot. Default: wA4
#' @param h Height of the plot. Default: hA4
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
    list.of.genes # Save multiple FeaturePlots, as jpeg, on A4 for each gene, which are stored as a list of gene names.
    , obj = combined.obj,
    foldername = substitute(list.of.genes), plot.reduction = "umap",
    intersectionAssay = c("RNA", "integrated")[1],
    layout = c("tall", "wide", FALSE)[2],
    colors = c("grey", "red"), nr.Col = 2, nr.Row = 4, cex = round(0.1 / (nr.Col * nr.Row), digits = 2),
    gene.min.exp = "q01", gene.max.exp = "q99", subdir = TRUE,
    prefix = NULL, suffix = NULL,
    background_col = "white",
    aspect.ratio = c(FALSE, 0.6)[2],
    saveGeneList = FALSE,
    w = wA4, h = hA4, scaling = 1,
    format = c("jpg", "pdf", "png")[1],
    ...
    # , jpeg.res = 225, jpeg.q = 90
    ) {
  tictoc::tic()
  ParentDir <- OutDir
  if (is.null(foldername)) foldername <- "genes"
  if (subdir) create_set_SubDir(paste0(foldername, "-", plot.reduction), "/")
  list.of.genes.found <- check.genes(list.of.genes = list.of.genes, obj = obj, assay.slot = intersectionAssay, makeuppercase = FALSE)
  DefaultAssay(obj) <- intersectionAssay

  if (layout == "tall") {
    w <- wA4 * scaling
    h <- hA4 * scaling
    nr.Col <- 2
    nr.Row <- 4
    print("layout active, nr.Col ignored.")
  }
  if (layout == "wide") {
    w <- hA4 * scaling
    h <- wA4 * scaling
    nr.Col <- 2
    nr.Row <- 2
    print("layout active, nr.Col ignored.")
  }

  lsG <- iterBy.over(1:length(list.of.genes.found), by = nr.Row * nr.Col)
  for (i in 1:length(lsG)) {
    genes <- list.of.genes.found[lsG[[i]]]
    iprint(i, genes)
    plotname <- kpp(c(prefix, plot.reduction, i, genes, suffix, format))

    plot.list <- Seurat::FeaturePlot(
      object = obj, features = genes, reduction = plot.reduction, combine = FALSE,
      ncol = nr.Col, cols = colors,
      min.cutoff = gene.min.exp, max.cutoff = gene.max.exp,
      pt.size = cex, ...
    )

    for (i in 1:length(plot.list)) {
      plot.list[[i]] <- plot.list[[i]] + NoLegend() + NoAxes()
      if (aspect.ratio) plot.list[[i]] <- plot.list[[i]] + ggplot2::coord_fixed(ratio = aspect.ratio)
    }

    pltGrid <- cowplot::plot_grid(plotlist = plot.list, ncol = nr.Col, nrow = nr.Row)
    ggsave(filename = plotname, width = w, height = h, bg = background_col, plot = pltGrid)
  }

  if (subdir) MarkdownReports::create_set_OutDir(... = ParentDir)
  if (saveGeneList) {
    if (is.null(obj@misc$gene.lists)) obj@misc$gene.lists <- list()
    obj@misc$gene.lists[[substitute(list.of.genes)]] <- list.of.genes.found
    print("Genes saved under: obj@misc$gene.lists")
    return(obj)
  }
  tictoc::toc()
}
# _________________________________________________________________________________________________
# Save multiple FeatureHeatmaps from a list of genes on A4 jpeg
# code for quantile: https://github.com/satijalab/seurat/blob/master/R/plotting_internal.R

#' @title multiFeatureHeatmap.A4
#'
#' @description Save multiple FeatureHeatmaps from a list of genes on A4 jpeg.
#' @param obj Seurat object, Default: combined.obj
#' @param list.of.genes A list of genes to plot. No default.
#' @param gene.per.page Number of genes to plot per page. Default: 5
#' @param group.cells.by Cell grouping variable for the heatmap. Default: 'batch'
#' @param plot.reduction Dimension reduction technique to use for plots. Default: 'umap'
#' @param cex Point size in the plot. Default: iround(3/gene.per.page)
#' @param sep_scale Logical, whether to scale the features separately. Default: FALSE
#' @param gene.min.exp Minimum gene expression level for plotting. Default: 'q5'
#' @param gene.max.exp Maximum gene expression level for plotting. Default: 'q95'
#' @param jpeg.res Resolution of the jpeg output. Default: 225
#' @param jpeg.q Quality of the jpeg output. Default: 90
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @seealso
#'  \code{\link[tictoc]{tic}}
#' @importFrom tictoc tic toc
#'
#' @export
multiFeatureHeatmap.A4 <- function(
    obj = combined.obj # Save multiple FeatureHeatmaps from a list of genes on A4 jpeg
    , list.of.genes, gene.per.page = 5,
    group.cells.by = "batch", plot.reduction = "umap",
    cex = iround(3 / gene.per.page), sep_scale = FALSE,
    gene.min.exp = "q5", gene.max.exp = "q95",
    jpeg.res = 225, jpeg.q = 90, ...) {
  tictoc::tic()
  list.of.genes <- check.genes(list.of.genes, obj = obj)

  lsG <- iterBy.over(1:length(list.of.genes), by = gene.per.page)
  for (i in 1:length(lsG)) {
    print(i)
    genes <- list.of.genes[lsG[[i]]]
    plotname <- kpp(c("FeatureHeatmap", plot.reduction, i, genes, "jpg"))
    print(plotname)
    jjpegA4(plotname, r = jpeg.res, q = jpeg.q)
    try(
      FeatureHeatmap(obj,
        features.plot = genes, group.by = group.cells.by,
        reduction.use = plot.reduction, do.return = FALSE,
        sep.scale = sep_scale, min.exp = gene.min.exp, max.exp = gene.max.exp,
        pt.size = cex, key.position = "top", ...
      ),
      silent = FALSE
    )
    try.dev.off()
  }
  tictoc::toc()
}


# __________________________________________
#' @title plot.UMAP.tSNE.sidebyside
#'
#' @description Plot a UMAP and tSNE side by side.
#' @param obj Seurat object. Default: combined.obj
#' @param grouping Variable to group cells by. Default: 'res.0.6'
#' @param no_legend Logical, whether to display legend. Default: FALSE
#' @param do_return Logical, whether to return plot object. Default: TRUE
#' @param do_label Logical, whether to display labels. Default: TRUE
#' @param label_size Size of labels. Default: 10
#' @param vector_friendly Logical, whether to optimize for vector outputs. Default: TRUE
#' @param cells_use A vector of cell names to use for the plot. Default: NULL
#' @param no_axes Logical, whether to hide axes. Default: TRUE
#' @param pt_size Size of points. Default: 0.5
#' @param name.suffix Suffix to append to the plot's name. Default: NULL
#' @param width Width of the plot. Default: hA4
#' @param heigth Height of the plot. Default: 1.75 * wA4
#' @param filetype Filetype to save plot as. Default: 'pdf'
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @seealso
#'  \code{\link[cowplot]{save_plot}}
#' @importFrom cowplot save_plot plot_grid
#'
#' @export plot.UMAP.tSNE.sidebyside
plot.UMAP.tSNE.sidebyside <- function(obj = combined.obj, grouping = "res.0.6", # Plot a UMAP and tSNE sidebyside
                                      no_legend = FALSE,
                                      do_return = TRUE,
                                      do_label = TRUE,
                                      label_size = 10,
                                      vector_friendly = TRUE,
                                      cells_use = NULL,
                                      no_axes = TRUE,
                                      pt_size = 0.5,
                                      name.suffix = NULL,
                                      width = hA4, heigth = 1.75 * wA4, filetype = "pdf", ...) {
  p1 <- Seurat::DimPlot(
    object = obj, reduction.use = "tsne", no.axes = no_axes, cells.use = cells_use,
    no.legend = no_legend, do.return = do_return, do.label = do_label, label.size = label_size,
    group.by = grouping, vector.friendly = vector_friendly, pt.size = pt_size, ...
  ) +
    ggtitle("tSNE") + theme(plot.title = element_text(hjust = 0.5))

  p2 <- Seurat::DimPlot(
    object = obj, reduction.use = "umap", no.axes = no_axes, cells.use = cells_use,
    no.legend = TRUE, do.return = do_return, do.label = do_label, label.size = label_size,
    group.by = grouping, vector.friendly = vector_friendly, pt.size = pt_size, ...
  ) +
    ggtitle("UMAP") + theme(plot.title = element_text(hjust = 0.5))

  plots <- cowplot::plot_grid(p1, p2, labels = c("A", "B"), ncol = 2)
  plotname <- kpp("UMAP.tSNE", grouping, name.suffix, filetype)

  cowplot::save_plot(
    filename = plotname, plot = plots,
    ncol = 2 # we're saving a grid plot of 2 columns
    , nrow = 1 # and 2 rows
    , base_width = width,
    base_height = heigth
    # each individual subplot should have an aspect ratio of 1.3
    # , base_aspect_ratio = 1.5
  )
}

# _________________________________________________________________________________________________
#' @title PlotTopGenesPerCluster
#'
#' @description Plot the top N diff. exp. genes in each cluster.
#' @param obj Seurat object. Default: combined.obj
#' @param cl_res Resolution value to use for determining clusters. Default: res
#' @param nrGenes Number of top differentially expressed genes to plot per cluster. Default: p$'n.markers'
#' @param order.by Column name to sort output tibble by. Default: c("combined.score","avg_log2FC", "p_val_adj")[1]
#' @param df_markers Data frame containing marker gene data. Default: combined.obj@misc$df.markers[[paste0("res.", cl_res)]]
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   PlotTopGenesPerCluster(obj = combined.obj, cl_res = 0.5, nrGenes = p$"n.markers")
#' }
#' }
#' @export

PlotTopGenesPerCluster <- function(
    obj = combined.obj, cl_res = res, nrGenes = p$"n.markers",
    order.by = c("combined.score", "avg_log2FC", "p_val_adj")[1],
    df_markers = obj@misc$"df.markers"[[paste0("res.", cl_res)]]) {
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
#' @title qQC.plots.BrainOrg
#'
#' @description Quickly plot key QC markers in brain organoids
#' @param QC.Features Any numeric metadata columns
#' @param obj Seurat object, Default: combined.obj
#'
#' @examples qQC.plots.BrainOrg.RV()
#' @export

qQC.plots.BrainOrg <- function(
    obj = combined.obj, title = "Top 4 QC markers on UMAP",
    nrow = 2, ncol = 2,
    QC.Features = c("nFeature_RNA", "percent.ribo", "percent.mito", "log10.HGA_Markers")) {
  print(QC.Features)
  n.found <- setdiff(QC.Features, colnames(obj@meta.data))
  stopif(length(n.found), message = paste("n.found:", n.found))
  px <- list(
    "A" = qUMAP(QC.Features[1], save.plot = FALSE, obj = obj) + NoAxes(),
    "B" = qUMAP(QC.Features[2], save.plot = FALSE, obj = obj) + NoAxes(),
    "C" = qUMAP(QC.Features[3], save.plot = FALSE, obj = obj) + NoAxes(),
    "D" = qUMAP(QC.Features[4], save.plot = FALSE, obj = obj) + NoAxes()
  )
  qA4_grid_plot(
    plot_list = px,
    plotname = title,
    w = hA4, h = wA4,
    nrow = nrow, ncol = ncol
  )
}

# _________________________________________________________________________________________________
#' @title qMarkerCheck.BrainOrg
#'
#' @description Quickly plot key markers in brain organoids
#' @param obj Seurat object, Default: combined.obj
#' @param custom.genes Use custom gene set? Default: FALSE
#' @param suffix Folder name suffix, Default: ""
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   qMarkerCheck.BrainOrg(combined.obj)
#' }
#' }
#' @export
qMarkerCheck.BrainOrg <- function(obj = combined.obj, custom.genes = FALSE, suffix = "") {
  Signature.Genes.Top16 <- if (!isFALSE(custom.genes)) {
    custom.genes
  } else {
    Signature.Genes.Top16 <- c(
      `dl-EN` = "KAZN", `ul-EN` = "SATB2" # dl-EN = deep layer excitatory neuron
      , `Immature neurons` = "SLA",
      Interneurons = "DLX6-AS1", Interneurons = "ERBB4",
      `Intermediate progenitor` = "EOMES"
      # ,  `Intermediate progenitor1` = "TAC3"
      , `S-phase` = "TOP2A", `G2M-phase` = "HIST1H4C",
      `oRG` = "ID4", `oRG` = "HOPX" # oRG outer radial glia
      , Astroglia = "GFAP", Astrocyte = "S100B",
      `Hypoxia/Stress` = "DDIT4", Glycolytic = "PDK1",
      `Low-Quality` = "POLR2A", `Choroid.Plexus` = "DCN"
      # , `Choroid.Plexus` = "OTX2", `Choroid.Plexus` = "BMP4"
    )
    print(Signature.Genes.Top16)
  }
  print(as_tibble_from_namedVec(Signature.Genes.Top16))
  multiFeaturePlot.A4(
    obj = obj, list.of.genes = Signature.Genes.Top16, layout = "tall",
    foldername = sppp("Signature.Genes.Top16", suffix)
  )
}





# _________________________________________________________________________________________________
#' @title PlotTopGenes
#'
#' @description Plot the highest expressed genes on umaps, in a subfolder. Requires calling calc.q99.Expression.and.set.all.genes before. #
#' @param obj Seurat object, Default: combined.obj
#' @param n Number of genes to plot, Default: 32
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   PlotTopGenes()
#' }
#' }
#' @export
PlotTopGenes <- function(obj = combined.obj, n = 32, exp.slot = "expr.q99") { # Plot the highest expressed genes on umaps, in a subfolder. Requires calling calc.q99.Expression.and.set.all.genes before.
  Highest.Expressed.Genes <- names(head(sort(obj@misc[[exp.slot]], decreasing = TRUE), n = n))
  multiFeaturePlot.A4(list.of.genes = Highest.Expressed.Genes, foldername = "Highest.Expressed.Genes")
}




# _________________________________________________________________________________________________
#' @title DimPlot.ClusterNames
#'
#' @description Plot UMAP with Cluster names. #
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
DimPlot.ClusterNames <- function(obj = combined.obj # Plot UMAP with Cluster names.
                                 , ident = "cl.names.top.gene.res.0.5", reduction = "umap", title = ident, ...) {
  Seurat::DimPlot(object = obj, reduction = reduction, group.by = ident, label = TRUE, repel = TRUE, ...) + NoLegend() + ggtitle(title)
}

# _________________________________________________________________________________________________
# Manipulating UMAP and PCA  ______________________________ ----
# _________________________________________________________________________________________________



# _________________________________________________________________________________________________
#' @title FlipReductionCoordinates
#'
#' @description Flip reduction coordinates (like UMAP upside down).
#' @param obj Seurat object, Default: combined.obj
#' @param dim Numer of dimensions used, Default: 2
#' @param reduction UMAP, tSNE, or PCA (Dim. reduction to use), Default: 'umap'
#' @param flip The axis (axes) to flip around. Default: c("x", "y", "xy", NULL)[1]
#' @param FlipReductionBackupToo Flip coordinates in backup slot too? Default: TRUE
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   clUMAP()
#'   combined.obj <- FlipReductionCoordinates(combined.obj)
#'   clUMAP()
#' }
#' }
#' @export
FlipReductionCoordinates <- function(
    obj = combined.obj, dim = 2, reduction = "umap",
    flip = c("x", "y", "xy", NULL)[1], FlipReductionBackupToo = TRUE) { # Set active UMAP to `obj@reductions$umap` from `obj@misc$reductions.backup`.
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
#' @title AutoNumber.by.UMAP
#'
#' @description Relabel cluster numbers along a UMAP (or tSNE) axis #
#' @param obj Seurat object, Default: combined.obj
#' @param dim Which dimension? Default: 1
#' @param swap Swap direction? Default: FALSE
#' @param reduction UMAP, tSNE, or PCA (Dim. reduction to use), Default: 'umap'
#' @param res Clustering resoluton to use, Default: 'integrated_snn_res.0.5'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   combined.obj <- AutoNumber.by.UMAP(obj = combined.obj, dim = 1, reduction = "umap", res = "integrated_snn_res.0.5")
#'   DimPlot.ClusterNames(ident = "integrated_snn_res.0.5.ordered")
#' }
#' }
#' @export
AutoNumber.by.UMAP <- function(obj = combined.obj # Relabel cluster numbers along a UMAP (or tSNE) axis
                               , dim = 1, swap = FALSE, reduction = "umap", res = "RNA_snn_res.0.5") {
  dim_name <- kppu(toupper(reduction), dim)
  coord.umap <- as.named.vector.df(FetchData(object = obj, vars = dim_name))
  identX <- as.character(obj@meta.data[[res]])

  ls.perCl <- split(coord.umap, f = identX)
  MedianClusterCoordinate <- unlapply(ls.perCl, median)

  OldLabel <- names(sort(MedianClusterCoordinate, decreasing = swap))
  NewLabel <- as.character(0:(length(MedianClusterCoordinate) - 1))
  NewMeta <- translate(vec = identX, oldvalues = OldLabel, newvalues = NewLabel)
  NewMetaCol <- kpp(res, "ordered")
  iprint("NewMetaCol:", NewMetaCol)

  obj[[NewMetaCol]] <- NewMeta
  return(obj)
}



# _________________________________________________________________________________________________
#' @title AutoNumber.by.PrinCurve
#'
#' @description Relabel cluster numbers along the principal curve of 2 UMAP (or tSNE) dimensions. #
#' @param obj Seurat object, Default: combined.obj
#' @param dim Dimensions to use, Default: 1:2
#' @param plotit Plot results (& show it), Default: TRUE
#' @param swap Swap Lambda paramter (multiplied with this) , Default: -1
#' @param reduction UMAP, tSNE, or PCA (Dim. reduction to use), Default: 'umap'
#' @param res Clustering resoluton to use, Default: 'integrated_snn_res.0.5'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   DimPlot.ClusterNames(ident = "integrated_snn_res.0.5")
#'   combined.obj <- AutoNumber.by.PrinCurve(
#'     obj = combined.obj, dim = 1:2, reduction = "umap", plotit = TRUE,
#'     swap = -1, res = "integrated_snn_res.0.5"
#'   )
#'   DimPlot.ClusterNames(ident = "integrated_snn_res.0.5.prin.curve")
#' }
#' }
#' @seealso
#'  \code{\link[princurve]{principal_curve}}
#' @importFrom princurve principal_curve whiskers
#' @importFrom MarkdownReports wplot_save_this
#'
#' @export
AutoNumber.by.PrinCurve <- function(
    obj = combined.obj # Relabel cluster numbers along the principal curve of 2 UMAP (or tSNE) dimensions.
    , dim = 1:2, plotit = TRUE, swap = -1,
    reduction = "umap", res = "integrated_snn_res.0.5") {
  # require(princurve)
  dim_name <- ppu(toupper(reduction), dim)
  coord.umap <- FetchData(object = obj, vars = dim_name)
  fit <- princurve::principal_curve(x = as.matrix(coord.umap))
  if (plotit) {
    plot(fit,
      xlim = range(coord.umap[, 1]), ylim = range(coord.umap[, 2]),
      main = "principal_curve"
    )
    # points(fit)
    points(coord.umap, pch = 18, cex = .25)
    princurve::whiskers(coord.umap, fit$s, lwd = .1)
    MarkdownReports::wplot_save_this(plotname = "principal_curve")
  }

  ls.perCl <- split(swap * fit$lambda, f = obj[[res]])
  MedianClusterCoordinate <- unlapply(ls.perCl, median)
  OldLabel <- names(sort(MedianClusterCoordinate))
  NewLabel <- as.character(0:(length(MedianClusterCoordinate) - 1))
  NewMeta <- translate(vec = obj[[res]], oldvalues = OldLabel, newvalues = NewLabel)
  NewMetaCol <- kpp(res, "prin.curve")
  iprint("NewMetaCol:", NewMetaCol)
  obj[[NewMetaCol]] <- NewMeta
  return(obj)
}




# _________________________________________________________________________________________________
# Saving plots ______________________________ ----
# _________________________________________________________________________________________________



# _________________________________________________________________________________________________
#' @title save2umaps.A4
#'
#' @description Save 2 umaps on 1 A4
#' @param plot_list A list of plots to be saved on an A4 page.
#' @param pname Boolean to determine if name is generated automatically. If FALSE, name is generated based on plot_list and suffix. Default: FALSE
#' @param suffix A suffix added to the filename. Default: NULL
#' @param scale Scaling factor for the size of the plot. Default: 1
#' @param nrow Number of rows to arrange the plots in. Default: 2
#' @param ncol Number of columns to arrange the plots in. Default: 1
#' @param h Height of the plot, calculated as height of A4 page times the scale. Default: hA4 * scale
#' @param w Width of the plot, calculated as width of A4 page times the scale. Default: wA4 * scale
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @importFrom cowplot plot_grid
#'
#' @export
save2umaps.A4 <- function(
    plot_list, pname = FALSE, suffix = NULL, scale = 1,
    nrow = 2, ncol = 1,
    h = hA4 * scale, w = wA4 * scale, ...) { # Save 2 umaps on an A4 page.
  if (pname == FALSE) pname <- Stringendo::sppp(substitute(plot_list), suffix)
  p1 <- cowplot::plot_grid(plotlist = plot_list, nrow = nrow, ncol = ncol, labels = LETTERS[1:length(plot_list)], ...)
  save_plot(plot = p1, filename = extPNG(pname), base_height = h, base_width = w)
}

# _________________________________________________________________________________________________
#' @title save4umaps.A4
#'
#' @description Save 4 umaps on 1 A4
#' @param plot_list A list of ggplot objects, each of which is one panel.
#' @param pname Plotname, Default: FALSE
#' @param suffix A suffix added to the filename, Default: NULL
#' @param scale Scaling factor of the canvas, Default: 1
#' @param nrow number of rows for panelson the page, Default: 2
#' @param ncol number of columns for panelson the page, Default: 2
#' @param h height of the plot, Default: wA4 * scale
#' @param w width of the plot, Default: hA4 * scale
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @importFrom cowplot plot_grid
#'
#' @export
save4umaps.A4 <- function(
    plot_list, pname = FALSE, suffix = NULL, scale = 1,
    nrow = 2, ncol = 2,
    h = wA4 * scale, w = hA4 * scale, ...) { # Save 4 umaps on an A4 page.
  if (pname == FALSE) pname <- Stringendo::sppp(substitute(plot_list), suffix)
  p1 <- cowplot::plot_grid(plotlist = plot_list, nrow = nrow, ncol = ncol, labels = LETTERS[1:length(plot_list)], ...)
  save_plot(plot = p1, filename = extPNG(pname), base_height = h, base_width = w)
}



# _________________________________________________________________________________________________
#' @title qqSaveGridA4
#'
#' @description Saves a grid of 2 or 4 ggplot objects onto an A4 page.
#' @param plotlist A list of ggplot objects. Default: pl.
#' @param plots A numeric vector indicating the indices of the plots to save from the 'plotlist'. Default: 1:2.
#' @param NrPlots Number of plots to save. Default: length(plots).
#' @param height Height for the saved image. Default: hA4.
#' @param width Width for the saved image. Default: wA4.
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
    plots = 1:2, NrPlots = length(plots), height = hA4, width = wA4,
    fname = "Fractions.Organoid-to-organoid variation.png", ...) {
  stopifnot(NrPlots %in% c(2, 4))
  iprint(NrPlots, "plots found,", plots, "are saved.")
  pg.cf <- cowplot::plot_grid(plotlist = plotlist[plots], nrow = 2, ncol = NrPlots / 2, labels = LETTERS[1:NrPlots], ...)
  if (NrPlots == 4) list2env(list(height = width, width = height), envir = as.environment(environment()))
  save_plot(
    filename = fname,
    plot = pg.cf, base_height = height, base_width = width
  )
  ww.FnP_parser(fname)
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
#' @title ww.check.if.3D.reduction.exist
#'
#' @description ww.check.if.3D.reduction.exist in backup slot #
#' @param obj Seurat object, Default: obj
#' @export
ww.check.if.3D.reduction.exist <- function(obj = obj) {
  if (!("UMAP_3" %in% colnames(obj@reductions$"umap"))) {
    stopif2(
      is.null(obj@misc$reductions.backup$"umap3d"),
      "No 3D umap found in backup slot, @misc$reductions.backup. Run SetupReductionsNtoKdimensions() first."
    )
    RecallReduction(obj = obj, dim = 3, reduction = "umap")
  } else { # Reduction found in normal UMAP slot
    obj
  }
}

# _________________________________________________________________________________________________
#' @title ww.check.quantile.cutoff.and.clip.outliers
#'
#' @description Function to check a specified quantile cutoff and clip outliers from a given expression vector.
#' @param expr.vec A numeric vector representing gene expression data. Default: plotting.data[, gene]
#' @param quantileCutoffX A numeric value representing the quantile at which to clip outliers. Default: quantileCutoff
#' @param min.cells.expressing A numeric value representing the minimum number of cells expressing a gene that should remain after clipping outliers. Default: 10
#' @return The input expression vector with outliers clipped.
#' @examples
#' \dontrun{
#' ww.check.quantile.cutoff.and.clip.outliers(expr.vec = expr.data, quantileCutoffX = 0.99, min.cells.expressing = 10)
#' }
#' @importFrom CodeAndRoll2 clip.outliers.at.percentile
#'
#' @export
ww.check.quantile.cutoff.and.clip.outliers <-
  function(expr.vec = plotting.data[, gene], quantileCutoffX = quantileCutoff,
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
#' @param AutoAnnotBy The cluster or grouping to be used for automatic annotation. Default: First returned result from GetNamedClusteringRuns(obj) function.
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
#'
#' @export

plot3D.umap.gene <- function(
    gene = "TOP2A", obj = combined.obj # Plot a 3D umap with gene expression. Uses plotly. Based on github.com/Dragonmasterx87.
    , quantileCutoff = .99, def.assay = c("integrated", "RNA")[2],
    suffix = NULL, AutoAnnotBy = GetNamedClusteringRuns(obj)[1],
    alpha = .5, dotsize = 1.25, ...) {
  # stopifnot(AutoAnnotBy %in% colnames(obj@meta.data) | AutoAnnotBy = FALSE)

  obj <- ww.check.if.3D.reduction.exist(obj = obj)
  stopifnot((gene %in% rownames(obj) | gene %in% colnames(obj@meta.data)))
  DefaultAssay(object = obj) <- def.assay
  iprint(DefaultAssay(object = obj), "assay")

  plotting.data <- FetchData(object = obj, vars = c("UMAP_1", "UMAP_2", "UMAP_3", "Expression" = gene), slot = "data")

  plotting.data$"Expression" <- ww.check.quantile.cutoff.and.clip.outliers(expr.vec = plotting.data[, gene], quantileCutoffX = quantileCutoff, min.cells.expressing = 10)
  CodeAndRoll2::clip.outliers.at.percentile(plotting.data[, gene], probs = c(1 - quantileCutoff, quantileCutoff))
  plotting.data$"label" <- paste(rownames(plotting.data), " - ", plotting.data[, gene], sep = "")

  ls.ann.auto <- if (AutoAnnotBy != FALSE) {
    Annotate4Plotly3D(obj = obj, plotting.data. = plotting.data, AnnotCateg = AutoAnnotBy)
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
#' @param category The metadata column based on which the 3D UMAP will be plotted. Default: 'v.project'
#' @param obj The Seurat object for which the 3D umap plot will be generated. Default: combined.obj
#' @param suffix A suffix added to the filename. Default: NULL
#' @param AutoAnnotBy The cluster or grouping to be used for automatic annotation. Default: First returned result from GetNamedClusteringRuns(obj) function.
#' @param dotsize The size of the dots in the plot. Default: 1.25
#' @param ... Pass any other parameter to the internally called `plotly::plot_ly`.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   plot3D.umap(combined.obj, category = "Phase")
#' }
#' }
#' @importFrom plotly plot_ly layout
#'
#' @export

plot3D.umap <- function(
    category = "v.project", obj = combined.obj # Plot a 3D umap based on one of the metadata columns. Uses plotly. Based on github.com/Dragonmasterx87.
    , suffix = NULL, AutoAnnotBy = GetNamedClusteringRuns(obj)[1],
    dotsize = 1.25, ...) {
  stopifnot(category %in% colnames(obj@meta.data))
  obj <- ww.check.if.3D.reduction.exist(obj = obj)

  plotting.data <- FetchData(object = obj, vars = c("UMAP_1", "UMAP_2", "UMAP_3", category))
  colnames(plotting.data)[4] <- "category"
  plotting.data$label <- paste(rownames(plotting.data)) # Make a column of row name identities (these will be your cell/barcode names)

  ls.ann.auto <- if (AutoAnnotBy != FALSE) {
    Annotate4Plotly3D(obj = obj, plotting.data. = plotting.data, AnnotCateg = AutoAnnotBy)
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
#' @title BackupReduction
#'
#' @description Backup UMAP to `obj@misc$reductions.backup` from `obj@reductions$umap`. #
#' @param obj Seurat object, Default: combined.obj
#' @param dim Numer of dimensions used, Default: 2
#' @param reduction UMAP, tSNE, or PCA (Dim. reduction to use), Default: 'umap'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   obj <- BackupReduction(obj = obj, dim = 2, reduction = "umap")
#' }
#' }
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
#' @description Function to compute dimensionality reductions for a given Seurat object and backup the computed reductions.
#' @param obj A Seurat object. Default: combined.obj
#' @param nPCs A numeric value representing the number of principal components to use. Default: p$n.PC
#' @param dimensions A numeric vector specifying the dimensions to use for the dimensionality reductions. Default: 3:2
#' @param reduction A character string specifying the type of dimensionality reduction to perform. Can be "umap", "tsne", or "pca". Default: 'umap'
#' @return The input Seurat object with computed dimensionality reductions and backups of these reductions.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   combined.obj <- SetupReductionsNtoKdimensions(obj = combined.obj, nPCs = 10, dimensions = 2:3, reduction = "umap")
#' }
#' }
#' @export
SetupReductionsNtoKdimensions <- function(obj = combined.obj, nPCs = p$"n.PC", dimensions = 3:2, reduction = "umap", ...) { # Calculate N-to-K dimensional umaps (default = 2:3); and back them up UMAP to `obj@misc$reductions.backup` from @reductions$umap
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
#' @title RecallReduction
#'
#' @description Set active UMAP to `obj@reductions$umap` from `obj@misc$reductions.backup`. #
#' @param obj Seurat object, Default: combined.obj
#' @param dim Numer of dimensions used, Default: 2
#' @param reduction UMAP, tSNE, or PCA (Dim. reduction to use), Default: 'umap'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   combined.obj <- RecallReduction(obj = combined.obj, dim = 2, reduction = "umap")
#'   qUMAP()
#'   combined.obj <- RecallReduction(obj = combined.obj, dim = 3, reduction = "umap")
#'   qUMAP()
#' }
#' }
#' @export
RecallReduction <- function(obj = combined.obj, dim = 2, reduction = "umap") { # Set active UMAP to `obj@reductions$umap` from `obj@misc$reductions.backup`.
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
#' @title Annotate4Plotly3D
#'
#' @description Create annotation labels for 3D plots. Source https://plot.ly/r/text-and-annotations/#3d-annotations.
#' @param obj The Seurat object for which the 3D plot annotations will be generated. Default: combined.obj
#' @param plotting.data. The data frame containing plotting data. Default: plotting.data
#' @param AnnotCateg The category for which the annotation is generated. Default: AutoAnnotBy
#' @export

Annotate4Plotly3D <- function(
    obj = combined.obj # Create annotation labels for 3D plots. Source https://plot.ly/r/text-and-annotations/#3d-annotations
    , plotting.data. = plotting.data,
    AnnotCateg = AutoAnnotBy) {
  stopifnot(AnnotCateg %in% colnames(obj@meta.data))

  plotting.data.$"annot" <- FetchData(object = obj, vars = c(AnnotCateg))[, 1]
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
    plot3D.umap.gene(obj = obj., gene = g, AutoAnnotBy = annotate.by, alpha = opacity, def.assay = default.assay, dotsize = cex)
  }
  try(oo())
  try(create_set_Original_OutDir(NewOutDir = ParentDir))
}


# _________________________________________________________________________________________________
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
    plot3D.umap(obj = obj., category = categ, AutoAnnotBy = annotate.by, dotsize = cex)
  }
  try(oo())
  try(create_set_Original_OutDir(NewOutDir = ParentDir))
}
