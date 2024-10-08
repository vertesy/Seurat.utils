% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat.Utils.Visualization.R
\name{scBarplot.CellFractions}
\alias{scBarplot.CellFractions}
\title{Generate Barplot of Cell Fractions}
\usage{
scBarplot.CellFractions(
  fill.by,
  group.by = GetNamedClusteringRuns()[1],
  obj = combined.obj,
  downsample = FALSE,
  min.nr.sampled.cells = 200,
  plotname = kppws("Cell proportions of", fill.by, "by", group.by),
  suffix = NULL,
  prefix = NULL,
  sub_title = suffix,
  hlines = c(0.25, 0.5, 0.75),
  return_table = FALSE,
  save_table = TRUE,
  save_plot = TRUE,
  also.pdf = FALSE,
  seedNr = 1989,
  w = NULL,
  h = 6,
  draw_plot = TRUE,
  show_numbers = FALSE,
  min.pct = 0.05,
  cex.pct = 2.5,
  min_frequency = 0,
  custom_col_palette = FALSE,
  color_scale = getDiscretePaletteObj(ident.used = fill.by, obj = obj, palette.used =
    "glasbey"),
  rnd_colors = FALSE,
  show.total.cells = TRUE,
  cex.total = 2,
  xlab.angle = 45,
  show_plot = TRUE,
  ...
)
}
\arguments{
\item{fill.by}{The variable to fill by for the bar plot.}

\item{group.by}{The variable to group by for the bar plot.}

\item{downsample}{Logical indicating whether to downsample data to equalize group sizes.}

\item{min.nr.sampled.cells}{The minimal number of cells to sample from each identity class. Defaults to 200 cells.}

\item{plotname}{The title of the plot.}

\item{suffix}{Optional suffix for the plot title.}

\item{sub_title}{Optional subtitle for the plot.}

\item{hlines}{Numeric vector specifying y-intercepts of horizontal lines to add to the plot.}

\item{return_table}{Logical; if TRUE, returns a contingency table instead of plotting.}

\item{save_plot}{Logical; if TRUE, saves the generated plot.}

\item{also.pdf}{Save plot in both png and pdf formats.}

\item{seedNr}{Seed for random number generation to ensure reproducibility.}

\item{w}{Width of the plot in inches.}

\item{h}{Height of the plot in inches.}

\item{draw_plot}{Logical; if FALSE, suppresses plotting (useful if only the table is desired).}

\item{show_numbers}{Logical; if TRUE, adds count numbers on top of each bar in the plot.}

\item{min.pct}{Show \% Labels above this threshold. Default = 0.05, or above 5 pct.}

\item{cex.pct}{Font size of pct labels.}

\item{min_frequency}{Minimum fraction to display individually in the plot; smaller fractions
are aggregated into an "Other" category.}

\item{custom_col_palette}{Specifies whether to use a standard or custom color palette.}

\item{color_scale}{Defines the color scale to use for the plot if a custom palette is selected.}

\item{show.total.cells}{Show total cells}

\item{cex.total}{Label size for total cells}

\item{...}{Additional parameters passed to internally called functions.}

\item{dsample.to.repl.thr}{Logical indicating if sampling should be done with replacement. Defaults to FALSE.}
}
\value{
Depending on the value of \code{return_table}, either returns a ggplot object or a list
containing values and percentages tables.
}
\description{
This function generates a bar plot of cell fractions per cluster from a Seurat object.
It offers the option to downsample data, equalizing the number of cells in each group
to the number in the smallest group. The plot's bars are grouped by one variable and filled by another.
The function supports custom color palettes, drawing numerical values on bars, and saving the plot.
}
\examples{
\dontrun{
if (interactive()) {
  scBarplot.CellFractions(obj = combined.obj, group.by = "integrated_snn_res.0.1", fill.by = "Phase", downsample = TRUE)
  scBarplot.CellFractions(obj = combined.obj, group.by = "integrated_snn_res.0.1", fill.by = "Phase", downsample = FALSE)
}
}
}
\seealso{
\code{\link[tools]{toTitleCase}}, \code{\link[ggplot2]{ggplot}}, \code{\link[dplyr]{group_by}}, \code{\link[dplyr]{summarise}}
}
