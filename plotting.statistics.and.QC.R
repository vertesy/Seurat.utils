######################################################################
# plotting.statistics.and.QC.R
######################################################################
# source ('~/GitHub/Seurat.utils/plotting.statistics.and.QC.R')
# Source: self + web

# Requirements ------------------------
require(Seurat)
require(ggplot2)
require(SoupX)
# May also require
# try (source ('~/GitHub/TheCorvinas/R/CodeAndRoll.R'),silent= F) # generic utilities funtions
# require('MarkdownReportsDev') # require("devtools") # plotting related utilities functions # devtools::install_github(repo = "vertesy/MarkdownReportsDev")


# PCA percent of variation associated with each PC ------------------------------------------------------------
seu.PC.var.explained <- function(obj =  combined.obj) { # Determine percent of variation associated with each PC.
  pct <- obj@reductions$pca@stdev / sum(obj@reductions$pca@stdev) * 100
  names(pct) =1:length(obj@reductions$pca@stdev)
  return(pct)
}

# plot percent of variation associated with each PC ------------------------------------------------------------
seu.plot.PC.var.explained <- function(obj =  combined.obj) { # Plot the percent of variation associated with each PC.
  pct <- seu.PC.var.explained(obj)
  wbarplot(pct, ylab= "% of variation explained" , xlab="Principal Components")
  barplot_label(round(pct, digits = 2), barplotted_variable = pct, cex=.5 )
}


# BarplotCellsPerObject ------------------------------------------------------------

BarplotCellsPerObject <- function(ls.Seu = ls.Seurat, # Take a List of Seurat objects and draw a barplot for the number of cells per object.
  plotname="Nr.Cells.After.Filtering", names=F ) {
  cellCounts = unlapply(ls.Seu, ncol)
  names(cellCounts) = if (length(names) == length(ls.Seurat)) names else names(ls.Seurat)
  wbarplot(cellCounts, plotname = plotname,tilted_text = T, ylab="Cells")
  barplot_label(cellCounts, TopOffset = 500, w = 4)
}


# sgCellFractionsBarplot.Mseq ------------------------------------------------------------------------
sgCellFractionsBarplot.Mseq <- function(data # Cell fractions Barplot for MULTI-seq. sg stands for "seurat ggplot".
  , seedNr=1989, group_by = "genotype", plotname="Cell proportions") {
  set.seed(seedNr)
  data %>%
    group_by( genotype ) %>% #eval(substitute(group_by))
    sample_n(NrCellsInSmallerDataSet ) %>%
    ssgCellFractionsBarplot.CORE(plotname = plotname)
}

# ssgCellFractionsBarplot.CORE ------------------------------------------------------------------------
ssgCellFractionsBarplot.CORE <- function(data # Cell Fractions Barplots, basic. sg stands for "seurat ggplot".
  , plotname="Cell proportions per ...", ClLabelExists = p$'clusternames.are.defined', AltLabel = p$'res.MetaD.colname') {
  LabelExists = ww.variable.exists.and.true(var = eval(ClLabelExists))

  data %>%
    aes(x=if (LabelExists) Cl.names else quote(AltLabel)) + # ggplot(aes(fill= genotype)) +
    ggtitle(plotname) +
    geom_bar( position="fill" ) +
    geom_hline(yintercept=.5, color='darkgrey')  +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    geom_text(aes(label=..count..), stat='count',position = position_fill(vjust=0.5)) +
    labs(x = "Clusters", y = "Fraction")
}

# sgCellFractionsBarplot ------------------------------------------------------------------------
sgCellFractionsBarplot <- function(data  # Cell Fractions Barplots. sg stands for "seurat ggplot".
  , seedNr=1989, group_by = "orig.ident", fill_by="experiment",label_sample_count=T, plotname="Cell proportions per ...",
                                   ClLabelExists = p$'clusternames.are.defined', AltLabel =p$'res.MetaD.colname' ) {
  LabelExists = ww.variable.exists.and.true(var = eval(ClLabelExists))
  if (LabelExists) iprint("Cl Labels found")
  set.seed(seedNr)
  data %>%
    group_by( eval(substitute(group_by)) ) %>%
    sample_n(NrCellsInSmallerDataSet ) %>%

    ggplot(aes_string(fill= fill_by)) +
    aes(x=if (LabelExists) Cl.names else quote(AltLabel)) + # ggplot(aes(fill= genotype)) +
    # ggplot(aes(fill= genotype, x=Cl.names)) + #OLD way
    geom_hline(yintercept=.5)  +
    geom_bar( position="fill" ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    { if (label_sample_count) geom_text(aes(label=..count..), stat='count', position = position_fill(vjust=0.5)) } +
    ggtitle(plotname) +
    labs(x = "Clusters", y = "Fraction")
}



# plotTheSoup ------------------------------------------------------------------------
plotTheSoup <- function(CellR.OutputDir = "~/Dropbox/Abel.IMBA/Data/SoupX_pbmc4k_demo/") { # Plot the ambient RNA content of droplets without a cell (background droplets).
  # Read In ------------------------
  sc = load10X(CellR.OutputDir, keepDroplets = TRUE)
  # Profiling the soup ------------------------
  sc = estimateSoup(sc)

  # Plot top gene's expression ----------------------------------------------------------------
  soupProfile = head(sc$soupProfile[order(sc$soupProfile$est, decreasing = TRUE), ], n = 20)
  soupX.cellfree.RNA.profile = 100 * col2named.vector(soupProfile[,1,drop=F])
  wbarplot(soupX.cellfree.RNA.profile
           , ylab="% Reads in the Soup"
           , sub = paste("Within the", basename(CellR.OutputDir), "dataset")
           , tilted_text = T)
  barplot_label(barplotted_variable = soupX.cellfree.RNA.profile
                , labels = percentage_formatter(soupX.cellfree.RNA.profile/100, digitz = 2)
                , TopOffset = .4, srt = 90, cex=.75)

  # Plot summarize expression ----------------------------------------------------------------
  soup.RP.sum   <- colSums(soupProfile[grep('^RPL|^RPS', rownames(soupProfile)),])
  soup.RPL.sum   <- colSums(soupProfile[grep('^RPL', rownames(soupProfile)),])
  soup.RPS.sum   <- colSums(soupProfile[grep('^RPS', rownames(soupProfile)),])
  soup.mito.sum <- colSums(soupProfile[grep('^MT-', rownames(soupProfile)),])

  soupProfile.summarized <- rbind(
    'Ribosomal' = soup.RP.sum,
    'Ribosomal.L' = soup.RPL.sum,
    'Ribosomal.S' = soup.RPS.sum,
    'Mitochondial' = soup.mito.sum,
    soupProfile[grep('^RPL|^RPS|^MT-', rownames(soupProfile), invert = T),]
  )

  NrColumns2Show  = min(10, nrow(soupProfile.summarized))
  soupX.cellfree.RNA.profile.summarized = 100 * col2named.vector(soupProfile.summarized[1:NrColumns2Show,1,drop=F])
  wbarplot(soupX.cellfree.RNA.profile.summarized
           , ylab="% Reads in the Soup"
           , sub = paste("Within the", basename(CellR.OutputDir), "dataset")
           , tilted_text = T)
  barplot_label(barplotted_variable = soupX.cellfree.RNA.profile.summarized
                , labels = percentage_formatter(soupX.cellfree.RNA.profile.summarized/100, digitz = 2)
                , TopOffset = .5)
  remove("sc")
  detach(SoupX)
} # plotTheSoup
