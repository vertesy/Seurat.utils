######################################################################
# plotting.statistics.and.QC.R
######################################################################
# source ('~/GitHub/Packages/Seurat.utils/plotting.statistics.and.QC.R')
# Source: self + web

# Requirements ------------------------
require(Seurat)
require(ggplot2)
# May also require
# try (source ('/GitHub/Packages/CodeAndRoll/CodeAndRoll.R'),silent= F) # generic utilities funtions
# require('MarkdownReportsDev') # require("devtools") # plotting related utilities functions # devtools::install_github(repo = "vertesy/MarkdownReportsDev")

# sgCellFractionsBarplot.Mseq ------------------------------------------------------------------------
sgCellFractionsBarplot.Mseq <- function(data # Cell fractions Barplot for MULTI-seq. sg stands for "seurat ggplot".
                                        , seedNr=1989, group_by = "genotype", plotname="Cell proportions", downsample_to = NrCellsInSmallerDataSet) {
  set.seed(seedNr)
  data %>%
    group_by( genotype ) %>% #eval(substitute(group_by))
    sample_n(downsample_to ) %>%
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

