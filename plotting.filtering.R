######################################################################
# plotting.filtering.R
######################################################################
# source ('~/GitHub/Seurat.utils/plotting.filtering.R')


# updateHGNC helper ------------------------------------------------------------------------------------

PlotFilters <- function(obj = ls.Seurat[[i]], # Plot filtering threshold and distributions, using four panels to highlight the relation between Gene- and UMI-count, ribosomal- and mitochondrial-content.
  suffix = "org1", filetype='png' ) {
  p1 = FeatureScatter(object = obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
    geom_point(size=0.25, aes(colour =
                                nFeature_RNA > p$'thr.hp.nFeature_RNA' &
                                nFeature_RNA < p$'thr.lp.nFeature_RNA' )) +
    geom_hline(yintercept = c(p$thr.lp.nFeature_RNA, p$'thr.hp.nFeature_RNA')
               , linetype = "dashed", color = "black")

  p2 = FeatureScatter(object = obj, feature1 = "nFeature_RNA", feature2 = "percent.ribo") +
    geom_point(size=0.25, aes(colour =
                                nFeature_RNA > p$'thr.hp.nFeature_RNA' &
                                nFeature_RNA < p$'thr.lp.nFeature_RNA' &
                                percent.ribo < p$'thr.lp.ribo' )) +
    geom_vline(xintercept = c(p$'thr.lp.nFeature_RNA', p$'thr.hp.nFeature_RNA')
               , linetype = "dashed", color = "black") +
    geom_hline(yintercept = p$'thr.lp.ribo'
               , linetype = "dashed", color = "black")

  p3 = FeatureScatter(object = obj, feature1 = "nFeature_RNA", feature2 = "percent.mito") +
    geom_point(size=0.25, aes(colour =
                                nFeature_RNA > p$'thr.hp.nFeature_RNA' &
                                nFeature_RNA < p$'thr.lp.nFeature_RNA' &
                                percent.mito < p$'thr.lp.mito' )) +
    geom_vline(xintercept = c(p$'thr.lp.nFeature_RNA', p$'thr.hp.nFeature_RNA')
               , linetype = "dashed", color = "black") +
    geom_hline(yintercept = p$'thr.lp.mito'
               , linetype = "dashed", color = "black")

  p4 = FeatureScatter(object = obj, feature1 = "percent.ribo", feature2 = "percent.mito") +
    geom_point(size=0.25, aes(colour =
                                percent.ribo < p$'thr.lp.ribo' &
                                percent.mito < p$'thr.lp.mito' )) +
    geom_vline(xintercept = p$'thr.lp.ribo'
               , linetype = "dashed", color = "black") +
    geom_hline(yintercept = p$'thr.lp.mito'
               , linetype = "dashed", color = "black")

  pg = plot_grid(p1,p2,p3,p4, nrow = 2)
  ggsave(filename = ppp('Filters', suffix, filetype))
}


