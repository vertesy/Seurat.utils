######################################################################
# plotting.filtering.R
######################################################################
# source ('~/GitHub/Packages/Seurat.utils/plotting.filtering.R')



# PlotFilters ------------------------------------------------------------------------------------
PlotFilters <- function(ls.obj = ls.Seurat # Plot filtering threshold and distributions, using four panels to highlight the relation between Gene- and UMI-count, ribosomal- and mitochondrial-content.
                        , subdir="Filtering.plots"
                        , parentdir= OutDirOrig
                        , suffices = names(ls.Seurat)
                        , filetype='.png'
                        , below.mito = p$"thr.lp.mito"
                        , above.mito = p$"thr.hp.mito"
                        , below.ribo = p$"thr.lp.ribo"
                        , above.ribo = p$"thr.hp.ribo"
                        , below.nFeature_RNA = p$"thr.hp.nFeature_RNA"
                        , above.nFeature_RNA = p$"thr.lp.nFeature_RNA"
                        , transparency = 0.25
                        , cex = 0.75
                        , LabelDistFromTop = 200 # for barplot_label
) {

  create_set_OutDir(parentdir, subdir)
  require(ggplot2)

  for (i in 1:l(ls.Seurat)) {
    print(samples[i])

    mm =  ls.obj[[i]]@meta.data
    filt.nFeature_RNA = (mm$nFeature_RNA < below.nFeature_RNA & mm$nFeature_RNA > above.nFeature_RNA)
    filt.below.mito = (mm$percent.mito < below.mito & mm$percent.mito > above.mito)
    # filt.below.mito = (mm$percent.mito < below.mito)
    filt.below.ribo = (mm$percent.ribo < below.ribo & mm$percent.ribo > above.ribo)
    mm =  cbind(mm, filt.nFeature_RNA, filt.below.mito, filt.below.ribo)

    mm$colour.thr.nFeature <- cut(mm$nFeature_RNA,
                                  breaks = c(-Inf, above.nFeature_RNA, below.nFeature_RNA, Inf),
                                  labels = c(p0("LQ (<", above.nFeature_RNA,")"),
                                             p0("HQ (", above.nFeature_RNA,"< X <", below.nFeature_RNA,")"),
                                             p0("Dbl/Outlier (>", below.nFeature_RNA,")")
                                  )
    )

    A = ggplot(data = mm, aes(x=nFeature_RNA, fill=colour.thr.nFeature)) +
      geom_histogram(binwidth = 100) +
      ggtitle(paste("Cells between", above.nFeature_RNA,"and",below.nFeature_RNA, " UMIs are selected (", pc_TRUE(filt.nFeature_RNA), ")")) +
      geom_vline(xintercept = below.nFeature_RNA) +
      geom_vline(xintercept = above.nFeature_RNA);
    # A

    B = ggplot(mm, aes(x=nFeature_RNA, y = percent.mito))+
      ggtitle(paste("Cells below", percentage_formatter(below.mito),
                    "mito reads are selected (with A:", pc_TRUE(filt.nFeature_RNA & filt.below.mito), ")")) +
      geom_point(alpha = transparency, size= cex,  show.legend = FALSE,
                 aes(color = filt.nFeature_RNA & filt.below.mito)  ) +
      scale_x_log10() + # scale_y_log10() +
      # annotation_logticks() +
      geom_hline(yintercept = below.mito) +
      geom_hline(yintercept = above.mito) +
      geom_vline(xintercept = below.nFeature_RNA) +
      geom_vline(xintercept = above.nFeature_RNA);
    # B


    C = ggplot(mm, aes(x=nFeature_RNA, y = percent.ribo))+
      ggtitle(paste("Cells below", percentage_formatter(below.mito),
                    "ribo reads are selected (with A:", pc_TRUE(filt.nFeature_RNA & filt.below.ribo), ")")) +
      geom_point(alpha = transparency, size= cex,  show.legend = FALSE,
                 aes(color = filt.nFeature_RNA & filt.below.ribo)  ) +
      scale_x_log10() + # scale_y_log10() +
      # annotation_logticks() +
      geom_hline(yintercept = below.ribo) +
      geom_hline(yintercept = above.ribo) +
      geom_vline(xintercept = below.nFeature_RNA) +
      geom_vline(xintercept = above.nFeature_RNA);
    # C

    D = ggplot(mm, aes(x=percent.ribo, y = percent.mito))+
      ggtitle(paste("Cells below", percentage_formatter(below.ribo),
                    "ribo reads are selected (with A,B,C:", pc_TRUE(filt.nFeature_RNA & filt.below.mito & filt.below.ribo), ")")) +
      geom_point(alpha = transparency, size= cex,  show.legend = FALSE,
                 aes(color = filt.nFeature_RNA & filt.below.mito & filt.below.ribo)  ) +
      scale_x_log10() + scale_y_log10() +
      # annotation_logticks() +
      geom_hline(yintercept = below.mito) +
      geom_hline(yintercept = above.mito) +
      geom_vline(xintercept = below.ribo) +
      geom_vline(xintercept = above.ribo);
    # D


    plot_list = list(A,B,C,D)
    px = plot_grid(plotlist = plot_list, nrow = 2, ncol = 2, labels = LETTERS[1:4])

    fname = ppp("Filtering.thresholds", suffices[i], filetype)
    save_plot(filename = fname, plot = px, base_height=12, ncol=1, nrow=1) #Figure 2

  } # for
  print(11)

  Nr.Cells.Before.Filtering = unlapply(ls.Seurat, ncol); names(Nr.Cells.Before.Filtering) = suffices
  wbarplot(Nr.Cells.Before.Filtering, tilted_text = T)
  barplot_label(Nr.Cells.Before.Filtering, TopOffset = LabelDistFromTop, w = 4)

  # End ------------------------------------------------------------------------
  create_set_Original_OutDir()

}
# PlotFilters(ls.obj = ls.Seurat)
