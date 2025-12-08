# ____________________________________________________________________
# Seurat.utils ----
# ____________________________________________________________________
# file.edit("~/GitHub/Packages/Seurat.utils/R/Seurat.Utils.R")
# source("~/GitHub/Packages/Seurat.utils/R/Seurat.Utils.Metadata.R")
# source("~/GitHub/Packages/Seurat.utils/R/Seurat.Utils.Visualization.R")
# source("~/GitHub/Packages/Seurat.utils/R/Seurat.utils.less.used.R")

# devtools::check_man("~/GitHub/Packages/Seurat.utils")
# devtools::load_all("~/GitHub/Packages/Seurat.utils")
# devtools::document("~/GitHub/Packages/Seurat.utils"); devtools::load_all("~/GitHub/Packages/Seurat.utils")

# file.remove("~/GitHub/Packages/Seurat.utils/weight3.bar.png")



# _________________________________________________________________________________________________
# One-stop-shop functions for entire pipelines _____________________________ ------
# _________________________________________________________________________________________________

#' @title Process Seurat Objects in Parallel
#'
#' @description Applies a series of Seurat processing steps to each Seurat object in a list.
#'              The operations include scaling data, running PCA, UMAP, finding neighbors, and finding clusters.
#'              This is done in parallel using multiple cores.
#'
#' @param obj A Seurat object to be processed.
#' @param param.list A list of parameters used in the processing steps.
#' @param species_ A character string indicating the species ("human" or "mouse"). Default: "human".
#' @param update_gene_symbols A boolean indicating whether to update gene symbols from HGNC. Default: `FALSE`.
#' @param add.meta.fractions A boolean indicating whether to add metadata for fractions of cells in each cluster. Default: `FALSE`.
#' @param precompute A boolean indicating whether to compute steps: `FindVariableFeatures()`,
#' `calc.q99.Expression.and.set.all.genes()`,  `ScaleData()` and `RunPCA()` Default: `TRUE`.
#' @param compute A boolean indicating whether to compute the steps: `IntegrateLayers() / RunHarmony()`,
#' `RunUMAP()`, `FindNeighbors()`, and `FindClusters()`. Default: `TRUE`.
#' @param save A boolean indicating whether to save the results. Default: `TRUE`.
#' @param plot A boolean indicating whether to plot the results. Default: `TRUE`.
#' @param nfeatures The number of variable genes to use. Default: 2000.
#' @param variables.2.regress A list of variables to regress out. Default: NULL.
#' @param harmony.covariates A list of covariates to use for Harmony. Default: variables.2.regress.
#' @param n.PC The number of principal components to use. Default: 30.
#' @param resolutions A list of resolutions to use for clustering. Default: c(0.1, 0.2, 0.3, 0.4, 0.5).
#' @param reduction_input The reduction method to use as input for clustering & UMAP. Default: "pca".
#' @param WorkingDir The working directory to save the results. Default: getwd().
#' @param harmony.seurat.implementation A boolean indicating whether to use the Seurat implementation
#' of Harmony. Default: `FALSE`.
#' @param ... Additional parameters to be passed to `ScaleData()`.
#'
#' @return A Seurat object after applying scaling, PCA, UMAP, neighbor finding, and clustering.
#'
#' @examples
#' # Assuming ls.Seurat is a list of Seurat objects and params is a list of parameters
#' # results <- mclapply(ls.Seurat, processSeuratObject, params, mc.cores = 4)
#'
#' @details
#' Recommended to run the pipeline in sections in case some parts break.
#'
#' \preformatted{
#' # 1. Precompute
#' obj_Seu <- processSeuratObject(obj = obj_Seu, precompute = TRUE, compute = FALSE)
#'
#' # 2. Compute and save
#' obj_Seu <- processSeuratObject(obj = obj_Seu, precompute = FALSE, compute = TRUE)
#'
#' # 3. Plot
#' obj_Seu <- processSeuratObject(obj = obj_Seu, compute = TRUE, plot = TRUE)
#' }
#'
#' @importFrom Seurat ScaleData RunPCA RunUMAP FindNeighbors FindClusters
#' @importFrom tictoc tic toc
#' @importFrom harmony RunHarmony
#'
#' @export
processSeuratObject <- function(obj, param.list = p, species_ = "human",
                                update_gene_symbols = FALSE,
                                add.meta.fractions = FALSE,
                                precompute = TRUE,
                                compute = TRUE,
                                save = TRUE,
                                plot = TRUE,
                                nfeatures = param.list$"n.var.genes",
                                variables.2.regress = param.list$"variables.2.regress.combined",
                                harmony.covariates = variables.2.regress,
                                n.PC = param.list$"n.PC",
                                resolutions = param.list$"snn_res",
                                reduction_input = "pca",
                                WorkingDir = getwd(),
                                harmony.seurat.implementation = FALSE,
                                ...) {
  #
  use_harmony <- (reduction_input == "harmony")
  warning("Make sure you cleaned up the memory!", immediate. = TRUE)
  message("\nWorkingDir: ", WorkingDir)
  if (use_harmony) message("Harmony integration is attempted, but it is experimental.")
  stopifnot(require(tictoc))

  tictoc::tic("processSeuratObject")

  # Assertions to check input types _________________________________________________
  stopifnot(
    "Seurat" %in% class(obj),
    is.list(param.list),
    all(c("n.PC", "snn_res") %in% names(param.list)),
    is.numeric(n.PC), is.numeric(resolutions),
    is.character(variables.2.regress) || is.null(variables.2.regress)
  )

  if (!is.null(variables.2.regress)) {
    stopifnot("variables.2.regress is not found in @meta" = variables.2.regress %in% colnames(obj@meta.data))
  }

  iprint("nfeatures:", nfeatures)
  iprint("n.PC:", n.PC)
  iprint("snn_res:", resolutions)
  iprint("variables.2.regress (ScaleData):", variables.2.regress)
  if (use_harmony) iprint("variables.2.regress (Harmony):", harmony.covariates)

  # Save parameters _________________________________________________
  param.list$"n.var.genes" <- nfeatures
  param.list$"variables.2.regress.combined" <- variables.2.regress
  param.list$"n.PC" <- n.PC
  param.list$"snn_res" <- resolutions

  # .checkListElements(param_list = param.list, elements = c("variables.2.regress.combined", "n.PC", "snn_res"))

  obj@misc$"p" <- param.list # overwrite previous parameters
  gc()


  if (update_gene_symbols) {
    message("------------------- UpdateGenesSeurat -------------------")
    obj <- UpdateGenesSeurat(obj, ShowStats = TRUE, species_ = species_)
  }

  if (precompute | add.meta.fractions) {
    if ("data" %!in% Layers(obj)) {
      message("------------------- NormalizeData -------------------")
      obj <- Seurat::NormalizeData(object = obj)
    }
  }

  if (add.meta.fractions) {
    message("Adding metadata for gene-class fractions, e.g., percent.mito, etc.")
    obj <- addGeneClassFractions(obj, species = species_)
  } # end if add.meta.fractions


  if (precompute) {
    message("------------------- FindVariableFeatures -------------------")
    tic("FindVariableFeatures")
    obj <- FindVariableFeatures(obj,
      mean.function = "FastExpMean",
      dispersion.function = "FastLogVMR", nfeatures = nfeatures
    )
    toc()

    obj <- calc.q99.Expression.and.set.all.genes(obj = obj, quantileX = .99)

    message("------------------- ScaleData -------------------")
    tic(kpipe("ScaleData", kppc(variables.2.regress)))
    obj <- ScaleData(obj, assay = "RNA", verbose = TRUE, vars.to.regress = variables.2.regress, ...)
    toc()

    message("------------------- PCA /UMAP -------------------")
    tic("PCA")
    obj <- RunPCA(obj, npcs = n.PC, verbose = TRUE)
    toc()
  }


  if (compute) {
    if (use_harmony) {
      # Split ________________________________________________
      message("------------------- Split layers -------------------")

      m.REGR <- obj@meta.data[, harmony.covariates, drop = FALSE]
      stopif("Harmony cannot regress numeric variables" = any(sapply(m.REGR, is.numeric)))

      obj$"regress_out" <- xr <- apply(m.REGR, 1, kppu)
      nr_new_layers <- nr.unique(xr)
      cells_per_layer <- table(xr)

      warnif(
        "Too few (<5) cells in some regress_out categories:" = any(cells_per_layer < 5),
        "Too many (>25) regress_out categories" = nr_new_layers > 25
      )
      message("Number of regress_out categories:", nr_new_layers)
      message("Cells per regress_out category:", kppc(head(sort(cells_per_layer))), "...")
      if (T) hist(cells_per_layer)

      tic("Split layers by regress_out")
      obj[["RNA"]] <- split(obj[["RNA"]], f = xr, drop = FALSE)
      toc()


      message("------------------- Harmony - EXPERIMENTAL -------------------")
      if (harmony.seurat.implementation) {
        message("Using Seurat's IntegrateLayers(method = HarmonyIntegration)")
        tic("IntegrateLayers / Harmony")
        obj <- IntegrateLayers(
          object = obj, method = HarmonyIntegration, orig.reduction = "pca",
          new.reduction = "harmony", verbose = TRUE
        )
        toc()
        # No need to pass group.by.vars, it is defined by the split layers
      } else {
        tic("RunHarmony")
        # obj <- harmony::RunHarmony(object = obj, group.by.vars = "regress_out", dims.use = 1:n.PC, plot_convergence = FALSE)
        obj <- harmony::RunHarmony(object = obj, group.by.vars = harmony.covariates, dims.use = 1:n.PC, plot_convergence = FALSE)
        toc()
        # It is generally better to provide individual group.by.vars, not as a single column / string.
        # Source: https://github.com/immunogenomics/harmony/issues/246

        tic("JoinLayers")
        obj <- JoinLayers(obj, assay = "RNA")
        toc()
      }

      obj@misc$"harmony.params" <- c("n.PC" = n.PC, "regress" = harmony.covariates)
    } # end if use_harmony

    # Compute UMAP, FindNeighbors, FindClusters _________________________________________________
    message("------------------- UMAP -------------------")
    tic("UMAP")
    obj <- SetupReductionsNtoKdimensions(obj,
      nPCs = n.PC, reduction_output = "umap",
      reduction_input = reduction_input, dimensions = 3:2
    )
    toc()

    message("------------------- FindNeighbors & Clusters -------------------")
    tic("FindNeighbors")
    obj <- FindNeighbors(obj, reduction = reduction_input, dims = 1:n.PC)
    toc()

    tic("FindClusters")
    obj <- FindClusters(obj, resolution = resolutions)
    toc()
  } # END if (compute)


  if (save) {
    message("------------------- Saving -------------------")
    create_set_OutDir(WorkingDir)
    xsave(obj, suffix = "reprocessed", paramList = param.list)
  }

  if (plot) {
    message("------------------- Plotting -------------------")

    try(suPlotVariableFeatures(obj = obj, assay = "RNA"), silent = TRUE)

    try(scPlotPCAvarExplained(obj), silent = TRUE)

    try(qQC.plots.BrainOrg(obj = obj), silent = TRUE)

    # multi_clUMAP.A4(obj = obj)

    # res.ident <- paste0(DefaultAssay(obj), "_snn_res.", resolutions)[1:4]
    try(qClusteringUMAPS(obj = obj), silent = TRUE) # , idents = res.ident

    # if (ncol(obj) < 50000)
    try(qMarkerCheck.BrainOrg(obj = obj), silent = TRUE)

    Signature.Genes.Top20 <- c(
      `dl-EN` = "KAZN", `ul-EN` = "SATB2" # dl-EN = deep layer excitatory neuron
      , `Immature neurons` = "SLA", Interneurons = "DLX6-AS1",
      Interneurons = "ERBB4", Interneurons = "SCGN",
      `Intermediate progenitor` = "EOMES" # ,  `Intermediate progenitor1` = "TAC3"
      , `S-phase` = "TOP2A", `G2M-phase` = "H4C3" # formerly: HIST1H4C
      , `oRG` = "HOPX", `oRG` = "ID4" # oRG outer radial glia
      , Astroglia = "GFAP",
      Astrocyte = "S100B", `Hypoxia/Stress` = "DDIT4",
      `Choroid.Plexus` = "TTR", `Low-Quality` = "POLR2A",
      `Mesenchyme` = "DCN", Glycolytic = "PDK1",
      `Choroid.Plexus` = "OTX2", `Mesenchyme` = "DCN"
    )
    try(plotQUMAPsInAFolder(genes = Signature.Genes.Top20, obj = obj), silent = TRUE)
  }
  toc()

  return(obj)
}

# _________________________________________________________________________________________________
#' @title Run Differential Gene Expression Analysis (DGEA)
#'
#' @description Runs a differential gene expression analysis based on specified parameters,
#' reorders clusters if needed, and optionally saves results. Supports output and plotting configurations.
#'
#' @param obj Seurat object, assumed to be pre-configured with necessary data.
#' @param param.list List of parameters for DE analysis. Default: p.
#' @param res.analyzed.DE Vector of numeric values specifying the resolutions to analyze.
#'        Default: c(0.1).
#' @param ident Use this to specify a non-standard cluster identity, such as named clusters.
#' `runDGEA` will use this ident explicitly for the DE analysis. Default: NULL.
#' @param reorder.clusters Logical indicating whether to reorder clusters based on dimension.
#'        Default: `TRUE`.
#' @param reorder.dimension Integer specifying the dimension for reordering (1 for x, -1 for y).
#'        Default: 1.
#' @param add.combined.score Logical indicating whether to add a combined score to the markers.
#'        Default: `TRUE`.
#' @param save.obj Logical indicating whether to save the modified Seurat object.
#'        Default: `TRUE`.
#' @param directory Character string specifying the base directory for saving results.
#'        Default: OutDir
#' @param dir_suffix Character string specifying the suffix for the subdirectory.
#' @param subdirectory Character string specifying the subdirectory for saving outputs within
#'        the base directory. Default: "DGEA + date".
#' @param calculate.DGEA Logical determining if the DE analysis should be calculated.
#'        Default: `TRUE`.
#' @param plot.DGEA Logical determining if results should be plotted.
#'        Default: `TRUE`.
#' @param umap_caption Character string specifying the caption for the UMAP plot. Default: "".
#' @param plot.av.enrichment.hist Logical indicating whether to plot the average enrichment histogram.
#'       Default: `TRUE`.
#' @param plot.log.top.gene.stats Logical indicating whether to plot the log top gene statistics.
#' @param auto.cluster.naming Logical indicating automatic labeling of clusters.
#'        Default: `TRUE`.
#' @param clean.misc.slot Logical indicating whether to clean the misc slots of previous
#' clustering results. Default: `TRUE`.
#' @param clean.meta.data Logical indicating whether to clean the metadata slots of
#' previous clustering results. Default: `TRUE`.
#' @param n.cores Integer specifying the number of cores to use for parallel processing (multisession).
#'       Default: 1.
#' @param presto Logical indicating whether to use presto for DE analysis. Default: `TRUE`.
#' @param WorkingDir Character string specifying the working directory. Default: getwd().
#'
#' @importFrom future plan
#' @return Modified Seurat object and markers list.
#' @examples
#' runDGEA(obj = mySeuratObject, param.list = myListParams, directory = "Results/MyAnalysis")
#'
#' @export

runDGEA <- function(obj,
                    param.list = p,
                    ident = NULL,
                    res.analyzed.DE = if (is.null(ident)) c(.1) else ident, # param.list$'res.analyzed.DE'
                    reorder.clusters = if (is.null(ident)) TRUE else FALSE,
                    reorder.dimension = 1,
                    # ordering = if(any(!testNumericCompatible(res.analyzed.DE))) "no" else "ordered", # param.list$"cl.annotation"
                    # ordering = "ordered", # param.list$"cl.annotation"
                    directory,
                    dir_suffix,
                    subdirectory = ppp("DGEA_res", idate()),
                    add.combined.score = TRUE,
                    save.obj = TRUE,
                    calculate.DGEA = TRUE,
                    plot.DGEA = TRUE,
                    umap_caption = "",
                    plot.av.enrichment.hist = TRUE,
                    plot.log.top.gene.stats = TRUE,
                    auto.cluster.naming = TRUE,
                    clean.misc.slot = TRUE,
                    clean.meta.data = TRUE,
                    n.cores = 1,
                    presto = TRUE,
                    WorkingDir = getwd()) {
  if (presto) require(presto)
  message("\nWorkingDir: ", WorkingDir)

  # Assertions for input parameters
  stopifnot(
    is(obj, "Seurat"),
    is.list(param.list),
    "res.analyzed.DE should be numeric, explicit strings should be provided in: ident" =
      is.numeric(res.analyzed.DE) | !is.null(ident),
    dir.exists(directory)
  )

  create_set_OutDir(directory, subdirectory, newName = "dir_DGEA")
  dir_DGEA <- OutDir

  # Log utilized parameters from param.list
  {
    message("cl.annotation: ", if (reorder.clusters) paste("ordered:", reorder.dimension) else "no")
    message("test: ", param.list$"test")
    message("only.pos: ", param.list$"only.pos")
    message("---------------------------------")
    message("return.thresh: ", param.list$"return.thresh")
    message("logfc.threshold: ", param.list$"logfc.threshold")
    message("min.pct: ", param.list$"min.pct")
    message("min.diff.pct: ", param.list$"min.diff.pct")
    message("min.cells.group: ", param.list$"min.cells.group")
    message("max.cells.per.ident: ", param.list$"max.cells.per.ident")
  }

  # Record changes in @misc$p
  obj@misc$p$"res.analyzed.DE" <- if (is.null(ident)) res.analyzed.DE else ident
  obj@misc$p$"cl.annotation" <- if (is.null(ident)) {
    if (reorder.clusters) reorder.dimension else "no"
  } else {
    "character"
  }

  # Retrieve analyzed DE resolutions
  message("Resolutions analyzed:")
  df.markers.all <- Idents.for.DEG <- list.fromNames(x = res.analyzed.DE)


  if (clean.misc.slot) {
    message("Clearing the misc slot: df.markers and top.markers.resX")
    topMslots <- grepv("top.markers.res", names(obj@misc))
    obj@misc[topMslots] <- NULL
  }

  if (clean.meta.data) {
    message("Clearing the meta.data clustering columns.")
    topMslots <- grepv("top.markers.res", names(obj@meta.data))
    cl.ordered <- GetOrderedClusteringRuns(obj = obj)
    obj@meta.data[, cl.ordered] <- NULL
    # cl.names <- GetNamedClusteringRuns(obj = obj, pat = "^cl.names.*[0-1]\\.[0-9]",
    #                                    find.alternatives = FALSE)
    # obj@meta.data[, c(cl.ordered, cl.names)] <- NULL
  }

  # Loop through each resolution setting to find markers ________________________________________
  if (reorder.clusters) {
    message("Renumbering ----------------------------------------")
    for (i in 1:length(res.analyzed.DE)) {
      res <- res.analyzed.DE[i]
      create_set_OutDir(paste0(dir_DGEA, ppp("res", res)))
      message(i)

      # Reorder clusters based on average expression of markers
      message("Reordering clusters along dimension: ", sign(reorder.dimension), "*", if (abs(reorder.dimension) == 1) "x" else "y")
      obj <- AutoNumber.by.UMAP(
        obj = obj,
        ident = GetClusteringRuns(res = res, obj = obj)[1],
        dim = abs(reorder.dimension), reduction = "umap",
        swap = (reorder.dimension < 0), plot = TRUE
      )
    } # end for loop
  } # end if reorder.clusters

  # Set up clustering identity for DE analysis _______________________________________________
  for (i in 1:length(res.analyzed.DE)) {
    Idents.for.DEG[[i]] <-
      if (!is.null(ident)) {
        ident
      } else {
        if (reorder.clusters) {
          GetOrderedClusteringRuns(res = res.analyzed.DE[i], obj = obj)[1]
        } else {
          GetClusteringRuns(res = res.analyzed.DE[i], obj = obj)[1]
        }
      } # end if is.null(ident)
    stopifnot(Idents.for.DEG[[i]] %in% names(obj@meta.data))
  } # end for loop



  # Loop through each resolution setting to find markers ________________________________________
  if (n.cores > 1) future::plan("multisession", workers = n.cores)

  if (calculate.DGEA) {
    message("Calculating ----------------------------------------")
    for (i in 1:length(res.analyzed.DE)) {
      res <- res.analyzed.DE[i]
      tag.res <- ppp("res", res)
      df.slot <- if (!is.null(ident)) ident else tag.res

      message("Resolution: ", res, " -----------")
      create_set_OutDir(paste0(dir_DGEA, tag.res))

      message("Ident.for.DEG: ", Idents.for.DEG[[i]])
      Idents(obj) <- Idents.for.DEG[[i]]

      # Perform differential expression analysis
      tic("FindAllMarkers")
      df.markers <- Seurat::FindAllMarkers(obj,
        verbose = TRUE,
        test.use = param.list$"test",
        logfc.threshold = param.list$"logfc.threshold",
        return.thresh = param.list$"return.thresh",
        min.pct = param.list$"min.pct",
        min.diff.pct = param.list$"min.diff.pct",
        min.cells.group = param.list$"min.cells.group",
        max.cells.per.ident = param.list$"max.cells.per.ident",
        only.pos = param.list$"only.pos",
      )
      toc()


      Stringendo::stopif(is.null(df.markers))

      # order df.markers by logFC
      df.markers <- df.markers[order(df.markers$"avg_log2FC", decreasing = TRUE), ]

      if (add.combined.score) df.markers <- Add.DE.combined.score(df.markers)

      obj@misc$"df.markers"[[df.slot]] <- df.markers

      # Save results to disk
      fname <- ppp("df.markers", res)
      ReadWriter::write.simple.tsv(df.markers, filename = fname, v = FALSE)
      df.markers.all[[i]] <- df.markers
      xsave(df.markers, suffix = df.slot, v = FALSE)
    } # end for loop

    # Save final results to disk
    create_set_OutDir(directory, subdirectory)

    # Assign df.markers.all to global environment
    ReadWriter::write.simple.xlsx(named_list = df.markers.all, filename = kpp("df.markers.all", kppd(res.analyzed.DE), idate()))
    assign("df.markers.all", df.markers.all, envir = .GlobalEnv)


    if (save.obj) {
      create_set_OutDir(WorkingDir)
      tag <- if (is.null(ident)) kpp("res", res.analyzed.DE) else ident
      xsave(obj, suffix = kpp("w.DGEA", tag))
    }
  } # end if calculate.DGEA

  # Loop through each resolution setting to find markers ________________________________________
  if (plot.DGEA) {
    message("Plotting results -----------------")

    for (i in 1:length(res.analyzed.DE)) {
      res <- res.analyzed.DE[i]
      message("Resolution: ", res)
      tag.res <- ppp("res", res)
      df.slot <- if (!is.null(ident)) ident else tag.res

      create_set_OutDir(paste0(dir_DGEA, df.slot))

      df.markers <- obj@misc$"df.markers"[[df.slot]]
      Stringendo::stopif(is.null(df.markers))

      PlotTopGenesPerCluster(
        obj = obj,
        cl_res = res,
        df_markers = df.markers,
        nrGenes = param.list$"n.markers",
        order.by = param.list$"DEG.ranking"
      )

      # Automatic cluster labeling by top gene ________________________________________
      if (auto.cluster.naming) {
        message("Automatic cluster labeling by top gene.")

        obj <- StoreAllMarkers(df_markers = df.markers, res = res, obj = obj)
        obj <- AutoLabelTop.logFC(group.by = Idents.for.DEG[[i]], obj = obj, plot.top.genes = FALSE) # already plotted above

        clUMAP(ident = ppp("cl.names.top.gene", Idents.for.DEG[[i]]), obj = obj, caption = umap_caption)
      } # end if auto.cluster.naming

      # Plot per-cluster gene enrichment histogram ________________________________________
      if (plot.av.enrichment.hist) {
        message("Plotting per-cluster gene enrichment histogram.")
        # create_set_OutDir(directory, subdirectory)

        df.markers.tbl <- as_tibble(df.markers)
        df.markers.tbl$"cluster" <- as.character(df.markers.tbl$"cluster")
        p.deg.hist <- ggpubr::gghistogram(df.markers.tbl,
          x = "avg_log2FC",
          title = "Number of enriched genes per cluster",
          subtitle = "Binned by Log2(FC)",
          caption = paste(res, "| vertical line at FC of 2."),
          rug = TRUE,
          color = "cluster", fill = "cluster",
          facet.by = "cluster", xlim = c(0, 3),
          ylab = "Nr. D.E. Genes"
        ) +
          geom_vline(xintercept = 1) +
          theme_linedraw()

        qqSave(ggobj = p.deg.hist, w = 10, h = 6, title = ppp("Enrichment log2FC per cluster", res))
      }

      # Plot per-cluster enriched gene counts ________________________________________
      if (plot.log.top.gene.stats) {
        message("Plotting per-cluster enriched gene counts.")

        # Filter genes with avg_log2FC > 2
        lfc2_hiSig_genes <- df.markers |>
          dplyr::filter(avg_log2FC > 1, p_val_adj < 0.05) |>
          group_by(cluster) |>
          arrange(cluster, desc(avg_log2FC))

        top_genes <- lfc2_hiSig_genes |>
          top_n(1, avg_log2FC) |>
          pull(gene)

        # Get the number of genes per cluster
        (NrOfHighlySignLFC2_genes <- lfc2_hiSig_genes |>
          summarise(n = n()) |>
          deframe() |>
          sortbyitsnames())

        qbarplot(NrOfHighlySignLFC2_genes,
          label = NrOfHighlySignLFC2_genes,
          plotname = "Number of diff. genes per cluster",
          sub = "Genes with avg_log2FC > 1 and p_val_adj < 0.05",
          xlab = "Clusters", ylab = "Number of diff. genes"
        )

        # Write out gene lists per cluster ________________________________________
        {
          # Get the genes in a list per cluster
          genes_list <- lfc2_hiSig_genes |>
            group_by(cluster) |>
            summarise(genes = list(gene)) |>
            select(genes) |>
            deframe()

          names(genes_list) <- unique(lfc2_hiSig_genes$"cluster")
          genes_list <- sortbyitsnames(genes_list)
          names(genes_list) <- ppp("cl", names(genes_list), top_genes, "DGs")

          # write out the gene list, each element to a txt file.
          create_set_OutDir(paste0(dir_DGEA, ppp("res", res), "/top_genes"))
          for (i in 1:length(genes_list)) {
            write.simple.vec(input_vec = genes_list[[i]], filename = names(genes_list)[i], v = FALSE)
          } # for cluster
        }
      } # end if plot.log.top.gene.stats
    } # end for loop of resolutions
  } # end if plot.DGEA
  # create_set_OutDir(directory, subdirectory)

  # Return obj and df.markers.all to global environment
  return(obj)
  create_set_Original_OutDir()
} # end runDGEA



# _________________________________________________________________________________________________
# General ______________________________ ----
# _________________________________________________________________________________________________


# _________________________________________________________________________________________________
#' @title Update Seurat Object Properly, including Assays and DimReducs
#'
#' @description This function is an extension on `SeuratObject::UpdateSeuratObject()`. It
#' first calls `UpdateSeuratObject()`, to updates the class definitions of of a (v3) Seurat object,
#' then it updates its assays to the 'Assay5' class, and updates the UMAP DimReduc to keys.
#'
#' @param obj A Seurat object to be updated. Default: None.
#' @param update.gene.symbols Logical. If TRUE, gene symbols are updated to the latest version.
#'
#' @return An updated Seurat object.
#'
#' @importFrom SeuratObject UpdateSeuratObject
#' @examples
#' \dontrun{
#' combined.obj <- UpdateSeuratObjectProperly(combined.obj)
#' }
#'
#' @export
UpdateSeuratObjectProperly <- function(obj, update.gene.symbols = TRUE) {
  # Input assertions
  stopifnot(is(obj, "Seurat"))

  warning("This function is not yet fully tested. Use with caution.", immediate. = TRUE)
  message("Input obj. version: ", obj@version)

  # Update Object Structure (not Assays, etc.) _________
  if (obj@version < "5") {
    obj <- SeuratObject::UpdateSeuratObject(obj)
  } else {
    message("Object already updated to version 5. Skipping 'UpdateSeuratObject()'.")
  }

  # Update assays individually __________________
  existing_assays <- names(obj@assays)
  message("Updating assays to 'Assay5' class. Found: \n", existing_assays)
  for (assay in existing_assays) {
    obj[[assay]] <- as(obj[[assay]], Class = "Assay5")
  }

  # Update UMAP DimReduc manually __________________
  umap.exists <- !is.null(obj@reductions$umap)
  if (umap.exists) {
    message("Updating UMAP DimReduc to keys.")
    colnames(obj@reductions$umap@cell.embeddings) <- tolower(colnames(obj@reductions$umap@cell.embeddings))
    obj@reductions$umap@key <- tolower(obj@reductions$umap@key)
  } else {
    message("No UMAP DimReduc found. Skipping.")
  }

  if (update.gene.symbols) {
    message("Updating gene symbols to the latest version.")
    obj <- Seurat.utils::UpdateGenesSeurat(obj)
  }

  message("Output obj. version: ", obj@version)
  return(obj)
}


# _________________________________________________________________________________________________
#' @title parallel.computing.by.future
#'
#' @description Run gc(), load multi-session computing and extend memory limits.
#' @param cores Number of cores
#' @param maxMemSize memory limit
#'
#' @export
parallel.computing.by.future <- function(cores = 4, maxMemSize = 4000 * 1024^2) {
  # https://satijalab.org/seurat/v3.0/future_vignette.html
  cat(
    "1. If you load futures before you finished using foreach loops,
    NormalizeData inside a foreach loop fails (Error writing to connection)
    -> I assume 'future' and 'doMC' are not compatible

    2. If you setup computing on e.g. six cores, it runs 6 instances of R with the entire memory space copied.
    If you run out of memory, the system starts using the SSD as memory, and it slows you down extremely extremely extremely.
    -> Therefore it is important to clean up the memory space before setting up multicore computation.

    Loaded: library(future), workers set to 6 (def),set Max mem size to 2GB (def)."
  )

  gc(full = TRUE)
  try(memory.biggest.objects(), silent = TRUE)
  user_input <- readline(prompt = "Are you sure that memory should not be cleaned before parallelizing? (y/n)")

  if (user_input == "y") {
    iprint("N. cores", cores)
    library(future)
    # plan("multiprocess", workers = cores)
    plan("multisession", workers = cores)
    # So to set Max mem size to 2GB, you would run :
    options(future.globals.maxSize = maxMemSize)
  } else {
    print("No parallelization")
  }
}


# _________________________________________________________________________________________________
#' @title Intersect Genes with Seurat Object
#'
#' @description Intersects a set of gene names with those found in a Seurat object.
#' @param genes A vector of gene names to be intersected with the Seurat object.
#' @param obj A Seurat object containing gene expression data.
#' @param n_genes_shown Number of missing genes to be printed. Default: 10.
#' @param strict All genes to be present in the Seurat object?  Default: `TRUE`.
#' @param verbose verbose
#' @return A vector of gene names that are found both in the input 'genes' vector and the
#'         Seurat object.
#'
#' @export
IntersectGeneLsWithObject <- function(genes, obj = combined.obj, n_genes_shown = 10,
                                      species_ = "human", EnforceUnique = TRUE, ShowStats = TRUE,
                                      strict = TRUE, verbose = TRUE) {
  message(" > Running IntersectGeneLsWithObject()...")
  # "formerly IntersectWithExpressed(), which still exist in gruffi."

  stopifnot(
    is.character(genes),
    is(obj, "Seurat"),
    is.numeric(n_genes_shown) && n_genes_shown > 0,
    is.logical(strict)
  )
  stopifnot(length(genes) > 0, length(rownames(obj)) > 0)

  # Strict mode: Ensure all genes are present in the Seurat object
  all.genes.found <- all(genes %in% rownames(obj))
  if (!all.genes.found) {
    symbols.missing <- setdiff(genes, rownames(obj))
    iprint(length(symbols.missing), "symbols.missing:", symbols.missing)
    message(" > Running HGNChelper::checkGeneSymbols() to update symbols")

    HGNC.updated <- HGNChelper::checkGeneSymbols(genes, unmapped.as.na = FALSE, map = NULL, species = species_)
    if (ShowStats) {
      HGNC.updated
      print(GetUpdateStats(HGNC.updated))
    }

    if (EnforceUnique) HGNC.updated <- HGNC.EnforceUnique(HGNC.updated)
    genes <- HGNC.updated$Suggested.Symbol

    # UpdateSymbolList(symbols.missing) # Does not catch CTIP2 !!!
    if (strict) stopifnot(all(genes %in% rownames(obj)))
  }

  # Finding genes that are missing in the Seurat object
  missing_in_obj <- setdiff(genes, rownames(obj))
  if (verbose) {
    Stringendo::iprint(
      length(missing_in_obj), " (of ", length(genes),
      ") genes are MISSING from the Seurat object with (", length(rownames(obj)),
      ") genes. E.g.:", head(missing_in_obj, n_genes_shown)
    )
  }

  # Finding genes that are found in both the input list and the Seurat object
  g_found <- intersect(genes, rownames(obj))

  # Output argument assertion
  stopifnot(length(g_found) > 0)

  return(g_found)
}

# _________________________________________________________________________________________________

#' @title Intersect Genes with the List of Noticeably Expressed Genes
#'
#' @description Intersects a vector of gene names with a Seurat object to find genes that are both
#' in the input list and have expression levels in the top quantiles as defined by the object's
#' q99 expression data. It aims to filter genes based on their expression levels being above a
#' specified threshold. Additionally, it offers an option to sort the genes by their expression
#' levels in decreasing order.
#'
#' @param genes A vector of gene names to be intersected with the Seurat object.
#' @param obj A Seurat object containing gene expression data. Default: `combined.obj`.
#' @param above The expression level threshold above which genes are considered noticeably
#' expressed. Default: 0.
#' @param sort A logical flag indicating whether to sort the filtered genes by their expression
#' levels in decreasing order. Default: `FALSE`.
#' @return A vector of gene names that are found both in the input 'genes' vector and the Seurat
#' object, and have expression levels above the specified 'above' threshold. If `sort` is TRUE,
#' these genes are returned in decreasing order of their expression levels.
#'
#' @examples
#' # Assuming `genes` is a vector of gene names and `
#'
#' @export
SelectHighlyExpressedGenesq99 <- function(genes, obj = combined.obj,
                                          above = 0, sort = FALSE, strict = FALSE) {
  message(" > Running SelectHighlyExpressedGenesq99()...")
  stopifnot(is.character(genes), is(obj, "Seurat"), is.numeric(above))

  genes.expr <- IntersectGeneLsWithObject(genes = genes, obj = obj, verbose = FALSE, strict = strict)
  if (length(genes.expr) < length(genes)) message("Some genes not expressed. Recommend to IntersectGeneLsWithObject() first.")

  q99.expression <- obj@misc$expr.q99
  print(pc_TRUE(q99.expression == 0, suffix = "of genes at q99.expression are zero"))
  genes.expr.high <- q99.expression[genes.expr]
  if (sort) genes.expr.high <- sort.decreasing(genes.expr.high)
  print(genes.expr.high)
  genes.filt <- names(genes.expr.high)[genes.expr.high > above]

  SFX <- kppws("of the genes are above min. q99 expression of:", above)
  print(pc_TRUE(genes.expr %in% genes.filt, suffix = SFX))

  return(genes.filt)
}




# _________________________________________________________________________________________________
#' @title SmallestNonAboveX
#'
#' @description replace small values with the next smallest value found, which is >X.
#' @param vec Numeric input vector
#' @param X Threshold, Default: 0
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   SmallestNonZero(vec = df.markers$"p_val")
#' }
#' }
#' @export
SmallestNonAboveX <- function(vec, X = 0) {
  newmin <- min(vec[vec > X])
  vec[vec <= X] <- newmin
  vec
}


# _________________________________________________________________________________________________
#' @title AreTheseCellNamesTheSame
#'
#' @description Assert and compare two character vectors (e.g.: cell IDs) how much they overlap and
#' plot a Venn diagram. The function aborts with an error if overlap is too small.
#' @param vec1 Character vector, eg. with cell names
#' @param vec2 Character vector, eg. with cell names
#' @param names Names for plotting
#' @param min.overlap Threshold below there is no there is no meaningful overlap between the two vectors.
#'
#' @export
#' @examples # reTheseCellNamesTheSame()
AreTheseCellNamesTheSame <- function(
    vec1 = names(UVI.annot),
    vec2 = names(nr_UVI),
    names = c("Cells in Targ.Ampl", "Cells in GEX"),
    min.overlap = 0.33) {
  Cellname.Overlap <- list(vec1, vec2)
  names(Cellname.Overlap) <- if (!isFALSE(names)) names else c(substitute_deparse(vec1), substitute_deparse(vec2))

  cells.in.both <- intersect(vec1, vec2)
  sbb <- percentage_formatter(length(cells.in.both) / length(vec2), suffix = "of cells (GEX) in have a UVI assigned")
  ggExpress::qvenn(Cellname.Overlap, subt = sbb)
  iprint("Venn diagram saved.")
  iprint(sbb)

  Nr.overlapping <- length(intersect(vec1, vec2))
  Nr.total <- length(union(vec1, vec2))
  Percent_Overlapping <- Nr.overlapping / Nr.total
  print("")
  report <- percentage_formatter(Percent_Overlapping,
    prefix = "In total,",
    suffix = paste("of the cellIDs overlap across", names(Cellname.Overlap)[1], "and", names(Cellname.Overlap)[2])
  )
  print(report[1])
  stopifnot(Percent_Overlapping > min.overlap)
}


# _________________________________________________________________________________________________
#' @title Add to Misc or Tools Slot
#'
#' @description This function creates and adds a sub-slot to either the 'misc' or 'tools' slot of a
#' Seurat object. If the sub-slot already exists, it can either be overwritten or a warning will be issued.
#'
#' @param obj A Seurat object.
#' @param pocket_name Which main pocket to use: 'misc' or 'tools'. Default: 'misc'.
#' @param slot_value The value to be assigned to the sub-slot.
#' @param slot_name The name of the sub-slot. Automatically derived from 'sub_slot_value' if not provided.
#' @param sub_slot_value The value to be assigned to the sub-slot.
#' @param sub_slot_name The name of the sub-slot. Automatically derived from 'sub_slot_value' if not provided.
#' @param overwrite A boolean indicating whether to overwrite an existing sub-slot with the same name.
#'
#' @return The modified Seurat object with the new or updated sub-slot.
#'
#' @export
addToMiscOrToolsSlot <- function(obj, pocket_name = "misc",
                                 slot_value = NULL,
                                 slot_name = substitute_deparse(slot_value),
                                 sub_slot_value = NULL,
                                 sub_slot_name = substitute_deparse(sub_slot_value),
                                 overwrite = FALSE) {
  message("Running addToMiscOrToolsSlot()...")

  stopifnot(is(obj, "Seurat"),
    pocket_name %in% c("misc", "tools"),
    is.character(slot_name), length(slot_name) == 1,
    is.character(sub_slot_name), length(sub_slot_name) == 1,
    "slot name or value is provided" = is.null(slot_value) || !is.null(slot_name),
    "sub_slot name or value is provided" = is.null(sub_slot_value) || !is.null(sub_slot_name)
  )

  # Accessing the specified slot
  pocket <- slot(object = obj, name = pocket_name)

  # Creating new sub_slot or reporting if it exists
  if (slot_name %in% names(pocket) && !overwrite) {
    warning(paste(slot_name, "in", pocket_name, "already exists. Not overwritten."), immediate. = TRUE)
  } else {
    pocket[[slot_name]] <- slot_value
  }

  # Creating new sub_sub_slot or reporting if it exists
  if (sub_slot_name %in% names(pocket[[slot_name]]) && !overwrite) {
    warning(paste(sub_slot_name, "in", pocket_name, "@", slot_name, "already exists. Not overwritten."), immediate. = TRUE)
  } else {
    pocket[[slot_name]][[sub_slot_name]] <- sub_slot_value
  }


  # Assigning the modified slot back to the object
  slot(object = obj, name = pocket_name) <- pocket

  return(obj)
}

# _________________________________________________________________________________________________
#' @title Display Slots in the @tools of an Seurat Object
#'
#' @description `showToolsSlots` prints the names of slots in the `@tools` of a given object.
#' It specifically targets list elements, skipping over data frames and other non-list objects.
#'
#' @param obj An object whose `@tools` slot needs to be examined.
#' @param max.level The maximum level of nesting to print.
#' @param subslot The name of a sub-slot within the `@tools` slot to examine.
#' @param ... Additional arguments to be passed to `str`.
#'
#' @details
#' The function iterates over the slots in the `@tools` of `obj`. If a slot is a list
#' (and not a data frame), it prints the names of elements within this list. If the slot
#' is not a list or is a data frame, it skips printing the names. The function currently
#' does not use the `indent` parameter but it could be incorporated in future enhancements
#' to control the formatting of the output.
#'
#' @examples showToolsSlots(obj)
#'
#' @export
showToolsSlots <- function(obj, max.level = 1, subslot = NULL, ...) {
  slotX <- if (is.null(subslot)) obj@tools else obj@tools[[subslot]]
  str(slotX, max.level = max.level, ...)

  # tools_slot <- names(obj@tools)
  # # i=4
  # for (i in seq(tools_slot)) {
  #   cat("", fill = TRUE)
  #   message("obj@tools$", tools_slot[i])
  #
  #   x <- obj@tools[[tools_slot[i]]]
  #   if (!is.data.frame(x) & is.list(x)) {
  #     print(paste("  ", names(x)), width = 12)
  #   } else {
  #     return(idim(x))
  #   }
  # }
}


# _________________________________________________________________________________________________
#' @title Display Slots in the @misc of an Seurat Object
#'
#' @description See `showToolsSlots` for details. Prints the names of slots in the `@misc` of a given object.
#' It specifically targets list elements, skipping over data frames and other non-list objects.
#'
#' @param obj An object whose `@misc` slot needs to be examined. Default: `combined.obj`
#' @param max.level Max depth to dive into sub-elements.
#' @param subslot A subslot within `@misc`.
#' @param ... ...
#'
#' @examples showToolsSlots(obj)
#'
#' @export
showMiscSlots <- function(obj = combined.obj, max.level = 1, subslot = NULL,
                          ...) {
  slotX <- if (is.null(subslot)) obj@misc else obj@misc[[subslot]]
  str(slotX, max.level = max.level, ...)

  # Path to slot
  msg <- paste0(substitute_deparse(obj), "@misc")
  if (!is.null(subslot)) msg <- paste0(msg, "$", substitute_deparse(subslot))
  message(msg)
}




# _________________________________________________________________________________________________
#' @title calc.q99.Expression.and.set.all.genes

#' @description Calculate the gene expression of the e.g.: 99th quantile (expression in the top 1% cells).
#' @param obj Seurat object, Default: `combined.obj`
#' @param quantileX Quantile level, Default: 0.9
#' @param max.cells Max number of cells to do the calculation on. Downsample if excdeeded. Default: 1e+05
#' @param slot slot in the Seurat object. Default: 'data'
#' @param assay RNA or integrated assay, Default: c("RNA", "integrated")[1]
#' @param set.misc Create the "all.genes" variable in @misc? Default: `TRUE`.
#' @param assign_to_global_env Create the "all.genes" variable in the global env?, Default: `TRUE`.
#' @param plot Plot the expression distribution? Default: `TRUE`.
#' @param show Show the distribution plot? Default: `TRUE`.
#' @param obj.version Manuallyoverride the Version of the Seurat object. Useful when you used the
#' problematic `SeuratObject::UpdateSeuratObject()`.Default: no override `obj@version`.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   combined.obj <- calc.q99.Expression.and.set.all.genes(
#'     obj = combined.obj, quantileX = 0.9,
#'     max.cells = 25000
#'   )
#'   head(sort(as.numeric.wNames(obj@misc$expr.q90), decreasing = TRUE))
#'   combined.obj <- calc.q99.Expression.and.set.all.genes(
#'     obj = combined.obj, quantileX = 0.95,
#'     max.cells = 25000, set.all.genes = FALSE
#'   )
#' }
#' }
#' @seealso
#'  \code{\link[sparseMatrixStats]{character(0)}}
#' @importFrom tictoc tic toc
#' @importFrom sparseMatrixStats rowQuantiles
#'
#' @export

calc.q99.Expression.and.set.all.genes <- function(
    obj = combined.obj,
    quantileX = 0.99, max.cells = 1e5,
    slot = "data",
    assay = c("RNA", "integrated", "SCT")[1],
    set.misc = TRUE,
    assign_to_global_env = TRUE,
    suffix = substitute_deparse(obj),
    plot = TRUE,
    show = TRUE,
    obj.version = obj@version) {
  top.quant <- (1 - quantileX)
  message("\nCalculating the gene expression level at the the top", percentage_formatter(top.quant), " of cells. | q: ", quantileX)
  message("slot: ", slot, " assay: ", assay, ".\n")
  nr.total.cells <- ncol(obj)

  n.cells.in.top.quantile <- floor(nr.total.cells * top.quant) # number of cells in the top quantileX

  tictoc::tic("calc.q99.Expression.and.set.all.genes")
  stopifnot(
    is(obj, "Seurat"),
    quantileX > 0 & quantileX < 1,
    max.cells > 1e3, max.cells < 1e6,
    is.logical(set.misc), is.logical(assign_to_global_env), is.logical(plot), is.logical(show),
    is.character(suffix)
  )

  warnifnot(
    slot %in% c("data", "scale.data", "counts"),
    assay %in% c("RNA", "integrated"),
    ">1000 cells in the top the quantileX (with >0 expression). Increase quantileX not to miss genes expressed in small populations!" = n.cells.in.top.quantile > 1000,
    "<50 cells in the top the quantileX (with >0 expression). Decrease quantileX for robustness!" = n.cells.in.top.quantile > 50
  )

  # Get the data matrix ____________________________________________________________
  # browser()
  assay_data <- obj@assays[[assay]]
  if (obj.version >= "5") {
    if (assay == "RNA") {
      layers <- assay_data@layers
      message(length(layers), " layers in RNA assay")
      stopifnot(slot %in% names(layers))
      data_mtx <- layers[[slot]]
    } else {
      data_mtx <- slot(assay_data, slot)
    }
  } else {
    data_mtx <- slot(assay_data, slot)
  }

  # Downsample if the number of cells is too high _________________________________________________
  if (ncol(data_mtx) > max.cells) {
    dsampled <- sample(x = 1:ncol(data_mtx), size = max.cells)
    data_mtx <- data_mtx[, dsampled]
    message("Downsampled from ", ncol(obj), " to ", max.cells, " cells")
  }

  # Calculate the number of cells in the top quantile (e.g.: 99th quantile) that is
  # required to for gene expression to be >0
  message(
    "Each gene has to be expressed in min. ", n.cells.in.top.quantile, " cells, to have >0 quantile-expression\n",
    "quantileX: ", quantileX, " max.cells: ", max.cells
  )

  # Prepare for plotting ____________________________________________________________
  qname <- paste0("q", quantileX * 100)
  slot_name <- kpp("expr", qname)

  print("Calculating Gene Quantiles")
  expr.q99.df <- sparseMatrixStats::rowQuantiles(data_mtx, probs = quantileX)
  expr.q99 <- iround(expr.q99.df)

  log2.gene.expr.of.the.Xth.quantile <- as.numeric(log2(expr.q99 + 1)) # strip names
  qnameP <- paste0(100 * quantileX, "th quantile")

  # Plot the distribution of gene expression in the 99th quantile _________________________________
  if (plot) {
    pobj <- ggExpress::qhistogram(log2.gene.expr.of.the.Xth.quantile,
      plotname = paste("Gene expression in the", qnameP, " in", suffix),
      ext = "pdf", breaks = 30,
      subtitle = kollapse(pc_TRUE(expr.q99 > 0, NumberAndPC = TRUE), " genes have ", qname, " expr. > 0 (in ", nr.total.cells, " cells)."),
      caption = paste(nr.total.cells, "cells in", qnameP, "from", ncol(data_mtx), "cells in (downsampled) object."),
      suffix = suffix,
      xlab = paste0("log2(expr. in the ", qnameP, "quantile+1) [UMI]"),
      ylab = "Nr. of genes",
      plot = TRUE, save = TRUE,
      vline = .15,
      filtercol = TRUE,
      palette_use = "npg"
    )
    tictoc::toc()
    if (show) print(pobj)
  }


  # Compute a percentile rank for each genes expression level in the 99th quantile
  all.genes <- dplyr::percent_rank(expr.q99)

  # Add gene names
  genes <- rownames(obj)
  stopifnot(length(genes) == length(expr.q99))
  names(all.genes) <- genes
  all.genes <- as.list(sort(all.genes, decreasing = TRUE))

  if (assign_to_global_env) assign("all.genes", all.genes, envir = as.environment(1))

  # if (set.all.genes) obj@misc$'all.genes' = all.genes
  if (set.misc) obj@misc[[slot_name]] <- expr.q99

  iprint(
    "Quantile", quantileX, "is now stored under obj@misc$", slot_name,
    " Please execute all.genes <- obj@misc$all.genes."
  )
  return(obj)
}



# _________________________________________________________________________________________________
#' @title Filter Coding Gene Symbols (or any matching input Patterns)
#'
#' @description This function filters out gene names that match specified patterns. It reports
#' the original and final number of gene symbols and the percentage remaining after filtering.
#' It filters out non-coding gene symbols by default.
#'
#' @param genes A character vector of gene symbols.
#' @param pattern_NC A character vector of patterns to filter out non-coding gene symbols.
#' Default: c("^AC.", "^AL.", "^c[1-9]orf", "\\.AS[1-9]$", ... etc.
#' @param v "verbose" Whether to print the number of genes before and after filtering.
#' @param unique Whether to return unique gene symbols. Default: `TRUE`.
#' @param ... Additional arguments to pass to \code{\link[stringr]{str_detect}}.
#'
#' @return A character vector of filtered gene symbols.
#'
#' @examples
#' genes <- c("AC123", "AL456", "c1orf7", "TP53", "BRCA1", "X1.AS1", "MYC")
#' genes_kept <- filterCodingGenes(genes)
#' print(genes_kept)
#'
#' @importFrom stringr str_detect
#' @export
#'
filterCodingGenes <- function(
    genes, pattern_NC = c(
      "^A[CFLP][0-9]{6}", "^Z[0-9]{5}",
      "^LINC0[0-9]{4}",
      "^C[1-9]+orf[1-9]+", "^C[1-9][0-9]+orf[1-9]+", "^CXorf[1-9]+",
      "[-|\\.]AS[1-9]*$", "[-|\\.]DT[1-9]*$",
      "^MIR[1-9]", "^SNHG[1-9]",
      "^CU[0-9]{6}", "^BX[0-9]{6}",
      "^FP[0-9]{6}", "^AC[0-9]{6}"
    ),
    v = TRUE, unique = TRUE, ...) {
  # Input assertions
  stopifnot(
    is.character(genes), length(genes) > 0,
    is.character(pattern_NC), length(pattern_NC) > 0
  )

  # Filter the genes
  combined_pattern <- paste(pattern_NC, collapse = "|")
  genes_discarded <- genes[stringr::str_detect(genes, combined_pattern)]
  iprint("Example discarded", CodeAndRoll2::trail(genes_discarded))

  genes_kept <- genes[stringr::str_detect(genes, combined_pattern, negate = TRUE)]

  # Report original and final list sizes and percentage remaining
  if (v) {
    original_length <- length(genes)
    filtered_length <- length(genes_kept)
    percentage_remaining <- (filtered_length / original_length) * 100

    message("Original number of gene symbols: ", original_length)
    message("Filtered number of gene symbols: ", filtered_length)
    message("Percentage remaining: ", round(percentage_remaining, 2), "%")
  }

  # Output assertions
  stopifnot(is.character(genes_kept), length(genes_kept) <= original_length)

  if (unique) genes_kept <- unique(genes_kept)

  return(genes_kept)
}

# _________________________________________________________________________________________________
#' @title Filter and Sort Gene Expression List Based on Specified Genes and Expression Threshold
#'
#' @description This function takes a named list of gene expression values and a character vector of gene
#' symbols. It identifies the intersection of gene symbols with names in the list, filters genes based on a
#' specified expression threshold, and returns a character vector of genes that meet the criteria, sorted
#' by expression in descending order.
#'
#' @param genes Character vector of gene symbols to search for in the gene list. Default: NULL.
#' @param gene_list A named list of gene expression values where names are gene symbols, and values are
#' expression levels. Default:  all.genes
#' @param sort_by_expr Logical value specifying whether to sort the resulting gene list by expression level.
#' @param threshold Numeric value specifying the minimum expression level for filtering. Genes with
#' expression values below this threshold will be excluded. Default: 0.1.
#'
#' @return A character vector of gene symbols that match the specified list, meet the expression threshold,
#' and are sorted in descending order by expression level.
#'
#' @examples
#' # Example usage:
#' gene_list <- list(ROBO2 = 0.9982406, CDH18 = 0.9981755, DCC = 0.9981755, AL589740.1 = 0.9981103)
#' genes <- c("ROBO2", "DCC", "AL589740.1", "UNKNOWN")
#' filterExpressedGenes(gene_list, genes, threshold = 0.9981)
#'
#' @export
filterExpressedGenes <- function(
    genes, gene_list = all.genes,
    sort_by_expr = TRUE, threshold = 0.1) {
  # Assertions
  stopifnot(
    is.list(gene_list),
    is.character(genes),
    is.numeric(threshold), length(threshold) == 1
  )
  stopif(is.null(gene_list))

  # Step 1: Intersect the gene symbols with the names in the list and report statistics
  matching_genes <- CodeAndRoll2::intersect.wNames(x = genes, y = names(gene_list), names = "x")
  message(
    "Number of matching genes: ", length(matching_genes), " from ", length(genes),
    ". Missing: ", head(setdiff(genes, names(gene_list))), " ..."
  )

  # Step 2: Filter out genes below the expression threshold
  filtered_genes <- matching_genes[sapply(matching_genes, function(g) gene_list[[g]] >= threshold)]
  message("Number of genes above the threshold: ", length(filtered_genes), " from ", length(matching_genes))

  # Step 3: Conditionally sort genes according to their expression in descending order
  if (sort_by_expr) {
    order_of_expr <- order(unlist(gene_list[filtered_genes]), decreasing = TRUE)
    filtered_genes <- filtered_genes[order_of_expr]
  }

  # sorted_genes <- filtered_genes[order(sapply(filtered_genes, function(g) gene_list[[g]]), decreasing = TRUE)]
  # sorted_genes <- names(sort(unlist(gene_list[filtered_genes]), decreasing = TRUE))

  # Step 4: Return the character vector
  return(filtered_genes)
}

# r$CodeAndRoll2()
# filterExpressedGenes(AstrocyteMarkers)





# _________________________________________________________________________________________________
# Clustering ______________________________ ----
# _________________________________________________________________________________________________



# _________________________________________________________________________________________________
#' @title RenameClustering
#'
#' @description Rename clustering in a Seurat object.
#' @param namedVector named vector, where values = new, names(vec) = old
#' @param orig.ident meta.data colname original
#' @param suffix.new.ident How to name (suffix) the new identity. Default: "ManualNames"
#' @param new.ident meta.data colname new
#' @param suffix.plot Suffix description (short string) to be added to the umap plots.
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @param obj Seurat object
#'
#' @export

RenameClustering <- function(
    namedVector = ManualNames,
    orig.ident = "RNA_snn_res.0.3",
    suffix.new.ident = "ManualNames",
    new.ident = ppp(orig.ident, suffix.new.ident),
    obj = combined.obj,
    suffix.plot = "",
    plot_umaps = TRUE,
    ...) {
  NewX <- CodeAndRoll2::translate(
    vec = as.character(obj@meta.data[, orig.ident]),
    old = names(namedVector),
    new = namedVector
  )

  obj <- AddMetaData(object = obj, metadata = NewX, col.name = new.ident)

  iprint("new.ident is", new.ident, "created from", orig.ident)
  print("")

  if (plot_umaps) {
    stopifnot(is.character(suffix.plot))
    suffix.plot <- if (nchar(suffix.plot)) make.names(suffix.plot)
    print(clUMAP(orig.ident, suffix = suffix.plot, sub = suffix.plot, obj = obj, ...))
    print(clUMAP(new.ident, suffix = suffix.plot, sub = suffix.plot, obj = obj, ...))
    clUMAP(new.ident, suffix = suffix.plot, sub = suffix.plot, label = FALSE, obj = obj, ...)
  } else {
    iprint("New ident:", new.ident)
  }

  return(obj)
}


# _________________________________________________________________________________________________
#' @title Shorten Clustering Names
#'
#' @description This function takes in a string representing a clustering name,
#' and shortens it according to specific rules. It replaces "snn_res." with "",
#' "best.matching.names" with "bmatch", "ordered" with "ord",
#' "ManualNames" with "mNames", and ".long" at the end of the string with ".L".
#'
#' @param str A character string representing the clustering name to be shortened.
#'
#' @return A character string representing the shortened clustering name.
#'
#' @examples
#' \dontrun{
#' shorten_clustering_names("RNA_snn_res.0.5.ordered.ManualNames") # Returns 'RNA.0.5.ord.mNames'
#' shorten_clustering_names("RNA_snn_res.0.3.best.matching.names.ManualNames.long") # Returns 'RNA.0.3.bmatch.mNames.L'
#' shorten_clustering_names("RNA_snn_res.1.7.ordered.ManualNames.Simplest") # Returns 'RNA.1.7.ord.mNames.Simplest'
#' shorten_clustering_names("RNA_snn_res.0.5.ordered.ManualNames.Simpler") # Returns 'RNA.0.5.ord.mNames.Simpler'
#' }
#'
#' @export
shorten_clustering_names <- function(str) {
  # Replace 'snn_res' with nothing
  str <- gsub("snn_res.", "", str)
  # Replace 'best.matching.names' with 'bmatch'
  str <- gsub("best.matching.names", "bmatch", str)
  # Replace 'ordered' with 'ord'
  str <- gsub("ordered", "ord", str)
  # Replace 'ManualNames' with 'mNames'
  str <- gsub("ManualNames", "mNames", str)
  # Replace 'long' with 'L'
  str <- gsub(".long$", ".L", str)
  return(str)
}



# _________________________________________________________________________________________________
#' @title Retrieve Cluster Names
#'
#' @description Extracts cluster names based on a specified identity class from a Seurat object.
#'
#' @param obj A Seurat object. Default: `combined.obj`.
#' @param ident The identity class from which to retrieve cluster names.
#' Default uses the second clustering run from `GetClusteringRuns(obj)`.
#' @examples
#' \dontrun{
#' getClusterNames(obj = combined.obj, ident = GetClusteringRuns(obj)[2])
#' }
#' @return Prints and returns the sorted unique cluster names as a character vector.
#' @export
getClusterNames <- function(obj = combined.obj, ident = GetClusteringRuns(obj)[2]) {
  iprint("ident used:", ident)
  clz <- as.character(sort(deframe(unique(obj[[ident]]))))
  cat(dput(clz))
}



# _________________________________________________________________________________________________
#' @title GetClusteringRuns
#'
#' @description The `GetClusteringRuns` function retrieves metadata column names associated with
#'  clustering runs, based on a pattern to match, `"*snn_res.[0-9].[0-9]$"`, by default.
#' @param obj Seurat object, Default: `combined.obj`
#' @param res Clustering resolution to use, Default: `FALSE`.
#' @param pat Pattern to match, Default: `*snn_res.*[0-9]$`
#' @param v verbose, Default: `TRUE`.
#'
#' @return Prints and returns the sorted unique cluster names as a character vector.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   GetClusteringRuns(obj = combined.obj, pat = "*snn_res.*[0-9]$")
#' }
#' }
#' @export
GetClusteringRuns <- function(obj = combined.obj,
                              res = FALSE, pat = "*snn_res.[0-9].[0-9]+$",
                              v = TRUE) {
  # Get clustering results
  clustering.results <- sort(CodeAndRoll2::grepv(x = colnames(obj@meta.data), pattern = pat))

  # Check if no clustering results were found
  if (identical(clustering.results, character(0))) if (v) warning("No matching (simple) clustering column found!", immediate. = TRUE)

  if (!isFALSE(res)) {
    # Extract numeric values from clustering.results
    clustering.res.found.numeric <- as.numeric(sub(".+_snn_res.", "", clustering.results))

    # Filter clustering.results based on the numeric vector res
    clustering.results <- clustering.results[clustering.res.found.numeric %in% res]
    if (length(clustering.results) == 0) warning("No clustering matches `res`!", immediate. = TRUE)
  }

  if (v) {
    message("Clustering runs found:")
    dput(clustering.results)
  }

  return(clustering.results)
}

# _________________________________________________________________________________________________
#' @title GetNamedClusteringRuns
#'
#' @description The `GetNamedClusteringRuns` function retrieves metadata column names associated with
#'  non-numeric ("named") clustering runs, based on a pattern to match, `"Name|name"`, by default.
#' @param obj Seurat object, Default: `combined.obj`
#' @param res Clustering resolution to use, Default: c(FALSE, 0.5)[1]
#' @param topgene Match clustering named after top expressed gene (see vertesy/Seurat.pipeline/~Diff gene expr.), Default: `FALSE`.
#' @param pat Pattern to match, Default: '^cl.names.Known.*[0,1]\.[0-9]$'
#' @param find.alternatives If TRUE, tries to find alternative clustering runs with
#' the same resolution, Default: `TRUE`.
#' @param v Verbose output, Default: `TRUE`.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   GetNamedClusteringRuns()
#' }
#' }
#' @export
GetNamedClusteringRuns <- function(
    obj = combined.obj,
    res = list(FALSE, 0.5)[[1]], topgene = FALSE,
    pat = c("^cl.names.top.gene.+[0-9]\\.[0-9]", "Name|name")[2],
    find.alternatives = TRUE,
    v = TRUE) {
  if (res) pat <- gsub(x = pat, pattern = "\\[.*\\]", replacement = res)
  if (topgene) pat <- gsub(x = pat, pattern = "Known", replacement = "top")
  clustering.results <- CodeAndRoll2::grepv(x = colnames(obj@meta.data), pattern = pat)

  if (identical(clustering.results, character(0))) {
    if (v) warning("No matching (named) clustering column found! Trying GetClusteringRuns(..., pat = '*_res.*[0,1]\\.[0-9]$)", immediate. = TRUE)
    if (find.alternatives) {
      clustering.results <-
        GetClusteringRuns(obj = obj, res = FALSE, pat = "*_res.*[0,1]\\.[0-9]$", v = FALSE)
    }
  }

  if (v) dput(clustering.results)

  return(clustering.results)
}



# _________________________________________________________________________________________________
#' @title GetOrderedClusteringRuns
#'
#' @description Get Clustering Runs: metadata column names.
#' @param obj Seurat object, Default: `combined.obj`.
#' @param res Clustering resolution to use, Default: `FALSE`.
#' @param pat Pattern to match, Default: '*snn_res.*[0,1]\.[0-9]\.ordered$'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   GetOrderedClusteringRuns()
#'   GetOrderedClusteringRuns(res = 0.5)
#' }
#' }
#' @export
GetOrderedClusteringRuns <- function(obj = combined.obj, res = FALSE,
                                     pat = "*snn_res.*[0,1]\\.[0-9]\\.ordered$") {
  if (res) pat <- gsub(x = pat, pattern = "\\[.*\\]", replacement = res)
  clustering.results <- CodeAndRoll2::grepv(x = colnames(obj@meta.data), pattern = pat)
  if (identical(clustering.results, character(0))) warning("No matching (ordered) clustering column found!", immediate. = TRUE)
  return(clustering.results)
}



# _________________________________________________________________________________________________
#' @title GetNumberOfClusters
#'
#' @description Print the number of clusters for each stored clustering run.
#' @param obj Seurat object, Default: `combined.obj`
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   GetNumberOfClusters()
#' }
#' }
#' @export
GetNumberOfClusters <- function(obj = combined.obj) { # Get Number Of Clusters
  clustering.results <- GetClusteringRuns(obj)
  print("## Number of clusters: ---------")
  for (cc in clustering.results) {
    NrCl <- length(unique(obj@meta.data[[cc]]))
    iprint(cc, "   ", NrCl)
  }
}


# _________________________________________________________________________________________________
#' @title calc.cluster.averages
#'
#' @description Calculates the average of a metadata column (numeric) per cluster.
#' @param col_name The name of the column for which the average is calculated. Default: 'Score.GO.0006096'.
#' @param plot.UMAP.too Whether to plot a UMAP as well. Default: `TRUE`.
#' @param return.plot Whether to return the plot. Default: `FALSE`.
#' @param obj The main Seurat object used for calculations. Default: `combined.obj`.
#' @param split_by Cluster to split by. Default: First entry of GetNamedClusteringRuns().
#' @param scale.zscore Whether to scale z-scores. Default: `FALSE`.
#' @param simplify Whether to simplify the result. Default: `TRUE`.
#' @param plotit Whether to plot the results. Default: `TRUE`.
#' @param histogram Whether to produce a histogram. Default: `FALSE`.
#' @param nbins The number of bins for the histogram. Default: 50.
#' @param suffix Suffix added to the filename. Default: NULL.
#' @param stat Statistical method applied, "mean" or "median". Default: "median".
#' @param quantile.thr The threshold for quantiles. Default: 0.9.
#' @param absolute.thr Absolute threshold used in computations. Default: `FALSE`.
#' @param filter The filter mode: 'above', 'below', or FALSE. Default: `FALSE`.
#' @param ylab.text Text for the y-axis label. Default: "Cluster" followed by the statistical method and "score".
#' @param title Title for the plot. Default: "Cluster" followed by the statistical method and column name.
#' @param subtitle The subtitle for the plot. Default: NULL.
#' @param width The width of the plot. Default: 8.
#' @param height The height of the plot. Default: 6.
#' @param ... Additional parameters passed to the internally called functions.
#' @param xlb The label for the x-axis. Default depends on the 'absolute.thr' parameter.
#' @param fname The filename for the plot. Default: based on column name and split_by value.
#' @export
#' @importFrom Stringendo percentage_formatter

calc.cluster.averages <- function(
    col_name = "Score.GO.0006096",
    plot.UMAP.too = TRUE,
    return.plot = FALSE,
    obj = combined.obj,
    split_by = GetNamedClusteringRuns()[1],
    scale.zscore = FALSE,
    simplify = TRUE, plotit = TRUE,
    histogram = FALSE, nbins = 50,
    suffix = NULL,
    stat = c("mean", "median")[2],
    quantile.thr = 0.9,
    absolute.thr = FALSE,
    filter = c(FALSE, "above", "below")[1],
    ylab.text = paste("Cluster", stat, "score"),
    title = paste("Cluster", stat, col_name),
    prefix.cl.names = FALSE,
    report = TRUE,
    subtitle = NULL,
    width = 8, height = 6,
    ...
    # , ylb = paste(ylab.text, col_name)
    # , xlb = paste("Clusters >",Stringendo::percentage_formatter(quantile.thr),"quantile are highlighted. |", split_by)
    , xlb = if (absolute.thr) {
      paste("Threshold at", absolute.thr)
    } else {
      paste(
        "Black lines: ", kppd(Stringendo::percentage_formatter(c(1 - quantile.thr, quantile.thr))), "quantiles |",
        "Cl. >", Stringendo::percentage_formatter(quantile.thr), "are highlighted. |", split_by
      )
    },
    fname = ppp(col_name, split_by, "cluster.average.barplot.pdf", ...)) { # calc.cluster.averages of a m
  iprint(substitute_deparse(obj), "split by", split_by)
  if (absolute.thr) iprint("In case of the absolute threshold, only the returned values are correct, the plot annotations are not!")

  if (plot.UMAP.too) qUMAP(obj = obj, feature = col_name)

  df.summary <-
    obj@meta.data |>
    select_at(c(col_name, split_by)) |>
    group_by_at(split_by) |>
    summarize(
      "nr.cells" = n(),
      "mean" = mean(!!sym(col_name), na.rm = TRUE),
      "SEM" = sem(!!sym(col_name), na.rm = TRUE),
      "median" = median(!!sym(col_name), na.rm = TRUE),
      "SE.median" = 1.2533 * sem(!!sym(col_name), na.rm = TRUE)
    )

  if (simplify) {
    av.score <- df.summary[[stat]]
    names(av.score) <- if (!isFALSE(prefix.cl.names)) ppp("cl", df.summary[[1]]) else df.summary[[1]]
    av.score <- sortbyitsnames(av.score)
    if (scale.zscore) av.score <- (scale(av.score)[, 1])

    cutoff <- if (absolute.thr) absolute.thr else quantile(av.score, quantile.thr)
    cutoff.low <- if (absolute.thr) NULL else quantile(av.score, (1 - quantile.thr))

    iprint("quantile.thr:", quantile.thr)
    if (plotit) {
      if (histogram) {
        p <- ggExpress::qhistogram(
          vec = as.numeric(av.score), save = FALSE,
          vline = cutoff,
          plotname = ppp(title, quantile.thr),
          bins = nbins,
          subtitle = paste(subtitle, "| median in blue/dashed"),
          ylab = ylab.text,
          xlab = xlb # Abused
          , xlab.angle = 45
          # , ylim = c(-1,1)
          , ...
          # , ext = "png", w = 7, h = 5
        ) + geom_vline(xintercept = cutoff.low, lty = 2)
        print(p)
        title_ <- ppp(title, suffix, flag.nameiftrue(scale.zscore))
        ggExpress::qqSave(ggobj = p, title = title_, ext = "png", w = width, h = height)
      } else {
        p <- ggExpress::qbarplot(
          vec = av.score, save = FALSE,
          hline = cutoff,
          plotname = title,
          suffix = quantile.thr,
          subtitle = subtitle,
          ylab = ylab.text,
          xlab = xlb # Abused
          , xlab.angle = 45
          # , ylim = c(-1,1)
          , ...
          # , ext = "png", w = 7, h = 5
        ) + geom_hline(yintercept = cutoff.low, lty = 2)

        print(p)
        title_ <- ppp(title, suffix, flag.nameiftrue(scale.zscore))
        qqSave(ggobj = p, title = title_, fname = ppp(title_, split_by, "png"), w = width, h = height)
      }
    }

    if (report) print(paste0(col_name, ": ", paste(iround(av.score), collapse = " vs. ")))
    if (filter == "below") {
      return(filter_LP(av.score, threshold = cutoff, plot.hist = FALSE))
    } else if (filter == "above") {
      return(filter_HP(av.score, threshold = cutoff, plot.hist = FALSE))
    } else {
      return(av.score)
    }
  } else if (return.plot) { # if /not/ simplify
    return(p)
  } else {
    return(df.summary)
  }
}






# _________________________________________________________________________________________________
#' @title plot.expression.rank.q90
#'
#' @description Plot gene expression based on the expression at the 90th quantile
#' (so you will not lose genes expressed in few cells).
#' @param obj Seurat object, Default: `combined.obj`
#' @param gene gene of interest, Default: 'ACTB'
#' @param filterZero Remove genes whose quantile-90 expression in 0? Default: `TRUE`.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   plot.expression.rank.q90(gene = "SATB2")
#' }
#' }
#' @importFrom Stringendo percentage_formatter
#' @importFrom MarkdownReports whist
#'
#' @export plot.expression.rank.q90
plot.expression.rank.q90 <- function(obj = combined.obj, gene = "ACTB", filterZero = TRUE) {
  expr.GOI <- obj@misc$expr.q90[gene]
  expr.all <- unlist(obj@misc$expr.q90)
  gene.found <- gene %in% names(expr.all)
  stopifnot(gene.found)

  if (expr.GOI == 0) iprint(gene, "is not expressed. q90-av.exp:", expr.GOI) else if (expr.GOI < 0.05) iprint(gene, "is lowly expressed. q90-av.exp:", expr.GOI)
  if (filterZero) {
    iprint("Zero 'q90 expression' genes (", pc_TRUE(expr.all == 0), ") are removed.")
    expr.all <- expr.all[expr.all > 0]
  }
  counts <- sum(obj@assays$RNA@counts[gene, ])
  if (expr.GOI == 0) {
    quantile.GOI <- 0
    title <- paste(gene, "is too lowly expressed: q90-av.exp is zero. \n There are", counts, "counts.")
  } else {
    pos.GOI <- which(names(expr.all) == gene)
    quantile.GOI <- ecdf(expr.all)(expr.all)[pos.GOI]
    title <- paste(
      gene, "is in the", Stringendo::percentage_formatter(quantile.GOI),
      "quantile of 'q90-av' expression. \n There are", counts, "counts"
    )
  }
  suppressWarnings(
    MarkdownReports::whist(expr.all,
      vline = expr.GOI, breaks = 100, main = title, plotname = make.names(title),
      ylab = "Genes",
      xlab = "Av. mRNA in the 10% top expressing cells (q90 av.exp.)"
    )
  )
}




# _________________________________________________________________________________________________
# Interacting with the environment ______________________________ ----
# _________________________________________________________________________________________________
# Subsetting, downsampling and manipulating the Seurat object



# _________________________________________________________________________________________________
#' @title set.mm
#'
#' @description Helps to find metadata columns. It creates a list with the names of of 'obj@meta.data'.
#' @param obj Seurat object, Default: `combined.obj`
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   set.mm()
#'   mm
#' }
#' }
#' @export
set.mm <- function(obj = combined.obj) {
  mm <- CodeAndRoll2::list.fromNames(colnames(obj@meta.data))
  assign(x = "mm", value = mm, envir = as.environment(1))
}


# _________________________________________________________________________________________________
#' @title Get the First Seurat Object from a List of Seurat Objects
#'
#' @description Return the first Seurat object from a list, or the object itself
#'   if a single Seurat object is supplied.
#'
#' @param obj A Seurat object, a list of Seurat objects, or any other list.
#'
#' @return The first Seurat object from the list or the Seurat object itself.
#' If the input is not a Seurat object or a list containing at least one Seurat
#' object, the function will throw an error.
#' @export
ww.get.1st.Seur.element <- function(obj) {
  if (is(obj)[1] == "list") {
    iprint("A list of objects is provided, taking the 1st from", length(obj), "elements.")
    obj <- obj[[1]]
  }
  stopifnot(is(obj) == "Seurat")
  return(obj)
}
# ww.get.1st.Seur.element(ls.Seurat[[1]])


# _________________________________________________________________________________________________
#' @title Recall all.genes global variable from a Seurat object
#'
#' @description Recall \code{all.genes} from a Seurat object's \code{misc} slot,
#'   which is stored by \code{calc.q99.Expression.and.set.all.genes()}, and optionally reset the global variable.
#' @param obj Seurat object, Default: `combined.obj`
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   recall.all.genes()
#'   all.genes
#' }
#' }
#' @importFrom MarkdownHelpers ww.assign_to_global
#'
#' @export
recall.all.genes <- function(obj = combined.obj, overwrite = FALSE) {
  obj <- ww.get.1st.Seur.element(obj)

  if ("all.genes" %in% names(obj@misc)) {
    if (!exists("all.genes") | overwrite) {
      all.genes <- obj@misc$all.genes
      print(head(unlist(all.genes)))
      MarkdownHelpers::ww.assign_to_global(name = "all.genes", value = all.genes, verbose = FALSE)
      message("all.genes is now (re)defined in the global environment.")
    } else {
      message("  ->   Variable 'all.genes' exists in the global namespace, and overwrite is: FALSE")
    }
  } else {
    message("  ->   Slot 'all.genes' does not exist in obj@misc.")
    hits <- grepv(pattern = "expr.", names(obj@misc))
    if (!is.null(hits)) {
      message("Found instead (", hits, "). Returning 1st element: ", hits[1])
      all.genes <- obj@misc[[hits[1]]]
      MarkdownHelpers::ww.assign_to_global(name = "all.genes", value = as.list(all.genes), verbose = FALSE)
    }
  }
}


# _________________________________________________________________________________________________
#' @title recall.meta.tags.n.datasets
#'
#' @description Recall  meta.tags from obj@misc to "meta.tags" in the global environment.
#' @param obj Seurat object, Default: `combined.obj`
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   recall.n.datasets()
#'   n.datasets
#' }
#' }
#' @importFrom MarkdownHelpers ww.assign_to_global
#'
#' @export
recall.meta.tags.n.datasets <- function(obj = combined.obj) {
  obj <- ww.get.1st.Seur.element(obj)

  if ("n.datasets" %in% names(obj@misc)) {
    if (!exists("n.datasets")) {
      n.datasets <- obj@misc$n.datasets
      print(head(unlist(n.datasets)))
      MarkdownHelpers::ww.assign_to_global(name = "n.datasets", value = n.datasets)
    } else {
      print("  ->   Variable 'n.datasets' already exists in the global namespace.")
    }
  } else {
    print("  ->   Slot 'n.datasets' does not exist in obj@misc.")
  }


  if ("meta.tags" %in% names(obj@misc)) {
    if (!exists("meta.tags")) {
      meta.tags <- obj@misc$meta.tags
      print(head(unlist(meta.tags)))
      MarkdownHelpers::ww.assign_to_global(name = "meta.tags", value = meta.tags)
    } else {
      iprint("  ->   Variable 'meta.tags' already exists in the global namespace.")
    }
  } else {
    print("  ->   Slot 'meta.tags' does not exist in obj@misc.")
  }
}


# _________________________________________________________________________________________________
#' @title recall.parameters
#'
#' @description Recall parameters from obj@misc to "p" in the global environment.
#' @param obj Seurat object, Default: `combined.obj`
#' @param overwrite Overwrite already existing in environment? Default: `FALSE`.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   recall.parameters()
#'   p
#' }
#' }
#' @importFrom MarkdownHelpers ww.assign_to_global
#'
#' @export
recall.parameters <- function(obj = combined.obj, overwrite = FALSE) {
  obj <- ww.get.1st.Seur.element(obj)

  if ("p" %in% names(obj@misc)) {
    p_found <- exists("p", envir = .GlobalEnv)
    if (p_found) message("  ->   Variable 'p' exists in the global namespace.")

    if (!p_found | (p_found & overwrite == TRUE)) {
      MarkdownHelpers::ww.assign_to_global(name = "p", value = obj@misc$"p", verbose = FALSE)
      message("p is now (re)defined in the global environment.")
    } else {
      message("p not overwritten.")
    }
  } else {
    message("  ->   Slot 'p' does not exist in obj@misc.")
  }
}



# _________________________________________________________________________________________________
#' @title recall.genes.ls
#'
#' @description Recall genes.ls from obj@misc to "genes.ls" in the global environment.
#' @param obj Seurat object, Default: `combined.obj`
#' @param overwrite Overwrite already existing in environment? Default: `FALSE`.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   recall.genes.ls()
#'   genes.ls
#' }
#' }
#' @importFrom MarkdownHelpers ww.assign_to_global
#'
#' @export

recall.genes.ls <- function(obj = combined.obj, overwrite = FALSE) { # genes.ls
  obj <- ww.get.1st.Seur.element(obj)

  if ("genes.ls" %in% names(obj@misc)) {
    if (!exists("genes.ls")) message("variable 'genes.ls' exists in the global namespace: ", head(p))

    if (!exists("genes.ls") | (exists("genes.ls") & overwrite == TRUE)) {
      MarkdownHelpers::ww.assign_to_global(name = "genes.ls", value = obj@misc$"genes.ls")
      message("Overwritten.")
    } else {
      message("Not overwritten.")
    }
  } else {
    message("  ->   Slot 'genes.ls' does not exist in obj@misc.")
  }
}


# _________________________________________________________________________________________________
#' @title Save Parameters to Seurat Object
#'
#' @description Stores a list of parameters within the `@misc$p` slot of a Seurat object,
#' allowing for easy reference and tracking of analysis parameters used.
#'
#' @param obj Seurat object to update. Default: `combined.obj`.
#' @param params List of parameters to save. Default: `p`.
#' @param overwrite Logical indicating if existing parameters should be overwritten. Default: `TRUE`.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   save.parameters(obj = combined.obj, params = p)
#' }
#' }
#'
#' @export
save.parameters <- function(obj = combined.obj, params = p, overwrite = TRUE) {
  obj <- ww.get.1st.Seur.element(obj)

  if (!is.null(obj@misc$"p") && overwrite) {
    print("Overwriting already existing obj@misc$p. Old version:")
    print(head(unlist(obj@misc$"p")))
    obj@misc$p <- params
  } else if (is.null(obj@misc$"p")) {
    obj@misc$p <- params
  }
}



# _________________________________________________________________________________________________
# List level metadata for ______________________________ ----
# _________________________________________________________________________________________________


#' @title Create Single-Cell Metadata Object for a collection of Seurat Objects
#'
#' @description This function creates a metadata object to correspond to a list of
#'   single-cell experiments, for storing parent level information.
#'   It initializes the object with the experiment and project name, and the
#'   creation date. The created object is of class 'scMetadata_class'.
#' @param experiment The name of the experiment for which metadata is being created.
#' @param project_ The project information to be associated with the metadata object.
#'   This defaults to the current project obtained using Seurat.utils::getProject().
#'
#' @return An 'scCollectionMetadata_class' object containing the metadata for a collection of experiment.
#' @export
#'
#' @examples
#' sc_meta <- create_scCombinedMeta(experiment = "Experiment1")
create_scCombinedMeta <- function(experiment, project_ = getProject()) {
  x <- list(
    experiment.corresponding = experiment,
    initialized = format(Sys.time(), format = "%Y.%m.%d | %H:%M:%S"),
    project = project_
  )
  class(x) <- "scCollectionMetadata_class"
  print(x)
  return(x)
}


# _________________________________________________________________________________________________
# Merging objects and @misc ______________________________ ----
# _________________________________________________________________________________________________


#' @title Copy Specified Elements from One Seurat Object's @misc to Another's
#'
#' @description Copies specified elements from the `@misc` slot of one Seurat object to the `@misc` slot
#' of another. It warns if some specified elements are missing in the source object or if elements are
#' overwritten in the destination object, depending on the `overwrite` argument.
#'
#' @param obj.from The source Seurat object from which elements in the `@misc` slot are to be copied.
#' @param obj.to The destination Seurat object to which elements in the `@misc` slot are to be copied.
#' @param elements.needed A vector of strings specifying the names of the elements in the `@misc` slot of
#' `obj.from` that should be copied to `obj.to`.
#' @param overwrite Logical indicating whether to overwrite elements in `obj.to` that already exist.
#' If `TRUE`, existing elements will be overwritten with a warning; if `FALSE`, the function will
#' stop with an error if it tries to copy an element that already exists in `obj.to`.
#' @return Returns the modified destination Seurat object (`obj.to`) with the specified elements
#' added to or updated in its `@misc` slot.
#' @examples
#' # Assuming `obj1` and `obj2` are Seurat objects and you wish to copy specific elements
#' # from obj1 to obj2, possibly overwriting existing elements in obj2
#' obj2 <- copyMiscElements(obj1, obj2, c("element1", "element2"), overwrite = TRUE)
#'
#' @export
copyMiscElements <- function(obj.from, obj.to, elements.needed, overwrite = TRUE) {
  obj.from <- ww.get.1st.Seur.element(obj.from)

  stopifnot(
    inherits(obj.from, "Seurat"),
    inherits(obj.to, "Seurat")
  )

  # Check for missing elements in obj.to@misc
  elements.from <- names(obj.from@misc)
  missing <- setdiff(elements.needed, elements.from)
  if (length(missing) > 0) {
    warning("Missing elements in obj.from@misc: ", paste(missing, collapse = ", "), immediate. = TRUE)
  }

  # Check for existing elements in obj.to@misc
  elements.already.existing <- intersect(elements.needed, names(obj.to@misc))
  if (length(elements.already.existing) > 0) {
    if (!overwrite) {
      stop(
        "The following elements already exist in obj.to@misc and 'overwrite' is FALSE: ",
        paste(elements.already.existing, collapse = ", ")
      )
    } else {
      warning("Overwriting the following elements in obj.to@misc: ",
        paste(elements.already.existing, collapse = ", "),
        immediate. = TRUE
      )
    }
  }

  # Copy specified elements from obj.from to obj.to
  existingElementsFrom <- intersect(elements.needed, names(obj.from@misc))
  for (element in existingElementsFrom) {
    obj.to@misc[[element]] <- obj.from@misc[[element]]
  }
  iprint("@misc contains: ", names(obj.to@misc))

  return(obj.to)
}


# _________________________________________________________________________________________________
#' @title Copy Tools Slots from Multiple Seurat Objects
#'
#' @description This function copies the `@tools` slots from a list of Seurat objects into a new slot
#' of a target Seurat object. This allows for the aggregation of tools information from multiple
#' experiments or datasets into a single, consolidated Seurat object.
#'
#' @param ls.obj A list of Seurat objects from which the `@tools` slots will be copied.
#' @param obj.to The target Seurat object to which the `@tools` slots will be added.
#' @param overwrite A logical parameter that is kept for compatibility but not used in this version.
#' Its presence does not affect the function's behavior.
#' @param new.slot The name of the new slot within `obj.to@tools` where the copied `@tools` information
#' will be stored. This allows for the organization of copied tools under a specific label, facilitating
#' easy access and interpretation.
#' @return Returns the modified target Seurat object (`obj.to`) with a new `@tools` slot containing
#' the copied information from the list of Seurat objects.
#' @examples
#' # Assuming `ls.obj` is a list of Seurat objects and `obj.to` is a target Seurat object
#' obj.to <- copyCompleteToolsSlots(ls.obj, obj.to, overwrite = TRUE, new.slot = "per.experiment")
#' @export
copyCompleteToolsSlots <- function(ls.obj, obj.to, overwrite = TRUE, new.slot = "per.experiment") {
  stopifnot(
    inherits(obj.to, "Seurat"),
    all(sapply(ls.obj, inherits, "Seurat"))
  )

  ls.tools <- lapply(ls.obj, function(x) x@tools)
  obj.to@tools[[new.slot]] <- ls.tools

  return(obj.to)
}




# _________________________________________________________________________________________________
# Subsetting the Seurat object ______________________________ ----
# _________________________________________________________________________________________________


#' @title Subset a Seurat Object by Identity
#'
#' @description Subsets a Seurat object based on a specified identity column and values. It allows
#'   for an optional inversion of the selection.
#'
#' @param obj A Seurat object. Default: `NULL`.
#' @param ident The name of the identity column to use for subsetting. It is recommended to
#'   specify this explicitly. Default: First entry from the result of `GetClusteringRuns()`.
#' @param identGroupKeep A vector of cluster values for which cells should be matched and retained.
#'   This parameter does not have a default value and must be specified.
#' @param invert A logical indicating whether to invert the selection, keeping cells that do
#'   not match the specified identity values. Default: `FALSE`.
#'
#' @return A Seurat object subsetted based on the specified identity and identity values.
#'
#' @examples
#' # Assuming `seurat_obj` is your Seurat object and you want to subset based on cluster 1
#' subsetted_obj <- subsetSeuObjByIdent(
#'   obj = seurat_obj, ident = "your_ident_column",
#'   identGroupKeep = c(1), invert = FALSE
#' )
#'
#' @importFrom tictoc tic toc
#' @export
subsetSeuObjByIdent <- function(
    obj = combined.obj,
    ident = GetClusteringRuns()[1],
    identGroupKeep,
    invert = FALSE) {
  tic("subsetSeuObjByIdent")
  # Input checks
  stopifnot(
    "obj must be a Seurat object" = inherits(obj, "Seurat"),
    "ident must be a character and exist in obj@meta.data" = is.character(ident) && ident %in% colnames(obj@meta.data),
    "identGroupKeep must exist in ident" = all(identGroupKeep %in% unique(obj@meta.data[[ident]]))
  )

  identGroupKeep <- if (invert) {
    setdiff(unique(obj@meta.data[[ident]]), identGroupKeep)
  } else {
    identGroupKeep
  }
  message(
    "ident: ", ident, " | ", length(identGroupKeep), " ID-groups selected: ", kppc(head(identGroupKeep)),
    "... | invert: ", invert, "\n"
  )

  idx.cells.pass <- obj@meta.data[[ident]] %in% identGroupKeep
  cellz <- colnames(obj)[idx.cells.pass]

  PCT <- percentage_formatter(length(cellz) / ncol(obj))
  message(
    PCT, " or ", length(cellz), " cells are selected from ", ncol(obj),
    ", using values (max 20): ", kppc(head(identGroupKeep, 20)), ", from ", ident, "."
  )

  x <- subset(x = obj, cells = cellz)
  toc()
  return(x)
}


# _________________________________________________________________________________________________
#' @title downsampleSeuObj
#'
#' @description Subset a compressed Seurat object and save it in the working directory.
#' @param obj A Seurat object to subset. Default: the i-th element of the list 'ls.Seurat'.
#' @param fractionCells The fraction of the object's data to keep. Default: 0.25.
#' @param nCells If set to a number greater than 1, indicates the absolute number of cells to keep.
#' If FALSE, the function uses 'fractionCells' to determine the number of cells. Default: `FALSE`.
#' @param seed A seed for random number generation to ensure reproducible results. Default: 1989.
#' @export
#' @importFrom Stringendo percentage_formatter

downsampleSeuObj <- function(obj = ls.Seurat[[i]], fractionCells = 0.25, nCells = FALSE,
                             seed = 1989) {
  set.seed(seed)
  if (isFALSE(nCells)) {
    cellIDs.keep <- sampleNpc(metaDF = obj@meta.data, pc = fractionCells)
    iprint(
      length(cellIDs.keep), "or", Stringendo::percentage_formatter(fractionCells),
      "of the cells are kept. Seed:", seed
    )
  } else if (nCells > 1) {
    nKeep <- min(ncol(obj), nCells)
    # print(nKeep)
    cellIDs.keep <- sample(colnames(obj), size = nKeep, replace = FALSE)
    if (nKeep < nCells) {
      iprint(
        "Only", nCells,
        "cells were found in the object, so downsampling is not possible."
      )
    }
  }
  obj <- subset(x = obj, cells = cellIDs.keep) # downsample
  return(obj)
}

# _________________________________________________________________________________________________
#' @title downsampleSeuObj.and.Save
#'
#' @description Downsample a Seurat object to a target fraction and save it.
#' @param obj Seurat object, Default: ORC
#' @param fraction Fractional size to downsample to. Default: 0.25
#' @param seed random seed used, Default: 1989
#' @param min.features Minimum features
#' @param dir Directory to save to. Default: OutDir
#' @param suffix A suffix added to the filename, Default: ''
#' @export
downsampleSeuObj.and.Save <- function(
    obj = ORC, fraction = 0.25, seed = 1989, dir = OutDir,
    min.features = p$"min.features", suffix = fraction,
    nthreads = .getNrCores()) {
  obj_Xpc <- downsampleSeuObj(obj = obj, fractionCells = fraction, seed = seed)
  nr.cells.kept <- ncol(obj_Xpc)

  # Seurat.utils:::.saveRDS.compress.in.BG(obj = obj_Xpc, fname = ppp(paste0(dir, substitute_deparse(obj),
  # suffix, nr.cells.kept, 'cells.with.min.features', min.features,"Rds" ) )
  xsave(obj_Xpc,
    suffix = ppp(suffix, nr.cells.kept, "cells.with.min.features", min.features),
    nthreads = nthreads, project = getProject(), showMemObject = TRUE, saveParams = FALSE
  )
}



# _________________________________________________________________________________________________
#' @title Sample max number of Cells From each identity in a Seurat Object
#'
#' @description This function samples a specified maximum number of cells from each identity class
#' in a Seurat object, in the meta.data. It ensures that the sampling does not exceed the total
#' number of cells available per identity.
#'
#' @param obj A Seurat object from which cells are to be sampled.
#' @param ident A character vector specifying the identity class from which cells are to be sampled.
#' @param max.cells A positive integer indicating the maximum number of cells to sample from each identity class.
#' @param verbose Logical indicating if messages about the sampling process should be printed to the console. Defaults to TRUE.
#' @param replacement.thr A numeric value between 0 and 1 indicating the percentage of cells to sample from each identity class. Defaults to 0.05.
#' @param dsample.to.repl.thr Logical indicating if sampling should be done with replacement. Defaults to FALSE.
#' @param plot_stats Logical indicating to plot a barplot.
#' @param seed An integer to set the seed for reproducibility.
#'
#'
#' @return Returns a Seurat object containing only the sampled cells.
#'
#' @details This function checks for the presence of the specified identity class within the object's metadata.
#' If the number of cells within any identity class is less than or equal to the `max.cells` parameter,
#' all cells from that class are retained. Otherwise, a random sample of `max.cells` is taken from the class.
#' The function updates the identity of the cells in the returned Seurat object to reflect the sampled cells.
#' If `verbose` is TRUE, it prints the total number of cells sampled and provides a visual summary of the fraction
#' of cells retained per identity class.
#'
#' @examples
#' # Assuming `seuratObj` is a Seurat object with identities stored in its metadata
#' sampledSeuratObj <- downsampleSeuObjByIdentAndMaxcells(obj = seuratObj, ident = "cellType", max.cells = 100)
#'
#' @importFrom CodeAndRoll2 df.col.2.named.vector
#'
#' @export
#'
downsampleSeuObjByIdentAndMaxcells <- function(obj,
                                               ident = GetNamedClusteringRuns()[1],
                                               max.cells = min(table(obj[[ident]])),
                                               verbose = TRUE,
                                               replacement.thr = 0.05,
                                               dsample.to.repl.thr = (max.cells / ncol(obj)) < replacement.thr, # if less than 5% of cells are sampled, sample with replacement
                                               plot_stats = TRUE,
                                               seed = 1989) {
  stopifnot(
    "obj must be a Seurat object" = inherits(obj, "Seurat"),
    "ident must be a character and exist in obj@meta.data" = is.character(ident) && ident %in% colnames(obj@meta.data),
    "max.cells must be a positive integer" = is.numeric(max.cells) && max.cells > 0,
    max.cells < ncol(obj)
  )

  data <- CodeAndRoll2::df.col.2.named.vector(obj[[ident]])
  uniqueCategories <- unique(data)

  set.seed(seed)
  if (dsample.to.repl.thr) {
    max.cells <- round(ncol(obj) * replacement.thr)
    msg <- percentage_formatter(replacement.thr,
      suffix = paste("of the data or", max.cells, "of cells."),
      prefix = "Sampling with replacement to:"
    )
    message(msg)
  }

  # Sample cells from each identity class
  sampledNames <- lapply(uniqueCategories, function(category) {
    namesInCategory <- names(data[data == category])
    if (length(namesInCategory) <= max.cells) {
      # If the number of cells in the category is less than or equal to max.cells, return all cells
      return(namesInCategory)
    } else {
      return(sample(namesInCategory, max.cells))
    }
  })

  sampledCells <- unlist(sampledNames)

  Idents(obj) <- ident
  obj2 <- subset(x = obj, cells = sampledCells)

  subb <- paste0("From ", ncol(obj), " reduced to ", ncol(obj2), " cells.")
  message(subb)

  if (verbose) {
    message("Total cells sampled: ", length(sampledCells))
    print(table(data))

    nr_remaining_cells <- orig_cells <- table(data)
    nr_remaining_cells[nr_remaining_cells > max.cells] <- max.cells
    fr_remaining_per_cluster <- iround(nr_remaining_cells / orig_cells)
    print(fr_remaining_per_cluster)
  }
  if (plot_stats) {
    pobj <- qbarplot(
      vec = fr_remaining_per_cluster, subtitle = subb, label = fr_remaining_per_cluster,
      ylab = "fr. of cells", save = FALSE
    )
    print(pobj)
  }
  return(obj2)
}

# _________________________________________________________________________________________________
#' @title Relabel Small Categories / Clusters
#'
#' @description
#' Relabels small categories in a specified metadata column of a Seurat object. Categories with
#' cell counts less than a minimum count are relabeled to a specified label. The function adds
#' a new metadata column with the updated labels.
#'
#' @param obj Seurat object. The Seurat object containing the metadata.
#' @param col_in Character string. Name of the metadata column to process.
#' @param backup_col_name Character string. Name of the new metadata column where to backup the original values.
#'   Default: `ppp(col_in, "orig")`.
#' @param min_count Numeric. Minimum number of cells required for a category to be retained.
#'   Categories with counts less than this number will be relabeled. Default: 100
#' @param small_label Character string. Label to assign to small categories. Default: "Other".
#' @param v Logical. If `TRUE`, prints verbose output. Default: `TRUE`.
#'
#' @return Seurat object. The modified Seurat object with the new metadata column added.
#'
#' @examples
#' # Assuming 'seurat_obj' is your Seurat object
#' seurat_obj <- RelabelSmallCategories(
#'   obj = seurat_obj,
#'   col_in = "cell_type",
#'   min_count = 50,
#'   small_label = "MinorType",
#'   v = TRUE
#' )
#'
RelabelSmallCategories <- function(
    obj, col_in, backup_col_name = ppp(col_in, "orig"),
    min_count = 100, small_label = "Other", v = TRUE) {
  # Input assertions
  stopifnot(
    inherits(obj, "Seurat"), # Check if obj is a Seurat object
    is.character(col_in), length(col_in) == 1, # col_in is a single string
    col_in %in% colnames(obj@meta.data), # col_in exists in metadata
    is.character(backup_col_name), length(backup_col_name) == 1, # backup_col_name is a single string
    is.numeric(min_count), min_count > 0, # min_count is a positive number
    is.character(small_label), length(small_label) == 1, # small_label is a single string
    is.logical(v), length(v) == 1 # v is a single logical value
  )

  message("backup_col_name: ", backup_col_name)

  # Extract the specified metadata column
  categories <- obj@meta.data[[backup_col_name]] <- obj@meta.data[[col_in]]

  # Count occurrences of each category
  category_counts <- table(categories)

  # Identify small categories
  small_categories <- names(category_counts[category_counts < min_count])

  new_categories <- as.character(categories) # Copy original categories

  # Relabel small categories
  new_categories[new_categories %in% small_categories] <- small_label

  # Add new column to metadata
  obj@meta.data[[col_in]] <- new_categories

  if (v) { # Verbose output
    total_cells <- length(categories)
    num_small_categories <- length(small_categories)
    num_large_categories <- length(category_counts) - num_small_categories
    percent_small_cells <- sum(category_counts[small_categories]) / total_cells * 100
    message(total_cells, " Total cells.")
    message(num_small_categories, " Categories relabeled to ", small_label)
    message(num_large_categories, " of ", length(category_counts), " Categories retained.")
    message(sprintf("Cells in relabeled categories: %d (%.2f%% of total)", sum(category_counts[small_categories]), percent_small_cells))
  }

  return(obj) # Return the modified Seurat object
}


# _________________________________________________________________________________________________
#' @title Remove Residual Small Clusters from a Seurat Object
#'
#' @description Removes clusters containing fewer cells than specified by `max.cells`
#' from a Seurat object. This function is particularly useful after subsetting a dataset,
#' where small, possibly unrepresentative clusters may remain.
#'
#' @param obj Seurat object from which small clusters will be removed. Default: `combined.obj`.
#' @param identities Vector of clustering identities to examine for small clusters;
#' Default: `GetClusteringRuns(obj)`.
#' @param max.cells Maximum number of cells a cluster can contain to still be considered for removal.
#' Default: The lesser of 0.5% of the dataset or 10 cells.
#' @param plot.removed Logical indicating if a umap of the cells removed should be plotted.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   combined.obj <- removeResidualSmallClusters(obj = combined.obj)
#' }
#' }
#'
#' @export
removeResidualSmallClusters <- function(
    obj = combined.obj,
    identities = GetClusteringRuns(obj, pat = "*snn_res.[0-9].[0-9]$")[1:5],
    max.cells = max(round((ncol(obj)) / 2000), 5),
    plot.removed = TRUE) {
  #
  stopifnot(
    inherits(obj, "Seurat"),
    is.character(identities), length(identities) > 0,
    all(identities %in% colnames(obj@meta.data)),
    is.numeric(max.cells), max.cells > 0,
    is.logical(plot.removed)
  )
  message("First, consider using RelabelSmallCategories() instead, to relabel small clusters.")

  META <- obj@meta.data
  all.cells <- rownames(META)


  message("max.cells: ", max.cells, " | Scanning over these identities:")
  small.clusters <- cells.to.remove <- CodeAndRoll2::list.fromNames(identities)

  for (i in 1:length(identities)) {
    colX <- identities[i]
    print(colX)
    tbl <- table(META[[colX]])

    small.clusters[[i]] <- which_names(tbl <= max.cells)
    cells.to.remove[[i]] <- all.cells[which(META[[colX]] %in% small.clusters[[i]])]
    if (length(cells.to.remove[[i]])) {
      iprint(
        length(cells.to.remove[[i]]), "cells in small clusters:", small.clusters[[i]],
        "| Cell counts:", tbl[small.clusters[[i]]]
      )
    }

    all.cells.2.remove <- unique(unlist(cells.to.remove))
    if (plot.removed) {
      SBT <- paste(length(all.cells.2.remove), "cells removed from small clusters across", length(identities), "identities.")
      pobj <- clUMAP(
        obj = obj, ident = identities[i],
        sub = SBT, caption = NULL,
        cells.highlight = all.cells.2.remove
      )
      print(pobj)
    } # if plot.removed


    if (length(all.cells.2.remove)) {
      iprint(
        ">>> a total of", length(all.cells.2.remove),
        "cells are removed which belonged to a small cluster in any of the identities."
      )
    } else {
      iprint(">>> No cells are removed because belonging to small cluster.")
    }

    cells.2.keep <- setdiff(all.cells, all.cells.2.remove)
    obj <- subset(x = obj, cells = cells.2.keep)
  } # for list of identities (meta columns)


  return(obj)
}


# _________________________________________________________________________________________________
#' @title dropLevelsSeurat
#'
#' @description Drop unused levels from `factor` variables in a Seurat object's meta.data.
#' @param obj A Seurat object.
#' @param verbose Logical. Whether to print a message indicating which levels are being dropped.
#' @param only Character vector. Explicit list of columns to only in the operation.
#' @param exclude Character vector. Names of columns to exclude from the operation.#'
#'
#' @export
dropLevelsSeurat <- function(obj = combined.obj, verbose = TRUE, also.character = FALSE,
                             only = NULL, exclude = NULL) {
  stopifnot(is(obj, "Seurat"))
  META <- obj@meta.data
  names.meta <- colnames(obj@meta.data)

  stopifnot(
    is.logical(verbose),
    is.logical(also.character),
    is.null(only) | only %in% names.meta,
    is.null(exclude) | exclude %in% names.meta
  )

  colclasses <- sapply(META, class)
  drop_in_these <- names(colclasses[colclasses %in% "factor"])

  if (!is.null(only)) drop_in_these <- only
  if (!is.null(exclude)) drop_in_these <- setdiff(drop_in_these, exclude)

  if (verbose) {
    message(
      "Dropping levels in ", length(drop_in_these), " identities:\n",
      kppc(drop_in_these)
    )
  }

  for (i in 1:length(drop_in_these)) {
    colX <- drop_in_these[i]
    META[[colX]] <- droplevels(META[[colX]])
  }

  obj@meta.data <- META
  return(obj)
}


# ____________________________________________________________________
#' @title Remove Clusters and Drop Levels from a List of Seurat Objects
#'
#' @description This function removes residual small clusters from specified Seurat objects and
#' drops levels in factor-like metadata.
#' @param ls_obj A list of Seurat objects.
#' @param object_names A character vector containing the names of the Seurat objects to process.
#' Default: names of all objects in the `ls_obj`.
#' @param indices A numeric vector indicating which datasets to process by their position in
#' the `object_names` vector. By default, it processes the second and third datasets.
#' @param ... Additional parameters passed to the `removeResidualSmallClusters` function.
#'
#' @details This function applies `removeResidualSmallClusters` and `dropLevelsSeurat` to
#' the Seurat objects specified by the `indices` in the `object_names`.
#' It operates in place, modifying the input `ls_obj` list.
#'
#' @return The function returns the modified list of Seurat objects.
#' @examples
#' \dontrun{
#' # Process the 2nd and 3rd datasets
#' removeClustersAndDropLevels(ls_obj, indices = c(2, 3))
#' }
#'
#' @export
removeClustersAndDropLevels <- function(ls_obj,
                                        object_names = names(ls_obj),
                                        indices = 2:3, ...) {
  #
  for (index in indices) {
    dataset_name <- object_names[index]
    obj <- ls_obj[[dataset_name]]
    obj <- removeResidualSmallClusters(obj = obj, identities = GetClusteringRuns(obj), ...)
    obj <- dropLevelsSeurat(obj)
    ls_obj[[dataset_name]] <- obj
  }
  return(ls_obj)
}




# _________________________________________________________________________________________________
#' @title Remove Cells by Dimension Reduction
#'
#' @description This function applies a cutoff in the specified dimension of a given
#' dimension reduction (UMAP, PCA, or t-SNE) to remove cells.
#' @param reduction A string specifying the dimension reduction technique to be used
#' ('umap', 'pca', or 'tsne'). Default: 'umap'.
#' @param umap_dim An integer specifying which dimension (axis) to apply the cutoff. Default: 1.
#' @param obj A Seurat object. Default: 'combined.obj'.
#' @param cutoff A numerical value indicating the cutoff value for the specified dimension. Default: 0.
#' @param cut_below A logical value indicating whether to remove cells below (TRUE) or
#' above (FALSE) the cutoff line. Default: `TRUE`.
#' @param only_plot_cutoff Simulate and plot cutoff only.
#' @param ... Any other parameters to be passed to internally called functions.
#' @return A Seurat object with cells removed according to the specified cutoff.
#' @export
removeCellsByUmap <- function(
    reduction = "umap",
    umap_dim = 1,
    obj = combined.obj,
    cutoff = 0,
    cut_below = TRUE,
    only_plot_cutoff = FALSE,
    ...) {
  # Plot cells
  sfx <- if (cut_below) "below" else "above"
  p <- clUMAP(obj = obj, save.plot = FALSE, sub = paste0("cutoff ", sfx, ": ", cutoff), ...)

  # Add cutoff line to plot
  if (umap_dim == 1) {
    p <- p + geom_vline(xintercept = cutoff)
  } else if (umap_dim == 2) {
    p <- p + geom_hline(yintercept = cutoff)
  }
  print(p)
  qqSave(p, fname = kpp("UMAP.with.cutoff", umap_dim, sfx, cutoff, "png"), h = 7, w = 7)

  if (!only_plot_cutoff) {
    # Retrieve cell embeddings
    cell_embedding <- obj@reductions[[reduction]]@cell.embeddings
    all_cells <- rownames(cell_embedding)
    stopifnot(ncol(cell_embedding) > 0)
    embedding_dim_x <- cell_embedding[, umap_dim]

    # Determine cells to remove based on cutoff
    cells_to_remove <- if (cut_below) which_names(embedding_dim_x < cutoff) else which_names(embedding_dim_x >= cutoff)

    # Report on cells removed
    if (length(cells_to_remove)) {
      iprint(">>> A total of", length(cells_to_remove), "cells are removed which fell on UMAP aside cutoff:", cutoff)
    } else {
      iprint(">>> No cells are removed because of the UMAP dimension cutoff.")
    }

    # Subset object to include only cells not removed
    cells_to_keep <- setdiff(all_cells, cells_to_remove)
    obj <- subset(x = obj, cells = cells_to_keep)
  } # only_plot_cutoff
  return(obj)
}




# _________________________________________________________________________________________________
# Downsampling Lists of Seurat objects ______________________________ ----
# _________________________________________________________________________________________________


#' @title Downsample a List of Seurat Objects to a Specific Number of Cells
#'
#' @description Downsampling each Seurat object in a list to a specified number of cells. This function is
#' particularly useful for creating smaller, more manageable subsets of large single-cell datasets for
#' preliminary analyses or testing.
#'
#' @param ls.obj List of Seurat objects to be downsampled. Default: `ls.Seurat`.
#' @param NrCells Target number of cells to downsample each Seurat object to.
#' @param save_object Logical indicating whether to save the downsampled Seurat objects using `isaveRDS`
#' or to return them. Default: `FALSE`.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   downsampledSeuratList <- downsampleListSeuObjsNCells(
#'     ls.obj =
#'       list(yourSeuratObj1, yourSeuratObj2), NrCells = 2000
#'   )
#'   downsampledSeuratList <- downsampleListSeuObjsNCells(NrCells = 200)
#' }
#' }
#'
#' @export
#' @importFrom tictoc tic toc
#' @importFrom Stringendo percentage_formatter
#' @importFrom foreach foreach %dopar% getDoParRegistered

downsampleListSeuObjsNCells <- function(
    ls.obj = ls.Seurat, NrCells = p$"dSample.Organoids",
    save_object = FALSE) {
  # Check if 'ls_obj' is a list of Seurat objects and 'obj_IDs' is a character vector of the same length
  if (!is.list(ls.obj) & inherits(ls.obj, "Seurat")) ls.obj <- list(ls.obj)
  stopifnot(is.list(ls.obj) & all(sapply(ls.obj, function(x) inherits(x, "Seurat"))))

  names.ls <- names(ls.obj)
  n.datasets <- length(ls.obj)
  iprint(NrCells, "cells")

  tictoc::tic("downsampleListSeuObjsNCells")
  if (foreach::getDoParRegistered()) {
    ls.obj.downsampled <- foreach::foreach(i = 1:n.datasets) %dopar% {
      iprint(names(ls.obj)[i], Stringendo::percentage_formatter(i / n.datasets, digitz = 2))
      downsampleSeuObj(obj = ls.obj[[i]], nCells = NrCells)
    }
    names(ls.obj.downsampled) <- names.ls
  } else {
    ls.obj.downsampled <- list.fromNames(names.ls)
    for (i in 1:n.datasets) {
      iprint(names(ls.obj)[i], Stringendo::percentage_formatter(i / n.datasets, digitz = 2))
      ls.obj.downsampled[[i]] <- downsampleSeuObj(obj = ls.obj[[i]], nCells = NrCells)
    }
  } # else
  tictoc::toc()

  print(head(sapply(ls.obj, ncol)))
  print(head(sapply(ls.obj.downsampled, ncol)))

  if (save_object) {
    isave.RDS(obj = ls.obj.downsampled, suffix = ppp(NrCells, "cells"), inOutDir = TRUE)
  } else {
    return(ls.obj.downsampled)
  }
}



# _________________________________________________________________________________________________
#' @title Downsample a List of Seurat Objects to a Fraction
#'
#' @description Downsampling a list of Seurat objects to a specified fraction of their original size.
#' This is useful for reducing dataset size for quicker processing or testing workflows.
#'
#' @param ls.obj List of Seurat objects to be downsampled. Default: `ls.Seurat`.
#' @param fraction Fraction of cells to retain in each Seurat object. Default: 0.1.
#' @param save_object Logical indicating whether to save the downsampled Seurat objects using
#' `isaveRDS` or return them. Default: `FALSE`.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   downsampled_objs <- downsampleListSeuObjsPercent(ls.obj = yourListOfSeuratObjects, fraction = 0.1)
#' }
#' }
#'
#' @export
#' @importFrom tictoc tic toc
#' @importFrom Stringendo percentage_formatter
#' @importFrom foreach foreach %dopar% getDoParRegistered
#'
downsampleListSeuObjsPercent <- function(
    ls.obj = ls.Seurat,
    fraction = 0.1,
    seed = 1989,
    save_object = FALSE) {
  # Check if 'ls_obj' is a list of Seurat objects and 'obj_IDs' is a character vector of the same length
  if (!is.list(ls.obj) & inherits(ls.obj, "Seurat")) ls.obj <- list(ls.obj)
  stopifnot(is.list(ls.obj) & all(sapply(ls.obj, function(x) inherits(x, "Seurat"))))

  names.ls <- names(ls.obj)
  n.datasets <- length(ls.obj)
  iprint(fraction, "fraction")

  tictoc::tic("downsampleListSeuObjsPercent")
  if (foreach::getDoParRegistered()) {
    ls.obj.downsampled <- foreach::foreach(i = 1:n.datasets) %dopar% {
      downsampleSeuObj(obj = ls.obj[[i]], fractionCells = fraction)
    }
    names(ls.obj.downsampled) <- names.ls
  } else {
    ls.obj.downsampled <- list.fromNames(names.ls)
    for (i in 1:n.datasets) {
      cells <- round(ncol(ls.obj[[1]]) * fraction)
      iprint(names(ls.obj)[i], cells, "cells=", Stringendo::percentage_formatter(i / n.datasets, digitz = 2))
      ls.obj.downsampled[[i]] <- downsampleSeuObj(obj = ls.obj[[i]], fractionCells = fraction, seed = seed)
    }
  }
  tictoc::toc() # else

  NrCells <- sum(sapply(ls.obj, ncol))

  print(head(sapply(ls.obj, ncol)))
  print(head(sapply(ls.obj.downsampled, ncol)))
  if (save_object) {
    isave.RDS(obj = ls.obj.downsampled, suffix = ppp(NrCells, "cells"), inOutDir = TRUE)
  } else {
    return(ls.obj.downsampled)
  }
}


# _________________________________________________________________________________________________
# DGEA ______________________________ ----
# _________________________________________________________________________________________________


# _________________________________________________________________________________________________
#' @title Add.DE.combined.score
#'
#' @description Add a combined score to differential expression (DE) results. The score is
#' calculated as log-fold change (LFC) times negative logarithm of scaled
#' p-value (LFC * -log10( p_cutoff / pval_scaling )).
#' @param df A data frame that holds the result of a differential gene expression analysis,
#' typically obtained via the 'FindAllMarkers' function. Default: df.markers.
#' @param p_val_min The minimum p-value considered. All values below this threshold are set to
#' this value. Default: 1e-25.
#' @param pval_scaling The value to scale p-values by in the calculation of the combined score. Default: 0.001.
#' @param colP The name of the column in the input data frame that holds p-values. Default: 'p_val'.
#' @param colLFC The name of the column in the input data frame that holds log-fold change values.
#' By default, it selects the first column not named "avg_logFC" or "avg_log2FC".
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   df.markers <- Add.DE.combined.score(df.markers)
#' }
#' }
#' @export
Add.DE.combined.score <- function(
    df = df.markers, p_val_min = 1e-25, pval_scaling = 0.001, colP = "p_val",
    colLFC = CodeAndRoll2::grepv(pattern = c("avg_logFC|avg_log2FC"), x = colnames(df), perl = TRUE)
    # , colLFC = "avg_log2FC"
    ) { # Score = -LOG10(p_val) * avg_log2FC
  p_cutoff <- SmallestNonAboveX(vec = df[[colP]], X = p_val_min)
  df$"combined.score" <- round(df[[colLFC]] * -log10(p_cutoff / pval_scaling))
  return(df)
}




# _________________________________________________________________________________________________
#' @title Save Top 25 Markers per Cluster
#'
#' @description Stores the top 25 markers for each cluster identified in a Seurat object, based on
#' the `avg_log2FC` from the output table of `FindAllMarkers()`. The result is saved under `@misc$df.markers$res...`,
#' rounding insignificant digits to three decimal places.
#'
#' @param obj Seurat object to update with top 25 markers information. Default: `combined.obj`.
#' @param df_markers Data frame containing results from differential gene expression analysis
#' via `FindAllMarkers()`, specifying significant markers across clusters. Default: `df.markers`.
#' @param res Clustering resolution at which the markers were identified. Default: 0.5.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   combined.obj <- StoreTop25Markers(obj = combined.obj, df_markers = df.markers, res = 0.5)
#' }
#' }
#'
#' @seealso \code{\link[Seurat]{FindAllMarkers}}, \code{\link[dplyr]{top_n}}
#'
#' @export
#' @importFrom Seurat FindAllMarkers
#' @importFrom dplyr group_by top_n select arrange

StoreTop25Markers <- function(
    obj = combined.obj,
    df_markers = df.markers, res = 0.5) {
  top25.markers <-
    df_markers |>
    group_by(cluster) |>
    top_n(n = 25, wt = avg_2logFC) |>
    dplyr::select(gene) |>
    col2named.vec.tbl() |>
    splitbyitsnames()

  obj@misc$"top25.markers"[[ppp("res", res)]] <- top25.markers
  return(obj)
}


# _________________________________________________________________________________________________
#' @title Store All Differential Expression Markers
#'
#' @description Saves the complete output table from `FindAllMarkers()` to a Seurat object, facilitating
#' easy access to differential expression analysis results. This function rounds numerical values to a
#' specified number of digits to maintain readability and manage file sizes.
#'
#' @param obj Seurat object to update with differential expression markers. Default: `combined.obj`.
#' @param df_markers Data frame containing the results from differential gene expression analysis
#' (`FindAllMarkers()` output). Default: `df.markers`.
#' @param res Clustering resolution identifier for storing and referencing the markers. Default: 0.5.
#' @param digit Number of significant digits to retain in numerical values. Default: 3.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   combined.obj <- StoreAllMarkers(obj = combined.obj, df_markers = df.markers, res = 0.5)
#' }
#' }
#'
#' @export
StoreAllMarkers <- function(
    obj = combined.obj,
    df_markers = df.markers, res = 0.5, digit = c(0, 3)[2]) {
  if (digit) df_markers[, 1:5] <- signif(df_markers[, 1:5], digits = digit)
  obj@misc$"df.markers"[[ppp("res", res)]] <- df_markers
  iprint("DF markers are stored under:", "obj@misc$df.markers$", ppp("res", res))
  return(obj)
}


# _________________________________________________________________________________________________
#' @title Get Top Differential Expression Genes Data Frame
#'
#' @description Retrieves a data frame of the top N differentially expressed genes from
#' differential gene expression analysis results, offering an option to exclude certain genes
#' based on patterns.
#'
#' @param dfDE Data frame containing the results of differential gene expression analysis
#' (e.g., output from `FindAllMarkers()`). Default: `df.markers`.
#' @param n Number of top markers to retrieve per cluster. Default: `p$n.markers`.
#' @param order.by Priority column for sorting markers before selection, such as `"avg_log2FC"`;
#' Default: `"avg_log2FC"`.
#' @param exclude Vector of regex patterns to exclude genes from the top markers list;
#' Default: `c("^AL*|^AC*|^LINC*")`.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   topMarkersDF <- GetTopMarkersDF(dfDE = df.markers, n = 3)
#' }
#' }
#'
#' @seealso \code{\link[Seurat]{FindAllMarkers}}, \code{\link[dplyr]{arrange}},
#'  \code{\link[dplyr]{filter}}, \code{\link[dplyr]{group_by}}
#'
#' @export
#' @importFrom dplyr arrange group_by slice select filter

GetTopMarkersDF <- function(
    dfDE = df.markers,
    n = p$"n.markers", order.by = c("avg_log2FC", "p_val_adj")[1],
    exclude = c(
      "^A[CFLP][0-9]{6}", "^Z[0-9]{5}",
      "^LINC0[0-9]{4}", "^C[1-9]+orf[1-9]+",
      "[-|\\.]AS[1-9]*$", "[-|\\.]DT[1-9]*$",
      "^MIR[1-9]", "^SNHG[1-9]"
    )) {
  "Works on active Idents() -> thus we call cluster"
  combined_pattern <- paste(exclude, collapse = "|")

  TopMarkers <- dfDE |>
    dplyr::filter(!grepl(combined_pattern, gene, perl = TRUE)) |>
    arrange(desc(!!as.name(order.by))) |>
    dplyr::group_by(cluster) |>
    dplyr::slice(1:n) |>
    dplyr::select(cluster, gene, avg_log2FC)

  return(TopMarkers)
}


# _________________________________________________________________________________________________
#' @title Get Top Differential Expression Markers from DGEA Results
#'
#' @description Retrieves the top N differentially expressed genes from the results of a differential
#' gene expression analysis, such as that provided by `FindAllMarkers()`.
#'
#' @param dfDE Data frame containing differential expression analysis results. Default: `df.markers`.
#' @param n Number of top markers to retrieve for each cluster. Default: `p$n.markers`.
#' @param order.by Column by which to sort the markers before selection, typically prioritizing
#' markers by significance or effect size. Default: `"avg_log2FC"`.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   topMarkers <- GetTopMarkers(df = df.markers, n = 3)
#' }
#' }
#'
#' @seealso \code{\link[Seurat]{FindAllMarkers}}, \code{\link[dplyr]{arrange}}, \code{\link[dplyr]{group_by}}
#'
#' @export
#' @importFrom dplyr arrange group_by slice select

GetTopMarkers <- function(dfDE = df.markers,
                          n = p$"n.markers",
                          order.by = c("combined.score", "avg_log2FC", "p_val_adj")[2]) {
  message("Works on active Idents()") # thus we call cluster
  TopMarkers <- dfDE |>
    dplyr::arrange(desc(!!as.name(order.by))) |>
    dplyr::group_by(cluster) |>
    dplyr::slice(1:n) |>
    dplyr::select(gene) |>
    CodeAndRoll2::col2named.vec.tbl()

  return(TopMarkers)
}




# _________________________________________________________________________________________________
#' @title AutoLabelTop.logFC
#'
#' @description Create a new "named identity" column in the metadata of a Seurat object,
#' with `Ident` set to a clustering output matching the `res` parameter of the function.
#' It requires the output table of `FindAllMarkers()`.  If you used `StoreAllMarkers()`
#' is stored under `@misc$df.markers$res...`, which location is assumed by default.
#' @param obj A Seurat object, with default value `combined.obj`.
#' @param group.by The clustering group to be used, defaults to the first entry by
#' `GetClusteringRuns()`.
#' @param res Clustering resolution tag. Default: extracted from `group.by`.
#' @param plot.top.genes Logical indicating whether to show a plot, default is `TRUE`.
#' @param suffix Suffix for the naming, defaults to the value of `res`.
#' @param order.by Sorting criterion for the output tibble, defaults to the second element
#' of `c("combined.score", "avg_log2FC", "p_val_adj")`.
#' @param exclude A vector of regular expressions to specify genes to exclude, with
#' default value `c("^AL*|^AC*|^LINC*|^C[0-9]+orf[0-9]*")`.
#' @param df_markers Data frame resulting from DGEA analysis (`FindAllMarkers`). The default
#' is `combined.obj@misc$df.markers[[paste0("res.", res)]]`.
#' @param plotEnrichment Logical indicating whether to plot enrichment, default is `TRUE`.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   combined.obj <- AutoLabelTop.logFC()
#'   combined.obj$"cl.names.top.gene.res.0.5"
#' }
#' }

#' @export
AutoLabelTop.logFC <- function(
    obj = combined.obj,
    group.by,
    res = stringr::str_extract(group.by, "\\d+\\.\\d+"),
    plot.top.genes = TRUE,
    suffix = res,
    order.by = c("combined.score", "avg_log2FC", "p_val_adj")[2],
    exclude = c(
      "^A[CFLP][0-9]{6}", "^Z[0-9]{5}",
      "^LINC0[0-9]{4}", "^C[1-9]+orf[1-9]+",
      "[-|\\.]AS[1-9]*$", "[-|\\.]DT[1-9]*$",
      "^MIR[1-9]", "^SNHG[1-9]"
    ),
    df_markers = obj@misc$"df.markers"[[paste0("res.", res)]],
    plotEnrichment = TRUE) {
  message(group.by)
  message(" > Running AutoLabelTop.logFC...")

  stopifnot(
    !is.null("df_markers"),
    order.by %in% colnames(df_markers)
  )

  df.top.markers <- GetTopMarkersDF(dfDE = df_markers, order.by = order.by, n = 1, exclude = exclude)

  # Enrichment plot ______________________________________________________________
  if (plotEnrichment) {
    top_log2FC <- df.top.markers$"avg_log2FC"
    names(top_log2FC) <- ppp(df.top.markers$"cluster", df.top.markers$"gene")
    ggExpress::qbarplot(top_log2FC,
      plotname = "The strongest fold change by cluster",
      label = iround(top_log2FC),
      subtitle = group.by,
      ylab = "avg_log2FC", xlab = "clusters",
      hline = 2,
      suffix = group.by
    )
  }

  top.markers <- col2named.vec.tbl(df.top.markers[, 1:2])

  obj@misc[[ppp("top.markers.res", res)]] <- top.markers

  ids_CBC <- deframe(obj[[group.by]])
  ids <- unique(ids_CBC)

  # Check if all clusters have DE-genes ____________________________________________________
  if (length(ids) != length(top.markers)) {
    warning("Not all clusters returned DE-genes!", immediate. = TRUE)
    missing <- setdiff(ids, names(top.markers))
    names(missing) <- missing
    iprint("missing:", missing)
    top.markers <- sortbyitsnames(c(top.markers, missing))
  }

  top.markers.ID <- ppp(names(top.markers), top.markers)
  names(top.markers.ID) <- names(top.markers)
  named.group.by <- top.markers.ID[ids_CBC]

  # Check if the clustering was ordered _____________________________________________________
  sfx.ord <- ifelse(grepl("ordered", group.by), group.by, "")
  namedIDslot <- sppp("cl.names.top.gene.", sfx.ord)

  obj <- addMetaDataSafe(obj = obj, metadata = as.character(named.group.by), col.name = namedIDslot, overwrite = TRUE)
  if (plot.top.genes) multiFeaturePlot.A4(list.of.genes = top.markers, suffix = suffix, obj = obj)

  return(obj)
}






# _________________________________________________________________________________________________
#' @title AutoLabel.KnownMarkers
#'
#' @description Creates a new "named identity" column in the metadata of a Seurat object,
#'  setting 'Ident' to a clustering output matching the 'res' parameter.
#'  This function requires the output table of `FindAllMarkers()`.
#' If you used `StoreAllMarkers()`, the output is stored under `@misc$df.markers$res...`,
#' which is the default location.
#' @param obj A Seurat object to work with. Default: `combined.obj`.
#' @param topN The top 'N' genes to consider. Default: 1.
#' @param res The clustering resolution to use. Default: 0.5.
#' @param KnownMarkers A character vector containing known marker genes to be used for annotation.
#' Default: `c("TOP2A", "EOMES", "SLA", "HOPX", "S100B", "DLX6-AS1", "POU5F1", "SALL4", "DDIT4",`
#' `"PDK1", "SATB2", "FEZF2")`.
#' @param order.by Specifies the column to sort the output tibble by.
#' Default: 'combined.score' (First among "combined.score", "avg_log2FC", "p_val_adj").
#' @param df_markers The data frame of markers. By default, it is stored under
#'  `@misc$df.markers$res...` in the provided Seurat object.
#'  Default: `combined.obj@misc$df.markers[[paste0("res.", res)]]`.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   combined.obj <- AutoLabel.KnownMarkers()
#'   DimPlot.ClusterNames(ident = "cl.names.KnownMarkers.0.5")
#' }
#' }
#' @seealso
#'  \code{\link[dplyr]{select}}, \code{\link[dplyr]{slice}}
#' @export
#' @importFrom dplyr select slice
AutoLabel.KnownMarkers <- function(
    obj = combined.obj, topN = 1, res = 0.5,
    KnownMarkers = c(
      `dl-EN` = "KAZN", `ul-EN` = "SATB2", `Immature neurons` = "SLA",
      Interneurons = "DLX6-AS1", Interneurons = "ERBB4", InterN_CGE = "SCGN",
      `Intermediate progenitor` = "EOMES",
      `S-phase` = "TOP2A", `G2M-phase` = "H4C3" # formerly: HIST1H4C
      , `oRG` = "HOPX", Astrocyte = "S100B",
      `Hypoxia/Stress` = "DDIT4", Glycolytic = "PDK1",
      `Choroid.Plexus` = "TTR", `Low-Quality` = "POLR2A",
      `Mesenchyme` = "DCN", `Choroid.Plexus` = "TTR"
    ),
    order.by = c("combined.score", "avg_log2FC", "p_val_adj")[1],
    df_markers = obj@misc$"df.markers"[[paste0("res.", res)]]) {
  stopifnot(!is.null("df_markers"))

  lfcCOL <- CodeAndRoll2::grepv(pattern = c("avg_logFC|avg_log2FC"), x = colnames(df_markers), perl = TRUE)
  keep <- unique(c(lfcCOL, "p_val_adj", "cluster", order.by, "gene"))


  matching.clusters <-
    df_markers |>
    dplyr::select(keep) |>
    arrange(desc(!!as.name(order.by))) |>
    dplyr::filter(gene %in% KnownMarkers) |>
    group_by(gene) |>
    dplyr::slice(1:topN) |>
    arrange(desc(!!as.name(order.by))) |>
    # top_n(n = 1, wt = avg_log2FC) |> # Select the top cluster for each gene
    arrange(cluster)

  print(matching.clusters)

  unique.matches <-
    matching.clusters |>
    group_by(cluster) |> # Select rows with unique values based on column "cluster"
    distinct(cluster, .keep_all = TRUE) |>
    dplyr::select(gene)

  print("Best matches:")
  print(unique.matches)

  "Error Here"
  "Error Here"
  "Error Here"
  "Error Here"
  "Error Here"
  "Error Here"

  top.markers.df <- GetTopMarkersDF(dfDE = df_markers, order.by = lfcCOL, n = 1)
  top.markers <- top.markers.df |> col2named.vec.tbl()

  missing.annotations <-
    top.markers.df |>
    dplyr::filter(!cluster %in% unique.matches$cluster) # filter for clusters that do not have a unique label already

  named.annotations <-
    rbind(unique.matches, missing.annotations) |> # merge the 2 df's
    arrange(cluster) |>
    CodeAndRoll2::col2named.vec.tbl()

  (top.markers.ID <- ppp(names(named.annotations), named.annotations))
  names(top.markers.ID) <- names(top.markers)
  named.ident <- top.markers.ID[Idents(object = obj)]

  namedIDslot <- ppp("cl.names.KnownMarkers", res)
  obj[[namedIDslot]] <- named.ident
  return(obj)
}



# _________________________________________________________________________________________________
# Correlations _________________________ ----
# _________________________________________________________________________________________________

#' @title Calculate Sparse Correlation Matrix
#'
#' @description Computes a sparse correlation matrix from a given sparse matrix input. This function is
#' useful for efficiently handling large datasets where most values are zero, facilitating the calculation
#' of both covariance and correlation matrices without converting to a dense format.
#'
#' @param smat A sparse matrix object, typically of class Matrix from the Matrix package.
#' @return A list with two elements:
#'   * `cov`: The covariance matrix derived from the input sparse matrix.
#'   * `cor`: The correlation matrix derived from the covariance matrix.
#'
#' @examples
#' \dontrun{
#' library(Matrix)
#' smat <- Matrix(rnorm(1000), nrow = 100, sparse = TRUE)
#' cor_res <- sparse.cor(smat)
#' print(cor_res$cor)
#' }
#'
#' @export
#' @importFrom Matrix colMeans crossprod tcrossprod
#' @importFrom stats sd
sparse.cor <- function(smat) {
  n <- nrow(smat)
  cMeans <- colMeans(smat)
  covmat <- (as.matrix(crossprod(smat)) - n * tcrossprod(cMeans)) / (n - 1)
  sdvec <- sqrt(diag(covmat))
  cormat <- covmat / tcrossprod(sdvec)
  list(cov = covmat, cor = cormat)
}


# _________________________________________________________________________________________________
#' @title Calc.Cor.Seurat
#'
#' @description Calculate gene correlation on a Seurat object.
#' @param assay.use The assay to use from the Seurat object. Default: 'RNA'
#' @param slot.use The slot to use from the assay in the Seurat object. Default: 'data'
#' @param quantileX The quantile level for the calculation. Default: `0.95`
#' @param max.cells Maximum number of cells to be used in the calculation. Default: `40000`
#' @param seed The random seed used for the calculation. Default: `p$seed`
#' @param digits The number of decimal places to round the correlation and covariance values. Default: `2`
#' @param obj The Seurat object to perform calculations on. Default: `combined.obj`
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   combined.obj <- calc.q99.Expression.and.set.all.genes(combined.obj, quantileX = 0.99, max.cells = 400000, set.all.genes = FALSE)
#'   combined.obj <- Calc.Cor.Seurat(assay.use = "RNA", slot.use = "data", digits = 2, obj = combined.obj, quantile = 0.99, max.cells = 40000)
#' }
#' }
#' @importFrom tictoc tic toc
#'
#' @export
Calc.Cor.Seurat <- function(
    assay.use = "RNA",
    slot.use = "data",
    quantileX = 0.95,
    max.cells = 40000,
    seed = p$"seed",
    digits = 2, obj = combined.obj) {
  expr.mat <- GetAssayData(slot = slot.use, assay = assay.use, object = obj)
  if (ncol(expr.mat) > max.cells) {
    set.seed(seed = seed)
    cells.use <- sample(x = colnames(expr.mat), size = max.cells)
  } else {
    cells.use <- ncol(expr.mat)
  }

  qname <- paste0("q", quantileX * 100)
  quantile_name <- kpp("expr", qname)

  if (is.null(obj@misc[[quantile_name]])) {
    iprint(
      "Call: combined.obj <- calc.q99.Expression.and.set.all.genes(combined.obj, quantileX =",
      quantileX, " first )"
    )
  }
  genes.HE <- which_names(obj@misc[[quantile_name]] > 0)
  iprint("Pearson correlation is calculated for", length(genes.HE), "HE genes with expr.", qname, ": > 0.")
  tictoc::tic("sparse.cor")
  ls.cor <- sparse.cor(smat = t(expr.mat[genes.HE, cells.use]))
  tictoc::toc()
  ls.cor <- lapply(ls.cor, round, digits = 2)

  slot__name <- kpp(slot.use, assay.use, quantile_name)
  obj@misc[[kpp("cor", slot__name)]] <- ls.cor$"cor"
  obj@misc[[kpp("cov", slot__name)]] <- ls.cor$"cov"
  iprint("Stored under obj@misc$", kpp("cor", slot.use, assay.use), "or cov... .")
  return(obj)
}


# _________________________________________________________________________________________________
#' @title Plot Gene Correlation Heatmap
#'
#' @description Generates a heatmap visualization of gene correlations based on expression data.
#' Useful for identifying groups of genes that exhibit similar expression patterns across different conditions
#' or cell types in a Seurat object.
#'
#' @param genes Vector of gene symbols to include in the correlation analysis and heatmap.
#' @param assay.use Assay from which to retrieve expression data within the Seurat object. Default: 'RNA'.
#' @param slot.use Specifies which slot of the assay to use for expression data `('data', 'scale.data', 'data.imputed')`;
#' Default: first item `('data')`.
#' @param quantileX Quantile level for calculating expression thresholds. Default: `0.95`.
#' @param min.g.cor Minimum absolute gene correlation value for inclusion in the heatmap. Default: `0.3`.
#' @param calc.COR Logical flag to calculate correlation matrix if not found in `@misc`. Default: `FALSE.`
#' @param cutRows Height at which to cut the dendrogram for rows, determining cluster formation. Default: `NULL.`
#' @param cutCols Height at which to cut the dendrogram for columns, determining cluster formation.
#' Default: same as `cutRows`.
#' @param obj Seurat object containing the data. Default: `combined.obj`.
#' @param ... Additional parameters passed to the internally called functions.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   plot.Gene.Cor.Heatmap(genes = c("Gene1", "Gene2", "Gene3"), obj = combined.obj)
#' }
#' }
#'
#' @importFrom Seurat GetAssayData
#' @importFrom pheatmap pheatmap
#' @importFrom MarkdownReports wplot_save_pheatmap
#'
#' @export plot.Gene.Cor.Heatmap
#'
plot.Gene.Cor.Heatmap <- function(
    genes,
    assay.use = "RNA", slot.use = c("data", "scale.data", "data.imputed")[1], quantileX = 0.95,
    min.g.cor = 0.3, calc.COR = FALSE,
    cutRows = NULL, cutCols = cutRows,
    obj = combined.obj, ...) {
  expr.mat <- GetAssayData(slot = slot.use, assay = assay.use, object = obj)
  if (slot.use == c("data.imputed")) {
    "WIP"
  }
  expr.mat <- GetAssayData(slot = slot.use, assay = assay.use, object = obj)

  qname <- paste0("expr.q", quantileX * 100)
  slotname_cor.mat <- kpp("cor", slot.use, assay.use, qname)
  cor.mat <- obj@misc[[slotname_cor.mat]]

  if (is.null(cor.mat)) {
    iprint(slotname_cor.mat, " not found in @misc.")
    iprint("Correlation slots present in @misc:", CodeAndRoll2::grepv(names(obj@misc), pattern = "^cor"))

    # Calculate --- --- --- --- ---
    if (calc.COR) {
      message("Calculating correlation now.")
      genes.found <- check.genes(genes = genes)
      message(length(genes.found), " genes are found in the object.")

      if (length(genes.found) > 200) iprint("Too many genes found in data, cor will be slow: ", length(genes.found))
      ls.cor <- sparse.cor(t(expr.mat[genes.found, ]))
      cor.mat <- ls.cor$cor
    } else {
      stop()
    }
  } else {
    print("Correlation is pre-calculated")
    genes.found <- intersect(genes, rownames(cor.mat))
    iprint(length(genes.found), "genes are found in the correlation matrix.")
    cor.mat <- cor.mat[genes.found, genes.found]
  }


  # Filter --- --- --- --- --- ---
  diag(cor.mat) <- NaN
  corgene.names <- union(
    which_names(rowMax(cor.mat) >= min.g.cor),
    which_names(rowMin(cor.mat) <= -min.g.cor)
  )
  iprint(length(corgene.names), "genes are more (anti-)correlated than +/-:", min.g.cor)

  pname <- paste0("Pearson correlations of ", substitute_deparse(genes), "\n min.cor:", min.g.cor, " | ", assay.use, ".", slot.use)
  o.heatmap <- pheatmap::pheatmap(cor.mat[corgene.names, corgene.names], main = pname, cutree_rows = cutRows, cutree_cols = cutCols, ...)
  MarkdownReports::wplot_save_pheatmap(o.heatmap, plotname = make.names(pname))

  # return values
  maxCorrz <- rowMax(cor.mat)[corgene.names]
  names(maxCorrz) <- corgene.names
  dput(maxCorrz)
}





# _________________________________________________________________________________________________
# Seurat.object.manipulations.etc.R ______________________________ ----
# _________________________________________________________________________________________________
# source('~/GitHub/Packages/Seurat.utils/Functions/Seurat.object.manipulations.etc.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Seurat.object.manipulations.etc.R"))



# _________________________________________________________________________________________________
#' @title Add Prefixes to Cell Names in Seurat Objects
#'
#' @description Adds prefixes derived from a vector of identifiers to cell names in a list of Seurat objects.
#' This is useful for ensuring unique cell names across multiple samples or conditions when combining or comparing datasets.
#'
#' @param ls_obj List of Seurat S4 objects to which prefixes will be added. Each object should correspond
#' to a different sample or condition.
#' @param obj_IDs Character vector of identifiers that will be used as prefixes. Each identifier in the vector
#' corresponds to a Seurat object in `ls_obj`. The length of `obj_IDs` must match the length of `ls_obj`.
#'
#' @examples
#' \dontrun{
#' # Assuming seurat_obj1 and seurat_obj2 are Seurat objects
#' ls_obj <- list(seurat_obj1, seurat_obj2)
#' obj_IDs <- c("sample1", "sample2")
#' ls_obj_prefixed <- prefix_cells_seurat(ls_obj = ls_obj, obj_IDs = obj_IDs)
#' # Now each cell name in seurat_obj1 and seurat_obj2 will be prefixed with 'sample1_' and 'sample2_', respectively.
#' }
#'
#' @return A list of Seurat objects with updated cell names, incorporating the specified prefixes.
#'
#' @export
#' @importFrom Seurat RenameCells
prefix_cells_seurat <- function(ls_obj, obj_IDs) {
  # Check if 'ls_obj' is a list of Seurat objects and 'obj_IDs' is a character vector of the same length
  if (!is.list(ls_obj) & inherits(ls_obj, "Seurat")) ls_obj <- list(ls_obj)
  stopifnot(is.list(ls_obj) & all(sapply(ls_obj, function(x) inherits(x, "Seurat"))))
  stopifnot(is.character(obj_IDs) & length(ls_obj) == length(obj_IDs))

  names_orig <- names(ls_obj)

  # Iterate over Seurat objects
  ls_obj_prefixed <- lapply(seq_along(ls_obj), function(i) {
    # Get the Seurat object and corresponding prefix
    obj <- ls_obj[[i]]
    prefix <- obj_IDs[i]

    # Add prefix to cell names
    new_cell_names <- paste0(prefix, "_", colnames(obj))

    # Rename cells in the Seurat object
    obj <- RenameCells(obj, new.names = new_cell_names)

    return(obj)
  })
  print(lapply(lapply(ls_obj_prefixed, colnames), head))

  names(ls_obj_prefixed) <- names_orig
  return(ls_obj_prefixed)
}


# _________________________________________________________________________________________________
#' @title Check Prefix in Seurat Object Cell IDs
#'
#' @description This function checks if a prefix has been added to the standard
#' cell-IDs (16 characters of A,TRUE,C,G) in a Seurat object. If so, it prints the number of unique prefixes found,
#' issues a warning if more than one unique prefix is found, and returns the identified prefix(es).
#'
#' @param obj A Seurat object with cell IDs possibly prefixed.
#' @param cell_ID_pattern Pattern to match cellIDs (with any suffix).
#' @return A character vector of the identified prefix(es).
#'
#' @examples
#' # Assuming 'obj' is your Seurat object
#' # prefix <- find_prefix_in_cell_IDs(obj)
#'
#' @export
find_prefix_in_cell_IDs <- function(obj, cell_ID_pattern = "[ATCG]{16}.*$") {
  stopifnot(inherits(obj, "Seurat"))

  # Extract cell IDs
  cell_IDs <- colnames(obj)

  # Remove the standard 16-character cell-IDs
  potential_prefixes <- gsub(pattern = cell_ID_pattern, replacement = "", x = cell_IDs)

  # Check if there is no prefix
  if (all(potential_prefixes == "")) {
    print("No prefix found in cell IDs.")
    return(NULL)
  }

  # Identify unique prefixes
  unique_prefixes <- unique(potential_prefixes)

  # Print the number of unique prefixes
  print(paste(length(unique_prefixes), "unique prefix(es) found:", head(unique_prefixes)))

  # Issue a warning if more than one unique prefix is found
  if (length(unique_prefixes) > 1) {
    warning("Multiple unique prefixes identified in cell IDs:", head(unique_prefixes), immediate. = TRUE)
  }

  # Return the identified prefix(es)
  return(unique_prefixes)
}




# _________________________________________________________________________________________________
#' @title Create Cluster Labels for Each Cell
#'
#' @description Generates labels for each cell by combining gene names and cluster IDs. This function
#' takes a named vector, typically representing top genes for clusters (values) and their corresponding
#' cluster IDs (names), along with a vector of cell IDs. It then creates a new vector where each cell
#' is labeled with its top gene and cluster ID in the format "GeneName.ClusterID".
#'
#' @param TopGenes A named vector with gene names as values and cluster IDs as names,
#' representing the top or defining gene for each cluster.
#' @param clID.per.cell A vector of cluster IDs for each cell, used to match each cell with its
#' corresponding top gene from `TopGenes`.
#'
#' @return A vector where each element corresponds to a cell labeled with both its defining gene
#' name and cluster ID, in the format "GeneName.ClusterID".
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Assuming `TopGenes.Classic` is a named vector of top genes and cluster IDs,
#'   # and `metaD.CL.colname` is a column in metadata with cluster IDs per cell
#'   cellLabels <- seu.Make.Cl.Label.per.cell(
#'     TopGenes = TopGenes.Classic,
#'     clID.per.cell = getMetadataColumn(ColName.metadata = metaD.CL.colname)
#'   )
#'   # `cellLabels` now contains labels for each cell in the format "GeneName.ClusterID"
#' }
#' }
#'
#' @export
seu.Make.Cl.Label.per.cell <- function(TopGenes, clID.per.cell) {
  Cl.names_class <- TopGenes[clID.per.cell]
  Cl.names_wNr <- paste0(Cl.names_class, " (", names(Cl.names_class), ")")
  return(Cl.names_wNr)
}


# _________________________________________________________________________________________________
#' @title Retrieve the Top Variable Genes from a Seurat Object
#'
#' @description Retrieves the names of the most variable genes from a Seurat object,
#' typically used to focus subsequent analyses on genes with the greatest variation across cells.
#'
#' @param obj A Seurat object containing gene expression data and,
#' pre-computed highly variable gene information.
#' @param nGenes The number of most variable genes to retrieve. Default: `p$nVarGenes`.
#'
#' @return A vector containing the names of the most variable genes.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Assuming `combined.obj` is a Seurat object with computed variable genes
#'   varGenes <- GetMostVarGenes(obj = combined.obj, nGenes = 100)
#' }
#' }
#'
#' @export
#' @importFrom Seurat FindVariableFeatures
#'
GetMostVarGenes <- function(obj, nGenes = p$nVarGenes) {
  head(rownames(slot(object = obj, name = "hvg.info")), n = nGenes)
}

# _________________________________________________________________________________________________
#' @title Check Gene Names in Seurat Object
#'
#' @description Examines gene names in a Seurat object for specific naming conventions,
#' such as the presence of hyphens (-) or dots (.) often found in mitochondrial gene names.
#' This function is useful for ensuring gene names conform to expected patterns,
#' especially when preparing data for compatibility with other tools or databases.
#'
#' @param Seu.obj A Seurat object containing gene expression data.
#'
#' @details This function prints out examples of gene names that contain specific characters
#' of interest (e.g., '-', '_', '.', '.AS[1-9]'). It is primarily used for data inspection
#' and cleaning before further analysis or data export.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Assuming `combined.obj` is your Seurat object
#'   gene.name.check(Seu.obj = combined.obj)
#'   # This will print examples of gene names containing '-', '_', '.', and '.AS[1-9]'
#' }
#' }
#'
#' @seealso \code{\link[Seurat]{GetAssayData}}
#'
#' @importFrom Seurat GetAssayData
#' @importFrom CodeAndRoll2 grepv
#' @importFrom MarkdownHelpers llprint llogit
#'
#' @export
gene.name.check <- function(Seu.obj) {
  rn <- rownames(GetAssayData(object = Seu.obj, slot = "counts"))
  MarkdownHelpers::llprint("### Gene name pattern")

  MarkdownHelpers::llogit('`rn = rownames(GetAssayData(object = ls.Seurat[[1]], slot = "counts"))`')
  MarkdownHelpers::llogit('`head(CodeAndRoll2::grepv(rn, pattern = "-"), 10)`')
  print("pattern = -")
  MarkdownHelpers::llprint(head(CodeAndRoll2::grepv(rn, pattern = "-"), 10))

  MarkdownHelpers::llogit('`head(CodeAndRoll2::grepv(rn, pattern = "_"), 10)`')
  print("pattern = _")
  MarkdownHelpers::llprint(head(CodeAndRoll2::grepv(rn, pattern = "_"), 10))

  MarkdownHelpers::llogit('`head(CodeAndRoll2::grepv(rn, pattern = "\\."), 10)`')
  print("pattern = \\.")
  MarkdownHelpers::llprint(head(CodeAndRoll2::grepv(rn, pattern = "\\."), 10))

  MarkdownHelpers::llogit('`head(CodeAndRoll2::grepv(rn, pattern = "\\.AS[1-9]"), 10)`')
  print("pattern = \\.AS[1-9]")
  MarkdownHelpers::llprint(head(CodeAndRoll2::grepv(rn, pattern = "\\.AS[1-9]"), 10))
}


# _________________________________________________________________________________________________
#' @title Check if Gene Names exist in Seurat Object or HGNC Database
#'
#' @description Verifies the presence of specified gene names within a Seurat object or
#' queries them against the HGNC database. This function is useful for ensuring gene names are
#' correctly formatted and exist within the dataset or are recognized gene symbols.
#'
#' @param genes A vector of gene names to be checked.
#' @param makeuppercase If `TRUE`, converts all gene names to uppercase before checking. Default: `FALSE`.
#' @param verbose If `TRUE`, prints information about any missing genes. Default: `TRUE`.
#' @param HGNC.lookup If `TRUE`, attempts to look up any missing genes in the HGNC database to
#' verify their existence. Default: `FALSE`.
#' @param obj The Seurat object against which the gene names will be checked.
#' @param assay.slot Assay slot of the Seurat object to check for gene names. Default: `'RNA'`.
#' @param data.slot Data slot of the assay to check for gene names. Default: `'data'`.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Check for the presence of a gene name in uppercase
#'   check.genes(genes = "top2a", makeuppercase = TRUE, obj = combined.obj)
#'
#'   # Check for a gene name with verbose output and HGNC lookup
#'   check.genes(genes = "VGLUT2", verbose = TRUE, HGNC.lookup = TRUE, obj = combined.obj)
#' }
#' }
#'
#' @seealso \code{\link[Seurat]{GetAssayData}}, \code{\link[DatabaseLinke.R]{qHGNC}}
#'
#' @export
#' @importFrom Seurat GetAssayData
#' @importFrom Stringendo percentage_formatter
#'
check.genes <- function(
    genes, makeuppercase = FALSE, HGNC.lookup = FALSE,
    obj,
    assay.slot = c("RNA", "integrated")[1],
    data.slot = c("counts", "data")[2],
    verbose = TRUE,
    ...) {
  tictoc::tic("check.genes")
  message(" > Running check.genes...")
  message("assay: ", assay.slot, ", data.slot: ", data.slot)

  if (makeuppercase) genes <- toupper(genes)

  all_genes <-
    if (obj@version < "5") {
      rownames(GetAssayData(object = obj, assay = assay.slot, slot = data.slot))
    } else {
      rownames(GetAssayData(object = obj, layer = data.slot))
    }

  missingGenes <- setdiff(genes, all_genes)
  if (length(missingGenes) > 0) {
    if (verbose) {
      message(
        "\n", length(missingGenes), " or ",
        Stringendo::percentage_formatter(length(missingGenes) / length(genes)),
        " genes not found in the data, e.g: ", kppc(head(missingGenes, n = 10))
      )
    }

    if (HGNC.lookup) {
      stopifnot("Package 'DatabaseLinke.R' must be installed to use the 'HGNC.lookup' option." = require("DatabaseLinke.R"))
      DatabaseLinke.R::qHGNC(missingGenes, Open = FALSE)
    }
  } else {
    message("All genes found.")
  }

  tictoc::toc()
  intersect_genes <- intersect(genes, all_genes)

  # Using logical indexing to return genes with names (if they had any)
  genes[intersect_genes %in% genes]
}



# _________________________________________________________________________________________________
#' @title Fix Zero Indexing in Seurat Clustering
#'
#' @description Adjusts Seurat object metadata to fix zero-based cluster indexing, converting it to one-based indexing.
#' This function modifies a specified metadata column in the Seurat object to replace zero-indexed cluster names with one-based indexing.
#'
#' @param ColName.metadata The name of the metadata column containing zero-based cluster indices. Default: `'res.0.6'`.
#' @param obj The Seurat object to be modified. Default: `org`.
#'
#' @return The Seurat object with the specified metadata column's cluster indices adjusted to one-based indexing.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Assuming `org` is a Seurat object with zero-based cluster indexing
#'   org <- fixZeroIndexing.seurat(ColName.metadata = "res.0.6", obj = org)
#'   # Now, `org` has its cluster indices in the 'res.0.6' metadata column adjusted to one-based indexing
#' }
#' }
#'
#' @export
fixZeroIndexing.seurat <- function(ColName.metadata = "res.0.6", obj = org) {
  obj@meta.data[, ColName.metadata] <- as.numeric(obj@meta.data[, ColName.metadata]) + 1
  print(obj@meta.data[, ColName.metadata])
  return(obj)
}


# _________________________________________________________________________________________________
#' @title Calculate Fraction of Genes in Transcriptome
#'
#' @description Calculates the fraction of specified genes within the entire transcriptome of
#' each cell in a Seurat object.
#' This function is useful for assessing the relative abundance of a set of genes across cells,
#' such as identifying cells with high expression of marker genes.
#'
#' @param geneset A character vector of gene symbols for which the fraction in the transcriptome will be calculated.
#' Default: `c("MALAT1")`. The function will check for the existence of these genes in the Seurat object.
#' @param obj A Seurat object containing gene expression data. Default: `combined.obj`.
#' The function extracts gene expression data from this object to calculate fractions.
#' @param data.slot The data slot from which to extract expression data. This can be `"counts"`
#' for raw counts or `"data"` for normalized data. Default: second element (`"data"`).
#'
#' @return A numeric vector where each element represents the fraction of the specified geneset's expression
#' relative to the total transcriptome of a cell, expressed as a percentage. The names of the vector correspond to cell IDs.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Assuming `combined.obj` is your Seurat object
#'   fractionInTranscriptome <- CalculateFractionInTranscriptome(geneset = c("MALAT1", "GAPDH"), obj = combined.obj)
#'   # This will return the fraction of MALAT1 and GAPDH in the transcriptome of each cell
#' }
#' }
#'
#' @note This function calls `check.genes` to verify the existence of the specified genes within the Seurat object.
#' If genes are not found, it will return a warning.
#'
#' @seealso \code{\link[Seurat]{GetAssayData}} for retrieving expression data from a Seurat object.
#'
#' @export
#'
CalculateFractionInTrome <- function(
    geneset = c("MALAT1"),
    obj = combined.obj,
    data.slot = c("counts", "data")[2]) {
  warning("    >>>> Use addMetaFraction() <<<<", immediate. = TRUE)
  geneset <- check.genes(genes = geneset)
  stopifnot(length(geneset) > 0)

  mat <- as.matrix(slot(obj@assays$RNA, name = data.slot))
  mat.sub <- mat[geneset, , drop = FALSE]
  RC.per.cell.geneset <- colSums(mat.sub)

  RC.per.cell <- colSums(mat)
  gene.fraction.per.cell <- 100 * RC.per.cell.geneset / RC.per.cell
  return(gene.fraction.per.cell)
}

# _________________________________________________________________________________________________
#' @title AddNewAnnotation
#'
#' @description This function creates a new metadata column based on an existing metadata column
#' and a list of mappings (name <- IDs).
#' @param obj A Seurat object for which the new annotation is to be created. Default: 'obj'.
#' @param source A character string specifying the existing metadata column to be used as the
#' basis for the new annotation. Default: 'RNA_snn_res.0.5'.
#' @param named.list.of.identities A named list providing the mappings for the new annotation.
#' Default: 'ls.Subset.ClusterLists'.
#' @return A character vector representing the new metadata column.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   ls.Subset.ClusterLists <- list("hESC.h9" = c("4", "10", "14"), "hESC.176" = c("0", "1", "2"))
#'   AddNewAnnotation()
#' }
#' }
#' @export
AddNewAnnotation <- function(
    obj = obj,
    source = "RNA_snn_res.0.5", named.list.of.identities = ls.Subset.ClusterLists) {
  NewID <- df.col.2.named.vector(obj[[source]])

  for (i in 1:length(named.list.of.identities)) {
    lx <- as.character(named.list.of.identities[[i]])
    name.lx <- names(named.list.of.identities)[i]
    NewID <- CodeAndRoll2::translate(vec = NewID, old = lx, new = name.lx)
  }
  print(table(NewID))
  return(NewID)
}


# _________________________________________________________________________________________________
#' @title whitelist.subset.ls.Seurat
#'
#' @description Subsets cells in a list of Seurat objects based on an externally provided list of cell IDs.
#' @param ls.obj A list of Seurat objects. Default: ls.Seurat.
#' @param metadir Directory for the metadata. Default: p$cellWhiteList.
#' @param whitelist.file Filename of the whitelist containing cell IDs. Default: "NonStressedCellIDs.2020.10.21_18h.tsv".
#' @return A list of Seurat objects containing only the cells specified in the whitelist.
#' @details The function first validates the presence of all identities from the metadata in the
#' Seurat objects. If all identities are present, the function subsets each Seurat object based on
#' the whitelist of cell IDs.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   ls.Seurat.subset <- whitelist.subset.ls.Seurat(
#'     ls.obj = ls.Seurat, metadir = p$"cellWhiteList",
#'     whitelist.file = "NonStressedCellIDs.2020.10.21_18h.tsv"
#'   )
#' }
#' }
#' @seealso
#' \code{\link[Seurat]{subset}}
#' @importFrom ReadWriter read.simple.tsv
#'
#' @export
whitelist.subset.ls.Seurat <- function(
    ls.obj = ls.Seurat,
    metadir = p$"cellWhiteList" #  '~/Dropbox/Abel.IMBA/MetadataD/POL.meta/cell.lists/'
    , whitelist.file = "NonStressedCellIDs.2020.10.21_18h.tsv") {
  cells.before <- sapply(ls.obj, ncol)
  # Find file
  df.cell.whitelist <- ReadWriter::read.simple.tsv(metadir, whitelist.file)
  dsets <- table(df.cell.whitelist[, 1])

  ls.orig.idents <- lapply(lapply(ls.Seurat, getMetadataColumn, ColName.metadata = "orig.ident"), unique)
  stopif(any(sapply(ls.orig.idents, l) == length(ls.Seurat)), message = "Some ls.Seurat objects have 1+ orig identity.")

  dsets.in.lsSeu <- unlist(ls.orig.idents)
  isMathced <- all(dsets.in.lsSeu == names(dsets)) # Stop if either ls.Seurat OR the metadata has identities not found in the other, in the same order.
  stopif(!isMathced, message = paste(
    "either ls.Seurat OR the metadata has identities not found in the other, or they are not in same order.",
    kpps(dsets.in.lsSeu), "vs.", kpps(names(dsets))
  ))

  # identX <- ls.orig.idents[[1]]
  for (i in 1:length(ls.orig.idents)) {
    identX <- ls.orig.idents[[i]]
    print(identX)

    # Extract and process cellIDs
    idx.match <- which(df.cell.whitelist[, 1] == identX)
    cell.whitelist <- rownames(df.cell.whitelist)[idx.match]
    cell.whitelist <- substr(
      x = cell.whitelist,
      start = 1, stop = nchar(cell.whitelist) - 2
    )

    # Extract and process cellIDs
    ls.obj[[i]] <- subset(x = ls.obj[[i]], cells = cell.whitelist)
  }
  cells.after <- sapply(ls.obj, ncol)
  iprint("cells.before", cells.before, "cells.after", cells.after)
  return(ls.obj)
}

# _________________________________________________________________________________________________
#' @title FindCorrelatedGenes
#'
#' @description Find correlated genes in a Seurat object
#' @param gene Gene of interest. Default: 'TOP2A'
#' @param obj Seurat object to find the correlated genes from. Default: `combined.obj`
#' @param assay Assay to be used from the Seurat object. Default: 'RNA'
#' @param slot Slot to be used from the specified assay in the Seurat object. Default: 'data'
#' @param HEonly Logical, if TRUE, filters matrix to high-expressing genes only. Default: `FALSE`.
#' @param minExpr Minimum expression level for a gene to be considered. Default: 1
#' @param minCells Minimum number of cells expressing a gene for the gene to be considered. Default: 1000
#' @param trailingNgenes Number of top genes to consider based on their correlation. Default: 1000
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   FindCorrelatedGenes(gene = "TOP2A", obj = combined.obj)
#'   write_clip(names(head(topGenes[-(1:6)], n = 50)))
#' }
#' }
#' @seealso
#'  \code{\link[matrixStats]{rowSums2}}
#' @importFrom matrixStats rowSums2
#' @importFrom tictoc tic toc
#' @importFrom MarkdownReports wbarplot
#'
#' @export
FindCorrelatedGenes <- function(
    gene = "TOP2A", obj = combined.obj, assay = "RNA", slot = "data",
    HEonly = FALSE, minExpr = 1, minCells = 1000,
    trailingNgenes = 1000) {
  tictoc::tic("FindCorrelatedGenes")
  AssayData <- GetAssayData(object = obj, assay = assay, slot = slot)
  matrix_mod <- iround(as.matrix(AssayData))
  if (HEonly) {
    idx.pass <- (matrixStats::rowSums2(matrix_mod > minExpr) > minCells)
    pc_TRUE(idx.pass)
    matrix_mod <- matrix_mod[which(idx.pass), ]
  }
  geneExpr <- as.numeric(matrix_mod[gene, ])
  correlations <- apply(matrix_mod, 1, cor, geneExpr)
  topGenes <- trail(sort(correlations, decreasing = TRUE), N = trailingNgenes)
  tictoc::toc()
  MarkdownReports::wbarplot(head(topGenes, n = 25))
  topGenes
}



# _________________________________________________________________________________________________
# _________________________________________________________________________________________________



# _________________________________________________________________________________________________
# Seurat.update.gene.symbols.HGNC.R ______________________________ ----
# _________________________________________________________________________________________________
# source('~/GitHub/Packages/Seurat.utils/Functions/Seurat.update.gene.symbols.HGNC.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Seurat.update.gene.symbols.HGNC.R"))
# require(HGNChelper)



#' @title Update Gene Symbols in a Seurat Object
#'
#' @description This function updates gene symbols in a Seurat object based on current gene
#' nomenclature guidelines, using HGNChelper(). It checks and updates gene symbols to their
#' latest approved versions,ensuring that gene annotations are current and consistent.
#' The function optionally enforces unique gene symbols and provides statistics on the update process.
#'
#' @param obj A Seurat object containing gene expression data. Default: `ls.Seurat[[i]]`
#' (ensure to replace `i` with the actual index or variable referencing your Seurat object).
#' @param species_ The species for which the gene symbols are checked and updated,
#' used to ensure the correct gene nomenclature is applied. Default: `'human'`,
#' Supports `'human'`, `'mouse'`, as specified in the `HGNChelper` package.
#' @param EnforceUnique Logical flag indicating whether to enforce unique gene symbols
#' within the Seurat object. When set to `TRUE`, it resolves issues with duplicated gene symbols
#' by appending unique identifiers. Default: `TRUE`.
#' @param ShowStats Logical flag indicating whether to display statistics about the gene
#' symbol update process. When set to `TRUE`, it prints detailed information on the console
#' about the changes made. Default: `FALSE`.
#'
#' @return A modified Seurat object with updated gene symbols. The function directly modifies
#' the input Seurat object, ensuring that gene symbols adhere to the latest nomenclature.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Assuming `mySeuratObject` is your Seurat object
#'   updatedSeuratObject <- UpdateGenesSeurat(
#'     obj = mySeuratObject, species_ = "human",
#'     EnforceUnique = TRUE, ShowStats = TRUE
#'   )
#'   # `updatedSeuratObject` now has updated gene symbols
#' }
#' }
#'
#' @seealso
#' \code{\link[HGNChelper]{checkGeneSymbols}} for details on checking and updating gene symbols.
#'
#' @export
#' @importFrom HGNChelper checkGeneSymbols
#'
UpdateGenesSeurat <- function(obj = ls.Seurat[[i]], species_ = "human", # assay = "RNA",
                              EnforceUnique = TRUE, ShowStats = F) {
  assays.present <- Assays(obj)
  for (assay in assays.present) {
    message("Renaming in assay: ", assay, "...")

    all_genes <- Features(obj, assay = assay)

    if( species_ %in% c("human", "mouse") ) {
      HGNC.updated <- HGNChelper::checkGeneSymbols(all_genes, unmapped.as.na = FALSE, map = NULL, species = species_)
    } else {
      message(species_)
      warning("Species not supported by HGNChelper. Skipping gene symbol update.", immediate. = TRUE)
      next
    }

    if (EnforceUnique) HGNC.updated <- HGNC.EnforceUnique(HGNC.updated)

    if (ShowStats) {
      print(HGNC.updated)
      print(GetUpdateStats(HGNC.updated))
    }

    obj <- RenameGenesSeurat(obj, newnames = HGNC.updated$"Suggested.Symbol", assay = assay)
  }
  return(obj)
}


# _________________________________________________________________________________________________
#' @title Rename Gene Symbols in a Seurat Object
#'
#' @description This function replaces gene names across various slots within a specified assay
#' of a Seurat object. It is designed to be run prior to any data integration or downstream analysis
#' processes. The function targets the `@counts`, `@data`, and `@meta.features` slots within
#' the specified assay, ensuring consistency in gene nomenclature across the object.
#'
#' @param obj A Seurat object containing the assay and slots to be updated. Default: `ls.Seurat[[i]]`
#' (replace `i` with the appropriate index).
#' @param newnames A character vector containing the new gene names intended to replace the
#' existing ones. Default: `HGNC.updated[[i]]$Suggested.Symbol`. Ensure this matches the order
#' and length of the genes in the specified assay.
#' @param assay The name of the assay within the Seurat object where gene names will be updated;
#' Default: `"RNA"`. This function assumes simple objects containing only an RNA assay.
#' @param slots A character vector specifying which slots within the assay to update. Possible
#' values include `"data"`, `"counts"`, and `"meta.features"`; other layers can be specified if present.
#'
#' @details It is crucial to run this function before any data integration or further analysis
#' to ensure gene symbol consistency. The function does not support complex objects with multiple
#' assays where dependencies between assays might lead to inconsistencies. Use with caution and
#' verify the results.
#'
#' @note This function modifies the Seurat object in place, changing gene symbols directly within
#' the specified slots. Be sure to have a backup of your Seurat object if needed before applying
#' this function.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Assuming `SeuratObj` is your Seurat object
#'   # and `HGNC.updated.genes` contains the updated gene symbols
#'   SeuratObj <- RenameGenesSeurat(
#'     obj = SeuratObj,
#'     newnames = HGNC.updated.genes$Suggested.Symbol
#'   )
#'   # `SeuratObj` now has updated gene symbols in the specified assay and slots
#' }
#' }
#'
#' @export
RenameGenesSeurat <- function(obj = ls.Seurat[[i]],
                              newnames = HGNC.updated[[i]]$Suggested.Symbol,
                              assay = "RNA",
                              slots = c("data", "counts", "meta.features")) {
  #
  # browser()
  message("RenameGenesSeurat, assay: ", assay)
  message("Slots expected: ", slots)
  message("Slots found: ", SeuratObject::Layers(obj@assays[[assay]]))
  warning("Run this before integration and downstream processing. It only attempts to change
          @counts, @data, and @meta.features in obj@assays$YOUR_ASSAY.", immediate. = TRUE)


  stopifnot(
    "Unequal gene name sets: nrow(assayobj) != nrow(newnames):" =
      length(Features(obj, assay = assay)) == length(newnames)
  )

  if (obj@version < "5") {
    warning("obj@version < 5. Old versions are not supported. Update the obj!", immediate. = TRUE)
  } else {
    layers <- Layers(obj@assays[[assay]])
    slots_found <- slots %in% layers
    stopifnot(any(slots_found))
    warnifnot("Not all slots present in the object - all(slots_found)" = all(slots_found))
    slots <- intersect(slots, layers)
  }


  if ("scale.data" %in% slots) {
    n_genes_sc_dta <- nrow(obj@assays[[assay]]$"scale.data")
    stopifnot(
      "scale.data does has different number of genes than newnames!" =
        n_genes_sc_dta == length(newnames)
    )
  }

  LayersFound <- SeuratObject::Layers(obj@assays[[assay]])
  iprint("Present: ", sort(LayersFound))

  slots <- sort(intersect(slots, LayersFound))
  iprint("Replaced: ", slots)

  for (slotX in slots) {
    print(slotX)
    nrO <- nrow(SeuratObject::GetAssayData(object = obj, assay = assay, layer = slotX))
    obj <- .check_and_rename(obj, assay, newnames = newnames, layer.name = slotX)
    nrN <- nrow(SeuratObject::GetAssayData(object = obj, assay = assay, layer = slotX))
    stopifnot(nrN == nrO)
  }
  return(obj)
}


# _________________________________________________________________________________________________
#' @title Check and Rename Gene Names in Seurat Assay Object
#'
#' @description This function renames rows (genes) in a specified slot of a Seurat assay object.
#' It supports slots storing data as either a dense or a sparse matrix (dgCMatrix) or data.frame.
#'
#' @param obj A Seurat object.
#' @param assay An Assay name in a Seurat object.
#' @param newnames A character vector of new gene names to be assigned.
#' @param layer.name A string specifying the slot in the Assay object to be updated.
#'                 Valid options typically include 'counts', 'data', or 'scale.data'.
#'
#' @return An Assay object with updated gene names in the specified slot.
#' @examples
#' \dontrun{
#' # Assuming 'seurat_obj' is a Seurat object and 'new_gene_names' is a vector of gene names
#' updated_assay <- check_and_rename(
#'   assayobj = seurat_obj[["RNA"]],
#'   newnames = new_gene_names,
#'   layer.name = "counts"
#' )
#' }
.check_and_rename <- function(obj, assay, newnames, layer.name) {
  cat(layer.name, fill = TRUE)

  length_newnames <- length(newnames)
  length_orig_names <- length(Features(obj, assay = assay))

  stopifnot(
    is(obj, "Seurat"),
    is.character(assay),
    is.character(layer.name),
    is.character(newnames),
    length_orig_names == length_newnames
  )

  assayobj <- obj@assays[[assay]]
  feature.list <- rownames(assayobj@features@.Data)

  if (length(feature.list) == length(newnames)) {
    rownames(assayobj@features@.Data) <- newnames
    nrX <- length(rownames(assayobj@features@.Data))
  } else {
    iprint("length feature.list", length(feature.list), "length newnames", length(newnames))
    stop()
  }

  if (layer.name %in% SeuratObject::Layers(assayobj)) {
    matrix_n <- SeuratObject::LayerData(assayobj, layer = layer.name)
    nr1 <- nrow(matrix_n)

    if (all(dim(matrix_n)) > 0) {
      stopifnot(nrow(matrix_n) == length(newnames))

      if ("dgCMatrix" %in% class(matrix_n)) {
        message(assay, "@", layer.name, " is of type dgeCMatrix!")
        matrix_n@Dimnames[[1]] <- newnames
      } else if ("matrix" %in% class(matrix_n)) {
        message(assay, "@", layer.name, " is of type Matrix!")
        rownames(matrix_n) <- newnames
      } else if ("data.frame" %in% class(matrix_n)) {
        message(assay, "@", layer.name, " is of type data.frame!")
        rownames(matrix_n) <- newnames
      } else {
        warning(">>> No renaming: ", assay, "@", layer.name,
          " not of type dgeCMatrix / Matrix / data.frame.",
          immediate. = TRUE
        )
      }
      stopifnot(nr1 == nrow(matrix_n))

      SeuratObject::LayerData(assayobj, layer = layer.name) <- matrix_n
      nr3 <- nrow(SeuratObject::LayerData(assayobj, layer = layer.name))
      stopifnot(nr3 == nrX)
    }
  } else {
    warning(paste(">>>", assay, "@", layer.name, "does not exist!"), immediate. = TRUE)
  }
  # obj <- SetAssayData(obj, layer = layer.name, new.data = matrix_n)
  obj@assays[[assay]] <- assayobj
  return(obj)
}

# _________________________________________________________________________________________________
#' @title Remove Specific Genes from a Seurat Object
#'
#' @description Removes specified genes from the metadata, counts, data, and scale.data slots of a Seurat object.
#' This operation is typically performed prior to data integration to ensure that gene sets are consistent
#' across multiple datasets. The function modifies the Seurat object in place.
#'
#' @param obj A Seurat object. Default: `ls.Seurat[[i]]` (please ensure to replace `i` with the actual index or variable).
#' @param symbols2remove A character vector specifying the genes to be removed from the Seurat object;
#' Default: `c("TOP2A")`.
#'
#' @details This function directly modifies the `@counts`, `@data`, and `@scale.data` slots within
#' the RNA assay of the provided Seurat object, as well as the `@meta.data` slot. It's important to run
#' this function as one of the initial steps after creating the Seurat object and before proceeding
#' with downstream analyses or integration processes.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Assuming `SeuratObj` is your Seurat object and you want to remove the gene "TOP2A"
#'   updatedSeuratObj <- RemoveGenesSeurat(obj = SeuratObj, symbols2remove = "TOP2A")
#'   # Now `updatedSeuratObj` does not contain "TOP2A" in the specified slots
#' }
#' }
#'
#' @return A Seurat object with the specified genes removed from the mentioned slots.
#'
#' @export
RemoveGenesSeurat <- function(obj = ls.Seurat[[i]], symbols2remove = c("TOP2A")) {
  print("Run this as the first thing after creating the Seurat object.
        It only removes genes from: metadata; obj@assays$RNA@counts, @data and @scale.data.")
  RNA <- obj@assays$RNA

  if (length(RNA@counts)) {
    NotFound <- setdiff(symbols2remove, RNA@counts@Dimnames[[1]])
    if (length(NotFound) == 0) {
      RNA@counts@Dimnames[[1]] <- symbols2remove
      print("Genes removed from RNA@counts")
    } else {
      print("Not All Genes Found in RNA@counts. Missing:")
      print(NotFound)
    }
  }
  if (length(RNA@data)) {
    if (length(setdiff(symbols2remove, RNA@data@Dimnames[[1]])) == 0) {
      RNA@data@Dimnames[[1]] <- symbols2remove
      print("Genes removed from RNA@data.")
    } else {
      print("Not All Genes Found in RNA@data")
    }
  }
  if (length(RNA@scale.data)) {
    if (length(setdiff(symbols2remove, RNA@scale.data@Dimnames[[1]])) == 0) {
      RNA@scale.data@Dimnames[[1]] <- symbols2remove
      print("Genes removed from RNA@scale.data.")
    } else {
      print("Not All Genes Found in RNA@scale.data")
    }
  }
  if (length(obj@meta.data)) {
    if (length(setdiff(symbols2remove, rownames(obj@meta.data))) == 0) {
      rownames(obj@meta.data) <- symbols2remove
      print("Genes removed from @meta.data.")
    } else {
      print("Not All Genes Found in @metadata")
    }
  }
  obj@assays$RNA <- RNA
  return(obj)
}



# _________________________________________________________________________________________________
#' @title Enforce Unique HGNC Gene Symbols
#'
#' @description Ensures that gene symbols are unique after being updated with HGNC symbols. This function
#' applies a suffix to duplicate gene symbols to enforce uniqueness. While using `make.unique` might not
#' be the ideal solution due to potential mismatches, it significantly reduces the number of mismatching
#' genes in certain scenarios, making it a practical approach for data integration tasks.
#'
#' @param updatedSymbols A data frame or matrix containing gene symbols updated via `HGNChelper::checkGeneSymbols()`.
#' The third column should contain the updated gene symbols that are to be made unique.
#'
#' @return A modified version of the input data frame or matrix with unique gene symbols in the third column.
#' If duplicates were found, they are made unique by appending `.1`, `.2`, etc., to the repeated symbols.
#'
#' @details The function specifically targets the issue of duplicate gene symbols which can occur after
#' updating gene symbols to their latest HGNC-approved versions. Duplicate symbols can introduce
#' ambiguity in gene expression datasets, affecting downstream analyses like differential expression or
#' data integration. By ensuring each gene symbol is unique, this function helps maintain the integrity
#' of the dataset.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Assuming `SymUpd` is your data frame of updated symbols from HGNChelper::checkGeneSymbols()
#'   uniqueSymbols <- HGNC.EnforceUnique(updatedSymbols = SymUpd)
#'   # `uniqueSymbols` now contains unique gene symbols in its third column
#' }
#' }
#'
#' @note This function is a workaround for ensuring unique gene symbols and might not be suitable
#' for all datasets or analyses. It's important to review the results and ensure that the gene
#' symbols accurately represent your data.
#'
#' @export
HGNC.EnforceUnique <- function(updatedSymbols) {
  NGL <- updatedSymbols[, 3]
  if (any.duplicated(NGL)) {
    updatedSymbols[, 3] <- make.unique(NGL)
    "Unique names are enforced by suffixing .1, .2, etc."
  }
  return(updatedSymbols)
}




# _________________________________________________________________________________________________
#' @title Gene Symbol Update Statistics
#'
#' @description Generates statistics on the gene symbol updates performed by `UpdateGenesSeurat()`.
#' This function analyzes the data frame of gene symbols before and after the update process,
#' providing insights into the proportion and total number of genes that were updated.
#'
#' @param genes A data frame of gene symbols before and after update, typically the output of
#' `UpdateGenesSeurat()`. Default: `HGNC.updated[[i]]` (where `i` is the index of the desired
#' Seurat object in a list).
#'
#' @return A named vector with statistics on gene updates, including the percentage of updated genes,
#' the absolute number of updated genes, and the total number of genes processed.
#'
#' @details The function examines the `Approved` column of the input data frame to identify
#' gene symbols marked for update and compares the original and suggested symbols to determine
#' actual updates. The statistics highlight the efficiency and impact of the gene symbol
#' updating process, aiding in the assessment of data preprocessing steps.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Assuming `HGNC.updated.genes` is your data frame containing the original and
#'   # suggested gene symbols, as returned by `UpdateGenesSeurat()`
#'   updateStats <- GetUpdateStats(genes = HGNC.updated.genes)
#'   # `updateStats` now contains the update statistics, including percentage and count of updated genes
#' }
#' }
#'
#' @note The function requires the input data frame to have specific columns as produced by
#' `HGNChelper::checkGeneSymbols()` and subsequently processed by `UpdateGenesSeurat()`.
#' Ensure that the input adheres to this format for accurate statistics.
#'
#' @seealso \code{\link{UpdateGenesSeurat}}, for the function that updates gene symbols and produces
#' the input data frame for this function.
#'
#' @importFrom Stringendo percentage_formatter
#'
#' @export
GetUpdateStats <- function(genes = HGNC.updated[[i]]) {
  MarkedAsUpdated <- genes[genes$Approved == FALSE, ]
  AcutallyUpdated <- sum(MarkedAsUpdated[, 1] != MarkedAsUpdated[, 3])
  UpdateStats <- c(
    "Updated (%)" = Stringendo::percentage_formatter(AcutallyUpdated / nrow(genes)),
    "Updated Genes" = floor(AcutallyUpdated), "Total Genes" = floor(nrow(genes))
  )
  return(UpdateStats)
}


# _________________________________________________________________________________________________
#' @title PlotUpdateStats
#'
#' @description Creates a scatter plot of update statistics.
#' @param mat A matrix containing update statistics. Default: UpdateStatMat.
#' @param column.names A character vector of column names in the mat parameter. Default: c("Updated (%)", "Updated (Nr.)").
#' @return A scatter plot displaying update statistics.
#' @details This function takes a matrix containing update statistics and column names to plot
#' the corresponding statistics. It colorizes the genes and plots the percentage of total genes
#' updated against the number of genes updated.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   PlotUpdateStats(mat = result.of.GetUpdateStats)
#' }
#' }
#' @seealso
#' \code{\link[wplot]{wplot}}, \code{\link[wcolorize]{wcolorize}}
#' @importFrom MarkdownReports wplot wlegend
#'
#' @export
PlotUpdateStats <- function(mat = UpdateStatMat, column.names = c("Updated (%)", "Updated (Nr.)")) { # Scatter plot of update stats.
  stopifnot(column.names %in% colnames(UpdateStatMat))
  HGNC.UpdateStatistics <- mat[, column.names]
  HGNC.UpdateStatistics[, "Updated (%)"] <- 100 * HGNC.UpdateStatistics[, "Updated (%)"]
  colnames(HGNC.UpdateStatistics) <- c("Gene Symbols updated (% of Total Genes)", "Number of Gene Symbols updated")
  lll <- wcolorize(vector = rownames(HGNC.UpdateStatistics))
  MarkdownReports::wplot(HGNC.UpdateStatistics,
    col = lll,
    xlim = c(0, max(HGNC.UpdateStatistics[, 1])),
    ylim = c(0, max(HGNC.UpdateStatistics[, 2]))
  )
  MarkdownReports::wlegend(NamedColorVec = lll, poz = 1)
}



# _________________________________________________________________________________________________
# Handling SNP demux table results coming from SoupOrCell ______________________________ ----
# _________________________________________________________________________________________________






# _________________________________________________________________________________________________


# _________________________________________________________________________________________________
# Read.Write.Save.Load.functions.R ______________________________ ----
# _________________________________________________________________________________________________
# source('~/GitHub/Packages/Seurat.utils/Functions/Read.Write.Save.Load.functions.R')
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Read.Write.Save.Load.functions.R"))

"Multicore read / write (I/O) functions are https://github.com/vertesy/Seurat.multicore"
"Single core read / write (I/O) functions are in https://github.com/vertesy/Seurat.utils/"


# _________________________________________________________________________________________________
#' @title Convert10Xfolders
#'
#' @description This function takes a parent directory with a number of subfolders, each
#' containing the standard output of 10X Cell Ranger. It (1) loads the (filtered) data matrices,
#' (2) converts them to Seurat objects, and (3) saves them as .qs files
#'
#' @param InputDir A character string specifying the input directory.
#' @param regex A logical value. If TRUE, the folderPattern is treated as a regular expression. Default: `FALSE`.
#' @param folderPattern A character vector specifying the pattern of folder names to be searched. Default: 'filtered_feature'.
#' @param suffix A character string specifying the suffix of the files saved.
#' @param min.cells An integer value specifying the minimum number of cells. Default: 5.
#' @param min.features An integer value specifying the minimum number of features. Default: 200.
#' @param normalize_data Add normalized "data" layer?. Default: `TRUE`.
#' @param updateHGNC A logical value indicating whether to update the HGNC. Default: `TRUE`.
#' @param save Save .qs object? Default: `TRUE`.
#' @param ShowStats A logical value indicating whether to show statistics. Default: `TRUE`.
#' @param writeCBCtable A logical value indicating whether to write out a list of cell barcodes (CBC) as a tsv file. Default: `TRUE`.
#' @param depth An integer value specifying the depth of scan (i.e., how many levels below the InputDir). Default: 2.
#' @param sort_alphanumeric sort files alphanumeric? Default: `TRUE`.
#' @param save_empty_droplets save empty droplets? Default: `TRUE`.
#'
#' @examples
#' \dontrun{
#' if (interactive()) Convert10Xfolders(InputDir)
#' }
#' @export
Convert10Xfolders <- function(
    InputDir,
    regex = FALSE,
    folderPattern = c("filtered_feature", "raw_feature", "SoupX_decont")[1],
    suffix = strsplit(folderPattern, "_")[[1]][1],
    depth = 4,
    min.cells = 5, min.features = 200,
    normalize_data = TRUE,
    updateHGNC = TRUE, ShowStats = TRUE,
    writeCBCtable = TRUE,
    nthreads = .getNrCores(),
    preset = "high",
    ext = "qs",
    sort_alphanumeric = TRUE,
    save_empty_droplets = TRUE,
    ...) {
  stopifnot(
    is.character(InputDir), dir.exists(InputDir),
    is.logical(regex), is.character(folderPattern), is.character(suffix), is.numeric(depth),
    is.numeric(min.cells), is.numeric(min.features), is.logical(updateHGNC), is.logical(ShowStats), is.logical(writeCBCtable),
    is.logical(sort_alphanumeric)
  )

  finOrig <- ReplaceRepeatedSlashes(list.dirs.depth.n(InputDir, depth = depth))
  fin <- CodeAndRoll2::grepv(x = finOrig, pattern = folderPattern, perl = regex)

  message(length(fin), " samples found.")

  samples <- basename(list.dirs(InputDir, recursive = FALSE))
  if (sort_alphanumeric) samples <- gtools::mixedsort(samples)
  iprint("Samples:", samples)

  if (!length(fin) > 0) {
    stop(paste("No subfolders found with pattern", folderPattern, "in dirs like: ", finOrig[1:3]))
  }

  for (i in 1:length(fin)) {
    print(i)
    pathIN <- Stringendo::FixPath(fin[i])
    message(pathIN)
    fnameIN <- basename(dirname(dirname(pathIN)))
    message(fnameIN)

    count_matrix <- Read10X(pathIN)
    if (!is.list(count_matrix) | length(count_matrix) == 1) {
      seu <- CreateSeuratObject(
        counts = count_matrix, project = fnameIN,
        min.cells = min.cells, min.features = min.features
      )
    } else {
      (stop("length(count_matrix) != 1"))
    }

    ncells <- ncol(seu)
    fname_X <- Stringendo::sppp(
      fnameIN, suffix, "min.cells", min.cells, "min.features", min.features,
      "cells", ncells
    )
    print(fname_X)

    f.path.out <- Stringendo::ParseFullFilePath(path = InputDir, file_name = fname_X, extension = ext)
    message(f.path.out)

    # update --- --- ---
    if (updateHGNC) seu <- UpdateGenesSeurat(seu, EnforceUnique = TRUE, ShowStats = TRUE)

    # NormalizeData --- --- ---
    if (normalize_data) seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000, verbose = TRUE)

    # write out --- --- ---
    if (save) qs::qsave(x = seu, file = f.path.out, nthreads = nthreads, preset = preset)

    # write cellIDs ---  --- ---
    if (writeCBCtable) {
      CBCs <- t(t(colnames(seu)))
      colnames(CBCs) <- "CBC"
      ReadWriter::write.simple.tsv(input_df = CBCs, manual_file_name = sppp(fnameIN, suffix, "CBC"), manual_directory = InputDir)
    }

    if (save_empty_droplets & suffix == "raw") {
      # Select and save empty droplets (the Soup)

      path_filtered <- gsub(x = pathIN, pattern = "/raw_feature_", replacement = "/filtered_feature_")
      fnp_filtered <- spps(path_filtered, "barcodes.tsv.gz")


      if (file.exists(fnp_filtered)) {
        SoupDir <- spps(InputDir, "Soup")
        dir.create(SoupDir)

        CBCs_HQ <- read.simple.vec(fnp_filtered)

        CBC_empty_drops <- setdiff(colnames(seu), CBCs_HQ)
        nr.empty.droplets <- length(CBC_empty_drops)
        umi_per_CBC <- colSums(seu@assays$RNA@layers$counts)
        pct.empty.droplets.max10umis <- pc_TRUE(umi_per_CBC < 11)
        message("We have ", nr.empty.droplets, " empty droplets, ", pct.empty.droplets.max10umis, " of which have max 10 umis.")
        FNM <- sppp("nr.empty.droplets", fnameIN, nr.empty.droplets)
        ReadWriter::write.simple.vec(nr.empty.droplets, manual_file_name = FNM, manual_directory = SoupDir)

        obj_empty_drops <- subset(seu, cells = CBC_empty_drops)

        f_path_out_ED <- Stringendo::ParseFullFilePath(path = SoupDir, file_name = sppp("obj.empty.droplets", fnameIN, nr.empty.droplets), extension = ext)
        qs::qsave(x = obj_empty_drops, file = f_path_out_ED, nthreads = nthreads, preset = preset)

        # save the bulk RNA counts of the empty droplets
        Soup.Bulk.RNA <- rowSums(count_matrix[, CBC_empty_drops])
        f_path_out_Bulk <- Stringendo::ParseFullFilePath(path = SoupDir, file_name = sppp("Soup.Bulk.RNA", fnameIN), extension = "qs")
        qs::qsave(x = Soup.Bulk.RNA, file = f_path_out_Bulk, nthreads = nthreads, preset = preset)
        ReadWriter::write.simple.tsv(Soup.Bulk.RNA, suffix = fnameIN, manual_directory = SoupDir)
      }
    } else {
      message("No empty droplets saved. suffix ", suffix)
    }
  } # for
}



# _________________________________________________________________________________________________
#' @title ConvertDropSeqfolders
#'
#' @description This function takes a parent directory with a number of subfolders, each
#' containing the standard output of 10X Cell Ranger. It (1) loads the filtered data matrices,
#' (2) converts them to Seurat objects, and (3) saves them as .RDS files.
#' @param InputDir A character string specifying the input directory.
#' @param folderPattern A character string specifying the pattern of folder names to be searched. Default: 'SRR*'.
#' @param filePattern A character string specifying the pattern of file names to be searched. Default: 'expression.tsv.gz'.
#' @param useVroom A logical value indicating whether to use vroom. Default: `TRUE`.
#' @param col_types.vroom A list defining column types for vroom. Default: list("GENE" = "c", .default = "d").
#' @param min.cells An integer value specifying the minimum number of cells. Default: 10.
#' @param min.features An integer value specifying the minimum number of features. Default: 200.
#' @param updateHGNC A logical value indicating whether to update the HGNC. Default: `TRUE`.
#' @param ShowStats A logical value indicating whether to show statistics. Default: `TRUE`.
#' @param minDimension An integer value specifying the minimum dimension. Default: 10.
#' @param overwrite A logical value indicating whether to overwrite files. Default: `FALSE`.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   ConvertDropSeqfolders(InputDir = InputDir)
#' }
#' }
#' @seealso
#'  \code{\link[vroom]{vroom}}
#'  \code{\link[readr]{read_delim}}
#' @importFrom readr read_tsv
#'
#' @export
ConvertDropSeqfolders <- function(
    InputDir,
    folderPattern = "SRR*", filePattern = "expression.tsv.gz",
    useVroom = TRUE, col_types.vroom = list("GENE" = "c", .default = "d"),
    min.cells = 10, min.features = 200, updateHGNC = TRUE, ShowStats = TRUE, minDimension = 10, overwrite = FALSE) {
  InputDir <- FixPath(InputDir)
  fin <- list.dirs(InputDir, recursive = FALSE)
  fin <- CodeAndRoll2::grepv(x = fin, pattern = folderPattern, perl = FALSE)

  for (i in 1:length(fin)) {
    print(i)
    pathIN <- FixPath(fin[i])
    print(pathIN)
    fnameIN <- basename(fin[i])
    subdir <- paste0(InputDir, fnameIN)
    fnameOUT <- ppp(subdir, "min.cells", min.cells, "min.features", min.features, "Rds")
    print(fnameOUT)
    if (!overwrite) {
      OutFile <- list.files(InputDir, pattern = basename(fnameOUT), recursive = TRUE)
      if (length(OutFile) > 0) {
        if (grepl(pattern = ".Rds$", OutFile, perl = TRUE)) {
          iprint("      RDS OBJECT ALREADY EXISTS.")
          next
        }
      } # if length
    }
    CountTable <- list.files(subdir, pattern = filePattern, recursive = FALSE)
    stopifnot(length(CountTable) == 1)
    count_matrix <- if (useVroom) {
      stopifnot("Package 'vroom' must be installed to use this function." = require("vroom"))
      vroom::vroom(file = kpps(subdir, CountTable), col_types = col_types.vroom)
    } else {
      readr::read_tsv(file = kpps(subdir, CountTable))
    }

    if (nrow(count_matrix) < minDimension | ncol(count_matrix) < minDimension) {
      iprint("")
      iprint("      EXPRESSION MATRIX TOO SMALL.", nrow(count_matrix), "x", ncol(count_matrix), ". Not processed.")
    } else {
      count_matrix <- FirstCol2RowNames(count_matrix)[, -1] # remove 1st "Cell column" # https://github.com/vertesy/SEO/issues/63
      seu <- CreateSeuratObject(
        counts = count_matrix, project = fnameIN,
        min.cells = min.cells, min.features = min.features
      )
      if (ncol(seu) < 1000) print("Only", ncol(seu), "cells survived filtering in the Seurat obj!")
      if (nrow(seu) < 1000) print("Only", nrow(seu), "genes found in the Seurat obj!")

      # update HGNC --- --- --- --- ---
      Sys.setenv("R_MAX_VSIZE" = 32000000000)
      if (updateHGNC) seu <- UpdateGenesSeurat(seu, EnforceUnique = TRUE, ShowStats = TRUE)
      saveRDS(seu, file = fnameOUT)
    }
  }
}


# _________________________________________________________________________________________________
#' @title LoadAllSeurats
#'
#' @description This function loads all Seurat objects found in a directory. It also works with
#' symbolic links (but not with aliases).
#' @param InputDir A character string specifying the input directory.
#' @param file.pattern A character string specifying the pattern of file names to be searched.
#' Default: '^filtered.+Rds$'.
#' @param string.remove1 A character string or FALSE. If a string is provided, it is removed from
#' file names. Default: "filtered_feature_bc_matrix.".
#' @param string.replace1 A character string of the new text instead of "string.remove1".
#' @param string.remove2 A character string or FALSE. If a string is provided, it is removed from
#' file names. Default: ".min.cells.10.min.features.200.Rds".
#' @param sort_alphanumeric sort files alphanumeric? Default: `TRUE`.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   ls.Seurat <- LoadAllSeurats(InputDir = InputDir)
#' }
#' }
#' @export
#' @importFrom tictoc tic toc
LoadAllSeurats <- function(
    InputDir,
    file.pattern = "^filtered.+Rds$",
    string.remove1 = list(FALSE, "filtered_feature_bc_matrix.", "raw_feature_bc_matrix.")[[2]],
    string.replace1 = "",
    string.remove2 = list(FALSE, ".min.cells.10.min.features.200.Rds")[[2]],
    sort_alphanumeric = TRUE) {
  tictoc::tic("LoadAllSeurats")
  InputDir <- FixPath(InputDir)

  print(file.pattern)
  use_rds <- grepl(pattern = "Rds", x = file.pattern) && !grepl(pattern = "qs", x = file.pattern)
  print(use_rds)

  fin.orig <- list.files(InputDir, include.dirs = FALSE, pattern = file.pattern)
  print(fin.orig)
  print(length(fin.orig))
  stopifnot(length(fin.orig) > 0)
  fin <- if (!isFALSE(string.remove1)) sapply(fin.orig, gsub, pattern = string.remove1, replacement = string.replace1) else fin.orig
  fin <- if (!isFALSE(string.remove2)) sapply(fin, gsub, pattern = string.remove2, replacement = "") else fin
  if (sort_alphanumeric) fin <- gtools::mixedsort(fin)


  ls.Seu <- list.fromNames(fin)
  for (i in 1:length(fin)) {
    print(fin[i])
    FNP <- paste0(InputDir, fin.orig[i])
    # print(paste("Attempting to load file:", FNP))  # Debug print

    if (use_rds) {
      ls.Seu[[i]] <- readRDS(FNP)
    } else if (!use_rds) {
      ls.Seu[[i]] <- qs::qread(file = FNP)
    } else {
      warning("File pattern ambiguous. Use either qs or rds:", file.pattern, immediate. = TRUE)
    }
  } # for
  print(tictoc::toc())
  return(ls.Seu)
}




# _________________________________________________________________________________________________
#' @title Load 10X Genomics Data as Seurat Object
#'
#' @description Reads 10X Genomics dataset files (gzipped) including matrix, features, and barcodes,
#' to a single expression matrix. This function handles the unzipping of these files, reads the data,
#' and re-compresses the files back to their original gzipped format.
#'
#' @param dir A character string specifying the path to the directory containing the 10X dataset files.
#' This directory should contain `matrix.mtx.gz`, `features.tsv.gz`, and `barcodes.tsv.gz` files.
#'
#' @return A Seurat object containing the single-cell RNA-seq data extracted from the provided 10X
#' Genomics dataset.
#'
#' @details This function facilitates the loading of 10X Genomics datasets into R for analysis with
#' the Seurat package. It specifically caters to gzipped versions of the `matrix.mtx`, `features.tsv`,
#' and `barcodes.tsv` files, automating their decompression, reading, and subsequent recompression.
#' The function relies on Seurat's `Read10X` function for data reading and object construction.
#'
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   # Replace `path_to_10x_data` with the path to your 10X data directory
#'   seuratObject <- read10x(dir = "path_to_10x_data")
#'   # `seuratObject` is now a Seurat object containing the loaded 10X data
#' }
#' }
#'
#' @note Ensure that the specified directory contains the required gzipped files.
#' If the `features.tsv.gz` file is named differently (e.g., `genes.tsv.gz`), please rename it
#' accordingly before running this function.
#'
#' @seealso \code{\link[Seurat]{Read10X}} for the underlying function used to read the 10X data.
#'
#' @importFrom tictoc tic toc
#' @importFrom R.utils gunzip gzip
#' @importFrom Seurat Read10X
#'
#' @export
read10x <- function(dir) {
  tictoc::tic("read10x")
  names <- c("barcodes.tsv", "features.tsv", "matrix.mtx")
  for (i in 1:length(names)) {
    R.utils::gunzip(paste0(dir, "/", names[i], ".gz"))
  }
  file.copy(paste0(dir, "/features.tsv"), paste0(dir, "/genes.tsv"))
  mat <- Seurat::Read10X(dir)
  file.remove(paste0(dir, "/genes.tsv"))
  for (i in 1:length(names)) {
    R.utils::gzip(paste0(dir, "/", names[i]))
  }
  tictoc::toc()
  mat
}



# _________________________________________________________________________________________________
#' @title .saveRDS.compress.in.BG
#'
#' @description Save and RDS object and compress resulting file in the background using system(gzip). OS X or unix.
#' @param obj Seurat object.
#' @param compress_internally Compress by R? Default: `FALSE`. (still compressed in background via CLI).
#' @param compr Compress at all? Default: `TRUE`.
#' @param fname File name
#' @param ... Additional parameters passed to saveRDS() function.
#' @seealso
#'  \code{\link[tictoc]{tic}}
#' @importFrom tictoc tic toc
.saveRDS.compress.in.BG <- function(obj, compr = FALSE, fname, compress_internally = FALSE, ...) {
  try(tictoc::tic(".saveRDS.compress.in.BG"), silent = TRUE)
  saveRDS(object = obj, compress = compress_internally, file = fname, ...)
  try(tictoc::toc(), silent = TRUE)
  if (compr) system(command = paste0("gzip '", fname, "'"), wait = FALSE) # execute in the background
  print(paste("Saved, optionally being .gz compressed", fname))
  try(say(), silent = TRUE)
}




# _________________________________________________________________________________________________
#' @title isave.RDS
#'
#' @description Save an RDS object, using a faster and efficient compression method that runs in the background.
#' @param obj The object to be saved, typically a Seurat object.
#' @param prefix A string prefix added to the filename. Default: NULL.
#' @param suffix A string suffix added to the filename. Default: NULL.
#' @param inOutDir A boolean flag, if TRUE the OutDir is used as save directory, if FALSE the
#' alternative_path_rdata is used. Default: `TRUE`.
#' @param project A string representing the project code. This is appended to the saved file name.
#' Default: the active project determined by getProject().
#' @param alternative_path_rdata A string that specifies the alternative path for storing the
#' RDS file if inOutDir is FALSE. Default: "~/Dropbox (VBC)/Abel.IMBA/AnalysisD/_RDS.files/"
#' appended with the basename of OutDir.
#' @param homepath A string representing the homepath. Will be replaced by '~' in the file path. Default: '~/'.
#' @param showMemObject A boolean flag, if TRUE the function will print out the memory size of the
#' largest objects in the workspace. Default: `TRUE`.
#' @param saveParams A boolean flag, if TRUE the parameters 'p' and 'all.genes' are added to the
#' 'misc' slot of the Seurat object if the object is of class Seurat. Default: `TRUE`.
#' @param compress Compress .Rds file after writing? Default: `TRUE`.
#' @param test_read Provide command to test validity by reading in the object just written.
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   isave.RDS(my.R.object)
#' }
#' }
#' @export
isave.RDS <- function(
    obj, prefix = NULL, suffix = NULL, inOutDir = TRUE,
    project = getProject(),
    alternative_path_rdata = paste0("~/Dropbox (VBC)/Abel.IMBA/AnalysisD/_RDS.files/", basename(OutDir)),
    homepath = if (Sys.info()[1] == "Darwin") "~/" else "/users/abel.vertesy/",
    showMemObject = TRUE, saveParams = TRUE,
    compress = TRUE,
    test_read = FALSE) {
  warning("isave.RDS() is deprecated. Use xsave() to save in .qs format.", immediate. = TRUE)
  path_rdata <- if (inOutDir) OutDir else alternative_path_rdata
  dir.create(path_rdata)

  if (showMemObject) {
    try(memory.biggest.objects(), silent = TRUE)
  }
  if ("Seurat" %in% is(obj) & saveParams) {
    try(obj@misc$p <- p, silent = TRUE)
    try(obj@misc$all.genes <- all.genes, silent = TRUE)
  }
  fnameBase <- kppu(prefix, substitute_deparse(obj), project, suffix, idate(Format = "%Y.%m.%d_%H.%M"))
  fnameBase <- trimws(fnameBase, whitespace = "_")
  FNN <- paste0(path_rdata, fnameBase, ".Rds")
  FNN <- gsub(pattern = "~/", replacement = homepath, x = FNN)
  print(FNN)
  if (test_read) {
    print(paste0('xx5 <- read_rds(\\"', FNN, '\\")'))
  } else {
    Seurat.utils:::.saveRDS.compress.in.BG(obj = obj, fname = FNN, compr = compress, compress_internally = FALSE)
  }
}

# _________________________________________________________________________________________________
#' @title Save an R Object Using 'qs' Package for Fast Compressed Saving
#'
#' @description This function saves an R object to a file in a quick and efficient format using the 'qs' package.
#' It constructs the file name based on various inputs and stores additional metadata if the object is a Seurat object.
#' The saving path can be adjusted by the presence of 'OutDir' in the global environment or defaults to the working directory.
#'
#' @param obj The R object to be saved.
#' @param suffix Optional; a suffix to add to the filename.
#' @param prefix Optional; a prefix to add to the filename.
#' @param nthreads Number of threads to use when saving, defaults to 12.
#' @param preset Compression preset, defaults to 'high'.
#' @param project The project name to be included in the filename, defaults to the result of `getProject()`.
#' @param dir Output Directory
#' @param showMemObject Logical; if TRUE, displays the memory size of the largest objects.
#' @param saveParams Logical; if TRUE and if the object is a Seurat object, additional parameters
#' are saved within it.
#' @param paramList Optional; a list of parameters to save within the Seurat object.
#' @param allGenes Optional; a list of all genes to save within the Seurat object.
#' @param saveLocation Logical; if TRUE and if the object is a Seurat object, file location is saved
#' into misc slot.
# #' @param backgroundJob NOT IMPLEMENTED. Logical; if TRUE, the compression is done in the background.
#' @param v Verbose output.
#'
#' @return Invisible; The function is called for its side effects (saving a file) and does not return anything.
#'
#' @note The function uses the 'qs' package for quick and efficient serialization of objects and
#' includes a timing feature from the 'tictoc' package.
#' @seealso \code{\link[qs]{qsave}} for the underlying save function used.
#' @importFrom qs qsave
#' @importFrom tictoc tic toc
#' @importFrom rstudioapi isAvailable
#'
#' @export
xsave <- function(
    obj,
    suffix = NULL,
    prefix = NULL,
    nthreads = if (object.size(obj) < 1e7) 1 else .getNrCores(12),
    preset = "high",
    project = getProject(),
    dir = if (exists("OutDir")) OutDir else getwd(),
    showMemObject = TRUE,
    saveParams = if (exists("p")) TRUE else FALSE, # save allGenes and paramList
    paramList = if (exists("p")) p else NULL,
    allGenes = if (exists("all.genes")) all.genes else NULL,
    saveLocation = TRUE,
    # backgroundJob = FALSE,
    v = TRUE) {
  #
  if (v) message(nthreads, " threads.\n-----------")
  if (v) message("project: ", project)

  # check if the object is a Seurat object
  obj_is_seurat <- inherits(obj, "Seurat")
  if (obj_is_seurat) {
    annot.suffix <- kpp(ncol(obj), "cells")
  } else {
    saveParams <- FALSE
    annot.suffix <- if (is.list(obj)) kppd("ls", length(obj)) else NULL
  }

  if (!isFALSE(saveParams)) message("paramList: ", if (exists("paramList")) paste(substitute_deparse(paramList), length(paramList), " elements.") else " not provided.")
  if (!isFALSE(saveParams)) message("allGenes: ", if (exists("allGenes")) " found as global variable." else " not provided.")

  try(tictoc::tic("xsave"), silent = TRUE)
  if (showMemObject & v) try(memory.biggest.objects(), silent = TRUE)

  fnameBase <- trimws(kppu(
    prefix, as.character(substitute(obj)), annot.suffix, suffix, project,
    idate(Format = "%Y.%m.%d_%H.%M")
  ), whitespace = "_")

  FNN <- paste0(dir, fnameBase, ".qs")
  CMND <- paste0(substitute(obj), " <- xread('", FNN, "')")
  if (v) message(CMND)

  if ("Seurat" %in% is(obj)) {
    if (saveParams) {
      if (exists("paramList")) try(obj@misc$"p" <- paramList, silent = TRUE)
      if (exists("allGenes")) try(obj@misc$"all.genes" <- allGenes, silent = TRUE)
    }
    if (saveLocation) try(obj@misc$"file.location" <- CMND, silent = TRUE)
  }

  qs::qsave(x = obj, file = FNN, nthreads = nthreads, preset = preset)

  try(tictoc::toc(), silent = TRUE)
}

# _________________________________________________________________________________________________
#' @title Read an R Object Using 'qs' Package for Fast Decompression
#'
#' @description This function reads an R object from a file saved in a format specific to the 'qs' package,
#' which is designed for quick and efficient compression and decompression of R objects.
#' It also times the read operation, providing feedback on the duration of the operation.
#'
#' @param file A character string specifying the path to the file where the R object is saved.
#' @param nthreads The number of threads to use when reading the object, defaults to 4.
#' @param loadParamsAndAllGenes Logical; if TRUE and if the object is a Seurat object, additional parameters
#' are loaded from within it.
#' @param overwriteParams Logical; if TRUE and if the object is a Seurat object, the parameters are overwritten.
#' @param overwriteAllGenes Logical; if TRUE and if the object is a Seurat object, the all genes are overwritten.
#' @param set_m Logical; if TRUE, the variable 'm', a list of @meta.data colnames, is assigned to
#' the global environment.
#' @param ... Further arguments passed on to the 'qs::qread' function.
#'
#' @return The R object that was saved in the specified file.
#' @note The function uses the 'qs' package for fast and efficient deserialization of objects
#' and includes a timing feature from the 'tictoc' package.
#'
#' @seealso \code{\link[qs]{qread}} for the underlying read function used.
#' @importFrom qs qread
#' @importFrom tictoc tic toc
#' @importFrom rstudioapi isAvailable
#'
#' @export
xread <- function(file,
                  nthreads = if (file.size(file) < 1e7) 1 else 4,
                  loadParamsAndAllGenes = TRUE,
                  overwriteParams = FALSE,
                  overwriteAllGenes = FALSE,
                  set_m = TRUE,
                  ...) {
  stopifnot(file.exists(file))

  message(nthreads, " threads.")
  try(tictoc::tic("xread"), silent = TRUE)

  obj <- qs::qread(file = file, nthreads = nthreads, ...)

  report <- if (is(obj, "Seurat")) {
    kppws("with", ncol(obj), "cells &", ncol(obj@meta.data), "meta columns.")
  } else if (is.list(obj)) {
    kppws("is a list of:", length(obj))
  } else {
    kppws("of length:", length(obj))
  }


  if ("Seurat" %in% is(obj)) {
    if (loadParamsAndAllGenes) {
      p_local <- obj@misc$"p"
      all.genes_local <- obj@misc$"all.genes"

      if (is.null(p_local)) {
        message("No parameter list 'p' found in object@misc.")
      } else {
        recall.parameters(obj = obj, overwrite = overwriteParams)
      }

      if (is.null(all.genes_local)) {
        message("No gene list 'all.genes' found in object@misc.")
      } else {
        recall.all.genes(obj = obj, overwrite = overwriteAllGenes)
      }
    } # loadParamsAndAllGenes

    if (set_m) {
      # if (!exists("m")) {
      # m <- list.fromNames(colnames(obj@meta.data))
      m <- lapply(data.frame(obj@meta.data), function(x) head(unique(x), 50))
      assign("m", m, envir = .GlobalEnv)
      message("Variable 'm', a list of @meta.data colnames and first 50 uq values, is now defined in the global environment.")
      # } else {
      #   message("Variable 'm' already exists in the global environment, not overwritten")
      # } # exists("m")
    } # set_m
  } # Seurat


  iprint(is(obj)[1], report)
  try(tictoc::toc(), silent = TRUE)
  invisible(obj)
}

# _________________________________________________________________________________________________
#' @title Load a .qs object with optional SLURM-based safe-memory check (CBE only)
#'
#' @description
#' Loads a `.qs` serialized object (e.g., Seurat object or list) with an optional memory-safety
#' check that prevents loading objects larger than the SLURM job's allocated memory. The memory
#' check runs **only** when:
#' 1) `safe_load = TRUE`
#' 2) the global variable `onCBE` exists and is `TRUE`
#' 3) the session is running inside a SLURM job with a defined memory limit
#'
#' If no SLURM memory limit is detected, the safety check is skipped and a warning is shown.
#'
#' @param path Path to the `.qs` file to load. Must exist. Default: none.
#' @param nthreads Number of threads for `qs::qread()`. Uses 1 if file < 1e7 bytes, else 4.
#'   Default: `if (file.size(path) < 1e7) 1 else 4`.
#' @param loadParamsAndAllGenes Logical; if `TRUE`, recall stored parameters and gene lists
#'   from `obj@misc`. Default: `TRUE`.
#' @param overwriteParams Logical; overwrite existing parameters when recalling. Default: `FALSE`.
#' @param overwriteAllGenes Logical; overwrite existing `all.genes` list. Default: `FALSE`.
#' @param set_m Logical; if `TRUE`, create variable `m` in global environment with metadata values.
#'   Default: `TRUE`.
#' @param safe_load Logical; enable SLURM-based memory safety check. Default: `TRUE`.
#' @param disk2mem_size_inflation Estimated expansion factor of `.qs` file once in memory.
#'   Default: `3`.
#' @param ... Additional arguments passed to `qs::qread()`. Default: none.
#'
#' @return Invisibly returns the loaded object.
#'
#' @examples
#' \dontrun{
#' obj <- xread2("/path/to/object.qs")
#' }
#'
#' @importFrom qs qread
#' @importFrom tictoc tic toc
#' @importFrom Stringendo ifExistsAndTrue
#' @export
xread2 <- function(path,
                   nthreads = if (file.size(path) < 1e7) 1 else 4,
                   loadParamsAndAllGenes = TRUE,
                   overwriteParams = FALSE,
                   overwriteAllGenes = FALSE,
                   set_m = TRUE,
                   safe_load = TRUE,
                   disk2mem_size_inflation = 3,
                   ...) {
  stopifnot(file.exists(path))

  # Pretty print bytes
  bytes <- function(x) format(structure(x, class = "object_size"), units = "auto")

  # Current R memory usage (RSS)
  get_rss <- function() {
    kb <- as.numeric(gsub(
      "\\D", "",
      grep("^VmRSS:", readLines("/proc/self/status"),
        value = TRUE
      )
    ))
    kb * 1024
  }

  # SLURM memory (in bytes), or NA if not inside SLURM job
  get_slurm_limit <- function() {
    m <- suppressWarnings(as.numeric(Sys.getenv("SLURM_MEM_PER_NODE", "")))
    if (!is.na(m) && m > 0) {
      return(m * 1024^2)
    }
    NA_real_
  }

  # ---- MEMORY SAFETY CHECK ----
  if (safe_load && Stringendo::ifExistsAndTrue("onCBE")) {
    lim <- get_slurm_limit()

    if (is.na(lim)) {
      warning("Memory safety check skipped: no SLURM memory limit detected.")
    } else {
      sz <- file.size(path)
      rss <- get_rss()

      # Leave 10% or 3GB free, whichever is larger
      safe_cap <- min(0.9 * lim, lim - 3 * 1024^3)
      safe_cap <- max(safe_cap, rss + 1 * 1024^3) # never below current usage
      headroom <- safe_cap - rss
      need <- disk2mem_size_inflation * sz

      message(
        "Need: ", bytes(need),
        " | Headroom: ", bytes(headroom),
        " | RSS (used memory): ", bytes(rss)
      )

      if (need > headroom) {
        message("⚠ Estimated need exceeds SLURM job memory. Loading may crash the session.")
        ans <- readline("Do you still want to load the file? [y/N]: ")
        if (tolower(ans) != "y") {
          return(invisible(NULL))
        }
      }
    }
  }

  # ---- ORIGINAL LOGIC BELOW ----
  message(nthreads, " threads.")
  try(tictoc::tic("xread2"), silent = TRUE)

  obj <- qs::qread(file = path, nthreads = nthreads, ...)

  report <- if (is(obj, "Seurat")) {
    kppws("with", ncol(obj), "cells &", ncol(obj@meta.data), "metadata columns.")
  } else if (is.list(obj)) {
    kppws("is a list of:", length(obj))
  } else {
    kppws("of length:", length(obj))
  }

  if ("Seurat" %in% is(obj)) {
    if (loadParamsAndAllGenes) {
      p_local <- obj@misc$"p"
      all.genes_local <- obj@misc$"all.genes"

      if (is.null(p_local)) {
        message("No parameter list 'p' found in object@misc.")
      } else {
        recall.parameters(obj = obj, overwrite = overwriteParams)
      }

      if (is.null(all.genes_local)) {
        message("No gene list 'all.genes' found in object@misc.")
      } else {
        recall.all.genes(obj = obj, overwrite = overwriteAllGenes)
      }
    }

    if (set_m) {
      m <- lapply(data.frame(obj@meta.data), function(x) head(unique(x), 50))
      assign("m", m, envir = .GlobalEnv)
      message("Variable 'm' created in global environment with metadata preview.")
    }
  }

  iprint(is(obj)[1], report)
  try(tictoc::toc(), silent = TRUE)
  invisible(obj)
}



# _________________________________________________________________________________________________
# Save workspace
# requires MarkdownReports (github) and defining OutDir
# requires github/vertesy/CodeAndRoll.r

#' @title isave.image
#'
#' @description Save an image of the current workspace using a faster and efficient compression
#' method that runs in the background.
#' @param ... Additional parameters passed to the idate() function in the creation of the file name.
#' @param path_rdata A string that specifies the path for storing the image of the workspace.
#' Default: "~/Dropbox/Abel.IMBA/AnalysisD/_Rdata.files/" appended with the basename of OutDir.
#' @param showMemObject A boolean flag, if TRUE the function will print out the memory size of the
#' largest objects in the workspace. Default: `TRUE`.
#' @param options A string for gzip options. Default: "--force".
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   isave.image(my.R.image)
#' }
#' }
#' @export
#' @importFrom Stringendo kollapse iprint
isave.image <- function(
    ..., path_rdata = paste0("~/Dropbox/Abel.IMBA/AnalysisD/_Rdata.files/", basename(OutDir)),
    showMemObject = TRUE, options = c("--force", NULL)[1]) {
  dir.create(path_rdata)

  if (showMemObject) {
    try(memory.biggest.objects(), silent = TRUE)
  }
  fname <- Stringendo::kollapse(path_rdata, "/", idate(), ..., ".Rdata")
  print(fname)
  if (nchar(fname) > 2000) stop()
  save.image(file = fname, compress = FALSE)
  iprint("Saved, being compressed", fname)
  system(paste("gzip", options, fname), wait = FALSE) # execute in the background
}


# _________________________________________________________________________________________________
#' @title Save workspace - qsave.image
#'
#' @description Save the workspace with external gzip compression (CPU intensive).
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @param showMemObject Logical; if TRUE, the function will print out the memory size of the largest
#' objects in the workspace. Default: `TRUE`.
#' @param options Options passed on to gzip, via CLI. Default: `c("--force", NULL)[1]`
#' @seealso
#'  \code{\link[Stringendo]{kollapse}}, \code{\link[function]{iprint}}
#' @export
#' @importFrom Stringendo kollapse iprint
#' @importFrom tictoc tic toc
qsave.image <- function(
    ..., showMemObject = TRUE,
    options = c("--force", NULL)[1]) {
  tictoc::tic("qsave.image")

  fname <- Stringendo::kollapse(getwd(), "/", basename(OutDir), idate(), ..., ".Rdata")
  print(fname)
  if (nchar(fname) > 2000) stop()
  save.image(file = fname, compress = FALSE)
  iprint("Saved, being compressed", fname)
  system(paste("gzip", options, fname), wait = FALSE) # execute in the background
  cat(tictoc::toc)
}


# _________________________________________________________________________________________________
#' @title Find 'Outs' Subdirectories in Specified Subdirectories
#'
#' @description This function searches through specified subdirectories within a root directory
#' to find all subdirectories named 'outs' and returns a character vector with their full paths.
#'
#' @param root_dir The root directory.
#' @param subdir A character vector of subdirectory names within the root directory to be scanned.
#' @param recursive Boolean indicating whether to search recursively within subdirectories.
#' @return A character vector containing the full paths to the 'outs' subdirectories.
#' @importFrom fs dir_ls
#' @export
find10XoutputFolders <- function(root_dir, subdir, recursive = TRUE) {
  stopifnot(
    is.character(root_dir), length(root_dir) == 1, dir.exists(root_dir),
    is.character(subdir), all(dir.exists(file.path(root_dir, subdir))),
    is.logical(recursive)
  )

  outs_dirs <- c()
  for (i in seq_along(subdir)) {
    path <- file.path(root_dir, subdir[i])
    printProgress(i, length(subdir), "Processing subdirectory")

    iprint("Searching in:", path)
    found_dirs <- fs::dir_ls(path, recurse = recursive, glob = "*/outs", type = "directory")
    iprint(length(found_dirs), "output folders found.")
    outs_dirs <- c(outs_dirs, found_dirs)
  }

  # Replace root_dir in the paths with an empty string for printing
  outs_print <- gsub(paste0("^", root_dir, "/?"), "", outs_dirs)
  iprint(length(outs_dirs), outs_print)

  return(outs_dirs)
}


# _________________________________________________________________________________________________
#' @title Clip Suffixes from 10X Cell Names
#'
#' @description Removes suffixes from cell names that are added by 10X technology and Seurat during data processing.
#'
#' @param cellnames A vector of cell names with potential suffixes.
#' @return A vector of cell names with suffixes removed.
#' @examples
#' cellnames <- c("cell1_1", "cell2_2")
#' clip10Xcellname(cellnames)
#' @export
#' @importFrom stringr str_split_fixed
clip10Xcellname <- function(cellnames) {
  stringr::str_split_fixed(cellnames, "_", n = 2)[, 1]
}

# _________________________________________________________________________________________________
#' @title Add Suffix to Cell Names (e.g. lane suffix: _1)
#'
#' @description Appends a specified suffix to cell names to mimic lane suffixes used in 10X datasets.
#'
#' @param cellnames A vector of cell names without numeric suffixes.
#' @param suffix The suffix to add to each cell name. Default: '_1'.
#' @return A vector of cell names with the specified suffix appended.
#' @examples
#' cellnames <- c("cell1", "cell2")
#' make10Xcellname(cellnames)
#' @export
make10Xcellname <- function(cellnames, suffix = "_1") {
  paste0(cellnames, suffix)
}



# _________________________________________________________________________________________________
# Soup.Analysis.of.ambient.RNA.R ______________________________ ----
# _________________________________________________________________________________________________
# source('~/GitHub/Packages/Seurat.utils/Functions/Soup.Analysis.of.ambient.RNA.R')
# try (source('https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Soup.Analysis.of.ambient.RNA.R'))
# Source: self + web



# _________________________________________________________________________________________________
#' @title plotTheSoup
#'
#' @description Plot stats about the ambient RNA content in a 10X experiment.
#'
#' @param CellRanger_outs_Dir CellRanger 'outs' (output) directory, Default: '~/Data/114593/114593'
#' @param library_name Aka SampleName (the folder above 'outs;).
#' @param out_dir_prefix Prefix for the output directory. Default: 'SoupStatistics'
#' @param add_custom_class Add a custom class of genes, matched by apattern in gene symbol. Default: `TRUE`.
#' @param pattern_custom The pattern to match in gene symbol. Default: `NA`.
#' @param ls.Alpha The alpha value for the label text. Default: 0.5.
#'
#' @seealso
#'  \code{\link[Matrix]{colSums}}
#'  \code{\link[tibble]{rownames}}
#'  \code{\link[ggrepel]{geom_label_repel}}
#'
#' @importFrom Matrix rowSums
#' @importFrom tibble rownames_to_column
#' @importFrom ggrepel geom_text_repel
#' @importFrom Stringendo percentage_formatter
#' @importFrom MarkdownReports wbarplot create_set_OutDir
#' @importFrom MarkdownHelpers ww.assign_to_global
#' @importFrom dplyr as_tibble
#'
#' @export
plotTheSoup <- function(CellRanger_outs_Dir = "~/Data/114593/114593",
                        # library_name = str_extract(CellRanger_outs_Dir, "[[:alnum:]_]+(?=/outs/)"),
                        library_name = basename(gsub("/outs", "", CellRanger_outs_Dir)),
                        out_dir_prefix = "SoupStatistics",
                        add_custom_class = FALSE, pattern_custom = "\\.RabV$",
                        ls.Alpha = 1) {
  iprint("library_name:", library_name)

  stopifnot( # Check input
    is.character(CellRanger_outs_Dir), dir.exists(CellRanger_outs_Dir),
    nchar(library_name) > 4, is.character(out_dir_prefix), nchar(out_dir_prefix) > 0,
    is.numeric(ls.Alpha)
  )

  if (add_custom_class) iprint("pattern_custom:", pattern_custom)

  # The regular expression `[[:alnum:]_]+(?=/outs/)` matches one or more alphanumeric characters or
  # underscores that are followed by the `/outs/` portion in the string. It ensures that the desired
  # substring is captured, but it does not include the `/outs/` in the matched result.
  # `[[:alnum:]_]+` matches one or more alphanumeric characters or underscores. The `[:alnum:]`
  # character class represents all alphabetic characters (both uppercase and lowercase) and digits.
  # The underscore character `_` is included as well. The `+` quantifier specifies that there should
  # be one or more occurrences of these characters in a row.
  # `(?=/outs/)` This is a positive lookahead assertion. It matches a position in the string where
  # `/outs/` is present immediately after. It doesn't consume any characters from the string; it
  # just checks for the presence of `/outs/` after the matched substring.


  Subfolders_10X_outs <- list.dirs(CellRanger_outs_Dir, full.names = FALSE, recursive = FALSE)
  stopifnot(length(Subfolders_10X_outs) > 0)

  # Identify raw and filtered files ___________________________________
  path.raw <- file.path(CellRanger_outs_Dir, grep(x = Subfolders_10X_outs, pattern = "^raw_*", value = TRUE))
  path.filt <- file.path(CellRanger_outs_Dir, grep(x = Subfolders_10X_outs, pattern = "^filt_*", value = TRUE))
  CR.matrices <- list.fromNames(c("raw", "filt"))

  # Adapter for Markdownreports background variable "OutDir"
  OutDirBac <- if (exists("OutDir")) OutDir else getwd()
  OutDir <- file.path(CellRanger_outs_Dir, paste0(kpp(out_dir_prefix, library_name)))

  MarkdownReports::create_set_OutDir(OutDir)
  MarkdownHelpers::ww.assign_to_global("OutDir", OutDir, 1)

  # Read raw and filtered data ___________________________________
  print("Reading raw CellRanger output matrices")
  CR.matrices$"raw" <- Seurat::Read10X(path.raw)
  if (length(CR.matrices$"raw") == 2) {
    CR.matrices$"raw" <- CR.matrices$"raw"[[1]]
  } # Maybe AB table is present too at slot 2!


  print("Reading filtered CellRanger output matrices")
  CR.matrices$"filt" <- Seurat::Read10X(path.filt)
  if (length(CR.matrices$"filt") == 2) {
    CR.matrices$"filt" <- CR.matrices$"filt"[[1]]
  } # Maybe AB table is present too at slot 2!


  # Profiling the soup ___________________________________
  print("Profiling the soup")
  GEMs.all <- CR.matrices$"raw"@Dimnames[[2]]
  GEMs.cells <- CR.matrices$"filt"@Dimnames[[2]]
  iprint("There are", length(GEMs.all), "GEMs sequenced, and", length(GEMs.cells), "are cells among those.")
  EmptyDroplets.and.Cells <- c("EmptyDroplets" = length(GEMs.all) - length(GEMs.cells), "Cells" = length(GEMs.cells))
  ggExpress::qbarplot(EmptyDroplets.and.Cells, label = EmptyDroplets.and.Cells, palette_use = "npg", col = 1:2, ylab = "GEMs")

  GEMs.soup <- setdiff(GEMs.all, GEMs.cells)
  CR.matrices$"soup" <- CR.matrices$"raw"[, GEMs.soup]
  CR.matrices$"soup.total.RC" <- Matrix::rowSums(CR.matrices$"soup")
  CR.matrices$"soup.total.sum" <- sum(CR.matrices$"soup")
  CR.matrices$"cells.total.sum" <- sum(CR.matrices$"filt")

  CR.matrices$"soup.rel.RC" <- CR.matrices$"soup.total.RC" / CR.matrices$"soup.total.sum"

  # Diff Exp ___________________________________
  Soup.VS.Cells.Av.Exp <- cbind(
    "Soup" = Matrix::rowSums(CR.matrices$"soup"),
    "Cells" = Matrix::rowSums(CR.matrices$"filt")
  )
  colnames(Soup.VS.Cells.Av.Exp)
  idx.HE <- rowSums(Soup.VS.Cells.Av.Exp) > 10
  pc_TRUE(idx.HE)
  Soup.VS.Cells.Av.Exp <- Soup.VS.Cells.Av.Exp[idx.HE, ]
  idim(Soup.VS.Cells.Av.Exp)
  Soup.VS.Cells.Av.Exp.log10 <- log10(Soup.VS.Cells.Av.Exp + 1)

  # ggplot prepare ___________________________________
  Soup.VS.Cells.Av.Exp.gg <- tibble::rownames_to_column(as.data.frame(Soup.VS.Cells.Av.Exp.log10), "gene")
  (Soup.VS.Cells.Av.Exp.gg <- dplyr::as_tibble(Soup.VS.Cells.Av.Exp.gg))
  soup.rate <- Soup.VS.Cells.Av.Exp.gg$Soup / (Soup.VS.Cells.Av.Exp.gg$Cells + Soup.VS.Cells.Av.Exp.gg$Soup)
  cell.rate <- Soup.VS.Cells.Av.Exp.gg$Cells / (Soup.VS.Cells.Av.Exp.gg$Cells + Soup.VS.Cells.Av.Exp.gg$Soup)

  axl.pfx <- "Total Expression in"
  axl.sfx <- "[log10(mRNA+1)]"

  HGNC <- Soup.VS.Cells.Av.Exp.gg$gene
  Class <- rep("Other", times = nrow(Soup.VS.Cells.Av.Exp.gg))
  Class[grep("^RPL|^RPS", HGNC)] <- "RP"
  Class[grep("^MT-", HGNC)] <- "MT"
  Class[grep("^LINC", HGNC)] <- "LINC"
  Class[grep("^AC", HGNC)] <- "AC"
  Class[grep("^AL", HGNC)] <- "AL"
  if (add_custom_class) Class[grep(pattern_custom, HGNC)] <- ReplaceSpecialCharacters(pattern_custom, remove_dots = TRUE)
  Nr.of.Genes.per.Class <- table(Class)


  ggExpress::qpie(Nr.of.Genes.per.Class)
  Soup.VS.Cells.Av.Exp.gg$Class <- Class

  fname <- kpp("Soup.VS.Cells.Av.Exp.GeneClasses", library_name, "pdf")
  pgg <-
    ggplot(
      Soup.VS.Cells.Av.Exp.gg |>
        arrange(-nchar(Class)), aes(x = Soup, y = Cells, label = gene, col = Class)
    ) +
    geom_abline(slope = 1, col = "darkgrey") +
    geom_point() +
    scale_alpha_manual(guide = "none", values = ls.Alpha) +
    xlab(paste(axl.pfx, "Soup", axl.sfx)) +
    ylab(paste(axl.pfx, "Cells", axl.sfx)) +
    ggtitle("Soup VS. Cells | gene classes")

  ggsave(pgg, filename = file.path(OutDir, fname), width = 7, height = 7)

  # ggplot ___________________________________
  quantiles <- c(0.025, 0.01, 0.0025)

  i <- 1
  for (i in 1:length(quantiles)) {
    pr <- quantiles[i]
    print(pr)
    HP.thr <- 200 * pr / quantiles[2]
    idx.HE2 <- rowSums(Soup.VS.Cells.Av.Exp) > HP.thr
    pc_TRUE(idx.HE2)

    fname <- kpp("Soup.VS.Cells.Av.Exp.quantile", pr, library_name, "pdf")

    Outlier <- idx.HE2 &
      (cell.rate < quantile(cell.rate, probs = pr) |
        soup.rate < quantile(soup.rate, probs = pr))

    pc_TRUE(Outlier)
    sum(Outlier)
    HP.thr.mod <- HP.thr
    while (sum(Outlier) > 40) {
      HP.thr.mod <- HP.thr.mod * 2
      Outlier <- Outlier & rowSums(Soup.VS.Cells.Av.Exp) > HP.thr.mod
    }

    pgg <-
      ggplot(Soup.VS.Cells.Av.Exp.gg, aes(
        x = Soup, y = Cells, label = gene,
        col = Outlier
      )) +
      geom_point() +
      theme(legend.position = "none") +
      xlab(paste(axl.pfx, "Soup", axl.sfx)) +
      ylab(paste(axl.pfx, "Cells", axl.sfx)) +
      ggtitle("Soup VS. Cells", subtitle = pr) +
      ggrepel::geom_text_repel(aes(label = ifelse(Outlier,
        as.character(gene), ""
      )))
    ggsave(pgg, filename = file.path(OutDir, fname), width = 7, height = 7)
  } # for


  # Per Gene ___________________________________
  PC.mRNA.in.Soup <- sum(CR.matrices$"soup") / sum(CR.matrices$"raw")
  PC.mRNA.in.Cells <- 100 * sum(CR.matrices$"filt") / sum(CR.matrices$"raw")
  MarkdownReports::wbarplot(
    variable = PC.mRNA.in.Cells, col = "seagreen", plotname = kppd("PC.mRNA.in.Cells", library_name),
    ylim = c(0, 100), ylab = "% mRNA in cells",
    sub = "% mRNA is more meaningful than % reads reported by CR"
  )
  barplot_label(
    barplotted_variable = PC.mRNA.in.Cells,
    labels = Stringendo::percentage_formatter(PC.mRNA.in.Cells / 100, digitz = 2),
    TopOffset = 10
  )


  # Plot top gene's expression ___________________________________
  Soup.GEMs.top.Genes <- 100 * head(sort(CR.matrices$"soup.rel.RC", decreasing = TRUE), n = 20)

  MarkdownReports::wbarplot(Soup.GEMs.top.Genes,
    plotname = kppd("Soup.GEMs.top.Genes", library_name),
    ylab = "% mRNA in the Soup",
    sub = paste("Within the", library_name, "dataset"),
    tilted_text = TRUE,
    ylim = c(0, max(Soup.GEMs.top.Genes) * 1.5)
  )
  barplot_label(
    barplotted_variable = Soup.GEMs.top.Genes,
    labels = Stringendo::percentage_formatter(Soup.GEMs.top.Genes / 100, digitz = 2),
    TopOffset = -.5, srt = 90, cex = .75
  )

  # Plot summarize expression ___________________________________
  soupProfile <- CR.matrices$"soup.total.RC"
  {
    soup.RP.sum <- sum(soupProfile[grep("^RPL|^RPS", names(soupProfile))])
    soup.RPL.sum <- sum(soupProfile[grep("^RPL", names(soupProfile))])
    soup.RPS.sum <- sum(soupProfile[grep("^RPS", names(soupProfile))])
    soup.mito.sum <- sum(soupProfile[grep("^MT-", names(soupProfile))])
    soup.LINC.sum <- sum(soupProfile[grep("^LINC", names(soupProfile))])
    soup.AC.sum <- sum(soupProfile[grep("^AC", names(soupProfile))])
    soup.AL.sum <- sum(soupProfile[grep("^AL", names(soupProfile))])
    genes.non.Above <- soupProfile[CodeAndRoll2::grepv("^RPL|^RPS|^MT-|^LINC|^AC|^AL", names(soupProfile), invert = TRUE)]
  }
  head(sort(genes.non.Above), n = 50)


  soupProfile.summarized <- c(
    "Mitochondial" = soup.mito.sum,
    "Ribosomal" = soup.RP.sum,
    "Ribosomal.L" = soup.RPL.sum,
    "Ribosomal.S" = soup.RPS.sum,
    "GenBank (AC)" = soup.AC.sum,
    "EMBL (AL)" = soup.AL.sum,
    "LINC" = soup.LINC.sum,
    sort(genes.non.Above, decreasing = TRUE)
  )
  NrColumns2Show <- min(10, nrow(soupProfile.summarized))
  ccc <- c("#FF4E00", "#778B04", "#8ea604", "#8ea604", "#F5BB00", "#F5BB00", "#EC9F05", rep(x = "#BF3100", times = NrColumns2Show - 6))


  Soup.GEMs.top.Genes.summarized <- 100 * soupProfile.summarized[1:NrColumns2Show] / CR.matrices$"soup.total.sum"
  maxx <- max(Soup.GEMs.top.Genes.summarized)
  MarkdownReports::wbarplot(Soup.GEMs.top.Genes.summarized,
    plotname = kppd("Soup.GEMs.top.Genes.summarized", library_name),
    ylab = "% mRNA in the Soup", ylim = c(0, maxx + 3),
    sub = paste("Within the", library_name, "dataset"),
    tilted_text = TRUE, col = ccc
  )
  barplot_label(
    barplotted_variable = Soup.GEMs.top.Genes.summarized,
    srt = 45, labels = Stringendo::percentage_formatter(Soup.GEMs.top.Genes.summarized / 100, digitz = 2),
    TopOffset = -1.5
  )

  # Absolute.fraction ___________________________________
  Absolute.fraction.soupProfile.summarized <- Soup.GEMs.top.Genes.summarized * PC.mRNA.in.Soup

  maxx <- max(Absolute.fraction.soupProfile.summarized)
  MarkdownReports::wbarplot(Absolute.fraction.soupProfile.summarized,
    plotname = kppd("Absolute.fraction.soupProfile.summarized", library_name),
    ylab = "% of mRNA in cells", ylim = c(0, maxx * 1.33),
    sub = paste(Stringendo::percentage_formatter(PC.mRNA.in.Soup), "of mRNA counts are in the Soup, in the dataset ", library_name),
    tilted_text = TRUE, col = ccc
  )
  barplot_label(
    barplotted_variable = Absolute.fraction.soupProfile.summarized,
    srt = 45, labels = Stringendo::percentage_formatter(Absolute.fraction.soupProfile.summarized / 100, digitz = 2)
    # formatC(Absolute.fraction.soupProfile.summarized, format="f", big.mark = " ", digits = 0)
    , TopOffset = -maxx * 0.15
  )

  # ___________________________________
  Soup.GEMs.top.Genes.non.summarized <- 100 * sort(genes.non.Above, decreasing = TRUE)[1:20] / CR.matrices$"soup.total.sum"
  maxx <- max(Soup.GEMs.top.Genes.non.summarized)
  MarkdownReports::wbarplot(Soup.GEMs.top.Genes.non.summarized,
    plotname = kppd("Soup.GEMs.top.Genes.non.summarized", library_name),
    ylab = "% mRNA in the Soup",
    sub = paste("Within the", library_name, "dataset"),
    tilted_text = TRUE, col = "#BF3100",
    ylim = c(0, maxx * 1.5)
  )
  barplot_label(
    barplotted_variable = Soup.GEMs.top.Genes.non.summarized,
    labels = Stringendo::percentage_formatter(Soup.GEMs.top.Genes.non.summarized / 100, digitz = 2),
    TopOffset = -maxx * 0.2, srt = 90, cex = .75
  )

  if (exists("OutDirBac")) MarkdownHelpers::ww.assign_to_global("OutDir", OutDirBac, 1)
} # plotTheSoup





# _________________________________________________________________________________________________
# Jaccard.toolkit _____________________________ ----
# _________________________________________________________________________________________________
# try(source('~/GitHub/Packages/Seurat.utils/Functions/Jaccard.toolkit.R'))
# try(source('https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Jaccard.toolkit.R'))


#  __________________________________________
# Fast direct calculation from a list


# _________________________________________________________________________________________________
#' @title jJaccardIndexVec
#'
#' @description Calculate jaccard similarity for 2 vectors. Helper to jPairwiseJaccardIndexList.
#' @param A Set A, Default: 1:3
#' @param B Set B, Default: 2:4
#' @export
jJaccardIndexVec <- function(A = 1:3, B = 2:4) length(intersect(A, B)) / length(union(A, B))

# _________________________________________________________________________________________________
#' @title jPairwiseJaccardIndexList
#'
#' @description Create a pairwise jaccard similarity matrix across all combinations of columns in
#' binary.presence.matrix. Modified from:
#' https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/
#' @param lsG List of genes, Default: ls_genes
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   jPairwiseJaccardIndexList(lsG = ls_genes)
#' }
#' }
#' @export
#' @importFrom Stringendo percentage_formatter
jPairwiseJaccardIndexList <- function(lsG = ls_genes) {
  if (length(names(lsG)) < length(lsG)) {
    iprint("Gene lists were not (all) named, now renamed as:")
    names(lsG) <- ppp("dataset", 1:length(lsG))
    print(names(lsG))
  }
  m <- matrix.fromNames(rowname_vec = names(lsG), colname_vec = names(lsG))
  n.sets <- length(lsG)
  for (r in 1:n.sets) {
    # print(Stringendo::percentage_formatter(r/n.sets))
    for (c in 1:n.sets) {
      if (c == r) {
        m[r, c] <- 1
      } else {
        m[r, c] <- signif(jJaccardIndexVec(lsG[[r]], lsG[[c]]), digits = 2)
      }
    }
  }
  return(m)
}


# Much slower Indirect calculation via PresenceMatrix
# _________________________________________________________________________________________________

# _________________________________________________________________________________________________
#' @title jPresenceMatrix
#'
#' @description Make a binary presence matrix from a list of vectors. Source:
#'   https://stackoverflow.com/questions/56155707/r-how-to-create-a-binary-relation-matrix-from-a-list-of-strings
#' @param string_list List of strings to compare overlapping entries.
#' Default: lst(a = 1:3, b = 2:5, c = 4:9, d = -1:4)
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   df.presence <- jPresenceMatrix(string_list = lst(a = 1:3, b = 2:5, c = 4:9, d = -1:4))
#' }
#' }
#' @export
jPresenceMatrix <- function(string_list = lst(a = 1:3, b = 2:5, c = 4:9, d = -1:4)) {
  df.presence <- string_list |>
    enframe() |>
    unnest(cols = "value") |>
    count(name, value) |>
    spread(value, n, fill = 0)
  df.presence2 <- FirstCol2RowNames(df.presence)
  return(t(df.presence2))
}


# _________________________________________________________________________________________________
#' @title jJaccardIndexBinary
#'
#' @description Calculate the Jaccard Index for two binary vectors. Modified from:
#'   https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/
#' @param x Set X
#' @param y Set Y
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   JaccardSimilarity <- jJaccardIndexBinary(
#'     x = sample(x = 0:1, size = 100, replace = TRUE),
#'     y = sample(x = 0:1, size = 100, replace = TRUE)
#'   )
#' }
#' }
#' @export
jJaccardIndexBinary <- function(x, y) {
  elements.found <- sort(unique(union(x, y)))
  stopifnot(length(elements.found) == 2) # check if you only have [0,1]
  stopifnot(as.numeric(elements.found) == 0:1) # check if you only have [0,1]

  M.11 <- sum(x == 1 & y == 1)
  M.10 <- sum(x == 1 & y == 0)
  M.01 <- sum(x == 0 & y == 1)
  return(M.11 / (M.11 + M.10 + M.01))
}



# _________________________________________________________________________________________________
#' @title jPairwiseJaccardIndex
#'
#' @description Create a pairwise jaccard similarity matrix across all combinations of columns in
#' binary.presence.matrix. Modified from:
#' https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/
#' @param binary.presence.matrix A boolean matrix. Default: df.presence
#' @examples
#' \dontrun{
#' if (interactive()) {
#'   PairwiseJaccardIndices <- jPairwiseJaccardIndex(binary.presence.matrix = df.presence)
#' }
#' }
#' @export
#' @importFrom Stringendo percentage_formatter
jPairwiseJaccardIndex <- function(binary.presence.matrix = df.presence) {
  m <- matrix.fromNames(rowname_vec = colnames(binary.presence.matrix), colname_vec = colnames(binary.presence.matrix))
  n.sets <- ncol(binary.presence.matrix)
  for (r in 1:n.sets) {
    print(Stringendo::percentage_formatter(r / n.sets))
    for (c in 1:n.sets) {
      if (c == r) {
        m[r, c] <- 1
      } else {
        m[r, c] <- signif(jJaccardIndexBinary(binary.presence.matrix[, r], binary.presence.matrix[, c]), digits = 2)
      }
    }
  }
  return(m)
}


# _________________________________________________________________________________________________
# Variable Features _____________________________ ------
# _________________________________________________________________________________________________

#' @title Compare variable features and their ranks in two Seurat objects.
#'
#' @description This function compares variable features (genes) between two Seurat objects,
#'   reporting the number of genes in each, the percentage of common genes, the percentage
#'   of unique genes in each object, and the similarity in the ranking of overlapping genes
#'   using Spearman's rank correlation coefficient. Optionally, it can also generate a scatterplot
#'   of the ranks of common genes using ggpubr's ggscatter. The function returns the common genes
#'   and the Spearman's rank correlation coefficient.
#'
#' @param obj1 The first Seurat object for comparison. Default: NULL.
#' @param obj2 The second Seurat object for comparison. Default: NULL.
#' @param cor.plot An optional boolean indicating whether to generate a scatterplot of the ranks
#'   of common genes. Default: `FALSE`.
#' @param plot_venn plot_venn
#' @param suffix suffix
#' @param save.plot save.plot
#' @return A list containing the common genes and Spearman's rank correlation coefficient.
#'   If cor.plot is TRUE, a scatterplot is also generated.
#' @importFrom Seurat VariableFeatures
#' @importFrom stats cor
#' @importFrom ggExpress qvenn qscatter
#' @examples
#' # Assuming obj1 and obj2 are Seurat objects
#' result <- compareVarFeaturesAndRanks(obj1, obj2, cor.plot = TRUE)
#' @export
compareVarFeaturesAndRanks <- function(
    obj1 = NULL, obj2 = NULL, cor.plot = TRUE, save.plot = TRUE,
    plot_venn = TRUE,
    suffix = NULL,
    ...) {
  stopifnot(!is.null(obj1), !is.null(obj2))
  stopifnot(is(obj1, "Seurat"), is(obj2, "Seurat"))

  name1 <- substitute_deparse(obj1)
  name2 <- substitute_deparse(obj2)

  var.genes1 <- Seurat::VariableFeatures(obj1)
  var.genes2 <- Seurat::VariableFeatures(obj2)

  if (plot_venn) {
    variable.genes.overlap <- list(var.genes1, var.genes2)
    names(variable.genes.overlap) <- c(name1, name2)
    ggExpress::qvenn(list = variable.genes.overlap, suffix = sppp(suffix, c(name1, name2)))
  }

  nr_genes1 <- length(var.genes1)
  nr_genes2 <- length(var.genes2)
  common_genes <- intersect(var.genes1, var.genes2)
  percent_common <- length(common_genes) / max(nr_genes1, nr_genes2) * 100
  percent_uniq1 <- (nr_genes1 - length(common_genes)) / nr_genes1 * 100
  percent_uniq2 <- (nr_genes2 - length(common_genes)) / nr_genes2 * 100

  ranks1 <- match(common_genes, var.genes1)
  ranks2 <- match(common_genes, var.genes2)

  spearman_correlation <- cor(ranks1, ranks2, method = "spearman")

  stopifnot(is.numeric(spearman_correlation))

  cat(sprintf("Nr of genes in obj1: %d\n", nr_genes1))
  cat(sprintf("Nr of genes in obj2: %d\n", nr_genes2))
  cat(sprintf("%% Common genes: %.2f%%\n", percent_common))
  cat(sprintf("%% Unique genes in obj1: %.2f%%\n", percent_uniq1))
  cat(sprintf("%% Unique genes in obj2: %.2f%%\n", percent_uniq2))
  cat(sprintf("Spearman's rank correlation: %.2f\n", spearman_correlation))

  if (cor.plot) {
    plot_data <- data.frame(ranks1, ranks2)
    colnames(plot_data) <- paste("Rank in", c(name1, name2))
    TTL <- paste("Spearman Rank Correlation of Shared Variable Genes")

    SUB <- paste(
      "between objects:", name1, "&", name2, "\n",
      length(common_genes), "or", percent_common, "% overlap from objects:",
      nr_genes1, "&", nr_genes2, "genes."
    )
    CPT <- paste("median ranks:", median(ranks1), "/", median(ranks2))
    file_name <- paste0(
      "Spearman_Rank_Correlation_of_",
      name1, "_and_", name2,
      "_", sprintf("%.2f", spearman_correlation), ".png"
    )
    print(head(plot_data))
    plt <- ggExpress::qscatter(
      df_XYcol = plot_data,
      plotname = TTL,
      subtitle = SUB,
      caption = CPT,
      # abline = c(0,1),
      save = save.plot,
      filename = file_name,
      correlation_r2 = TRUE,
      also.pdf = FALSE,
      cor.coef = TRUE, cor.method = "spearman",
      ...
    )
    print(plt)
  }

  unique.genes <- symdiff(var.genes1, var.genes2)
  names(unique.genes) <- paste0("Unique.", c(name1, name2))
  return(list(
    "common_genes" = common_genes,
    "unique.genes" = unique.genes,
    "spearman_correlation" = spearman_correlation
  ))
}


# _________________________________________________________________________________________________
# Helper Functions _____________________________ ------
# _________________________________________________________________________________________________



# _________________________________________________________________________________________________
#' @title Get the number of CPUs to use for CBE processing
#'
#' @description This function checks for the presence of a global `CBE.params` list and,
#' if found and contains a `cpus` entry, returns the number of CPUs specified by `cpus` minus one.
#' Otherwise, it returns a default number of CPUs.
#'
#' @param n.cpus.def The default number of CPUs to return if `CBE.params` does not exist
#' or does not contain a `cpus` entry. Defaults to 8.
#'
#' @return The number of CPUs to use for CBE processing. If `CBE.params$cpus` is set,
#' returns `CBE.params$cpus - 1`, ensuring at least 1 CPU is returned. Otherwise, returns `n.cpus.def`.
#'
#' @examples
#' # Assuming CBE.params does not exist or does not have a `cpus` entry
#' getCPUsCBE() # returns 8 by default
#'
#' # Assuming CBE.params exists and has a `cpus` entry of 4
#' getCPUsCBE() # returns 3
#' @export
#'
.getNrCores <- function(n.cpus.def = 8) {
  # n_cores_detected <- as.numeric(system("nproc", intern = TRUE))
  n_cores_detected <- as.numeric(system("echo $SLURM_CPUS_PER_TASK", intern = TRUE))
  n_cores_avail <- min(n_cores_detected, n.cpus.def)
  return(max(n_cores_avail, 1))
}


# _________________________________________________________________________________________________
#' @title Check List Elements
#'
#' @description Tests if list elements are defined and reports the value or issues a warning.
#'
#' @param param_list A list containing variables to be checked. Default: `NULL`.
#' @param elements A character vector of element names in `param_list` to check.
#' Default: `character(0)`.
#'
#' @return A message for each element that is defined, and a warning for elements that are not.
#' @examples
#' param_list <- list(a = 1, b = NULL)
#' elements <- c("a", "b", "c")
#' .checkListElements(param_list, elements)
.checkListElements <- function(param_list = NULL, elements = character(0)) {
  stopifnot(
    is.list(param_list),
    is.character(elements)
  )

  sapply(elements, function(element) {
    if (is.null(param_list[[element]])) {
      warning(sprintf("`%s` is not defined", element), immediate. = TRUE, call. = FALSE)
    } else {
      message(sprintf("`%s` is: %s", element, param_list[[element]]))
    }
  }, USE.NAMES = FALSE)

  invisible(NULL)
}


# _________________________________________________________________________________________________
#' @title Get number of scaled features
#'
#' @param obj A Seurat object containing scaled data in  `obj@assays$RNA@scale.data`.
#' @param assay The name of the assay to search for scaled data. Default: `RNA`.
#' @param v Verbose? Default: `TRUE`.
#'
#' @return Integer representing the number of scaled features
.getNrScaledFeatures <- function(obj, assay = Seurat::DefaultAssay(obj),
                                 obj.version = obj@version, v = TRUE) {
  res <- NA
  if (v) message(" > Running .getNrScaledFeatures...")
  if (v) message("Seurat version: ", obj@version, " | Assay searched: ", assay)


  if (obj.version >= "5") { # Check if Seurat version is 5 or higher
    if ("scale.data" %in% names(obj@assays[[assay]]@layers)) {
      res <- nrow(obj@assays[[assay]]@layers[["scale.data"]])
    } else {
      if (v) warning("No scaled data found in object.", immediate. = TRUE)
    }
  } else { # For Seurat versions below 5
    if ("scale.data" %in% names(obj@assays[[assay]])) {
      res <- nrow(obj@assays[[assay]][["scale.data"]])
    } else {
      if (v) warning("No scaled data found in object.", immediate. = TRUE)
    }
  }
  return(res)
}



# _________________________________________________________________________________________________
#' @title Get number of principal components
#'
#' @param obj A Seurat object containing PCA cell embeddings in `reductions$pca@cell.embeddings`
#' @param v Verbose? Default: `TRUE`.
#' @return Integer representing the number of principal components
#'
.getNrPCs <- function(obj, v = TRUE, reduc = "pca") {
  if ("pca" %in% names(obj@reductions)) {
    ncol(obj@reductions[[reduc]]@"cell.embeddings")
  } else {
    if (v) warning("No PCA cell embeddings found in object.", immediate. = TRUE)
    NA
  }
}

# _________________________________________________________________________________________________
#' @title Parse regression variables for name tagging
#'
#' @description This function extracts the regression variables from the `@commands` slot of a Seurat object.
#' If no regression variables are found, a message is printed.
#' @param obj A Seurat object
#' @param assay The name of the assay to search for scaled data. Default: `DefaultAssay()`.
#' @param v Verbose? Default: `TRUE`.
#'
#' @return Integer representing the number of principal components
.getRegressionVariablesForScaleData <- function(obj, assay = Seurat::DefaultAssay(obj), v = TRUE, ...) {
  if (v) message(" > Running .getRegressionVariablesForScaleData...")

  # Input assertions
  stopifnot(is(obj, "Seurat"), is.character(assay))

  # Check if the "commands" slot exists in the object
  if (!"commands" %in% slotNames(obj)) {
    if (v) warning("No commands slot found in object.", immediate. = TRUE)
    return(NULL)
  }

  # Find the ScaleData command using the helper function
  func_slot <- .FindCommandInObject(obj, pattern = paste0("^ScaleData.", assay))

  if (is.null(func_slot)) {
    if (v) message("No ScaleData command found in @commands.")
    return(NULL)
  }

  # Extract regression variables
  regressionVariables <- func_slot$"vars.to.regress"
  if (is.null(regressionVariables)) {
    if (v) message("No regression variables found in @commands")
  } else {
    if (v) message("Regression variables found in @commands: ", paste(regressionVariables, collapse = ", "))
  }

  return(regressionVariables)
}




# _________________________________________________________________________________________________
#' @title Parse key parameters from an object and format as a string
#'
#' @description This function extracts the number of scaled features, the number of principal components,
#' and formats additional information including regression variables.
#' @param obj An object to extract information from.
#' @param regressionVariables A list or vector containing variables for regression. Default: NULL.
#' If NULL, the function will attempt to extract the variables from the `object@commands$ScaleData`.
#' @param nrVarFeatures You can provide this number manually. Default: NULL.
#' @param return.as.name If TRUE, returns the name of the object. Default: `FALSE`.
#' @param assay The assay to extract scaled features from. Default: "RNA".
#' @param suffix A suffix string to add.
#' @param v Verbose? Default: `TRUE`.
#' @return A character string summarizing the key parameters.

.parseKeyParams <- function(obj,
                            regressionVariables = NULL,
                            nrVarFeatures = NULL,
                            return.as.name = FALSE,
                            assay = Seurat::DefaultAssay(obj),
                            suffix = NULL,
                            v = TRUE) {
  tictoc::tic(".parseKeyParams")

  if (v) message(" > Running .parseKeyParams...")
  scaledFeatures <- .getNrScaledFeatures(obj, assay, v = FALSE)

  if (is.null(regressionVariables)) regressionVariables <- .getRegressionVariablesForScaleData(obj = obj, assay = assay, v = FALSE)

  if (!is.null(nrVarFeatures)) {
    if (nrVarFeatures != scaledFeatures) {
      warning("nrVarFeatures !=  scaledFeatures. Reporting nrVarFeatures: ", nrVarFeatures, immediate. = TRUE)
    }
    scaledFeatures <- nrVarFeatures
  } # else use scaledFeatures

  pcs <- .getNrPCs(obj)
  regressionInfo <- kppc(regressionVariables)
  reg <- if (!is.null(regressionVariables)) paste0(regressionInfo, " regressed out") else "no regression"
  if (return.as.name) {
    reg <- ReplaceSpecialCharacters(RemoveWhitespaces(reg, replacement = "."))
    tag <- kpp(scaledFeatures, "ScaledFeatures", pcs, "PCs", reg, suffix)
  } else {
    tag <- paste0(scaledFeatures, " ScaledFeatures | ", pcs, " PCs | ", reg, " ", suffix)
  }
  tictoc::toc()
  return(tag)
}


# _________________________________________________________________________________________________
#' @title Find Command in Seurat Object by Partial Match
#'
#' @description
#' This function searches for commands in a list within a Seurat object using a partial match
#' (e.g., pattern matching) on the command names. It returns the content of the first match if only
#' one match is found. If multiple matches are found, it outputs the number of hits and their names.
#'
#' @param obj A Seurat object. **Default:** None.
#' @param pattern A character string representing the pattern to match command names. **Default:** None.
#'
#' @return If exactly one match is found, the function returns the content of the first match. If
#' multiple matches are found, it returns `NULL` after displaying the number of matches and their names.
#'
#' @examples
#' # Assuming 'combined.obj' is your Seurat object
#' result <- FindCommandInObject(combined.obj, "^FindVariable")
#'
#' @importFrom checkmate assert_class assert_character assert_string

.FindCommandInObject <- function(obj, pattern, perl = TRUE) {
  command_names <- names(obj@commands) # Get all command names

  # Find matches using partial pattern matching
  matches <- grep(pattern, command_names, value = TRUE, perl = perl)

  # Check the number of matches
  if (length(matches) == 0) {
    stop("No matching commands found.")
  } else {
    if (length(matches) > 1) {
      # Multiple matches found, print the number of hits and their names
      message(length(matches), " matches found: ", paste(matches, collapse = ", "))
    }
    # Return the content of the last match
    return(obj@commands[[matches[length(matches)]]])
  }
}




# _________________________________________________________________________________________________
#' @title Parse basic obj stats
#'
#' @description Parse cell and feature count to a string.
#' @param obj An object to extract information from.
#' @return A character string summarizing the key parameters.
#'
.parseBasicObjStats <- function(obj, sep = " ", assay = DefaultAssay(obj),
                                simple = FALSE, suffix = NULL) {
  n.cells <- format(length(Cells(obj)), big.mark = sep, scientific = FALSE)
  n.feat <- format(length(Features(obj, assay = assay)), big.mark = sep, scientific = FALSE)
  if (simple) {
    return(paste(n.cells, "cells.", suffix))
  } else {
    return(paste(n.cells, "cells,", n.feat, assay, "features.", suffix))
  }
}



# _________________________________________________________________________________________________
# New additions,  categorized _____________________________ ------
# _________________________________________________________________________________________________













# _________________________________________________________________________________________________
# Temp _____________________________ ------
# _________________________________________________________________________________________________
