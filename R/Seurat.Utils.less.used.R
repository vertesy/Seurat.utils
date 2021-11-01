

# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------------------------------
#' @title cell_pal
#' @description Generate paletter
#' @param cell_vars PARAM_DESCRIPTION
#' @param pal_fun PARAM_DESCRIPTION
#' @param ... Pass any other parameter to the internally called functions (most of them should work).
#' @examples
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export
cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# Duplicates ------------------------------------------------------------------------------------------------

# ------------------------------------------------------------------------
# #' @title Calc.Cor.Seurat
# #' @description
# #' @param assay.use PARAM_DESCRIPTION, Default: 'RNA'
# #' @param slot.use PARAM_DESCRIPTION, Default: 'data'
# #' @param geneset PARAM_DESCRIPTION, Default: FALSE
# #' @param quantileX Quantile level, Default: 0.95
# #' @param max.cells PARAM_DESCRIPTION, Default: 10000
# #' @param seed random seed used, Default: p$seed
# #' @param digits PARAM_DESCRIPTION, Default: 2
# #' @param obj Seurat object, Default: combined.obj
# #' @examples
# #' \dontrun{
# #' if(interactive()){
# #'  #EXAMPLE1
# #'  }
# #' }
# #' @export
# Calc.Cor.Seurat <- function(assay.use = "RNA", slot.use = "data", geneset = FALSE
#                             , quantileX = 0.95, max.cells =  10000, seed = p$"seed"
#                             , digits = 2, obj = combined.obj) {
#   expr.mat <- GetAssayData(slot = slot.use, assay = assay.use, object = obj)
#   if (ncol(expr.mat) > max.cells) {
#     set.seed(seed = seed)
#     cells.use <- sample(x = colnames(expr.mat), size = max.cells)
#   } else {
#     cells.use <- colnames(obj)
#   }

#   qname = p0("q", quantileX * 100)
#   quantile_name = kpp("expr", qname)
#   if (is.null(obj@misc[[quantile_name]])) { iprint("Quantile data missing! Call: combined.obj <- calc.q90.Expression.and.set.all.genes(combined.obj, quantileX =",quantileX,") first!"); stop()}

#   genes.HE  <- if (isFALSE(geneset)) {  which_names(obj@misc[[quantile_name]] > 0) } else {
#     check.genes(geneset)  }
#   iprint("Pearson correlation is calculated for", l(genes.HE), "HE genes with expr."
#          , qname,": > 0 on a sample of", max.cells, " cells.")
#   tic(); ls.cor <- sparse.cor(smat = t(expr.mat[genes.HE, cells.use])); toc()

#   ls.cor <- lapply(ls.cor, round, digits = 2)

#   slot__name <- kpp(slot.use, assay.use, quantile_name)
#   obj@misc[[kpp('cor', slot__name)]] <- ls.cor$'cor'
#   obj@misc[[kpp('cov', slot__name)]] <- ls.cor$'cov'
#   iprint("Stored under obj@misc$", kpp('cor', slot.use, assay.use), "or cov... ." )
#   return(obj)
# }

