######################################################################
# 00.Load.Seurat.Utils.LOCAL.R
######################################################################
# source('~/GitHub/Packages/Seurat.utils/00.Load.Seurat.Utils.LOCAL.R')
# try(source('https://raw.githubusercontent.com/vertesy/Seurat.utils/master/00.Load.Seurat.Utils.LOCAL.R'), silent =   T)

try(source("~/GitHub/Packages/Seurat.utils/Functions/Seurat.update.gene.symbols.HGNC.R"))
try(source("~/GitHub/Packages/Seurat.utils/Functions/Metadata.manipulation.R"))
try(source("~/GitHub/Packages/Seurat.utils/Functions/Plotting.dim.reduction.2D.R"))
try(source("~/GitHub/Packages/Seurat.utils/Functions/Plotting.dim.reduction.3D.R"))
try(source("~/GitHub/Packages/Seurat.utils/Functions/Plotting.filtering.R"))
try(source("~/GitHub/Packages/Seurat.utils/Functions/Plotting.statistics.and.QC.R"))
try(source("~/GitHub/Packages/Seurat.utils/Functions/Read.Write.Save.Load.functions.R"))
try(source("~/GitHub/Packages/Seurat.utils/Functions/Seurat.object.manipulations.etc.R"))
try(source("~/GitHub/Packages/Seurat.utils/Functions/Cluster.Auto-naming.DE.R"))
try(source('~/GitHub/Packages/Seurat.utils/Functions/MULTI-seq.functions.R'))
try(source('~/GitHub/Packages/Seurat.utils/Functions/Jaccard.toolkit.R'))

# try(source('~/GitHub/Packages/Seurat.utils/Functions/Soup.Analysis.of.ambient.RNA.R'))
# try (source("https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Custom.Functions.for.Slingshot.R"))

# source('https://raw.githubusercontent.com/vertesy/DatabaseLinke.R/master/DatabaseLinke.R')


# try(source('~/GitHub/Projects/SEO/GO-scoring/Seurat.gene.sets.and.GO.terms.R'))
try(source('~/GitHub/TheCorvinas/R/GO-scoring/Seurat.gene.sets.and.GO.terms.R'))

# Requirements ------------------------
# May also require
# try (source('/GitHub/Packages/CodeAndRoll/CodeAndRoll.R'),silent= F) # generic utilities funtions
# require('MarkdownReportsDev') # require("devtools") # plotting related utilities functions # devtools::install_github(repo = "vertesy/MarkdownReportsDev")
