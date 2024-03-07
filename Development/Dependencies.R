Dependency file generated on Mon Nov 27 17:45:02 2023 

#################################################################################################### 
Seurat.utils.less.used.R
#################################################################################################### 
$`c("package:Seurat", "package:SeuratObject")`
[1] "CreateSeuratObject"

$`package:base`
[1] "basename"  "c"         "is.list"   "length"    "list.dirs" "paste0"    "print"    
[8] "saveRDS"  

$`package:CodeAndRoll2`
[1] "grepv"

$`package:Seurat`
[1] "Read10X"

$`package:Seurat.utils`
[1] "UpdateGenesSeurat"

$`package:Stringendo`
[1] "ppp"

c("Seurat", "SeuratObject")
base
CodeAndRoll2
Seurat
Seurat.utils
Stringendo
#################################################################################################### 
Seurat.Utils.Metadata.R
#################################################################################################### 
$`c(".GlobalEnv", "package:Seurat.utils")`
[1] "clUMAP"

$`c("package:Seurat", "package:SeuratObject")`
[1] "AddMetaData"  "GetAssayData" "Idents"      

$`c("package:SeuratObject", "package:base")`
[1] "colSums"   "intersect"

$`character(0)`
 [1] "arrange"                             "as.numeric.wNames"                  
 [3] "deframe"                             "desc"                               
 [5] "FirstCol2RowNames"                   "ggcorrplot"                         
 [7] "group_by_at"                         "l"                                  
 [9] "melt"                                "n"                                  
[11] "pheatmap"                            "read_rds"                           
[13] "replace_by_most_frequent_categories" "gtools::mixedsort"                       
[15] "stopif"                             "summarise"                          
[17] "summarize_all"                       "summarize_at"                       
[19] "wplot_save_pheatmap"                

$`package:base`
 [1] "all"           "as.character"  "as.data.frame" "basename"      "c"            
 [6] "cat"           "cbind"         "ceiling"       "colnames"      "dim"          
[11] "dput"          "duplicated"    "file.exists"   "floor"         "gsub"         
[16] "inherits"      "is.data.frame" "is.null"       "isFALSE"       "length"       
[21] "levels"        "make.names"    "make.unique"   "match"         "message"      
[26] "min"           "names"         "nrow"          "paste"         "paste0"       
[31] "print"         "rep"           "return"        "rownames"      "sample"       
[36] "scale"         "setdiff"       "stopifnot"     "sub"           "substitute"   
[41] "table"         "unique"        "which"        

$`package:CodeAndRoll2`
[1] "as.named.vector.df" "colMax"             "grepv"              "trail"             
[5] "translate"          "vec.fromNames"      "which_names"       

$`package:ggExpress`
[1] "qbarplot" "qpie"     "qqSave"  

$`package:ggplot2`
[1] "position_stack" "sym"            "vars"          

$`package:ggpubr`
[1] "ggbarplot" "group_by" 

$`package:MarkdownHelpers`
[1] "stopif"

$`package:Seurat`
[1] "FindTransferAnchors" "TransferData"       

$`package:Seurat.utils`
[1] "check.genes"       "GetClusteringRuns" "isave.RDS"        

$`package:stats`
[1] "cor"    "filter" "median"

$`package:Stringendo`
[1] "FixPlotName"          "iprint"               "kollapse"            
[4] "kpp"                  "kpps"                 "percentage_formatter"
[7] "ppp"                 

$`package:utils`
[1] "head"

c(".GlobalEnv", "Seurat.utils")
c("Seurat", "SeuratObject")
c("SeuratObject", "base")
character(0)
base
CodeAndRoll2
ggExpress
ggplot2
ggpubr
MarkdownHelpers
Seurat
Seurat.utils
stats
Stringendo
utils
#################################################################################################### 
Seurat.Utils.R
#################################################################################################### 
$`c(".GlobalEnv", "package:Seurat.utils")`
[1] "clUMAP"                        "multiFeaturePlot.A4"          
[3] "qqSaveGridA4"                  "qUMAP"                        
[5] "SetupReductionsNtoKdimensions"

$`c("package:Seurat", "package:SeuratObject")`
[1] "AddMetaData"        "CreateSeuratObject" "GetAssayData"       "Idents"            
[5] "RenameCells"        "WhichCells"        

$`c("package:SeuratObject", "package:base")`
[1] "colMeans"  "colSums"   "intersect" "rowSums"  

$`character(0)`
 [1] ".saveRDS.compress.in.BG" "AddTrailingSlash"        "arrange"                
 [4] "as_tibble"               "barplot_label"           "check_and_rename"       
 [7] "checkGeneSymbols"        "count"                   "create_set_OutDir"      
[10] "deframe"                 "desc"                    "distinct"               
[13] "enframe"                 "EnhancedVolcano"         "eucl.dist.pairwise"     
[16] "FirstCol2RowNames"       "foreach"                 "geom_text_repel"        
[19] "getActiveProject"        "getDoParRegistered"      "group_by_"              
[22] "group_by_at"             "gunzip"                  "gzip"                   
[25] "isAvailable"             "CodeAndRoll2::split_vec_to_list_by_N"             "job"                    
[28] "left_join"               "list.dirs.depth.n"       "lst"                    
[31] "memory.biggest.objects"  "n"                       "percent_format"         
[34] "percent_rank"            "pheatmap"                "qHGNC"                  
[37] "qread"                   "qsave"                   "read_tsv"               
[40] "read.simple.tsv"         "rownames_to_column"      "rowQuantiles"           
[43] "rowSums2"                "sample_n"                "say"                    
[46] "select"                  "select_at"               "slice"                  
[49] "gtools::mixedsort"            "SoupChannel"             "spread"                 
[52] "str_extract"             "str_split_fixed"         "summarise"              
[55] "summarize"               "tic"                     "toc"                    
[58] "top_n"                   "toTitleCase"             "unnest"                 
[61] "vroom"                   "wbarplot"                "whist"                  
[64] "wlegend"                 "wplot"                   "wplot_save_pheatmap"    
[67] "write.simple.tsv"       

$`package:base`
  [1] "abs"              "all"              "any"              "apply"           
  [5] "as.character"     "as.data.frame"    "as.environment"   "as.list"         
  [9] "as.matrix"        "as.name"          "as.numeric"       "assign"          
 [13] "basename"         "c"                "cat"              "cbind"           
 [17] "ceiling"          "character"        "class"            "colnames"        
 [21] "crossprod"        "data.frame"       "diag"             "dim"             
 [25] "dir.create"       "dir.exists"       "dirname"          "do.call"         
 [29] "dput"             "droplevels"       "exists"           "file.copy"       
 [33] "file.exists"      "file.path"        "file.remove"      "floor"           
 [37] "gc"               "getwd"            "grep"             "grepl"           
 [41] "gsub"             "identical"        "ifelse"           "inherits"        
 [45] "is.character"     "is.list"          "is.null"          "isFALSE"         
 [49] "lapply"           "length"           "library"          "list"            
 [53] "list.dirs"        "list.files"       "log10"            "log2"            
 [57] "make.names"       "make.unique"      "match"            "max"             
 [61] "mean"             "message"          "min"              "names"           
 [65] "nchar"            "ncol"             "nrow"             "options"         
 [69] "paste"            "paste0"           "print"            "rbind"           
 [73] "readline"         "readRDS"          "rep"              "return"          
 [77] "round"            "rownames"         "sample"           "sapply"          
 [81] "save.image"       "saveRDS"          "scale"            "seq_along"       
 [85] "set.seed"         "setdiff"          "setwd"            "signif"          
 [89] "sort"             "sprintf"          "sqrt"             "stop"            
 [93] "stopifnot"        "strsplit"         "subset"           "substitute"      
 [97] "substr"           "sum"              "suppressWarnings" "Sys.Date"        
[101] "Sys.info"         "Sys.setenv"       "system"           "t"               
[105] "table"            "tcrossprod"       "toupper"          "trimws"          
[109] "try"              "tryCatch"         "union"            "unique"          
[113] "unlist"           "warning"          "which"           

$`package:CodeAndRoll2`
 [1] "any.duplicated"     "as.named.vector.df" "col2named.vec.tbl"  "cv"                
 [5] "grepv"              "idim"               "intersect.ls"       "iround"            
 [9] "list.fromNames"     "matrix.fromNames"   "pc_TRUE"            "rowDivide"         
[13] "rowMax"             "rowMin"             "sem"                "sortbyitsnames"    
[17] "splitbyitsnames"    "trail"              "translate"          "unlapply"          
[21] "which_names"       

$`package:ggExpress`
[1] "qbarplot"   "qhistogram" "qpie"       "qqSave"     "qvenn"     

$`package:ggplot2`
 [1] "aes"                 "annotation_logticks" "element_text"        "geom_abline"        
 [5] "geom_bar"            "geom_hline"          "geom_point"          "geom_text"          
 [9] "geom_vline"          "ggplot"              "ggsave"              "ggtitle"            
[13] "labs"                "position_fill"       "scale_alpha_manual"  "scale_fill_manual"  
[17] "scale_x_log10"       "scale_y_discrete"    "scale_y_log10"       "sym"                
[21] "theme"               "theme_classic"       "xlab"                "ylab"               

$`package:ggpubr`
[1] "ggbarplot" "group_by"  "mutate"   

$`package:grDevices`
[1] "rgb"

$`package:MarkdownHelpers`
[1] "filter_HP"           "filter_LP"           "llogit"              "llprint"            
[5] "stopif"              "wcolorize"           "ww.assign_to_global"

$`package:methods`
[1] "is"        "slot"      "slotNames"

$`package:Seurat`
 [1] "AverageExpression"    "DiscretePalette"      "FindClusters"        
 [4] "FindNeighbors"        "FindVariableFeatures" "FontSize"            
 [7] "LabelPoints"          "Read10X"              "RunPCA"              
[10] "RunTSNE"              "ScaleData"           

$`package:Seurat.utils`
 [1] "calc.q99.Expression.and.set.all.genes" "check.genes"                          
 [3] "dropLevelsSeurat"                      "GetClusteringRuns"                    
 [5] "GetNamedClusteringRuns"                "GetOrderedClusteringRuns"             
 [7] "getProject"                            "GetTopMarkersDF"                      
 [9] "GetUpdateStats"                        "HGNC.EnforceUnique"                   
[11] "isave.RDS"                             "jJaccardIndexBinary"                  
[13] "jJaccardIndexVec"                      "remove.residual.small.clusters"       
[15] "RenameGenesSeurat"                     "sampleNpc"                            
[17] "shorten_clustering_names"              "SmallestNonAboveX"                    
[19] "sparse.cor"                            "downsampleSeuObj"                         
[21] "UpdateGenesSeurat"                     "ww.get.1st.Seur.element"              
[23] "xsave"                                

$`package:SeuratObject`
[1] "plan"

$`package:stats`
[1] "ecdf"     "filter"   "median"   "quantile"

$`package:Stringendo`
 [1] "FixPath"              "FixPlotName"          "flag.nameiftrue"     
 [4] "idate"                "iprint"               "kollapse"            
 [7] "kpp"                  "kppd"                 "kpps"                
[10] "kppu"                 "percentage_formatter" "ppp"                 
[13] "sppp"                

$`package:utils`
[1] "head"     "read.csv"

c(".GlobalEnv", "Seurat.utils")
c("Seurat", "SeuratObject")
c("SeuratObject", "base")
character(0)
base
CodeAndRoll2
ggExpress
ggplot2
ggpubr
grDevices
MarkdownHelpers
methods
Seurat
Seurat.utils
SeuratObject
stats
Stringendo
utils
#################################################################################################### 
Seurat.Utils.Visualization.R
#################################################################################################### 
$`c(".GlobalEnv", "package:Seurat.utils")`
 [1] "Annotate4Plotly3D"                         
 [2] "BackupReduction"                           
 [3] "getDiscretePalette"                        
 [4] "gg_color_hue"                              
 [5] "multiFeaturePlot.A4"                       
 [6] "plot3D.umap"                               
 [7] "plot3D.umap.gene"                          
 [8] "qUMAP"                                     
 [9] "RecallReduction"                           
[10] "SavePlotlyAsHtml"                          
[11] "seu.PC.var.explained"                      
[12] "ww.check.if.3D.reduction.exist"            
[13] "ww.check.quantile.cutoff.and.clip.outliers"

$`c("package:Seurat", "package:SeuratObject")`
[1] "DefaultAssay" "Embeddings"   "FetchData"    "GetAssayData"

$`c("package:SeuratObject", "package:base")`
[1] "colSums"   "intersect" "rowSums"  

$`c("package:sp", "package:base")`
[1] "split"

$`c("package:sp", "package:graphics", "package:base")`
[1] "plot"

$`character(0)`
 [1] "barplot_label"              "MarkdownHelpers::color_check"                "create_set_Original_OutDir"
 [4] "create_set_OutDir"          "create_set_SubDir"          "deframe"                   
 [7] "FeatureHeatmap"             "FirstCol2RowNames"          "group_by_"                 
[10] "hue_pal"                    "CodeAndRoll2::split_vec_to_list_by_N"                "n"                         
[13] "oo"                         "plot_grid"                  "plot_ly"                   
[16] "principal_curve"            "save_plot"                  "saveWidget"                
[19] "select"                     "stopif"                    "summarise"                 
[22] "summarize"                  "tic"                        "toc"                       
[25] "wbarplot"                   "whiskers"                   "wplot_save_this"           

$`package:base`
 [1] "all"            "anyNA"          "apply"          "as.character"   "as.environment"
 [6] "as.factor"      "as.matrix"      "as.name"        "c"              "cbind"         
[11] "ceiling"        "colnames"       "cut"            "data.frame"     "do.call"       
[16] "environment"    "exists"         "floor"          "getwd"          "is.null"       
[21] "isFALSE"        "lapply"         "length"         "levels"         "list"          
[26] "list2env"       "make.names"     "mean"           "missing"        "names"         
[31] "ncol"           "nrow"           "paste"          "paste0"         "print"         
[36] "range"          "rep"            "return"         "rm"             "round"         
[41] "row.names"      "rownames"       "seq"            "setdiff"        "sort"          
[46] "stop"           "stopifnot"      "substitute"     "sum"            "table"         
[51] "toupper"        "try"            "unique"         "unlist"         "warning"       

$`package:CodeAndRoll2`
 [1] "as_tibble_from_namedVec"     "as.named.vector.df"         
 [3] "clip.at.fixed.value"         "clip.outliers.at.percentile"
 [5] "iround"                      "na.omit.strip"              
 [7] "pc_TRUE"                     "splitbyitsnames"            
 [9] "translate"                   "unlapply"                   

$`package:ggExpress`
[1] "qA4_grid_plot" "qbarplot"      "qhistogram"    "qqSave"       

$`package:ggplot2`
 [1] "aes"                 "annotation_logticks" "coord_fixed"         "element_blank"      
 [5] "element_text"        "geom_histogram"      "geom_hline"          "geom_point"         
 [9] "geom_vline"          "ggplot"              "ggsave"              "ggtitle"            
[13] "labs"                "scale_x_log10"       "scale_y_log10"       "theme"              
[17] "theme_bw"            "theme_linedraw"      "theme_set"          

$`package:ggpubr`
[1] "group_by"

$`package:graphics`
[1] "layout" "points"

$`package:grDevices`
[1] "hcl"

$`package:MarkdownHelpers`
[1] "color_check"   "filter_HP"     "jjpegA4"       "llprint"       "stopif"       
[6] "TRUE.unless"   "try.dev.off"   "ww.FnP_parser"

$`package:Seurat`
 [1] "DimPlot"         "DiscretePalette" "FeaturePlot"     "FeatureScatter" 
 [5] "NoAxes"          "NoLegend"        "RunPCA"          "RunTSNE"        
 [9] "RunUMAP"         "SplitObject"     "VlnPlot"        

$`package:Seurat.utils`
[1] "check.genes"            "getCellIDs.from.meta"   "GetClusteringRuns"     
[4] "GetNamedClusteringRuns" "GetTopMarkers"         

$`package:stats`
[1] "quantile"

$`package:Stringendo`
 [1] "extPNG"               "FixPlotName"          "flag.nameiftrue"     
 [4] "idate"                "iprint"               "kollapse"            
 [7] "kpp"                  "kpps"                 "kppu"                
[10] "percentage_formatter" "ppp"                  "ppu"                 
[13] "sppp"                

$`package:utils`
[1] "head"

c(".GlobalEnv", "Seurat.utils")
c("Seurat", "SeuratObject")
c("SeuratObject", "base")
c("sp", "base")
c("sp", "graphics", "base")
character(0)
base
CodeAndRoll2
ggExpress
ggplot2
ggpubr
graphics
grDevices
MarkdownHelpers
Seurat
Seurat.utils
stats
Stringendo
utils
