######################################################################
# Jaccard.toolkit.R
######################################################################
# try(source('~/GitHub/Packages/Seurat.utils/Functions/Jaccard.toolkit.R'))
# try(source('https://raw.githubusercontent.com/vertesy/Seurat.utils/master/Functions/Jaccard.toolkit.R'))


# Functions ------------------------
try(source('~/GitHub/Packages/CodeAndRoll/CodeAndRoll.R'),silent= F)
try(require('MarkdownReportsDev'),  silent = T)
try(require('tidyverse'),  silent = T)
# source('~/Github/TheCorvinas/R/DatabaseLinke.r')


# --------------------------------------------------------------------------------------------
# Fast direct calculation from a list --------------------------------------------------------
# --------------------------------------------------------------------------------------------


# jJaccardIndexVec ----------------------------------------
jJaccardIndexVec <- function(A = 1:3, B = 2:4) length(intersect(A,B)) / length(union(A,B))

# jPairwiseJaccardIndexList ----------------------------------------

jPairwiseJaccardIndexList <- function(lsG = ls_genes) { # Create a pairwise jaccard similarity matrix across all combinations of columns in binary.presence.matrix. Modified from: https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/
  if (l(names(lsG)) < l(lsG)) {
    iprint("Gene lists were not (all) named, now renamed as:")
    names(lsG) <- ppp("dataset", 1:l(lsG))
    print(names(lsG))
  }
  m = matrix.fromNames(rowname_vec = names(lsG), colname_vec = names(lsG))
  n.sets <- length(lsG)
  for (r in 1:n.sets) {
    # print(percentage_formatter(r/n.sets))
    for (c in 1:n.sets) {
      if (c == r) {
        m[r,c] = 1
      } else {
        m[r,c] =signif(jJaccardIndexVec(lsG[[r]], lsG[[c]]), digits = 2)
      }
    }
  }
  return(m)
}
# jPairwiseJaccardIndexList(lsG = ls_genes)






# --------------------------------------------------------------------------------------------
# Much slower Indirect calculation via PresenceMatrix ----------------------------------------
# --------------------------------------------------------------------------------------------


# jPresenceMatrix ----------------------------------------
jPresenceMatrix <- function(string_list = lst(a=1:3, b=2:5,c=4:9, d=-1:4) ) { # Make a binary presence matrix from a list. Source: https://stackoverflow.com/questions/56155707/r-how-to-create-a-binary-relation-matrix-from-a-list-of-strings
  df.presence <- string_list %>%
    enframe %>%
    unnest(cols = "value") %>%
    count(name, value) %>%
    spread(value, n, fill = 0)
  df.presence2 <- FirstCol2RowNames(df.presence)
  return(t(df.presence2))
}
# df.presence <- jPresenceMatrix(string_list = lst(a=1:3, b=2:5,c=4:9, d=-1:4))

# jJaccardIndexBinary ----------------------------------------
jJaccardIndexBinary = function (x, y) { # Calculate Jaccard Index. Modified from: https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/
  elements.found <- sort(unique(union(x, y)))
  stopifnot(l(elements.found) == 2) # check if you only have [0,1]
  stopifnot(as.numeric(elements.found) == 0:1) # check if you only have [0,1]

  M.11 = sum(x == 1 & y == 1)
  M.10 = sum(x == 1 & y == 0)
  M.01 = sum(x == 0 & y == 1)
  return (M.11 / (M.11 + M.10 + M.01))
}
# JaccardSimilarity <- jJaccardIndexBinary(  x=sample(x = 0:1, size = 100, replace = T)
#               , y=sample(x = 0:1, size = 100, replace = T))



# jPairwiseJaccardIndex ----------------------------------------
jPairwiseJaccardIndex <- function(binary.presence.matrix = df.presence) { # Create a pairwise jaccard similarity matrix across all combinations of columns in binary.presence.matrix. Modified from: https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/
  m = matrix.fromNames(rowname_vec = colnames(binary.presence.matrix), colname_vec = colnames(binary.presence.matrix) )
  n.sets <- ncol(binary.presence.matrix)
  for (r in 1:n.sets) {
    print(percentage_formatter(r/n.sets))
    for (c in 1:n.sets) {
      if (c == r) {
        m[r,c] = 1
      } else {
        m[r,c] = signif(jJaccardIndexBinary(binary.presence.matrix[,r], binary.presence.matrix[,c]), digits = 2)
      }
    }
  }
  return(m)
}
# PairwiseJaccardIndices <- jPairwiseJaccardIndex(binary.presence.matrix = df.presence)

