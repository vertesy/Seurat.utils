######################################################################
# Jaccard.toolkit.R
######################################################################
# source('Jaccard.toolkit.R')
# rm(list = ls(all.names = TRUE)); try(dev.off(), silent = T)

# Functions ------------------------
try(source('~/GitHub/CodeAndRoll/CodeAndRoll.R'),silent= F)
require('MarkdownReportsDev')
require('tidyverse')
# source('~/Github/TheCorvinas/R/DatabaseLinke.r')

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

jJaccardIndex = function (x, y) { # Calculate Jaccard Index. Modified from: https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/
  elements.found <- sort(unique(union(x, y)))
  stopifnot(l(elements.found) == 2) # check if you only have [0,1]
  stopifnot(as.numeric(elements.found) == 0:1) # check if you only have [0,1]

  M.11 = sum(x == 1 & y == 1)
  M.10 = sum(x == 1 & y == 0)
  M.01 = sum(x == 0 & y == 1)
  return (M.11 / (M.11 + M.10 + M.01))
}
# JaccardSimilarity <- jJaccardIndex(  x=sample(x = 0:1, size = 100, replace = T)
#               , y=sample(x = 0:1, size = 100, replace = T))


jPairwiseJaccardIndex <- function(binary.presence.matrix = df.presence) { # Create a pairwise jaccard similarity matrix across all combinations of columns in binary.presence.matrix. Modified from: https://www.displayr.com/how-to-calculate-jaccard-coefficients-in-displayr-using-r/
  m = matrix.fromNames(rowname_vec = colnames(binary.presence.matrix), colname_vec = colnames(binary.presence.matrix) )
  n.sets <- ncol(binary.presence.matrix)
  for (r in 1:n.sets) {
    print(percentage_formatter(r/n.sets))
    for (c in 1:n.sets) {
      if (c == r) {
        m[r,c] = 1
      } else {
        m[r,c] = signif(jJaccardIndex(binary.presence.matrix[,r], binary.presence.matrix[,c]), digits = 2)
      }
    }
  }
  return(m)
}
# PairwiseJaccardIndices <- jPairwiseJaccardIndex(binary.presence.matrix = df.presence)
