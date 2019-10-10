#' Hierarchical Clustering of Persistence Landscapes by Energy Distance
#' 
#' 
#' @export
landscape.eclust <- function(dlist){
  #############################################
  # Preprocessing : checkers
  if (!dlist_check(dlist)){
    stop("* landscape.eclust : an input 'dlist' should be a list of landscapes as 'kit.landscape' objects. Consult with 'd2landscape' function.")
  }
  
  #############################################
  # Computation
  #   1. 2-norm only.
  dmat = landscape.normpair(dlist, p=2, as.dist=TRUE)
  #   2. apply the function
  output = energy::energy.hclust(dmat)
  
  #############################################
  # Report
  return(output)
}


# # personal test -----------------------------------------------------------
# library(TDA)
# library(TDAkit)
# dlist = list()
# ntest = 25
# for (i in 1:(2*ntest)){
#   x1 = TDA::circleUnif(50)
#   x2 = TDA::circleUnif(50)*0.5
#   x2[,1] = x2[,1] + rep(1,nrow(x2))
#   X <- rbind(x1,x2)
#   X <- X + matrix(rnorm(nrow(X)*ncol(X),sd=0.05), ncol=ncol(X))
#   diagx = ripsDiag(X, maxdimension = 1, maxscale = Inf)$diagram
#   dlist[[i]] = d2landscape(diagx, dimension=0, k=10, inf.replace = FALSE)
# }
# output = landscape.eclust(dlist)