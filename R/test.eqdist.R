#' Energy-based Testing Equality of Distributions for Persistence Landscapes
#' 
#' 
#' @export
test.eqdist <- function(dlist, label, method=c("original","discoB","discoF"), mc.iter=496){
  #############################################
  # Preprocessing : checkers
  if (!dlist_check(dlist)){ # dlist : list of landscapes
    stop("* test.eqdist : an input 'dlist' should be a list of landscapes as 'kit.landscape' objects. Consult with 'd2landscape' function.")
  }
  nlist = length(dlist)
  label = round(label)      # label vector
  if ((!is.vector(label))||(length(label)!=nlist)||(any(is.infinite(label)))||(any(is.na(label)))||(missing(label))){
    stop("* test.eqdist : 'label' vector of same length as 'dlist' must be provided without any NaN or Infs.")
  }
  mymethod = match.arg(method)
  mc.iter  = round(mc.iter)
  
  
  #############################################
  # Computation : Preliminary
  #   pairwise L2 distance
  dmat = as.matrix(landscape.normpair(dlist, p=2)) # (nlist x nlist) L2 distance matrix
  #   let's take care of label
  ulabels = sort(unique(label), decreasing = FALSE)
  nunique = length(ulabels)
  #   re-arrange distance matrix and label vector as sizes 
  orderlab = base::order(label)
  mat.dist = dmat[orderlab,orderlab]
  vec.size = rep(0,nunique)
  for (i in 1:nunique){
    vec.size[i] = length(which(label==ulabels[i]))
  }
  
  #############################################
  # Computation : Main Case
  if (all(mymethod=="original")){
    output = energy::eqdist.etest(mat.dist, vec.size, distance=TRUE, method="original", R=mc.iter)
  } else if (all(mymethod=="discoB")){
    output = energy::eqdist.etest(mat.dist, vec.size, distance=TRUE, method="discoB", R=mc.iter)
  } else {
    output = energy::eqdist.etest(mat.dist, vec.size, distance=TRUE, method="discoF", R=mc.iter)
  }
  
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
# test.eqdist(dlist, c(rep(1,ntest), rep(2,ntest)), method="original")
# test.eqdist(dlist, c(rep(1,ntest), rep(2,ntest)), method="discoB")
# test.eqdist(dlist, c(rep(1,ntest), rep(2,ntest)), method="discoF")
