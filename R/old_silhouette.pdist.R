#' Pairwise \eqn{L_p} Distance of Multiple Persistence Silhouettes
#' 
#' @export
silhouette.pdist <- function(slist, p=2, as.dist=FALSE){
  #############################################
  # Preprocessing : checkers
  if (!check_list_silhouette(slist)){
    stop("* landscape.pdist : 'dlist' is invalid. Consult with 'd2landscape' for each element in the list.")
  }
  if (p < 1){
    stop("* landscape.pdist : 'p'-norm is defined for the range [1,Inf]. ")
  }
  nlist = length(slist)
  p     = as.double(p)
  
  #############################################
  # Main Computation
  #   1. rearrange as a list
  list2d = adjust_list_silhouette(slist, as.list=FALSE)
  mytseq = list2d$tseq
  #   2. iterate
  myeps  = 10*.Machine$double.eps
  output = array(0,c(nlist,nlist))
  for (i in 1:(nlist-1)){
    tgti = list2d$array[,i]
    for (j in (i+1):nlist){
      tgtj = list2d$array[,j]
      if (is.infinite(p)){
        output[i,j] <- output[j,i] <- max(abs(tgti-tgtj), na.rm=TRUE)
      } else {
        output[i,j] <- output[j,i] <- (simple_integral_1d((abs(tgti-tgtj)^p), mytseq)^(1/p))
      }
    }
  }
  
  #############################################
  # Report
  if (as.dist){
    return(stats::as.dist(output))
  } else {
    return(output)
  }
}

# # # personal test -----------------------------------------------------------
# library(TDA)
# library(TDAkit)
# slist = list()
# ntest = 10
# for (i in 1:ntest){
#   X = TDA::circleUnif(100)
#   diagx = ripsDiag(X, maxdimension = 1, maxscale = Inf)$diagram
#   slist[[i]] = d2silhouette(diagx, dimension=1)
# }
# for (i in (ntest+1):(2*ntest)){
#   X = TDAkit::gen.2circles(n=100)
#   diagx = ripsDiag(X, maxdimension = 1, maxscale = Inf)$diagram
#   slist[[i]] = d2silhouette(diagx, dimension=1)
# }
# image(silhouette.pdist(slist))

