#' Mean of Multiple Persistence Landscapes
#' 
#' @export
landscape.mean <- function(dlist){
  #############################################
  # Preprocessing : checkers
  weight = rep(1,length(dlist))/length(dist)
  if (!dlist_check(dlist)){
    stop("* landscape.mean : an input 'dlist' should be a list of landscapes as 'kit.landscape' objects. Consult with 'd2landscape' function.")
  }
  if (length(weight)!=length(dlist)){
    stop("* landscape.mean : length of 'weight' vector should equal to the number of landscapes in 'dlist'.")
  }
  
  #############################################
  # Main Computation
  mainout = dlist_adjust(dlist, as.list=FALSE)
  maincom = compute_slicewsum(mainout$array3d, weight)
  
  #############################################
  # Report the results
  output = list()
  output$lambda = maincom
  output$tseq   = mainout$tseq
  output$dimension = dlist[[1]]$dimension
  class(output) = "kit.landscape"
  return(output)
}



# # personal test -----------------------------------------------------------
# library(TDA)
# dlist = list()
# for (i in 1:10){
#   x1 = TDA::circleUnif(30)
#   x2 = TDA::circleUnif(30)*0.5
#   x2[,1] = x2[,1] + rep(1,nrow(x2))
#   X <- rbind(x1,x2)
#   X <- X + matrix(rnorm(nrow(X)*ncol(X),sd=0.1), ncol=ncol(X))
#   diagx = ripsDiag(X, maxdimension = 1, maxscale = Inf)$diagram
#   dlist[[i]] = d2landscape(diagx, dimension=1, k=0, inf.replace = FALSE)
# }
# output = landscape.wsum(dlist)
# plot(output$tseq, output$lambda[,1], "l")
# lines(output$tseq, output$lambda[,2], col="red")