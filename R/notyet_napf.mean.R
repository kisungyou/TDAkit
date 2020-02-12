#' Mean of Multiple Multiple Normalized Accumulative Persistence Functions
#' 
#' @export
napf.mean <- function(flist){
  #############################################
  # Preprocessing : checkers
  myweight = rep(1,length(flist))/length(flist)
  if (!check_list_napf(flist)){
    stop("* napf.mean : an input 'flist' should be a list of NAPFs as 'kit.napf' objects. Consult with 'd2napf' function.")
  }
  
  #############################################
  # Main Computation using 'napf.wsum'
  return(napf.wsum(flist, myweight))
}


# # personal test -----------------------------------------------------------
# library(TDA)
# library(TDAkit)
# 
# # data generation
# alist = list()
# ntest = 10
# for (i in 1:ntest){
#   X = TDAkit::gen.2circles(n=100)
#   diagx = ripsDiag(X, maxdimension = 1, maxscale = Inf)$diagram
#   alist[[i]] = d2napf(diagx, dimension=1)
# }
# 
# # mean
# amean <- napf.mean(alist)
# 
# # visualize
# opar  <- par(mfrow=c(2,3))
# for (i in 1:5){
#   plot(alist[[i]]$tseq, alist[[i]]$lambda, "l", main=paste0("mean ",i))
# }
# plot(amean$tseq, amean$lambda, "l", main="mean NAPF")