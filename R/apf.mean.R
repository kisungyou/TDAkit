#' Mean of Multiple Accumulative Persistence Functions
#' 
#' @export
apf.mean <- function(alist){
  #############################################
  # Preprocessing : checkers
  myweight = rep(1,length(alist))/length(alist)
  if (!check_list_apf(alist)){
    stop("* apf.mean : an input 'alist' should be a list of APFs as 'kit.apf' objects. Consult with 'd2apf' function.")
  }
  
  #############################################
  # Main Computation using 'apf.wsum'
  return(apf.wsum(alist, myweight))
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
#   alist[[i]] = d2apf(diagx, dimension=1)
# }
# 
# # mean
# amean <- apf.mean(alist)
# 
# # visualize
# opar  <- par(mfrow=c(2,3))
# for (i in 1:5){
#   plot(alist[[i]]$tseq, alist[[i]]$lambda, "l", main=paste0("mean ",i))
# }
# plot(amean$tseq, amean$lambda, "l", main="mean APF")