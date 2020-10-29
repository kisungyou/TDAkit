#' Weighted Sum of Multiple Persistence Landscapes
#' 
#' Given multiple landscape functions \eqn{\Lambda_1 (t), \Lambda_2 (t), \ldots, \Lambda_N (t)} of same dimension, 
#' compute the weighted sum 
#' \deqn{\bar{\Lambda} (t) = \sum_{n=1}^N w_n \Lambda_n (t)}
#' given weights \eqn{w_1,w_2,\ldots,w_N}.
#' 
#' @param dlist a length-\eqn{N} list of \code{'landscape'} objects from \code{\link{diag2landscape}}.
#' @param weight a weight vector of length \eqn{N}. If \code{NULL} (default), weights are automatically set as \eqn{w_1=\cdots=\w_N = 1/N}.
#' 
#' @return a list object of \code{"landscape"} class containing\describe{
#' \item{lambda}{an \eqn{(\code{nseq} \times k)} landscape functions.}
#' \item{tseq}{a length-\code{nseq} vector of domain grid.}
#' \item{dimension}{dimension of features considered.}
#' }
#' 
#' @examples 
#' # ---------------------------------------------------------------------------
#' #     Weighted Average of 10 Persistence Landscapes from '2holes' data
#' # ---------------------------------------------------------------------------
#' ## Generate 10 Diagrams with 'gen2holes()' function
#' list_rips = list()
#' for (i in 1:10){
#'   list_rips[[i]] = diagRips(gen2holes(n=100, sd=2)$data, maxdim=1)
#' }
#' 
#' ## Compute Persistence Landscapes from Each Diagram with k=5 Functions
#' list_land = list()
#' for (i in 1:10){
#'   list_land[[i]] = diag2landscape(list_rips[[i]], dimension=0, k=5)
#' }
#' 
#' ## Some Random Weights
#' wrand = abs(stats::rnorm(10))
#' wrand = wrand/sum(wrand)
#' 
#' ## Compute Weighted Sum of Landscapes
#' ldsum = plsum(list_land, weight=wrand)
#' 
#' ## Visualize
#' sam5  <- sort(sample(1:10, 5, replace=FALSE))
#' opar  <- par(no.readonly=TRUE)
#' par(mfrow=c(2,3), pty="s")
#' for (i in 1:5){
#'   tgt = list_land[[sam5[i]]]
#'   matplot(tgt$tseq, tgt$lambda[,1:5], type="l", lwd=3, main=paste("landscape no.",sam5[i]))
#' }
#' matplot(ldsum$tseq, ldsum$lambda[,1:5], type="l", lwd=3, main="weighted sum")
#' par(opar)
#' 
#' @references 
#' \insertRef{bubenik_persistence_2018}{TDAkit}
#' 
#' @concept landscape
#' @export
plsum <- function(dlist, weight=NULL){
  ## PREPROCESSING
  if (!check_list_landscape(dlist)){
    stop("* plsum : 'dlist' should be a list of 'landscape' objects.")
  }
  if ((length(weight)<1)&&(is.null(weight))){
    nn     = length(dlist)
    weight = rep(1/nn, nn)
  }
  if (length(weight)!=length(dlist)){
    stop("* plsum : length of 'dlist' should match to that of 'weight'.")
  }
  
  ## MAIN COMPUTATION
  mainout = adjust_list_landscapes(dlist, as.list=FALSE)
  maincom = compute_slicewsum(mainout$array3d, weight)
    
  ## REPORT
  output = list()
  output$lambda = maincom
  output$tseq   = as.vector(mainout$tseq)
  output$dimension = dlist[[1]]$dimension
  class(output) = "landscape"
  return(output)
}
