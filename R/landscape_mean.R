#' Mean of Multiple Persistence Landscapes
#' 
#' Given \eqn{N} landscape functions \eqn{\Lambda_1 (t), \Lambda_2 (t), \ldots, \Lambda_N (t)} of same dimension, 
#' compute the weighted sum 
#' \deqn{\bar{\Lambda} (t) = \frac{1}{N} \sum_{n=1}^N \Lambda_n (t)}.
#' 
#' @param dlist a length-\eqn{N} list of \code{'landscape'} objects from \code{\link{diag2landscape}}.
#' 
#' @return a list object of \code{"landscape"} class containing\describe{
#' \item{lambda}{an \eqn{(\code{nseq} \times k)} landscape functions.}
#' \item{tseq}{a length-\code{nseq} vector of domain grid.}
#' \item{dimension}{dimension of features considered.}
#' }
#' 
#' @examples 
#' # ---------------------------------------------------------------------------
#' #         Mean of 10 Persistence Landscapes from '2holes' data
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
#' ## Compute Weighted Sum of Landscapes
#' ldsum = plsum(list_land)
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
plmean <- function(dlist){
  ## PREPROCESSING
  if (!check_list_landscape(dlist)){
    stop("* plmean : 'dlist' should be a list of 'landscape' objects.")
  }
  weight = rep(1,length(dlist))/length(dlist)
  
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
