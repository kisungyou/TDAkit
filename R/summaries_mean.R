#' Mean of Multiple Functional Summaries
#' 
#' Given multiple functional summaries \eqn{\Lambda_1 (t), \Lambda_2 (t), \ldots, \Lambda_N (t)},
#' compute the mean 
#' \deqn{\bar{\Lambda} (t) = \frac{1}{N} \sum_{n=1}^N \Lambda_n (t)}.
#' 
#' @param fslist a length-\eqn{N} list of functional summaries of persistent diagrams.
#' 
#' @return a functional summary object.
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
#' ldsum = fsmean(list_land)
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
#' @concept summaries
#' @export
fsmean <- function(fslist){
  ## PREPROCESSING
  dtype  = check_list_summaries("fsmean", fslist)
  nlist  = length(fslist)
  weight = rep(1/nlist, nlist)
  
  ## COMPUTE AND RETURN VIA FSSUM
  output = TDAkit::fssum(fslist, weight)
  return(output)
}