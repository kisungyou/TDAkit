#' \eqn{k}-Groups Clustering of Multiple Functional Summaries by Energy Distance
#' 
#' Given \eqn{N} functional summaries \eqn{\Lambda_1 (t), \Lambda_2 (t), \ldots, \Lambda_N (t)}, 
#' perform \eqn{k}-groups clustering by energy distance using \eqn{L_2} metric.
#' 
#' @param fslist a length-\eqn{N} list of functional summaries of persistent diagrams.
#' @param k the number of clusters.
#' @param ... extra parameters including \describe{
#'   \item{maxiter}{the number of iterations (default: 50).}
#'   \item{nstart}{the number of restarts (default: 2).}
#' }
#' 
#' @return a length-\eqn{N} vector of class labels (from \eqn{1:k}).
#' 
#' @examples 
#' \donttest{
#' # ---------------------------------------------------------------------------
#' #           K-Groups Clustering via Energy Distance
#' #
#' # We will cluster dim=0 under top-5 landscape functions with 
#' # - Class 1 : 'iris' dataset with noise
#' # - Class 2 : samples from 'gen2holes()'
#' # - Class 3 : samples from 'gen2circles()'
#' # ---------------------------------------------------------------------------
#' ## Generate Data and Diagram from VR Filtration
#' ndata     = 10
#' list_rips = list()
#' for (i in 1:ndata){
#'   dat1 = as.matrix(iris[,1:4]) + matrix(rnorm(150*4), ncol=4)
#'   dat2 = gen2holes(n=100, sd=1)$data
#'   dat3 = gen2circles(n=100, sd=1)$data
#'   
#'   list_rips[[i]] = diagRips(dat1, maxdim=1)
#'   list_rips[[i+ndata]] = diagRips(dat2, maxdim=1)
#'   list_rips[[i+(2*ndata)]] = diagRips(dat3, maxdim=1)
#' }
#' 
#' ## Compute Persistence Landscapes from Each Diagram with k=5 Functions
#' list_land0 = list()
#' for (i in 1:(3*ndata)){
#'   list_land0[[i]] = diag2landscape(list_rips[[i]], dimension=0, k=5)
#' }
#' 
#' ## Run K-Groups Clustering with different K's
#' label2  = fskgroups(list_land0, k=2)
#' label3  = fskgroups(list_land0, k=3)
#' label4  = fskgroups(list_land0, k=4)
#' truelab = rep(c(1,2,3), each=ndata)
#' 
#' ## Run MDS & Visualization
#' embed = fsmds(list_land0, ndim=2)
#' opar  = par(no.readonly=TRUE)
#' par(mfrow=c(2,2), pty="s")
#' plot(embed, col=truelab, pch=19, main="true label")
#' plot(embed, col=label2,  pch=19, main="k=2 label")
#' plot(embed, col=label3,  pch=19, main="k=3 label")
#' plot(embed, col=label4,  pch=19, main="k=4 label")
#' par(opar)
#' }
#' 
#' @concept summaries
#' @export
fskgroups <- function(fslist, k=2, ...){
  ## PREPROCESSING
  dtype    = check_list_summaries("fskgroups", fslist)
  myk      = max(1, round(k))
  params   = list(...)
  pnames   = names(params)
  myiter   = max(20, ifelse(("maxiter"%in%pnames), params$maxiter, 50))
  mynstart = max(2, ifelse(("nstart"%in%pnames), params$nstart, 2))
  
  ## COMPUTE 2-NORM
  dmat = TDAkit::fsdist(fslist, p=2, as.dist = TRUE)
  
  ## APPLY ENERGY'S K-GROUPS
  output = energy::kgroups(dmat, myk, nstart=mynstart, iter.max=myiter)
  return(as.vector(output$cluster))
}