#' Spectral Clustering by Zelnik-Manor and Perona (2005)
#' 
#' Given \eqn{N} functional summaries \eqn{\Lambda_1 (t), \Lambda_2 (t), \ldots, \Lambda_N (t)}, 
#' perform spectral clustering proposed by Zelnik-Manor and Perona using a set of 
#' data-driven bandwidth parameters.
#' 
#' @param fslist a length-\eqn{N} list of functional summaries of persistent diagrams.
#' @param k the number of cluster (default: 2).
#' @param nnbd neighborhood size to define data-driven bandwidth parameter (default: 5).
#' 
#' @return a length-\eqn{N} vector of class labels (from \eqn{1:k}).
#' 
#' @examples 
#' \donttest{
#' # ---------------------------------------------------------------------------
#' #           Spectral Clustering Clustering via Energy Distance
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
#' ## Run Spectral Clustering using Different K's.
#' label2  = fssc05Z(list_land0, k=2)
#' label3  = fssc05Z(list_land0, k=3)
#' label4  = fssc05Z(list_land0, k=4)
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
#' @references  
#' Zelnik-manor L, Perona P (2005). ``Self-Tuning Spectral Clustering.'' In Saul LK, Weiss Y, Bottou L (eds.), \emph{Advances in Neural Information Processing Systems 17}, 1601–1608. MIT Press.
#' 
#' @concept summaries
#' @export
fssc05Z <- function(fslist, k=2, nnbd=5){
  ## PREPROCESSING
  myk    = max(1, round(k))
  mynnbd = max(1, round(nnbd))
  
  ## PAIRWISE DISTANCE
  if (inherits(fslist, "dist")){
    pdmat = fslist
  } else {
    dtype = check_list_summaries("fssc05Z", fslist)
    pdmat = TDAkit::fsdist(fslist, p=2, as.dist=TRUE)
  }
  
  ## RUN CLUSTERING
  run_sc = T4cluster::sc05Z(pdmat, k=myk, nnbd=mynnbd)
  output = as.vector(run_sc$cluster)
  return(output)
}