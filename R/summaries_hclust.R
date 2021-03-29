#' Hierarchical Agglomerative Clustering
#' 
#' Given multiple functional summaries \eqn{\Lambda_1 (t), \Lambda_2 (t), \ldots, \Lambda_N (t)}, 
#' perform hierarchical agglomerative clustering with \eqn{L_2} distance.
#' 
#' @param fslist a length-\eqn{N} list of functional summaries of persistent diagrams.
#' @param method agglomeration method to be used. This must be one of \code{"single"}, \code{"complete"}, \code{"average"}, \code{"mcquitty"}, \code{"ward.D"}, \code{"ward.D2"}, \code{"centroid"} or \code{"median"}.
#' @param members \code{NULL} or a vector whose length equals the number of observations. See \code{\link[stats]{hclust}} for details.
#' 
#' @return an object of class \code{hclust}. See \code{\link[stats]{hclust}} for details. 
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
#' list_lab = c(rep(1,ndata), rep(2,ndata), rep(3,ndata))
#' 
#' ## Compute Persistence Landscapes from Each Diagram with k=5 Functions
#' list_land0 = list()
#' for (i in 1:(3*ndata)){
#'   list_land0[[i]] = diag2landscape(list_rips[[i]], dimension=0, k=5)
#' }
#' 
#' ## Run MDS for Visualization
#' embed = fsmds(list_land0, ndim=2)
#' 
#' ## Clustering with 'single' and 'complete' linkage
#' hc.sing <- fshclust(list_land0, method="single")
#' hc.comp <- fshclust(list_land0, method="complete")
#' 
#' ## Visualize
#' opar  = par(no.readonly=TRUE)
#' par(mfrow=c(1,3))
#' plot(embed, pch=19, col=list_lab, main="2-dim embedding")
#' plot(hc.sing, main="single linkage")
#' plot(hc.comp, main="complete linkage")
#' par(opar)
#' }
#' 
#' @concept summaries
#' @export
fshclust <- function(fslist, method = c("single", "complete", "average", "mcquitty", "ward.D", "ward.D2",
                                        "centroid", "median"), members=NULL){
  ## PREPROCESSING
  mymethod  = match.arg(method)
  mymembers = members
  
  ## PAIRWISE DISTANCE
  if (inherits(fslist, "dist")){
    pdmat = fslist
  } else {
    dtype = check_list_summaries("fseqdist", fslist)
    pdmat = TDAkit::fsdist(fslist, p=2, as.dist=TRUE)
  }
  
  
  ## RUN VIA MAOTAI
  fimport = utils::getFromNamespace("hidden_hclust", "maotai")
  hcout   = fimport(pdmat, mymethod, mymembers)
  return(hcout)
}