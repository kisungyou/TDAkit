#' t-distributed Stochastic Neighbor Embedding
#' 
#' Given \eqn{N}  functional summaries \eqn{\Lambda_1 (t), \Lambda_2 (t), \ldots, \Lambda_N (t)}, 
#' t-SNE mimicks the pattern of probability distributions over pairs of Banach-valued 
#' objects on low-dimensional target embedding space by minimizing Kullback-Leibler divergence.
#' 
#' @param fslist a length-\eqn{N} list of functional summaries of persistent diagrams.
#' @param ndim an integer-valued target dimension.
#' @param ... extra parameters for \code{\link[Rtsne]{Rtsne}} algorithm, such as perplexity, momentum, and others.
#' 
#' @return a named list containing \describe{
#' \item{embed}{an \eqn{(N\times ndim)} matrix whose rows are embedded observations.}
#' \item{stress}{discrepancy between embedded and original distances as a measure of error.}
#' }
#' 
#' @examples
#' \donttest{
#' # ---------------------------------------------------------------------------
#' #     Multidimensional Scaling for Multiple Landscapes and Silhouettes
#' #
#' # We will compare dim=0 with top-5 landscape and silhouette functions with 
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
#' ## Compute Landscape and Silhouettes of Dimension 0
#' list_land = list()
#' list_sils = list()
#' for (i in 1:(3*ndata)){
#'   list_land[[i]] = diag2landscape(list_rips[[i]], dimension=0)
#'   list_sils[[i]] = diag2silhouette(list_rips[[i]], dimension=0)
#' }
#' list_lab = rep(c(1,2,3), each=ndata)
#' 
#' ## Run t-SNE and Classical/Metric MDS
#' land_cmds = fsmds(list_land, method="classical")
#' land_mmds = fsmds(list_land, method="metric")
#' land_tsne = fstsne(list_land, perplexity=5)$embed
#' sils_cmds = fsmds(list_sils, method="classical")
#' sils_mmds = fsmds(list_sils, method="metric")
#' sils_tsne = fstsne(list_land, perplexity=5)$embed
#' 
#' ## Visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,3))
#' plot(land_cmds, pch=19, col=list_lab, main="Landscape+CMDS")
#' plot(land_mmds, pch=19, col=list_lab, main="Landscape+MMDS")
#' plot(land_tsne, pch=19, col=list_lab, main="Landscape+tSNE")
#' plot(sils_cmds, pch=19, col=list_lab, main="Silhouette+CMDS")
#' plot(sils_mmds, pch=19, col=list_lab, main="Silhouette+MMDS")
#' plot(sils_tsne, pch=19, col=list_lab, main="Silhouette+tSNE")
#' par(opar)
#' }
#' 
#' @seealso \code{\link[Rtsne]{Rtsne}}
#' 
#' @concept summaries
#' @export
fstsne <- function(fslist, ndim=2, ...){
  ## PREPROCESSING
  dtype  = check_list_summaries("fstsne", fslist)
  myndim = max(2, round(ndim))
  
  ## PAIRWISE DISTANCE
  distobj = TDAkit::fsdist(fslist, p=2, as.dist=TRUE)
  
  ## RUN VIA MAOTAI
  func.import = utils::getFromNamespace("hidden_tsne", "maotai")
  out.tsne    = func.import(distobj, ndim=myndim, ...)
  return(out.tsne)
}