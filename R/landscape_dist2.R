#' Pairwise \eqn{L_p} Distances for Two Sets of Persistence Landscapes
#' 
#' Given two sets of landscape functions \eqn{\Lambda_1 (t),\ldots, \Lambda_M (t)} and 
#' \eqn{\tilde{\Lambda}_1 (t), \ldots, \tilde{\Lambda}_N (t)}, compute \eqn{L_p} distance 
#' across pairs.
#' 
#' @param dlist1 a length-\eqn{M} list of \code{'landscape'} objects from \code{\link{diag2landscape}}.
#' @param dlist2 a length-\eqn{N} list of \code{'landscape'} objects from \code{\link{diag2landscape}}.
#' @param p an exponent in \eqn{[1,\infty)} (default: 2).
#' 
#' @return an \eqn{(M\times N)} distance matrix.
#' 
#' @examples 
#' # ---------------------------------------------------------------------------
#' #         Compute L1 and L2 Distance for Two Sets of Landscapes
#' #
#' # First  set consists of {Class 1, Class 2}, while
#' # Second set consists of {Class 1, Class 3} where
#' #
#' # - Class 1 : 'iris' dataset with noise
#' # - Class 2 : samples from 'gen2holes()'
#' # - Class 3 : samples from 'gen2circles()'
#' # ---------------------------------------------------------------------------
#' ## Generate Data and Diagram from VR Filtration
#' ndata      = 10
#' list_rips1 = list()
#' list_rips2 = list()
#' for (i in 1:ndata){
#'   dat1 = as.matrix(iris[,1:4]) + matrix(rnorm(150*4, sd=4), ncol=4)
#'   dat2 = gen2holes(n=100, sd=1)$data
#'   dat3 = as.matrix(iris[,1:4]) + matrix(rnorm(150*4, sd=4), ncol=4)
#'   dat4 = gen2circles(n=100, sd=1)$data
#'   
#'   list_rips1[[i]]       = diagRips(dat1, maxdim=1)
#'   list_rips1[[i+ndata]] = diagRips(dat2, maxdim=1)
#'   
#'   list_rips2[[i]]       = diagRips(dat3, maxdim=1)
#'   list_rips2[[i+ndata]] = diagRips(dat4, maxdim=1)
#' }
#' 
#' ## Compute Persistence Landscapes from Each Diagram with k=10 Functions
#' #  We try to get distance in dimension 1 only for faster comparison.
#' list_pset1 = list()
#' list_pset2 = list()
#' for (i in 1:(2*ndata)){
#'   list_pset1[[i]] = diag2landscape(list_rips1[[i]], dimension=1, k=10)
#'   list_pset2[[i]] = diag2landscape(list_rips2[[i]], dimension=1, k=10)
#' }
#' 
#' ## Compute L1 and L2 Distance Matrix
#' dmat1 = pldist2(list_pset1, list_pset2, p=1)
#' dmat2 = pldist2(list_pset1, list_pset2, p=2)
#' 
#' ## Visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' image(dmat1[,(2*ndata):1], axes=FALSE, main="distance for p=1")
#' image(dmat2[,(2*ndata):1], axes=FALSE, main="distance for p=2")
#' par(opar)
#' 
#' @concept landscape
#' @export
pldist2 <- function(dlist1, dlist2, p=2){
  ## PREPROCESSING
  if (!check_list_landscape(dlist1)){
    stop("* pldist2 : 'dlist1' should be a list of 'landscape' objects.")
  }
  if (!check_list_landscape(dlist2)){
    stop("* pldist2 : 'dlist2' should be a list of 'landscape' objects.")
  }
  myp    = max(1, as.double(p))
  dlist3 = c(dlist1, dlist2)
  m = length(dlist1)
  n = length(dlist2)
  
  ## MAIN COMPUTATION
  #  Rearrange as a 3d array
  list3d = adjust_list_landscapes(dlist3, as.list=FALSE)
  mytseq = as.vector(list3d$tseq)
  #  Iterate
  output = array(0,c(m,n))
  for (i in 1:m){
    tgti = list3d$array3d[,,i]
    for (j in 1:n){
      tgtj = list3d$array3d[,,(m+j)]
      output[i,j] = landscape.norm.internal(tgti-tgtj, mytseq, myp)
    }
  }
  return(output)
}