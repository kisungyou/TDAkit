#' Persistence Landscape Kernel 
#' 
#' Given multiple persistence landscapes \eqn{\Lambda_1 (t), \Lambda_2 (t), \ldots, \Lambda_N (t)}, compute 
#' the persistence landscape kernel under the \eqn{L_2} sense. 
#' 
#' @param landlist a length-\eqn{N} list of \code{"landscape"} objects, which can be obtained from \code{\link{diag2landscape}} function.
#' 
#' @return an \eqn{(N\times N)} kernel matrix.
#' 
#' @examples 
#' \donttest{
#' # ---------------------------------------------------------------------------
#' #      Persistence Landscape Kernel in Dimension 0 and 1
#' #
#' # We will compare dim=0,1 with top-20 landscape functions with 
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
#' #  We try to get distance in dimensions 0 and 1.
#' list_land0 = list()
#' list_land1 = list()
#' for (i in 1:(3*ndata)){
#'   list_land0[[i]] = diag2landscape(list_rips[[i]], dimension=0, k=5)
#'   list_land1[[i]] = diag2landscape(list_rips[[i]], dimension=1, k=5)
#' }
#' 
#' ## Compute Persistence Landscape Kernel Matrix
#' plk0 <- plkernel(list_land0)
#' plk1 <- plkernel(list_land1)
#' 
#' ## Visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' image(plk0[,(3*(ndata)):1], axes=FALSE, main="Kernel : dim=0")
#' image(plk1[,(3*(ndata)):1], axes=FALSE, main="Kernel : dim=1")
#' par(opar)
#' }
#' 
#' @references 
#' Jan Reininghaus, Stefan Huber, Ulrich Bauer, and Roland Kwitt (2015). ``A stable multi-scale kernel for topological machine learning.'' 
#' \emph{Proc. 2015 IEEE Conf. Comp. Vision & Pat. Rec. (CVPR ’15)}.
#' 
#' @concept landscape
#' @export
plkernel <- function(landlist){
  ## PREPROCESSING
  type.obsolete = check_list_summaries("plkernel", landlist)
  if (!inherits(landlist[[1]], "landscape")){
    stop("* plkernel : input 'landlist' should be a list of 'landscape' objects.")
  }
  
  ## MAIN COMPUTATION
  #  Rearrange as a 3d array (T,K,nlist) & list$tseq
  nlist  = length(landlist)
  list3d = adjust_list_landscapes(landlist, as.list=FALSE)
  mytseq = as.vector(list3d$tseq)
  #  Iterate
  output = array(0,c(nlist, nlist))
  for (i in 1:(nlist-1)){
    tgti = as.matrix(list3d$array3[,,i])
    for (j in (i+1):nlist){
      tgtj = as.matrix(list3d$array3[,,j])
      output[i,j] <- plkernel_single(mytseq, tgti, tgtj)
      output[j,i] <- output[i,j]
    }
  }
  return(output)
}

#' @keywords internal
#' @noRd
plkernel_single <- function(tseq, mat1, mat2){
  K = base::ncol(mat1)
  output = 0
  for (k in 1:K){
    yseq   = as.vector(mat1[,k])*as.vector(mat2[,k])
    output = output + simple_integral_1d(yseq, tseq)
  }
  return(output)
}