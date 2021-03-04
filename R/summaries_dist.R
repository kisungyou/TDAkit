#' Pairwise \eqn{L_p} Distance of Multiple Functional Summaries
#' 
#' Given multiple functional summaries \eqn{\Lambda_1 (t), \Lambda_2 (t), \ldots, \Lambda_N (t)}, 
#' compute \eqn{L_p} distance in a pairwise sense.
#' 
#' @param fslist a length-\eqn{N} list of functional summaries of persistent diagrams.
#' @param p an exponent in \eqn{[1,\infty)} (default: 2).
#' @param as.dist logical; if TRUE, it returns \code{dist} object, else it returns an \eqn{(N\times N)} symmetric matrix.
#' 
#' @return a S3 \code{dist} object or \eqn{(N\times N)} symmetric matrix of pairwise distances according to \code{as.dist} parameter.
#' 
#' @examples 
#' \donttest{
#' # ---------------------------------------------------------------------------
#' #      Compute L_2 Distance for 3 Types of Landscapes and Silhouettes
#' #
#' # We will compare dim=0,1 with top-5 landscape functions with 
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
#' ## Compute Silhouettes
#' list_sil0 = list()
#' list_sil1 = list()
#' for (i in 1:(3*ndata)){
#'   list_sil0[[i]] = diag2silhouette(list_rips[[i]], dimension=0)
#'   list_sil1[[i]] = diag2silhouette(list_rips[[i]], dimension=1)
#' }
#' 
#' ## Compute L2 Distance Matrices
#' ldmat0 = fsdist(list_land0, p=2, as.dist=FALSE)
#' ldmat1 = fsdist(list_land1, p=2, as.dist=FALSE)
#' sdmat0 = fsdist(list_sil0, p=2, as.dist=FALSE)
#' sdmat1 = fsdist(list_sil1, p=2, as.dist=FALSE)
#' 
#' ## Visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,2), pty="s")
#' image(ldmat0[,(3*(ndata)):1], axes=FALSE, main="Landscape : dim=0")
#' image(ldmat1[,(3*(ndata)):1], axes=FALSE, main="Landscape : dim=1")
#' image(sdmat0[,(3*(ndata)):1], axes=FALSE, main="Silhouette : dim=0")
#' image(sdmat1[,(3*(ndata)):1], axes=FALSE, main="Silhouette : dim=1")
#' par(opar)
#' }
#' 
#' @concept summaries
#' @export
fsdist <- function(fslist, p=2, as.dist=TRUE){
  ## PREPROCESSING
  dtype  = check_list_summaries("fsdist", fslist)
  myp    = max(1, as.double(p))
  mydist = as.logical(as.dist)
  
  ## CASE BRANCHING
  output = switch(dtype,
                  "landscape"  = fsdist_land(fslist, myp),
                  "silhouette" = fsdist_sils(fslist, myp))
  
  ## RETURN
  if (mydist){
    return(stats::as.dist(output))
  } else {
    return(output)
  }
}


# individual functions ----------------------------------------------------
#' @keywords internal
#' @noRd
fsdist_land <- function(fslist, myp){
  ## MAIN COMPUTATION
  #  Rearrange as a 3d array
  nlist  = length(fslist)
  list3d = adjust_list_landscapes(fslist, as.list=FALSE)
  mytseq = as.vector(list3d$tseq)
  #  Iterate
  output = array(0,c(nlist, nlist))
  for (i in 1:(nlist-1)){
    tgti = as.matrix(list3d$array3[,,i])
    for (j in (i+1):nlist){
      tgtj = as.matrix(list3d$array3[,,j])
      output[i,j] <- landscape.norm.internal(tgti-tgtj, mytseq, myp)
      output[j,i] <- output[i,j]
    }
  }
  return(output)
}
#' @keywords internal
#' @noRd
fsdist_sils <- function(slist, myp){
  nlist = length(slist)
  # Main Computation
  #   1. rearrange as a list
  list2d = adjust_list_silhouette(slist, as.list=FALSE)
  mytseq = list2d$tseq
  #   2. iterate
  myeps  = 10*.Machine$double.eps
  output = array(0,c(nlist,nlist))
  for (i in 1:(nlist-1)){
    tgti = list2d$array[,i]
    for (j in (i+1):nlist){
      tgtj = list2d$array[,j]
      if (is.infinite(myp)){
        output[i,j] <- output[j,i] <- max(abs(tgti-tgtj), na.rm=TRUE)
      } else {
        output[i,j] <- output[j,i] <- (simple_integral_1d((abs(tgti-tgtj)^myp), mytseq)^(1/myp))
      }
    }
  }
  return(output)
}