#' Pairwise \eqn{L_p} Distance for Two Sets of Functional Summaries
#' 
#' Given two sets of functional summaries \eqn{\Lambda_1 (t), \ldots, \Lambda_M (t)} and 
#' \eqn{\Omega_1 (t), \ldots, \Omega_N (t)}, compute \eqn{L_p} distance across pairs.
#' 
#' @param fslist1 a length-\eqn{M} list of functional summaries of persistent diagrams.
#' @param fslist2 a length-\eqn{N} list of functional summaries of persistent diagrams.
#' @param p an exponent in \eqn{[1,\infty)} (default: 2).
#' 
#' @return an \eqn{(M\times N)} distance matrix.
#' 
#' @examples 
#' \donttest{
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
#' dmat1 = fsdist2(list_pset1, list_pset2, p=1)
#' dmat2 = fsdist2(list_pset1, list_pset2, p=2)
#' 
#' ## Visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' image(dmat1[,(2*ndata):1], axes=FALSE, main="distance for p=1")
#' image(dmat2[,(2*ndata):1], axes=FALSE, main="distance for p=2")
#' par(opar)
#' }
#' 
#' @concept summaries
#' @export
fsdist2 <- function(fslist1, fslist2, p=2){
  ## PREPROCESSING
  dtype1  = check_list_summaries("fsdist2", fslist1)
  dtype2  = check_list_summaries("fsdist2", fslist2)
  myp    = max(1, as.double(p))
  if (!all(dtype1==dtype2)){
    stop("* fsdist2 : two inputs 'fslist1' and 'fslist2' are of different types of functional summaries.")
  }
  
  ## SWITCH CASE
  output = switch(dtype1,
                  "landscape"  = fsdist2_land(fslist1, fslist2, myp),
                  "silhouette" = fsdist2_sils(fslist1, fslist2, myp))
  return(output)
}


# individual functions ----------------------------------------------------
#' @keywords internal
#' @noRd
fsdist2_land <- function(fslist1, fslist2, myp){
  # preprocessing
  m = length(fslist1)
  n = length(fslist2)
  dlist3 = c(fslist1, fslist2)
  
  # transform
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
#' @keywords internal
#' @noRd
fsdist2_sils <- function(fslist1, fslist2, myp){
  # preprocess
  m = length(fslist1)
  n = length(fslist2)
  slists = c(fslist1, fslist2)
  
  # transform
  list2d = adjust_list_silhouette(slists, as.list=FALSE)
  funcs1 = list2d$array[,1:m]
  funcs2 = list2d$array[,(m+1):(m+n)]
  mytseq = list2d$tseq
  
  # iterate (be careful since we have merged two lists)
  output = array(0,c(m,n))
  for (i in 1:m){
    tgti = funcs1[,i]
    for (j in 1:n){
      tgtj = funcs2[,j]
      if (is.infinite(myp)){
        output[i,j] <- max(abs(tgti-tgtj))
      } else {
        trimmed = matrix(abs(trim.the.norm(tgti-tgtj))^myp, ncol=1)
        output[i,j] <- (simple_integral(trimmed, mytseq)^(1/myp))
      }
    }
  }
  
  # return
  return(output)
}