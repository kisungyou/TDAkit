#' Convert Persistence Diagram into Persistent Silhouette
#' 
#' Persistence Silhouette (PS) is a functional summary of persistent homology 
#' that is constructed given a \code{homology} object. PS is a weighted average of 
#' landscape functions so that it becomes a uni-dimensional function.
#' 
#' @param homology an object of S3 class \code{"homology"} generated from \code{diagRips} or other diagram-generating functions.
#' @param dimension dimension of features to be considered (default: 1).
#' @param p an exponent for the weight function of form \eqn{|a-b|^p} (default: 2).
#' @param nseq grid size for which the landscape function is evaluated. 
#' 
#' @return a list object of \code{"silhouette"} class containing\describe{
#' \item{lambda}{an \eqn{(\code{nseq} \times k)} landscape functions.}
#' \item{tseq}{a length-\code{nseq} vector of domain grid.}
#' \item{dimension}{dimension of features considered.}
#' }
#' 
#' @examples 
#' # ---------------------------------------------------------------------------
#' #              Persistence Silhouette of 'iris' Dataset
#' #
#' # We will extract silhouettes of dimensions 0, 1, and 2.
#' # ---------------------------------------------------------------------------
#' ## Prepare 'iris' data
#' XX = as.matrix(iris[,1:4])
#' 
#' ## Compute Persistence Diagram 
#' pdrips = diagRips(XX, maxdim=2)
#' 
#' ## Convert to Silhouettes of Each Dimension
#' sil0 <- diag2silhouette(pdrips, dimension=0)
#' sil1 <- diag2silhouette(pdrips, dimension=1)
#' sil2 <- diag2silhouette(pdrips, dimension=2)
#' 
#' ## Visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,2))
#' plot(pdrips$Birth, pdrips$Death, col=as.factor(pdrips$Dimension),
#'      pch=19, main="persistence diagram", xlab="Birth", ylab="Death")
#' plot(sil0$tseq, sil0$lambda, type="l", lwd=3, main="dimension 0", xlab="t")
#' plot(sil1$tseq, sil1$lambda, type="l", lwd=3, main="dimension 1", xlab="t")
#' plot(sil2$tseq, sil2$lambda, type="l", lwd=3, main="dimension 2", xlab="t")
#' par(opar)
#' 
#' @concept convert
#' @export
diag2silhouette <- function(homology, dimension=1, p=2, nseq=100){
  ## PREPROCESSING
  #  homology
  if (!inherits(homology,"homology")){
    stop("* diag2silhouette : input 'homology' is not a valid homology object. Please use an output from 'diagRips' or other construction algorithms.")
  }
  #  dimension
  dimension = round(dimension) # target dimension : cannot be multiple
  if (!(dimension %in% homology$Dimension)){
    stop("* diag2silhouette : input 'dimension' does not have corresponding information in the given 'homology'.")
  }
  idin      = which(homology$Dimension==dimension)
  dat.dim   = round(homology$Dimension[idin])
  dat.birth = homology$Birth[idin]
  dat.death = homology$Death[idin]
  #  weight function
  myp = max(1, as.double(p))
  myt = seq(from=0, to=max(dat.death), length.out=max(10, round(nseq))) # time sequence part  
  
  ## MAIN COMPUTATION
  #  We use weight function |a-b|^p
  out.numerator   = rep(0,length(myt))
  out.denominator = 0
  for (i in 1:length(dat.birth)){
    # select birth and death pair
    bj = dat.birth[i]
    dj = dat.death[i]
    
    # update denominator
    weightj = (base::abs(bj-dj)^myp)
    out.denominator = out.denominator + weightj
    # update numerator
    out.numerator = out.numerator + (weightj*base::pmax(base::pmin(myt-bj, dj-myt), rep(0,length(myt))))
  }
  # finalize
  lbdfun = as.vector((out.numerator/out.denominator))
  
  ## RETURN
  res = list(lambda=lbdfun, tseq=myt, dimension=dimension)
  class(res) = "silhouette"
  return(res)
}