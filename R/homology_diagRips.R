#' Compute Vietoris-Rips Complex for Persistent Homology
#' 
#' \code{diagRips} computes the persistent diagram of the Vietoris-Rips filtration 
#' constructed on a point cloud represented as \code{matrix} or \code{dist} object. 
#' This function is a second-hand wrapper to \pkg{TDAstats}'s wrapping for \code{Ripser} library.
#' 
#' @param data a \code{'matrix'} or a S3 \code{'dist'} object.
#' @param maxdim maximum dimension of the computed homological features (default: 1).
#' @param threshold maximum value of the filtration (default: \code{Inf}). 
#' 
#' @return a dataframe object of S3 class \code{"homology"} with following columns\describe{
#' \item{Dimension}{dimension corresponding to a feature.}
#' \item{Birth}{birth of a feature.}
#' \item{Death}{death of a feature.}
#' }
#' 
#' @seealso \code{\link[TDAstats]{calculate_homology}}
#' 
#' @examples 
#' # ---------------------------------------------------------------------------
#' # Check consistency of two types of inputs : 'matrix' and 'dist' objects
#' # ---------------------------------------------------------------------------
#' # Use 'iris' data and compute its distance matrix
#' XX = as.matrix(iris[,1:4])
#' DX = stats::dist(XX)
#' 
#' # Compute VR Diagram with two inputs
#' vr.mat = diagRips(XX)
#' vr.dis = diagRips(DX)
#' 
#' col1 = as.factor(vr.mat$Dimension)
#' col2 = as.factor(vr.dis$Dimension)
#' 
#' # Visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' plot(vr.mat$Birth, vr.mat$Death, pch=19, col=col1, main="from 'matrix'")
#' plot(vr.dis$Birth, vr.dis$Death, pch=19, col=col2, main="from 'dist'")
#' par(opar)
#' 
#' @references
#' \insertRef{wadhwa_tdastats_2018}{TDAkit}
#' 
#' Ulrich Bauer (2019). ``Ripser: Efficient Computation of Vietoris-Rips Persistence Barcodes.'' \emph{arXiv:1908.02518}.
#' 
#' @concept diagram
#' @export
diagRips <- function(data, maxdim=1, threshold=Inf){
  ## PREPARATION
  if (is.vector(data)){
    data = matrix(data, ncol=1)
  }
  mydim = max(0, round(maxdim))
  mythr = as.double(threshold)
  if ((mythr <= 0)||(is.infinite(mythr))){
    mythr = -1
  }
  
  ## COMPUTATION WITH RIPSER
  if (inherits(data, "dist")){
    ripser_obj = TDAstats::calculate_homology(data, dim=mydim, threshold=mythr, format="distmat", return_df = TRUE)
  } else {
    ripser_obj = TDAstats::calculate_homology(data, dim=mydim, threshold=mythr, format="cloud", return_df = TRUE)
  }

  ## REARRANGE THE TERMS
  colnames(ripser_obj) = c("Dimension","Birth","Death")
  class(ripser_obj) = "homology"
  return(ripser_obj)
}