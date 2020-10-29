#' Compute Alpha Shape Persistence Diagram in 3D
#' 
#' \code{diagAS} computes the persistent diagram of the Alpha Shape filtration 
#' constructed on a point cloud in \eqn{\mathbf{R}^3}. This function is a 
#' wrapper to \pkg{TDA}'s implementation. It is noted that AS filtration is \emph{only} 
#' available in 3-dimensional point cloud.
#' 
#' @param data an \eqn{(n\times 3)} data matrix.
#' @param maxdim maximum dimension of the computed homological features (default: 1).
#' @param threshold maximum value of the filtration (default: \code{Inf}). 
#' 
#' @return a data frame with following columns\describe{
#' \item{Dimension}{dimension corresponding to a feature.}
#' \item{Birth}{birth of a feature.}
#' \item{Death}{death of a feature.}
#' }
#' 
#' @examples 
#' # ---------------------------------------------------------------------------
#' #                 Compare VR and AS Persistence Diagrams
#' # ---------------------------------------------------------------------------
#' # Use 'iris' data for the first 3 columns
#' XX = as.matrix(iris[,1:3])
#' 
#' # Compute VR and AS Diagram
#' run.vr = diagRips(XX, maxdim=1)
#' run.as = diagAS(XX, maxdim=1)
#' 
#' col1 = as.factor(run.vr$Dimension)
#' col2 = as.factor(run.as$Dimension)
#' 
#' # Visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2), pty="s")
#' plot(run.vr$Birth, run.vr$Death, pch=19, col=col1, main="VR diagram")
#' plot(run.as$Birth, run.as$Death, pch=19, col=col2, main="AS diagram")
#' par(opar)
#' 
#' @seealso \code{\link[TDA]{alphaShapeDiag}}
#' 
#' @concept diagram
#' @export
diagAS <- function(data, maxdim=1){
  ## PREPARATION
  if (!is.matrix(data)){
    stop(paste0("* diagAS : input 'data' must be a matrix."))
  }
  if (ncol(data)!=3){
    stop("* diagAS : input 'data' must have 3 columns.")
  }
  mydim = max(0, round(maxdim))
  mydim = min(2, mydim)
  
  ## COMPUTATION
  tdarun = TDA::alphaShapeDiag(data, maxdimension = mydim, library = c("GUDHI", "Dionysus"))$diagram
  delete = which(is.infinite(tdarun[,3]))
  tdarun = tdarun[-delete,]
  
  ## WRAP
  output = data.frame(Dimension = round(tdarun[,1]),
                      Birth = tdarun[,2],
                      Death = tdarun[,3])
  return(output)
}
