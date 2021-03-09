#' Convert Persistence Diagram into Persistence Landscape
#' 
#' Persistence Landscape (PL) is a functional summary of persistent homology 
#' that is constructed given a \code{homology} object. 
#' 
#' @param homology an object of S3 class \code{"homology"} generated from \code{diagRips} or other homology-generating functions.
#' @param dimension dimension of features to be considered (default: 1).
#' @param k the number of top landscape functions to be used (default: 0). When \code{k=0} is set, it gives all relevant landscape functions that are non-zero.
#' @param nseq grid size for which the landscape function is evaluated (default: 1000).
#' 
#' @return a list object of \code{"landscape"} class containing\describe{
#' \item{lambda}{an \eqn{(\code{nseq} \times k)} landscape functions.}
#' \item{tseq}{a length-\code{nseq} vector of domain grid.}
#' \item{dimension}{dimension of features considered.}
#' }
#' 
#' @examples 
#' \donttest{
#' # ---------------------------------------------------------------------------
#' #              Persistence Landscape of 'iris' Dataset
#' #
#' # We will extract landscapes of dimensions 0, 1, and 2.
#' # For each feature, only the top 5 landscape functions are plotted.
#' # ---------------------------------------------------------------------------
#' ## Prepare 'iris' data
#' XX = as.matrix(iris[,1:4])
#' 
#' ## Compute Persistence Diagram 
#' pdrips = diagRips(XX, maxdim=2)
#' 
#' ## Convert to Landscapes of Each Dimension
#' land0 <- diag2landscape(pdrips, dimension=0, k=5)
#' land1 <- diag2landscape(pdrips, dimension=1, k=5)
#' land2 <- diag2landscape(pdrips, dimension=2, k=5)
#' 
#' ## Visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,2))
#' plot(pdrips$Birth, pdrips$Death, col=as.factor(pdrips$Dimension),
#'      pch=19, main="persistence diagram", xlab="Birth", ylab="Death")
#' matplot(land0$tseq, land0$lambda, type="l", lwd=3, main="dimension 0", xlab="t")
#' matplot(land1$tseq, land1$lambda, type="l", lwd=3, main="dimension 1", xlab="t")
#' matplot(land2$tseq, land2$lambda, type="l", lwd=3, main="dimension 2", xlab="t")
#' par(opar)
#' }
#' @references 
#' Peter Bubenik (2018). ``The Persistence Landscape and Some of Its Properties.'' \emph{arXiv:1810.04963}.
#' 
#' @concept convert
#' @export
diag2landscape <- function(homology, dimension=1, k=0, nseq=1000){
  ## PREPROCESSING
  #  homology
  if (!inherits(homology,"homology")){
    stop("* diag2landscape : input 'homology' is not a valid homology object. Please use an output from 'diagRips' or other construction algorithms.")
  }
  #  dimension
  dimension = round(dimension) # target dimension : cannot be multiple
  if (!(dimension %in% homology$Dimension)){
    stop("* diag2landscape : input 'dimension' does not have corresponding information in the given 'homology'.")
  }
  idin      = which(homology$Dimension==dimension)
  dat.dim   = round(homology$Dimension[idin])
  dat.birth = homology$Birth[idin]
  dat.death = homology$Death[idin]
  #  others
  myk = ifelse((round(k)<1), length(dat.birth), round(k))               # number of landscape functions
  myt = seq(from=0, to=max(dat.death), length.out=max(10, round(nseq))) # time sequence part  

  ## MAIN COMPUTATION
  #  if too many features are expected, reduce the number
  if (length(dat.birth) < myk){
    myk = length(dat.birth)
  }
  #  main compute !
  lambdas = compute_lambdas(myt, dat.birth, dat.death, myk)
  
  ## WRAP AND RETURN
  res = list(lambda=lambdas, tseq=myt, dimension=dimension)
  class(res) = "landscape"
  return(res)
}

#' @keywords internal
#' @noRd
compute_lambdas <- function(tseq, births, deaths, maxK){
  ntest = length(births)
  ntime = length(tseq)
  
  output = array(0,c(ntime,ntest))
  for (i in 1:ntest){
    b = births[i]
    d = deaths[i]
    output[,i] = base::pmax(base::pmin(tseq-b, d-tseq), rep(0,ntime))
  }
  
  newout = array(0,c(ntime,ntest))
  for (i in 1:ntime){
    newout[i,] = sort(output[i,], decreasing = TRUE)
  }
  return((newout[,1:maxK]))
}