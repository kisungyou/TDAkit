#' \eqn{L_p} norm of a single Persistence Landscape
#' 
#' Given a persistence landscape \eqn{\Lambda (t)}, compute the \eqn{p}-norm 
#' \deqn{\| \Lambda (t) \|_p = \left( \sum_{k=1}^K \| \Lambda_{k} (t) \|_p^p \right)^{1/p}} 
#' for top \eqn{K} landscape functions.
#' 
#' @param landscape a \code{'landscape'} object from \code{\link{diag2landscape}}.
#' @param p an exponent in \eqn{[1,\infty)} (default: 2).
#' 
#' @examples 
#' ## Generate Toy Data from 'gen2circles()'
#' dat = gen2circles(n=100)$data
#' 
#' ## Compute PD and PL
#' myPD  = diagRips(dat, maxdim=1)
#' myPL0 = diag2landscape(myPD, dim=0)
#' myPL1 = diag2landscape(myPD, dim=1)
#' 
#' ## Compute 2-norm
#' plnorm(myPL0, p=2)
#' plnorm(myPL1, p=2)
#' 
#' @references 
#' \insertRef{bubenik_persistence_2018}{TDAkit}
#' 
#' @concept landscape
#' @export
plnorm <- function(landscape, p=2){
  ## PREPROCESSING
  if (!inherits(landscape, "landscape")){
    stop("* plnorm : provide a compatiable landscape object from 'diag2landscape'.")
  }
  myp = max(1, as.double(p))
  
  ## MAIN COMPUTATION
  if (is.infinite(myp)){
    return(max(landscape$lambda))
  } else {
    trimmed = as.matrix(abs(trim.the.norm(landscape$lambda)))^myp
    return(simple_integral(trimmed, landscape$tseq)^(1/myp))
  }
}