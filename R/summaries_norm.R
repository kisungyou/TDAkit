#' \eqn{L_p} Norm of a Single Functional Summary
#' 
#' Given a functional summary \eqn{\Lambda (t)}, compute the \eqn{p}-norm.
#' 
#' @param fsobj a functional summary object.
#' @param p an exponent in \eqn{[1,\infty)} (default: 2).
#' 
#' @return an \eqn{L_p}-norm value.
#' 
#' @examples 
#' \donttest{
#' ## Generate Toy Data from 'gen2circles()'
#' dat = gen2circles(n=100)$data
#' 
#' ## Compute PD, Landscapes, and Silhouettes
#' myPD  = diagRips(dat, maxdim=1)
#' myPL0 = diag2landscape(myPD, dimension=0)
#' myPL1 = diag2landscape(myPD, dimension=1)
#' myPS0 = diag2silhouette(myPD, dimension=0)
#' myPS1 = diag2silhouette(myPD, dimension=1)
#' 
#' ## Compute 2-norm
#' fsnorm(myPL0, p=2)
#' fsnorm(myPL1, p=2)
#' fsnorm(myPS0, p=2)
#' fsnorm(myPS1, p=2)
#' }
#' 
#' @concept summaries
#' @export
fsnorm <- function(fsobj, p=2){
  ## PREPROCESS
  myp   = max(1, as.double(p))
  clobj = class(fsobj)
  
  ## BRANCHING CASE
  output = switch(clobj, 
                  "silhouette" = fsnorm_sils(fsobj, myp),
                  "landscape"  = fsnorm_land(fsobj, myp))
  return(output)
}


# individual functions ----------------------------------------------------
#' @keywords internal
#' @noRd
fsnorm_sils <- function(fsobj, p){
  if (is.infinite(p)){
    return(max(fsobj$lambda))
  } else {
    myeps   = 10*.Machine$double.eps
    trimmed = matrix(abs(trim.the.norm(fsobj$lambda))^p, ncol=1)
    return(simple_integral(trimmed, fsobj$tseq)^(1/p))
  }
}
#' @keywords internal
#' @noRd
fsnorm_land <- function(fsobj, myp){
  ## MAIN COMPUTATION
  if (is.infinite(myp)){
    return(max(fsobj$lambda))
  } else {
    trimmed = as.matrix(abs(trim.the.norm(fsobj$lambda)))^myp
    return(simple_integral(trimmed, fsobj$tseq)^(1/myp))
  }
}
