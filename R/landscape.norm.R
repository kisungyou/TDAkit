#' \eqn{L_p} norm of a single Persistence Landscape
#' 
#' @export
landscape.norm <- function(landscape, p=2){
  #############################################
  # Preprocessing : checkers
  if (!inherits(landscape, "kit.landscape")){
    stop("* landscape.norm : provide a compatiable landscape object from 'd2landscape'.")
  }
  if (p < 1){
    stop("* landscape.norm : 'p'-norm is defined for the range [1,Inf]. ")
  }
  
  #############################################
  # Main Computation with Branching
  if (is.infinite(p)){
    return(max(landscape$lambda))
  } else {
    myeps   = 10*.Machine$double.eps
    trimmed = as.matrix(abs(trim.the.norm(landscape$lambda)))^p
    return(simple_integral(trimmed, landscape$tseq)^(1/p))
  }
}

# simply, without kit, internal use only ----------------------------------
#' @keywords internal
#' @noRd
landscape.norm.internal <- function(lmat, tseq, p=2){
  if (is.infinite(p)){
    return(max(lmat))
  } else {
    myeps   = .Machine$double.eps
    trimmed = as.matrix(abs(trim.the.norm(lmat)))^p
    return(simple_integral(trimmed, tseq)^(1/p))
  }
}


# trim the zero ones ------------------------------------------------------
#' @keywords internal
#' @noRd
trim.the.norm <- function(lmat){
  if ((is.vector(lmat))||(ncol(lmat)==1)){
    return(as.matrix(lmat))
  } else {
    check.small <- function(x){
      return(all(x < 10*.Machine$double.eps))
    }
    lastid = suppressWarnings(min(which(apply(lmat, 2, check.small) == TRUE)) - 1)
    if (is.infinite(lastid)){
      lastid = 0
    }
    return(lmat[,1:max(lastid,1)])  
  }
}
