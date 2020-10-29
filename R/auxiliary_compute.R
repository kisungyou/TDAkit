## AUXILIARY : Computation
#  (01) trim.the.norm          : trim the zero ones
#  (02) landscape.norm.internal : simply internal use only


# (01) trim.the.norm ------------------------------------------------------
#' @keywords internal
#' @noRd
trim.the.norm <- function(lmat){
  if ((is.vector(lmat))||(ncol(lmat)==1)){
    return(matrix(lmat, ncol=1))
  } else {
    p = ncol(lmat)
    vecnorm = rep(0,p)
    for (i in 1:p){
      tgt = as.vector(lmat[,i])
      vecnorm[i] = sqrt(sum(tgt^2))
    }
    whew = which(vecnorm > sqrt(.Machine$double.eps))
    if (length(whew) < 1){
      idlarge = 1
    } else {
      idlarge = max(whew, 1)
    }
    return(lmat[,1:idlarge])  
  }
}

# (02) landscape.norm.internal --------------------------------------------
#' @keywords internal
#' @noRd
landscape.norm.internal <- function(lmat, tseq, p=2){
  if (is.infinite(p)){
    return(max(lmat))
  } else {
    trimmed = as.matrix(abs(trim.the.norm(lmat)))^p
    return(simple_integral(trimmed, tseq)^(1/p))
  }
}
