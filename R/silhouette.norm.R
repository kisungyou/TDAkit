#' \eqn{L_p} norm of a single Persistence Silhouette
#' 
#' 
#' @export
silhouette.norm <- function(silhouette, p=2){
  #############################################
  # Preprocessing : checkers
  if (!inherits(silhouette, "kit.silhouette")){
    stop("* silhouette.norm : provide a compatiable silhouette object from 'd2silhouette'.")
  }
  if (p < 1){
    stop("* silhouette.norm : 'p'-norm is defined for the range [1,Inf]. ")
  }
  
  #############################################
  # Main Computation with Branching
  if (is.infinite(p)){
    return(max(silhouette$lambda))
  } else {
    myeps   = 10*.Machine$double.eps
    trimmed = matrix(abs(trim.the.norm(silhouette$lambda))^p, ncol=1)
    return(simple_integral(trimmed, silhouette$tseq)^(1/p))
  }
}