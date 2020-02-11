#' Pairwise \eqn{L_p} Distance of Multiple Persistence Landscapes
#' 
#' @export
landscape.pdist <- function(dlist, p=2, as.dist=FALSE){
  #############################################
  # Preprocessing : checkers
  if (!check_list_landscape(dlist)){
    stop("* landscape.pdist : 'dlist' is invalid. Consult with 'd2landscape' for each element in the list.")
  }
  if (p < 1){
    stop("* landscape.pdist : 'p'-norm is defined for the range [1,Inf]. ")
  }
  nlist = length(dlist)
  
  #############################################
  # Main Computation
  #   1. rearrange as a 3d array
  list3d = dlist_adjust(dlist, as.list=FALSE)
  #   2. iterate
  output = array(0,c(nlist,nlist))
  for (i in 1:(nlist-1)){
    tgti = list3d$array3d[,,i]
    for (j in (i+1):nlist){
      tgtj = list3d$array3d[,,j]
      output[i,j] <- landscape.norm.internal(tgti-tgtj, list3d$tseq, p)
      output[j,i] <- output[i,j]
    }
  }
 
  #############################################
  # Report
  if (as.dist){
    return(stats::as.dist(output))
  } else {
    return(output)
  }
}


# internal usage for older dependency -------------------------------------
#' @keywords internal
#' @noRd
landscape.normpair <- function(dlist, p=2, as.dist=FALSE){
  myp    = round(p)
  mydist = as.dist
  return(landscape.pdist(dlist, p=myp, as.dist=mydist))
}