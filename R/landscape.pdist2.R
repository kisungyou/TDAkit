#' Pairwise \eqn{L_p} Distances for Two Sets of Persistence Landscapes
#' 
#' @export
landscape.pdist2 <- function(dlist1, dlist2, p=2){
  #############################################
  # Preprocessing : checkers
  if (!dlist_check(dlist1)){    stop("* landscape.pdist2 : 'dlist1' is invalid. Consult with 'd2landscape' for each element in the list.")  }
  if (!dlist_check(dlist2)){    stop("* landscape.pdist2 : 'dlist2' is invalid. Consult with 'd2landscape' for each element in the list.")  }
  if (p < 1){
    stop("* landscape.pdist2 : 'p'-norm is defined for the range [1,Inf]. ")
  }
  m = length(dlist1)
  n = length(dlist2)
  
  #############################################
  # Main Computation
  #   1. merge two lists
  dlists = c(dlist1, dlist2)
  #   2. rearrange as 3d array
  list3d = dlist_adjust(dlists, as.list=FALSE)
  #   3. iterate (be careful since we have merged two lists)
  output = array(0,c(m,n))
  for (i in 1:m){
    tgti = list3d$array3d[,,i]
    for (j in 1:n){
      tgtj = list3d$array3d[,,(m+j)]
      output[i,j] = landscape.norm.internal(tgti-tgtj, list3d$tseq, p)
    }
  }
  
  #############################################
  # Report
  return(output)
}