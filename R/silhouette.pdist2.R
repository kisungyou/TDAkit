#' Pairwise \eqn{L_p} Distances for Two Sets of Persistence Silhouettes
#' 
#' @export
silhouette.pdist2 <- function(slist1, slist2, p=2){
  #############################################
  # Preprocessing : checkers
  if (!check_list_silhouette(slist1)){    stop("* silhouette.pdist2 : 'slist1' is invalid. Consult with 'd2silhouette' for each element in the list.")  }
  if (!check_list_silhouette(slist2)){    stop("* silhouette.pdist2 : 'slist2' is invalid. Consult with 'd2silhouette' for each element in the list.")  }
  if (p < 1){
    stop("* silhouette.pdist2 : 'p'-norm is defined for the range [1,Inf]. ")
  }
  p = as.double(p)
  m = length(slist1)
  n = length(slist2)
  
  #############################################
  # Main Computation
  #   1. merge two lists
  slists = c(slist1, slist2)
  #   2. extract some information
  list2d = adjust_list_silhouette(slists, as.list=FALSE)
  funcs1 = list2d$array[,1:m]
  funcs2 = list2d$array[,(m+1):(m+n)]
  mytseq = list2d$tseq
  #   3. iterate (be careful since we have merged two lists)
  output = array(0,c(m,n))
  for (i in 1:m){
    tgti = funcs1[,i]
    for (j in 1:n){
      tgtj = funcs2[,j]
      if (is.infinite(p)){
        output[i,j] <- max(abs(tgti-tgtj))
      } else {
        trimmed = matrix(abs(trim.the.norm(tgti-tgtj))^p, ncol=1)
        output[i,j] <- (simple_integral(trimmed, mytseq)^(1/p))
      }
    }
  }
  
  #############################################
  # Report
  return(output)
}