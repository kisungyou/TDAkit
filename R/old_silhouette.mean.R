#' Mean of Multiple Persistence Silhouettes
#' 
#' @export
silhouette.mean <- function(slist){
  #############################################
  # Preprocessing : checkers
  myweight = rep(1,length(slist))/length(slist)
  if (!check_list_silhouette(slist)){
    stop("* silhouette.mean : an input 'slist' should be a list of silhouettes as 'kit.silhouette' objects. Consult with 'd2silhouette' function.")
  }

  #############################################
  # Main Computation using 'silhouette.wsum'
  return(silhouette.wsum(slist, myweight))
}
