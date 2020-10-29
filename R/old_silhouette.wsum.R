#' Weighted Sum of Multiple Persistence Silhouettes
#' 
#' 
#' @export
silhouette.wsum <- function(slist, weight){
  #############################################
  # Preprocessing : checkers
  if (missing(weight)){
    warning("* silhouette.wsum : since 'weight' is missing, we apply arithmetic mean of landscapes.")
    weight = rep(1,length(slist))/length(slist)
  }
  if (!check_list_silhouette(slist)){
    stop("* silhouette.wsum : an input 'slist' should be a list of silhouettes as 'kit.silhouette' objects. Consult with 'd2silhouette' function.")
  }
  if (length(weight)!=length(slist)){
    stop("* silhouette.wsum : length of 'weight' vector should equal to the number of silhouettes in 'slist'.")
  }
  
  #############################################
  # Main Computation
  mainout = adjust_list_silhouette(slist, as.list=FALSE)
  dat.tseq = mainout$tseq
  dat.func = mainout$array # each column is silhouette function
  vec.output = rep(0,length(dat.tseq))
  for (i in 1:length(weight)){
    vec.output = vec.output + (weight[i])*as.vector(dat.func[,i])
  }
  
  #############################################
  # Report the results
  output = list()
  output$lambda = vec.output
  output$tseq   = dat.tseq
  output$dimension = slist[[1]]$dimension
  class(output) = "kit.silhouette"
  return(output)
}
