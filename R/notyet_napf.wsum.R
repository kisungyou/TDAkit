#' Weighted Sum of Multiple Normalized Accumulative Persistence Functions
#' 
#' 
#' @export
napf.wsum <- function(flist, weight){
  #############################################
  # Preprocessing : checkers
  if (missing(weight)){
    warning("* napf.wsum : since 'weight' is missing, we apply arithmetic mean of NAPFs.")
    weight = rep(1,length(flist))/length(flist)
  }
  if (!check_list_napf(flist)){
    stop("* napf.wsum : an input 'flist' should be a list of silhouettes as 'kit.napf' objects. Consult with 'd2napf' function.")
  }
  if (length(weight)!=length(flist)){
    stop("* napf.wsum : length of 'weight' vector should equal to the number of NAPFs in 'flist'.")
  }
  
  #############################################
  # Main Computation
  mainout = adjust_list_napf(flist, as.list=FALSE)
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
  output$dimension = flist[[1]]$dimension
  class(output) = "kit.napf"
  return(output)
}