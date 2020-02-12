#' Weighted Sum of Multiple Accumulative Persistence Functions
#' 
#' 
#' @export
apf.wsum <- function(alist, weight){
  #############################################
  # Preprocessing : checkers
  if (missing(weight)){
    warning("* apf.wsum : since 'weight' is missing, we apply arithmetic mean of landscapes.")
    weight = rep(1,length(alist))/length(alist)
  }
  if (!check_list_apf(alist)){
    stop("* apf.wsum : an input 'alist' should be a list of APFs as 'kit.apf' objects. Consult with 'd2apf' function.")
  }
  if (length(weight)!=length(alist)){
    stop("* apf.wsum : length of 'weight' vector should equal to the number of APFs in 'slist'.")
  }
  
  #############################################
  # Main Computation
  mainout  = adjust_list_apf(alist, as.list=FALSE) # columns are APFs
  dat.tseq = mainout$tseq
  dat.func = mainout$array 
  vec.output = rep(0,length(dat.tseq))
  for (i in 1:length(weight)){
    vec.output = vec.output + (weight[i])*as.vector(dat.func[,i])
  }
  
  #############################################
  # Report the results
  output = list()
  output$lambda = vec.output
  output$tseq   = dat.tseq
  output$dimension = alist[[1]]$dimension
  class(output) = "kit.apf"
  return(output)
}