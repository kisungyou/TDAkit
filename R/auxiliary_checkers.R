## AUXILIARY : CHECKERS
#  (01) check_diagram        : (T/F) for Persistence Diagram
#  (02) check_list_landscape : check a list of landscapes of same dimension



# (01) check_diagram ------------------------------------------------------
#' @keywords internal
#' @noRd
check_diagram <- function(diagram){
  cond1 = is.data.frame(diagram)
  cond2 = (ncol(diagram)==3)
  cond3 = all(as.vector(diagram[,3]) >= as.vector(diagram[,2]))
  cond4 = all(diagram >= 0)
  
  if (cond1&&cond2&&cond3&&cond4){
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# (02) check_list_landscape -----------------------------------------------
check_list_landscape <- function(dlist){
  cond1 = (is.list(dlist)&&(length(dlist)>1))
  cond2 = all(unlist(lapply(dlist, inherits, "landscape"))==TRUE)
  
  vd = rep(0, length(dlist))
  for (i in 1:length(dlist)){
    vd[i] = round(dlist[[i]]$dimension)
  }
  cond3 = (length(unique(vd))==1)
  if (cond1&&cond2&&cond3){
    return(TRUE)
  } else {
    return(FALSE)
  }
}