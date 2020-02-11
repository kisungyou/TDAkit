# Auxiliary Functions -----------------------------------------------------
# (1) check_list_landscape  : see if it is a list of landscapes
#     check_list_silhouette 
# (2) dlist_adjust : adjust and return to 
#         - as.list = FALSE; list$array3d 3d array (T,K,nlist) & list$tseq
#         - as.list = TRUE;  dlist by logical
# (3) adjust_list_silhouette 
#         - as.list = TRUE;  return as list with common tseq
#         - as.list = FALSE; list$array 2d array (T,nlist) & list$tseq

# (1) check_list_landscape ------------------------------------------------
#' @keywords internal
#' @noRd
check_list_landscape <- function(dlist){
  cond1 = (is.list(dlist)&&(length(dlist)>1))
  cond2 = all(unlist(lapply(dlist, inherits, "kit.landscape"))==TRUE)
  if (cond1&&cond2){
    return(TRUE)
  } else {
    return(FALSE)
  }
}
#' @keywords internal
#' @noRd
check_list_silhouette <- function(slist){
  cond1 = (is.list(slist)&&(length(slist)>1))
  cond2 = all(unlist(lapply(slist, inherits, "kit.silhouette"))==TRUE)
  if (cond1&&cond2){
    return(TRUE)
  } else {
    return(FALSE)
  }
}


# (2) dlist_adjust --------------------------------------------------------
#' @keywords internal
#' @noRd
dlist_adjust <- function(dlist, as.list=TRUE){
  # 1. check if it is a list of landscapes
  if (!check_list_landscape(dlist)){
    stop("* dlist_adjust : the given input is not a proper list of landscapes.")
  }
  # 2. check same dimension or not
  nlist   = length(dlist)
  vec.dim = rep(0,nlist)
  for (i in 1:nlist){
    vec.dim[i] = dlist[[i]]$dimension
  }
  if (length(unique(vec.dim))!=1){
    stop("* dlist_adjust : all landscapes should be computed for the same dimension.")
  }
  # 3. max.k
  vec.k = rep(0,nlist)
  for (i in 1:nlist){
    if (is.vector(dlist[[i]]$lambda)){
      vec.k[i] = 1
    } else {
      vec.k[i] = ncol(dlist[[i]]$lambda)
    }
  }
  maxk = unique(max(vec.k))
  # 4. vect
  vect = c()
  for (i in 1:nlist){
    vect = c(vect, dlist[[i]]$tseq)
  }
  vect = sort(unique(vect), decreasing = FALSE)
  
  # 5. compute into 3d array
  output = array(0,c(length(vect), maxk, nlist))
  for (i in 1:nlist){
    tgt = dlist[[i]]
    if ((length(tgt$tseq)==length(vect))&&(all(tgt$tseq==vect))){
      output[,1:ncol(tgt$lambda),i] = tgt$lambda
    } else {
      output[,1:ncol(tgt$lambda),i] = dlist_adjust_single(tgt$lambda, tgt$tseq, vect)
    }
  }
  output[is.na(output)] = 0
  output[(output<0)] = 0
  
  # 6. return the output
  if (as.list){
    outlist = list()
    for (i in 1:nlist){
      tmplist = list(lambda=output[,,i], tseq=vect, dimension=tgt$dimension)
      class(tmplist) = "kit.landscape"
      outlist[[i]] = tmplist
    }
    return(outlist)
  } else {
    outlist = list()
    outlist$array3d = output
    outlist$tseq    = vect
    return(outlist)
  }
}
#' @keywords internal
#' @noRd
dlist_adjust_single <- function(lmat, tseq, newtseq){
  if (is.vector(lmat)){
    lmat = matrix(lmat)
  }
  TT = length(newtseq)
  KK = ncol(lmat)
  output = array(0,c(TT,KK))
  for (i in 1:KK){
    ytmp = stats::approx(tseq, as.vector(lmat[,i]), xout=newtseq)$y
    ytmp[is.na(ytmp)] = 0
    output[,i] = ytmp
  }
  return(output)
}

# (3) adjust_list_silhouette ----------------------------------------------
#' @keywords internal
#' @noRd
adjust_list_silhouette <- function(slist, as.list=TRUE){
  # 1. check if it is a list of landscapes
  if (!check_list_silhouette(slist)){
    stop("* adjust_list_silhouette : the given input is not a proper list of silhouettes.")
  }
  # 2. check same dimension or not
  nlist   = length(slist)
  vec.dim = rep(0,nlist)
  for (i in 1:nlist){
    vec.dim[i] = slist[[i]]$dimension
  }
  if (length(unique(vec.dim))!=1){
    stop("* adjust_list_silhouette : all landscapes should be computed for the same dimension.")
  }
  # 3. extract common tseq
  mytseq = c()
  for (i in 1:nlist){
    cseq = slist[[i]]$tseq
    mytseq = sort(unique(c(mytseq, cseq)), decreasing = FALSE)
  }
  # 4. transform
  arrayout = array(0,c(length(mytseq), nlist))
  listout  = list()
  for (i in 1:nlist){
    tgti = slist[[i]]
    ytmp = stats::approx(tgti$tseq, tgti$lambda, xout=mytseq)$y
    ytmp[is.na(ytmp)] = 0.0
    ytmp[(ytmp<0.0)]  = 0.0
    if (as.list){
      tmplist = list(lambda=ytmp, tseq=mytseq, dimension=tgti$dimension)
      class(tmplist) = "kit.silhouette"
      listout[[i]]   = tmplist
    } else {
      arrayout[,i] = ytmp
    }
  }
  arrayout[is.na(arrayout)]=0.0
  arrayout[(arrayout < 0)] =0.0
  # 5. return
  if (as.list){
    return(listout)
  } else {
    output = list()
    output$array = arrayout
    output$tseq  = mytseq
    return(output)
  }
}