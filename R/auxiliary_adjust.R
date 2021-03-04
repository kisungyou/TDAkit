## AUXILIARY : ADJUSTING THE INPUTS
#  (01) adjust_list_landscapes
#         - as.list = FALSE; list$array3d 3d array (T,K,nlist) & list$tseq
#         - as.list = TRUE;  dlist by logical
#  (02) adjust_list_silhouette
#         - as.list = FALSE; list$array 2d array (T,nlist) & list$tseq
#         - as.list = TRUE;  dlist by logical



# (02) adjust_list_silhouette ---------------------------------------------
#' @keywords internal
#' @noRd
adjust_list_silhouette <- function(slist, as.list=FALSE){
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
    stop("* adjust_list_silhouette : all silhouettes should be computed for the same dimension.")
  }
  # 3. extract common tseq
  mytmin = 0
  mytmax = 0
  mytlength = 0
  for (i in 1:nlist){
    cseq = slist[[i]]$tseq
    if (min(cseq) <= mytmin){      mytmin = min(cseq)    }
    if (max(cseq) >= mytmax){      mytmax = max(cseq)    }
    if (length(cseq) > mytlength){      mytlength = length(cseq)    }
  }
  mytseq = base::seq(from=mytmin, to=mytmax, length=mytlength)
  
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


# (01) adjust_list_landscapes ---------------------------------------------
#' @keywords internal
#' @noRd
adjust_list_landscapes <- function(dlist, as.list=FALSE){
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
  mytmin = 0
  mytmax = 0
  mytlength = 0
  for (i in 1:nlist){
    cseq = dlist[[i]]$tseq
    if (min(cseq) <= mytmin){      mytmin = min(cseq)    }
    if (max(cseq) >= mytmax){      mytmax = max(cseq)    }
    if (length(cseq) > mytlength){      mytlength = length(cseq)    }
  }
  vect = base::seq(from=mytmin, to=mytmax, length=mytlength)
  
  
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
      class(tmplist) = "landscape"
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
