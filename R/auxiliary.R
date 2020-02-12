# Auxiliary Functions -----------------------------------------------------
# (1) check_list_landscape  : see if it is a list of landscapes
#     check_list_silhouette 
#     check_list_apf
#     check_list_napf
# (2) dlist_adjust : adjust and return to 
#         - as.list = FALSE; list$array3d 3d array (T,K,nlist) & list$tseq
#         - as.list = TRUE;  dlist by logical
# (3) adjust_list_silhouette 
#         - as.list = TRUE;  return as list with common tseq
#         - as.list = FALSE; list$array 2d array (T,nlist) & list$tseq
# (4) adjust_list_apf
#         - as.list = TRUE;  return as list with common tseq
#         - as.list = FALSE; list$array 2d array (T,nlist) & list$tseq
#         - since it's ECDF-like, we need ECDF manipulation
# (5) adjust_list_napf
#         - as.list = TRUE;  return as list with common tseq
#         - as.list = FALSE; list$array 2d array (T,nlist) & list$tseq
#         - since it's ECDF-like, we need ECDF manipulation
# (6) compute_median

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
#' @keywords internal
#' @noRd
check_list_apf <- function(alist){
  cond1 = (is.list(alist)&&(length(alist)>1))
  cond2 = all(unlist(lapply(alist, inherits, "kit.apf"))==TRUE)
  if (cond1&&cond2){
    return(TRUE)
  } else {
    return(FALSE)
  }
}
#' @keywords internal
#' @noRd
check_list_napf <- function(alist){
  cond1 = (is.list(alist)&&(length(alist)>1))
  cond2 = all(unlist(lapply(alist, inherits, "kit.napf"))==TRUE)
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

# (4) adjust_list_apf -----------------------------------------------------
#' @keywords internal
#' @noRd
adjust_list_apf_single <- function(xval, yval, xvalnew){
  ## my manual padding of maximal values
  rxval    = as.double(utils::tail(xval,    n=1))
  rxvalnew = as.double(utils::tail(xvalnew, n=1))

  if (rxvalnew > rxval){
    xinc  = max((100*.Machine$double.eps), (rxvalnew-rxval)/10)
    xvaladd = base::seq(from=rxval+xinc, to=rxvalnew, length.out=10)
    yvaladd = rep(max(yval)[1], 10)
    
    xval = c(xval, xvaladd)
    yval = c(yval, yvaladd)
  }

    # approximation
  ytmp = stats::approx(xval, yval, xout=xvalnew, method = "constant")$y
  return(ytmp)
}
#' @keywords internal
#' @noRd
adjust_list_apf <- function(slist, as.list=TRUE){
  # 1. check if it is a list of APFs
  if (!check_list_apf(slist)){
    stop("* adjust_list_apf : the given input is not a proper list of APFs.")
  }
  # 2. check same dimension or not
  nlist   = length(slist)
  vec.dim = rep(0,nlist)
  for (i in 1:nlist){
    vec.dim[i] = slist[[i]]$dimension
  }
  if (length(unique(vec.dim))!=1){
    stop("* adjust_list_apf : all APFs should be computed for the same dimension.")
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
    
    ytmp = adjust_list_apf_single(tgti$tseq, tgti$lambda, mytseq)
    ytmp[is.na(ytmp)] = 0.0
    ytmp[(ytmp<0.0)]  = 0.0
    if (as.list){
      tmplist = list(lambda=ytmp, tseq=mytseq, dimension=tgti$dimension)
      class(tmplist) = "kit.apf"
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


# (5) adjust_list_napf ----------------------------------------------------
#' @keywords internal
#' @noRd
adjust_list_napf <- function(slist, as.list=TRUE){
  # 1. check if it is a list of APFs
  if (!check_list_napf(slist)){
    stop("* adjust_list_napf : the given input is not a proper list of NAPFs.")
  }
  # 2. check same dimension or not
  nlist   = length(slist)
  vec.dim = rep(0,nlist)
  for (i in 1:nlist){
    vec.dim[i] = slist[[i]]$dimension
  }
  if (length(unique(vec.dim))!=1){
    stop("* adjust_list_napf : all NAPFs should be computed for the same dimension.")
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
    
    ytmp = adjust_list_apf_single(tgti$tseq, tgti$lambda, mytseq)
    ytmp[is.na(ytmp)] = 0.0
    ytmp[(ytmp<0.0)]  = 0.0
    ytmp[(ytmp>1.0)]  = 1.0
    if (as.list){
      tmplist = list(lambda=ytmp, tseq=mytseq, dimension=tgti$dimension)
      class(tmplist) = "kit.napf"
      listout[[i]]   = tmplist
    } else {
      arrayout[,i] = ytmp
    }
  }
  arrayout[is.na(arrayout)]=0.0
  arrayout[(arrayout < 0)] =0.0
  arrayout[(arrayout > 1)] =1.0
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

# (6) compute_median ------------------------------------------------------
# (6-1) Weiszfeld
#       For smaller ones, do the padding with sqrt(.Machine$double.eps)
#' @keywords internal
#' @noRd
compute_median_weiszfeld <- function(cols, tseq, weights, maxiter, abstol, fname, print.progress=FALSE){
  # parameter
  nsummary = ncol(cols)
  
  # prepare
  sol.old = base::rowMeans(cols%*%diag(weights))
  epsmall = 100*(.Machine$double.eps)
  
  # iterate
  wnt = rep(0,nsummary)
  for (i in 1:round(maxiter)){
    # compute |Z(t) - x_n|
    for (j in 1:nsummary){
      vecj   = (sol.old - as.vector(cols[,j]))^2
      wnt[j] = sqrt(simple_integral_1d(vecj, tseq))
    }
    # zero padding
    if (any(wnt < epsmall)){
      wnt = wnt + epsmall
    }
    # now, invert
    wnt = weights/wnt
    # try to update
    sol.new = base::rowSums(cols%*%diag(wnt))/base::sum(wnt)
    
    # updating
    sol.inc = sqrt(sum(abs(sol.new-sol.old)^2))
    sol.old = sol.new
    if (sol.inc < abstol){
      print(paste0("* ",fname," : iteration ",i,"/",maxiter," terminates for small incremental change..", sep=""))
      break
    }
    if (i%%5 == 0){
      print(paste0("* ",fname," : iteration ",i,"/",maxiter," complete..", sep=""))
    }
  }
  
  # return
  return(sol.old)
}