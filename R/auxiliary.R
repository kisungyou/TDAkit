# Auxiliary Functions -----------------------------------------------------
# (1) check_list_landscape  : see if it is a list of landscapes
#     check_list_silhouette 
#     check_list_apf
#     check_list_napf

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
# (6) compute_median_weiszfeld
#     compute_median_sgd
#     compute_median_asgd

# (1) check_list_landscape ------------------------------------------------

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
      if (print.progress){
        print(paste0("* ",fname," : iteration ",i,"/",maxiter," terminates for small incremental change..", sep=""))  
      }
      break
    }
    if (print.progress){
      print(paste0("* ",fname," : iteration ",i,"/",maxiter," complete..", sep=""))
    }
  }
  
  # return
  return(sol.old)
}
# (6-2) SGD
#       Stochastic Gradient by Cardot (Bernoulli, 2013)
#' @keywords internal
#' @noRd
compute_median_sgd <- function(cols, tseq, weights, maxiter, abstol, fname, print.progress=FALSE,
                               par.C, par.alpha){
  # parameter & shuffling
  nsummary = ncol(cols)
  cols     = cols[,base::sample(1:nsummary)]
  epsmall  = 100*(.Machine$double.eps)
  
  # setup for iterative update
  z0     = as.vector(cols[,1])
  zold   = z0
  zstack = c() # columns would be
  for (i in 1:(nsummary-1)){
    # prepare
    xnow   = cols[,(i+1)] # we should start from the second one.
    xzdiff = xnow-zold    # x_{n+1} - z_{n}
    xznorm = sqrt(simple_integral_1d( (xzdiff^2), tseq))
    if (xznorm < epsmall){
      xznorm = xznorm + epsmall
    }
    gamma.sgd = ifelse(i < 2, par.C, par.C*(i^(-par.alpha)))
    
    # update
    znew   = zold + gamma.sgd*(xzdiff/xznorm)
    zstack = cbind(zstack, znew)
    zinc   = sqrt(sum(znew-zold)^2)
    zold   = znew
    
    # stopping criterion
    if (zinc < abstol){
      if (print.progress){
        print(paste0("* ",fname," : iteration ",i,"/",nsummary-1," terminates for small incremental change..", sep=""))  
      }
      break
    }
    if (print.progress){
      print(paste0("* ",fname," : iteration ",i,"/",nsummary-1," complete..", sep=""))
    }
  }
  
  # return
  return(zold)
}
# (6-3) ASGD
#       Aveeraged Stochastic Gradient by Cardot (Bernoulli, 2013)
#' @keywords internal
#' @noRd
compute_median_asgd <- function(cols, tseq, weights, maxiter, abstol, fname, print.progress=FALSE,
                               par.C, par.alpha){
  # parameter & shuffling
  nsummary = ncol(cols)
  cols     = cols[,base::sample(1:nsummary)]
  epsmall  = 100*(.Machine$double.eps)
  
  # setup for iterative update
  z0     = as.vector(cols[,1])
  zold   = z0
  zstack = c() # columns would be
  for (i in 1:(nsummary-1)){
    # prepare
    xnow   = cols[,(i+1)] # we should start from the second one.
    xzdiff = xnow-zold    # x_{n+1} - z_{n}
    xznorm = sqrt(simple_integral_1d( (xzdiff^2), tseq))
    if (xznorm < epsmall){
      xznorm = xznorm + epsmall
    }
    gamma.sgd = ifelse(i < 2, par.C, par.C*(i^(-par.alpha)))
    
    # update
    znew   = zold + gamma.sgd*(xzdiff/xznorm)
    zstack = cbind(zstack, znew)
    zinc   = sqrt(sum(znew-zold)^2)
    zold   = znew
    
    # stopping criterion
    if (zinc < abstol){
      if (print.progress){
        print(paste0("* ",fname," : iteration ",i,"/",nsummary-1," terminates for small incremental change..", sep=""))  
      }
      break
    }
    if (print.progress){
      print(paste0("* ",fname," : iteration ",i,"/",nsummary-1," complete..", sep=""))
    }
  }
  
  # We need averaging
  z0bar = rep(0,length(tseq))
  for (i in 1:ncol(zstack)){
    z0bar = z0bar + (1/i)*(as.vector(zstack[,i])-z0bar)
  }
  return(z0bar)
}
