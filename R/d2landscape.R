#' Convert Persistence Diagram to Persistence Landscape
#' 
#' 
#' @export
d2landscape <- function(diagram, dimension=1, k=0, tseq, 
                        inf.replace=FALSE, inf.repval=2*max(diagram[!is.infinite(diagram[,3]),3])){
  #############################################
  # Preprocessing : Check
  # input diagram
  if (!inherits(diagram,"diagram")){
    cond1 = is.matrix(diagram)
    cond2 = (ncol(diagram)==3)
    cond3 = all(as.vector(diagram[,3]) >= as.vector(diagram[,2]))
    if (!(cond1&&cond2&&cond3)){
      stop("* d2landscape : for a matrix-valued input 'diagram', please match the format as 'diagram' object from TDA package.")
    }
  }
  dat.dimension = round(as.vector(diagram[,1]))
  dat.birth     = as.vector(diagram[,2])
  dat.death     = as.vector(diagram[,3])
  # dimension
  dimension = round(dimension)
  if (!(dimension %in% dat.dimension)){
    stop("* d2landscape : 'diagram' does not have information pertaining to the given 'dimension'.")
  }
  if (inf.replace){
    diagram = diagram[is.infinite(diagram[,3]),3] = inf.repval
  }
  idin = base::intersect(which((!is.infinite(dat.dimension))), which(dat.dimension==dimension))
  dat.dimension = dat.dimension[idin] # separate out the ones
  dat.birth     = dat.birth[idin]
  dat.death     = dat.death[idin]
  # k : number of landscape functions
  kk = round(k)
  if (kk < 1){
    myKK = length(dat.birth)
  } else {
    myKK = kk
  }
  if (length(dat.birth) < myKK){
    stop("* d2landscape : for the matching dimension, we have smaller number of topological features than 'kk'. Try a smaller 'k' value.")
  }
  # tseq
  if (missing(tseq)){
    mytseq = seq(from=0, to=max(dat.death), length.out=200)
  } else {
    mytseq = sort(tseq, decreasing = FALSE)
  }
  
  #############################################
  # Main Computation
  # wrap the data into my diagram : columns are each lambda function; k=1 -> column matrix
  mydiag = cbind(dat.dimension, dat.birth, dat.death)
  colnames(mydiag) = c("dimension","Birth","Death")
  # # computation
  # if (myKK == 1){
  #   lambdas = (TDA::landscape(mydiag, dimension=unique(dat.dimension), 
  #                           tseq=mytseq, KK=1))
  # } else {
  #   lambdas = TDA::landscape(mydiag, dimension=unique(dat.dimension),
  #                            tseq=mytseq, KK=1:myKK)
  # }
  lambdas = compute_lambdas(mytseq, dat.birth, dat.death, myKK)
  
  
    
  
  
  #############################################
  # Return the output
  res = list(lambda=lambdas, tseq=mytseq, dimension=dimension)
  class(res) = "kit.landscape"
  return(res)
}

# compute lambdas ---------------------------------------------------------
#' @keywords internal
#' @noRd
compute_lambdas <- function(tseq, births, deaths, maxK){
  ntest = length(births)
  ntime = length(tseq)
  
  output = array(0,c(ntime,ntest))
  for (i in 1:ntest){
    b = births[i]
    d = deaths[i]
    output[,i] = base::pmax(base::pmin(tseq-b, d-tseq), rep(0,ntime))
  }
  
  newout = array(0,c(ntime,ntest))
  for (i in 1:ntime){
    newout[i,] = sort(output[i,], decreasing = TRUE)
  }
  return((newout[,1:maxK]))
}


# # Personal Test
# library(TDA)
# x1 = TDA::circleUnif(30)
# x2 = TDA::circleUnif(30)*0.5
# x2[,1] = x2[,1] + rep(1,nrow(x2))
# X <- rbind(x1,x2)
# diagx = ripsDiag(X, maxdimension=1, maxscale=Inf)$diagram
# mytseq  = seq(from=0, to=11, length.out=100)
# x1 = TDA::landscape(diagx, dimension=1, KK=1:2, tseq=mytseq)
# x2 = TDAkit::d2landscape(diagx, dimension=1, k=2, tseq = mytseq)
# 
# par(mfrow=c(1,2))
# plot(mytseq, x1[,1],"l")
# lines(mytseq,x1[,2],"l",col="red")
# plot(mytseq, x2$lambda[,1], "l")
# lines(mytseq, x2$lambda[,2], col="red")