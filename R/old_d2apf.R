#' Convert Persistence Diagram to Accumulated Persistence Function
#' 
#' @export
d2apf <- function(diagram, dimension=1, tseq){
  #############################################
  # Preprocessing : Check
  # input diagram
  if (!inherits(diagram,"diagram")){
    cond1 = is.matrix(diagram)
    cond2 = (ncol(diagram)==3)
    cond3 = all(as.vector(diagram[,3]) >= as.vector(diagram[,2]))
    if (!(cond1&&cond2&&cond3)){
      stop("* d2apf : for a matrix-valued input 'diagram', please match the format as 'diagram' object from TDA package.")
    }
  }
  dat.dimension = round(as.vector(diagram[,1]))
  dat.birth     = as.vector(diagram[,2])
  dat.death     = as.vector(diagram[,3])
  
  # dimension
  dimension = round(dimension)
  if (!(dimension %in% dat.dimension)){
    stop("* d2apf : 'diagram' does not have information pertaining to the given 'dimension'.")
  }
  
  idin = base::intersect(which((!is.infinite(dat.death))), which(dat.dimension==dimension))
  dat.dimension = dat.dimension[idin] # separate out the ones
  dat.birth     = dat.birth[idin]
  dat.death     = dat.death[idin]
  
  # tseq
  if (missing(tseq)){
    mytseq = seq(from=0, to=(1.1*max(dat.birth+dat.death)/2), length.out=200)
  } else {
    mytseq = sort(tseq, decreasing = FALSE)
  }
  
  #############################################
  # Main Computation
  # unlike Silhouette, there is no weight now
  out.vec = rep(0,length(mytseq))
  tmp.den = 0
  for (i in 1:length(dat.birth)){
    # select birth and death pair
    bj = dat.birth[i]
    dj = dat.death[i]
    
    # update
    out.vec = out.vec + (abs(bj-dj))*(1.0*((bj+dj)<=2*mytseq))
    tmp.den = tmp.den + abs(bj-dj)
  }
  ## if normalize, out.vec/tmp.den
  
  #############################################
  # Return the output
  res = list(lambda=out.vec, tseq=mytseq, dimension=dimension)
  class(res) = "kit.apf"
  return(res)
}


# # personal test : visualize -----------------------------------------------
# library(TDA)
# library(TDAkit)
# 
# # data generation
# X = TDAkit::gen.2circles(n=100)
# diagx = ripsDiag(X, maxdimension = 1, maxscale = Inf)$diagram
# 
# # construct apf function
# mytseq = seq(from=0, to=0.8, length.out=1234)
# apf0 = d2apf(diagx, dimension=0, tseq = mytseq)
# apf1 = d2apf(diagx, dimension=1, tseq = mytseq)
# 
# # visualize
# par(mfrow=c(1,2))
# plot(apf0$tseq, apf0$lambda, "l", main="dimension 0")
# plot(apf1$tseq, apf1$lambda, "l", main="dimension 1")
