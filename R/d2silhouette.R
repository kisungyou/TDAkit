#' Convert Persistence Diagram to Persistence Silhouette
#' 
#' @export
d2silhouette <- function(diagram, dimension=1, p=2, tseq){
  #############################################
  # Preprocessing : Check
  # input diagram
  if (!inherits(diagram,"diagram")){
    cond1 = is.matrix(diagram)
    cond2 = (ncol(diagram)==3)
    cond3 = all(as.vector(diagram[,3]) >= as.vector(diagram[,2]))
    if (!(cond1&&cond2&&cond3)){
      stop("* d2silhouette : for a matrix-valued input 'diagram', please match the format as 'diagram' object from TDA package.")
    }
  }
  dat.dimension = round(as.vector(diagram[,1]))
  dat.birth     = as.vector(diagram[,2])
  dat.death     = as.vector(diagram[,3])
  
  # weight function
  myp = as.double(p)
  if ((length(myp)>1)||(myp <= 0)||(is.infinite(myp))){
    stop("* d2silhouette : weight exponent parameter 'p' should be a positive real number.")
  }
  
  # dimension
  dimension = round(dimension)
  if (!(dimension %in% dat.dimension)){
    stop("* d2silhouette : 'diagram' does not have information pertaining to the given 'dimension'.")
  }
  
  idin = base::intersect(which((!is.infinite(dat.death))), which(dat.dimension==dimension))
  dat.dimension = dat.dimension[idin] # separate out the ones
  dat.birth     = dat.birth[idin]
  dat.death     = dat.death[idin]
  
  # tseq
  if (missing(tseq)){
    mytseq = seq(from=0, to=max(dat.death), length.out=200)
  } else {
    mytseq = sort(tseq, decreasing = FALSE)
  }
  
  #############################################
  # Main Computation
  # We use weight function |a-b|^p
  out.numerator   = rep(0,length(mytseq))
  out.denominator = 0
  for (i in 1:length(dat.birth)){
    # select birth and death pair
    bj = dat.birth[i]
    dj = dat.death[i]
    
    # update denominator
    weightj = (base::abs(bj-dj)^myp)
    out.denominator = out.denominator + weightj
    # update numerator
    out.numerator = out.numerator + (weightj*base::pmax(base::pmin(mytseq-bj, dj-mytseq), rep(0,length(mytseq))))
  }
  # finalize
  lbdfun = (out.numerator/out.denominator)
  
  #############################################
  # Return the output
  res = list(lambda=lbdfun, tseq=mytseq, dimension=dimension)
  class(res) = "kit.silhouette"
  return(res)
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
# x3 = TDAkit::d2silhouette(diagx, dimension=1, tseq=mytseq)
# 
# par(mfrow=c(1,3))
# plot(mytseq, x1[,1],"l")
# lines(mytseq,x1[,2],"l",col="red")
# plot(mytseq, x2$lambda[,1], "l")
# lines(mytseq, x2$lambda[,2], col="red")
# plot(mytseq, x3$lambda, "l")
# 
# x11() # test different values
# par(mfrow=c(2,3))
# for (i in 1:6){
#   xx = TDAkit::d2silhouette(diagx, dimension=1, p=(i/2))
#   plot(xx$tseq, xx$lambda, "l")
# }