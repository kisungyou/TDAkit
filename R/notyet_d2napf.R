#' Convert Persistence Diagram to Normalized Accumulated Persistence Function
#' 
#' @export
d2napf <- function(diagram, dimension=1, tseq){
  #############################################
  # Preprocessing : Check
  # input diagram
  if (!inherits(diagram,"diagram")){
    cond1 = is.matrix(diagram)
    cond2 = (ncol(diagram)==3)
    cond3 = all(as.vector(diagram[,3]) >= as.vector(diagram[,2]))
    if (!(cond1&&cond2&&cond3)){
      stop("* d2napf : for a matrix-valued input 'diagram', please match the format as 'diagram' object from TDA package.")
    }
  }
  dat.dimension = round(as.vector(diagram[,1]))
  dat.birth     = as.vector(diagram[,2])
  dat.death     = as.vector(diagram[,3])
  
  # dimension
  dimension = round(dimension)
  if (!(dimension %in% dat.dimension)){
    stop("* d2napf : 'diagram' does not have information pertaining to the given 'dimension'.")
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
  napfvec = out.vec/tmp.den
  
  #############################################
  # Return the output
  res = list(lambda=napfvec, tseq=mytseq, dimension=dimension)
  class(res) = "kit.napf"
  return(res)
}
