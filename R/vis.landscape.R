#' Visualize Persistence Landscape with ggplot2
#' 
#' if show.k = 0; remove all zeros,
#' 
#' @export
vis.landscape <- function(landscape, show.k=0){
  #############################################
  # Preprocessing : checkers
  if (!inherits(landscape, "kit.landscape")){
    stop("* vis.landscape : provide a compatiable landscape object from 'd2landscape'.")
  }
  show.k = round(show.k)

  #############################################
  # Create an object
  #   1. remove nulls
  if (show.k < 1){
    lmat = trim.the.norm(landscape$lambda)
  } else {
    lmat = landscape$lambda
    if (ncol(lmat) >= show.k){
      lmat = lmat[,1:show.k]
    } else {
      lmat = cbind(lmat, array(0,c(nrow(lmat), k-ncol(lmat))))
    }
  }
  if (is.vector(lmat)){
    lmat = as.matrix(lmat)
  }
  KK   = ncol(lmat)
  seqt = landscape$tseq
  
  #   2. create a weird data.frame
  df.x = c()
  df.y = c()
  df.lab = c()
  for (i in 1:KK){
    df.x = c(df.x, seqt)                # update t
    df.y = c(df.y, as.vector(lmat[,i])) # update y values
    df.lab = c(df.lab, rep(i, length(seqt)))
  }
  df.lab = as.factor(df.lab)
  mydf = data.frame(x=df.x, y=df.y, lab=df.lab)
  
  #   3. time to ggplot : theme_{bw, minimal, classic, void}
  ggp <- ggplot(mydf, aes(y=y, x=x, color=lab)) + geom_line() + 
    xlab("t") + ylab("lambda") + 
    scale_color_discrete(name="k") 

  
  #############################################
  # Return the ggplot2 object
  return(ggp)
}

# # personal test -----------------------------------------------------------
# library(TDA)
# x1 = TDA::circleUnif(30)
# x2 = TDA::circleUnif(30)*0.5
# x2[,1] = x2[,1] + rep(1,nrow(x2))
# X <- rbind(x1,x2)
# X <- X + matrix(rnorm(nrow(X)*ncol(X),sd=0.1), ncol=ncol(X))
# diagx  = ripsDiag(X, maxdimension = 1, maxscale = Inf)$diagram
# myland0 =  d2landscape(diagx, dimension=0, k=0, inf.replace = FALSE)
# myland1 =  d2landscape(diagx, dimension=1, k=0, inf.replace = FALSE)
# 
# vis00 = vis.landscape(myland0, show.k = 0)
# vis01 = vis.landscape(myland0, show.k = 1)
# vis02 = vis.landscape(myland0, show.k = 2)
# 
# 
# vis10 = vis.landscape(myland1, show.k = 0)
# vis11 = vis.landscape(myland1, show.k = 1)
# vis12 = vis.landscape(myland1, show.k = 2)
