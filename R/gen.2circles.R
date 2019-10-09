#' Generate Two Intertwined Circles
#' 
#' 
#' @examples 
#' ## see two sets of intertwined circles
#' nn = 200
#' x1 = gen.2circles(n=nn)
#' x2 = gen.2circles(n=nn, sd=0.05)
#' x3 = gen.2circles(n=nn, sd=0.1)
#' cols = as.factor(c(rep(1,nn/2),rep(2,nn/2)))
#' 
#' ## visualize
#' opar = par(mfrow=c(1,3), pty="s")
#' plot(x1[,1],x1[,2],pch=18,main="sd=0.00",col=cols)
#' plot(x2[,1],x2[,2],pch=18,main="sd=0.05",col=cols)
#' plot(x3[,1],x3[,2],pch=18,main="sd=0.10",col=cols)
#' on.exit(par(opar))
#' 
#' @export
gen.2circles <- function(n=100, sd=0){
  #############################################
  # Preprocessing : checkers
  n = round(n)
  if (n < 2){
    stop("* gen.2circles : select a larger 'n' please.")
  }
  if ((length(sd)>1)||(sd < 0)){
    stop("* gen.2circles : 'sd' should be a nonnegative real number.")
  }
  
  #############################################
  # generate
  n1 = round(n/2)
  n2 = round(n-n1)
  x1 = TDA::circleUnif(n1)
  x2 = TDA::circleUnif(n2)
  x2[,1] = x2[,1] + 1.25
  xx = rbind(x1,x2)
  if (sd > 0){
    xx = xx + matrix(rnorm(n*2, sd=sd), ncol=2)
  }
  
  #############################################
  # report
  return(xx)
}
