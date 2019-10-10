#' Generate Two Intersecting Filled Circles
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
gen.2circles <- function(n=496, sd=0){
  #############################################
  # Preprocessing : checkers
  n = round(n)
  if (n < 2){
    stop("* gen.2circles : select a larger 'n' please.")
  }
  if ((length(sd)>1)||(sd < 0)){
    stop("* gen.2circles : 'sd' should be a nonnegative real number.")
  }
  mysd = as.double(sd)
  
  #############################################
  # generate
  vec.r = stats::runif(n, min=0, max=1)
  vec.t = stats::runif(n, min=0, max=(2*pi))
  
  tmp = cbind(vec.r*cos(vec.t), vec.r*sin(vec.t))
  if (sd > 0){
    tmp = tmp + matrix(rnorm(n*2, sd=mysd),ncol=2)
  }
  n1 = round(n/2)
  tmp[1:n1,1] = tmp[1:n1,1] + 1.25
  
  #############################################
  # report
  return(tmp)
}
