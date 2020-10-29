#' Generate Two Intersecting Circles
#' 
#' It generates data from two intersecting circles.
#' 
#' @param n the total number of observations to be generated.
#' @param sd level of additive white noise.
#' 
#' @return a list containing\describe{
#' \item{data}{an \eqn{(n\times 2)} data matrix for row-stacked observations.}
#' \item{label}{a length-\eqn{n} vector for class label.}
#' }
#' 
#' @examples 
#' ## Generate Data with Different Noise Levels
#' nn = 200
#' x1 = gen2circles(n=nn, sd=0)
#' x2 = gen2circles(n=nn, sd=0.1)
#' x3 = gen2circles(n=nn, sd=0.25)
#' 
#' ## Visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' plot(x1$data, pch=19, main="sd=0.00", col=x1$label)
#' plot(x2$data, pch=19, main="sd=0.10", col=x2$label)
#' plot(x3$data, pch=19, main="sd=0.25", col=x3$label)
#' par(opar)
#' 
#' @concept data
#' @export
gen2circles <- function(n=496, sd=0){
  ## PREPROCESSING
  myn  = max(5, round(n))
  mysd = max(0, as.double(sd))
  
  ## GENERATE
  vec.r = stats::runif(n, min=0, max=1)
  vec.t = stats::runif(n, min=0, max=(2*pi))
  
  tmp = cbind(vec.r*cos(vec.t), vec.r*sin(vec.t))
  if (sd > 0){
    tmp = tmp + matrix(rnorm(n*2, sd=mysd),ncol=2)
  }
  n1 = round(n/2)
  tmp[1:n1,1] = tmp[1:n1,1] + 1.25
  
  
  ## RETURN
  output = list()
  output$data  = tmp
  output$label = c(rep(1,n1), rep(2,myn-n1))
  return(output)
}