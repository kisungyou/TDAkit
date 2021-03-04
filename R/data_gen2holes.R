#' Generate Two Intertwined Holes
#' 
#' It generates data from two intertwine circles with empty interiors(holes).
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
#' x1 = gen2holes(n=nn, sd=0)
#' x2 = gen2holes(n=nn, sd=0.1)
#' x3 = gen2holes(n=nn, sd=0.25)
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
gen2holes <- function(n=496, sd=0){
  ## PREPROCESSING
  myn  = max(5, round(n))
  mysd = max(0, as.double(sd))
  
  ## GENERATE
  n1 = round(myn/2)
  n2 = round(myn-n1)
  
  theta1 = stats::runif(n1, min=0, max=2*pi)
  theta2 = stats::runif(n2, min=0, max=2*pi)
  
  x1 = cbind(base::cos(theta1), base::sin(theta1))
  x2 = cbind(base::cos(theta2), base::sin(theta2))
  x2[,1] = x2[,1] + 1.25
  xx = rbind(x1,x2)
  if (sd > 0){
    xx = xx + matrix(rnorm(n*2, sd=sd), ncol=2)
  }
  
  ## RETURN
  output = list()
  output$data  = xx
  output$label = c(rep(1,n1), rep(2,n2))
  return(output)
}
