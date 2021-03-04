#' Multi-sample Energy Test of Equal Distributions
#' 
#' Also known as \eqn{k}-sample problem, it tests whether multiple functional summaries 
#' are equally distributed or not via Energy statistics. 
#'
#' @param fslist a length-\eqn{N} list of functional summaries of persistent diagrams.
#' @param label a length-\eqn{N} vector of class labels.
#' @param method (case-sensitive) name of methods; one of \code{"original"} or \code{"disco"}.
#' @param mc.iter number of bootstrap replicates.
#' 
#' @return a (list) object of S3 class \code{htest} containing:\describe{
#'   \item{method}{name of the test.}
#'   \item{statistic}{a test statistic.}
#'   \item{p.value}{\eqn{p}-value under \eqn{H_0} of equal distributions.}
#' }
#' 
#' @examples 
#' \donttest{
#' # ---------------------------------------------------------------------------
#' #         Test for Equality of Distributions via Energy Statistics
#' #
#' # We will compare dim=0's top-5 landscape functions with 
#' # - Class 1 : 'iris' dataset with noise
#' # - Class 2 : samples from 'gen2holes()'
#' # - Class 3 : samples from 'gen2circles()'
#' # ---------------------------------------------------------------------------
#' ## Generate Data and Diagram from VR Filtration
#' ndata     = 10
#' list_rips = list()
#' for (i in 1:ndata){
#'   dat1 = as.matrix(iris[,1:4]) + matrix(rnorm(150*4), ncol=4)
#'   dat2 = gen2holes(n=100, sd=1)$data
#'   dat3 = gen2circles(n=100, sd=1)$data
#'   
#'   list_rips[[i]] = diagRips(dat1, maxdim=1)
#'   list_rips[[i+ndata]] = diagRips(dat2, maxdim=1)
#'   list_rips[[i+(2*ndata)]] = diagRips(dat3, maxdim=1)
#' }
#' 
#' ## Compute Persistence Landscapes from Each Diagram with k=5 Functions
#' list_land0 = list()
#' for (i in 1:(3*ndata)){
#'   list_land0[[i]] = diag2landscape(list_rips[[i]], dimension=0, k=5)
#' }
#' 
#' ## Create Label and Run the Test with Different Options
#' list_lab = c(rep(1,ndata), rep(2,ndata), rep(3,ndata))
#' fseqdist(list_land0, list_lab, method="original")
#' fseqdist(list_land0, list_lab, method="disco")
#' }
#' 
#' @concept summaries
#' @export
fseqdist <- function(fslist, label, method=c("original","disco"), mc.iter=999){
  ## PREPROCESSING
  dtype = check_list_summaries("fseqdist", fslist)
  label = as.factor(label)
  if (length(label)!=length(fslist)){
    stop("* fseqdist : length of 'label' must be equal to length of 'fslist'.")
  }
  mymethod = match.arg(method)
  myiter   = max(99, round(mc.iter))
  
  ## L2 DISTANCE
  dmat = TDAkit::fsdist(fslist, p=2, as.dist=FALSE)
  dobj = stats::as.dist(dmat)
  
  ## LABEL
  ulabels = sort(unique(label), decreasing = FALSE)
  nunique = length(ulabels)
  orderlab = base::order(label)
  mat.dist = dmat[orderlab,orderlab]
  vec.size = rep(0,nunique)
  for (i in 1:nunique){
    vec.size[i] = length(which(label==ulabels[i]))
  }
  
  # COMPUTATION : MAIN CASE
  if (all(mymethod=="original")){
    output = energy::eqdist.etest(dmat, vec.size, distance=TRUE, method="original", R=myiter)
    output$data.name = deparse(substitute(fslist))
  } else {
    output = energy::eqdist.etest(dmat, vec.size, distance=TRUE, method="discoB", R=myiter)
    output$method    = "DISCO : between-sample component approach"
    output$data.name = deparse(substitute(fslist))
  } 
  
  # Report
  return(output)
}