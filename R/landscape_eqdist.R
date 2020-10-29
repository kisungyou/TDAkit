#' Multisample Energy Test of Equal Distributions
#' 
#' Also known as \eqn{k}-sample problem, it tests whether multiple landscapes 
#' are equally distributed or not via Energy statistics. 
#'
#' @param dlist a length-\eqn{N} list of \code{'landscape'} objects from \code{\link{diag2landscape}}.
#' @param label a length-\eqn{N} vector of class labels.
#' @param method (case-insensitive) name of methods; one of \code{"original"} or \code{"discoB"}.
#' @param mc.iter number of bootstrap replicates.
#' 
#' @return a (list) object of S3 class \code{htest} containing:\describe{
#'   \item{method}{name of the test.}
#'   \item{statistic}{a test statistic.}
#'   \item{p.value}{\eqn{p}-value under \eqn{H_0}.}
#' }
#' 
#' @examples 
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
#' pleqdists(list_land0, list_lab, method="original")
#' pleqdists(list_land0, list_lab, method="DiscoB")
#' 
#' @concept landscape
#' @export
pleqdists <- function(dlist, label, method=c("original","discoB"), mc.iter=999){
  ## PREPROCESSING
  #  dlist input
  if (!check_list_landscape(dlist)){
    stop("* pleqdist : 'dlist' should be a list of 'landscape' objects.")
  }
  #  label
  label = as.factor(label)
  if (length(label)!=length(dlist)){
    stop("* pleqdist : length of 'label' must be equal to length of 'dlist'.")
  }
  #  others
  mymethod = match.arg(tolower(method), tolower(c("original","discoB","discoF")))
  myiter   = max(99, round(mc.iter))
  
  ## COMPUTATION : PRELIMINARY
  #  Pairwise L2 Distance
  dmat = pldist(dlist, p=2, as.dist = FALSE)
  dobj = stats::as.dist(dmat)
  #  let's take care of label
  ulabels = sort(unique(label), decreasing = FALSE)
  nunique = length(ulabels)
  #  re-arrange distance matrix and label vector as sizes 
  orderlab = base::order(label)
  mat.dist = dmat[orderlab,orderlab]
  vec.size = rep(0,nunique)
  for (i in 1:nunique){
    vec.size[i] = length(which(label==ulabels[i]))
  }
  
  # COMPUTATION : MAIN CASE
  if (all(mymethod=="original")){
    output = energy::eqdist.etest(dmat, vec.size, distance=TRUE, method="original", R=myiter)
    output$data.name = deparse(substitute(dlist))
  } else if (all(mymethod=="discob")){
    output = energy::eqdist.etest(dmat, vec.size, distance=TRUE, method="discoB", R=myiter)
    output$method    = "DISCO : between-sample component approach"
    output$data.name = deparse(substitute(dlist))
  } 
  
  #############################################
  # Report
  return(output)
}
