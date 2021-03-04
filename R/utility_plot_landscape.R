#' Plot Persistence Landscape
#' 
#' Given a persistence landscape object in S3 class \code{landscape}, visualize the 
#' landscapes using \pkg{ggplot2}.
#' 
#' @param x a \code{landscape} object. 
#' @param ... extra parameters including \describe{
#' \item{top.k}{the number of landscapes to be plotted (default: 5).}
#' \item{colored}{a logical; \code{TRUE} to assign different colors for landscapes, or \code{FALSE} to use grey color for all landscapes.}
#' }
#' 
#' @return a \pkg{ggplot2} object.
#' 
#' @examples 
#' \donttest{
#' # Use 'iris' data
#' XX = as.matrix(iris[,1:4])
#' 
#' # Compute Persistence diagram and landscape of order 0 
#' homology  = diagRips(XX)
#' landscape = diag2landscape(homology, dimension=0)
#' 
#' # Plot with 'barcode'
#' opar <- par(no.readonly=TRUE)
#' plot(landscape)
#' par(opar)
#' }
#' 
#' @concept utility
#' @export
plot.landscape <- function(x, ...){
  ## preprocess
  object = x
  if (!inherits(object, "landscape")){
    stop("* landscape : input 'x' should be a 'landscape' object.")
  }
  params    = list(...)
  pnames    = names(params)
  
  top.k     = ifelse("top.k"%in%pnames, round(params$top.k), 5)
  colored   = ifelse("colored"%in%pnames, as.logical(params$colored), FALSE)
  numtoshow = min(max(1, round(top.k)), ncol(object$lambda))
  
  ## prepare for inputs
  df_tseq <- df_lbds <- df_show <- NULL
  tseq <- lambda <- group <- NULL
  df_tseq = rep(as.vector(object$tseq), times=numtoshow)
  df_lbds = c()
  for (i in 1:numtoshow){
    df_lbds = c(df_lbds, as.vector(object$lambda[,i]))
  }
  df_nums = rep(1:numtoshow, each=length(as.vector(object$tseq)))
  df_show = data.frame(tseq=df_tseq, lambda=df_lbds, group=as.factor(df_nums))

  ## visualize with ggplot2
  if (colored){
    ggout <- ggplot2::ggplot(data=df_show) + 
      ggplot2::geom_line(ggplot2::aes_string(x="tseq", y="lambda", group="group", colour="group", linetype="group")) + 
      ggplot2::scale_color_discrete() + 
      ggplot2::theme_minimal() + 
      ggplot2::theme(legend.position = "none")
  } else {
    ggout <- ggplot2::ggplot(data=df_show) + 
      ggplot2::geom_line(ggplot2::aes_string(x="tseq", y="lambda", group="group", linetype="group"),  color="grey") + 
      ggplot2::theme_minimal() + 
      ggplot2::theme(legend.position = "none")
  }
  
  ## post-processing of the figure
  ggout <- ggout + 
    ggplot2::xlab("") + 
    ggplot2::ylab("") + 
    ggplot2::theme(panel.grid         = ggplot2::element_blank(),
                   axis.line.x.bottom = ggplot2::element_line(colour="black"),
                   axis.line.y.left   = ggplot2::element_line(colour="black")) 
  return(ggout)
}