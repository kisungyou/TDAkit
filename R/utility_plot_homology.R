#' Plot Persistent Homology via Barcode or Diagram
#' 
#' Given a persistent homology of the data represented by a reconstructed 
#' complex in S3 class \code{homology} object, visualize it as either a barcode 
#' or a persistence diagram using \pkg{ggplot2}.
#' 
#' @param x a \code{homology} object. 
#' @param ... extra parameters including\describe{
#' \item{method}{type of visualization; either \code{"barcode"} or \code{"diagram"}.}
#' }
#' 
#' @return a \pkg{ggplot2} object.
#' 
#' @examples 
#' \donttest{
#' # Use 'iris' data
#' XX = as.matrix(iris[,1:4])
#' 
#' # Compute VR Diagram 
#' homology = diagRips(XX)
#' 
#' # Plot with 'barcode'
#' opar <- par(no.readonly=TRUE)
#' plot(homology, method="barcode")
#' par(opar)
#' }
#' 
#' @concept utility
#' @export
plot.homology <- function(x, ...){
  ## preprocess
  object = x
  if (!inherits(object, "homology")){
    stop("* plot.homology : input 'x' should be a 'homology' object.")
  }

  params = list(...)
  pnames = names(params)
  
  mymethod  = ifelse("method"%in%pnames, params$method, "barcode")
  vismethod = match.arg(mymethod,c("barcode","diagram"))
  
  ## case branching
  outplot <- switch(vismethod,
                    "barcode" = vis_homology_barcode(object),
                    "diagram" = vis_homology_diagram(object))
  return(outplot)
}


# individual functions ----------------------------------------------------
#' @keywords internal
#' @noRd
vis_homology_barcode <- function(object){
  # extract information
  dat_vertical <- dat_dims <- dat_birth <- dat_death <- NULL
  birth <- death <- vertical <- dimension <- NULL
  dat_dims     = as.factor(object$Dimension)
  dat_birth    = object$Birth
  dat_death    = object$Death
  
  # prepare data for plotting
  dat_vertical = (1:length(dat_birth))/length(dat_birth)
  df.plot = data.frame(birth=dat_birth, death=dat_death, 
                       vertical=dat_vertical, dimension=dat_dims)
  
  # ggplot object
  ggout <- ggplot2::ggplot(data=df.plot) + 
    ggplot2::geom_segment(ggplot2::aes_string(x="birth",y="vertical",
                                              xend="death",yend="vertical",
                                              colour="dimension")) + 
    ggplot2::theme_minimal() + 
    ggplot2::xlab("Radius") + 
    ggplot2::ylab("") + 
    ggplot2::theme(axis.text.y  = ggplot2::element_blank(),               # remove y-axis text
                   axis.ticks.y = ggplot2::element_blank(),               # remove y-tick
                   panel.grid   = ggplot2::element_blank(),               # remove grid line
                   axis.line.x  = ggplot2::element_line(colour="black")) + # add x-axis line) 
    ggplot2::coord_fixed()
  return(ggout)
}
#' @keywords internal
#' @noRd
vis_homology_diagram <- function(object){
  # extract information
  dat_dims <- dat_birth <- dat_death <- NULL
  birth <- death <- dimension <- NULL
  dat_dims  = as.factor(object$Dimension)
  dat_birth = object$Birth
  dat_death = object$Death
  df.plot   = data.frame(birth=dat_birth, death=dat_death, dimension=dat_dims)
  
  range.max = max(max(dat_birth), max(dat_death))
  
  # ggplot object
  ggout <- ggplot2::ggplot(df.plot) + 
    ggplot2::geom_point(ggplot2::aes_string(x="birth", y="death", 
                                            colour="dimension",
                                            shape="dimension")) + 
    ggplot2::theme_minimal() + 
    ggplot2::xlab("Birth") + 
    ggplot2::ylab("Death") + 
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.line  = ggplot2::element_line(colour="black")) + 
    ggplot2::geom_abline(intercept=0, slope=1) + 
    ggplot2::xlim(0, range.max) + 
    ggplot2::ylim(0, range.max) + 
    ggplot2::coord_fixed()
  
  return(ggout)
}

## issue with no visible binding
#  help 1 : https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
#  help 2 : https://github.com/STAT545-UBC/Discussion/issues/451
#
# my solutions
#