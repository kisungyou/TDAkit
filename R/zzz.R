.pkgenv <- new.env(parent = emptyenv())

.onLoad <- function(lib, pkg){
  Rdpack::Rdpack_bibstyles(package=pkg, authors="LongNames")
  invisible(NULL)
}

.onAttach <- function(...){
  ## Retrieve Year Information
  date <- date()
  x <- regexpr("[0-9]{4}", date)
  this.year <- substr(date, x[1], x[1] + attr(x, "match.length") - 1)
  
  # Retrieve Current Version
  this.version = utils::packageVersion("TDAkit")
  
  ## Print on Screen
  packageStartupMessage("** --------------------------------------------------------- **")
  packageStartupMessage("**          Toolkit for Topological Data Analysis")
  packageStartupMessage("**")
  packageStartupMessage("** Version    : ",this.version,"      (",this.year,")",sep="")
  packageStartupMessage("** Maintainer : Kisung You (kisung.you@outlook.com)")
  packageStartupMessage("**")
  packageStartupMessage("** Please share any bugs or suggestions to the maintainer.")
  packageStartupMessage("** --------------------------------------------------------- **")
}

.onUnload <- function(libpath) {
  library.dynam.unload("TDAkit", libpath)
}
