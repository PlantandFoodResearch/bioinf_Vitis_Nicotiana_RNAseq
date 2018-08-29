panel.points <- function(x, y, col=NULL, fitLoess = FALSE,  maTransform = FALSE, ...){
  if(is.null(col)) {
    dcols <- densCols(x,y)
    if(maTransform) {
      points((x+y)/2, y-x, col=dcols, ...)
    } else {
      points(x, y, col=dcols, ...)
    }
  } else {
    if(maTransform) {
      points((x+y)/2, y-x, col=col, ...)
    } else {
      points(x, y, col=col, ...)
    }
  }
  if(fitLoess) {
    loessObj <- loess(y~x)
    yhat <- fitted(loessObj)
    o <- order(yhat)
    lines(x[o], yhat[o], col="red")
  }
  if(maTransform) {
    abline(h=0, cex=0.01, col="darkgrey")
  } else {
    abline(0,1, cex=0.01, col="darkgrey")
  }
}


panel.corr <- function(x,y, maTransform = FALSE, ...){
  NAs  <- (is.na(x) | is.na(y))
  infs <- (is.infinite(x) | is.infinite(y))
  correlation <- round(cor(x[!(NAs | infs)], y[!(NAs | infs)]), 2)
  if(maTransform) {
    text(x=10, y=0, labels=substitute(rho==correlation, list(correlation=correlation)))
  } else {
    text(x=7.5, y=7.5, labels=substitute(rho==correlation, list(correlation=correlation)))
  }
}


openDevice <- function (device = "png", name = "Rplot", ...){
  dev <- get(device)
  Interactive <- c("X11", "quartz", "windows")
  if (any(device == Interactive)) {
    dev(...)
  }
  else {
    extMap <- list(cairo_pdf = "pdf", cairo_ps = "ps", jpeg = "jpg", 
                   postscript = "ps", tiff = "tif")
    if (!is.null(extMap[[device]])) {
      filename <- paste(name, extMap[[device]], sep = ".")
    }
    else {
      filename <- paste(name, device, sep = ".")
    }
    dev(file = filename, ...)
    return(filename)
  }
  return(invisible())
}


gitRevision <- function(allcommits=FALSE){
  ## Fetch all commits for a specific file, knitr::current_input() can be NULL
  ## If the filename is omitted, the latest commits will be returned
  shasum <- system(paste("git log --pretty=format:'%h'", knitr::current_input()), intern=TRUE)
  if(!allcommits) {
    shasum <- shasum[1]
  }
  return(shasum)
}
