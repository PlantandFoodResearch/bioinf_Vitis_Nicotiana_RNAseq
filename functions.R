panel.points <- function(x, y, col=NULL, loess = FALSE,  maTransform = TRUE, ...){
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
  if(loess) {
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


panel.corr <- function(x,y, maTransform = TRUE, ...){
  NAs  <- (is.na(x) | is.na(y))
  infs <- (is.infinite(x) | is.infinite(y))
  correlation <- round(cor(x[!(NAs | infs)], y[!(NAs | infs)]), 2)
  if(maTransform) {
    text(x=10, y=0, labels=substitute(rho==correlation, list(correlation=correlation)))
  } else {
    text(x=7.5, y=7.5, labels=substitute(rho==correlation, list(correlation=correlation)))
  }
}
