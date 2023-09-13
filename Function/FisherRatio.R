#================ Fisher ratio discriminant measure between two vector ============
FisherRatio <- function(xy){
  ## Convert type and remove NA from input vectors
  x <- xy[NameSamCan]
  y <- xy[NameSamNor]
  x <- x[which(!is.na(x))]
  y <- y[which(!is.na(y))]
  ## Claculate mean x and y
  meanx <- xy["CanMean"]
  meany <- xy["NorMean"]
  ## Claculate standard deviation x and y
  sdx = xy["CanVar"]
  sdy = xy["NorVar"]
  ## Calculate fisher ratio
  FR = ((meanx-meany)^2) / (sdx + sdy + .Machine$double.eps)
  return(FR)
}