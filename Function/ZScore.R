#====================== Z Score discriminant between two vector ===================
ZScore <- function(xy){
  ## Convert type and remove NA from input vectors
  x <- xy[NameSamCan]
  y <- xy[NameSamNor]
  x <- x[which(!is.na(x))]
  y <- y[which(!is.na(y))]
  ## Claculate mean x and y
  meanx <- xy["CanMean"]
  meany <- xy["NorMean"]
  ## Calculate stndars deviation for c(x,y)
  library(fishmethods)
  sdxy <- sqrt(combinevar(c(meanx, meany) , c(xy["CanVar"],xy["NorVar"]) , c(length(x),length(y)))[2])
  ## Calculate ZScore
  S = abs(meanx-meany) / (sdxy + .Machine$double.eps)
  return(S)
}