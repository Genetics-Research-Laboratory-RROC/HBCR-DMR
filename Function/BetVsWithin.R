#========== Between versus Within discriminant measure between two vector ==========
BetVSWit <- function(xy){
  ## Convert type and remove NA from input vectors
  x <- xy[NameSamCan]
  y <- xy[NameSamNor]
  x <- x[which(!is.na(x))]
  y <- y[which(!is.na(y))]
  ## Claculate mean x, y and concatenate these
  meanx <- xy["CanMean"]
  meany <- xy["NorMean"]
  library(fishmethods)
  meanxy <- combinevar(c(meanx, meany) , c(xy["CanVar"],xy["NorVar"]) , c(length(x),length(y)))[1]
  ## Calculate sw and sb
  SW <- sum((x-meanx)^2 ,(y-meany)^2)
  SB <- sum((c(meanx,meany) - meanxy)^2)
  ## Calculate S and return this
  S = SB / (SW + .Machine$double.eps)
  
  return(S)
}