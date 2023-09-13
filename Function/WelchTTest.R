#===================== Welch test discriminant between two vector ================
WelchTTest <- function(xy){
  ## Convert type and remove NA from input vectors
  x <- xy[NameSamCan]
  y <- xy[NameSamNor]
  x <- x[which(!is.na(x))]
  y <- y[which(!is.na(y))]
  ## Claculate mean x and y
  meanx <- xy["CanMean"]
  meany <- xy["NorMean"]
  ## Claculate length x and y
  NOfx <- length(x)
  NOfy <- length(y)
  ## Claculate standard deviation x and y
  sdx <- sqrt(xy["CanVar"])
  sdy <- sqrt(xy["NorVar"])
  ## Calculate welch test
  S = abs(meanx-meany) / (sqrt((sdx/NOfx)+(sdy/NOfy))+.Machine$double.eps)
  return(S)
}


