#================ Information gain discriminant between two vector ================
InfoGain <- function(x, y){
  ## Convert type and remove NA from input vectors
  x <- as.numeric(x) 
  y <- as.numeric(y)
  x <- x[which(!is.na(x))]
  y <- y[which(!is.na(y))]
  ## CORElearn report error if all entire of x and y is same
  if (length(unique(c(x,y)))==1) return(0)
  Data <- data.frame(MInfo=c(x,y),Labels=factor(c(rep(-1,length(x)),rep(1,length(y)))))
  ## Calculate information gain with CORElearn package
  # Data <- CORElearn::attrEval(Labels~., Data, estimator = "InfGain")
  # Data <- FSelector::information.gain(Labels~., Data)
  Data <- RWeka::InfoGainAttributeEval(Labels~., data = Data)

  return(Data)
}