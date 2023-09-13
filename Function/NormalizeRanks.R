#== Normalization each ranking method result in [0,1] range for each chromosome ===
NormalizeRanks <- function(RankData){
  ## Ranking methods loop
  for (i in 1:dim(RankData)[2]) {
    temp <- RankData[,i]
    ## Normalization
    temp <- (temp - min(temp)) / (max(temp) - min(temp))
    RankData[,i] <- temp
  }
  return(RankData)
}