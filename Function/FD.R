FD <- function(MeanCpGInfoGroup){
  Res <- abs((max(MeanCpGInfoGroup$CanMean) - max(MeanCpGInfoGroup$NorMean)) / (min(MeanCpGInfoGroup$CanMean) - min(MeanCpGInfoGroup$NorMean)))
  return(Res)
}



