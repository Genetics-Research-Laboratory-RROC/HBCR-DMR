#======= Calculate ranking measures for each CpG site with apply function ========
RankingMethods <- function(CpGInfoDataMod, DSSResi){
  ## Calculate cancer mean and variance methylation info
  CpGInfoDataMod <- cbind(CpGInfoDataMod, "CanMean" = DSSResi$mu1)
  CpGInfoDataMod <- cbind(CpGInfoDataMod, "CanVar" = DSSResi$var1)
  ## Calculate normal mean and variance methylation info
  CpGInfoDataMod <- cbind(CpGInfoDataMod, "NorMean" = DSSResi$mu2)
  CpGInfoDataMod <- cbind(CpGInfoDataMod, "NorVar" = DSSResi$var2)
  ## Calculate hyper or hypo state
  CpGInfoDataMod <- cbind(CpGInfoDataMod, HyperFlag = apply(CpGInfoDataMod[,c("CanMean","NorMean"),drop=FALSE], 1, function(x) (x["CanMean"]>x["NorMean"])))
  ## Calculate Between versus Within discriminant measure
  if(CriFlags[1] == 1){
    CpGInfoDataMod <- cbind(CpGInfoDataMod, apply(CpGInfoDataMod, 1, function(x) BetVSWit(x)))
    colnames(CpGInfoDataMod)[dim(CpGInfoDataMod)[2]] <- CriNames[1]
  }
  ## Calculate Information Gain discriminant measure
  if (CriFlags[2] == 1) {
    CpGInfoDataMod <- cbind(CpGInfoDataMod, apply(CpGInfoDataMod[,c(SampleName),drop=FALSE], 1, function(x) InfoGain(x[NameSamCan],x[NameSamNor])))
    colnames(CpGInfoDataMod)[dim(CpGInfoDataMod)[2]] <- CriNames[2]
  }
  ## Calculate fisher ratio discriminant measure
  if (CriFlags[3] == 1){ 
    CpGInfoDataMod <- cbind(CpGInfoDataMod, apply(CpGInfoDataMod, 1, function(x) FisherRatio(x)))
    colnames(CpGInfoDataMod)[dim(CpGInfoDataMod)[2]] <- CriNames[3]
  }
  ## Calculate ZScore discriminant measure
  if (CriFlags[4] == 1){
    CpGInfoDataMod <- cbind(CpGInfoDataMod, apply(CpGInfoDataMod, 1, function(x) ZScore(x)))
    colnames(CpGInfoDataMod)[dim(CpGInfoDataMod)[2]] <- CriNames[4]
  }
  
  ## Calculate welch test discriminant measure
  if (CriFlags[5] == 1){
    CpGInfoDataMod <- cbind(CpGInfoDataMod, apply(CpGInfoDataMod, 1, function(x) WelchTTest(x)))
    colnames(CpGInfoDataMod)[dim(CpGInfoDataMod)[2]] <- CriNames[5]
  }
  return(CpGInfoDataMod)
}



