#================== Calculate ranks and normalization these =======================
Ranking <- function(DSSResi){
  CpGInfoDataMod <- data.frame(reshape::cast(CpGInfoData, SP~NameSample, value = if(MCMRFlag) "MC" else "MRatio"))
  CpGInfoDataMod <- CpGInfoDataMod[order(CpGInfoDataMod$SP),]
  DSSResi <- DSSResi[order(DSSResi$pos),]
  # cat(sum(CpGInfoDataMod$SP != DSSResi$pos)) # if print zero then data is matched correct
  ## Calling RankingMethods function
  CpGInfoDataMod <- RankingMethods(CpGInfoDataMod, DSSResi)
  
  ## Normalization
  CpGInfoDataMod[,CriNames[which(as.logical(CriFlags))]] <- NormalizeRanks(CpGInfoDataMod[,CriNames[which(as.logical(CriFlags))]])
  
  ## Calculate high or low methylation state (1 mean's high and 0 mean's low)
  CpGInfoDataMod <- cbind(CpGInfoDataMod, apply(CpGInfoDataMod[,CriNames[min(which(as.logical(CriFlags)))],drop=FALSE],1, function(x) x >= HighHypThr))
  colnames(CpGInfoDataMod)[dim(CpGInfoDataMod)[2]] <- "HighHyp"
  
  ## Determine DMC according to user-defined threshold
  # CpGInfoDataMod <- cbind(CpGInfoDataMod, apply(CpGInfoDataMod[,"RRA",drop=FALSE],1, function(x) x <= SigRank))
  CpGInfoDataMod <- cbind(CpGInfoDataMod, apply(CpGInfoDataMod[,CriNames[which(as.logical(CriFlags))],drop=FALSE],1, function(x) sum(x>=SigRank)>=NOfRankMethodSig))
  colnames(CpGInfoDataMod)[dim(CpGInfoDataMod)[2]] <- "DMCFlag"
  
  CpGInfoDataMod <- cbind(CpGInfoDataMod, "pvalue" = apply(CpGInfoDataMod[,c(SampleName),drop=FALSE], 1, function(x) t.test.mod(x[NameSamCan],x[NameSamNor])))
  
  return(CpGInfoDataMod)
}




