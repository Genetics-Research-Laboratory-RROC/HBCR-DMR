#========================= Extract DMR Infoemation ================================
ExtractDMRInfo <- function(SubCluster, DMRRanksCpGSite){
  ## Calculate region length
  DMRInfo <- cbind(SubCluster, RegLen = apply(SubCluster, 1, function(x) x[2]-x[1]+1))
  ## Calculate number of CpG in each region
  DMRInfo <- cbind(DMRInfo, NOfCpG = apply(SubCluster, 1, function(x) length(which(x[1] <= DMRRanksCpGSite$SP & x[2] >= DMRRanksCpGSite$SP))))
  ## Calculate fold difference
  DMRInfo <- cbind(DMRInfo, FD = apply(SubCluster, 1, function(x) 
    FD(DMRRanksCpGSite[which(x[1] <= DMRRanksCpGSite$SP & x[2] >= DMRRanksCpGSite$SP), c("CanMean","NorMean")])))
  ## Calculate hyper or hypo state (1 mean's hyper and 0 mean's hypo)
  DMRInfo <- cbind(DMRInfo, HyperFlag = apply(SubCluster, 1, function(x) 
    MaxVote(DMRRanksCpGSite[which(x[1] <= DMRRanksCpGSite$SP & x[2] >= DMRRanksCpGSite$SP), "HyperFlag"])))
  ## Calculate high or low methylation state (1 mean's high and 0 mean's low) with MaxVote function
  DMRInfo <- cbind(DMRInfo, HighHypFlag = apply(SubCluster, 1, function(x) 
    MaxVote(DMRRanksCpGSite[which(x[1] <= DMRRanksCpGSite$SP & x[2] >= DMRRanksCpGSite$SP), "HighHyp"])))
  ## Calculate variance
  DMRInfo <- cbind(DMRInfo, VarRank = apply(SubCluster, 1, function(x) 
    var(DMRRanksCpGSite[ which(x[1] <= DMRRanksCpGSite$SP & x[2] >= DMRRanksCpGSite$SP) , CriNames[min(which(as.logical(CriFlags)))] ])))
  DMRInfo <- as.data.frame(DMRInfo)
  # Change type column's data from factor to numeric
  DMRInfo[,c("RStart","REnd","RegLen","NOfCpG","FD","VarRank")] <- apply(DMRInfo[,c("RStart","REnd","RegLen","NOfCpG","FD","VarRank")], 2, function(x) as.numeric(as.character(x)))
  DMRInfo$VarRank[is.na(DMRInfo$VarRank)] <- 0
  
  return(DMRInfo)
}


