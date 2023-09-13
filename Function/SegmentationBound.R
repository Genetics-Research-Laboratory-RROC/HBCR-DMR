#======================= Divide cluster to sub-cluster (region) ==================
SegmentationBound <- function(DataCluster){
  SP <- DataCluster$SP
  DMCFlag <- as.numeric(DataCluster$DMCFlag)
  ## Add 0 to first and last of DMCFlag
  DMCFlag <- c(0, DMCFlag, 0)
  ## different between DMCFlag and left shifted DMCFlag
  Difforward  <- DMCFlag[1:(length(DMCFlag)-1)] - DMCFlag[2:length(DMCFlag)]
  ## Remove last entire
  Difforward  <- Difforward[-length(Difforward)]
  ## different between left shifted DMCFlag and DMCFlag
  DifBackward <- DMCFlag[2:length(DMCFlag)] - DMCFlag[1:(length(DMCFlag)-1)]
  ## Remove first entire+
  DifBackward <- DifBackward[-1]
  ## When middle region
  Res <- NULL
  if (sum(Difforward==-1) >= 1 & sum(DifBackward==-1)>=1) {
    Res <- rbind(Res, cbind(SP[which(Difforward==-1)], SP[which(DifBackward==-1)]))
  }

  return(Res)
}


