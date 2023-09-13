#============= Cluster CpG site according to inner CpG distances ==================
ClusterFinding <- function(SP, DisThr){
  SP   <- unique(SP)
  ## Claculate distance between all pair CpG sites (with diff function)
  Dist <- diff(SP)
  ## Distance thresholding
  Dist <- which(Dist >= DisThr | Dist<0)
  ## Calculate cluster
  FirstCluster <- c(SP[1], SP[Dist[1]])                            # First region
  MidCluster  <- cbind(SP[(Dist[-length(Dist)]+1)],SP[Dist[-1]])   # All region except for first and last region
  LastIdx     <- Dist[length(Dist)]+1                             
  LastCluster <- c(SP[LastIdx], SP[length(SP)])                    # Last region
  Clusters <- as.data.frame(rbind(FirstCluster, MidCluster, LastCluster))
  colnames(Clusters) <- c("SPStart","SPEnd")
  row.names(Clusters) <- NULL
  return(Clusters)
}



