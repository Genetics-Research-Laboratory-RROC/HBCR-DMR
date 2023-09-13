#==================== Find subcluster (region) in all clusters ====================
FindingSubCluster <- function(Data, Clusters){
  SubCluster <- NULL
  ## Cluster loop
  for (j in 1:dim(Clusters)[1]) 
  {
    STemp <- Clusters[j,1]
    ETemp <- Clusters[j,2]
    Idx <- which(STemp <= Data$SP & ETemp >= Data$SP)
    # Calculate subcluster in each cluster with SegmentationBound function
    SubCluster <- rbind(SubCluster,SegmentationBound(Data[Idx,]))
  }
  colnames(SubCluster) <- c('RStart','REnd')
  return(SubCluster)
}



