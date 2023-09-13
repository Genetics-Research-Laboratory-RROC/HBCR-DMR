RegPValCal <- function(DMRInfo, CpGInfo){
  #==================================================================
  library(metap)
  colnames(DMRInfo) <- c("chr","start","end")
  colnames(CpGInfo) <- c("chr","pos","pval")
  
  DMRInfo$chr <- as.character(DMRInfo$chr)
  CpGInfo$chr <- as.character(CpGInfo$chr)
  
  RPValCal <- rep(0, dim(DMRInfo)[1])
  for (i in 1:dim(DMRInfo)[1]) {
    Start <- DMRInfo$start[i]
    End   <- DMRInfo$end[i]
    Chr   <- DMRInfo$chr[i]
    Ind <- which(CpGInfo$chr == Chr & CpGInfo$pos >= Start & CpGInfo$pos <= End)
    if (length(Ind) < 5){
      RPValCal[i] <- median(CpGInfo$pval[Ind], na.rm = TRUE)
    } else {
      RPValCal[i] <- max(unlist(metap::allmetap(CpGInfo$pval[Ind], method = "all")$p)[c("minimump","sumlog","sump","sumz")],na.rm = TRUE)
    }
  }
  return(RPValCal)
}