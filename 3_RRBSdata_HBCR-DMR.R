# DMLtest_Function_Check.R  -> check My_DMLtest doing same orginal DMLtest
# Plot_Mean&Var.R           -> Visually difference between mean and variance normal and HBB
# firt run RRBSdata_Preparation.R
# then run Extract_Mean&Var_DSS.R
# now run below code
#======================================================================================
rm(list = ls())
gc()

#======================================================================================
#============================ Select requirment directory =============================
MainDir <- paste(tcltk::tk_choose.dir(default = getwd(), caption = "Select Main Directory"),"/",sep = "")

DirLoad <- paste(MainDir,"Data/CpGInfo/"   ,sep = "")
DirSave <- paste(MainDir,"Data/CpGDMRInfo/",sep = "")
DirFunc <- paste(MainDir,"Function/"       ,sep = "")

LFiles <- list.files(DirFunc)
LFiles <- LFiles[! LFiles %in% c("My_DMLtest.R","DMLtest.R","My_makeBSseqData.R")]
for (FuncName in LFiles) {
  source(paste(DirFunc,FuncName,sep = ""))
}

#=============== Extract sample information from defined data folder =============
library(RRBSdata)
data("rrbs")
data("islands")
library(GenomicRanges)

OBorOT <- unique(as.character(strand(rrbs@rowRanges)))
ChrNames <- unique(as.character(seqnames(rrbs@rowRanges)))
NOfAllSample <- dim(colData(rrbs))[1]
NameSamCan <- rownames(colData(rrbs))[which(colData(rrbs)$group == "cancer")]
NameSamNor <- rownames(colData(rrbs))[which(colData(rrbs)$group == "normal")]
SampleName <- c(NameSamCan, NameSamNor)
rm(rrbs) 

#============================ Define parameters ==================================
CriNames <- c("BetVsWithin","InfoGain","FisherRatio","ZScore","WelchTTest")
CriFlags <- c(1, 1, 1, 1, 1)
CovThr <- 1
DisThr <- 100
MCMRFlag <- 1   # 1=MC    0=MR
SigRank  <- 0.04
HighHypThr <- 0.5
RepetationThr <- 0.75
NOfRankMethodSig <- 3

#============================== Main Section =====================================
DMRs <- GenomicRanges::GRanges()
CpGInfo <- data.frame()
ElapsedTime <- 0
CreateFolder(DirSave)

for (StrandInd in OBorOT) {
  cat(paste(StrandInd, "processing start ...\n"))
  ## define and create result folder
  SubDir <- paste(DirSave,StrandInd,"/",sep = "")
  CreateFolder(SubDir)
  
  SubDir <- paste(DirSave,StrandInd,"/CpGSiteRank/",sep = "")
  CreateFolder(SubDir)
  SubDir <- paste(DirSave,StrandInd,"/DMRInfo/",sep = "")
  CreateFolder(SubDir)
  SubDir <- paste(DirSave,StrandInd,"/DMRInfo/Hyper/",sep = "")
  CreateFolder(SubDir)
  SubDir <- paste(DirSave,StrandInd,"/DMRInfo/Hypo/",sep = "")
  CreateFolder(SubDir)
  ## Chromosome loop
  DSSRes <- read.table(paste(MainDir,"Result/DSS_Output.txt",sep = ""),header = TRUE)
  
  for (ChrInd in ChrNames) {
    tictoc::tic()
    cat(paste(" ",ChrInd," ...\n",sep = ""))
    ## Read data from all samples
    CpGInfoData <- data.frame()
    TypeSample <- list.files(DirLoad)
    ## Sample loop
    for (TypeIndex in TypeSample){
      NameSample <- list.files(paste(DirLoad,TypeIndex,"/",sep = ""))
      for (NameIndex in NameSample){
        SubCpGInfoData <- paste(DirLoad,TypeIndex,"/",NameIndex,"/",StrandInd,"/",NameIndex,"_CPGInfo",StrandInd,"_",ChrInd,".txt",sep = "")
        SubCpGInfoData <- read.table(SubCpGInfoData,header = TRUE)
        CpGInfoData <- rbind(CpGInfoData, cbind(SubCpGInfoData,TypeSample = rep(TypeIndex,dim(SubCpGInfoData)[1]),NameSample = rep(NameIndex,dim(SubCpGInfoData)[1])))
      }
    }
    rm(SubCpGInfoData)
    ## Filtering according to coverage threshold
    CpGInfoData <- CpGInfoData[CpGInfoData$Cov >= CovThr,]
    ## Filtering according to frequency in samples
    SPRep <- as.data.frame(table(CpGInfoData$SP))
    SPRep[,1] <- as.numeric(levels(SPRep[,1]))[SPRep[,1]]
    SPRep <- SPRep[SPRep[,2] >= round(NOfAllSample*RepetationThr),1]
    CpGInfoData <- CpGInfoData[which(CpGInfoData$SP %in% SPRep),]
    rm(SPRep)
    ## Sort data according to CpG position
    CpGInfoData <- CpGInfoData[order(CpGInfoData$SP),]
    
    DSSResi <- DSSRes[which(DSSRes$chr == ChrInd),-1]
    CpGInfoData <- CpGInfoData[which(CpGInfoData$SP %in% DSSResi$pos),]
    DSSResi <- DSSResi[which(DSSResi$pos %in% unique(CpGInfoData$SP)),]
    
    ## Find cluster with according to inner CpG distance
    Clusters <- ClusterFinding(CpGInfoData$SP, DisThr)
    
    ## Ranking each CpG site with multiple discriminant measure methods
    DMRRanksCpGSite <- Ranking(DSSResi)
    CpGInfo <-rbind(CpGInfo, cbind("Chr" = rep(ChrInd,dim(DMRRanksCpGSite)[1]),DMRRanksCpGSite[,c("SP", "pvalue")]))
    DMRRanksCpGSite <- DMRRanksCpGSite[,which(! names(DMRRanksCpGSite) %in% "pvalue")]
    ## Extract subcluster (Region) in each cluster according to consecutive Significant CpG ()
    SubCluster <- FindingSubCluster(DMRRanksCpGSite[,c("SP","DMCFlag")], Clusters)
    
    ## Extract DMR Information
    SubClusterInfo <- ExtractDMRInfo(SubCluster, DMRRanksCpGSite)
    # SubClusterInfo <- SubClusterInfo[SubClusterInfo$NOfCpG > 4,]
    ##
    DMRs <- c(DMRs, GenomicRanges::GRanges(seqnames = ChrInd, 
                                           ranges = IRanges::IRanges(start = SubClusterInfo$RStart, end = SubClusterInfo$REnd)))
    ETi     <- tictoc::toc()
    ElapsedTime <- ElapsedTime + (ETi$toc - ETi$tic)
    
    ## Save result as text file
    write.table(DMRRanksCpGSite, paste(DirSave,StrandInd,"/CpGSiteRank/CpGSiteRank_",toString(ChrInd),".txt",sep = ""),
                col.names = TRUE, row.names = FALSE)
    write.table(SubClusterInfo[SubClusterInfo$HyperFlag == TRUE ,!(names(SubClusterInfo) %in% "HyperFlag")],
                paste(DirSave,StrandInd,"/DMRInfo/Hyper/DMRInfo_",ChrInd,".txt",sep = ""),
                col.names = TRUE, row.names = FALSE)
    write.table(SubClusterInfo[SubClusterInfo$HyperFlag == FALSE,!(names(SubClusterInfo) %in% "HyperFlag")],
                paste(DirSave,StrandInd,"/DMRInfo/Hypo/DMRInfo_" ,ChrInd,".txt",sep = ""),
                col.names = TRUE, row.names = FALSE)
  }
}
save.image(file = paste(MainDir, "/Result/EnvironmentP1.RData", sep = ""))

#======================================================================================
# Evaluation section
Method_DMR <- DMRs
Method_DMR <- GenomicRanges::intersect(Method_DMR, islands)

IndRRBS  <- !is.na(islands$dmr.meth.diff)
RRBS_DMR <- islands[IndRRBS,]
start(RRBS_DMR) <- RRBS_DMR$dmr.start;  end(RRBS_DMR) <- RRBS_DMR$dmr.end;

#======================================================================================
Method_DMR$pvalue <- RegPValCal(data.frame(seqnames(Method_DMR),start(Method_DMR),end(Method_DMR)), CpGInfo[,c("Chr","SP","pvalue")])
Method_DMR$status <- as.numeric(as.logical(GenomicRanges::countOverlaps(Method_DMR, RRBS_DMR, ignore.strand=TRUE))) #?

Res <- Evaluation_Criteria()
Res[[1]] <- c(Res[[1]],"RunningTime" = ElapsedTime)
Method_DMR <- islands[Res[[4]]]
RRBS_DMR   <- islands[Res[[3]]]

#======================================================================================
# Save result (venn diagram, ROC, Evaluation measures and DMRs)
CreateFolder(paste(MainDir,"Result/",sep = ""))
png(paste(MainDir,"Result/VD.png",sep = ""),width = 1200, height = 800,  res = 200)
ChIPpeakAnno::makeVennDiagram(list(islands, RRBS_DMR, Method_DMR), NameOfPeaks = c("Islands", "ActualDMR","PredictedDMR"),
                              fill = c("gray90","gray80","gray30"), col = "black", cat.col = "black",
                              cat.pos = c(-10,0,5), cex = 1.25, cat.cex = 1.25, cat.fontface = "bold")
dev.off()

png(paste(MainDir,"Result/ROC.png",sep = ""),width = 1000, height = 1000,  res = 200)
plot(Res[[2]], legend = FALSE, color = FALSE, lwd = 1)
dev.off()

Res[[1]][4:length(Res[[1]])] <- round(Res[[1]][4:length(Res[[1]])],4)
png(paste(MainDir,"Result/Result.png",sep = ""),width = 2000, height = 200,  res = 150)
gridExtra::grid.table(t(Res[[1]]))
dev.off()
#======================================================================================
saveRDS(Res,paste(MainDir,"/Result/DSSRes.rds",sep = ""))
save.image(file = paste(MainDir, "/Result/EnvironmentP2.RData", sep = ""))



