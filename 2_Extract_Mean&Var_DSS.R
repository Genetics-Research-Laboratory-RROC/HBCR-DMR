# Calculate transformed mean and variance with modified function from DSS package
# and save result for next step
#======================================================================================
rm(list = ls())
gc()

#======================================================================================
# Load rrbs data and neccesary packages
data(rrbs, package = "RRBSdata")

library(BiSeq)
library(bsseq)
library(GenomicRanges)

#======================================================================================
# Cjoose main directory
MainDir <- paste(tcltk::tk_choose.dir(default = getwd(), 
                                      caption = "Select HBCR-DMR directory"),"/",sep = "")

#======================================================================================
# load data from two group and save in list for next step 
ChrLoc <- rowRanges(rrbs)
ChrLoc <- data.frame("chr" = seqnames(ChrLoc), "pos" = start(ChrLoc))

N <- totalReads(rrbs)
X <- methReads(rrbs)
SampleType <- as.character(colData(rrbs)$group)
SampleNames <- colnames(rrbs)
# note: totalReads, methReads and rowRanges have same order and identical rownames.
Normal <- list()
Cancer <- list()
for (i in 1:length(SampleNames)) {
  temp <- cbind(ChrLoc, N[,SampleNames[i]], X[,SampleNames[i]])
  colnames(temp) <- c("chr","pos","N","X" )
  if (SampleType[i] == "normal") {
    Normal[[length(Normal)+1]] <- temp
  } else {
    Cancer[[length(Cancer)+1]] <- temp
  }
}

rm(list = c('rrbs','N','X','i','ChrLoc','temp'))

#======================================================================================
# Load and run My_DMLtest function for claculation transformed mean and variance 
source(paste(MainDir,"Function/My_DMLtest.R",sep = ""))

BSobj <- DSS::makeBSseqData(c(Cancer, Normal),SampleNames)
ModCpGTest <- My_DMLtest(BSobj, group1=SampleNames[SampleType=="cancer"], group2=SampleNames[SampleType=="normal"], smoothing=TRUE)

#======================================================================================
# Save result as text file in Result subdirectory
source(paste(MainDir,"Function/CreateFolder.R",sep = ""))
CreateFolder(paste(MainDir,"Result/",sep = ""))
write.table(ModCpGTest, paste(MainDir,"Result/DSS_Output.txt",sep = ""), 
            col.names = TRUE, row.names = FALSE)





