# Load rrbs data from RRBSdata package and change format for input DSS package
#======================================================================================
rm(list = ls())
gc()

#======================================================================================
# Load package and data (rrbs from RRBSdata package), choose main directory and source function
library(RRBSdata)
data("rrbs")

MainDir <- paste(tcltk::tk_choose.dir(default = getwd(), 
                                      caption = "Select HBCR-DMR directory"),"/",sep = "")
source(paste(MainDir,"Function/CreateFolder.R",sep = ""))

#======================================================================================
# extract initial info of data from data object
TypeSample  <- as.character(rrbs$group)
NameSample  <- colnames(rrbs)
Strand <- as.character(unique(strand(rrbs@rowRanges)))

#======================================================================================
# Create necessry directory for saving data
CreateFolder(paste(MainDir, "Data/", sep = ""))
CurrentDir <- paste(MainDir, "Data/CpGInfo/", sep = "")
CreateFolder(CurrentDir)

for (TypeIndex in unique(TypeSample)) {
  SubDir <- paste(CurrentDir,TypeIndex,"/",sep = "")
  CreateFolder(SubDir)
  for (NameIndex in NameSample[which(TypeSample == TypeIndex)]) {
    SubDir<- paste(CurrentDir,TypeIndex,"/",NameIndex,"/",sep = "")
    CreateFolder(SubDir)
    for (StrandIndex in Strand) {
      SubDir<- paste(CurrentDir,TypeIndex,"/",NameIndex,"/",StrandIndex,"/",sep = "")
      CreateFolder(SubDir)
    }
  }
}

#======================================================================================
# Extract necessary fields from rrbs data for DSS package and save these
Chr    <- seqnames(rrbs@rowRanges)
SP     <- start(rrbs@rowRanges)
N      <- totalReads(rrbs)
X      <- methReads(rrbs)
MRatio <- X / N  # Methylation percentage (call)

for (TypeIndex in unique(TypeSample)) {
  for (NameIndex in NameSample[which(TypeSample == TypeIndex)]) {
    for (StrandIndex in Strand){
      SaveDir <- paste(CurrentDir,TypeIndex,"/",NameIndex,"/",StrandIndex,"/",sep = "")
      cat(paste(TypeIndex,NameIndex,"\n\n",sep = " "))
      for (chri in unique(as.character(Chr))) {
        cat(chri,"\n")
        temp   <- which(Chr == chri)
        SPi    <- SP[temp]
        MCi    <- MRatio[temp,NameIndex]
        COVi   <- N[temp,NameIndex]
        MRatioi <- MCi  # methylation ratio smoothed from methylation call (remove own smoothing)
        temp <- data.frame("SP" = SPi, "MC" = MCi, "Cov" = COVi, "MRatio" = MRatioi)
        temp <- temp[complete.cases(temp),]  # remove cpg with cov = 0
        write.table(temp,paste(SaveDir,NameIndex,"_CPGInfo",StrandIndex,"_",chri,".txt",sep = ""), 
                    col.names = TRUE, row.names = FALSE)
      }
    }
  }
}


