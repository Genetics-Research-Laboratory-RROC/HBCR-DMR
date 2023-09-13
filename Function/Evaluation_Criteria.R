Evaluation_Criteria <- function() {
    #===================================================================================
    IndMethod  <- as.logical(GenomicRanges::countOverlaps(islands, Method_DMR, ignore.strand=TRUE))
    
    PPIndex <- IndMethod == TRUE  & IndRRBS == TRUE
    PNIndex <- IndMethod == FALSE & IndRRBS == TRUE
    PPRRBS_DMR <- islands[PPIndex];   start(PPRRBS_DMR) <- PPRRBS_DMR$dmr.start;  end(PPRRBS_DMR) <- PPRRBS_DMR$dmr.end;
    PPIndex[PPIndex == TRUE] <- as.logical(GenomicRanges::countOverlaps(PPRRBS_DMR, Method_DMR))
    IndRRBS <- PPIndex | PNIndex
    
    Res <- NULL
    TP <- Res["TP"] <- as.numeric(sum(IndMethod == TRUE & IndRRBS == TRUE ))
    FP <- Res["FP"] <- as.numeric(sum(IndMethod == TRUE & IndRRBS == FALSE))
    
    FN <- Res["FN"] <- as.numeric(sum(IndMethod == FALSE & IndRRBS == TRUE ))
    TN <- Res["TN"] <- as.numeric(sum(IndMethod == FALSE & IndRRBS == FALSE))
    
    #============= Recall ===================
    Res["Recall"] <- TP / (TP + FN)
    #============= Precision ================
    Res["Precision"] <- TP / (TP + FP)
    #============= AUC ======================
    PVal <- Method_DMR$pvalue
    roccurve <- PRROC::roc.curve(scores.class0 = 1-PVal, weights.class0 = Method_DMR$status, curve = TRUE)
    Res["AUC"] <- roccurve$auc
    #============= MCC ======================
    Res["MCC"] <- ((TP*TN)-(FP*FN)) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    #============= F1-score =================
    Res["F1"]  <- 2 * ((Res["Precision"] * Res["Recall"]) / (Res["Precision"] + Res["Recall"]))
    #============= Specificity ==============
    Res["Specificity"] <- TN / (TN + FP)
    #============= Accuracy =================
    Res["Acc"] <- (TP + TN) / (TP + TN + FP + FN)
    #============= NPV ======================
    Res["NPV"] <- TN / (TN + FN)
    #============= PPV ======================
    Res["PPV"] <- TP / (TP + FP)
    
    return(list(Res, roccurve, IndRRBS, IndMethod))
}



