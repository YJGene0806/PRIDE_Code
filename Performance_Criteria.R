###############################################
## Performance criteria
###############################################
## function: ACC.Criterion()
## input XX is a 2 by 1 matrix 
## first row represents predictive differential edge (0/1)
## second row represents true differential edge (0/1)
###############################################

ACC.Criterion <- function(XX){

## predict + true
T1 <- apply(XX, 2, sum)

## predict - true
T2 <- apply(XX, 2, function(input) input[1]-input[2])

TP <- length(which(T1 == 2))

TN <- length(which(T1 == 0))

FP <- length(which(T2 == 1))

FN <- length(which(T2 == -1))

SEN <- TP/(TP+FN)

SPE <- TN/(TN+FP)

ACC <- (TP + TN)/(dim(XX)[2])

FDR <- FP/(FP+TP)

denominator <- as.double((TP+FP))*as.double((TP+FN))*as.double((TN+FP))*as.double((TN+FN))

MCC <- (TP*TN-FP*FN)/sqrt(denominator)

Recall <- TP/(TP+FN)

Precision <- TP/(TP+FP)

F1 <- (2*TP)/(2*TP+FP+FN)

Result <- c(TP, TN, FP, FN, SEN, SPE, ACC, FDR, MCC, Recall, Precision, F1)

names(Result) <- c("TP", "TN", "FP", "FN", "SEN", "SPE", "ACC", "FDR", "MCC", "Recall", "Precision", "F1")

return(Result)

}



