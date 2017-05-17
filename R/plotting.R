#File for ploting function
#'
#'
#' Hierarchical Reciever Operating Characteristic
#'
#' Function for creating HROC plots.  These plots show clustering sensitivity
#' on the y-axis and 1 - clustering specificity on the x-axis.  Changes in the
#' values correspond to vertical increments along the hierarchical clustering
#' dengrogram.  Hierarchical clustering is done via the agnes function from
#' the cluster package.  All of the agnes function parameters can be passed
#' throught the hroc function.
#'
#' @export
#'
#' @param dat An NxP data matrix where the P columns are going to be clustered.
#'
#' @param labelID A vector of length P giving integer or text labels for each
#'  column.  NA's are allowed which results in the ROC plot being made
#'  based on the clustering of the full data but only using the non-NA labels
#'  in computing sensitivity and specificity.
#'
#'
hroc <- function(dat, labelID, ...){
  #Columns will be clustered.  Accepts options for the Agnes function.
  #First reduce the data for efficient processing
  if(nrow(dat) > ncol(dat)){
  dat <- ldr(dat)
  }

  hc <- agnes(t(dat), ...)

  mergeTable <- hc$merge
  distVec <- hc$height[order(hc$height)]

  #create vector for observations with known labels
  if(sum(is.na(labelID)) > 0){
    realObs <- (1:length(labelID))[-which(is.na(labelID))]
    realLabs <- labelID[-which(is.na(labelID))]
  }else{
    realObs <- (1:length(labelID))
    realLabs <- labelID
  }
  #create matrix with all possible pairs and a linkage vector
  pairMat <- t(combn(realObs, 2))
  linked <- labelID[pairMat[ , 1]] == labelID[pairMat[ , 2]]

  #Finally get the distance at which each point merged
  pairDist <- apply(pairMat, 1, function(x)
    getDist(x[1], x[2], mergeTable, distVec))

  pred <- ROCR::prediction(-pairDist, linked)
  perf <- ROCR::performance(pred, measure = "tpr", x.measure = "fpr")

  plot(perf, col="red", main = paste("HROC with", hc$method, "linkage"))
  df <- data.frame(pairMat, linked, pairDist)
  df
}


getDist <- function(x1, x2, mergeTable, distVec){
  #Function to find the distance at which two points (x1, x2) are merged
  i <- 1
  loc1 <- -1*x1
  loc2 <- -1*x2
  while(loc1 != loc2){
    if(mergeTable[i, 1] == loc1 | mergeTable[i, 2] == loc1){loc1 <- i}
    if(mergeTable[i, 1] == loc2 | mergeTable[i, 2] == loc2){loc2 <- i}
    if(loc1 == loc2){return(distVec[i])}
    i <- i + 1
  }
}

