#File for ploting function
#'
#'
#' Hierarchical Reciever Operating Characteristic
#'
#' Function for creating HROC plots.  These plots show the reciever operating
#' characteristic for a given set of labels (possibly incomplete) and a type
#' of hierarchical clustering algorithm.  The distance at which any two
#' points are merged along a dengrogram is used as a predictor of linkage
#' while true pair-wise linkage is based on the labels.
#' Hierarchical clustering is done via the \code{\link{agnes}} function from
#' the \code{\link{cluster}} package.  All of the agnes function parameters
#' can be passed through the hroc function. In addition to creating a plot,
#' this function returns a table with all of the pairwise combinations,
#' along with their respective linkages and distances.
#'
#' @export
#'
#' @param dat An NxP data matrix where the P columns are going to be clustered.
#'
#' @param labelID A list of length P giving a vector of integer or text
#'  labels for each column that was clustered.  Vectors of labels are
#'  accepted to allow observations to belong to multiple categories.  In this
#'  case, pairs of observations will be considered linked if any of the
#'  labels in the vector match any of the labels corresponding to the other
#'  observation.  NA's are allowed which results in the ROC plot being made
#'  based on the clustering of the full data but only using the non-NA labels
#'  in computing sensitivity and specificity.
#'
#' @param methods A vector of strings that determine what sort of linkage
#'  function should be used in the hierarchical clustering.  Parameters
#'  can be passed to the agnes function through "..." but this is separate
#'  because it allows users to put multiple HROC curves on the same plot.
#'  By default single, complete and average linkage functions will be
#'  plotted together.
#'
#'
hroc <- function(dat,
                 labelID,
                 methods = c("single", "complete", "average"),
                 ...){
  if (is.list(labelID) == FALSE){
    labelID <- as.list(labelID)
  }


  if(max(table(unlist(labelID))) <= 1){
    stop("At least 2 observations must share a label")
    }

  if(max(table(unlist(labelID))) == length(labelID)){
    stop("Labels require at least 2 levels")
    }

  #First reduce the data for efficient processing
  if(nrow(dat) > ncol(dat)){
    dat <- ldr(dat)
  }

  hc <- list()
  mergeTable <- list()
  distVec <- list()
  for (i in 1:length(methods)){
    print("Creating clusters")
    hc[[i]] <- cluster::agnes(t(dat), method = methods[i])
    mergeTable[[i]] <- hc[[i]]$merge
    distVec[[i]] <- hc[[i]]$height[order(hc[[i]]$height)]
    plot(hc[[i]], which.plots = 2)
  }


  #create vector for observations with known labels
  anyMissing <- (sum(is.na(unlist(labelID))) > 0)
  if(anyMissing){
    missIndex <- which(sapply(labelID,
                        function(x) sum(is.na(x)) > 0))
    realLabs <- labelID[-missIndex]
  }else{
    realLabs <- labelID
  }
  realObs <- 1:length(realLabs)

  #create matrix with all possible pairs and a linkage vector
  pairMat <- t(combn(realObs, 2))
  linked <- Map('%in%', realLabs[pairMat[ , 1]], realLabs[pairMat[ , 2]])
  linked <- sapply(linked, function(x) sum(x) > 0)

  #Finally get the distances at which each point merged
  
  #Figure out which pairs are not used because of missing labels
  if(anyMissing){
  obsIndex <- (1:length(labelID))[-missIndex]
  nRemove <- length(missIndex) * length(obsIndex)
  rmIndex <- rep(0, nRemove)
  
  vecMiss <- rep(missIndex, each = length(obsIndex))
  vecObs <- rep(obsIndex, length(missIndex))
  bound <- cbind(vecMiss, vecObs)
  badPairs <- rbind(bound, t(combn(missIndex, 2)))
  
  rmIndex <- apply(badPairs, 1, function(x) 
    pairmatIndex(x[1], x[2], length(labelID)))
  }
  
  #now make the distance vector based on the complete dendrogram
  pairDist <- list()
  rocrDist <- list()
  rocrLab <- list()
  for (i in 1:length(methods)){
    pairDist[[i]] <- pushDist(mergeTable[[i]], distVec[[i]])
    if(anyMissing){
      rocrDist[[i]] <- pairDist[[i]][-rmIndex, ] * (-1)
    }else{
      rocrDist[[i]] <- pairDist[[i]] * (-1)
    }
    rocrLab[[i]] <- linked
  }

  pred <- ROCR::prediction(rocrDist, rocrLab)
  perf <- ROCR::performance(pred, measure = "tpr", x.measure = "fpr")
  auc <- ROCR::performance(pred, measure="auc")

  aucs <- signif(unlist(auc@y.values), 2)
  line2 <- paste(methods, "=", aucs, collapse = ", ")

  title <- paste("HROC \n", "linkage AUC's:", line2 )

  if(length(methods) > 1){
  plot(perf@x.values[[1]], perf@y.values[[1]], col=1, main = title,
       type = "l", xlab="FPR", ylab="TPR")

    for(i in 2:length(methods)){
      lines(perf@x.values[[i]], perf@y.values[[i]], col=i,
             main = title, lty = i,
             type = "l", xlab="FPR", ylab="TPR")
    }
  }else{
    plot(perf@x.values[[1]], perf@y.values[[1]], col=1, main = title,
         type = "l", xlab="FPR", ylab="TPR")
  }

  legend("bottomright", methods, lty=c(1:length(methods)),
         col = c(1:length(methods)))

  pairDf <- -1*do.call(cbind, rocrDist)
  colnames(pairDf) <- paste(methods, "Dist", sep = "")
  df <- data.frame(pairMat, linked, pairDf)
  df
}



