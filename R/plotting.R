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
#'  column.
#'
#'
hroc <- function(dat, labelID, ...){
  #First reduce the data for efficient processing
  smalldat <- ldr(dat)


}
