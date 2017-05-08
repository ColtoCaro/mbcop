#File for ploting function
#'
#'
#' Hierarchical Reciever Operating Characteristic
#'
#' Function for creating HROC plots.  These plots show clustering sensitivity
#' on the y-axis and 1 - clustering specificity on the x-axis.  Changes in the
#' values correspond to vertical increments along the hierarchical clustering
#' dengrogram.
#'
#' @export
#'
#' @param dat An NxP data matrix where the P columns are going to be clustered.
#'
#' @param labelID A vector of length P giving integer or text labels for each
#'  column.
#'
#' @param linkage Linkage function for deciding how to group similar vectors.
#'
hroc <- function(dat, labelID, linkage){
  #First reduce the data for efficient processing
  smalldat <- ldr(dat)


}
