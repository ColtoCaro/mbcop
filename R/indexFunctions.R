#File with code for creating indices to compare partitions
#'
#'
#' Creating Pairwise indices
#'
#' Function for creating indices based on the 2x2 table of pairwise
#' comparisons. This includes sensitivity, specificity, positive predictive
#' value, negative predictive value and a kappa statistic. The function
#' returns a list object containing all of these pariwise indices.
#'
#'
#' @export
#'
#' @param trueLabels Either a character or numeric vector denoting the known
#'  or hypothesized labels for each observation.
#'
#' @param clustLabels Either a character or numeric vector denoting the
#'  clustering assignments for each observation.
#'
#'
#'
pairwise = function(trueLabels, clustLabels){

  freq <- smry2x2(trueLabels, clustLabels)

  if(is.null(dim(freq))){
    sens = freq[4] / (freq[3] + freq[4])
    spec = freq[1] / (freq[1] + freq[2])
  }
  if(!(is.null(dim(freq)))){
    sens = freq[,4] / (freq[,3] + freq[,4])
    spec = freq[,1] / (freq[,1] + freq[,2])
  }

  labelfreq <- table(truelabels)/length(truelabels)
  prev = labelfreq %*% labelfreq

  ppv <- prev*sens / (prev * sens + (1-prev)*(1-spec))
  npv <- (1-prev)*spec / ((1-prev)*spec + prev*(1-sens))

  kappa <- get_kappa(sens, spec, prev)

  list(sensitivity = sens, specificity = spec,
              ppv = ppv, npv = npv, kappa = kappa)

}

#' Creating indices based on triplets of observations
#'
#' Function for creating indices based on the 5x5 table of triplet
#' comparisons. This includes kappa computed on the 5x5 table, along with
#' five other indices representing submatrix associations.
#'
#'
#' @export
#'
#' @param trueLabels Either a character or numeric vector denoting the known
#'  or hypothesized labels for each observation.
#'
#' @param clustLabels Either a character or numeric vector denoting the
#'  clustering assignments for each observation.
#'
#'
#'
triplets = function(trueLabels, clustLabels){

  xt3u = xtab_utriplets(trueLabels, clustLabels)

  kappa <- kappa_stat(xt3u)

  subKappas <- five_kappas(xt3u)

  list(kappa = kappa, subKappas = subKappas)

}
