

#' Fits the Maximum Likelihood method to estimate completeness of the fossil record
#' developed by Foote.
#'
#' @param fr a matrix, data.frame, or tbl with 2 columns of per-taxon first and
#'  last occurrence, relative to the modern. If it's a list, it'll try to
#'  extract the matrix from the second entry.
#'
#' @details Runs the analysis of the fossil record completeness developed by:
#'
#' @return a numeric vector with the DR estimated for each tip
#'
#' @author Mauro TC Sugawara, based on David W. Bapst make_durationFreq* functions
#'  from paleotree
#'
#' @references
#' Foote, M. 1997. Paleobiology 23(3):278–300.
#' Foote, M., and D. M. Raup. 1996. Paleobiology 22(2):121–140.
#'
#' @examples
#' set.seed(666)
#' record = paleotree::simFossilRecord(
#' p = 0.3,
#' q = 0.1,
#' nruns = 1,
#' r = 3,
#' nTotalTaxa = c(30, 40)
#' )
#' firstlast = sapply(record, function(x)
#' range(x[["sampling.times"]]))
#' FR = t(apply(firstlast, MARGIN = 2, function(x)
#' ifelse(is.finite(x), x, NA)))
#' res = extract_ML_completeness(FR)
#' res$par
extract_ML_completeness = function(fr) {
  if (inherits(fr, "list")) {
    fr = try(fr[[2]])
  }
  if (inherits(fr, "tbl") | inherits(fr, "data.frame")) {
    fr = as.matrix(fr)
  }
  if (!inherits(fr, "matrix")) {
    stop("'fr' must be a matrix, data.frame, or tbl.
         Or a list with a matrix as second item.")
  }
  fr = na.exclude(fr)
  if (dim(fr)[2] != 2) {
    stop(
      "matrix must have 2 columns of per-taxon first and last occurrence,
      relative to the modern"
    )
  }
  if (any(fr[, 2] >= fr[, 1])) {
    fr = cbind(fr[, 2], fr[, 1])
  }
  if (any(as.numeric(fr) %% 1 != 0)) {
    likFun = paleotree::make_durationFreqCont(fr)
  } else{
    likFun = paleotree::make_durationFreqDisc(fr)
  }
  bounds = attributes(likFun)$parbounds
  start = apply(do.call(rbind, bounds), MARGIN = 2, mean)
  out = optim(
    par = start,
    fn = likFun,
    lower = bounds[[1]],
    upper = bounds[[2]],
    method = "L-BFGS-B",
    control = list(maxit = 1000000)
  )
  if (out$convergence != 0) {
    warning("Optimization did NOT converge!")
  }

  return(out)
  }
