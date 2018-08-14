
#' Auxiliary functions for get_DR.
#'
#' @author Mauro TC Sugawara

Get_Path = function(tip, edge, root) {
  out = c()
  while (tip != root) {
    out = append(which(edge[, 2] == tip), out)
    tip = edge[out[1], 1]
  }

  return(rev(out))
}
.get_path = Vectorize(FUN = Get_Path, vectorize.args = "tip")
EstimateES = function(ind, edge.length) {
  lj = edge.length[ind]
  weight = 2 ^ ((1:length(ind)) - 1)
  out = sum(lj * (1 / weight))

  return(out)
}
.estimateES = Vectorize(FUN = EstimateES, vectorize.args = "ind")
