



#' Estimates DR statistics
#'
#' @param phy Ultrametric phylogeny
#'
#' @details From the SM of Jetz et al. 2012: The ES measure for a
#' focal tip on a rooted bifurcating tree is the sum of the edge lengths from the
#' species i to the root, with each consecutive edge discounted by a factor of
#' 1/2:
#'    ES_i = sum_{j = 1}^{Ni} l_j * ( 1/2^(j-1) ),
#'
#' where N_i is the number of edges on the path from species i and the root, and
#' l_j is the length of the edge j, with j=1 being the pendant edge leading to the
#' species and j=N_i being the edge nearest the root on that path. For a more
#' general equation for non-bifurcating trees, see Redding et al. (2008). The
#' inverse of this measure can be seen as a measure of the splitting rate of the
#' path to a tip: species in rapidly-diversifying clades will have short edge
#' lengths shared among many species and low ES values, while isolated species on
#' a tree have no evidence of recent diversification and large ES values. We term
#' the 1/ES metric for a species its species-level lineage diversification rate,
#' or DR.
#'
#' @return a numeric vector with the DR estimated for each tip
#'
#' @author Mauro TC Sugawara
#'
#' @references
#'  W. Jetz, G. H. Thomas, J. B. Joy, K. Hartmann, A. O. Mooers. The global diversity of birds in space and time. Nature, 491: 444-448. doi:10.1038/nature11631
#'
#' @examples
#' tree = rcoal(30)
#' get_DR(tree)
#'
#' @export
#' @import ape
get_DR = function(phy) {
  if (!inherits(phy, "phylo")) {
    stop("'phy' must be of class phylo.")
  }
  nt = Ntip(phy)
  if (nt != (phy$Nnode + 1)) {
    stop("'phy' must be fully resolved.")
  }
  if (!is.ultrametric(phy)) {
    stop("'phy' must be ultrametric.")
  }
  path = .get_path(1:nt, phy$edge, nt + 1)
  ES = .estimateES(path, edge.length = phy$edge.length)
  DR = 1 / ES
  names(DR) = phy$tip.label

  return(DR)
}
