#' Calculate the minkowski sum of two polygons
#'
#' The minkowski sum is the point-wise sum of two polygons. Practically it means
#' that you trace one polygon around the boundary of the other, expanding the
#' polygon into the covered area.
#'
#' @param P A `polyclid_polygon` vector that should be expanded
#' @param Q A `polyclid_polygon` vector to trace around `P`
#'
#' @return A `polyclid_polygon` vector
#'
#' @family minkowski sum functions
#'
#' @importFrom polyclid make_valid
#' @export
#'
#' @examples
#'
#' poly <- polyclid::denmark()
#'
#' # Expand Jutland by tracing the Island of Bornhold around it
#'
#' # First we should center Bornholm on the origin.
#' bornholm <- transform(poly[7], affine_translate(-vec(centroid(poly[7]))))
#' plot(minkowski_sum(poly[1], bornholm))
#'
minkowski_sum <- function(P, Q) {
  polygon_minkowski_sum(make_valid(P), make_valid(Q))
}
