#' Increase number of vertices along a boundary
#'
#' This function increases the number if vertices in a polygon or polyline by
#' either adding vertices so that no edge is longer than a certain threshold,
#' until the cardinality reaches a certain level, or by subdividing edges a
#' specific number of times. Vertices are always inserted _on_ the boundary so
#' no change in the actual geometry is induced by the operation.
#'
#' @param poly A `polyclid_polygon` or `polyclid_polyline` vector
#' @param max_dist Either `NULL` or the maximum distance between consecutive
#' vertices allowed after the densification
#' @param min_cardinality Either `NULL` or the lowest number of vertices in the
#' resulting polygon. If an input polygon has a higher cardinality it is
#' returned unchanged.
#' @param n_splits Either `NULL` or the number of new vertices to insert into
#' every edge.
#'
#' @return A vector of the same type and length as `poly`
#'
#' @family boundary resolution
#'
#' @importFrom polyclid is_polygon is_polyline
#' @export
#'
#' @examples
#' poly <- polygon(
#'   c(391, 240, 252, 374, 289, 134, 68, 154, 161, 435, 208, 295, 421, 441),
#'   c(374, 431, 340, 320, 214, 390, 186, 259, 107, 108, 148, 160, 212, 303)
#' )
#' plot(poly)
#' euclid_plot(vert(poly))
#'
#' # Add 1 verice to each edges
#' plot(poly)
#' euclid_plot(vert(densify_poly(poly, n_splits = 1)))
#'
#' # Increase cardinality to 40
#' plot(poly)
#' euclid_plot(vert(densify_poly(poly, min_cardinality = 40)))
#'
#' # Make sure edges have at most a length of 100
#' plot(poly)
#' euclid_plot(vert(densify_poly(poly, max_dist = 100)))
#'
densify_poly <- function(poly, max_dist = NULL, min_cardinality = NULL, n_splits = NULL) {
  if (!is_polygon(poly) && !is_polyline(poly)) {
    cli_abort("{.arg poly} must be either a polygon or polyline vector")
  }
  if (sum(c(is.null(max_dist), is.null(min_cardinality), is.null(n_splits))) != 2) {
    cli_abort("Only one of {.arg max_dist}, {.arg min_cardinality}, and {.arg n_splits} may be given")
  }
  if (!is.null(max_dist)) {
    max_dist <- as.numeric(max_dist)
    if (anyNA(max_dist) || any(max_dist <= 0)) {
      cli_abort("{.arg max_dist} must be a positive numeric vector")
    }
    poly_add_detail(poly, 0L, max_dist, 0L, TRUE, FALSE)
  } else if (!is.null(min_cardinality)) {
    min_cardinality <- as.integer(min_cardinality)
    if (anyNA(min_cardinality) || any(min_cardinality < 1)) {
      cli_abort("{.arg min_cardinality} must be a positive integer vector")
    }
    poly_add_detail(poly, min_cardinality, 0.0, 0L, FALSE, FALSE)
  } else {
    n_splits <- as.integer(n_splits)
    if (anyNA(n_splits) || any(n_splits < 0)) {
      cli_abort("{.arg n_splits} must be a positive integer vector")
    }
    poly_add_detail(poly, 0L, 0.0, n_splits, FALSE, TRUE)
  }
}
