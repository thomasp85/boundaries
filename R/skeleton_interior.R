#' Calculate the interior straight skeleton
#'
#' While straight skeletons are mostly used for offsetting, they can also in
#' themselves be useful, e.g. for finding the central spine of a polygon. While
#' the boundaries package does not provide a class for representing the straight
#' skeleton itself (ie. the internal CGAL representation), you can get a
#' representation of it as a polyline_set.
#'
#' @param polygon A `polyclid_polygon` vector
#' @param keep_boundary Should the boundary (ie. the input polygon) be returned
#' as part of the skeleton. This is useful if you are mainly interested in using
#' the straight skeleton for partitioning the polygon based on its skeleton
#' @param only_inner Should the returned skeleton only contain the inner
#' bisectors of the skeleton (ie the parts not connected to the vertices of the
#' input polygon). This effectively returns the spine of the polygon
#'
#' @return A `polyclid_polyline_set` vector
#'
#' @family straight skeleton functions
#'
#' @importFrom polyclid make_valid
#' @export
#'
#' @examples
#'
#' poly <- polyclid::denmark()[9]
#' plot(poly)
#'
#' # Only return the bisectors
#' plot(skeleton_interior(poly))
#'
#' # Only return the inner bisectors
#' plot(skeleton_interior(poly, only_inner = TRUE))
#'
#' # Return both the skeleton and boundary
#' plot(skeleton_interior(poly, keep_boundary = TRUE))
#'
skeleton_interior <- function(polygon, keep_boundary = FALSE, only_inner = FALSE) {
  if (!is_logical(only_inner, 1L)) {
    cli_abort("{.arg only_inner} must be a scalar logical")
  }
  if (!is_logical(keep_boundary, 1L)) {
    cli_abort("{.arg keep_boundary} must be a scalar logical")
  }
  if (keep_boundary && only_inner) {
    cli_warn("{.arg keep_boundary} is ignored when {.code only_inner = TRUE}")
  }
  polygon_skeleton_polylineset(make_valid(polygon), keep_boundary, only_inner)
}
