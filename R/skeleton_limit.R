#' Limit of straight skeletons
#'
#' The limit of a straight skeleton is the distance at which the inner offset
#' disappears completely. You can get the geometry left at the limit as a
#' polyline set consisting of isolated vertices and/or segments.
#'
#' @param polygon A `polyclid_polygon` vector
#'
#' @return `skeleton_limit` returns a `euclid_exact_numeric` vector and
#' `skeleton_limit_location` returns a `polyclid_polyline_set` vector
#'
#' @family straight skeleton functions
#'
#' @importFrom polyclid make_valid
#' @export
#'
#' @examples
#'
#' poly <- polyclid::denmark()[9]
#' # Complex polygons will often have a single limit point
#'
#' plot(poly)
#' euclid_plot(vert(skeleton_limit_location(poly)))
#'
#' # But certain geometries will result in segments rather than lines
#' poly <- polygon(
#'   c(1, 1, 2, 2, -2, -2, -1, -1),
#'   c(-2, 0, 0, 2, 2, 0, 0, -2)
#' )
#' euclid_plot(skeleton_limit_location(poly))
#'
#' # You can get the distance from the limit location to the boundary
#' skeleton_limit(poly)
#'
skeleton_limit <- function(polygon) {
  polygon_skeleton_limit(make_valid(polygon))
}

#' @rdname skeleton_limit
#' @importFrom polyclid make_valid
#' @export
skeleton_limit_location <- function(polygon) {
  polygon_skeleton_limit_location(make_valid(polygon))
}
