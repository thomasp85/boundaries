#' Offset polygon using minkowski sum with a disc
#'
#' The minkowski sum of two polygons is the polygon you get when you trace one
#' polygon around the boundary of the other and adding the covered area to the
#' boundary. As such it provides an offset with rounded corners for obtuse
#' angles. The algorithm provided here is an approximation with the allowed
#' error bound being controlled by the `eps` argument. Further, the introduced
#' corner arcs are estimations of the disc as they are converted to line
#' segments.
#'
#' @param polygon A `polyclid_polygon` vector. If shorter than `offset` it will
#' be recycled to the length of `offset`
#' @param offset An `euclid_exact_numeric` or numeric vector. If shorter than
#' `polygon` it will be recycled to the length of `polygon`
#' @param arc_segments The number of segments used to draw a full circle
#' @param eps The error bound of the approximation
#'
#' @return A `polyclid_polygon_set` vector
#'
#' @family polygon offsetting
#' @family minkowski sum functions
#'
#' @importFrom polyclid make_valid
#' @importFrom euclid as_exact_numeric
#' @export
#'
#' @examples
#' poly <- polyclid::denmark()
#'
#' plot(poly, col = "grey")
#'
#' # Negative offset (inset)
#' ins <- minkowski_offset(poly, -0.05, arc_segments = 10)
#' euclid_plot(ins, lty = 2)
#'
#' # Positive offset
#' off <- minkowski_offset(poly, 0.1, arc_segments = 10)
#' euclid_plot(off, lty = 3)
#'
minkowski_offset <- function(polygon, offset, arc_segments = 50, eps = 0.00001) {
  arc_segments <- as.integer(arc_segments)
  if (any(is.na(arc_segments) | arc_segments < 1)) {
    cli_abort("{.arg arc_segments} must be positive integers")
  }
  eps <- as.numeric(eps)
  if (length(eps) != 1 || !is.finite(eps) || eps < 0) {
    cli_abort("{.arg eps} must be a scalar positive numeric")
  }
  polygon_minkowski_offset(make_valid(polygon), as_exact_numeric(offset), arc_segments, eps)
}
