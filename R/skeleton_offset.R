#' Offset a polygon based on it's straight skeleton
#'
#' A straight skeleton is the structure made up of the angular bisectors of each
#' vertex in a polygon (see [skeleton_interior()] for how to get the skeleton).
#' It can be used to offset polygons either inward or outward. This offsetting
#' technique will keep the corners of the polygon _as-is_.
#'
#' @param polygon A `polyclid_polygon` vector. If shorter than `offset` it will
#' be recycled to the length of `offset`
#' @param offset An `euclid_exact_numeric` or numeric vector. If shorter than
#' `polygon` it will be recycled to the length of `polygon`
#'
#' @return A `polyclid_polygon_set` vector
#'
#' @family polygon offsetting
#' @family straight skeleton functions
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
#' ins <- skeleton_offset(poly, -0.05)
#' euclid_plot(ins, lty = 2)
#'
#' # Positive offset
#' off <- skeleton_offset(poly, 0.1)
#' euclid_plot(off, lty = 3)
#'
#' # As can be seen, small sharp angles can become severely exaggerated with
#' # this offset type
#' plot(skeleton_offset(poly[9], c(0.01, 0.02, 0.04, 0.06, 0.1)))
#'
skeleton_offset <- function(polygon, offset) {
  polygon_skeleton_offset(make_valid(polygon), as_exact_numeric(offset))
}

