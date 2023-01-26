#' Simplify polygons and polylines
#'
#' This function allows you to reduce the number of vertices used to encode the
#' geometry while attempting to keep the geometry as close to the original as
#' possible. You can choose between different cost functions and stop
#' conditions in order to fine tune the procedure.
#'
#' @param poly A `polyclid_polygon` or `polyclid_polyline`
#' @param cost The cost function to use. See details.
#' @param stop The stop condition to use for terminating the simplification. See
#' details.
#' @param cost_ratio The ratio to use if `cost = "hybrid squared`
#' @param stop_threshold The threshold to use with the stop condition
#'
#' @return A vector of the same type as `poly`
#'
#' @export
#'
#' @examples
#' poly <- polyclid::denmark()
#'
#' plot(poly)
#' poly_s <- simplify_poly(poly)
#'
#' # Number of original vertices
#' sum(cardinality(poly))
#'
#' # Number after simplification
#' sum(cardinality(poly_s))
#'
#' # Close to same visual
#' plot(poly_s)
#'
#' # Plot Jutland using only 20 vertices
#' plot(simplify_poly(poly[1], stop = "below count", stop_threshold = 20))
#' euclid_plot(poly[1], lty = 2, border = "red")
#'
simplify_poly <- function(poly, cost = "squared", stop = "below count ratio", cost_ratio = 0.5, stop_threshold = 0.5) {
  cost_fun <- c("hybrid squared", "scaled squared", "squared")
  cost <- arg_match0(cost, cost_fun)
  cost <- match(cost, cost_fun) - 1L

  stop_fun <- c("above cost", "below count ratio", "below count")
  stop <- arg_match0(stop, stop_fun)
  stop <- match(stop, stop_fun) - 1L

  cost_ratio <- as.numeric(cost_ratio)
  if (length(cost_ratio) != 1 || !is.finite(cost_ratio) || cost_ratio < 0) {
    cli_abort("{.arg cost_ratio} must be a scalar positive numeric")
  }

  stop_threshold <- as.numeric(stop_threshold)
  if (length(stop_threshold) != 1 || !is.finite(stop_threshold) || stop_threshold < 0) {
    cli_abort("{.arg stop_threshold} must be a scalar positive numeric")
  }
  if (stop == 2 && !is_integerish(stop_threshold)) {
    cli_abort("When {.arg stop} is {.val below count} {.arg stop_threshold} should be an integer")
  }
  poly_simplify(poly, cost, stop, cost_ratio, stop_threshold)
}
