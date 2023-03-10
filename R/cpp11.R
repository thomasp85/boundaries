# Generated by cpp11: do not edit by hand

poly_add_detail <- function(poly, c, l, n, use_l, use_n) {
  .Call(`_boundaries_poly_add_detail`, poly, c, l, n, use_l, use_n)
}

poly_corner_cutting <- function(poly, max_angle, max_cut, n_cut) {
  .Call(`_boundaries_poly_corner_cutting`, poly, max_angle, max_cut, n_cut)
}

polygon_minkowski_offset <- function(polygons, offset, n, eps) {
  .Call(`_boundaries_polygon_minkowski_offset`, polygons, offset, n, eps)
}

polygon_minkowski_sum <- function(p, q) {
  .Call(`_boundaries_polygon_minkowski_sum`, p, q)
}

poly_simplify <- function(poly, cost, stop, cost_par, stop_par) {
  .Call(`_boundaries_poly_simplify`, poly, cost, stop, cost_par, stop_par)
}

polygon_skeleton_offset <- function(polygons, offset) {
  .Call(`_boundaries_polygon_skeleton_offset`, polygons, offset)
}

polygon_skeleton_polylineset <- function(polygons, keep_boundary, only_inner) {
  .Call(`_boundaries_polygon_skeleton_polylineset`, polygons, keep_boundary, only_inner)
}

polygon_skeleton_limit <- function(polygons) {
  .Call(`_boundaries_polygon_skeleton_limit`, polygons)
}

polygon_skeleton_limit_location <- function(polygons) {
  .Call(`_boundaries_polygon_skeleton_limit_location`, polygons)
}
