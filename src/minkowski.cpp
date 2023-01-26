#include <vector>
#include <utility>
#include <algorithm>

#include <cpp11/declarations.hpp>
#include <cpp11/integers.hpp>

#include <CGAL/approximated_offset_2.h>
#include <CGAL/Gps_circle_segment_traits_2.h>
#include <CGAL/minkowski_sum_2.h>

#include <euclid.h>
#include <polyclid.h>

typedef CGAL::Gps_circle_segment_traits_2<Kernel> Traits;
typedef Traits::Polygon_2 Circ_polygon_2;
typedef Traits::Polygon_with_holes_2 Circ_polygon_with_holes_2;

Polyline circ_to_linear_polygon(Circ_polygon_2& poly, int n) {
  std::vector<Point_2> ring;
  std::vector< std::pair<double, double> > approx_segment;

  for (auto iter = poly.curves_begin(); iter != poly.curves_end(); iter++) {
    approx_segment.clear();
    int n_real = 2;
    if (iter->is_circular()) {
      Vector_2 vec1(iter->supporting_circle().center(), Point_2(CGAL::to_double(iter->source().x()), CGAL::to_double(iter->source().y())));
      Vector_2 vec2(iter->supporting_circle().center(), Point_2(CGAL::to_double(iter->target().x()), CGAL::to_double(iter->target().y())));
      double y = CGAL::sqrt(CGAL::to_double((vec1 - vec2).squared_length().exact()));
      double x = CGAL::sqrt(CGAL::to_double((vec1 + vec2).squared_length().exact()));
      n_real = 1.0 + n - n * (std::atan2(y, x) / 6.283185307179586);
    }
    iter->approximate(std::back_inserter(approx_segment), n_real);
    if (approx_segment.size() == 0) continue;
    if (ring.empty()) ring.emplace_back(approx_segment[0].first, approx_segment[0].second);

    for (size_t i = 1; i < approx_segment.size(); ++i) {
      ring.emplace_back(approx_segment[i].first, approx_segment[i].second);
    }
  }
  ring.pop_back();
  return {ring.begin(), ring.end()};
}

Polygon_set circ_to_linear_polygon(Circ_polygon_with_holes_2& poly, int n) {
  Polygon res;

  if (!poly.is_unbounded()) {
    res = Polygon(circ_to_linear_polygon(poly.outer_boundary(), n));
  }

  for (auto iter = poly.holes_begin(); iter != poly.holes_end(); iter++) {
    res.add_hole(circ_to_linear_polygon(*iter, n));
  }

  return Polygon_set(res);
}

[[cpp11::register]]
SEXP polygon_minkowski_offset(SEXP polygons, SEXP offset, cpp11::integers n, double eps) {
  std::vector<Polygon> poly = polyclid::get_polygon_vec(polygons);
  std::vector<Exact_number> os = euclid::get_exact_numeric_vec(offset);

  size_t max_size = std::max(poly.size(), os.size());
  std::vector<Polygon_set> res;
  res.reserve(max_size);

  for (size_t i = 0; i < max_size; ++i) {
    Polygon p = poly[i % poly.size()];
    Exact_number of = os[i % os.size()];
    if (!(p.get_flag(VALIDITY_CHECKED) && p.get_flag(IS_VALID))) {
      cpp11::stop("Input polygons must be valid");
    }
    if (p.is_na() || of.is_na()) {
      res.push_back(Polygon_set::NA_value());
      continue;
    }
    if (of > 0) {
      Circ_polygon_with_holes_2 offset = CGAL::approximated_offset_2(p, of, eps);
      res.push_back(circ_to_linear_polygon(offset, n[i % n.size()]));
    } else {
      std::vector<Circ_polygon_2> offset;
      res.push_back(Polygon_set());
      if (!p.is_unbounded()) {
        CGAL::approximated_inset_2(p.outer_boundary(), -of, eps, std::back_inserter(offset));
        for (auto iter = offset.begin(); iter != offset.end(); iter++) {
          res.back().insert(circ_to_linear_polygon(*iter, n[i % n.size()]));
        }
      }
      for (auto iter = p.holes_begin(); iter != p.holes_end(); iter++) {
        Segment_trait::Polygon_2 hole(*iter);
        hole.reverse_orientation();
        Circ_polygon_with_holes_2 hole_offset = CGAL::approximated_offset_2(hole, -of, eps);
        res.back().difference(circ_to_linear_polygon(hole_offset, n[i % n.size()]));
      }
    }
  }

  return polyclid::create_polygon_set_vec(res);
}

[[cpp11::register]]
SEXP polygon_minkowski_sum(SEXP p, SEXP q) {
  std::vector<Polygon> P = polyclid::get_polygon_vec(p);
  std::vector<Polygon> Q = polyclid::get_polygon_vec(q);
  size_t max_size = std::max(P.size(), Q.size());
  std::vector<Polygon> res;
  res.reserve(max_size);

  for (size_t i = 0; i < max_size; ++i) {
    if (P[i % P.size()].is_na() || Q[i % Q.size()].is_na()) {
      res.push_back(Polygon::NA_value());
      continue;
    }
    res.push_back(CGAL::minkowski_sum_2(P[i % P.size()], Q[i % Q.size()]));
  }

  return polyclid::create_polygon_vec(res);
}
