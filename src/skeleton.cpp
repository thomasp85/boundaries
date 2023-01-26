#include <vector>
#include <algorithm>

#include <cpp11/declarations.hpp>

#include <CGAL/create_offset_polygons_from_polygon_with_holes_2.h>
#include <CGAL/create_straight_skeleton_from_polygon_with_holes_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Cartesian_converter.h>

#include <euclid.h>
#include <polyclid.h>

#include <boost/shared_ptr.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel I_Kernel;
typedef CGAL::Cartesian_converter<I_Kernel,Kernel> I_to_E;
typedef CGAL::Cartesian_converter<Kernel,I_Kernel> E_to_I;
typedef CGAL::Polygon_with_holes_2<I_Kernel> I_Polygon;

I_Polygon exact_to_inexact_poly(const Polygon& poly) {
  static E_to_I converter;
  I_Polygon ipoly;
  std::vector<I_Kernel::Point_2> ring;
  if (!poly.is_unbounded()) {
    std::transform(poly.outer_boundary().vertices_begin(), poly.outer_boundary().vertices_end(),
                   std::back_inserter(ring),
                   [](const Point_2& p) { return converter(p); });
    ipoly = I_Polygon(CGAL::Polygon_2<I_Kernel>(ring.begin(), ring.end()));
  }
  for (auto iter = poly.holes_begin(); iter != poly.holes_end(); iter++) {
    ring.clear();
    std::transform(iter->vertices_begin(), iter->vertices_end(),
                   std::back_inserter(ring),
                   [](const Point_2& p) { return converter(p); });
    ipoly.add_hole(CGAL::Polygon_2<I_Kernel>(ring.begin(), ring.end()));
  }
  return ipoly;
}
Polygon inexact_to_exact_poly(const I_Polygon& ipoly) {
  static I_to_E converter;
  Polygon poly;
  std::vector<Point_2> ring;
  if (!ipoly.is_unbounded()) {
    std::transform(ipoly.outer_boundary().vertices_begin(), ipoly.outer_boundary().vertices_end(),
                   std::back_inserter(ring),
                   [](const I_Kernel::Point_2& p) { return converter(p); });
    poly = Polygon(Segment_trait::Polygon_2(ring.begin(), ring.end()));
  }
  for (auto iter = ipoly.holes_begin(); iter != ipoly.holes_end(); iter++) {
    ring.clear();
    std::transform(iter->vertices_begin(), iter->vertices_end(),
                   std::back_inserter(ring),
                   [](const I_Kernel::Point_2& p) { return converter(p); });
    poly.add_hole(Segment_trait::Polygon_2(ring.begin(), ring.end()));
  }
  return poly;
}

[[cpp11::register]]
SEXP polygon_skeleton_offset(SEXP polygons, SEXP offset) {
  std::vector<Polygon> poly = polyclid::get_polygon_vec(polygons);
  std::vector<Exact_number> os = euclid::get_exact_numeric_vec(offset);

  size_t max_size = std::max(poly.size(), os.size());
  std::vector<Polygon_set> res;
  res.reserve(max_size);

  E_to_I converter;

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
    std::vector< boost::shared_ptr<I_Polygon> > offset;
    if (of > 0) {
      offset = CGAL::create_exterior_skeleton_and_offset_polygons_with_holes_2(converter(of), exact_to_inexact_poly(p));
    } else {
      offset = CGAL::create_interior_skeleton_and_offset_polygons_with_holes_2(converter(-of), exact_to_inexact_poly(p));
    }
    Polygon_set offsetted_polygon;
    res.push_back(offsetted_polygon);
    for (auto iter = offset.begin(); iter != offset.end(); iter++) {
      res.back().insert(inexact_to_exact_poly(*iter->get()));
    }
  }

  return polyclid::create_polygon_set_vec(res);
}

[[cpp11::register]]
SEXP polygon_skeleton_polylineset(SEXP polygons, bool keep_boundary, bool only_inner) {
  std::vector<Polygon> poly = polyclid::get_polygon_vec(polygons);
  std::vector<Polyline_set> res;
  res.reserve(poly.size());

  I_to_E converter;

  for (size_t i = 0; i < poly.size(); ++i) {
    Polygon p = poly[i];
    if (!(p.get_flag(VALIDITY_CHECKED) && p.get_flag(IS_VALID))) {
      cpp11::stop("Input polygons must be valid");
    }
    if (p.is_na() || p.is_unbounded()) {
      res.push_back(Polyline_set::NA_value());
      continue;
    }

    boost::shared_ptr< CGAL::Straight_skeleton_2<I_Kernel> > skeleton = CGAL::create_interior_straight_skeleton_2(exact_to_inexact_poly(p));
    Polyline_set skeleton_lines;
    res.push_back(skeleton_lines);
    for (auto iter = skeleton->halfedges_begin(); iter != skeleton->halfedges_end(); iter++) {
      if ((iter->id()%2)==0) {
        if (!keep_boundary && !iter->is_bisector()) {
          continue;
        }
        if (only_inner && !iter->is_inner_bisector()) {
          continue;
        }
        res.back().insert_non_overlapping(Segment_2(converter(iter->vertex()->point()), converter(iter->opposite()->vertex()->point())));
      }
    }
  }

  return polyclid::create_polyline_set_vec(res);
}

[[cpp11::register]]
SEXP polygon_skeleton_limit(SEXP polygons) {
  std::vector<Polygon> poly = polyclid::get_polygon_vec(polygons);
  std::vector<Exact_number> res;
  res.reserve(poly.size());

  for (size_t i = 0; i < poly.size(); ++i) {
    Polygon p = poly[i];
    if (!(p.get_flag(VALIDITY_CHECKED) && p.get_flag(IS_VALID))) {
      cpp11::stop("Input polygons must be valid");
    }
    if (p.is_na() || p.is_unbounded()) {
      res.push_back(Exact_number::NA_value());
      continue;
    }

    boost::shared_ptr< CGAL::Straight_skeleton_2<I_Kernel> > skeleton = CGAL::create_interior_straight_skeleton_2(exact_to_inexact_poly(p));
    auto max_time = skeleton->vertices_begin()->time();
    for (auto iter = skeleton->vertices_begin(); iter != skeleton->vertices_end(); iter++) {
      if (iter->is_skeleton() && iter->time() > max_time) {
        max_time = iter->time();
      }
    }
    res.push_back(max_time);
  }

  return euclid::create_exact_numeric_vec(res);
}

[[cpp11::register]]
SEXP polygon_skeleton_limit_location(SEXP polygons) {
  std::vector<Polygon> poly = polyclid::get_polygon_vec(polygons);
  std::vector<Polyline_set> res;
  res.reserve(poly.size());

  I_to_E converter;

  for (size_t i = 0; i < poly.size(); ++i) {
    Polygon p = poly[i];
    if (!(p.get_flag(VALIDITY_CHECKED) && p.get_flag(IS_VALID))) {
      cpp11::stop("Input polygons must be valid");
    }
    if (p.is_na() || p.is_unbounded()) {
      res.push_back(Polyline_set::NA_value());
      continue;
    }

    boost::shared_ptr< CGAL::Straight_skeleton_2<I_Kernel> > skeleton = CGAL::create_interior_straight_skeleton_2(exact_to_inexact_poly(p));
    auto max_time = skeleton->vertices_begin()->time();
    for (auto iter = skeleton->vertices_begin(); iter != skeleton->vertices_end(); iter++) {
      if (iter->is_skeleton() && iter->time() > max_time) {
        max_time = iter->time();
      }
    }
    res.push_back(Polyline_set());
    for (auto iter = skeleton->halfedges_begin(); iter != skeleton->halfedges_end(); iter++) {
      if ((iter->id()%2)==0 && iter->is_bisector()) {
        bool prim_include = iter->vertex()->time() == max_time;
        bool sec_include = iter->opposite()->vertex()->time() == max_time;
        if (prim_include && sec_include) {
          res.back().insert_non_overlapping(Segment_2(converter(iter->vertex()->point()), converter(iter->opposite()->vertex()->point())));
        } else {
          if (prim_include) {
            CGAL::insert_point(res.back(), converter(iter->vertex()->point()));
          }
          if (sec_include) {
            CGAL::insert_point(res.back(), converter(iter->opposite()->vertex()->point()));
          }
        }
      }
    }
  }

  return polyclid::create_polyline_set_vec(res);
}
