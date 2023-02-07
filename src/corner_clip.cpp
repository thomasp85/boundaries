#include <vector>
#include <utility>

#include <cpp11/declarations.hpp>
#include <cpp11/integers.hpp>

#include <euclid.h>
#include <polyclid.h>

#include <list>

struct Edge {
  Point_2 source;
  double length;
  bool corner;
  double cut_dist_back;
  double cut_dist_forward;
};

std::vector<Point_2> cut_corners(std::list<Edge>& edges, size_t n_cut, bool wrap) {
  while (n_cut != 0) {
    auto iter = edges.begin();
    while (iter != edges.end()) {
      if (!iter->corner) {
        iter++;
        continue;
      }
      auto next_corner = std::next(iter, 1);
      while (next_corner != edges.end()) {
        if (next_corner->corner) break;
        next_corner++;
      }
      auto back_iter = std::prev(iter == edges.begin() ? edges.end() : iter, 1);
      double len = 0.0;
      while (len + back_iter->length <= iter->cut_dist_back) {
        len += back_iter->length;
        if (back_iter == edges.begin()) {
          if (!wrap) break;
          back_iter = edges.end();
        }
        back_iter--;
      }
      double length_reduction = 1.0 - (iter->cut_dist_back - len) / back_iter->length;
      auto beyond = std::next(back_iter, 1);
      if (wrap && beyond == edges.end()) beyond = edges.begin();
      Point_2 new_source = back_iter->source + Vector_2(back_iter->source, beyond->source) * Kernel::FT(length_reduction);
      back_iter->length *= length_reduction;

      auto front_iter = iter;
      len = 0.0;
      while (len + front_iter->length <= iter->cut_dist_forward) {
        len += front_iter->length;
        front_iter++;
        if (front_iter == edges.end()) {
          if (!wrap) {
            front_iter = std::prev(front_iter, 2);
            break;
          }
          front_iter = edges.begin();
        }
      }
      length_reduction = (iter->cut_dist_forward - len) / front_iter->length;
      beyond = std::next(front_iter, 1);
      if (wrap && beyond == edges.end()) beyond = edges.begin();
      Point_2 new_target = front_iter->source + Vector_2(front_iter->source, beyond->source) * Kernel::FT(length_reduction);
      front_iter->source = new_target;
      front_iter->length *= 1.0 - length_reduction;
      front_iter->corner = true;
      front_iter->cut_dist_forward = iter->cut_dist_forward * 0.5;
      double new_length = CGAL::sqrt(CGAL::to_double((new_target - new_source).squared_length().exact()));
      auto new_edge = edges.insert(std::next(back_iter, 1), {
        new_source,
        new_length,
        true,
        iter->cut_dist_back * 0.5,
        new_length * 0.25
      });
      front_iter->cut_dist_back = new_length * 0.25;
      auto edge_iter = std::next(new_edge, 1);
      if (wrap && edge_iter == edges.end()) edge_iter = edges.begin();
      while (edge_iter != front_iter) {
        auto remove_edge = edge_iter;
        edge_iter++;
        if (edge_iter == edges.end()) {
          if (!wrap) break;
          edge_iter = edges.begin();
        }
        edges.erase(remove_edge);
      }
      iter = next_corner;
    }
    n_cut--;
  }

  std::vector<Point_2> new_poly;
  std::transform(edges.begin(), edges.end(), std::back_inserter(new_poly), [](const Edge& e) { return e.source; });

  return new_poly;
}

std::list<Edge> create_edge_ring(const Segment_trait::Polygon_2& poly, double max_angle, double max_cut) {
  std::list<Edge> ring;
  auto segment_start = poly.edges_circulator();
  auto s_iter = segment_start;
  do {
    auto last_edge = std::prev(s_iter, 1);
    auto last_len = ring.size() == 0 ? CGAL::sqrt(CGAL::to_double(last_edge->squared_length().exact())) : ring.back().length;
    double cur_len = CGAL::sqrt(CGAL::to_double(s_iter->squared_length().exact()));
    Vector_2 a = Vector_2(*last_edge) * Kernel::FT(cur_len);
    Vector_2 b = Vector_2(*s_iter) * Kernel::FT(last_len);
    double angle = 3.14159265359 - std::atan2(CGAL::sqrt(CGAL::to_double((a - b).squared_length().exact())),
                                              CGAL::sqrt(CGAL::to_double((a + b).squared_length().exact())));
    ring.push_back({s_iter->source(), cur_len, angle < max_angle, 0.0, 0.0});
  } while (++s_iter != segment_start);

  double last_cut = 0.0;
  auto corner = ring.end();
  for (auto iter = ring.begin(); iter != ring.end(); iter++) {
    if (iter->corner) {
      if (corner == ring.end()) corner = iter;
      iter->cut_dist_back = last_cut;
      double dist_to_next = 0.0;
      auto forward_iter = iter;
      do {
        dist_to_next += forward_iter->length;
        forward_iter++;
        if (forward_iter == ring.end()) {
          forward_iter = ring.begin();
        }
      } while (!forward_iter->corner);
      last_cut = std::min(dist_to_next * 0.5, max_cut) * 0.5;
      iter->cut_dist_forward = last_cut;
    }
  }

  // No corners
  if (corner != ring.end()) corner->cut_dist_back = last_cut;

  return ring;
}

std::list<Edge> create_edge_line(const Polyline& poly, double max_angle, double max_cut) {
  std::list<Edge> line;
  for (auto iter = poly.edges_begin(); iter != poly.edges_end(); iter++) {
    double cur_len = CGAL::sqrt(CGAL::to_double(iter->squared_length().exact()));
    if (line.size() == 0) {
      line.push_back({iter->source(), cur_len, false, 0.0, 0.0});
    } else {
      auto last_edge = std::prev(iter, 1);
      auto last_len = line.back().length;
      Vector_2 a = Vector_2(*last_edge) * Kernel::FT(cur_len);
      Vector_2 b = Vector_2(*iter) * Kernel::FT(last_len);
      double angle = 3.14159265359 - std::atan2(CGAL::sqrt(CGAL::to_double((a - b).squared_length().exact())),
                                                CGAL::sqrt(CGAL::to_double((a + b).squared_length().exact())));
      line.push_back({iter->source(), cur_len, angle < max_angle, 0.0, 0.0});
    }
  }
  line.push_back({*std::prev(poly.vertices_end(), 1), 0.0, false, 0.0, 0.0});

  double last_cut = 0.0;
  auto corner = line.begin();
  while (corner != line.end()) {
    if (corner->corner) break;
    last_cut += corner->length;
    corner++;
  }
  last_cut = std::min(last_cut, max_cut) * 0.5;
  if (corner == line.end()) {
    return line;
  }
  for (auto iter = corner; iter != line.end(); iter++) {
    if (iter->corner) {
      iter->cut_dist_back = last_cut;
      double dist_to_next = 0.0;
      auto forward_iter = iter;
      do {
        dist_to_next += forward_iter->length;
        forward_iter++;
      } while (forward_iter != line.end() && !forward_iter->corner);
      last_cut = std::min(dist_to_next * (forward_iter == line.end() ? 1.0 : 0.5), max_cut) * 0.5;
      iter->cut_dist_forward = last_cut;
    }
  }

  return line;
}

Segment_trait::Polygon_2 clip_corner_ring(const Segment_trait::Polygon_2& poly, double max_angle, double max_cut, size_t n_cut) {
  std::list<Edge> ring = create_edge_ring(poly, max_angle, max_cut);

  std::vector<Point_2> new_poly = cut_corners(ring, n_cut, true);

  return {new_poly.begin(), new_poly.end()};
}

Polygon clip_corner_ring(const Polygon& poly, double max_angle, double max_cut, size_t n_cut) {
  Polygon new_poly;
  if (!poly.is_unbounded()) {
    new_poly = Polygon(clip_corner_ring(poly.outer_boundary(), max_angle, max_cut, n_cut));
  }
  for (auto iter = poly.holes_begin(); iter != poly.holes_end(); iter++) {
    new_poly.add_hole(clip_corner_ring(*iter, max_angle, max_cut, n_cut));
  }

  return new_poly;
}

Polyline clip_corner_line(const Polyline& poly, double max_angle, double max_cut, size_t n_cut) {
  std::list<Edge> line = create_edge_line(poly, max_angle, max_cut);

  std::vector<Point_2> new_poly = cut_corners(line, n_cut, false);

  return {new_poly.begin(), new_poly.end()};
}

[[cpp11::register]]
SEXP poly_corner_cutting(SEXP poly, cpp11::doubles max_angle, cpp11::doubles max_cut, cpp11::integers n_cut) {
  switch (polyclid::get_geometry_type(poly)) {
    case polyclid::POLYGON: {
      std::vector<Polygon> polygons = polyclid::get_polygon_vec(poly);
      std::vector<Polygon> res;
      res.reserve(polygons.size());
      for (size_t i = 0; i < polygons.size(); ++i) {
        if (polygons[i].is_na()) {
          res.push_back(Polygon::NA_value());
        } else {
          res.push_back(clip_corner_ring(polygons[i], max_angle[i % max_angle.size()], max_cut[i % max_cut.size()], n_cut[i % n_cut.size()]));
        }
      }
      return polyclid::create_polygon_vec(res);
    }
    case polyclid::POLYLINE: {
      std::vector<Polyline> polylines = polyclid::get_polyline_vec(poly);
      std::vector<Polyline> res;
      res.reserve(polylines.size());
      for (size_t i = 0; i < polylines.size(); ++i) {
        if (polylines[i].is_na()) {
          res.push_back(Polyline::NA_value());
        } else {
          res.push_back(clip_corner_line(polylines[i], max_angle[i % max_angle.size()], max_cut[i % max_cut.size()], n_cut[i % n_cut.size()]));
        }
      }
      return polyclid::create_polyline_vec(res);
    }
    default: cpp11::stop("Geometry not supported");
  }
}
