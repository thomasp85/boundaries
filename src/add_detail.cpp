#include <vector>
#include <utility>

#include <cpp11/declarations.hpp>
#include <cpp11/integers.hpp>

#include <euclid.h>
#include <polyclid.h>

template<typename Iter>
std::vector<double> get_segment_lengths(Iter begin, Iter end) {
  std::vector<double> res;
  for (auto iter = begin; iter != end; iter++) {
    res.push_back(CGAL::sqrt(CGAL::to_double(iter->squared_length().exact())));
  }
  return res;
}

template<typename Iter>
std::vector<Point_2> insert_vert_to_c(Iter begin, Iter end, unsigned int n,
                                      std::vector<double>& lengths,
                                      double total, bool wrap) {
  std::vector<Point_2> res;
  unsigned int cur_n = 0;
  for (auto iter = begin; iter != end; iter++) {
    cur_n++;
  }
  if (!wrap) cur_n++;
  if (cur_n >= n) {
    for (auto iter = begin; iter != end; iter++) {
      res.push_back(iter->source());
    }
    if (!wrap) res.push_back(std::prev(end, 1)->target());
    return res;
  }
  n -= cur_n;

  std::vector<unsigned int> splits;
  splits.reserve(cur_n);
  unsigned int n_remain = n;
  double max_seg = 0.0;
  size_t max_i = 0;
  for (size_t i = 0; i < lengths.size(); ++i) {
    unsigned int split = std::round(n * (lengths[i]) / total);
    split = std::min(split, n_remain);
    splits.push_back(split);
    lengths[i] /= split + 1;
    if (lengths[i] > max_seg) {
      max_seg = lengths[i];
      max_i = i;
    }
    n_remain -= split;
    if (n_remain <= 0) break;
  }
  splits[max_i] += n_remain;
  size_t i = 0;
  for (auto iter = begin; iter != end; iter++) {
    res.push_back(iter->source());
    unsigned int n_splits = splits[i];
    i++;
    if (n_splits == 0) continue;
    Vector_2 vec(*iter);
    vec /= n_splits + 1;
    while (n_splits--) {
      res.push_back(res.back() + vec);
    }
  }
  if (!wrap) res.push_back(std::prev(end, 1)->target());
  return res;
}

template<typename Iter>
std::vector<Point_2> insert_vert_to_l(Iter begin, Iter end, double max_l, bool wrap) {
  std::vector<Point_2> res;
  for (auto iter = begin; iter != end; iter++) {
    res.push_back(iter->source());
    Vector_2 vec(*iter);
    double l_seg = CGAL::sqrt(CGAL::to_double(vec.squared_length().exact()));
    unsigned int n = std::floor(l_seg / max_l);
    vec /= n + 1;
    while (n--) {
      res.push_back(res.back() + vec);
    }
  }
  if (!wrap) {
    res.push_back(std::prev(end, 1)->target());
  }
  return res;
}

template<typename Iter>
std::vector<Point_2> insert_vert_to_n(Iter begin, Iter end, unsigned int n, bool wrap) {
  std::vector<Point_2> res;
  for (auto iter = begin; iter != end; iter++) {
    res.push_back(iter->source());
    Vector_2 vec(*iter);
    vec /= n + 1;
    unsigned int n_temp = n;
    while (n_temp--) {
      res.push_back(res.back() + vec);
    }
  }
  if (!wrap) {
    res.push_back(std::prev(end, 1)->target());
  }
  return res;
}

[[cpp11::register]]
SEXP poly_add_detail(SEXP poly, cpp11::integers c, cpp11::doubles l, cpp11::integers n, bool use_l, bool use_n) {
  switch (polyclid::get_geometry_type(poly)) {
    case polyclid::POLYLINE: {
      std::vector<Polyline> lines = polyclid::get_polyline_vec(poly);
      std::vector<Polyline> res;
      res.reserve(lines.size());
      for (size_t i = 0; i < lines.size(); ++i) {
        if (lines[i].is_na() || lines[i].is_empty()) {
          res.push_back(lines[i]);
          continue;
        }
        std::vector<Point_2> line;
        if (use_l) {
          line = insert_vert_to_l(lines[i].edges_begin(), lines[i].edges_end(), l[i%l.size()], false);
        } else if (use_n) {
          line = insert_vert_to_n(lines[i].edges_begin(), lines[i].edges_end(), n[i%n.size()], false);
        } else {
          std::vector<double> lengths = get_segment_lengths(lines[i].edges_begin(), lines[i].edges_end());
          double full_length = std::accumulate(lengths.begin(), lengths.end(), 0.0);
          line = insert_vert_to_c(lines[i].edges_begin(), lines[i].edges_end(),
                                  c[i%c.size()], lengths, full_length, false);
        }
        res.emplace_back(line.begin(), line.end());
      }
      return polyclid::create_polyline_vec(res);
    }
    case polyclid::POLYGON: {
      std::vector<Polygon> polygons = polyclid::get_polygon_vec(poly);
      std::vector<Polygon> res;
      res.reserve(polygons.size());
      for (size_t i = 0; i < polygons.size(); ++i) {
        if (polygons[i].is_na() || (polygons[i].is_unbounded() && polygons[i].number_of_holes() == 0)) {
          res.push_back(polygons[i]);
          continue;
        }
        std::vector<Point_2> line;
        Polygon p;
        if (use_l) {
          if (!polygons[i].is_unbounded()) {
            line = insert_vert_to_l(polygons[i].outer_boundary().edges_begin(), polygons[i].outer_boundary().edges_end(), l[i%l.size()], true);
            p = Polygon(Segment_trait::Polygon_2(line.begin(), line.end()));
          }
          for (auto h_iter = polygons[i].holes_begin(); h_iter != polygons[i].holes_end(); h_iter++) {
            line = insert_vert_to_l(h_iter->edges_begin(), h_iter->edges_end(), l[i%l.size()], true);
            p.add_hole(Segment_trait::Polygon_2(line.begin(), line.end()));
          }
        } else if (use_n) {
          if (!polygons[i].is_unbounded()) {
            line = insert_vert_to_n(polygons[i].outer_boundary().edges_begin(), polygons[i].outer_boundary().edges_end(), n[i%n.size()], true);
            p = Polygon(Segment_trait::Polygon_2(line.begin(), line.end()));
          }
          for (auto h_iter = polygons[i].holes_begin(); h_iter != polygons[i].holes_end(); h_iter++) {
            line = insert_vert_to_n(h_iter->edges_begin(), h_iter->edges_end(), n[i%n.size()], true);
            p.add_hole(Segment_trait::Polygon_2(line.begin(), line.end()));
          }
        } else {
          std::vector< std::vector<double> > lengths;
          std::vector<double> full_lengths;
          if (!polygons[i].is_unbounded()) {
            lengths.push_back(get_segment_lengths(polygons[i].outer_boundary().edges_begin(), polygons[i].outer_boundary().edges_end()));
            full_lengths.push_back(std::accumulate(lengths.back().begin(), lengths.back().end(), 0.0));
          }
          for (auto h_iter = polygons[i].holes_begin(); h_iter != polygons[i].holes_end(); h_iter++) {
            lengths.push_back(get_segment_lengths(h_iter->edges_begin(), h_iter->edges_end()));
            full_lengths.push_back(std::accumulate(lengths.back().begin(), lengths.back().end(), 0.0));
          }
          std::vector<double> mean_length;
          for (size_t j = 0; j < full_lengths.size(); ++j) {
            mean_length.push_back(full_lengths[j] / lengths[j].size());
          }
          double summed_mean = std::accumulate(mean_length.begin(), mean_length.end(), 0.0);
          std::vector<unsigned> ns;
          for (size_t j = 0; j < mean_length.size() - 1; ++j) {
            ns.push_back(std::round(c[i%c.size()] * mean_length[j] / summed_mean));
          }
          ns.push_back(c[i%c.size()] - std::accumulate(ns.begin(), ns.end(), 0));
          size_t h = 0;
          if (!polygons[i].is_unbounded()) {
            line = insert_vert_to_c(polygons[i].outer_boundary().edges_begin(), polygons[i].outer_boundary().edges_end(),
                                    ns[h], lengths[h], full_lengths[h], true);
            p = Polygon(Segment_trait::Polygon_2(line.begin(), line.end()));
            h++;
          }
          for (auto h_iter = polygons[i].holes_begin(); h_iter != polygons[i].holes_end(); h_iter++) {
            line = insert_vert_to_c(h_iter->edges_begin(), h_iter->edges_end(),
                                    ns[h], lengths[h], full_lengths[h], true);
            p.add_hole(Segment_trait::Polygon_2(line.begin(), line.end()));
            h++;
          }
        }
        res.push_back(p);
        return polyclid::create_polygon_vec(res);
      }
    }
    default: cpp11::stop("Unknown geometry type");
  }
  return R_NilValue;
}

