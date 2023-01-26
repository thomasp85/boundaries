#include <vector>
#include <utility>
#include <algorithm>

#include <cpp11/declarations.hpp>
#include <cpp11/integers.hpp>
#include <cpp11/strings.hpp>

#include <euclid.h>
#include <polyclid.h>

#include <CGAL/Polyline_simplification_2/simplify.h>

namespace PS = CGAL::Polyline_simplification_2;

#define HYBRID 0;
#define SCALED 1;
#define SQUARED 2;

#define ABOVE_COST 0;
#define BELOW_RATIO 1;
#define BELOW_COUNT 2;

template<typename T, typename C, typename S>
T poly_simplify_cost_stop_impl(T& poly, C cost, S stop);

template<typename C, typename S>
Polygon poly_simplify_cost_stop_impl(Polygon& poly, C cost, S stop) {
  return PS::simplify(poly, cost, stop);
}

template<typename C, typename S>
Polyline poly_simplify_cost_stop_impl(Polyline& poly, C cost, S stop) {
  std::vector<Point_2> res;
  PS::simplify(poly.begin(), poly.end(), cost, stop, std::back_inserter(res), false);
  return {res.begin(), res.end()};
}

template<typename T, typename C>
T poly_simplify_cost_impl(T& poly, C cost, int stop, double stop_par) {
  switch (stop) {
    case 0: {
      return poly_simplify_cost_stop_impl(poly, cost, PS::Stop_above_cost_threshold(stop_par));
    }
    case 1: {
      return poly_simplify_cost_stop_impl(poly, cost, PS::Stop_below_count_ratio_threshold(stop_par));
    }
    case 2: {
      return poly_simplify_cost_stop_impl(poly, cost, PS::Stop_below_count_threshold(stop_par));
    }
  }
  return poly;
}

template<typename T>
T poly_simplify_impl(T& poly, int cost, int stop, double cost_par, double stop_par) {
  switch (cost) {
    case 0: {
      return poly_simplify_cost_impl(poly, PS::Hybrid_squared_distance_cost<Kernel::FT>(cost_par), stop, stop_par);
    }
    case 1: {
      return poly_simplify_cost_impl(poly, PS::Scaled_squared_distance_cost(), stop, stop_par);
    }
    case 2: {
      return poly_simplify_cost_impl(poly, PS::Squared_distance_cost(), stop, stop_par);
    }
  }
  return poly;
}



[[cpp11::register]]
SEXP poly_simplify(SEXP poly, int cost, int stop, double cost_par, double stop_par) {
  switch(polyclid::get_geometry_type(poly)) {
    case polyclid::POLYGON: {
      std::vector<Polygon> p = polyclid::get_polygon_vec(poly);
      std::vector<Polygon> res;
      res.reserve(p.size());
      for (size_t i = 0; i < p.size(); ++i) {
        if (p[i].is_na()) {
          res.push_back(Polygon::NA_value());
        } else {
          res.push_back(poly_simplify_impl(p[i], cost, stop, cost_par, stop_par));
        }
      }
      return polyclid::create_polygon_vec(res);
    }
    case polyclid::POLYLINE: {
      std::vector<Polyline> p = polyclid::get_polyline_vec(poly);
      std::vector<Polyline> res;
      res.reserve(p.size());
      for (size_t i = 0; i < p.size(); ++i) {
        if (p[i].is_na()) {
          res.push_back(Polyline::NA_value());
        } else {
          res.push_back(poly_simplify_impl(p[i], cost, stop, cost_par, stop_par));
        }
      }
      return polyclid::create_polyline_vec(res);

    }
    default: {
      cpp11::stop("Don't know how to simplify the provided geometry");
    }
  }
  return R_NilValue;
}
