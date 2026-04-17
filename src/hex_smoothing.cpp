// hexsmoothR: Hexagonal Grid Smoothing for Satellite Data
// Author: Max M. Lang (max.lang@stx.ox.ac.uk)
//
// C++ implementation of the spatial smoothing algorithms used by the package.
//
// Outputs per variable, per processed cell:
//   - raw                : the original (unsmoothed) centre-cell value
//   - weighted_combined  : weighted mean of centre + all-order neighbours
//   - neighbors_<N>{st|nd|rd|th} : mean of valid neighbours at order N
//
// The implementations are optimised for:
//   - Multi-variable processing (all variables done together per cell)
//   - Cache-friendly access via contiguous std::vector storage
//   - Single-pass weighted-sum accumulation
//   - Periodic R interrupt checks for long runs

#include <Rcpp.h>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <cstring>
using namespace Rcpp;

namespace {
inline std::string ordinal_suffix(int n) {
  int mod100 = n % 100;
  if (mod100 >= 11 && mod100 <= 13) return "th";
  switch (n % 10) {
    case 1: return "st";
    case 2: return "nd";
    case 3: return "rd";
    default: return "th";
  }
}

inline std::string fmt_invalid_hex(int hex_idx_one_based, int max_one_based) {
  return "Invalid hex_idx: " + std::to_string(hex_idx_one_based) +
    " (must be between 1 and " + std::to_string(max_one_based) + ")";
}
}  // namespace

// [[Rcpp::export]]
List process_district_all_vars_n_orders(
    List variable_values,      // List of NumericVectors (one per variable)
    List neighbors,            // List of length n_orders, each a list per cell
    List weights,              // List with center_weight + neighbor_weights
    IntegerVector hex_indices, // 1-based cell indices to process
    CharacterVector var_names
) {
  if (variable_values.size() == 0) stop("variable_values cannot be empty");
  if (hex_indices.size() == 0)     stop("hex_indices cannot be empty");
  if (neighbors.size() == 0)       stop("neighbors cannot be empty");

  double center_weight = as<double>(weights["center_weight"]);
  NumericVector neighbor_weights = weights["neighbor_weights"];
  const int n_orders = neighbor_weights.size();

  if (neighbors.size() != n_orders) {
    stop(
      "Number of neighbour orders (" + std::to_string(neighbors.size()) +
      ") must match neighbor_weights length (" + std::to_string(n_orders) + ")"
    );
  }
  if (var_names.size() != variable_values.size()) {
    stop("var_names size must match variable_values size");
  }

  const int n_vars    = variable_values.size();
  const int n_process = hex_indices.size();

  std::vector<std::vector<double>> var_data(n_vars);
  std::vector<std::vector<unsigned char>> is_valid(n_vars);
  int n_cells = -1;
  for (int v = 0; v < n_vars; v++) {
    NumericVector values = variable_values[v];
    if (n_cells < 0) {
      n_cells = values.size();
    } else if (values.size() != n_cells) {
      stop(
        "All variables must have the same length. Variable 0: " +
        std::to_string(n_cells) + ", variable " + std::to_string(v) +
        ": " + std::to_string(values.size())
      );
    }
    var_data[v].assign(values.begin(), values.end());
    is_valid[v].resize(values.size());
    for (int i = 0; i < values.size(); i++) {
      is_valid[v][i] = !NumericVector::is_na(values[i]) ? 1u : 0u;
    }
  }

  List first_order_neighbors = neighbors[0];
  const int n_neighbor_cells = first_order_neighbors.size();

  std::vector<std::vector<std::vector<double>>> order_means(
      n_vars,
      std::vector<std::vector<double>>(
          n_orders, std::vector<double>(n_process, NA_REAL)));
  std::vector<std::vector<double>> raw_values(
      n_vars, std::vector<double>(n_process, NA_REAL));
  std::vector<std::vector<double>> weighted_combined(
      n_vars, std::vector<double>(n_process, NA_REAL));

  std::vector<std::vector<int>> valid_neighbor_indices(n_orders);
  for (int order = 0; order < n_orders; order++) {
    valid_neighbor_indices[order].reserve(64);
  }

  for (int i = 0; i < n_process; i++) {
    if ((i & 4095) == 0) {
      R_CheckUserInterrupt();
      if (i > 0) Rprintf("Processing cell %d of %d\r", i, n_process);
    }

    int hex_idx = hex_indices[i] - 1;
    if (hex_idx < 0 || hex_idx >= n_neighbor_cells) {
      stop(fmt_invalid_hex(hex_indices[i], n_neighbor_cells));
    }

    for (int order = 0; order < n_orders; order++) {
      valid_neighbor_indices[order].clear();
      List order_neighbors = neighbors[order];
      IntegerVector neighs = order_neighbors[hex_idx];
      for (int j = 0; j < neighs.size(); j++) {
        int idx = neighs[j] - 1;
        if (idx >= 0 && idx < n_cells) {
          valid_neighbor_indices[order].push_back(idx);
        }
      }
    }

    for (int v = 0; v < n_vars; v++) {
      const std::vector<double>& values   = var_data[v];
      const std::vector<unsigned char>& valid = is_valid[v];

      raw_values[v][i] = values[hex_idx];

      for (int order = 0; order < n_orders; order++) {
        if (!valid_neighbor_indices[order].empty()) {
          double sum = 0.0;
          int count = 0;
          for (int idx : valid_neighbor_indices[order]) {
            if (valid[idx]) {
              sum += values[idx];
              count++;
            }
          }
          if (count > 0) order_means[v][order][i] = sum / count;
        }
      }

      double weighted_sum = 0.0;
      double weight_sum   = 0.0;

      if (valid[hex_idx]) {
        weighted_sum += values[hex_idx] * center_weight;
        weight_sum   += center_weight;
      }

      for (int order = 0; order < n_orders; order++) {
        const double w = neighbor_weights[order];
        for (int idx : valid_neighbor_indices[order]) {
          if (valid[idx]) {
            weighted_sum += values[idx] * w;
            weight_sum   += w;
          }
        }
      }

      if (weight_sum > 0) weighted_combined[v][i] = weighted_sum / weight_sum;
    }
  }

  if (n_process > 5000) Rprintf("\n");

  List all_results;
  for (int v = 0; v < n_vars; v++) {
    std::string var_name = as<std::string>(var_names[v]);
    List var_results;

    var_results["raw"] =
        NumericVector(raw_values[v].begin(), raw_values[v].end());
    var_results["weighted_combined"] =
        NumericVector(weighted_combined[v].begin(),
                      weighted_combined[v].end());

    for (int order = 0; order < n_orders; order++) {
      int order_num = order + 1;
      std::string key = "neighbors_" + std::to_string(order_num) +
                        ordinal_suffix(order_num);
      var_results[key] =
          NumericVector(order_means[v][order].begin(),
                        order_means[v][order].end());
    }

    all_results[var_name] = var_results;
  }

  return all_results;
}

// [[Rcpp::export]]
List process_district_all_vars(
    List variable_values,      // List of NumericVectors
    List first_neighbors,
    List second_neighbors,
    NumericVector weights,     // c(center, first_order, second_order)
    IntegerVector hex_indices,
    CharacterVector var_names
) {
  if (variable_values.size() == 0) stop("variable_values cannot be empty");
  if (hex_indices.size() == 0)     stop("hex_indices cannot be empty");
  if (weights.size() != 3) {
    stop("weights must have exactly 3 elements (center, first_order, second_order)");
  }
  if (var_names.size() != variable_values.size()) {
    stop("var_names size must match variable_values size");
  }

  const int n_vars    = variable_values.size();
  const int n_process = hex_indices.size();

  const double center_weight       = weights[0];
  const double first_order_weight  = weights[1];
  const double second_order_weight = weights[2];

  std::vector<std::vector<double>> var_data(n_vars);
  std::vector<std::vector<unsigned char>> is_valid(n_vars);
  int n_cells = -1;
  for (int v = 0; v < n_vars; v++) {
    NumericVector values = variable_values[v];
    if (n_cells < 0) {
      n_cells = values.size();
    } else if (values.size() != n_cells) {
      stop(
        "All variables must have the same length. Variable 0: " +
        std::to_string(n_cells) + ", variable " + std::to_string(v) +
        ": " + std::to_string(values.size())
      );
    }
    var_data[v].assign(values.begin(), values.end());
    is_valid[v].resize(values.size());
    for (int i = 0; i < values.size(); i++) {
      is_valid[v][i] = !NumericVector::is_na(values[i]) ? 1u : 0u;
    }
  }

  const int n_neighbor_cells = first_neighbors.size();

  std::vector<std::vector<double>> raw_values(
      n_vars, std::vector<double>(n_process, NA_REAL));
  std::vector<std::vector<double>> first_means(
      n_vars, std::vector<double>(n_process, NA_REAL));
  std::vector<std::vector<double>> second_means(
      n_vars, std::vector<double>(n_process, NA_REAL));
  std::vector<std::vector<double>> weighted_first(
      n_vars, std::vector<double>(n_process, NA_REAL));
  std::vector<std::vector<double>> weighted_combined(
      n_vars, std::vector<double>(n_process, NA_REAL));

  std::vector<int> valid_first_indices;
  std::vector<int> valid_second_indices;
  valid_first_indices.reserve(20);
  valid_second_indices.reserve(40);

  for (int i = 0; i < n_process; i++) {
    if ((i & 4095) == 0) {
      R_CheckUserInterrupt();
      if (i > 0) Rprintf("Processing cell %d of %d\r", i, n_process);
    }

    int hex_idx = hex_indices[i] - 1;
    if (hex_idx < 0 || hex_idx >= n_neighbor_cells) {
      stop(fmt_invalid_hex(hex_indices[i], n_neighbor_cells));
    }

    IntegerVector first_neighs  = first_neighbors[hex_idx];
    IntegerVector second_neighs = second_neighbors[hex_idx];

    valid_first_indices.clear();
    valid_second_indices.clear();

    for (int j = 0; j < first_neighs.size(); j++) {
      int idx = first_neighs[j] - 1;
      if (idx >= 0 && idx < n_cells) valid_first_indices.push_back(idx);
    }
    for (int j = 0; j < second_neighs.size(); j++) {
      int idx = second_neighs[j] - 1;
      if (idx >= 0 && idx < n_cells) valid_second_indices.push_back(idx);
    }

    for (int v = 0; v < n_vars; v++) {
      const std::vector<double>& values   = var_data[v];
      const std::vector<unsigned char>& valid = is_valid[v];

      raw_values[v][i] = values[hex_idx];

      if (!valid_first_indices.empty()) {
        double sum = 0.0;
        int count = 0;
        for (int idx : valid_first_indices) {
          if (valid[idx]) { sum += values[idx]; count++; }
        }
        if (count > 0) first_means[v][i] = sum / count;
      }

      if (!valid_second_indices.empty()) {
        double sum = 0.0;
        int count = 0;
        for (int idx : valid_second_indices) {
          if (valid[idx]) { sum += values[idx]; count++; }
        }
        if (count > 0) second_means[v][i] = sum / count;
      }

      double weighted_sum_first = 0.0;
      double weight_sum_first   = 0.0;
      double weighted_sum_all   = 0.0;
      double weight_sum_all     = 0.0;

      if (valid[hex_idx]) {
        double val = values[hex_idx];
        weighted_sum_first += val * center_weight;
        weight_sum_first   += center_weight;
        weighted_sum_all   += val * center_weight;
        weight_sum_all     += center_weight;
      }

      for (int idx : valid_first_indices) {
        if (valid[idx]) {
          double val = values[idx];
          weighted_sum_first += val * first_order_weight;
          weight_sum_first   += first_order_weight;
          weighted_sum_all   += val * first_order_weight;
          weight_sum_all     += first_order_weight;
        }
      }
      if (weight_sum_first > 0) {
        weighted_first[v][i] = weighted_sum_first / weight_sum_first;
      }

      for (int idx : valid_second_indices) {
        if (valid[idx]) {
          double val = values[idx];
          weighted_sum_all += val * second_order_weight;
          weight_sum_all   += second_order_weight;
        }
      }
      if (weight_sum_all > 0) {
        weighted_combined[v][i] = weighted_sum_all / weight_sum_all;
      }
    }
  }

  if (n_process > 5000) Rprintf("\n");

  List all_results;
  for (int v = 0; v < n_vars; v++) {
    std::string var_name = as<std::string>(var_names[v]);

    List var_results = List::create(
      Named("raw")               = NumericVector(raw_values[v].begin(),
                                                 raw_values[v].end()),
      Named("neighbors_1st")     = NumericVector(first_means[v].begin(),
                                                 first_means[v].end()),
      Named("neighbors_2nd")     = NumericVector(second_means[v].begin(),
                                                 second_means[v].end()),
      Named("smoothed_1st")      = NumericVector(weighted_first[v].begin(),
                                                 weighted_first[v].end()),
      Named("smoothed_all")      = NumericVector(weighted_combined[v].begin(),
                                                 weighted_combined[v].end()),
      Named("weighted_combined") = NumericVector(weighted_combined[v].begin(),
                                                 weighted_combined[v].end())
    );

    all_results[var_name] = var_results;
  }

  return all_results;
}

// [[Rcpp::export]]
List simple_smoothing_single_var(
    NumericVector values,
    List first_neighbors,
    List second_neighbors,
    IntegerVector hex_indices
) {
  List variable_values   = List::create(values);
  NumericVector weights  = NumericVector::create(1.0, 1.0, 1.0);
  CharacterVector var_names = CharacterVector::create("var1");

  List results = process_district_all_vars(
    variable_values, first_neighbors, second_neighbors,
    weights, hex_indices, var_names
  );
  List var_results = results["var1"];

  return List::create(
    Named("first_order_means")  = var_results["neighbors_1st"],
    Named("second_order_means") = var_results["neighbors_2nd"],
    Named("combined_means")     = var_results["weighted_combined"]
  );
}

// [[Rcpp::export]]
List gaussian_smoothing_single_var(
    NumericVector values,
    List first_neighbors,
    List second_neighbors,
    double center_weight,
    double first_order_weight,
    double second_order_weight,
    IntegerVector hex_indices
) {
  List variable_values  = List::create(values);
  NumericVector weights = NumericVector::create(
    center_weight, first_order_weight, second_order_weight
  );
  CharacterVector var_names = CharacterVector::create("var1");

  List results = process_district_all_vars(
    variable_values, first_neighbors, second_neighbors,
    weights, hex_indices, var_names
  );
  List var_results = results["var1"];

  return List::create(
    Named("first_order_weighted") = var_results["smoothed_1st"],
    Named("combined_weighted")    = var_results["smoothed_all"]
  );
}
