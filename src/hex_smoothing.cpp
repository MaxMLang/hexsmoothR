// hexsmoothR: Hexagonal Grid Smoothing for Satellite Data
// Author: Max M. Lang (max.lang@stx.ox.ac.uk)
// C++ implementation of spatial smoothing algorithms
//
// This file implements the core spatial smoothing algorithm for hexagonal grids.
// The algorithm computes multiple types of spatial smoothing:
// - First-order neighbor means (simple average of adjacent cells)
// - Second-order neighbor means (simple average of second-order neighbors)
// - Combined means (average of center + all neighbors)
// - Weighted first-order means (weighted average of center + first-order neighbors)
// - Weighted combined means (weighted average of center + all neighbors)
// - N-order neighbor smoothing (configurable number of neighbor orders)
//
// The implementation is optimized for:
// - Multi-variable processing (all variables processed simultaneously)
// - Memory efficiency (pre-allocated arrays, contiguous memory layout)
// - Performance (manual loop unrolling, single-pass calculations)
// - Large datasets (chunking support, progress reporting)

#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstring>
using namespace Rcpp;

// [[Rcpp::export]]
List process_district_all_vars_n_orders(
    List variable_values,      // List of NumericVectors
    List neighbors,            // List of neighbor lists for each order
    List weights,              // List with center_weight and neighbor_weights
    IntegerVector hex_indices,
    CharacterVector var_names
) {
  // Input validation
  if (variable_values.size() == 0) {
    stop("variable_values cannot be empty");
  }
  if (hex_indices.size() == 0) {
    stop("hex_indices cannot be empty");
  }
  if (neighbors.size() == 0) {
    stop("neighbors cannot be empty");
  }
  
  // Extract weights
  double center_weight = weights["center_weight"];
  NumericVector neighbor_weights = weights["neighbor_weights"];
  int n_orders = neighbor_weights.size();
  
  if (neighbors.size() != n_orders) {
    stop("Number of neighbor orders (%d) must match neighbor_weights length (%d)", 
         neighbors.size(), n_orders);
  }
  
  int n_vars = variable_values.size();
  int n_process = hex_indices.size();
  
  // Pre-allocate all output arrays
  std::vector<std::vector<std::vector<double>>> order_means(n_vars, 
    std::vector<std::vector<double>>(n_orders, std::vector<double>(n_process, NA_REAL)));
  std::vector<std::vector<double>> weighted_combined(n_vars, std::vector<double>(n_process, NA_REAL));
  
  // Extract all variable data into contiguous arrays for better cache performance
  std::vector<std::vector<double>> var_data(n_vars);
  for (int v = 0; v < n_vars; v++) {
    NumericVector values = variable_values[v];
    var_data[v].resize(values.size());
    std::memcpy(var_data[v].data(), values.begin(), values.size() * sizeof(double));
  }
  
  // Pre-compute NA masks for all variables (faster than repeated NA checks)
  std::vector<std::vector<bool>> is_valid(n_vars);
  for (int v = 0; v < n_vars; v++) {
    NumericVector values = variable_values[v];
    is_valid[v].resize(values.size());
    for (int i = 0; i < values.size(); i++) {
      is_valid[v][i] = !NumericVector::is_na(values[i]);
    }
  }
  
  // Pre-allocate workspace vectors to avoid repeated allocations
  std::vector<std::vector<int>> valid_neighbor_indices(n_orders);
  for (int order = 0; order < n_orders; order++) {
    valid_neighbor_indices[order].reserve(50);  // Most hexagons have reasonable neighbor counts
  }
  
  // Process each hexagon once, computing all variables together
  for (int i = 0; i < n_process; i++) {
    // Check for R interrupt periodically
    if (i % 5000 == 0) {
      R_CheckUserInterrupt();
      if (i > 0) {
        Rprintf("Processing hexagon %d of %d\r", i, n_process);
      }
    }
    
    int hex_idx = hex_indices[i] - 1;  // Convert to 0-based indexing
    
    // Validate hex_idx bounds
    List first_order_neighbors = neighbors[0];
    if (hex_idx < 0 || hex_idx >= first_order_neighbors.size()) {
      stop("Invalid hex_idx: %d (must be between 1 and %d)", hex_indices[i], first_order_neighbors.size());
    }
    
    // Clear and reuse vectors
    for (int order = 0; order < n_orders; order++) {
      valid_neighbor_indices[order].clear();
    }
    
    // Pre-compute valid neighbor indices for all orders
    int max_idx = var_data[0].size();
    for (int v = 1; v < n_vars; v++) {
      if (var_data[v].size() != max_idx) {
        stop("All variables must have the same size. Variable 0: %d, Variable %d: %d", 
             max_idx, v, var_data[v].size());
      }
    }
    
    // Process neighbors for each order
    for (int order = 0; order < n_orders; order++) {
      List order_neighbors = neighbors[order];
      IntegerVector neighs = order_neighbors[hex_idx];
      
      for (int j = 0; j < neighs.size(); j++) {
        int idx = neighs[j] - 1;
        if (idx >= 0 && idx < max_idx) {
          valid_neighbor_indices[order].push_back(idx);
        }
      }
    }
    
    // Process all variables for this hexagon
    for (int v = 0; v < n_vars; v++) {
      const std::vector<double>& values = var_data[v];
      const std::vector<bool>& valid = is_valid[v];
      
      // Calculate means for each order
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
          
          if (count > 0) {
            order_means[v][order][i] = sum / count;
          }
        }
      }
      
      // Calculate weighted combined result
      double weighted_sum = 0.0;
      double weight_sum = 0.0;
      
      // Add center cell
      if (valid[hex_idx]) {
        double val = values[hex_idx];
        weighted_sum += val * center_weight;
        weight_sum += center_weight;
      }
      
      // Add neighbors from all orders
      for (int order = 0; order < n_orders; order++) {
        for (int idx : valid_neighbor_indices[order]) {
          if (valid[idx]) {
            double val = values[idx];
            weighted_sum += val * neighbor_weights[order];
            weight_sum += neighbor_weights[order];
          }
        }
      }
      
      // Store weighted combined result
      if (weight_sum > 0) {
        weighted_combined[v][i] = weighted_sum / weight_sum;
      }
    }
  }
  
  // Clear progress message
  if (n_process > 5000) {
    Rprintf("\n");
  }
  
  // Convert results back to R format
  List all_results;
  
  for (int v = 0; v < n_vars; v++) {
    std::string var_name = as<std::string>(var_names[v]);
    
    // Create result structure
    List var_results;
    var_results["raw"] = NumericVector(weighted_combined[v].begin(), weighted_combined[v].end());
    var_results["weighted_combined"] = NumericVector(weighted_combined[v].begin(), weighted_combined[v].end());
    
    // Add order-specific results
    for (int order = 0; order < n_orders; order++) {
      int order_num = order + 1;
      std::string suffix = (order_num == 1) ? "st" : (order_num == 2) ? "nd" : (order_num == 3) ? "rd" : "th";
      std::string order_name = "neighbors_" + std::to_string(order_num) + suffix;
      var_results[order_name] = NumericVector(order_means[v][order].begin(), order_means[v][order].end());
    }
    
    // Maintain backward compatibility for 2-order system
    if (n_orders == 2) {
      var_results["neighbors_1st"] = NumericVector(order_means[v][0].begin(), order_means[v][0].end());
      var_results["neighbors_2nd"] = NumericVector(order_means[v][1].begin(), order_means[v][1].end());
      var_results["smoothed_1st"] = NumericVector(weighted_combined[v].begin(), weighted_combined[v].end());
      var_results["smoothed_all"] = NumericVector(weighted_combined[v].begin(), weighted_combined[v].end());
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
  // Input validation
  if (variable_values.size() == 0) {
    stop("variable_values cannot be empty");
  }
  if (hex_indices.size() == 0) {
    stop("hex_indices cannot be empty");
  }
  if (weights.size() != 3) {
    stop("weights must have exactly 3 elements (center, first_order, second_order)");
  }
  if (var_names.size() != variable_values.size()) {
    stop("var_names size must match variable_values size");
  }
  
  int n_vars = variable_values.size();
  int n_process = hex_indices.size();
  
  double center_weight = weights[0];
  double first_order_weight = weights[1];
  double second_order_weight = weights[2];
  
  // Pre-allocate all output arrays
  std::vector<std::vector<double>> first_means(n_vars, std::vector<double>(n_process, NA_REAL));
  std::vector<std::vector<double>> second_means(n_vars, std::vector<double>(n_process, NA_REAL));
  std::vector<std::vector<double>> combined_means(n_vars, std::vector<double>(n_process, NA_REAL));
  std::vector<std::vector<double>> weighted_first(n_vars, std::vector<double>(n_process, NA_REAL));
  std::vector<std::vector<double>> weighted_combined(n_vars, std::vector<double>(n_process, NA_REAL));
  
  // Extract all variable data into contiguous arrays for better cache performance
  std::vector<std::vector<double>> var_data(n_vars);
  for (int v = 0; v < n_vars; v++) {
    NumericVector values = variable_values[v];
    var_data[v].resize(values.size());
    std::memcpy(var_data[v].data(), values.begin(), values.size() * sizeof(double));
  }
  
  // Pre-compute NA masks for all variables (faster than repeated NA checks)
  std::vector<std::vector<bool>> is_valid(n_vars);
  for (int v = 0; v < n_vars; v++) {
    NumericVector values = variable_values[v];
    is_valid[v].resize(values.size());
    for (int i = 0; i < values.size(); i++) {
      is_valid[v][i] = !NumericVector::is_na(values[i]);
    }
  }
  
  // Pre-allocate workspace vectors to avoid repeated allocations
  std::vector<int> valid_first_indices;
  std::vector<int> valid_second_indices;
  valid_first_indices.reserve(20);   // Most hexagons have ~6 first neighbors
  valid_second_indices.reserve(30);  // Most hexagons have ~12 second neighbors
  
  // Process each hexagon once, computing all variables together
  for (int i = 0; i < n_process; i++) {
    // Check for R interrupt periodically
    if (i % 5000 == 0) {
      R_CheckUserInterrupt();
      if (i > 0) {
        Rprintf("Processing hexagon %d of %d\r", i, n_process);
      }
    }
    
    int hex_idx = hex_indices[i] - 1;  // Convert to 0-based indexing
    
    // Validate hex_idx bounds
    if (hex_idx < 0 || hex_idx >= first_neighbors.size()) {
      stop("Invalid hex_idx: %d (must be between 1 and %d)", hex_indices[i], first_neighbors.size());
    }
    
    // Get neighbors once for all variables
    IntegerVector first_neighs = first_neighbors[hex_idx];
    IntegerVector second_neighs = second_neighbors[hex_idx];
    
    // Clear and reuse vectors
    valid_first_indices.clear();
    valid_second_indices.clear();
    
    // Pre-compute valid neighbor indices to avoid repeated bounds checking
    // Validate that all variables have the same size
    int max_idx = var_data[0].size();
    for (int v = 1; v < n_vars; v++) {
      if (var_data[v].size() != max_idx) {
        stop("All variables must have the same size. Variable 0: %d, Variable %d: %d", 
             max_idx, v, var_data[v].size());
      }
    }
    
    for (int j = 0; j < first_neighs.size(); j++) {
      int idx = first_neighs[j] - 1;
      if (idx >= 0 && idx < max_idx) {
        valid_first_indices.push_back(idx);
      }
    }
    
    for (int j = 0; j < second_neighs.size(); j++) {
      int idx = second_neighs[j] - 1;
      if (idx >= 0 && idx < max_idx) {
        valid_second_indices.push_back(idx);
      }
    }
    
    // Process all variables for this hexagon
    for (int v = 0; v < n_vars; v++) {
      const std::vector<double>& values = var_data[v];
      const std::vector<bool>& valid = is_valid[v];
      
      // Calculate first-order mean
      if (!valid_first_indices.empty()) {
        double sum = 0.0;
        int count = 0;
        
        // Manual unrolling for better performance
        int n_first = valid_first_indices.size();
        int j = 0;
        
        // Process 4 at a time when possible
        for (; j + 3 < n_first; j += 4) {
          int idx0 = valid_first_indices[j];
          int idx1 = valid_first_indices[j+1];
          int idx2 = valid_first_indices[j+2];
          int idx3 = valid_first_indices[j+3];
          
          if (valid[idx0]) { sum += values[idx0]; count++; }
          if (valid[idx1]) { sum += values[idx1]; count++; }
          if (valid[idx2]) { sum += values[idx2]; count++; }
          if (valid[idx3]) { sum += values[idx3]; count++; }
        }
        
        // Handle remaining elements
        for (; j < n_first; j++) {
          int idx = valid_first_indices[j];
          if (valid[idx]) {
            sum += values[idx];
            count++;
          }
        }
        
        if (count > 0) {
          first_means[v][i] = sum / count;
        }
      }
      
      // Calculate second-order mean
      if (!valid_second_indices.empty()) {
        double sum = 0.0;
        int count = 0;
        
        for (int idx : valid_second_indices) {
          if (valid[idx]) {
            sum += values[idx];
            count++;
          }
        }
        
        if (count > 0) {
          second_means[v][i] = sum / count;
        }
      }
      
      // Calculate combined mean and weighted means in single pass
      double combined_sum = 0.0;
      int combined_count = 0;
      double weighted_sum_first = 0.0;
      double weight_sum_first = 0.0;
      double weighted_sum_all = 0.0;
      double weight_sum_all = 0.0;
      
      // Add self value
      if (valid[hex_idx]) {
        double val = values[hex_idx];
        combined_sum += val;
        combined_count++;
        weighted_sum_first += val * center_weight;
        weight_sum_first += center_weight;
        weighted_sum_all += val * center_weight;
        weight_sum_all += center_weight;
      }
      
      // Add first-order neighbors (single pass for all calculations)
      for (int idx : valid_first_indices) {
        if (valid[idx]) {
          double val = values[idx];
          combined_sum += val;
          combined_count++;
          weighted_sum_first += val * first_order_weight;
          weight_sum_first += first_order_weight;
          weighted_sum_all += val * first_order_weight;
          weight_sum_all += first_order_weight;
        }
      }
      
      // Store first-order weighted result
      if (weight_sum_first > 0) {
        weighted_first[v][i] = weighted_sum_first / weight_sum_first;
      }
      
      // Add second-order neighbors
      for (int idx : valid_second_indices) {
        if (valid[idx]) {
          double val = values[idx];
          combined_sum += val;
          combined_count++;
          weighted_sum_all += val * second_order_weight;
          weight_sum_all += second_order_weight;
        }
      }
      
      // Store combined results
      if (combined_count > 0) {
        combined_means[v][i] = combined_sum / combined_count;
      }
      
      if (weight_sum_all > 0) {
        weighted_combined[v][i] = weighted_sum_all / weight_sum_all;
      }
    }
  }
  
  // Clear progress message
  if (n_process > 5000) {
    Rprintf("\n");
  }
  
  // Convert results back to R format
  List all_results;
  
  for (int v = 0; v < n_vars; v++) {
    std::string var_name = as<std::string>(var_names[v]);
    
    // Convert std::vectors to NumericVectors
    NumericVector first_mean_r(first_means[v].begin(), first_means[v].end());
    NumericVector second_mean_r(second_means[v].begin(), second_means[v].end());
    NumericVector combined_mean_r(combined_means[v].begin(), combined_means[v].end());
    NumericVector weighted_first_r(weighted_first[v].begin(), weighted_first[v].end());
    NumericVector weighted_combined_r(weighted_combined[v].begin(), weighted_combined[v].end());
    
    List var_results = List::create(
      Named("raw") = combined_mean_r,
      Named("neighbors_1st") = first_mean_r,
      Named("neighbors_2nd") = second_mean_r,
      Named("smoothed_1st") = weighted_first_r,
      Named("smoothed_all") = weighted_combined_r
    );
    
    all_results[var_name] = var_results;
  }
  
  return all_results;
}

// Keep the original single-variable functions for compatibility
// [[Rcpp::export]]
List simple_smoothing_single_var(
    NumericVector values,
    List first_neighbors,
    List second_neighbors,
    IntegerVector hex_indices
) {
  // Create dummy inputs for the optimized function
  List variable_values = List::create(values);
  NumericVector weights = NumericVector::create(1.0, 1.0, 1.0);
  CharacterVector var_names = CharacterVector::create("var1");
  
  // Call the correct function name
  List results = process_district_all_vars(
    variable_values, first_neighbors, second_neighbors,
    weights, hex_indices, var_names
  );
  
  List var_results = results["var1"];
  
  return List::create(
    Named("first_order_means") = var_results["neighbors_1st"],
                                            Named("second_order_means") = var_results["neighbors_2nd"],
                                                                                     Named("combined_means") = var_results["raw"]
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
  // Create inputs for the optimized function
  List variable_values = List::create(values);
  NumericVector weights = NumericVector::create(center_weight, first_order_weight, second_order_weight);
  CharacterVector var_names = CharacterVector::create("var1");
  
  // Call the correct function name
  List results = process_district_all_vars(
    variable_values, first_neighbors, second_neighbors,
    weights, hex_indices, var_names
  );
  
  List var_results = results["var1"];
  
  return List::create(
    Named("first_order_weighted") = var_results["smoothed_1st"],
                                               Named("combined_weighted") = var_results["smoothed_all"]
  );
}
