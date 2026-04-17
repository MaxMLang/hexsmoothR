# R fallback for N-order spatial smoothing

Pure-R implementation matching the C++ algorithm in
\`process_district_all_vars_n_orders()\`. Used automatically by
\[\`smooth_variables()\`\] when the compiled library is unavailable.

## Usage

``` r
smooth_variables_r_fallback(
  variable_values,
  neighbors,
  weights,
  hex_indices,
  var_names
)
```

## Arguments

- variable_values:

  List of numeric vectors (one per variable).

- neighbors:

  List of length \`n_orders\`; each element is a list of length
  \`n_cells\` of integer neighbour indices.

- weights:

  List with \`center_weight\` and \`neighbor_weights\`.

- hex_indices:

  Integer vector of cell indices (1-based) to process.

- var_names:

  Character vector of variable names.

## Value

Named list with one entry per variable, each containing \`raw\` (the
centre cell value), \`weighted_combined\` (the smoothed value), and
\`neighbors\_\<N\>\<suffix\>\` (per-order neighbour means, e.g.
\`neighbors_1st\`, \`neighbors_2nd\`, ...).
