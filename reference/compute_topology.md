# Compute spatial topology and weights for hexagonal smoothing

Finds spatial neighbours and computes Gaussian-based weights for
smoothing. Neighbours are computed up to \`neighbor_orders\` orders
(1st-order = touching, 2nd-order = neighbours of neighbours, etc.).
Weights decay with order using a Gaussian kernel.

## Usage

``` r
compute_topology(
  grid,
  projection_crs = NULL,
  neighbor_orders = 2,
  sigma = NULL,
  center_weight = 1,
  neighbor_weights_param = NULL,
  adaptive_sigma_factor = 0.5,
  sample_size = 100
)
```

## Arguments

- grid:

  sf object with polygonal geometries (or a list of sf objects).

- projection_crs:

  CRS used for distance calculations. Default \`NULL\` selects an
  appropriate UTM zone via \[\`get_utm_crs()\`\].

- neighbor_orders:

  Number of neighbour orders (positive integer).

- sigma:

  Gaussian bandwidth. \`NULL\` (default) auto-computes from
  \`avg_distance\` and \`adaptive_sigma_factor\`.

- center_weight:

  Weight for the centre cell.

- neighbor_weights_param:

  Optional list of length \`neighbor_orders\` with per-order weights
  (overrides Gaussian computation).

- adaptive_sigma_factor:

  Scaling factor for auto-sigma.

- sample_size:

  Cells to sample for the average-distance estimate.

## Value

List (or named list of lists if \`grid\` was a list) containing
\`neighbors\`, \`weights\`, \`avg_distance\`, \`sigma\`, \`grid_ids\`,
\`grid_indices\`, \`neighbor_orders\`.

## Weight computation

\- \`sigma = avg_distance \* adaptive_sigma_factor\` (auto-bandwidth)
when \`sigma\` is \`NULL\`. - For order N: weight ~ \`exp(-N^2 / (2 \*
sigma^2))\`. - All weights (including \`center_weight\`) are normalised
to sum to 1.
