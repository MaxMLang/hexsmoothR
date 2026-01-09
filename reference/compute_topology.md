# Compute spatial topology and weights for hexagonal smoothing

Finds spatial neighbors and computes Gaussian-based weights for
smoothing. This function determines which grid cells are neighbors at
specified orders (e.g., first-order neighbors are touching, second-order
are neighbors of neighbors), then calculates appropriate weights based
on spatial distance using Gaussian decay.

## Usage

``` r
compute_topology(
  grid,
  projection_crs = 32636,
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

  sf object with polygonal geometries (hexagons or squares) or list of
  sf objects. Must contain columns \`grid_id\` and \`grid_index\`, or
  they will be created automatically.

- projection_crs:

  CRS for distance calculations (default UTM 36N). Should be a projected
  coordinate system (like UTM) for accurate distance measurements in
  meters.

- neighbor_orders:

  Number of neighbor orders to compute (default 2). Higher orders
  provide more extensive smoothing but increase computation time and
  memory usage.

- sigma:

  Gaussian bandwidth parameter in meters (NULL = auto-computed as
  avg_distance \* adaptive_sigma_factor). Controls how quickly weights
  decay with distance.

- center_weight:

  Weight for the center cell (default 1.0). Higher values give more
  influence to the original cell value.

- neighbor_weights_param:

  Optional list of weights for each neighbor order (overrides Gaussian
  computation). If specified, must have length equal to
  \`neighbor_orders\`. Example: \`list(0.5, 0.25)\` for 2 orders.

- adaptive_sigma_factor:

  Scaling factor for automatic sigma computation (default 0.5). Smaller
  values = less smoothing, larger values = more smoothing.

- sample_size:

  Number of cells to sample for distance calculations (default 100).
  Used to estimate average neighbor distance efficiently.

## Value

List containing:

- neighbors:

  List of lists, each containing neighbor indices for each order

- weights:

  List with center_weight and neighbor weights for each order (all
  normalized)

- avg_distance:

  Average distance between neighboring cells in meters

- sigma:

  Gaussian bandwidth parameter used

- grid_ids:

  Character vector of grid cell identifiers

- grid_indices:

  Integer vector of grid cell indices

- neighbor_orders:

  Number of neighbor orders computed

## Details

\*\*Weight Computation:\*\*

The function uses Gaussian weights based on spatial distance: - \`sigma
= avg_distance \* adaptive_sigma_factor\` (adaptive bandwidth) - Weights
decay exponentially with neighbor order: \`exp(-order^2 / (2 \*
sigma^2))\` - All weights are normalized so they sum to 1

\*\*How weights are applied in smoothing:\*\* - Each individual neighbor
gets the full weight (not averaged) - Center cell: \`value \*
center_weight\` - Each neighbor at order N: \`value \*
neighbor_weight_order_N\` - Final value = sum of all weighted values /
sum of all weights used

This means cells with more neighbors get more influence from their
neighborhood, which is realistic for spatial analysis.
