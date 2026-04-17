# Apply spatial smoothing to variables on a hexagonal grid

Smooths variables using the topology produced by
\[\`compute_topology()\`\]. Uses the compiled C++ implementation when
available; if the C++ call fails for an unexpected reason, falls back to
the pure-R implementation and re-raises clear validation errors.

## Usage

``` r
smooth_variables(
  variable_values,
  neighbors,
  weights,
  hex_indices = NULL,
  var_names = NULL
)
```

## Arguments

- variable_values:

  Named list of numeric vectors (one per variable).

- neighbors:

  List of neighbour lists (per order). Either the \`neighbors\` element
  of a topology, or the topology object itself.

- weights:

  List with \`center_weight\` and \`neighbor_weights\`. Either the
  \`weights\` element of a topology, or the topology object itself.

- hex_indices:

  Integer vector of cell indices to process. Defaults to all cells.

- var_names:

  Character vector of variable names. Defaults to
  \`names(variable_values)\`.

## Value

Named list with one entry per variable. Each entry contains:

- \`raw\`: the original (unsmoothed) centre-cell values

- \`weighted_combined\`: weighted average of centre + all neighbours

- \`neighbors\_\<N\>\<suffix\>\`: mean of neighbours at order N (e.g.
  \`neighbors_1st\`, \`neighbors_2nd\`, \`neighbors_3rd\`)
