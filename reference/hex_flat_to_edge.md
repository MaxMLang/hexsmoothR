# Convert between hexagon measurements

Helper functions to convert between different hexagon measurements. For
a regular hexagon, the circumradius equals the edge length, so the
"edge" and "circumradius" helpers are mathematically identical (provided
for clarity at the call site).

## Usage

``` r
hex_flat_to_edge(flat_to_flat)

hex_flat_to_circumradius(flat_to_flat)

hex_edge_to_flat(edge_length)

hex_circumradius_to_flat(circumradius)
```

## Arguments

- flat_to_flat:

  Flat-to-flat distance (between opposite edges).

- edge_length:

  Edge length of the hexagon.

- circumradius:

  Circumradius (centre to vertex).

## Value

Numeric value in the same units as the input.

## Details

\- \`hex_flat_to_edge()\`: flat-to-flat distance to edge length -
\`hex_flat_to_circumradius()\`: flat-to-flat distance to circumradius -
\`hex_edge_to_flat()\`: edge length to flat-to-flat distance -
\`hex_circumradius_to_flat()\`: circumradius to flat-to-flat distance

## Examples

``` r
hex_flat_to_edge(1000)          # ~577.35
#> [1] 577.3503
hex_flat_to_circumradius(1000)  # ~577.35
#> [1] 577.3503
hex_edge_to_flat(577.35)        # ~1000
#> [1] 999.9995
hex_circumradius_to_flat(577.35)# ~1000
#> [1] 999.9995
```
