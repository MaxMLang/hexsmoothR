# Convert between hexagon measurements

Helper functions to convert between different hexagon measurements:
flat-to-flat distance, edge length, and circumradius (center to vertex).

## Usage

``` r
hex_flat_to_edge(flat_to_flat)

hex_flat_to_circumradius(flat_to_flat)

hex_edge_to_flat(edge_length)

hex_circumradius_to_flat(circumradius)
```

## Arguments

- flat_to_flat:

  Flat-to-flat distance (distance between opposite edges)

- edge_length:

  Edge length of the hexagon

- circumradius:

  Circumradius of the hexagon

## Value

Numeric value in the same units as input

## Details

Convert between hexagon measurements

These functions convert between different hexagon measurements: -
\`hex_flat_to_edge\`: Convert flat-to-flat distance to edge length -
\`hex_flat_to_circumradius\`: Convert flat-to-flat distance to
circumradius - \`hex_edge_to_flat\`: Convert edge length to flat-to-flat
distance - \`hex_circumradius_to_flat\`: Convert circumradius to
flat-to-flat distance

## Examples

``` r
# Convert 1000m flat-to-flat to edge length
hex_flat_to_edge(1000)  # Returns ~577.35m
#> [1] 577.3503

# Convert 1000m flat-to-flat to circumradius  
hex_flat_to_circumradius(1000)  # Returns ~577.35m
#> [1] 577.3503

# Convert 577.35m edge length to flat-to-flat
hex_edge_to_flat(577.35)  # Returns ~1000m
#> [1] 999.9995

# Convert 577.35m circumradius to flat-to-flat
hex_circumradius_to_flat(577.35)  # Returns ~1000m
#> [1] 999.9995
```
