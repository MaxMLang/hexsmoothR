# hexsmoothR: Hexagonal Grid Smoothing for Satellite Data

<div align="center">
  <img src="man/figures/hexsmoothR_sticker.png" alt="hexsmoothR hex sticker" width="200">
</div>

<!-- badges: start -->
[![R-CMD-check](https://github.com/MaxMLang/hexsmoothR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/MaxMLang/hexsmoothR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

[![R](https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white)](https://www.r-project.org/)
[![Rcpp](https://img.shields.io/badge/Rcpp-276DC3?style=for-the-badge&logo=r&logoColor=white)](https://cran.r-project.org/package=Rcpp)
[![C++](https://img.shields.io/badge/C++-00599C?style=for-the-badge&logo=c%2B%2B&logoColor=white)](https://isocpp.org/)
[![GIS](https://img.shields.io/badge/GIS-4A90E2?style=for-the-badge&logo=arcgis&logoColor=white)](https://en.wikipedia.org/wiki/Geographic_information_system)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=for-the-badge)](https://opensource.org/licenses/MIT)


A comprehensive R package for creating hexagonal grids and applying spatial smoothing to satellite raster data. The package is specifically designed for hexagonal grid analysis, providing tools for extracting environmental variables from TIF files and applying Gaussian-weighted spatial smoothing using C++ optimization with chunking support for large datasets.

## Features

- **Hexagonal Grid Creation**: Create hexagonal grids from study areas (primary focus)
- **Spatial Topology**: Compute neighbor relationships and distances for hexagonal grids
- **Raster Processing**: Efficient extraction of values from TIF files
- **C++ Optimization**: High-performance spatial smoothing algorithms

- **Custom Weights**: Flexible Gaussian weight parameters
- **Multiple Formats**: Support for various raster formats and projections
- **Flexible Grid Input**: Accepts any sf POLYGON grid with proper structure

## Installation

### Prerequisites

**Windows users:** You need Rtools to compile C++ code. Install from [https://cran.r-project.org/bin/windows/Rtools/](https://cran.r-project.org/bin/windows/Rtools/)

**macOS users:** Install Xcode Command Line Tools: `xcode-select --install`

**Linux users:** Install build-essential: `sudo apt-get install build-essential` (Ubuntu/Debian) or equivalent for your distribution.

### Install hexsmoothR

```r
# Install from GitHub
devtools::install_github("MaxMLang/hexsmoothR")

# Or install from local source
install.packages("hexsmoothR_0.1.0.tar.gz", repos = NULL, type = "source")
```

### Troubleshooting

If you get compilation errors:

1. **Windows**: Make sure Rtools is installed and in your PATH
2. **macOS**: Ensure Xcode Command Line Tools are installed
3. **Linux**: Install build tools for your distribution
4. **All platforms**: Make sure you have the latest version of Rcpp: `install.packages("Rcpp")`

## Quick Start

```r
library(hexsmoothR)
library(sf)
library(terra)

# Create hexagonal grid
hex_grid <- create_grid(study_area, cell_size = 20000, type = "hexagonal")

# Compute spatial topology
topology <- compute_topology(hex_grid)

# Extract raster data
raster_files <- c(ndvi = "path/to/ndvi.tif", elevation = "path/to/elevation.tif")
extracted_data <- extract_raster_data(raster_files, study_area, cell_size = 20000)

# Apply spatial smoothing
smoothed_results <- smooth_variables(
  variable_values = list(ndvi = extracted_data$data$ndvi, elevation = extracted_data$data$elevation),
  neighbors = topology$neighbors,
  weights = topology$weights,
  var_names = c("ndvi", "elevation")
)
```

## Exported Functions

hexsmoothR exports the following public functions:

- **`create_grid`** - Create hexagonal or square grids from study areas
- **`compute_topology`** - Compute neighbor relationships and spatial weights
- **`extract_raster_data`** - Extract values from raster files to grid cells
- **`smooth_variables`** - Apply spatial smoothing with C++ optimization
- **`find_hex_cell_size_for_target_cells`** - Calculate optimal cell size for target number of cells
- **`get_utm_crs`** - Get appropriate UTM CRS for a study area
- **`hex_flat_to_edge`** - Convert between hexagon measurements
- **`hex_edge_to_flat`** - Edge length to flat-to-flat distance
- **`hex_flat_to_circumradius`** - Flat-to-flat to circumradius conversion
- **`hex_circumradius_to_flat`** - Circumradius to flat-to-flat conversion

**Note:** Internal C++ and R fallback functions are not exported and should not be called directly.

## Grid Compatibility

hexsmoothR works with **any hexagonal grid** that has the required structure:

- **Uber H3 grids** from Python/R
- **Custom hexagonal grids** 
- **Square grids** (for comparison)

Required grid structure:
```r
grid <- st_sf(
  geometry = st_sfc(polygons, crs = your_crs),
  grid_id = unique_identifiers,      # Character vector
  grid_index = 1:length(polygons)    # Numeric sequence
)
```

## Output

The smoothing results contain:
- `raw` - Weighted average of center cell and all neighbors
- `neighbors_1st`, `neighbors_2nd`, etc. - Mean of neighbors at each order
- `weighted_combined` - Final weighted average (recommended for analysis)

## Advanced Usage

### Multiple Variables
```r
raster_files <- c(ndvi = "path/to/ndvi.tif", elevation = "path/to/elevation.tif")
extracted_data <- extract_raster_data(raster_files, study_area, cell_size = 20000)

smoothed_results <- smooth_variables(
  variable_values = list(ndvi = extracted_data$data$ndvi, elevation = extracted_data$data$elevation),
  neighbors = topology$neighbors,
  weights = topology$weights,
  var_names = c("ndvi", "elevation")
)
```

### Custom Weights
```r
topology <- compute_topology(hex_grid, neighbor_orders = 3, center_weight = 1.0)
```

## Grid Requirements

hexsmoothR works with any hexagonal grid that follows the required structure. You can use grids created with:

- **Uber H3 grids** from Python, R, or any programming language
- **Custom hexagonal grids** created with any GIS software
- **Grids from other packages** (as long as they follow the structure below)

**Required grid structure:**
```r
grid <- st_sf(
  geometry = st_sfc(polygons, crs = your_crs),
  grid_id = unique_identifiers,      # Character vector
  grid_index = 1:length(polygons)    # Numeric sequence
)
```

**Example with H3 grids from Python:**
```python
import h3
import geopandas as gpd
from shapely.geometry import Polygon

# Create H3 hexagons in Python
h3_hexes = h3.polygon_to_cells(polygon, resolution=8)
geometries = [Polygon(h3.cell_to_boundary(hex_id)) for hex_id in h3_hexes]

# Convert to GeoDataFrame with required structure
gdf = gpd.GeoDataFrame({
    'grid_id': h3_hexes,
    'grid_index': range(len(h3_hexes))
}, geometry=geometries, crs='EPSG:4326')
```

**Tip:** Use projected CRS (UTM) for real-world analysis with cell sizes in meters.

## Performance

- **C++ Optimization**: High-performance spatial smoothing
- **Memory Efficient**: Chunking support for large datasets
- **Scalable**: From small areas to continental scales

## Implementation Details

hexsmoothR uses Rcpp for high-performance spatial smoothing algorithms with automatic fallback to R implementation if needed.

## Dependencies

- **Rcpp**: C++ integration
- **sf**: Spatial data handling
- **terra**: Raster processing
- **exactextractr**: Efficient raster extraction
- **data.table**: Fast data manipulation

## Windows-Specific Notes

**Important for Windows users:** hexsmoothR contains C++ code that must be compiled during installation. This requires:

1. **Rtools**: Download and install from [CRAN Rtools page](https://cran.r-project.org/bin/windows/Rtools/)
2. **Restart R/RStudio** after installing Rtools
3. **Verify installation**: Run `Sys.which("make")` - it should return a path, not `""`

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests
5. Submit a pull request

## License

MIT License - see LICENSE file for details.

## Citation

```bibtex
@software{hexsmoothR2025,
  title={hexsmoothR: Hexagonal Grid Smoothing for Satellite Data},
  author={Max M. Lang},
  year={2025},
  url={https://github.com/maxmlang/hexsmoothR}
}
```

## Support

For issues and questions:
- Check the vignettes: `vignette("hexsmoothR-complete-guide", package = "hexsmoothR")`
- Run examples: `example(create_grid)`
- Report bugs on GitHub
