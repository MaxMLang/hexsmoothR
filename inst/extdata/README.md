# Sample Data Documentation

## sentinel2_ndvi.tif

This file contains real Sentinel-2 L2A NDVI (Normalized Difference Vegetation Index) data for demonstration purposes.

### Data Source
- **Satellite**: Sentinel-2 L2A (Level-2A Atmospheric Correction)
- **Platform**: Copernicus Hub Ecosystem
- **Source**: [Copernicus Data Space Ecosystem](https://browser.dataspace.copernicus.eu/)
- **Processing**: Level-2A atmospheric correction applied
- **Region**: Coastal area with complex shorelines
- **Resolution**: 10m (native Sentinel-2 resolution)
- **Date**: July 30, 2025
- **Cloud Filter**: < 20% cloud cover

### Sentinel-2 L2A Details
The data was obtained from the Copernicus Hub Ecosystem, which provides:
- **Atmospheric Correction**: Level-2A processing removes atmospheric effects
- **Surface Reflectance**: Calibrated reflectance values for accurate analysis
- **Quality Control**: Automated quality assessment and filtering
- **Open Access**: Free and open data under Copernicus license

### NDVI Calculation
NDVI was calculated using Sentinel-2 bands:
- **NIR Band**: B8A (Near-Infrared, 865nm)
- **Red Band**: B4 (Red, 665nm)
- **Formula**: (B8A - B4) / (B8A + B4)

### Citation
When using this sample data in publications, please cite:

```
European Space Agency (2021). Sentinel-2 Level-2A Product. 
Copernicus Open Access Hub. 
https://browser.dataspace.copernicus.eu/
```

### Usage in R
```r
# Access the Sentinel-2 sample data
sentinel_file <- system.file("extdata", "sentinel2_ndvi.tif", package = "hexsmoothR")
sentinel_raster <- terra::rast(sentinel_file)

# View basic information
print(sentinel_raster)
```

### Data Characteristics
- **Format**: GeoTIFF
- **CRS**: WGS84 (EPSG:4326)
- **Extent**: Coastal region with complex shorelines
- **Values**: NDVI ranging from -1 to 1
- **Purpose**: Demonstration of hexsmoothR package with real coastal data
- **Special Features**: Ideal for demonstrating hexagonal smoothing on complex shorelines

### License
This sample data is provided under the Copernicus license for educational and testing purposes. 
For research applications, please obtain your own Sentinel-2 data through 
the [Copernicus Hub Ecosystem](https://browser.dataspace.copernicus.eu/).

---

## default.tif

This file contains sample NDVI (Normalized Difference Vegetation Index) data for demonstration purposes.

### Data Source
- **Satellite**: Landsat 8 Collection 2 Tier 1 Surface Reflectance
- **Processing**: Google Earth Engine (GEE)
- **Region**: Mediterranean area (approximately 50,000 km²)
- **Resolution**: 1km
- **Date Range**: June-August 2023
- **Cloud Filter**: < 20% cloud cover

### GEE Processing Details
The data was generated using the following Google Earth Engine workflow:

1. **Region of Interest**: Rectangle covering 2.0° longitude × 2.25° latitude
2. **Data Collection**: Landsat 8 Surface Reflectance (LANDSAT/LC08/C02/T1_L2)
3. **Temporal Filter**: Summer 2023 (June 1 - August 31)
4. **Cloud Filter**: Images with < 20% cloud cover
5. **Compositing**: Median composite for cloud-free representation
6. **NDVI Calculation**: (NIR - Red) / (NIR + Red) using bands SR_B5 and SR_B4
7. **Export**: 1km resolution GeoTIFF

### Citation
When using this sample data in publications, please cite:

```
USGS Earth Resources Observation and Science Center (2021). 
Landsat 8 Collection 2 Tier 1 Surface Reflectance. 
Google Earth Engine Data Catalog. 
https://developers.google.com/earth-engine/datasets/catalog/LANDSAT_LC08_C02_T1_L2
```

### Usage in R
```r
# Access the sample data
sample_file <- system.file("extdata", "default.tif", package = "hexsmoothR")
sample_raster <- terra::rast(sample_file)

# View basic information
print(sample_raster)
```

### Data Characteristics
- **Format**: GeoTIFF
- **CRS**: WGS84 (EPSG:4326)
- **Extent**: Mediterranean region
- **Values**: NDVI ranging from -1 to 1
- **Purpose**: Demonstration and testing of hexsmoothR package functionality

### License
This sample data is provided for educational and testing purposes only. 
For research applications, please obtain your own satellite data through 
appropriate channels (Google Earth Engine, USGS Earth Explorer, etc.). 