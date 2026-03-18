# 🌍 Climate Hotspot Analysis

A geospatial analysis project that identifies and visualises **urban heat islands and climate stress hotspots** using satellite-derived temperature data and land cover classification.

---

## 📌 Project Overview

This project analyses land surface temperature (LST) data alongside land use/land cover (LULC) datasets to detect areas under climate stress. It generates choropleth maps, hotspot clusters, and statistical summaries to support environmental research and urban planning decisions.

**Key Questions Answered:**
- Where are the most intense urban heat islands?
- How does land cover type correlate with surface temperature?
- Which areas have worsened over a 10-year period?

---

## 🗂️ Project Structure

```
climate-hotspot-analysis/
├── data/
│   ├── raw/                  # Raw GeoTIFF / shapefiles (not tracked by git)
│   └── processed/            # Cleaned and reprojected datasets
├── notebooks/
│   ├── 01_data_exploration.ipynb
│   ├── 02_hotspot_detection.ipynb
│   └── 03_visualisation.ipynb
├── src/
│   ├── __init__.py
│   ├── data_loader.py        # Load and validate spatial datasets
│   ├── hotspot_detector.py   # Getis-Ord Gi* hotspot detection
│   ├── land_cover.py         # LULC classification helpers
│   └── visualiser.py         # Map generation and plotting
├── outputs/
│   ├── maps/                 # Exported PNG/HTML maps
│   └── reports/              # Summary statistics CSV
├── tests/
│   └── test_hotspot_detector.py
├── requirements.txt
├── .gitignore
└── README.md
```

---

## 🔬 Methodology

1. **Data Acquisition** — MODIS LST (MOD11A1) and ESA CCI Land Cover data
2. **Preprocessing** — Reprojection to EPSG:4326, cloud masking, temporal averaging
3. **Hotspot Detection** — Getis-Ord Gi* spatial autocorrelation statistic
4. **Land Cover Correlation** — Zonal statistics per LULC class
5. **Temporal Analysis** — Decadal comparison (2013 vs 2023)
6. **Visualisation** — Interactive Folium maps + static Matplotlib figures

---

## 🛠️ Tech Stack

| Tool | Purpose |
|------|---------|
| `geopandas` | Vector spatial data manipulation |
| `rasterio` | Raster data I/O and processing |
| `numpy` | Numerical computation |
| `scipy` | Spatial statistics |
| `matplotlib` / `seaborn` | Static visualisation |
| `folium` | Interactive web maps |
| `pysal` | Spatial autocorrelation (Gi*) |
| `pandas` | Tabular data handling |

---

## 🚀 Getting Started

### Prerequisites
- Python 3.10+
- pip

### Installation

```bash
git clone https://github.com/yourusername/climate-hotspot-analysis.git
cd climate-hotspot-analysis
pip install -r requirements.txt
```

### Running the Analysis

```bash
# Run the full pipeline
python src/hotspot_detector.py --input data/raw/lst_2023.tif --output outputs/maps/

# Or step through the notebooks
jupyter notebook notebooks/01_data_exploration.ipynb
```

---

## 📊 Sample Output

The analysis produces:
- **Hotspot maps** highlighting statistically significant heat clusters (p < 0.05)
- **Land cover overlays** showing temperature by urban/vegetation/water class
- **Trend charts** comparing 2013 vs 2023 mean temperatures per district

---

## 📁 Data Sources

| Dataset | Source | Resolution |
|---------|--------|------------|
| MODIS LST (MOD11A1) | [NASA EarthData](https://earthdata.nasa.gov/) | 1 km / daily |
| ESA CCI Land Cover | [ESA Climate Office](https://www.esa-landcover-cci.org/) | 300 m / annual |
| Administrative Boundaries | [GADM](https://gadm.org/) | Vector |

> **Note:** Raw data files are not included in this repository due to size. See `data/README.md` for download instructions.

---

## 🧪 Running Tests

```bash
pytest tests/
```


---

## 🙋 Author

** Patrick Nkwachukwu Ezeh **\ 
·[Email](mailto:patrickezeh02@gmail.com)
