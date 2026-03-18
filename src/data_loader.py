"""
data_loader.py
--------------
Handles loading, validation, and preprocessing of spatial datasets
for the Climate Hotspot Analysis pipeline.
"""

import os
import numpy as np
import geopandas as gpd
import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.mask import mask
from shapely.geometry import mapping
import pandas as pd
from pathlib import Path


TARGET_CRS = "EPSG:4326"


def load_raster(filepath: str) -> tuple:
    """
    Load a raster file and return the data array and metadata.

    Parameters
    ----------
    filepath : str
        Path to a GeoTIFF or compatible raster file.

    Returns
    -------
    data : np.ndarray
        2D array of raster values with nodata replaced by np.nan.
    meta : dict
        Rasterio metadata dictionary.
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"Raster file not found: {filepath}")

    with rasterio.open(filepath) as src:
        data = src.read(1).astype(np.float32)
        meta = src.meta.copy()
        nodata = src.nodata

        if nodata is not None:
            data[data == nodata] = np.nan

        # Reproject to target CRS if needed
        if src.crs and src.crs.to_string() != TARGET_CRS:
            data, meta = _reproject_raster(src, TARGET_CRS)

    print(f"[✓] Loaded raster: {path.name} | Shape: {data.shape} | CRS: {meta['crs']}")
    return data, meta


def _reproject_raster(src, target_crs: str) -> tuple:
    """Reproject an open rasterio dataset to a target CRS."""
    transform, width, height = calculate_default_transform(
        src.crs, target_crs, src.width, src.height, *src.bounds
    )
    meta = src.meta.copy()
    meta.update({"crs": target_crs, "transform": transform, "width": width, "height": height})

    data = np.empty((height, width), dtype=np.float32)
    reproject(
        source=rasterio.band(src, 1),
        destination=data,
        src_transform=src.transform,
        src_crs=src.crs,
        dst_transform=transform,
        dst_crs=target_crs,
        resampling=Resampling.bilinear,
    )
    return data, meta


def load_vector(filepath: str) -> gpd.GeoDataFrame:
    """
    Load a vector file (shapefile, GeoJSON, GeoPackage) as a GeoDataFrame.

    Parameters
    ----------
    filepath : str
        Path to the vector file.

    Returns
    -------
    gdf : gpd.GeoDataFrame
        Reprojected to TARGET_CRS.
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"Vector file not found: {filepath}")

    gdf = gpd.read_file(filepath)

    if gdf.crs is None:
        raise ValueError(f"Vector file has no CRS defined: {filepath}")

    if gdf.crs.to_string() != TARGET_CRS:
        gdf = gdf.to_crs(TARGET_CRS)

    print(f"[✓] Loaded vector: {path.name} | Features: {len(gdf)} | CRS: {gdf.crs}")
    return gdf


def clip_raster_to_boundary(raster_path: str, boundary_gdf: gpd.GeoDataFrame) -> tuple:
    """
    Clip a raster to the extent of a vector boundary.

    Parameters
    ----------
    raster_path : str
        Path to the raster file.
    boundary_gdf : gpd.GeoDataFrame
        Boundary polygon(s) to clip to.

    Returns
    -------
    clipped_data : np.ndarray
    clipped_meta : dict
    """
    shapes = [mapping(geom) for geom in boundary_gdf.geometry]

    with rasterio.open(raster_path) as src:
        clipped_data, clipped_transform = mask(src, shapes, crop=True, nodata=np.nan)
        clipped_meta = src.meta.copy()
        clipped_meta.update({
            "height": clipped_data.shape[1],
            "width": clipped_data.shape[2],
            "transform": clipped_transform,
        })

    clipped_data = clipped_data[0].astype(np.float32)
    clipped_data[clipped_data == src.nodata] = np.nan

    print(f"[✓] Raster clipped | New shape: {clipped_data.shape}")
    return clipped_data, clipped_meta


def compute_zonal_statistics(raster_data: np.ndarray, meta: dict, zones_gdf: gpd.GeoDataFrame,
                              stat: str = "mean") -> gpd.GeoDataFrame:
    """
    Compute zonal statistics for each polygon in zones_gdf over a raster.

    Parameters
    ----------
    raster_data : np.ndarray
        2D raster array.
    meta : dict
        Rasterio metadata with transform and CRS.
    zones_gdf : gpd.GeoDataFrame
        Polygon zones for aggregation.
    stat : str
        Statistic to compute: 'mean', 'max', 'min', 'std'.

    Returns
    -------
    result_gdf : gpd.GeoDataFrame
        Input GDF with an additional column for the computed statistic.
    """
    from rasterio.features import geometry_mask

    stat_funcs = {"mean": np.nanmean, "max": np.nanmax, "min": np.nanmin, "std": np.nanstd}
    if stat not in stat_funcs:
        raise ValueError(f"Unsupported stat '{stat}'. Choose from: {list(stat_funcs.keys())}")

    results = []
    transform = meta["transform"]

    for _, row in zones_gdf.iterrows():
        geom_mask = geometry_mask(
            [mapping(row.geometry)],
            transform=transform,
            invert=True,
            out_shape=raster_data.shape,
        )
        masked_vals = raster_data[geom_mask]
        value = stat_funcs[stat](masked_vals) if masked_vals.size > 0 else np.nan
        results.append(value)

    result_gdf = zones_gdf.copy()
    result_gdf[f"lst_{stat}"] = results
    print(f"[✓] Zonal statistics computed | Stat: {stat} | Zones: {len(result_gdf)}")
    return result_gdf


def kelvin_to_celsius(data: np.ndarray, scale_factor: float = 0.02) -> np.ndarray:
    """
    Convert MODIS LST values to Celsius.

    MODIS LST is stored as scaled integers: actual_temp_K = value * 0.02
    Then convert K → °C by subtracting 273.15.

    Parameters
    ----------
    data : np.ndarray
        Raw MODIS LST values.
    scale_factor : float
        MODIS scale factor (default 0.02).

    Returns
    -------
    np.ndarray of temperatures in °C.
    """
    celsius = (data * scale_factor) - 273.15
    celsius[celsius < -100] = np.nan  # mask invalid
    print(f"[✓] Converted LST to Celsius | Range: {np.nanmin(celsius):.1f}°C – {np.nanmax(celsius):.1f}°C")
    return celsius


if __name__ == "__main__":
    print("data_loader.py — run via the analysis pipeline or notebooks.")
