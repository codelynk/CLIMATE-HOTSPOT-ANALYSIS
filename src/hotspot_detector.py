"""
hotspot_detector.py
--------------------
Implements Getis-Ord Gi* spatial autocorrelation to identify
statistically significant climate stress hotspots and coldspots.
"""

import numpy as np
import geopandas as gpd
import pandas as pd
from scipy import stats
from scipy.spatial import cKDTree
import warnings


def build_spatial_weights(gdf: gpd.GeoDataFrame, method: str = "knn", k: int = 8,
                           distance_threshold: float = None) -> np.ndarray:
    """
    Build a spatial weights matrix for Gi* computation.

    Parameters
    ----------
    gdf : gpd.GeoDataFrame
        GeoDataFrame with geometry column.
    method : str
        'knn' for k-nearest neighbours or 'distance' for distance band.
    k : int
        Number of nearest neighbours (used when method='knn').
    distance_threshold : float
        Distance threshold in CRS units (used when method='distance').

    Returns
    -------
    W : np.ndarray
        Row-standardised spatial weights matrix (n x n).
    """
    coords = np.array([(geom.centroid.x, geom.centroid.y) for geom in gdf.geometry])
    n = len(coords)
    W = np.zeros((n, n))

    tree = cKDTree(coords)

    if method == "knn":
        _, indices = tree.query(coords, k=k + 1)  # +1 because query includes self
        for i, neighbours in enumerate(indices):
            for j in neighbours[1:]:  # skip self
                W[i, j] = 1.0

    elif method == "distance":
        if distance_threshold is None:
            raise ValueError("distance_threshold must be set when method='distance'")
        pairs = tree.query_pairs(distance_threshold)
        for i, j in pairs:
            W[i, j] = 1.0
            W[j, i] = 1.0

    else:
        raise ValueError(f"Unknown method '{method}'. Use 'knn' or 'distance'.")

    # Row standardise
    row_sums = W.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1  # avoid division by zero
    W = W / row_sums

    print(f"[✓] Spatial weights matrix built | Method: {method} | Shape: {W.shape}")
    return W


def getis_ord_gi_star(values: np.ndarray, W: np.ndarray) -> tuple:
    """
    Compute the Getis-Ord Gi* statistic for hotspot analysis.

    A high positive Gi* z-score indicates a hotspot (cluster of high values).
    A high negative Gi* z-score indicates a coldspot (cluster of low values).

    Parameters
    ----------
    values : np.ndarray
        1D array of attribute values (e.g., mean LST per zone).
    W : np.ndarray
        Row-standardised spatial weights matrix (n x n).

    Returns
    -------
    gi_star : np.ndarray
        Gi* z-scores for each location.
    p_values : np.ndarray
        Two-tailed p-values for each z-score.
    """
    n = len(values)
    x_bar = np.mean(values)
    s = np.std(values)

    if s == 0:
        warnings.warn("Standard deviation is zero — all values are identical. Gi* will be zero.")
        return np.zeros(n), np.ones(n)

    gi_star = np.zeros(n)

    for i in range(n):
        w_i = W[i]  # weights for location i
        sum_w = w_i.sum()
        sum_w2 = (w_i ** 2).sum()

        numerator = np.dot(w_i, values) - x_bar * sum_w
        denominator = s * np.sqrt(
            (n * sum_w2 - sum_w ** 2) / (n - 1)
        )

        gi_star[i] = numerator / denominator if denominator != 0 else 0.0

    # Two-tailed p-values from standard normal distribution
    p_values = 2 * (1 - stats.norm.cdf(np.abs(gi_star)))

    print(f"[✓] Gi* computed | Hotspots (p<0.05): {(p_values < 0.05).sum()} / {n} zones")
    return gi_star, p_values


def classify_hotspots(gi_star: np.ndarray, p_values: np.ndarray,
                       confidence_levels: dict = None) -> np.ndarray:
    """
    Classify zones into hotspot/coldspot categories based on Gi* and p-values.

    Parameters
    ----------
    gi_star : np.ndarray
        Gi* z-scores.
    p_values : np.ndarray
        Corresponding p-values.
    confidence_levels : dict, optional
        Mapping of label → p-value threshold. Defaults to 90/95/99% confidence.

    Returns
    -------
    classifications : np.ndarray of str
        Labels for each zone.
    """
    if confidence_levels is None:
        confidence_levels = {
            "Hotspot 99%": 0.01,
            "Hotspot 95%": 0.05,
            "Hotspot 90%": 0.10,
            "Coldspot 90%": 0.10,
            "Coldspot 95%": 0.05,
            "Coldspot 99%": 0.01,
        }

    classifications = np.full(len(gi_star), "Not Significant", dtype=object)

    for i, (z, p) in enumerate(zip(gi_star, p_values)):
        if z > 0 and p <= 0.01:
            classifications[i] = "Hotspot 99%"
        elif z > 0 and p <= 0.05:
            classifications[i] = "Hotspot 95%"
        elif z > 0 and p <= 0.10:
            classifications[i] = "Hotspot 90%"
        elif z < 0 and p <= 0.01:
            classifications[i] = "Coldspot 99%"
        elif z < 0 and p <= 0.05:
            classifications[i] = "Coldspot 95%"
        elif z < 0 and p <= 0.10:
            classifications[i] = "Coldspot 90%"

    counts = pd.Series(classifications).value_counts().to_dict()
    print("[✓] Hotspot classification:")
    for label, count in counts.items():
        print(f"     {label}: {count}")

    return classifications


def run_hotspot_analysis(gdf: gpd.GeoDataFrame, value_col: str,
                          weights_method: str = "knn", k: int = 8) -> gpd.GeoDataFrame:
    """
    Full hotspot detection pipeline on a GeoDataFrame.

    Parameters
    ----------
    gdf : gpd.GeoDataFrame
        Input zones with a numeric value column.
    value_col : str
        Column name containing values to analyse.
    weights_method : str
        Spatial weights method ('knn' or 'distance').
    k : int
        Neighbours for KNN weights.

    Returns
    -------
    result_gdf : gpd.GeoDataFrame
        Input GDF with added columns: gi_star, p_value, hotspot_class.
    """
    # Drop rows with missing values
    valid_mask = gdf[value_col].notna()
    valid_gdf = gdf[valid_mask].copy().reset_index(drop=True)

    if len(valid_gdf) < 10:
        raise ValueError(f"Too few valid zones ({len(valid_gdf)}) for meaningful Gi* analysis.")

    values = valid_gdf[value_col].values
    W = build_spatial_weights(valid_gdf, method=weights_method, k=min(k, len(valid_gdf) - 1))
    gi_star, p_values = getis_ord_gi_star(values, W)
    classifications = classify_hotspots(gi_star, p_values)

    valid_gdf["gi_star"] = gi_star
    valid_gdf["p_value"] = p_values
    valid_gdf["hotspot_class"] = classifications

    return valid_gdf


def temporal_change_analysis(gdf_t1: gpd.GeoDataFrame, gdf_t2: gpd.GeoDataFrame,
                              value_col: str, id_col: str) -> gpd.GeoDataFrame:
    """
    Compare hotspot classifications between two time periods.

    Parameters
    ----------
    gdf_t1 : gpd.GeoDataFrame
        Hotspot results for time period 1 (with hotspot_class column).
    gdf_t2 : gpd.GeoDataFrame
        Hotspot results for time period 2 (with hotspot_class column).
    value_col : str
        Column used for the original values.
    id_col : str
        Column used to join the two GDFs.

    Returns
    -------
    change_gdf : gpd.GeoDataFrame
        Joined GDF with change analysis columns.
    """
    merged = gdf_t1[[id_col, value_col, "hotspot_class", "geometry"]].merge(
        gdf_t2[[id_col, value_col, "hotspot_class"]],
        on=id_col,
        suffixes=("_t1", "_t2"),
    )

    merged["value_change"] = merged[f"{value_col}_t2"] - merged[f"{value_col}_t1"]
    merged["class_change"] = merged["hotspot_class_t1"] + " → " + merged["hotspot_class_t2"]
    merged["worsened"] = (
        merged["hotspot_class_t2"].str.contains("Hotspot") &
        ~merged["hotspot_class_t1"].str.contains("Hotspot")
    )

    print(f"[✓] Temporal change analysis complete | Worsened zones: {merged['worsened'].sum()}")
    return gpd.GeoDataFrame(merged, geometry="geometry", crs=gdf_t1.crs)


if __name__ == "__main__":
    # Quick smoke test with synthetic data
    print("Running smoke test with synthetic data...")
    from shapely.geometry import Point

    np.random.seed(42)
    n = 50
    gdf = gpd.GeoDataFrame(
        {"zone_id": range(n), "lst_mean": np.random.normal(35, 5, n)},
        geometry=[Point(np.random.uniform(3.0, 4.5), np.random.uniform(6.0, 7.5)) for _ in range(n)],
        crs="EPSG:4326",
    )

    result = run_hotspot_analysis(gdf, value_col="lst_mean", k=6)
    print(result[["zone_id", "lst_mean", "gi_star", "p_value", "hotspot_class"]].head(10))
    print("\nSmoke test passed ✓")
