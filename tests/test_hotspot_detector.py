"""
test_hotspot_detector.py
------------------------
Unit tests for the hotspot_detector module.
"""

import pytest
import numpy as np
import geopandas as gpd
from shapely.geometry import Point
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from hotspot_detector import (
    build_spatial_weights,
    getis_ord_gi_star,
    classify_hotspots,
    run_hotspot_analysis,
)


@pytest.fixture
def synthetic_gdf():
    """Create a small synthetic GeoDataFrame for testing."""
    np.random.seed(42)
    n = 30
    return gpd.GeoDataFrame(
        {
            "zone_id": range(n),
            "lst_mean": np.concatenate([
                np.random.normal(45, 2, 10),   # high temps — expected hotspot
                np.random.normal(30, 2, 10),   # low temps — expected coldspot
                np.random.normal(37, 3, 10),   # mixed
            ]),
        },
        geometry=[
            Point(np.random.uniform(3.0, 4.5), np.random.uniform(6.0, 7.5))
            for _ in range(n)
        ],
        crs="EPSG:4326",
    )


class TestBuildSpatialWeights:

    def test_knn_shape(self, synthetic_gdf):
        W = build_spatial_weights(synthetic_gdf, method="knn", k=4)
        n = len(synthetic_gdf)
        assert W.shape == (n, n), "Weights matrix shape should be (n, n)"

    def test_knn_row_standardised(self, synthetic_gdf):
        W = build_spatial_weights(synthetic_gdf, method="knn", k=4)
        row_sums = W.sum(axis=1)
        # Each row should sum to 1 (or 0 for isolated points)
        assert np.allclose(row_sums[row_sums > 0], 1.0, atol=1e-6)

    def test_distance_method(self, synthetic_gdf):
        W = build_spatial_weights(synthetic_gdf, method="distance", distance_threshold=1.0)
        assert W.shape == (len(synthetic_gdf), len(synthetic_gdf))

    def test_unknown_method_raises(self, synthetic_gdf):
        with pytest.raises(ValueError, match="Unknown method"):
            build_spatial_weights(synthetic_gdf, method="invalid")

    def test_distance_missing_threshold_raises(self, synthetic_gdf):
        with pytest.raises(ValueError, match="distance_threshold"):
            build_spatial_weights(synthetic_gdf, method="distance")


class TestGetisOrdGiStar:

    def test_output_length(self, synthetic_gdf):
        W = build_spatial_weights(synthetic_gdf, method="knn", k=4)
        gi, p = getis_ord_gi_star(synthetic_gdf["lst_mean"].values, W)
        assert len(gi) == len(synthetic_gdf)
        assert len(p) == len(synthetic_gdf)

    def test_p_values_in_range(self, synthetic_gdf):
        W = build_spatial_weights(synthetic_gdf, method="knn", k=4)
        _, p = getis_ord_gi_star(synthetic_gdf["lst_mean"].values, W)
        assert np.all(p >= 0) and np.all(p <= 1), "P-values must be in [0, 1]"

    def test_uniform_values_returns_zeros(self, synthetic_gdf):
        W = build_spatial_weights(synthetic_gdf, method="knn", k=4)
        uniform_values = np.ones(len(synthetic_gdf)) * 35.0
        gi, _ = getis_ord_gi_star(uniform_values, W)
        assert np.allclose(gi, 0.0, atol=1e-6)


class TestClassifyHotspots:

    def test_output_length(self):
        gi = np.array([3.5, -3.5, 0.5, -0.5])
        p = np.array([0.001, 0.001, 0.5, 0.5])
        result = classify_hotspots(gi, p)
        assert len(result) == 4

    def test_hotspot_99_detected(self):
        gi = np.array([4.0])
        p = np.array([0.005])
        result = classify_hotspots(gi, p)
        assert result[0] == "Hotspot 99%"

    def test_coldspot_99_detected(self):
        gi = np.array([-4.0])
        p = np.array([0.005])
        result = classify_hotspots(gi, p)
        assert result[0] == "Coldspot 99%"

    def test_not_significant(self):
        gi = np.array([0.3])
        p = np.array([0.75])
        result = classify_hotspots(gi, p)
        assert result[0] == "Not Significant"


class TestRunHotspotsAnalysis:

    def test_returns_geodataframe(self, synthetic_gdf):
        result = run_hotspot_analysis(synthetic_gdf, value_col="lst_mean", k=5)
        assert isinstance(result, gpd.GeoDataFrame)

    def test_columns_present(self, synthetic_gdf):
        result = run_hotspot_analysis(synthetic_gdf, value_col="lst_mean", k=5)
        assert "gi_star" in result.columns
        assert "p_value" in result.columns
        assert "hotspot_class" in result.columns

    def test_all_classes_valid(self, synthetic_gdf):
        valid_classes = {
            "Hotspot 99%", "Hotspot 95%", "Hotspot 90%",
            "Not Significant",
            "Coldspot 90%", "Coldspot 95%", "Coldspot 99%",
        }
        result = run_hotspot_analysis(synthetic_gdf, value_col="lst_mean", k=5)
        assert set(result["hotspot_class"].unique()).issubset(valid_classes)

    def test_too_few_zones_raises(self):
        tiny_gdf = gpd.GeoDataFrame(
            {"zone_id": range(5), "lst_mean": [30, 35, 40, 38, 32]},
            geometry=[Point(i, i) for i in range(5)],
            crs="EPSG:4326",
        )
        with pytest.raises(ValueError, match="Too few valid zones"):
            run_hotspot_analysis(tiny_gdf, value_col="lst_mean", k=4)
