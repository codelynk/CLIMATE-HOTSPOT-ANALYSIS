"""
visualiser.py
-------------
Generates static and interactive maps for the Climate Hotspot Analysis project.
"""

import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import folium
from folium.plugins import HeatMap
import os


# ── Colour palette for hotspot classifications ──────────────────────────────
HOTSPOT_PALETTE = {
    "Hotspot 99%":      "#d7191c",
    "Hotspot 95%":      "#fdae61",
    "Hotspot 90%":      "#fee08b",
    "Not Significant":  "#ffffbf",
    "Coldspot 90%":     "#a6d96a",
    "Coldspot 95%":     "#1a9641",
    "Coldspot 99%":     "#006837",
}

ORDERED_CLASSES = list(HOTSPOT_PALETTE.keys())


def plot_hotspot_map(result_gdf: gpd.GeoDataFrame, title: str = "Climate Stress Hotspot Map",
                     figsize: tuple = (12, 10), output_path: str = None) -> plt.Figure:
    """
    Plot a static choropleth hotspot map using Matplotlib.

    Parameters
    ----------
    result_gdf : gpd.GeoDataFrame
        Must contain 'hotspot_class' and 'geometry' columns.
    title : str
        Map title.
    figsize : tuple
        Figure size in inches.
    output_path : str, optional
        If provided, saves the figure to this path.

    Returns
    -------
    fig : plt.Figure
    """
    fig, ax = plt.subplots(1, 1, figsize=figsize, facecolor="#1a1a2e")
    ax.set_facecolor("#1a1a2e")

    for cls in ORDERED_CLASSES:
        subset = result_gdf[result_gdf["hotspot_class"] == cls]
        if not subset.empty:
            subset.plot(ax=ax, color=HOTSPOT_PALETTE[cls], edgecolor="white",
                        linewidth=0.3, alpha=0.9)

    # Legend
    patches = [
        mpatches.Patch(color=HOTSPOT_PALETTE[cls], label=cls)
        for cls in ORDERED_CLASSES
        if cls in result_gdf["hotspot_class"].values
    ]
    legend = ax.legend(
        handles=patches, loc="lower left", frameon=True,
        framealpha=0.85, fontsize=9, title="Hotspot Classification",
        title_fontsize=10, facecolor="#2d2d44", labelcolor="white",
    )
    legend.get_title().set_color("white")

    ax.set_title(title, color="white", fontsize=15, fontweight="bold", pad=15)
    ax.set_xlabel("Longitude", color="#aaaacc", fontsize=9)
    ax.set_ylabel("Latitude", color="#aaaacc", fontsize=9)
    ax.tick_params(colors="#aaaacc")
    for spine in ax.spines.values():
        spine.set_edgecolor("#444466")

    plt.tight_layout()

    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        fig.savefig(output_path, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
        print(f"[✓] Static map saved to: {output_path}")

    return fig


def plot_lst_distribution(gdf: gpd.GeoDataFrame, value_col: str = "lst_mean",
                           group_col: str = "hotspot_class", output_path: str = None) -> plt.Figure:
    """
    Box plot of LST values grouped by hotspot classification.
    """
    import matplotlib.gridspec as gridspec

    fig = plt.figure(figsize=(13, 5), facecolor="#1a1a2e")
    gs = gridspec.GridSpec(1, 2, width_ratios=[2, 1])
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])

    for ax in [ax1, ax2]:
        ax.set_facecolor("#2d2d44")

    ordered = [c for c in ORDERED_CLASSES if c in gdf[group_col].values]
    data_by_class = [gdf[gdf[group_col] == cls][value_col].dropna().values for cls in ordered]
    colours = [HOTSPOT_PALETTE[cls] for cls in ordered]

    bp = ax1.boxplot(data_by_class, patch_artist=True, notch=False, vert=True,
                     medianprops=dict(color="white", linewidth=2))
    for patch, colour in zip(bp["boxes"], colours):
        patch.set_facecolor(colour)
        patch.set_alpha(0.85)

    ax1.set_xticklabels(ordered, rotation=30, ha="right", color="white", fontsize=8)
    ax1.set_ylabel("Land Surface Temperature (°C)", color="white")
    ax1.set_title("LST Distribution by Hotspot Class", color="white", fontweight="bold")
    ax1.tick_params(colors="white")
    ax1.spines["bottom"].set_color("#555577")
    ax1.spines["left"].set_color("#555577")
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.yaxis.label.set_color("white")

    # Count bar chart
    counts = gdf[group_col].value_counts().reindex(ordered, fill_value=0)
    bars = ax2.barh(ordered, counts.values, color=colours, alpha=0.85, edgecolor="white", linewidth=0.4)
    ax2.set_xlabel("Number of Zones", color="white")
    ax2.set_title("Zone Counts", color="white", fontweight="bold")
    ax2.tick_params(colors="white")
    ax2.spines["bottom"].set_color("#555577")
    ax2.spines["left"].set_color("#555577")
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.set_yticklabels(ordered, color="white", fontsize=8)

    for bar, val in zip(bars, counts.values):
        ax2.text(val + 0.3, bar.get_y() + bar.get_height() / 2,
                 str(val), va="center", color="white", fontsize=8)

    fig.patch.set_facecolor("#1a1a2e")
    plt.tight_layout()

    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        fig.savefig(output_path, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
        print(f"[✓] Distribution plot saved to: {output_path}")

    return fig


def create_interactive_map(result_gdf: gpd.GeoDataFrame, value_col: str = "lst_mean",
                            output_path: str = None) -> folium.Map:
    """
    Generate an interactive Folium map with hotspot classification popups.

    Parameters
    ----------
    result_gdf : gpd.GeoDataFrame
        Must contain 'hotspot_class', value_col, and geometry.
    value_col : str
        Column name for the temperature values.
    output_path : str, optional
        If provided, saves the map as an HTML file.

    Returns
    -------
    m : folium.Map
    """
    centroid = result_gdf.geometry.unary_union.centroid
    m = folium.Map(location=[centroid.y, centroid.x], zoom_start=9,
                   tiles="CartoDB dark_matter")

    def style_fn(feature):
        cls = feature["properties"].get("hotspot_class", "Not Significant")
        return {
            "fillColor": HOTSPOT_PALETTE.get(cls, "#ffffbf"),
            "color": "white",
            "weight": 0.5,
            "fillOpacity": 0.75,
        }

    def highlight_fn(feature):
        return {"weight": 2, "color": "#ffffff", "fillOpacity": 0.95}

    tooltip_cols = ["hotspot_class", value_col]
    if "gi_star" in result_gdf.columns:
        tooltip_cols.append("gi_star")
    if "p_value" in result_gdf.columns:
        tooltip_cols.append("p_value")

    folium.GeoJson(
        result_gdf[tooltip_cols + ["geometry"]].to_json(),
        style_function=style_fn,
        highlight_function=highlight_fn,
        tooltip=folium.GeoJsonTooltip(
            fields=tooltip_cols,
            aliases=["Class", "LST (°C)", "Gi* Z-score", "P-value"][:len(tooltip_cols)],
            localize=True,
            sticky=False,
        ),
    ).add_to(m)

    # Legend
    legend_html = """
    <div style="position: fixed; bottom: 30px; left: 30px; z-index: 9999;
         background-color: rgba(20,20,40,0.88); padding: 12px 16px;
         border-radius: 8px; border: 1px solid #444466; font-family: monospace;">
      <b style="color:white; font-size:12px;">Hotspot Classification</b><br>
    """
    for cls, colour in HOTSPOT_PALETTE.items():
        legend_html += f"""
        <div style="display:flex;align-items:center;margin-top:4px;">
          <div style="width:14px;height:14px;background:{colour};
               border-radius:3px;margin-right:8px;border:1px solid #888;"></div>
          <span style="color:white;font-size:11px;">{cls}</span>
        </div>"""
    legend_html += "</div>"

    m.get_root().html.add_child(folium.Element(legend_html))

    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        m.save(output_path)
        print(f"[✓] Interactive map saved to: {output_path}")

    return m


def plot_temporal_change(change_gdf: gpd.GeoDataFrame, output_path: str = None) -> plt.Figure:
    """
    Visualise temperature change between two time periods.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), facecolor="#1a1a2e")

    vmax = change_gdf["value_change"].abs().quantile(0.95)
    divnorm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    cmap = plt.cm.RdBu_r

    for ax in axes:
        ax.set_facecolor("#1a1a2e")

    change_gdf.plot(column="value_change", cmap=cmap, norm=divnorm,
                    ax=axes[0], edgecolor="none", linewidth=0)
    axes[0].set_title("Temperature Change (T2 − T1)", color="white", fontweight="bold")
    axes[0].tick_params(colors="#aaaacc")

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=divnorm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes[0], fraction=0.03, pad=0.04)
    cbar.set_label("ΔT (°C)", color="white")
    cbar.ax.yaxis.set_tick_params(color="white")
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color="white")

    change_gdf[change_gdf["worsened"]].plot(
        ax=axes[1], color="#ff4444", edgecolor="white", linewidth=0.4, alpha=0.85
    )
    change_gdf[~change_gdf["worsened"]].plot(
        ax=axes[1], color="#2d2d44", edgecolor="#444466", linewidth=0.3, alpha=0.6
    )
    axes[1].set_title(f"Newly Worsened Hotspots\n({change_gdf['worsened'].sum()} zones)",
                       color="white", fontweight="bold")
    axes[1].tick_params(colors="#aaaacc")

    fig.patch.set_facecolor("#1a1a2e")
    plt.tight_layout()

    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        fig.savefig(output_path, dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor())
        print(f"[✓] Temporal change map saved to: {output_path}")

    return fig
