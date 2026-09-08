"""Temperature station map colored by start year, with a horizontal colorbar
below the map instead of an inline legend.

Recreates the "ta" map from plots4paper.ipynb cell 25, styled with the
discrete colorbar sidebar from cell 27 (reoriented horizontal, below the map).
"""

import re

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.gridspec import GridSpec

MAP_DOMAIN = [-12, 61.5, 29.8, 73]
DATE_BLOCK_RE = re.compile(r"(\d{8})[-_](\d{8})")


def extract_keys(filename):
    match = DATE_BLOCK_RE.search(filename)
    if not match:
        return (None, None)
    left = filename[: match.start()].rstrip("_")
    right = filename[match.end() :].lstrip("_").replace(".tsv", "").replace(".TSV", "")
    var = right.split("_")[0]
    return left, var


def load_df_map():
    df_meta = pd.read_csv("metadata_summary.csv")
    df_v2 = pd.read_csv("sef_series_summary_v2.csv")

    keys = df_meta["filename"].apply(lambda x: pd.Series(extract_keys(x)))
    df_meta["station_source"] = keys[0]
    df_meta["variable"] = keys[1]

    df_coords = df_meta.dropna(subset=["station_source", "variable"]).drop_duplicates(
        subset=["station_source", "variable"]
    )
    df_coords = df_coords[["station_source", "variable", "Lat", "Lon"]]

    return df_v2.merge(df_coords, on=["station_source", "variable"], how="left")


def main():
    df_map = load_df_map()

    start_year = df_map["start"].str.split("-").str[0].astype(int)
    bins = [-9999, 1700, 1750, 1800, 1850, 1900]
    labels = ["<1700", "1700–1749", "1750–1799", "1800–1849", ">1850"]
    df_map["start_period"] = pd.cut(start_year, bins=bins, labels=labels)

    colors = sns.color_palette("plasma", len(labels))[::-1]
    cats = df_map["start_period"].cat.categories

    df_var = df_map[df_map["variable"] == "ta"].copy()
    df_var = df_var.sort_values(by="start", ascending=False)
    codes = df_var["start_period"].cat.codes
    point_colors = [colors[c] for c in codes]

    fig = plt.figure(figsize=(2.3, 3.2))
    gs = GridSpec(2, 1, height_ratios=[1, 0.05], hspace=0.25, figure=fig)

    ax = fig.add_subplot(gs[0], projection=ccrs.PlateCarree())
    ax.set_extent(MAP_DOMAIN, crs=ccrs.PlateCarree())

    ax.set_aspect('auto')  # Let Cartopy handle the aspect ratio for geographic accuracy
    ax.add_feature(cfeature.COASTLINE, linewidth=0.1)
    ax.add_feature(cfeature.LAND, facecolor='#f4f4f4', edgecolor='black')
    ax.spines["geo"].set_visible(False)

    gl = ax.gridlines(draw_labels=True, linewidth=0.3, alpha=0.3, color="gray", linestyle="--")
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {"size": 6}
    gl.ylabel_style = {"size": 6}

    ax.scatter(
        df_var["Lon"].astype(float),
        df_var["Lat"].astype(float),
        facecolor=point_colors,
        s=10,
        edgecolor="black",
        linewidth=0.2,
        marker="o",
        transform=ccrs.PlateCarree(),
        alpha=0.9,
        zorder=3,
    )

    ax.set_title("Temperature Series in HIST-DAILY", fontsize=7, fontweight="bold", pad=3)

    cmap_discrete = ListedColormap(colors)
    bounds = np.arange(-0.5, len(labels), 1)
    norm = BoundaryNorm(bounds, cmap_discrete.N)

    cbar_ax = fig.add_subplot(gs[1])
    sm = plt.cm.ScalarMappable(cmap=cmap_discrete, norm=norm)
    sm.set_array([])

    cbar = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal", extend="both", extendfrac=0.05)
    cbar.ax.xaxis.set_tick_params(which="both", length=0)
    cbar.set_ticks([0.5, 1.5, 2.5, 3.5])
    cbar.set_ticklabels(["1700", "1750", "1800", "1850"], fontsize=6)
    cbar.set_label("Start Year", fontsize=7, fontweight="bold", labelpad=3)

    plt.savefig("image/station_map_start_year_ta_horizontal.svg", bbox_inches="tight")
    plt.savefig("image/station_map_start_year_ta_horizontal.pdf", bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()
