# -*- coding: utf-8 -*-
"""
Seasonal 2D maps of extreme events
Seasons: DJF, MAM, JJA, SON
Variables: Rainfall, Tmax, Tmin — 95th and 99th percentiles
Rain 99th uses corrected file
"""

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
from shapely.vectorized import contains
import os

# --------------------------------------------------
# File paths
# --------------------------------------------------
base_dir = r"C:\Users\khvj92\OneDrive - Durham University\Met Office\95&99_Extremes"
shp_file = r"C:\Users\khvj92\OneDrive - Durham University\Met Office\NPg-full.shp"

files = {
    'rain_95' : os.path.join(base_dir, "Rain_95pct_Events.csv"),
    'rain_99' : os.path.join(base_dir, "Exceedance_99pct_Events_corrected.csv"),  # corrected
    'tmax_95' : os.path.join(base_dir, "Tmax_95th_Exceedance_Events.csv"),
    'tmax_99' : os.path.join(base_dir, "Tmax_99th_Exceedance_Events.csv"),
    'tmin_95' : os.path.join(base_dir, "Tmin_ColdExtreme_95th.csv"),
    'tmin_99' : os.path.join(base_dir, "Tmin_ColdExtreme_99th.csv"),
}

# --------------------------------------------------
# Seasons
# --------------------------------------------------
seasons = {
    'DJF' : [12, 1, 2],
    'MAM' : [3,  4, 5],
    'JJA' : [6,  7, 8],
    'SON' : [9, 10, 11],
}
season_labels = {
    'DJF' : 'December–January–February',
    'MAM' : 'March–April–May',
    'JJA' : 'June–July–August',
    'SON' : 'September–October–November',
}

# --------------------------------------------------
# Dataset definitions
# --------------------------------------------------
col_info = {
    'rain_95' : dict(
        val_col='rainfall_mm',  val_agg='max',
        excess_col='excess_mm', excess_agg='max',
        val_cbar='Rainfall (mm/day)',
        excess_cbar='Rainfall Excess (mm/day)',
        val_cmap='YlGnBu', excess_cmap='Blues',
        period='1990–2022', pct='95th', name='Rainfall',
    ),
    'rain_99' : dict(
        val_col='rainfall_mm',  val_agg='max',
        excess_col='excess_mm', excess_agg='max',
        val_cbar='Rainfall (mm/day)',
        excess_cbar='Rainfall Excess (mm/day)',
        val_cmap='YlOrRd', excess_cmap='OrRd',
        period='1990–2022', pct='99th', name='Rainfall',
    ),
    'tmax_95' : dict(
        val_col='tmax_c',       val_agg='max',
        excess_col='excess_c',  excess_agg='max',
        val_cbar='Temperature (°C)',
        excess_cbar='Temperature Excess (°C)',
        val_cmap='YlOrRd', excess_cmap='Reds',
        period='1990–2022', pct='95th', name='Tmax',
    ),
    'tmax_99' : dict(
        val_col='tmax_c',       val_agg='max',
        excess_col='excess_c',  excess_agg='max',
        val_cbar='Temperature (°C)',
        excess_cbar='Temperature Excess (°C)',
        val_cmap='OrRd', excess_cmap='Reds',
        period='1990–2022', pct='99th', name='Tmax',
    ),
    'tmin_95' : dict(
        val_col='tmin_c',               val_agg='min',
        excess_col='below_threshold_c', excess_agg='max',
        val_cbar='Temperature (°C)',
        excess_cbar='Temperature Below Threshold (°C)',
        val_cmap='YlGnBu', excess_cmap='Blues',
        period='1990–2022', pct='95th', name='Tmin',
    ),
    'tmin_99' : dict(
        val_col='tmin_c',               val_agg='min',
        excess_col='below_threshold_c', excess_agg='max',
        val_cbar='Temperature (°C)',
        excess_cbar='Temperature Below Threshold (°C)',
        val_cmap='PuBu', excess_cmap='Blues',
        period='1990–2022', pct='99th', name='Tmin',
    ),
}

# --------------------------------------------------
# Load shapefile
# --------------------------------------------------
print("Loading shapefile ...")
gdf = gpd.read_file(shp_file)
if gdf.crs is None or gdf.crs.to_epsg() != 4326:
    gdf = gdf.to_crs(epsg=4326)
minx, miny, maxx, maxy = gdf.total_bounds
try:
    region_wgs84 = gdf.geometry.union_all()
except AttributeError:
    region_wgs84 = gdf.geometry.unary_union
extent = [minx, maxx, miny, maxy]

# --------------------------------------------------
# Build regular WGS84 grid and mask
# --------------------------------------------------
n_lon = int((maxx - minx) / 0.009) + 1
n_lat = int((maxy - miny) / 0.009) + 1
grid_lon, grid_lat = np.meshgrid(
    np.linspace(minx, maxx, n_lon),
    np.linspace(miny, maxy, n_lat)
)

try:
    from shapely import contains_xy
    in_region = contains_xy(
        region_wgs84,
        grid_lon.ravel(),
        grid_lat.ravel()
    ).reshape(grid_lon.shape)
except ImportError:
    from shapely.vectorized import contains
    in_region = contains(
        region_wgs84,
        grid_lon.ravel(),
        grid_lat.ravel()
    ).reshape(grid_lon.shape)

# --------------------------------------------------
# Helper: reproject scattered points to regular grid
# --------------------------------------------------
def reproject(df_loc, val_col, sigma=1.0):
    grid_z = griddata(
        points = df_loc[['lon', 'lat']].values,
        values = df_loc[val_col].values,
        xi     = (grid_lon, grid_lat),
        method = 'linear'
    ).astype(np.float32)
    grid_z[~in_region] = np.nan
    smooth = gaussian_filter(
        np.where(np.isnan(grid_z), 0, grid_z), sigma=sigma
    )
    return np.where(np.isnan(grid_z), np.nan, smooth)

# --------------------------------------------------
# Degree formatters
# --------------------------------------------------
def deg_fmt_lat(x, pos): return f'{x:.1f}°'
def deg_fmt_lon(x, pos): return f'{x:.1f}°'

# --------------------------------------------------
# Plot settings
# --------------------------------------------------
plt.rcParams.update({
    'font.family'    : 'Arial',
    'font.size'      : 16,
    'axes.titlesize' : 20,
    'axes.labelsize' : 18,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
})

# --------------------------------------------------
# Helper: make and save one map
# --------------------------------------------------
def make_map(grid_data, cmap, title, cbar_label,
             out_path, gdf, extent):

    fig, ax = plt.subplots(figsize=(8, 10), dpi=600)

    vmin = np.nanpercentile(grid_data, 2)
    vmax = np.nanpercentile(grid_data, 98)

    im = ax.imshow(
        grid_data,
        origin='lower', cmap=cmap,
        vmin=vmin, vmax=vmax,
        extent=extent, interpolation='bilinear'
    )
    gdf.boundary.plot(ax=ax, edgecolor='black', linewidth=0.8)

    ax.set_title(title, fontweight='bold', fontsize=18, pad=14)
    ax.set_xlabel('Longitude', fontweight='bold', fontsize=18)
    ax.set_ylabel('Latitude',  fontweight='bold', fontsize=18)
    ax.set_aspect('equal')

    ax.xaxis.set_major_formatter(FuncFormatter(deg_fmt_lon))
    ax.yaxis.set_major_formatter(FuncFormatter(deg_fmt_lat))
    ax.tick_params(labelsize=14)
    for tick in ax.get_xticklabels() + ax.get_yticklabels():
        tick.set_fontweight('bold')

    for spine in ax.spines.values():
        spine.set_linewidth(1.2)

    cbar = plt.colorbar(im, ax=ax,
                        fraction=0.04, pad=0.035, shrink=0.78)
    cbar.set_label(cbar_label, fontsize=16, fontweight='bold')
    cbar.set_ticks(np.round(np.linspace(vmin, vmax, 6), 1))
    cbar.ax.tick_params(labelsize=13)
    for tick in cbar.ax.get_yticklabels():
        tick.set_fontweight('bold')

    plt.tight_layout()
    plt.savefig(out_path, dpi=600, bbox_inches='tight')
    plt.show()
    plt.close(fig)

# ==================================================
# MAIN LOOP
# ==================================================
for key, info in col_info.items():

    fpath = files[key]
    fname = os.path.basename(fpath)

    print(f"\n{'='*60}")
    print(f"Processing: {fname}")
    print(f"{'='*60}")

    print("  Loading CSV ...")
    df = pd.read_csv(fpath, parse_dates=['date'])
    df['month']  = df['date'].dt.month
    df['season'] = df['month'].map(
        {m: s for s, months in seasons.items() for m in months}
    )
    print(f"  Shape          : {df.shape}")
    print(f"  Columns        : {list(df.columns)}")
    print(f"  Unique locations: "
          f"{df[['lat','lon']].drop_duplicates().shape[0]:,}")
    print(f"  Lat range      : {df['lat'].min():.4f} → "
          f"{df['lat'].max():.4f}")
    print(f"  Lon range      : {df['lon'].min():.4f} → "
          f"{df['lon'].max():.4f}")

    # Output directory
    out_dir = os.path.join(base_dir, "Seasonal_Maps", key)
    os.makedirs(out_dir, exist_ok=True)

    for season_code in seasons:

        season_name = season_labels[season_code]
        df_s = df[df['season'] == season_code].copy()

        print(f"\n  Season: {season_code} — events: {len(df_s):,}")
        if len(df_s) == 0:
            print("    No events — skipping.")
            continue

        # One point per cell
        loc = df_s.groupby(['lat', 'lon']).agg(
            val    = (info['val_col'],    info['val_agg']),
            excess = (info['excess_col'], info['excess_agg'])
        ).reset_index()

        print(f"    Cells  : {len(loc):,}")
        print(f"    val    : {loc['val'].min():.2f} – "
              f"{loc['val'].max():.2f}")
        print(f"    excess : {loc['excess'].min():.2f} – "
              f"{loc['excess'].max():.2f}")

        # Reproject to regular grid
        grid_val    = reproject(loc, 'val')
        grid_excess = reproject(loc, 'excess')

        thresh_phrase = ('Below Threshold'
                         if info['val_agg'] == 'min'
                         else 'Above Threshold')

        # ── MAP 1: Main variable ───────────────────
        title1 = (f"{info['name']} {info['pct']} Percentile\n"
                  f"{season_code} ({season_name})\n"
                  f"{info['period']}")
        out1 = os.path.join(
            out_dir,
            f"{key}_{season_code}_{info['val_col']}.png"
        )
        make_map(grid_val, info['val_cmap'],
                 title1, info['val_cbar'],
                 out1, gdf, extent)
        print(f"    Saved: {os.path.basename(out1)}")

        # ── MAP 2: Excess ──────────────────────────
        title2 = (f"{info['name']} {info['pct']} Percentile "
                  f"{thresh_phrase}\n"
                  f"{season_code} ({season_name})\n"
                  f"{info['period']}")
        out2 = os.path.join(
            out_dir,
            f"{key}_{season_code}_{info['excess_col']}.png"
        )
        make_map(grid_excess, info['excess_cmap'],
                 title2, info['excess_cbar'],
                 out2, gdf, extent)
        print(f"    Saved: {os.path.basename(out2)}")

        del loc, df_s, grid_val, grid_excess

    del df

print("\n" + "=" * 60)
print("ALL SEASONAL MAPS SAVED")
print("=" * 60)