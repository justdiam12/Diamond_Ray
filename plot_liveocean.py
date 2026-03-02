import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.io import loadmat
from pyproj import Geod
from scipy.interpolate import RegularGridInterpolator, PchipInterpolator
from matplotlib.animation import FFMpegWriter
from datetime import datetime
import matplotlib.dates as mdates

# Read BTY
def read_bty(bty_file, lat_start, lat_end, lon_start, lon_end, num_points):

    bty_data = loadmat(bty_file)
    bty_map = bty_data['bath_map']
    lat_range = bty_data['lat_range'].flatten()
    lon_range = bty_data['lon_range'].flatten()

    interp_bty = RegularGridInterpolator(
        (lat_range, lon_range),
        bty_map,
        bounds_error=False,
        fill_value=np.nan
    )

    g = Geod(ellps="WGS84")
    _, _, total_dist_m = g.inv(lon_start, lat_start, lon_end, lat_end)
    pts = g.npts(lon_start, lat_start, lon_end, lat_end, num_points - 2)

    lons = np.array([lon_start] + [p[0] for p in pts] + [lon_end])
    lats = np.array([lat_start] + [p[1] for p in pts] + [lat_end])

    bty_ranges = np.linspace(0, total_dist_m, num_points)
    track_points = np.column_stack((lats, lons))
    bty_depths = interp_bty(track_points)

    return bty_ranges, bty_depths


# ==========================================================
# PATHS
# ==========================================================
liveocean_dir = '/Users/justindiamond/Documents/Documents/UW-APL/Research/ARMS_DATA_MASTER/LIVEOCEAN'
tide_file = '/Users/justindiamond/Documents/Documents/UW-APL/Research/ARMS_DATA_MASTER/Tidal_Cycle/whitney_point.txt'
bty_file = '/Users/justindiamond/Documents/Documents/UW-APL/Research/Diamond_Ray/data_files/bty.mat'

lon_start, lon_end = -122.8, -122.85
lat_start, lat_end = 47.78, 47.71
num_points = 1000

# ==========================================================
# LOAD BATHYMETRY
# ==========================================================
bty_ranges, bty_depths = read_bty(
    bty_file, lat_start, lat_end,
    lon_start, lon_end, num_points
)

bty_ranges_km = bty_ranges / 1000.0
bty_depths_plot = np.abs(bty_depths)

# ==========================================================
# SORT FOLDERS
# ==========================================================
folders = sorted([
    f for f in os.listdir(liveocean_dir)
    if os.path.isdir(os.path.join(liveocean_dir, f))
])

# folders = folders[0:20]

# ==========================================================
# LOAD FULL TIDE DATA (NO WINDOW CLIP)
# ==========================================================
tide_times = []
tide_heights = []

with open(tide_file, "r") as f:
    for line in f:
        if line.startswith("2020/"):
            parts = line.split()
            dt = datetime.strptime(
                f"{parts[0]} {parts[2]} {parts[3]}",
                "%Y/%m/%d %I:%M %p"
            )
            tide_times.append(dt)
            tide_heights.append(float(parts[4]) / 3.281)

tide_heights = tide_heights

tide_df = pd.DataFrame({
    "datetime": tide_times,
    "height": tide_heights
}).sort_values("datetime")

# Build spline using FULL dataset
t0_full = tide_df["datetime"].iloc[0]

t_hours_full = np.array([
    (t - t0_full).total_seconds() / 3600
    for t in tide_df["datetime"]
])

tide_spline = PchipInterpolator(t_hours_full, tide_df["height"].values)

# ==========================================================
# DEFINE DISPLAY WINDOW
# ==========================================================
tide_start = datetime(2020, 1, 21, 0)
tide_end   = datetime(2020, 1, 30, 0)

display_hours = np.linspace(
    (tide_start - t0_full).total_seconds()/3600,
    (tide_end - t0_full).total_seconds()/3600,
    3000
)

tide_smooth = tide_spline(display_hours)

tide_time_smooth = [
    t0_full + pd.Timedelta(hours=h)
    for h in display_hours
]

# ==========================================================
# FIGURE SETUP
# ==========================================================
fig, (ax1, ax2, ax3) = plt.subplots(
    3, 1, figsize=(14, 11), constrained_layout=True
)

writer = FFMpegWriter(fps=6)
video_path = os.path.join(liveocean_dir, "liveocean_PL_evolution.mp4")

# ==========================================================
# INITIALIZE PANELS
# ==========================================================
im1 = ax1.imshow(np.zeros((10,10)), aspect='auto', vmin=115, vmax=150)
fig.colorbar(im1, ax=ax1).set_label("Pressure Level (dB re 1 μPa)")

im2 = ax2.imshow(np.zeros((10,10)), aspect='auto', vmin=1472, vmax=1485)
fig.colorbar(im2, ax=ax2).set_label("Sound Speed (m/s)")
im2.cmap.set_bad(color='black')

global_levels = np.arange(1472, 1486, 1.5)
contours = None
bty_line = None

tide_line, = ax3.plot([], [], linewidth=2)
tide_dot, = ax3.plot([], [], 'ro')

ax3.set_xlim(tide_start, tide_end)
ax3.set_ylim(min(tide_smooth)-0.2, max(tide_smooth)+0.2)
ax3.set_ylabel("Tide Height (m)")
ax3.set_xlabel("Time")
ax3.grid(True)
ax3.xaxis.set_major_locator(mdates.DayLocator())
ax3.xaxis.set_major_formatter(mdates.DateFormatter('%b %d'))

# ==========================================================
# ANIMATION LOOP
# ==========================================================
with writer.saving(fig, video_path, dpi=100):

    for folder in folders:

        data_dir = os.path.join(liveocean_dir, folder)
        pyram_file = os.path.join(data_dir, "pyram.mat")
        ssp_file = os.path.join(data_dir, "ssp.mat")

        if not (os.path.isfile(pyram_file) and os.path.isfile(ssp_file)):
            continue

        print(f"Adding frame: {folder}")
        current_time = datetime.strptime(folder, "%Y-%m-%dT%H")

        # ---------------- PL ----------------
        pyram_data = loadmat(pyram_file)
        ranges = pyram_data["ranges"].flatten()
        depths = pyram_data["depths"].flatten()
        pl = pyram_data["pressure_level"]

        ranges_km = ranges / 1000.0

        im1.set_data(pl)
        im1.set_extent([
            ranges_km.min(),
            ranges_km.max(),
            depths.max(),
            depths.min()
        ])

        # ---------------- SSP ----------------
        ssp_data = loadmat(ssp_file)
        ssp_depths = ssp_data["ssp_depths"].flatten()
        ssp = ssp_data["ssp"]
        ssp_ranges = ssp_data["ssp_ranges"].flatten()

        ssp_interp = RegularGridInterpolator(
            (ssp_depths, ssp_ranges),
            ssp,
            bounds_error=False,
            fill_value=np.nan
        )

        R_high, D_high = np.meshgrid(bty_ranges_km, ssp_depths)
        points = np.column_stack((D_high.ravel(), R_high.ravel()))
        ssp_high = ssp_interp(points).reshape(D_high.shape)

        mask = D_high > bty_depths_plot[np.newaxis, :]
        ssp_masked = np.ma.masked_where(mask, ssp_high)

        im2.set_data(ssp_masked)
        im2.set_extent([
            bty_ranges_km.min(),
            bty_ranges_km.max(),
            ssp_depths.max(),
            ssp_depths.min()
        ])

        ax2.set_xlim(bty_ranges_km.min(), bty_ranges_km.max())
        ax2.set_ylim(ssp_depths.max(), ssp_depths.min())
        ax2.margins(x=0)

        if contours is not None:
            for coll in contours.collections:
                coll.remove()

        contours = ax2.contour(
            R_high, D_high, ssp_masked,
            levels=global_levels,
            colors='k',
            linewidths=1.8
        )

        if bty_line is not None:
            bty_line.remove()

        bty_line, = ax2.plot(
            bty_ranges_km,
            bty_depths_plot,
            color='k',
            linewidth=2.5
        )

        # ---------------- TIDE ----------------
        idx = np.searchsorted(tide_time_smooth, current_time)

        tide_line.set_data(
            tide_time_smooth[:idx],
            tide_smooth[:idx]
        )

        if idx > 0:
            tide_dot.set_data(
                tide_time_smooth[idx-1],
                tide_smooth[idx-1]
            )

        fig.suptitle(
            "LIVEOCEAN Pressure Level, Sound Speed, and Tide\n" + folder
        )

        writer.grab_frame()

print(f"\nVideo saved to:\n{video_path}")