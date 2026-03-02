import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
from scipy.interpolate import PchipInterpolator 

# Path
file_path = "/Users/justindiamond/Documents/Documents/UW-APL/Research/ARMS_DATA_MASTER/Tidal_Cycle"

stations = {
    "Whitney Point": "whitney_point.txt",
    "Zelatched Point": "zelatched_point.txt",
    "Quilcene": "quilcene.txt"
}

# Days to plot
target_days = [
    datetime(2020, 1, 21).date(),
    datetime(2020, 1, 24).date(),
    datetime(2020, 1, 28).date(),
    datetime(2020, 1, 29).date()
]

# --------------------------------------------
# Load all station data first
# --------------------------------------------
station_data = {}

for station_name, file_name in stations.items():

    dates = []
    heights = []
    high_low = []

    with open(os.path.join(file_path, file_name), "r") as f:
        for line in f:
            if line.startswith("2020/"):
                parts = line.split()

                dt = datetime.strptime(
                    f"{parts[0]} {parts[2]} {parts[3]}",
                    "%Y/%m/%d %I:%M %p"
                )

                dates.append(dt)
                heights.append(float(parts[4]))
                high_low.append(parts[5])

    df = pd.DataFrame({
        "datetime": dates,
        "height_m": heights,
        "type": high_low
    }).sort_values("datetime")

    station_data[station_name] = df

# --------------------------------------------
# Create 1x4 subplot layout
# --------------------------------------------
fig, axes = plt.subplots(1, 4, figsize=(16,8), sharey=True)

for ax, day in zip(axes, target_days):

    day_start = pd.Timestamp(day)
    day_end = day_start + pd.Timedelta(days=1)

    for station_name, df in station_data.items():

        window_start = day_start - pd.Timedelta(hours=12)
        window_end = day_end + pd.Timedelta(hours=12)

        window_df = df[
            (df["datetime"] >= window_start) &
            (df["datetime"] <= window_end)
        ].sort_values("datetime")

        if len(window_df) < 4:
            continue

        times_all = window_df["datetime"]
        heights_all = window_df["height_m"].values

        t0 = times_all.iloc[0]

        hours_all = np.array([
            (t - t0).total_seconds() / 3600
            for t in times_all
        ])

        spline = PchipInterpolator(hours_all, heights_all)

        hours_plot = np.linspace(
            (day_start - t0).total_seconds() / 3600,
            (day_end - t0).total_seconds() / 3600,
            400
        )

        heights_smooth = spline(hours_plot)

        time_smooth = [
            t0 + pd.Timedelta(hours=h)
            for h in hours_plot
        ]

        # Plot smooth curve
        ax.plot(time_smooth, heights_smooth, label=station_name)

        # Plot extrema inside day
        day_points = window_df[
            (window_df["datetime"] >= day_start) &
            (window_df["datetime"] <= day_end)
        ]

        ax.scatter(day_points["datetime"],
                   day_points["height_m"],
                   s=20)

    ax.set_xlim(day_start, day_end)
    ax.set_title(day.strftime("%d-%b-%Y"))
    ax.set_xlabel("Time")
    ax.grid(True)
    ax.tick_params(axis='x', rotation=45)

axes[0].set_ylabel("Tide Height (m)")
axes[3].legend(loc="lower right")

# Single legend for figure
handles, labels = axes[0].get_legend_handles_labels()
# fig.legend(handles, labels, loc="upper right")

plt.suptitle("Tide Height Comparison Across Stations", fontsize=14)
plt.tight_layout()
plt.savefig(os.path.join(file_path, "all_stations_tidal_cycle.png"), dpi=300)
plt.show()