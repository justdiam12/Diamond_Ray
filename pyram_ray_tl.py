import os
import numpy as np
import matplotlib.pyplot as plt
from diamond_pyram import Diamond_PyRAM
from diamond_ray import Diamond_Ray_Code
from scipy.io import loadmat
from scipy.interpolate import interp1d


data_dir = "/Users/justindiamond/Documents/Documents/UW-APL/Research/ARMS_DATA_MASTER/LIVEOCEAN_DAILY"


files = [
    f for f in os.listdir(data_dir)
    if os.path.isfile(os.path.join(data_dir, f)) and f != "__pycache__"
]

day_21 = loadmat("/Users/justindiamond/Documents/Documents/UW-APL/Research/ARMS_DATA_MASTER/array_depths/depths_20200121.mat")
depth_21 = day_21["unique_depths"]
day_24 = loadmat("/Users/justindiamond/Documents/Documents/UW-APL/Research/ARMS_DATA_MASTER/array_depths/depths_20200124.mat")
depth_24 = day_24["unique_depths"]
day_28 = loadmat("/Users/justindiamond/Documents/Documents/UW-APL/Research/ARMS_DATA_MASTER/array_depths/depths_20200128.mat")
depth_28 = day_28["unique_depths"]
day_29 = loadmat("/Users/justindiamond/Documents/Documents/UW-APL/Research/ARMS_DATA_MASTER/array_depths/depths_20200129.mat")
depth_29 = day_29["unique_depths"]

pyram_21_dir = "/Users/justindiamond/Documents/Documents/UW-APL/Research/ARMS_DATA_MASTER/LIVEOCEAN_DAILY/2020-01-21/pyram.mat"
pyram_24_dir = "/Users/justindiamond/Documents/Documents/UW-APL/Research/ARMS_DATA_MASTER/LIVEOCEAN_DAILY/2020-01-24/pyram.mat"
pyram_28_dir = "/Users/justindiamond/Documents/Documents/UW-APL/Research/ARMS_DATA_MASTER/LIVEOCEAN_DAILY/2020-01-28/pyram.mat"
pyram_29_dir = "/Users/justindiamond/Documents/Documents/UW-APL/Research/ARMS_DATA_MASTER/LIVEOCEAN_DAILY/2020-01-29/pyram.mat"

pyram_dirs = [pyram_21_dir, pyram_24_dir, pyram_28_dir, pyram_29_dir]

SL = 195.0  # Source level (dB)

depths_list = [depth_21, depth_24, depth_28, depth_29]
dates = ["2020-01-21", "2020-01-24", "2020-01-28", "2020-01-29"]

plt.figure(figsize=(6, 10))

for date, pyram_path, unique_depths in zip(dates, pyram_dirs, depths_list):

    data = loadmat(pyram_path)

    z_pyram = np.squeeze(data["depths"])
    PL = data["pressure_level"]

    ranges = np.squeeze(data["ranges"])   # (Nr,)

    r_max = ranges.max()

    # Mask for last 0.5 meters
    r_mask = ranges >= (r_max - 1)

    # If grid is coarse and nothing selected, fallback to last index
    if np.any(r_mask):
        PL_final = np.mean(PL[:, r_mask], axis=1)   # (Nz,)
        # PL_final = PL[:, -1]
    else:
        PL_final = PL[:, -1]

    unique_depths = np.squeeze(unique_depths)

    pl_mean = []
    pl_std = []

    for d in unique_depths:
        mask = (z_pyram >= d - 1) & (z_pyram <= d + 1)

        if np.any(mask):
            vals = PL_final[mask]
        else:
            idx = np.argmin(np.abs(z_pyram - d))
            vals = np.array([PL_final[idx]])

        pl_mean.append(np.mean(vals))
        pl_std.append(np.std(vals))

    pl_mean = np.array(pl_mean)
    pl_std = np.array(pl_std)

    # Convert to Transmission Loss
    tl_mean = SL - pl_mean
    tl_std = pl_std   # unchanged

    # Sort by depth
    idx_sort = np.argsort(unique_depths)
    unique_depths = unique_depths[idx_sort]
    tl_mean = tl_mean[idx_sort]
    tl_std = tl_std[idx_sort]

    # Scatter with error bars
    plt.errorbar(
        tl_mean,
        unique_depths,
        xerr=tl_std,
        fmt='o',
        linestyle='none',
        capsize=3,
        label=date
    )

# Formatting
plt.gca().invert_yaxis()
plt.xlim(55, 75)
plt.xlabel("Transmission Loss (dB)")
plt.ylabel("Depth (m)")
plt.title("Transmission Loss vs Depth (Mean ± Std, ±1 m window)")
plt.legend()
plt.grid()
plt.savefig(os.path.join(data_dir, "tl_vs_depth.png"), dpi=300)
plt.show()