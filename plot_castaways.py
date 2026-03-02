import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat, savemat
import os       

castaway_dir = '/Users/justindiamond/Documents/Documents/UW-APL/Research/CAAS_DA/CASTAWAY'
files = os.listdir(castaway_dir)
depths_dir = '/Users/justindiamond/Documents/Documents/UW-APL/Research/Diamond_Ray/array_depths'

depths_21 = loadmat(os.path.join(depths_dir, 'depths_20200121.mat'))['unique_depths']
depths_24 = loadmat(os.path.join(depths_dir, 'depths_20200124.mat'))['unique_depths']
depths_28 = loadmat(os.path.join(depths_dir, 'depths_20200128.mat'))['unique_depths']
depths_29 = loadmat(os.path.join(depths_dir, 'depths_20200129.mat')) ['unique_depths'] 

for i in range(len(files)):
    if os.path.isdir(os.path.join(castaway_dir, files[i])):
        print(f"Processing {files[i]}...")
        data_dir = os.path.join(castaway_dir, files[i])
        ssp_file = os.path.join(data_dir, 'ssp.mat')
        pyram_file = os.path.join(data_dir, 'pyram.mat')
        source_level = 195 # dB re 1 μPa @ 1 m

        ssp_data = loadmat(ssp_file)

        ssp_depths = ssp_data['ssp_depths']
        ssp = ssp_data['ssp']

        pyram_data = loadmat(pyram_file)
        pl_ranges = pyram_data['ranges']
        pl_depths = pyram_data['depths']
        pl = pyram_data['pressure_level']

         # Convert to TL
        TL = source_level - pl   # shape: (nz, nr)

        # Flatten arrays from MATLAB (they load as 2D)
        pl_depths = pl_depths.flatten()
        pl_ranges = pl_ranges.flatten()

        # ----------------------------
        # Choose correct depth set
        # ----------------------------
        folder_name = files[i]
        if "21-Jan-2020" in folder_name:
            target_depths = depths_21.flatten()

        elif "24-Jan-2020" in folder_name:
            target_depths = depths_24.flatten()

        elif "28-Jan-2020" in folder_name:
            target_depths = depths_28.flatten()

        elif "29-Jan-2020" in folder_name:
            target_depths = depths_29.flatten()

        else:
            print("No matching depth file — skipping")
            continue

        mean_tl = []
        std_tl = []
        valid_depths = []

        # --------------------------------------
        # Extract TL statistics at each depth
        # --------------------------------------
        for d in target_depths:

            # Find nearest depth index in PyRAM grid
            idx = np.argmin(np.abs(pl_depths - d))

            tl_slice = TL[idx, :]  # TL vs range at this depth

            mean_tl.append(np.nanmean(tl_slice))
            std_tl.append(np.nanstd(tl_slice))
            valid_depths.append(pl_depths[idx])

        mean_tl = np.array(mean_tl)
        std_tl = np.array(std_tl)
        valid_depths = np.array(valid_depths)

        # -----------------------------------------
        # Setup depth lookup and storage container
        # -----------------------------------------
        depth_lookup = {
            "21-Jan-2020": depths_21.flatten(),
            "24-Jan-2020": depths_24.flatten(),
            "28-Jan-2020": depths_28.flatten(),
            "29-Jan-2020": depths_29.flatten()
        }

        daily_profiles = {}

        # -----------------------------------------
        # Loop through all cast folders ONCE
        # -----------------------------------------
        for folder in files:

            full_path = os.path.join(castaway_dir, folder)

            if not os.path.isdir(full_path):
                continue

            date_key = folder.split("_")[0]

            if date_key not in depth_lookup:
                continue

            print(f"Processing {folder}")

            pyram_file = os.path.join(full_path, "pyram.mat")
            data = loadmat(pyram_file)

            pl_depths = data["depths"].flatten()
            pl = data["pressure_level"]

            source_level = 195
            TL = source_level - pl

            target_depths = depth_lookup[date_key]

            cast_profile = []

            for d in target_depths:
                idx = np.argmin(np.abs(pl_depths - d))
                tl_slice = TL[idx, :]
                cast_profile.append(np.nanmean(tl_slice))

            cast_profile = np.array(cast_profile)

            if date_key not in daily_profiles:
                daily_profiles[date_key] = []

            daily_profiles[date_key].append(cast_profile)

        # -----------------------------------------
        # PLOT AFTER LOOP IS COMPLETE
        # -----------------------------------------

        plt.figure(figsize=(7,9))

        colors = {
            "21-Jan-2020": "blue",
            "24-Jan-2020": "red",
            "28-Jan-2020": "green",
            "29-Jan-2020": "purple"
        }

        for date_key, profiles in daily_profiles.items():

            profiles = np.array(profiles)  # shape: (num_casts, num_depths)

            mean_profile = np.nanmean(profiles, axis=0)
            std_profile = np.nanstd(profiles, axis=0)

            depths_plot = depth_lookup[date_key]

            # Plot dots with vertical error bars
            plt.errorbar(
                mean_profile,
                depths_plot,
                xerr=std_profile,
                fmt='o',                 # dots only
                color=colors[date_key],
                label=date_key,
                capsize=3
            )

        plt.gca().invert_yaxis()
        plt.xlabel("Mean Transmission Loss (dB)")
        plt.ylabel("Depth (m)")
        plt.title("Daily Mean TL Profiles")
        plt.legend()
        plt.tight_layout()
        plt.show()