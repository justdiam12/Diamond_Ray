import os
import numpy as np
from pyram.pyram.PyRAM import PyRAM
import matplotlib.pyplot as plt
from pyproj import Geod
from scipy.io import loadmat, savemat
from datetime import datetime
from scipy.interpolate import RegularGridInterpolator


class Diamond_PyRAM():

    def __init__(self,
                 ssp_file=None,
                 freq=None,
                 source_level=None,
                 source_depth=None,
                 receiver_depth=None,
                 bottom_prop=None,
                 lon_start=None,
                 lon_end=None,
                 lat_start=None,
                 lat_end=None,
                 num_points=None,
                 bty_file=None,
                 save_dir=None):

        if ssp_file is not None:
            self.ssp_file = ssp_file
        if bty_file is not None:    
            self.bty_file = bty_file
        if save_dir is not None:
            self.save_dir = save_dir


        self.freq = freq
        self.source_level = source_level
        self.receiver_depth = receiver_depth

        if bottom_prop is not None:
            self.bottom_density = bottom_prop[0]
            self.bottom_ss = bottom_prop[1]
            self.bottom_atten = bottom_prop[2]

        self.lon_start, self.lon_end = lon_start, lon_end
        self.lat_start, self.lat_end = lat_start, lat_end

        self.num_points = num_points

        if ssp_file is not None:
            self.ssp, self.ssp_depths, self.ssp_ranges = self.read_ssp()
        if bty_file is not None:
            self.rbzb = self.read_bty()

        if source_depth is not None:
            self.source_depth = source_depth
        # else:
        #     self.source_depth = self.rbzb[0,1] - 4


    def read_ssp(self):
        ssp_data = loadmat(self.ssp_file)

        ssp = ssp_data['ssp']
        ssp_depths = ssp_data['ssp_depths'].flatten()

        if ssp.ndim == 1 or ssp.shape[0] == 1:
            ssp = ssp.flatten()
            ssp_ranges = None
        else:
            ssp_ranges = ssp_data['ssp_ranges'].flatten() * 1000.0  # km → m

        return ssp, ssp_depths, ssp_ranges

    def read_bty(self):

        if self.bty_file is None:
            return None

        bty_data = loadmat(self.bty_file)
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
        _, _, total_dist_m = g.inv(
            self.lon_start, self.lat_start,
            self.lon_end, self.lat_end
        )

        pts = g.npts(
            self.lon_start, self.lat_start,
            self.lon_end, self.lat_end,
            self.num_points - 2
        )

        lons = np.array([self.lon_start] + [p[0] for p in pts] + [self.lon_end])
        lats = np.array([self.lat_start] + [p[1] for p in pts] + [self.lat_end])

        ranges = np.linspace(0, total_dist_m, self.num_points)
        track_points = np.column_stack((lats, lons))
        depths = np.abs(interp_bty(track_points))

        rbzb = np.column_stack((ranges, depths))

        return rbzb
    
    def create_model(self, dr=None):
        
        if self.ssp.ndim == 1:
            z_ss = self.ssp_depths
            rp_ss = np.array([0.0])
            cw = self.ssp.reshape(-1, 1)

            z_sb = np.array([0.0])
            rp_sb = np.array([0.0])

            cb = np.array([[self.bottom_ss]])
            rhob = np.array([[self.bottom_density]])
            attn = np.array([[self.bottom_atten]])

            rmax = self.rbzb[-1, 0]

            model = PyRAM(
                self.freq,
                self.source_depth,
                self.receiver_depth,
                z_ss,
                rp_ss,
                cw,
                z_sb,
                rp_sb,
                cb,
                rhob,
                attn,
                self.rbzb,
                rmax=rmax
            )


            return model
        else: 
            z_ss = self.ssp_depths
            rp_ss = self.ssp_ranges
            cw = self.ssp

            z_sb = np.array([0.0])
            rp_sb = np.array([0.0])

            cb = np.array([[self.bottom_ss]])
            rhob = np.array([[self.bottom_density]])
            attn = np.array([[self.bottom_atten]])

            rmax = self.rbzb[-1, 0]

            model = PyRAM(
                self.freq,
                self.source_depth,
                self.receiver_depth,
                z_ss,
                rp_ss,
                cw,
                z_sb,
                rp_sb,
                cb,
                rhob,
                attn,
                self.rbzb,
                rmax=rmax
            )

            return model


# # Save LIVEOCEAN results in .mat file and plot the TL grid
# if __name__ == "__main__":

#     liveocean_dir = '/Users/justindiamond/Documents/Documents/UW-APL/Research/ARMS_DATA_MASTER/LIVEOCEAN_DAILY'
#     files = os.listdir(liveocean_dir)

#     for i in range(len(files)):

#         data_dir_full = os.path.join(liveocean_dir, files[i])

#         if os.path.isdir(data_dir_full) and os.path.isfile(os.path.join(data_dir_full, 'ssp.mat')):

#             print(f"Processing {files[i]}...")

#             ssp_file = os.path.join(data_dir_full, 'ssp.mat')
#             bty_file = os.path.join(data_dir_full, 'bty.mat')

#             source_level = 195  # dB re 1 μPa @ 1 m
#             freq = 3500  # Hz
#             # source_depth = 35
#             receiver_depth = 50

#             bottom_prop = (2000, 1600, 0.8)

#             transect_coords = {
#                 "2020-01-21": {"start": (47.7729, -122.80275), "end": (47.73, -122.84)},
#                 "2020-01-22": {"start": (47.7729, -122.80275), "end": (47.73, -122.84)},
#                 "2020-01-23": {"start": (47.77295, -122.802983), "end": (47.73, -122.84)},
#                 "2020-01-24": {"start": (47.77295, -122.802983),     "end": (47.73, -122.84)},
#                 "2020-01-25": {"start": (47.772967, -122.802917), "end": (47.73, -122.84)},
#                 "2020-01-26": {"start": (47.772967, -122.802917), "end": (47.73, -122.84)},
#                 "2020-01-27": {"start": (47.773000, -122.802717), "end": (47.73, -122.84)},
#                 "2020-01-28": {"start": (47.773000, -122.802717),     "end": (47.73, -122.84)},
#                 "2020-01-29": {"start": (47.773000, -122.802717),     "end": (47.73, -122.84)},
#                 "2020-01-30": {"start": (47.773000, -122.802717), "end": (47.73, -122.84)},
#             }
            
#             # --------------------------
#             # Get transect for this day
#             # --------------------------
#             date_str = files[i]  # folder name like "2020-01-21"

#             if date_str not in transect_coords:
#                 print(f"Skipping {date_str} (no transect defined)")
#                 continue

#             (start_lat, start_lon) = transect_coords[date_str]["start"]
#             (end_lat, end_lon)     = transect_coords[date_str]["end"]

#             lat_start, lon_start = start_lat, start_lon
#             lat_end, lon_end     = end_lat, end_lon

#             num_points = 1000

#             # --- Run PyRAM ---
#             pyram = Diamond_PyRAM(
#                 ssp_file=ssp_file,
#                 freq=freq,
#                 source_level=source_level,
#                 source_depth=None,
#                 receiver_depth=receiver_depth,
#                 bottom_prop=bottom_prop,
#                 lon_start=lon_start,
#                 lon_end=lon_end,
#                 lat_start=lat_start,
#                 lat_end=lat_end,
#                 num_points=num_points,
#                 bty_file=bty_file,
#                 save_dir=data_dir_full
#             )

#             model = pyram.create_model()
#             results = model.run()

#             ranges = results["Ranges"]              # (Nr,)
#             depths = results["Depths"]              # (Nz,)
#             TL = results["TL Grid"]                 # (Nz, Nr)

#             pl = pyram.source_level - TL            # convert to pressure level

#             # --- Save pyram output ---
#             save_path = os.path.join(data_dir_full, 'pyram.mat')
#             savemat(save_path, {
#                 "ranges": ranges,
#                 "depths": depths,
#                 "pressure_level": pl
#             })

#             # --- Convert to km ---
#             ranges_km = ranges / 1000.0

#             # --- SSP from pyram object ---
#             z_ssp = pyram.ssp_depths
#             r_ssp = pyram.ssp_ranges
#             c = pyram.ssp

#             # Ensure km consistency
#             r_ssp_km = r_ssp / 1000.0 if np.max(r_ssp) > 100 else r_ssp

#             # --- Bathymetry ---
#             rbzb = pyram.rbzb
#             bty_ranges = rbzb[:, 0]
#             bty_depths = rbzb[:, 1]

#             bty_ranges_km = bty_ranges / 1000.0 if np.max(bty_ranges) > 100 else bty_ranges

#             # --- Figure ---
#             fig, (ax1, ax2) = plt.subplots(
#                 2, 1,
#                 figsize=(12, 8),
#                 sharex=True,
#                 gridspec_kw={"height_ratios": [1, 1]}
#             )

#             # =========================
#             # Top: Pressure Level
#             # =========================
#             im1 = ax1.imshow(
#                 pl,
#                 extent=[ranges_km.min(), ranges_km.max(),
#                         depths.max(), depths.min()],
#                 aspect='auto',
#                 vmin=115,
#                 vmax=150
#             )

#             # Bathymetry overlay
#             ax1.plot(bty_ranges_km, bty_depths, color='black', linewidth=2)

#             # Mask below seafloor
#             ax1.fill_between(
#                 bty_ranges_km,
#                 bty_depths,
#                 depths.max(),
#                 color='black'
#             )

#             ax1.set_ylabel("Depth (m)")
#             ax1.set_title(f"{files[i]} - Pressure Level")
#             fig.colorbar(im1, ax=ax1, label="PL (dB re 1 μPa)")

#             # =========================
#             # Bottom: SSP
#             # =========================
#             im2 = ax2.imshow(
#                 c,
#                 extent=[r_ssp_km.min(), r_ssp_km.max(),
#                         z_ssp.max(), z_ssp.min()],
#                 aspect='auto',
#                 vmin=1470,
#                 vmax=1488
#             )

#             # Bathymetry overlay
#             ax2.plot(bty_ranges_km, bty_depths, color='black', linewidth=2)

#             # Mask below seafloor
#             ax2.fill_between(
#                 bty_ranges_km,
#                 bty_depths,
#                 z_ssp.max(),
#                 color='black'
#             )

#             ax2.set_xlabel("Range (km)")
#             ax2.set_ylabel("Depth (m)")
#             ax2.set_title("Sound Speed Profile (m/s)")
#             fig.colorbar(im2, ax=ax2, label="Sound Speed (m/s)")

#             # --- Layout ---
#             plt.tight_layout()

#             # --- Save figure ---
#             fig_path = os.path.join(data_dir_full, "pyram_ssp.png")
#             plt.savefig(fig_path, dpi=300)
#             plt.close()

#             print(f"Saved figure to {fig_path}")


# Run PyRAM with Castaways
if __name__ == "__main__":

    castaway_dir = '/Users/justindiamond/Documents/Documents/UW-APL/Research/ARMS_DATA_MASTER/CASTAWAY'
    files = sorted(os.listdir(castaway_dir))

    # Transects by day
    transect_coords = {
        "2020-01-21": {"start": (47.7729, -122.80275), "end": (47.73, -122.84)},
        "2020-01-22": {"start": (47.7729, -122.80275), "end": (47.73, -122.84)},
        "2020-01-23": {"start": (47.77295, -122.802983), "end": (47.73, -122.84)},
        "2020-01-24": {"start": (47.77295, -122.802983), "end": (47.73, -122.84)},
        "2020-01-25": {"start": (47.772967, -122.802917), "end": (47.73, -122.84)},
        "2020-01-26": {"start": (47.772967, -122.802917), "end": (47.73, -122.84)},
        "2020-01-27": {"start": (47.773000, -122.802717), "end": (47.73, -122.84)},
        "2020-01-28": {"start": (47.773000, -122.802717), "end": (47.73, -122.84)},
        "2020-01-29": {"start": (47.773000, -122.802717), "end": (47.73, -122.84)},
        "2020-01-30": {"start": (47.773000, -122.802717), "end": (47.73, -122.84)},
    }

    # Loop over CASTAWAY folders
    for folder in files:

        data_dir_full = os.path.join(castaway_dir, folder)

        if not (os.path.isdir(data_dir_full) and 
                os.path.isfile(os.path.join(data_dir_full, 'ssp.mat'))):
            continue

        print(f"Processing {folder}...")

        # CASTAWAY date
        try:
            date_part = folder.split('_')[0]
            date_obj = datetime.strptime(date_part, "%d-%b-%Y")
            date_str = date_obj.strftime("%Y-%m-%d")

        except ValueError:
            print(f"Skipping {folder} (invalid date format)")
            continue

        # Get transect
        if date_str not in transect_coords:
            print(f"Skipping {folder} (no transect defined)")
            continue

        (start_lat, start_lon) = transect_coords[date_str]["start"]
        (end_lat, end_lon)     = transect_coords[date_str]["end"]

        lat_start, lon_start = start_lat, start_lon
        lat_end, lon_end     = end_lat, end_lon

        # Files + params
        ssp_file = os.path.join(data_dir_full, 'ssp.mat')
        bty_file = os.path.join(data_dir_full, 'bty.mat')

        source_level = 195
        freq = 3500
        receiver_depth = 50

        bottom_prop = (2000, 1600, 0.8)
        num_points = 1000

        # Run PyRAM
        pyram = Diamond_PyRAM(
            ssp_file=ssp_file,
            freq=freq,
            source_level=source_level,
            source_depth=None,
            receiver_depth=receiver_depth,
            bottom_prop=bottom_prop,
            lon_start=lon_start,
            lon_end=lon_end,
            lat_start=lat_start,
            lat_end=lat_end,
            num_points=num_points,
            bty_file=bty_file,
            save_dir=data_dir_full
        )

        model = pyram.create_model()
        results = model.run()

        ranges = results["Ranges"]
        depths = results["Depths"]
        TL = results["TL Grid"]

        pl = pyram.source_level - TL

        # Save results
        savemat(os.path.join(data_dir_full, 'pyram.mat'), {
            "ranges": ranges,
            "depths": depths,
            "pressure_level": pl
        })

        # Prep plotting
        ranges_km = ranges / 1000.0

        z_ssp = pyram.ssp_depths
        r_ssp = pyram.ssp_ranges
        c = pyram.ssp

        r_ssp_km = r_ssp / 1000.0 if r_ssp is not None and np.max(r_ssp) > 100 else r_ssp

        rbzb = pyram.rbzb
        bty_ranges_km = rbzb[:, 0] / 1000.0
        bty_depths = rbzb[:, 1]

        # Plot
        fig = plt.figure(figsize=(10, 4))
        gs = fig.add_gridspec(1, 2, width_ratios=[3, 1])

        ax1 = fig.add_subplot(gs[0])  # TL
        ax2 = fig.add_subplot(gs[1])  # SSP profile

        # PL
        im1 = ax1.imshow(
            pl,
            extent=[ranges_km.min(), ranges_km.max(),
                    depths.max(), depths.min()],
            aspect='auto',
            vmin=115,
            vmax=150
        )

        # BTY
        ax1.plot(bty_ranges_km, bty_depths, 'k', lw=2)
        ax1.fill_between(bty_ranges_km, bty_depths, depths.max(), color='black')

        ax1.set_xlabel("Range (km)")
        ax1.set_ylabel("Depth (m)")
        ax1.set_title(f"{folder} - Pressure Level")
        fig.colorbar(im1, ax=ax1, label="PL (dB re 1 μPa)")

        # SSP
        # If SSP is 2D but constant in range, just take first column
        if c.ndim == 2:
            c_profile = c[:, 0]
        else:
            c_profile = c

        ax2.plot(c_profile, z_ssp, 'b', lw=2)

        ax2.set_ylim(z_ssp.max(), z_ssp.min())
        ax2.set_xlabel("Sound Speed (m/s)")
        ax2.set_ylabel("Depth (m)")
        ax2.set_title("SSP")
        ax2.grid(True)

        plt.tight_layout()

        # Save
        fig_path = os.path.join(data_dir_full, "pyram_ssp.png")
        plt.savefig(fig_path, dpi=300)
        plt.close()

        print(f"Saved figure to {fig_path}")


