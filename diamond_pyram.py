import os
import numpy as np
from pyram.pyram.PyRAM import PyRAM
import matplotlib.pyplot as plt
from pyproj import Geod
from scipy.io import loadmat, savemat
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
        self.source_depth = source_depth
        self.receiver_depth = receiver_depth

        if bottom_prop is not None:
            self.bottom_density = bottom_prop[0]
            self.bottom_ss = bottom_prop[1]
            self.bottom_atten = bottom_prop[2]

        self.lon_start, self.lon_end = lon_start, lon_end
        self.lat_start, self.lat_end = lat_start, lat_end

        self.num_points = num_points  # IMPORTANT

        if ssp_file is not None:
            self.ssp, self.ssp_depths, self.ssp_ranges = self.read_ssp()
        if bty_file is not None:
            self.rbzb = self.read_bty()


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

# # Save CASTAWAY results in .mat file and plot the TL grid
# if __name__ == "__main__":

#     castaway_dir = '/Users/justindiamond/Documents/Documents/UW-APL/Research/CAAS_DA/CASTAWAY'
#     files = os.listdir(castaway_dir)
#     for i in range(len(files)):
#         if os.path.isdir(os.path.join(castaway_dir, files[i])):
#             print(f"Processing {files[i]}...")
#             data_dir = os.path.join(castaway_dir, files[i])
#             ssp_file = os.path.join(data_dir, 'ssp.mat')
#             bty_file = os.path.join(data_dir, 'bty.mat')
#             source_level = 195 # dB re 1 μPa @ 1 m
#             freq = 3500 # Hz
#             source_depth = 25
#             receiver_depth = 50
#             bottom_prop = (2000, 1600, 0.8)  # density (kg/m^3), sound speed (m/s), attenuation (dB/m kHz)
#             lon_start, lon_end = -122.8, -122.85
#             lat_start, lat_end = 47.78, 47.71
#             num_points = 1000

#             # Run Ray Tracing
#             pyram = Diamond_PyRAM(ssp_file=ssp_file, 
#                                 freq=freq,
#                                 source_level=source_level,
#                                 source_depth=source_depth, 
#                                 receiver_depth=receiver_depth,
#                                 bottom_prop=bottom_prop,
#                                 lon_start=lon_start,
#                                 lon_end=lon_end,
#                                 lat_start=lat_start,
#                                 lat_end=lat_end,
#                                 num_points=num_points,
#                                 bty_file=bty_file, 
#                                 save_dir=data_dir)

#             model = pyram.create_model()
#             results = model.run()

#             ranges = results["Ranges"]
#             depths = results["Depths"]
#             pl = pyram.source_level - results["TL Grid"]

#             save_path = os.path.join(data_dir, 'pyram.mat')

#             savemat(save_path, {
#                 "ranges": ranges,
#                 "depths": depths,
#                 "pressure_level": pl
#             })

#             # Convert ranges to km for plotting
#             ranges_km = ranges / 1000.0

#             plt.figure(figsize=(10, 6))

#             im = plt.imshow(
#                 pl,
#                 extent=[ranges_km.min(), ranges_km.max(),
#                         depths.max(), depths.min()],
#                 aspect='auto',
#                 vmin=115,
#                 vmax=150
#             )

#             plt.colorbar(im, label="Pressure Level (dB re 1 μPa)")
#             plt.xlabel("Range (km)")
#             plt.ylabel("Depth (m)")
#             plt.title(f"{files[i]}")

#             plt.tight_layout()

#             fig_path = os.path.join(data_dir, "pyram.png")
#             plt.savefig(fig_path, dpi=300)
#             plt.close()

#             print(f"Saved figure to {fig_path}")



# Save LIVEOCEAN results in .mat file and plot the TL grid
if __name__ == "__main__":

    liveocean_dir = '/Users/justindiamond/Documents/Documents/UW-APL/Research/ARMS_DATA_MASTER/LIVEOCEAN'
    files = os.listdir(liveocean_dir)
    for i in range(len(files)):
        if os.path.isdir(os.path.join(liveocean_dir, files[i])) and os.path.isfile(os.path.join(liveocean_dir, files[i]) + '/ssp.mat'):
            print(f"Processing {files[i]}...")
            data_dir = os.path.join(liveocean_dir, files[i])
            ssp_file = os.path.join(data_dir, 'ssp.mat')
            bty_file = os.path.join(data_dir, 'bty.mat')
            source_level = 195 # dB re 1 μPa @ 1 m
            freq = 3500 # Hz
            source_depth = 25
            receiver_depth = 50
            bottom_prop = (2000, 1600, 0.8)  # density (kg/m^3), sound speed (m/s), attenuation (dB/m kHz)
            lon_start, lon_end = -122.8, -122.85
            lat_start, lat_end = 47.78, 47.71
            num_points = 1000

            # Run Ray Tracing
            pyram = Diamond_PyRAM(ssp_file=ssp_file, 
                                freq=freq,
                                source_level=source_level,
                                source_depth=source_depth, 
                                receiver_depth=receiver_depth,
                                bottom_prop=bottom_prop,
                                lon_start=lon_start,
                                lon_end=lon_end,
                                lat_start=lat_start,
                                lat_end=lat_end,
                                num_points=num_points,
                                bty_file=bty_file, 
                                save_dir=data_dir)

            model = pyram.create_model()
            results = model.run()

            ranges = results["Ranges"]
            depths = results["Depths"]
            pl = pyram.source_level - results["TL Grid"]

            save_path = os.path.join(data_dir, 'pyram.mat')

            savemat(save_path, {
                "ranges": ranges,
                "depths": depths,
                "pressure_level": pl
            })

            # Convert ranges to km for plotting
            ranges_km = ranges / 1000.0

            plt.figure(figsize=(10, 6))

            im = plt.imshow(
                pl,
                extent=[ranges_km.min(), ranges_km.max(),
                        depths.max(), depths.min()],
                aspect='auto',
                vmin=115,
                vmax=150
            )

            plt.colorbar(im, label="Pressure Level (dB re 1 μPa)")
            plt.xlabel("Range (km)")
            plt.ylabel("Depth (m)")
            plt.title(f"{files[i]}")

            plt.tight_layout()

            fig_path = os.path.join(data_dir, "pyram.png")
            plt.savefig(fig_path, dpi=300)
            plt.close()

            print(f"Saved figure to {fig_path}")