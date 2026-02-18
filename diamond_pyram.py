import os
import numpy as np
from pyram.pyram.PyRAM import PyRAM
import matplotlib.pyplot as plt
from pyproj import Geod
from scipy.io import loadmat
from scipy.interpolate import RegularGridInterpolator



class Diamond_PyRAM():

    def __init__(self,
                 ssp_file,
                 freq,
                 source_level,
                 source_depth,
                 receiver_depth,
                 bottom_prop,
                 lon_start,
                 lon_end,
                 lat_start,
                 lat_end,
                 num_points,
                 bty_file,
                 save_dir):

        self.ssp_file = ssp_file
        self.bty_file = bty_file
        self.save_dir = save_dir

        self.freq = freq
        self.source_level = source_level
        self.source_depth = source_depth
        self.receiver_depth = receiver_depth

        self.bottom_density = bottom_prop[0]
        self.bottom_ss = bottom_prop[1]
        self.bottom_atten = bottom_prop[2]

        self.lon_start, self.lon_end = lon_start, lon_end
        self.lat_start, self.lat_end = lat_start, lat_end

        self.num_points = num_points  # IMPORTANT

        self.ssp, self.ssp_depths, self.ssp_ranges = self.read_ssp()
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
    
    def create_model(self):
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

# # -----------------------------
# # Basic environment definition
# # -----------------------------

# freq = 3500.0        # Hz
# zs = 33.0            # source depth (m)
# zr = 50.0            # receiver depth (m)
# water_depth = 100.0  # total water depth (m)
# rmax = 8500.0        # 8.5 km

# # -----------------------------
# # Water column (constant SSP)
# # -----------------------------

# z_ss = np.linspace(0, water_depth, 101)      # water depths
# rp_ss = np.array([0.0])                      # range-independent
# cw = 1500.0 * np.ones((z_ss.size, 1))        # 1500 m/s everywhere

# # -----------------------------
# # Bottom properties (halfspace)
# # -----------------------------

# z_sb = np.array([0.0])                       # bottom layer thickness reference
# rp_sb = np.array([0.0])                      # range-independent

# cb = np.array([[1700.0]])                    # bottom sound speed (m/s)
# rhob = np.array([[1.8]])                     # density (g/cm^3)
# attn = np.array([[0.5]])                     # attenuation (dB/wavelength)

# # -----------------------------
# # Bathymetry (flat bottom)
# # -----------------------------

# rbzb = np.array([
#     [0.0, water_depth],
#     [rmax, water_depth]
# ])

# # -----------------------------
# # Create model
# # -----------------------------

# model = PyRAM(
#     freq, zs, zr,
#     z_ss, rp_ss, cw,
#     z_sb, rp_sb, cb, rhob,
#     attn, rbzb,
#     rmax=rmax
# )

# # -----------------------------
# # Run model
# # -----------------------------

# results = model.run()

# ranges = results["Ranges"]
# depths = results["Depths"]
# tl_grid = results["TL Grid"]

# # -----------------------------
# # Plot TL
# # -----------------------------

# plt.figure(figsize=(10,6))
# plt.imshow(
#     tl_grid,
#     extent=[ranges.min()/1000, ranges.max()/1000,
#             depths.max(), depths.min()],
#     aspect='auto',
#     vmin=0,
#     vmax=75
# )
# plt.colorbar(label="Transmission Loss (dB)")
# plt.xlabel("Range (km)")
# plt.ylabel("Depth (m)")
# plt.title("PyRAM TL (Flat 100m Waveguide)")
# plt.tight_layout()
# plt.show()

if __name__ == "__main__":
    data_dir = '/Users/justindiamond/Documents/Documents/UW-APL/Research/Diamond_Ray/data_files'
    ssp_file = os.path.join(data_dir, 'ssp.mat')
    bty_file = os.path.join(data_dir, 'bty.mat')
    source_level = 195 # dB re 1 μPa @ 1 m
    freq = 3500 # Hz
    source_depth = 40
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
    tl_grid = results["TL Grid"]

    # -----------------------------
    # Plot Pressure Level
    # -----------------------------

    plt.figure(figsize=(10,6))
    plt.imshow(
        195 -tl_grid,
        extent=[ranges.min()/1000, ranges.max()/1000,
                depths.max(), depths.min()],
        aspect='auto',
        vmin=110,
        vmax=150
    )
    plt.colorbar(label="Pressure Level (dB re 1 μPa)")
    plt.xlabel("Range (km)")
    plt.ylabel("Depth (m)")
    plt.title("Pressure Level Using ARMS SSP and Bathymetry")
    plt.tight_layout()
    plt.show()