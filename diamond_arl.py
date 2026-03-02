import os
import subprocess
import arlpy.uwapm as pm
from matplotlib import pyplot as plt
from pyproj import Geod
import numpy as np
import pandas as pd
from scipy.io import loadmat
from scipy.interpolate import interp1d, RegularGridInterpolator

class Diamond_ARL():
    def __init__(self, 
                 ssp_file=None, 
                 freq=None,
                 source_level=None,
                 angle_min=None,
                 angle_max=None,
                 source_depth=None,
                 receiver_depth=None,
                 bottom_prop=None,
                 surface_prop=None,
                 water_prop=None,
                 lon_start=None,
                 lon_end=None,
                 lat_start=None,
                 lat_end=None,
                 num_points=None,
                 num_beams=None,
                 range_points=None,
                 depth_points=None,
                 bty_file=None, 
                 ati_file=None,
                 save_dir=None):
        
        if ssp_file is not None:
            self.ssp_file = ssp_file
        if bty_file is not None:    
            self.bty_file = bty_file
        if ati_file is not None:
            self.ati_file = ati_file
        if save_dir is not None:
            self.save_dir = save_dir

        # Read data
        if ssp_file is not None:
            self.ssp = self.read_ssp()

        self.freq = freq
        self.source_level = source_level
        self.angle_min = angle_min
        self.angle_max = angle_max
        self.source_depth = source_depth
        self.receiver_depth = receiver_depth

        # Bottom, Surface, and Water Properties
        if bottom_prop is not None:
            self.bottom_density = bottom_prop[0] # kg/m^3
            self.bottom_ss = bottom_prop[1]      # m/s
            self.bottom_atten = bottom_prop[2]   # 
        if surface_prop is not None:
            self.surface_density = surface_prop[0] # kg/m^3
            self.surface_ss = surface_prop[1]      # m/s
            self.surface_atten = surface_prop[2]   # 
        if water_prop is not None:
            self.water_density = water_prop[0] # kg/m^3
            self.water_atten = water_prop[1]   # dB/m kHz

        self.lon_start, self.lon_end = lon_start, lon_end
        self.lat_start, self.lat_end = lat_start, lat_end
        self.num_points = num_points
        self.num_beams = num_beams
        self.range_points = range_points
        self.depth_points = depth_points

        if bty_file is not None:
            self.bty, self.max_range, self.max_depth = self.read_bty()
        if self.angle_min is not None and self.angle_max is not None:
            self.angles = np.arange(self.angle_min, self.angle_max + 1, 1)

    def read_ssp(self):
        ssp_data = loadmat(self.ssp_file)

        ssp = ssp_data['ssp']
        ssp_depths = ssp_data['ssp_depths'].flatten()

        # -------------------------
        # Single profile
        # -------------------------
        if ssp.ndim == 1 or ssp.shape[0] == 1:
            ssp = np.atleast_1d(ssp).flatten()
            return [[ssp_depths[i], ssp[i]] for i in range(len(ssp))]

        # -------------------------
        # Along-track profile
        # -------------------------
        ssp_ranges = ssp_data['ssp_ranges'].flatten() 

        # Ensure orientation = (depths, ranges)
        if ssp.shape[0] != len(ssp_depths):
            ssp = ssp.T

        ssp_df = pd.DataFrame(
            ssp,
            index=ssp_depths,
            columns=ssp_ranges
        )

        return ssp_df
    
    def read_bty(self):
        if self.bty_file is None:
            return None, None

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

        bty_ranges = np.linspace(0, total_dist_m, self.num_points)
        track_points = np.column_stack((lats, lons))
        bty_depths = interp_bty(track_points)

        bty_export = []

        for i in range(len(bty_ranges)):
            bty_export.append([bty_ranges[i], abs(bty_depths[i])])
        
        max_range = bty_ranges[-1]
        max_depth = np.nanmax(abs(bty_depths))

        return bty_export, max_range, max_depth
    
    def eigenray(self):
        rx_depths = np.array([self.receiver_depth, self.receiver_depth + 0.01], dtype=float)
        rx_ranges = np.array([self.max_range, self.max_range], dtype=float)

        env = pm.create_env2d(
            depth=self.bty,
            soundspeed=self.ssp,
            soundspeed_interp="spline",
            bottom_soundspeed=self.bottom_ss,
            bottom_density=self.bottom_density,
            bottom_absorption=self.bottom_atten,
            tx_depth=self.source_depth,
            rx_depth=rx_depths,
            rx_range=rx_ranges,
            min_angle=self.angle_min,
            max_angle=self.angle_max,
            frequency=self.freq
        )

        # bathy = [
        #     [0, 30],    # 30 m water depth at the transmitter
        #     [5000, 30],  # 20 m water depth 300 m away
        #     [10000, 30]  # 25 m water depth at 10 km
        # ]

        # ssp = [
        #     [ 0, 1540],  # 1540 m/s at the surface
        #     [10, 1530],  # 1530 m/s at 10 m depth
        #     [20, 1532],  # 1532 m/s at 20 m depth
        #     [25, 1533],  # 1533 m/s at 25 m depth
        #     [30, 1540]   # 1540 m/s at the seabed
        # ]

        # env = pm.create_env2d(
        #     depth=bathy,
        #     soundspeed=ssp,
        #     bottom_soundspeed=1450,
        #     bottom_density=1200,
        #     bottom_absorption=1.0,
        #     tx_depth=15,
        #     rx_depth=15,
        #     rx_range=10000
        # )

        # return env
    
    def plot_rays(self, rays, env=None, invert_colors=False, figsize=(10,6)):

        if rays is None or len(rays) == 0:
            print("No rays to plot.")
            return

        fig, ax = plt.subplots(figsize=figsize)

        # Sort so highest bounce rays plot first
        rays = rays.sort_values('bottom_bounces', ascending=False)

        max_bounce = np.max(np.abs(rays.bottom_bounces)) if len(rays.bottom_bounces) > 0 else 1
        if max_bounce == 0:
            max_bounce = 1

        # Determine if we should plot in km
        all_ranges = np.concatenate([row.ray[:,0] for _, row in rays.iterrows()])
        divisor = 1
        xlabel = 'Range (m)'
        if np.max(all_ranges) - np.min(all_ranges) > 10000:
            divisor = 1000
            xlabel = 'Range (km)'

        # Plot rays
        for _, row in rays.iterrows():
            intensity = abs(row.bottom_bounces) / max_bounce
            gray = 1 - intensity if invert_colors else intensity
            color = (gray, gray, gray)

            ax.plot(row.ray[:,0]/divisor,
                    row.ray[:,1],
                    color=color,
                    linewidth=1)

        # Overlay bathymetry if available
        if env is not None and 'depth' in env:
            depth = np.array(env['depth'])
            ranges = depth[:,0] / divisor
            depths = depth[:,1]
            ax.plot(ranges, depths, color='black', linewidth=2)

        ax.set_xlabel(xlabel)
        ax.set_ylabel('Depth (m)')
        ax.set_title('Ray Paths')
        ax.grid(True)
        ax.invert_yaxis()

        plt.tight_layout()
        plt.show()
    
    def coherent_tl(self):
        env = pm.create_env2d(
            depth=self.bty,
            soundspeed=self.ssp,
            bottom_soundspeed=self.bottom_ss,
            bottom_density=self.bottom_density,
            bottom_absorption=self.bottom_atten,
            tx_depth=self.source_depth,
            rx_depth=np.linspace(0, self.max_depth-0.01, self.depth_points),
            rx_range=np.linspace(0, self.max_range, self.range_points),
            min_angle=self.angle_min,
            max_angle=self.angle_max,
            frequency=self.freq,
            nbeams=self.num_beams)

        return env

        # bathy = [
        #     [0, 30],    # 30 m water depth at the transmitter
        #     [5000, 30],  # 20 m water depth 300 m away
        #     [10000, 30]  # 25 m water depth at 10 km
        # ]

        # ssp = [
        #     [ 0, 1540],  # 1540 m/s at the surface
        #     [10, 1530],  # 1530 m/s at 10 m depth
        #     [20, 1532],  # 1532 m/s at 20 m depth
        #     [25, 1533],  # 1533 m/s at 25 m depth
        #     [30, 1540]   # 1540 m/s at the seabed
        # ]

        # env = pm.create_env2d(
        #     depth=bathy,
        #     soundspeed=ssp,
        #     frequency=self.freq,
        #     bottom_soundspeed=1450,
        #     bottom_density=1200,
        #     bottom_absorption=1.0,
        #     tx_depth=15,
        #     rx_depth=15,
        #     rx_range=10000
        # )

        # env['rx_range'] = np.linspace(0, 10000, 10001)
        # env['rx_depth'] = np.linspace(0, 30, 301)

        # return env
    

    def plot_pressure_level(self, tloss, env=None):

        # Convert to dB properly
        PL = self.source_level + 20 * np.log10(np.abs(tloss.values) + 1e-12)

        ranges = env['rx_range']
        depths = env['rx_depth']

        # Convert to km if large
        if ranges.max() > 10000:
            ranges = ranges / 1000.0
            xlabel = "Range (km)"
        else:
            xlabel = "Range (m)"

        # Mask values below bathymetry
        PL_masked = PL.copy()
        if env is not None and 'depth' in env and isinstance(env['depth'], list):
            bty = np.array(env['depth'])
            bty_range = bty[:, 0]
            bty_depth = bty[:, 1]

            if xlabel == "Range (km)":
                bty_range = bty_range / 1000

            # Interpolate bathymetry to receiver ranges
            bty_interp = np.interp(ranges, bty_range, bty_depth)

            # Create mask: True where depth > bathymetry
            mask = np.zeros_like(PL_masked, dtype=bool)
            for i, r in enumerate(ranges):
                mask[:, i] = depths > bty_interp[i]

            PL_masked = np.ma.masked_array(PL_masked, mask=mask)

        # Setup extent
        extent = [
            ranges.min(),
            ranges.max(),
            depths.min(),  # top
            depths.max()   # bottom
        ]

        plt.figure(figsize=(10, 6))

        # Use a colormap and make masked values black
        cmap = plt.get_cmap('viridis').copy()
        cmap.set_bad('k')  # masked = black

        im = plt.imshow(
            PL_masked,
            extent=extent,
            aspect='auto',
            origin='upper',
            cmap=cmap
        )

        plt.colorbar(im, label="Pressure Level (dB re 1 µPa)")
        plt.xlabel(xlabel)
        plt.ylabel("Depth (m)")
        plt.title("Pressure Level")

        # Overlay bathymetry
        if env is not None and 'depth' in env and isinstance(env['depth'], list):
            plt.plot(bty_range, bty_depth, 'k', linewidth=2)

        plt.tight_layout()
        plt.show()
    
    

if __name__ == "__main__":

    # Set the Path to Include Bellhop
    models_path = "/Users/justindiamond/Documents/Documents/UW-APL/Research/Models"
    os.environ["PATH"] = f"{models_path}:{os.environ['PATH']}"

    # Data Files and Parameters
    data_dir = '/Users/justindiamond/Documents/Documents/UW-APL/Research/Diamond_Ray/data_files'
    ssp_file = os.path.join(data_dir, 'ssp.mat')
    bty_file = os.path.join(data_dir, 'bty.mat')
    ati_file = os.path.join(data_dir, 'ati.mat')
    source_level = 195 # dB re 1 μPa @ 1 m
    freq = 25000 # Hz
    angle_min, angle_max = -80, 80 # degrees
    source_depth = 33
    receiver_depth = 55
    water_prop = (1026, 0.1)  # density (kg/m^3), attenuation (dB/m kHz)
    bottom_prop = (2000, 1500, 0.5)  # density (kg/m^3), sound speed (m/s), attenuation (dB/m kHz)
    surface_prop = (1000, 350, 0.1) # density (kg/m^3), sound speed (m/s), attenuation (dB/m kHz)
    lon_start, lon_end = -122.8, -122.85
    lat_start, lat_end = 47.78, 47.71
    num_points = 50

    # Run Ray Tracing
    run = Diamond_ARL(ssp_file=ssp_file, 
                      freq=freq,
                      source_level=source_level,
                      angle_min=angle_min, 
                      angle_max=angle_max, 
                      source_depth=source_depth, 
                      receiver_depth=receiver_depth,
                      bottom_prop=bottom_prop,
                      surface_prop=surface_prop,
                      water_prop=water_prop,
                      lon_start=lon_start,
                      lon_end=lon_end,
                      lat_start=lat_start,
                      lat_end=lat_end,
                      num_points=num_points,
                      bty_file=bty_file, 
                      ati_file=ati_file,
                      save_dir=data_dir)
    
    # env = run.eigenray()
    # rays = pm.compute_eigenrays(env)
    # run.plot_rays(rays, env=env)

    env = run.coherent_tl()
    rays = pm.compute_transmission_loss(env)
    run.plot_pressure_level(rays, env=env)