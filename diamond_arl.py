import os
import subprocess
import arlpy.uwapm as pm
from matplotlib import pyplot as plt
from pyproj import Geod
import numpy as np
import pandas as pd
from scipy.io import loadmat
from scipy.interpolate import interp1d, RegularGridInterpolator
from transmit_sig import *

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

        # Range-Independent SSP
        if ssp.ndim == 1 or ssp.shape[0] == 1:
            ssp = np.atleast_1d(ssp).flatten()
            return [[ssp_depths[i], ssp[i]] for i in range(len(ssp))]

        # Range-Dependent SSP
        ssp_ranges = ssp_data['ssp_ranges'].flatten() 
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
    
    
    def compute_s_t(self, ray):
        ssp = self.ssp
        ssp_z = []
        ssp_c = []
        for i in range(len(ssp)):
            ssp_i = ssp[i]
            ssp_z.append(ssp_i[0])
            ssp_c.append(ssp_i[1])
        ssp_func = interp1d(ssp_z, ssp_c, fill_value="extrapolate")

        r = ray[:, 0]
        z = ray[:, 1]

        total_time = 0.0
        total_dist = 0.0

        # Loop through segments
        for i in range(len(r) - 1):
            dr = r[i+1] - r[i]
            dz = z[i+1] - z[i]

            ds = np.sqrt(dr**2 + dz**2)
            z_mid = 0.5 * (z[i] + z[i+1])

            c = float(ssp_func(z_mid))

            total_time += ds / c
            total_dist += ds

        # Attach to DataFrame
        return total_time, total_dist


    def get_bounce_sequence(self, ray, env, tol=1.0):
        r = ray[:, 0]
        z = ray[:, 1]

        # Bathymetry
        bty = self.bty
        bty_z = []
        bty_r = []
        for i in range(len(bty)):
            bty_i = bty[i]
            bty_r.append(bty_i[0])
            bty_z.append(bty_i[1])
        bathy_func = interp1d(bty_r, bty_z, fill_value="extrapolate")
        z_bottom = bathy_func(r)

        sequence = []
        reflection_coeffs = []

        in_surface = False
        in_bottom = False

        # Determine Reflection Coefficient
        rho_w = self.water_density
        c_w   = self.ssp[0][1]  
        rho_b = self.bottom_density
        c_b   = self.bottom_ss 
        Z_w = rho_w * c_w
        Z_b = rho_b * c_b
        R_bottom = (Z_b - Z_w) / (Z_b + Z_w)

        for i in range(len(z)):

            # Surface bounce (pressure release)
            if abs(z[i]) <= tol:
                if not in_surface:
                    sequence.append(0)
                    reflection_coeffs.append(-1.0)  # phase inversion
                    in_surface = True
            else:
                in_surface = False

            # Bottom bounce
            if abs(z[i] - z_bottom[i]) <= tol:
                if not in_bottom:
                    sequence.append(1)
                    reflection_coeffs.append(R_bottom)
                    in_bottom = True
            else:
                in_bottom = False

        return sequence, reflection_coeffs


    def eigenray(self):

        env = pm.create_env2d(
            depth=self.bty,
            soundspeed=self.ssp,
            soundspeed_interp="spline",
            bottom_soundspeed=self.bottom_ss,
            bottom_density=self.bottom_density,
            bottom_absorption=self.bottom_atten,
            tx_depth=self.source_depth,
            rx_depth=self.receiver_depth,
            rx_range=self.max_range,
            min_angle=self.angle_min,
            max_angle=self.angle_max,
            frequency=self.freq
        )
        rays = pm.compute_eigenrays(env)
        return env, rays
    

    def unique_eigenrays(self, rays, env, num_rays=10, tol=1):

        if rays is None or len(rays) == 0:
            print("No rays to plot.")
            return

        selected_ray_idxs = []
        
        for index, (_, row) in enumerate(rays.iterrows()):
            r = row.ray[:, 0]
            z = row.ray[:, 1]

            r_idx = np.argmin(np.abs(r - self.max_range))
            z_rr = z[r_idx]

            if np.abs(z_rr - self.receiver_depth) <= tol:
                selected_ray_idxs.append(index)

        if len(selected_ray_idxs) == 0:
            print("No rays hit receiver within tolerance.")
            return

        # Filter + sort
        selected = rays.iloc[selected_ray_idxs].copy()
        selected = selected.sort_values('bottom_bounces', ascending=True)

        unique_indices = []
        seen_sequences = []

        travel_times = []
        path_lengths = []
        r_totals = []

        for index, (_, row) in enumerate(selected.iterrows()):
            b_seq, r_coeffs = self.get_bounce_sequence(row.ray, env)
            r_total = np.prod(r_coeffs)

            total_time, total_dist = self.compute_s_t(row.ray)

            # Store values separately
            travel_times.append(total_time)
            path_lengths.append(total_dist)
            r_totals.append(r_total)

            if b_seq not in seen_sequences:
                seen_sequences.append(b_seq)
                unique_indices.append(index)

        # Attach computed values
        selected['travel_time'] = travel_times
        selected['path_length'] = path_lengths
        selected['r_total'] = r_totals

        # Select unique rays
        unique = selected.iloc[unique_indices]

        # Final sorting
        unique = unique.sort_values('bottom_bounces', ascending=True).head(num_rays)

        return unique


    def p_at_depth(self, unique, date=21, save_file=None):
        p_total = 0.0
        signal, t = generate_arms_waveform(date)
        dt = t[1] - t[0]

        # Arrival Times
        arrivals = unique['travel_time'].values
        t0 = np.min(arrivals)
        rel_arrivals = arrivals - t0
        max_delay = np.max(rel_arrivals)
        total_duration = max_delay + t[-1]
        n_samples = int(np.ceil(total_duration / dt)) + 1
        total_signal = np.zeros(n_samples)

        time_axis = np.arange(0, n_samples) * dt + t0
        for _, row in unique.iterrows():

            s = row['path_length']
            tau = row['travel_time'] - t0  # relative arrival time
            R = row['r_total']

            # Amplitude scaling
            A = R * (1 / s) * np.exp(-self.water_atten * s) * 10 ** (self.source_level / 20)
            p_total += np.abs(A)

            # Time shift index
            shift_idx = int(np.round(tau / dt))
            end_idx = shift_idx + len(signal)
            if end_idx > len(total_signal):
                end_idx = len(total_signal)

            total_signal[shift_idx:end_idx] += A * signal[:end_idx - shift_idx]

        # Filter the signal
        t_min = t0 + 0.5
        t_max = t0 + 0.95
        mask = (time_axis >= t_min) & (time_axis <= t_max)

        # Extract segment
        segment = total_signal[mask]

        # RMS calculation
        p_rms = np.sqrt(np.mean(np.abs(segment)**2))
        p_rms_dB = 20 * np.log10(p_rms)
        p_total_dB = 20 * np.log10(p_total)

        print("--------------------")
        print(f"Peak Pressure: {p_total} µPa, {p_total_dB} dB")
        print("--------------------")
        print(f"RMS Pressure: {p_rms} µPa, {p_rms_dB} dB")
        print("--------------------")

        return p_total, p_rms, total_signal, time_axis

    
    def plot_eigenrays(self, unique, save_file, figsize=(10,6)):

        fig, ax = plt.subplots(figsize=figsize)

        # Plot rays
        for _, row in unique.iterrows():
            ax.plot(row.ray[:,0] / 1000,
                    row.ray[:,1],
                    linewidth=1.5)
        
        # Plot Source and Receiver
        ax.plot(0, self.source_depth, 'bo', label="Source")
        ax.plot(self.max_range / 1000, self.receiver_depth, 'ro', label="Receiver")

        # Bathymetry
        if env is not None and 'depth' in env:
            depth = np.array(env['depth'])
            ranges = depth[:,0] / 1000
            depths = depth[:,1]
            ax.plot(ranges, depths, color='black', linewidth=2)

        ax.set_xlabel('Range (km)')
        ax.set_ylabel('Depth (m)')
        ax.set_title('Eigenrays')
        ax.grid(True)
        ax.invert_yaxis()
        ax.legend()
        plt.tight_layout()
        plt.savefig(save_file)
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
        # plt.show()
    
    

# if __name__ == "__main__":

#     # Set the Path to Include Bellhop
#     models_path = "/Users/justindiamond/Documents/Documents/UW-APL/Research/Models"
#     os.environ["PATH"] = f"{models_path}:{os.environ['PATH']}"

#     # Data Files and Parameters
#     data_dir = '/Users/justindiamond/Documents/Documents/UW-APL/Research/Diamond_Ray/data_files'
#     ssp_file = os.path.join(data_dir, 'ssp.mat')
#     bty_file = os.path.join(data_dir, 'bty.mat')
#     ati_file = os.path.join(data_dir, 'ati.mat')
#     source_level = 195 # dB re 1 μPa @ 1 m
#     freq = 25000 # Hz
#     angle_min, angle_max = -80, 80 # degrees
#     source_depth = 33
#     receiver_depth = 55
#     water_prop = (1026, 0.1)  # density (kg/m^3), attenuation (dB/m kHz)
#     bottom_prop = (2000, 1500, 0.5)  # density (kg/m^3), sound speed (m/s), attenuation (dB/m kHz)
#     surface_prop = (1000, 350, 0.1) # density (kg/m^3), sound speed (m/s), attenuation (dB/m kHz)
#     lon_start, lon_end = -122.8, -122.85
#     lat_start, lat_end = 47.78, 47.71
#     num_points = 50

#     # Run Ray Tracing
#     run = Diamond_ARL(ssp_file=ssp_file, 
#                       freq=freq,
#                       source_level=source_level,
#                       angle_min=angle_min, 
#                       angle_max=angle_max, 
#                       source_depth=source_depth, 
#                       receiver_depth=receiver_depth,
#                       bottom_prop=bottom_prop,
#                       surface_prop=surface_prop,
#                       water_prop=water_prop,
#                       lon_start=lon_start,
#                       lon_end=lon_end,
#                       lat_start=lat_start,
#                       lat_end=lat_end,
#                       num_points=num_points,
#                       bty_file=bty_file, 
#                       ati_file=ati_file,
#                       save_dir=data_dir)
    
#     # env = run.eigenray()
#     # rays = pm.compute_eigenrays(env)
#     # run.plot_rays(rays, env=env)

#     env = run.coherent_tl()
#     rays = pm.compute_transmission_loss(env)
#     run.plot_pressure_level(rays, env=env)

if __name__ == "__main__":
    # Paths
    models_path = "/Users/justindiamond/Documents/Documents/UW-APL/Research/Models"
    os.environ["PATH"] = f"{models_path}:{os.environ['PATH']}"

    data_dir = '/Users/justindiamond/Documents/Documents/UW-APL/Research/Diamond_Ray/data_files'
    bty_file = os.path.join(data_dir, 'bty.mat')
    ssp_file = os.path.join(data_dir, 'ssp_smooth.mat')
    save_file_dir = '/Users/justindiamond/Documents/Documents/UW-APL/Research/Diamond_Ray/Bellhop_Signal_Plots'
    save_file_1 = 'num_ray_5.png'
    save_file_2 = 'num_ray_5_signal.png'

    # Acoustic Properties
    ssp_depths = np.linspace(0, 200, 201)
    ssp = np.ones(len(ssp_depths)) * 1500
    num_points = 551
    bty_ranges = np.linspace(0, 5500, num_points)
    bty_depths = np.ones(len(bty_ranges)) * 200
    source_level = 195
    freq = 3500
    angle_min, angle_max = -40, 40
    source_depth = 43.5
    receiver_depth = 55
    water_prop = (1026, 0.0)
    bottom_prop = (2000, 2000, 0.0)
    lon_start, lon_end = -122.802983, -122.84
    lat_start, lat_end = 47.77295, 47.73

    # BELLHOP
    arlpy = Diamond_ARL()
    # arlpy.ssp = [[ssp_depths[i], ssp[i]] for i in range(len(ssp_depths))]
    arlpy.ssp_file = ssp_file
    arlpy.ssp = arlpy.read_ssp()
    arlpy.bty_file = bty_file
    arlpy.source_level = source_level
    arlpy.freq = freq
    arlpy.angle_min = angle_min
    arlpy.angle_max = angle_max
    arlpy.source_depth = source_depth
    arlpy.receiver_depth = receiver_depth
    arlpy.bottom_density = bottom_prop[0]
    arlpy.bottom_ss = bottom_prop[1]
    arlpy.bottom_atten = bottom_prop[2]
    arlpy.water_density = water_prop[0]
    arlpy.water_atten = water_prop[1]
    arlpy.lon_start = lon_start
    arlpy.lon_end = lon_end
    arlpy.lat_start = lat_start
    arlpy.lat_end = lat_end
    arlpy.range_points = 1000
    arlpy.depth_points = 5500
    arlpy.num_points = num_points
    # arlpy.bty, arlpy.max_range, arlpy.max_depth = [[bty_ranges[i], bty_depths[i]] for i in range(len(bty_ranges))], max(bty_ranges), max(bty_depths)
    arlpy.bty, arlpy.max_range, arlpy.max_depth = arlpy.read_bty()
    arlpy.angles = np.arange(angle_min, angle_max + 0.01, 0.01)
    arlpy.num_beams = len(arlpy.angles)

    env, rays = arlpy.eigenray()
    unique = arlpy.unique_eigenrays(rays, num_rays=5, env=env)
    p_total, p_rms, total_signal, time_axis = arlpy.p_at_depth(unique)
    arlpy.plot_eigenrays(unique, save_file=os.path.join(save_file_dir, save_file_1))


    plt.figure(figsize=(12,6))
    plt.plot(time_axis, total_signal)
    plt.xlabel("Time (s)")
    plt.ylabel("Pressure (µPa)")
    plt.title("Received Signal (All Eigenrays)")
    plt.savefig(os.path.join(save_file_dir, save_file_2))
    # plt.show()