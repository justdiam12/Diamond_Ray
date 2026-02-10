import os
import numpy as np
from matplotlib import pyplot as plt
from pyproj import Geod
from scipy.io import loadmat
from scipy.interpolate import interp1d, RegularGridInterpolator
from scipy.io import savemat


class Diamond_Ray():
    def __init__(self, 
                 ssp_file, 
                 freq,
                 angle_min,
                 angle_max,
                 source_depth,
                 bottom_prop,
                 surface_prop,
                 water_prop,
                 lon_start,
                 lon_end,
                 lat_start,
                 lat_end,
                 num_points,
                 bty_file, 
                 ati_file,
                 save_dir):

        self.ssp_file = ssp_file
        self.bty_file = bty_file
        self.ati_file = ati_file
        self.save_dir = save_dir

        # Read data
        self.ssp, self.ssp_depths, self.ssp_ranges = self.read_ssp()

        self.freq = freq
        self.angle_min = angle_min
        self.angle_max = angle_max
        self.source_depth = source_depth

        # Bottom, Surface, and Water Properties
        self.bottom_density = bottom_prop[0] # kg/m^3
        self.bottom_ss = bottom_prop[1]      # m/s
        self.bottom_atten = bottom_prop[2]   # 
        self.surface_density = surface_prop[0] # kg/m^3
        self.surface_ss = surface_prop[1]      # m/s
        self.surface_atten = surface_prop[2]   #
        self.water_density = water_prop[0] # kg/m^3
        self.water_atten = water_prop[1]   # dB/m kHz

        self.lon_start, self.lon_end = lon_start, lon_end
        self.lat_start, self.lat_end = lat_start, lat_end
        self.num_points = num_points

        self.bty_depths, self.bty_ranges = self.read_bty()
        self.angles = np.arange(self.angle_min, self.angle_max + 1, 1)

        # SSP Interpolator
        if self.ssp.ndim == 1:
            # Range-independent SSP c(z)
            self.range_dependent_ssp = False

            self.c_interp = interp1d(
                self.ssp_depths,
                self.ssp,
                bounds_error=False,
                fill_value="extrapolate"
            )

            dc_dz = np.gradient(self.ssp, self.ssp_depths)
            self.dc_dz_interp = interp1d(
                self.ssp_depths,
                dc_dz,
                bounds_error=False,
                fill_value="extrapolate"
            )

            self.z_bottom = self.ssp_depths.max()

        else:
            # Range-dependent SSP c(r, z)
            self.range_dependent_ssp = True

            # Ensure shape = (Nr, Nz)
            if self.ssp.shape == (len(self.ssp_depths), len(self.ssp_ranges)):
                self.ssp = self.ssp.T

            self.c_interp = RegularGridInterpolator(
                (self.ssp_ranges, self.ssp_depths),
                self.ssp,
                bounds_error=False,
                fill_value=None
            )

            dc_dz = np.gradient(self.ssp, self.ssp_depths, axis=1)
            self.dc_dz_interp = RegularGridInterpolator(
                (self.ssp_ranges, self.ssp_depths),
                dc_dz,
                bounds_error=False,
                fill_value=None
            )

            self.z_bottom = self.ssp_depths.max()

        # Surface
        self.z_surface = 0.0

        # Bathymetry Interpolator
        if self.bty_ranges is not None:
            self.bty_interp = interp1d(
                self.bty_ranges,
                self.bty_depths,
                bounds_error=False,
                fill_value="extrapolate"
            )

            dz_dr = np.gradient(self.bty_depths, self.bty_ranges)
            self.dbty_dr_interp = interp1d(
                self.bty_ranges,
                dz_dr,
                bounds_error=False,
                fill_value="extrapolate"
            )
        else:
            self.bty_interp = None
            self.dbty_dr_interp = None


    # Helper Functions
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

        return np.abs(bty_depths), bty_ranges


    def c(self, z, r):
        if self.range_dependent_ssp:
            return self.c_interp((r, z))
        else:
            return self.c_interp(z)


    def dc_dz(self, z, r):
        if self.range_dependent_ssp:
            return self.dc_dz_interp((r, z))
        else:
            return self.dc_dz_interp(z)


    def propagate_ray(self, theta0_deg, dr=10.0, r_max=100e3):

        r = self.bty_ranges[0] if self.bty_ranges is not None else 0.0
        z = self.source_depth
        theta = np.deg2rad(theta0_deg)

        amp = 1 + 0j  # complex amplitude

        r_hist = [r]
        z_hist = [z]
        theta_hist = [theta]
        amp_hist = [amp]

        while r < r_max:

            theta = np.clip(theta, -np.pi/2 + 1e-6, np.pi/2 - 1e-6)

            c = self.c(z, r)
            dc_dz = self.dc_dz(z, r)

            if c is None:
                break

            # Ray equations
            dz_dr = np.tan(theta)
            dtheta_dr = -(1.0 / c) * dc_dz / np.cos(theta)

            # Step forward
            r_new = r + dr
            z_new = z + dz_dr * dr
            theta_new = theta + dtheta_dr * dr

            # Path length increment
            ds = dr / np.cos(theta)

            # Water Absorption
            alpha_db = self.water_atten * (self.freq / 1000.0)  # dB/m
            alpha_np = alpha_db / 8.686  # nepers/m

            amp *= np.exp(-alpha_np * ds)

            # Spherical Spreading
            if r > 0.0:
                amp *= (r / r_new)
            else:
                # avoid singularity at source
                amp *= 1.0

            # Surface Reflection
            if z_new <= self.z_surface:

                theta_i = theta
                R = self.reflection_coefficient(theta_i, z, r, medium='surface')
                amp *= R

                z_new = -z_new
                theta_new = -theta_new

            # Bottom Reflection
            elif self.bty_interp is not None:

                z_b = self.bty_interp(r_new)

                if z_new >= z_b:

                    slope = self.dbty_dr_interp(r_new)
                    phi = np.arctan(slope)

                    theta_rel = theta - phi
                    R = self.reflection_coefficient(theta_rel, z, r, medium='bottom')

                    amp *= R

                    z_new = z_b - (z_new - z_b)
                    theta_new = 2 * phi - theta_new

            elif z_new >= self.z_bottom:
                break

            # Update state
            r = r_new
            z = z_new
            theta = theta_new

            r_hist.append(r)
            z_hist.append(z)
            theta_hist.append(theta)
            amp_hist.append(amp)

        return (
            np.array(r_hist),
            np.array(z_hist),
            np.array(theta_hist),
            np.array(amp_hist)
        )
    

    def reflection_coefficient(self, theta_i, z, r, medium='bottom'):

        # Water properties (use user-specified density + local c)
        rho1 = self.water_density
        c1 = self.c(z, r)

        if medium == 'bottom':
            rho2 = self.bottom_density
            c2 = self.bottom_ss
        else:
            rho2 = self.surface_density
            c2 = self.surface_ss

        Z1 = rho1 * c1
        Z2 = rho2 * c2

        sin_theta_t = (c2 / c1) * np.sin(theta_i)

        # Allow complex transmission angle (critical angle physics)
        cos_theta_t = np.sqrt(1 - sin_theta_t**2 + 0j)

        R = (Z2 * np.cos(theta_i) - Z1 * cos_theta_t) / \
            (Z2 * np.cos(theta_i) + Z1 * cos_theta_t)

        return R


    # Plot Rays and SSP
    def run(self):
        ray_data = {}

        # Filter for Range Independent SSP
        if self.ssp.ndim == 1 or self.ssp.shape[1] == 1:
            # Figure
            fig, (ax_ssp, ax_ray) = plt.subplots(
                ncols=2,
                figsize=(14, 6),
                gridspec_kw={'width_ratios': [1, 3]},
                sharey=True
            )

            # Range-independent SSP: just plot line
            ax_ssp.plot(self.ssp, self.ssp_depths, 'k', lw=2)
            ax_ssp.set_xlabel('c (m/s)')
            ax_ssp.set_title('Sound Speed Profile')
            ax_ssp.set_ylabel('Depth (m)')
            ax_ssp.invert_yaxis()
            ax_ssp.grid()

            for angle in self.angles:
                print(f"Propagating ray at {angle}°")

                r, z, theta, amp = self.propagate_ray(
                    angle,
                    r_max=self.bty_ranges[-1] if self.bty_ranges is not None else 100e3
                )

                # Store ray results
                ray_data[f"ray_{angle}deg"] = {
                    "r": r,
                    "z": z,
                    "theta": theta,
                    "amp": amp
                }

                ax_ray.plot(r / 1000, z)

            # Plot bathymetry
            if self.bty_ranges is not None:
                ax_ray.plot(self.bty_ranges / 1000, self.bty_depths, 'k', lw=2)

            ax_ray.set_xlabel('Range (km)')
            ax_ray.set_title('Ray Paths')
            ax_ray.grid()
            ax_ray.legend(loc='upper right', fontsize=8)

            plt.tight_layout()
            plt.savefig(os.path.join(self.save_dir, 'ray_paths.png'))

        # Range-Dependent SSP
        else:
            fig, (ax_ray, ax_ssp) = plt.subplots(
                nrows=2,
                figsize=(14, 8),
                sharex=True,
                gridspec_kw={'height_ratios': [1, 1]},
                constrained_layout=True
            )

            # Ray Paths
            for angle in self.angles:
                print(f"Propagating ray at {angle}°")

                r, z, theta, amp = self.propagate_ray(
                    angle,
                    r_max=self.bty_ranges[-1] if self.bty_ranges is not None else 100e3
                )

                # Store ray results
                ray_data[f"ray_{angle}deg"] = {
                    "r": r,
                    "z": z,
                    "theta": theta,
                    "amp": amp
                }

                ax_ray.plot(r / 1000, z)

            # Bathymetry
            if self.bty_ranges is not None:
                ax_ray.plot(self.bty_ranges / 1000, self.bty_depths, 'k', lw=2)

            ax_ray.set_ylabel('Depth (m)')
            ax_ray.set_title('Ray Paths')
            ax_ray.invert_yaxis()
            ax_ray.grid()

            # SSP
            ssp_plot = np.copy(self.ssp)

            # Ensure (depth, range)
            if ssp_plot.ndim == 2 and ssp_plot.shape[0] != len(self.ssp_depths):
                ssp_plot = ssp_plot.T

            # Mask below bathymetry
            if self.bty_ranges is not None and ssp_plot.ndim == 2:
                bty_interp_for_plot = np.interp(
                    self.ssp_ranges,
                    self.bty_ranges,
                    self.bty_depths
                )
                for i, z_b in enumerate(bty_interp_for_plot):
                    ssp_plot[self.ssp_depths > z_b, i] = np.nan

            im = ax_ssp.imshow(
                ssp_plot,
                extent=[
                    self.ssp_ranges[0] / 1000,
                    self.ssp_ranges[-1] / 1000,
                    self.ssp_depths[-1],
                    self.ssp_depths[0]
                ],
                aspect='auto',
                cmap='viridis'
            )

            # Bathymetry overlay
            if self.bty_ranges is not None:
                ax_ssp.plot(self.bty_ranges / 1000, self.bty_depths, 'k', lw=2)

            ax_ssp.set_xlabel('Range (km)')
            ax_ssp.set_ylabel('Depth (m)')
            ax_ssp.set_title('Sound Speed Profile')
            ax_ssp.grid()

            # Colorbar
            cbar = fig.colorbar(im, ax=ax_ssp, location='right')
            cbar.set_label('c (m/s)')
            plt.savefig(os.path.join(self.save_dir, 'ray_paths.png'))

        output_file = os.path.join(self.save_dir, "ray.mat")

        savemat(output_file, {
            "frequency": self.freq,
            "angles_deg": self.angles,
            "source_depth": self.source_depth,
            "ray_data": ray_data
        })

        print(f"Ray data saved to {output_file}")


if __name__ == "__main__":
    data_dir = '/Users/justindiamond/Documents/Documents/UW-APL/Research/Diamond_Ray/data_files_ssp_long'
    ssp_file = os.path.join(data_dir, 'ssp.mat')
    bty_file = os.path.join(data_dir, 'bty.mat')
    ati_file = os.path.join(data_dir, 'ati.mat')
    freq = 3500 # Hz
    angle_min, angle_max = -10, 10
    source_depth = 33
    water_prop = (1026, 0.1)  # density (kg/m^3), attenuation (dB/m kHz)
    bottom_prop = (2000, 1500, 0.5)  # density (kg/m^3), sound speed (m/s), attenuation (dB/m kHz)
    surface_prop = (1000, 350, 0.1) # density (kg/m^3), sound speed (m/s), attenuation (dB/m kHz)
    lon_start, lon_end = -122.8, -122.85
    lat_start, lat_end = 47.78, 47.71
    num_points = 1000

    # Run Ray Tracing
    ray = Diamond_Ray(ssp_file=ssp_file, 
                      freq=freq,
                      angle_min=angle_min, 
                      angle_max=angle_max, 
                      source_depth=33, 
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
    
    ray.run()