import os
import numpy as np
from matplotlib import pyplot as plt
from pyproj import Geod
from scipy.io import loadmat
from scipy.interpolate import interp1d, RegularGridInterpolator
from scipy.io import savemat


class Diamond_Ray_Code():
    def __init__(self, 
                 ssp_file, 
                 freq,
                 source_level,
                 angle_min,
                 angle_max,
                 angle_precision,
                 source_depth,
                 receiver_range,
                 receiver_depth,
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
        self.source_level = source_level
        self.angle_min = angle_min
        self.angle_max = angle_max
        self.angle_precision = angle_precision
        self.source_depth = source_depth
        self.receiver_range = receiver_range
        self.receiver_depth = receiver_depth

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


    def propagate_ray(self, theta0_deg, dr=5.0, r_max=100e3):
        
        r = 0.0
        z = self.source_depth
        s = 0.0
        theta = np.deg2rad(theta0_deg)
        travel_time = 0.0

        # Range, Depth, and Travel Time History
        r_hist = [r]
        z_hist = [z]
        theta_hist = [np.rad2deg(theta)]
        time_hist = [travel_time]
        s_hist = [s]
        R_hist = []
        R_coeff_hist = []

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

            # Along Track Distance and Travel Time
            ds = dr / np.cos(theta)
            travel_time += ds / c
            s += ds

            # Surface reflection
            if z_new <= self.z_surface and dz_dr < 0:
                dr_hit = (self.z_surface - z) / dz_dr
                r_new = r + dr_hit
                z_new = self.z_surface + 1e-6
                ds = dr_hit / np.cos(theta)
                travel_time -= (dr / np.cos(theta)) / c
                s -= (dr / np.cos(theta))
                travel_time += ds / c    
                s += ds
                theta_new = -theta

                # Refelction Coefficient
                theta_rel = theta
                R = self.reflection_coefficient(theta_rel, z, r, 'surface')
                R_hist.append('T')
                R_coeff_hist.append(R)

            # Bottom reflection 
            elif self.bty_interp is not None:
                z_b = self.bty_interp(r_new)

                if z_new >= z_b and dz_dr > 0:
                    dr_hit = (z_b - z) / dz_dr
                    r_new = r + dr_hit
                    z_new = z_b - 1e-6
                    travel_time -= (dr / np.cos(theta)) / c
                    s -= (dr / np.cos(theta))
                    ds = dr_hit / np.cos(theta)
                    s += ds
                    travel_time += ds / c

                    # Reflection
                    slope = self.dbty_dr_interp(r_new)
                    phi = np.arctan(slope)
                    theta_new = 2 * phi - theta

                    # Refelction Coefficient
                    theta_rel = theta - phi
                    R = self.reflection_coefficient(theta_rel, z, r, 'bottom')
                    R_hist.append('B')
                    R_coeff_hist.append(R)

            r = r_new
            z = z_new
            theta = theta_new

            r_hist.append(r)
            z_hist.append(z)
            theta_hist.append(np.rad2deg(theta))
            time_hist.append(travel_time)
            s_hist.append(s)

        return (
            np.array(r_hist),
            np.array(z_hist),
            np.array(theta_hist),
            np.array(time_hist),
            np.array(s_hist),
            np.array(R_hist),
            np.array(R_coeff_hist)
        )
    

    def reflection_coefficient(self, theta_i, z, r, medium='bottom'):

        # Local Seawater Properties
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

        # Snell's Law
        sin_theta_t = (c2 / c1) * np.sin(theta_i)

        # Complex transmission angle (critical angle physics)
        cos_theta_t = np.sqrt(1 - sin_theta_t**2 + 0j)

        R = (Z2 * np.cos(theta_i) - Z1 * cos_theta_t) / \
            (Z2 * np.cos(theta_i) + Z1 * cos_theta_t)

        return R
    

    def eigenrays(self, dtheta=1.0):

        receiver_r = self.receiver_range
        receiver_z = self.receiver_depth

        if self.bty_ranges is not None:
            bty_max = np.max(self.bty_ranges)
            r_max = min(1.01*self.receiver_range, bty_max)
        else:
            r_max = 1.01*self.receiver_range

        theta_vals = np.arange(self.angle_min, self.angle_max + dtheta, dtheta)

        eigen_list = []

        prev_theta = None
        prev_error = None

        for theta0 in theta_vals:

            print(f"Shooting {theta0:.2f}°")

            r_hist, z_hist, _, _, _, _, _ = self.propagate_ray(
                theta0,
                r_max=r_max
            )

            # If ray does not reach receiver range, skip
            if receiver_r > r_hist[-1]:
                prev_theta = theta0
                prev_error = None
                continue

            z_at_rec = np.interp(receiver_r, r_hist, z_hist)
            error = z_at_rec - receiver_z

            print(f"   depth at receiver = {z_at_rec:.3f}, error = {error:.3f}")

            # Detect sign change (root bracket)
            if prev_error is not None and error * prev_error < 0:

                print(f"   --> Bracket found between {prev_theta}° and {theta0}°")

                theta_root = self.refine_eigenray(
                    prev_theta,
                    theta0,
                    receiver_z
                )

                if theta_root is not None:

                    # Re-propagate refined ray
                    r_hist, z_hist, _, _, _, _, _ = self.propagate_ray(
                        theta_root,
                        r_max=r_max
                    )

                    if receiver_r <= r_hist[-1]:
                        z_final = np.interp(receiver_r, r_hist, z_hist)
                        miss_distance = abs(z_final - receiver_z)

                        print(f"   --> Refined depth = {z_final:.4f}, miss = {miss_distance:.4f} m")

                        # Acceptance criterion
                        if miss_distance < 0.1:
                            eigen_list.append(theta_root)
                            print(f"   --> Accepted eigenray at {theta_root:.6f}°")
                        else:
                            print("   --> Rejected (miss too large)")

            prev_theta = theta0
            prev_error = error

        return eigen_list


    def refine_eigenray(self,
                        theta_low,
                        theta_high,
                        z_rec,
                        tol=1e-4):

        r_rec = self.receiver_range
        if self.bty_ranges is not None:
            bty_max = np.max(self.bty_ranges)
            r_max = min(1.01*self.receiver_range, bty_max)
        else:
            r_max = 1.01*self.receiver_range

        for _ in range(100):

            theta_mid = 0.5 * (theta_low + theta_high)

            r_hist, z_hist, _, _, _, _, _ = self.propagate_ray(
                theta_mid,
                r_max=r_max
            )

            if r_rec > r_hist[-1]:
                return None

            z_mid = np.interp(r_rec, r_hist, z_hist)
            error_mid = z_mid - z_rec

            r_hist, z_hist, _, _, _, _, _ = self.propagate_ray(
                theta_low,
                r_max=r_max
            )

            z_low = np.interp(r_rec, r_hist, z_hist)
            error_low = z_low - z_rec

            if error_mid * error_low < 0:
                theta_high = theta_mid
            else:
                theta_low = theta_mid

            if abs(theta_high - theta_low) < tol:
                break

        return 0.5 * (theta_low + theta_high)


    # Plot Rays and SSP
    def plot_ray_fan(self):
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

                r, z, theta, time, s, r_hist, r_coeff_hist = self.propagate_ray(
                    angle,
                    r_max=self.bty_ranges[-1] if self.bty_ranges is not None else 100e3
                )

                # Store ray results
                ray_data[f"ray_{angle}deg"] = {
                    "r": r,
                    "z": z,
                    "theta": theta,
                    "time": time,
                    "distance": s,
                    "R_hist": r_hist,
                    "R_coeff_hist": r_coeff_hist
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

                r, z, theta, time, s, r_hist, r_coeff_hist = self.propagate_ray(
                    angle,
                    r_max=self.bty_ranges[-1] if self.bty_ranges is not None else 100e3
                )

                # Store ray results
                ray_data[f"ray_{angle}deg"] = {
                    "r": r,
                    "z": z,
                    "theta": theta,
                    "time": time,
                    "distance": s,
                    "R_hist": r_hist,
                    "R_coeff_hist": r_coeff_hist
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
            plt.savefig(os.path.join(self.save_dir, 'ray_paths.png'), Resolution=300)

        output_file = os.path.join(self.save_dir, "ray.mat")

        savemat(output_file, {
            "frequency": self.freq,
            "angles_deg": self.angles,
            "source_depth": self.source_depth,
            "ray_data": ray_data
        })

        print(f"Ray data saved to {output_file}")


    def pressure_at_zr(self, theta0):
        # Receiver Position
        receiver_r = self.receiver_range
        receiver_z = self.receiver_depth
        if self.bty_ranges is not None:
            bty_max = np.max(self.bty_ranges)
            r_max = min(1.01*self.receiver_range, bty_max)
        else:
            r_max = 1.01*self.receiver_range

        # Propagate Ray
        r, z, theta, time, s, r_hist, r_coeff_hist = self.propagate_ray(
            theta0,
            r_max=r_max)

        # Closest approach
        dr = r - receiver_r
        dz = z - receiver_z
        distance = np.sqrt(dr**2 + dz**2)
        idx = np.argmin(distance)

        r_closest = r[idx]
        z_closest = z[idx]
        theta_closest = theta[idx]
        time_closest = time[idx]
        s_closest = s[idx]
        
        # Reflection Product
        R_total = 1.0
        for R in r_coeff_hist:
            R_total *= R

        # # Will worry about phase later
        # c_local = self.c0  # or interpolate SSP if range-dependent
        # k = 2*np.pi*self.freq / c_local

        p = np.real(R_total * (1 / s_closest) * 10 ** (self.source_level / 20))

        return p


    def plot_eigenrays(self, eigen_list):

        receiver_r = self.receiver_range
        if self.bty_ranges is not None:
            bty_max = np.max(self.bty_ranges)
            r_max = min(1.01*self.receiver_range, bty_max)
        else:
            r_max = 1.01*self.receiver_range

        fig, ax = plt.subplots(figsize=(10, 6))

        ray_data = {}
        p_total = 0.0

        for theta0 in eigen_list:

            if theta0 is None:
                continue

            print(f"Propagating eigenray at {theta0:.6f}°")

            r, z, theta, time, s, r_hist, r_coeff_hist = self.propagate_ray(
                theta0,
                r_max=r_max
            )
            p = self.pressure_at_zr(theta0=theta0)
            p_total += np.abs(p)
            print(f"Pressure: {p} uPa, {20*np.log10(np.abs(p))} dB")

            # Store eigenray data
            ray_data[f"ray_{theta0:.6f}deg"] = {
                "r": r,
                "z": z,
                "theta": theta,
                "time": time,
                "distance": s,
                "R_hist": r_hist,
                "R_coeff_hist": r_coeff_hist
            }

            # Plot
            ax.plot(
                r / 1000,
                z,
                lw=2.5,
                label=f"{theta0:.3f}°"
            )

        # Bathymetry
        if self.bty_ranges is not None:
            ax.plot(
                self.bty_ranges / 1000,
                self.bty_depths,
                'k',
                lw=2
            )

        # Source
        ax.plot(
            0,
            self.source_depth,
            'ro',
            markersize=8,
            label="Source"
        )

        # Receiver
        ax.plot(
            receiver_r / 1000,
            self.receiver_depth,
            'ks',
            markersize=8,
            label="Receiver"
        )

        ax.set_xlabel("Range (km)")
        ax.set_ylabel("Depth (m)")
        ax.set_title("Eigenrays")
        ax.invert_yaxis()
        ax.grid()
        ax.legend()

        plt.tight_layout()
        plt.savefig(os.path.join(self.save_dir, "eigenrays.png"), dpi=300)

        print("----------")
        print(f"Total Pressure at Receiver: {p_total} uPa, {np.abs(20*np.log10(p_total))} dB")
        print("----------")

        # Save MAT file
        output_file = os.path.join(self.save_dir, "eigenrays.mat")

        savemat(output_file, {
            "frequency": self.freq,
            "angles_deg": np.array(eigen_list),
            "source_depth": self.source_depth,
            "receiver_depth": self.receiver_depth,
            "receiver_range": self.receiver_range,
            "ray_data": ray_data
        })

        print(f"Eigenray data saved to {output_file}")


if __name__ == "__main__":
    data_dir = '/Users/justindiamond/Documents/Documents/UW-APL/Research/Diamond_Ray/data_files_ssp_long'
    ssp_file = os.path.join(data_dir, 'ssp.mat')
    bty_file = os.path.join(data_dir, 'bty.mat')
    ati_file = os.path.join(data_dir, 'ati.mat')
    source_level = 195 # dB re 1 μPa @ 1 m
    freq = 3500 # Hz
    angle_min, angle_max = -50, 50
    angle_precision = 1
    source_depth = 33
    receiver_range = 8100 # Meters (m)
    receiver_depth = 55
    water_prop = (1026, 0.1)  # density (kg/m^3), attenuation (dB/m kHz)
    bottom_prop = (2000, 3000, 0)  # density (kg/m^3), sound speed (m/s), attenuation (dB/m kHz)
    surface_prop = (200, 350, 0) # density (kg/m^3), sound speed (m/s), attenuation (dB/m kHz)
    lon_start, lon_end = -122.8, -122.85
    lat_start, lat_end = 47.78, 47.71
    num_points = 1000

    # Run Ray Tracing
    ray = Diamond_Ray_Code(ssp_file=ssp_file, 
                           freq=freq,
                           source_level=source_level,
                           angle_min=angle_min, 
                           angle_max=angle_max, 
                           angle_precision=angle_precision,
                           source_depth=source_depth, 
                           receiver_range=receiver_range,
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
    
    # ray.ray_fan()
    eigen_list = ray.eigenrays()
    ray.plot_eigenrays(eigen_list)