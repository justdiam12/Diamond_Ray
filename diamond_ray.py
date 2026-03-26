import os
import numpy as np
from matplotlib import pyplot as plt
from pyproj import Geod
from scipy.io import loadmat
from scipy.interpolate import interp1d, RegularGridInterpolator
from scipy.io import savemat
from scipy.signal import chirp, windows
from transmit_sig import gen_waveform
import gc


class Diamond_Ray_Code():
    def __init__(self, 
                 ssp_file=None, 
                 freq=None,
                 source_level=None,
                 angle_min=None,
                 angle_max=None,
                 angle_precision=None,
                 source_depth=None,
                 receiver_range=None,
                 receiver_depth=None,
                 bottom_prop=None,
                 surface_prop=None,
                 water_prop=None,
                 lon_start=None,
                 lon_end=None,
                 lat_start=None,
                 lat_end=None,
                 num_points=None,
                 bty_file=None, 
                 ati_file=None,
                 signal=None,
                 signal_time=None,
                 save_dir=None):

        if ssp_file is not None:
            self.ssp_file = ssp_file
            self.ssp, self.ssp_depths, self.ssp_ranges = self.read_ssp()
        if bty_file is not None:
            self.bty_file = bty_file
        if ati_file is not None:
            self.ati_file = ati_file
        if save_dir is not None:
            self.save_dir = save_dir

        # Read SSP and BTY
        if lon_start is not None and lon_end is not None:
            self.lon_start, self.lon_end = lon_start, lon_end
        if lat_start is not None and lat_end is not None:
            self.lat_start, self.lat_end = lat_start, lat_end
        if num_points is not None:
            self.num_points = num_points
        if bty_file is not None and num_points is not None:
            self.bty_depths, self.bty_ranges = self.read_bty()
        else:
            self.bty_depths = None
            self.bty_ranges = None

        # Properties
        if freq is not None:
            self.freq = freq
        if source_level is not None:
            self.source_level = source_level
        if angle_min is not None:
            self.angle_min = angle_min
        if angle_max is not None:
            self.angle_max = angle_max
        if angle_precision is not None:
            self.angle_precision = angle_precision
        if source_depth is not None:
            self.source_depth = source_depth
        else:
            self.source_depth = self.bty_depths - 4

        if receiver_range is not None:
            self.receiver_range = receiver_range
        else:
            self.receiver_range = np.round(np.max(self.bty_ranges-1))
        self.receiver_depth = receiver_depth

        # Bottom, Surface, and Water Properties
        if bottom_prop is not None:
            self.bottom_density = bottom_prop[0] # kg/m^3
            self.bottom_ss = bottom_prop[1]      # m/s
            self.bottom_atten = bottom_prop[2] # Np/m
        if surface_prop is not None:
            self.surface_density = surface_prop[0] # kg/m^3
            self.surface_ss = surface_prop[1]      # m/s
            self.surface_atten = surface_prop[2] # Np/m
        if water_prop is not None:
            self.water_density = water_prop[0] # kg/m^3
            self.water_atten = water_prop[1] # Np/m

        # Angles
        if angle_min is not None and angle_max is not None and angle_precision is not None:
            self.angles = np.arange(self.angle_min, self.angle_max + 1, 1)

        # Signal to Append
        if signal is not None:
            self.signal = signal
        if signal_time is not None:
            self.signal_time = signal_time

        # SSP Interpolator
        if self.ssp.ndim == 1:
            # Range-independent SSP c(z)
            self.range_dependent_ssp = False
        else:
            # Range-dependent SSP c(r, z)
            self.range_dependent_ssp = True
        self.build_c_derivative()

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


    def build_c_derivative(self):
        # Range Independent SSP
        if self.range_dependent_ssp is False:

            z = self.ssp_depths
            c = self.ssp
            n = len(z)

            slopes = np.zeros(n)

            for i in range(n - 1):
                dz = z[i+1] - z[i]
                slopes[i] = (c[i+1] - c[i]) / dz

            slopes[-1] = slopes[-2]
            self.ssp_slopes = slopes

        # Range-Dependent SSP
        else:
            z = self.ssp_depths
            c = self.ssp   # shape (Nz, Nr)
            Nz, Nr = c.shape
            slopes = np.zeros((Nz, Nr))
            for i in range(Nz - 1):
                dz = z[i+1] - z[i]
                slopes[i, :] = (c[i+1, :] - c[i, :]) / dz
            slopes[-1, :] = slopes[-2, :]
            self.ssp_slopes = slopes

    def c(self, z, r):

        # Range Independent SSP
        if self.range_dependent_ssp is False:
            depths = self.ssp_depths
            cvals = self.ssp
            z = np.clip(z, depths[0], depths[-1])
            idx = np.searchsorted(depths, z) - 1
            idx = np.clip(idx, 0, len(depths)-2)
            z0 = depths[idx]
            z1 = depths[idx+1]
            c0 = cvals[idx]
            c1 = cvals[idx+1]
            w = (z - z0) / (z1 - z0)
            return c0 + w*(c1 - c0)

        # Range-Dependent SSP
        else:
            depths = self.ssp_depths
            ranges = self.ssp_ranges
            cgrid = self.ssp   # (Nz, Nr)

            # Get depth and range indices
            z = np.clip(z, depths[0], depths[-1])
            r = np.clip(r, ranges[0], ranges[-1])
            iz = np.searchsorted(depths, z) - 1
            ir = np.searchsorted(ranges, r) - 1
            iz = np.clip(iz, 0, len(depths)-2)
            ir = np.clip(ir, 0, len(ranges)-2)

            # Interpolation
            z0, z1 = depths[iz], depths[iz+1]
            r0, r1 = ranges[ir], ranges[ir+1]
            c00 = cgrid[iz,   ir]
            c01 = cgrid[iz+1, ir]
            c10 = cgrid[iz,   ir+1]
            c11 = cgrid[iz+1, ir+1]
            wz = (z - z0) / (z1 - z0)
            wr = (r - r0) / (r1 - r0)
            c0 = c00 + wz*(c01 - c00)
            c1 = c10 + wz*(c11 - c10)

            return c0 + wr*(c1 - c0)


    def dc_dz(self, z, r):
        
        # Range Independent SSP
        if self.range_dependent_ssp is False:
            depths = self.ssp_depths
            slopes = self.ssp_slopes
            z = np.clip(z, depths[0], depths[-1])
            idx = np.searchsorted(depths, z) - 1
            idx = np.clip(idx, 0, len(depths)-2)
            return slopes[idx]

        # Range-Dependent SSP
        else:
            depths = self.ssp_depths
            ranges = self.ssp_ranges
            slopes = self.ssp_slopes  # (Nz, Nr)

            z = np.clip(z, depths[0], depths[-1])
            r = np.clip(r, ranges[0], ranges[-1])

            iz = np.searchsorted(depths, z) - 1
            ir = np.searchsorted(ranges, r) - 1

            iz = np.clip(iz, 0, len(depths)-2)
            ir = np.clip(ir, 0, len(ranges)-2)

            r0, r1 = ranges[ir], ranges[ir+1]
            wr = (r - r0) / (r1 - r0)

            # interpolate along range
            dc0 = slopes[iz, ir]
            dc1 = slopes[iz, ir+1]

            return dc0 + wr*(dc1 - dc0)
    

    def propagate_ray(self, theta0_deg, dr=10, r_max=100e3):

        r = 0.0
        z = self.source_depth
        tau = np.tan(np.deg2rad(theta0_deg))

        travel_time = 0.0
        s = 0.0

        r_hist = [r]
        z_hist = [z]
        theta_hist = [theta0_deg]
        time_hist = [travel_time]
        s_hist = [s]
        R_hist = []
        R_coeff_hist = []

        eps = 1e-6

        while r < r_max:

            c = self.c(z, r)
            dc_dz = self.dc_dz(z, r)

            tang = -dc_dz / c

            # RK4 coefficients
            k1_tau = (1 + tau**2) * tang
            k1_z = tau

            z2 = z + k1_z * dr / 2
            tau2 = tau + k1_tau * dr / 2
            c2 = self.c(z2, r + dr/2)
            dc_dz2 = self.dc_dz(z2, r + dr/2)
            tang2 = -dc_dz2 / c2

            k2_tau = (1 + tau2**2) * tang2
            k2_z = tau2

            z3 = z + k2_z * dr / 2
            tau3 = tau + k2_tau * dr / 2
            c3 = self.c(z3, r + dr/2)
            dc_dz3 = self.dc_dz(z3, r + dr/2)
            tang3 = -dc_dz3 / c3

            k3_tau = (1 + tau3**2) * tang3
            k3_z = tau3

            z4 = z + k3_z * dr
            tau4 = tau + k3_tau * dr
            c4 = self.c(z4, r + dr)
            dc_dz4 = self.dc_dz(z4, r + dr)
            tang4 = -dc_dz4 / c4

            k4_tau = (1 + tau4**2) * tang4
            k4_z = tau4

            # RK4 update
            z_new = z + (dr/6) * (k1_z + 2*k2_z + 2*k3_z + k4_z)
            tau_new = tau + (dr/6) * (k1_tau + 2*k2_tau + 2*k3_tau + k4_tau)

            r_new = r + dr

            theta = np.arctan(tau)

            # path length
            ds = dr / np.cos(theta)
            s += ds
            travel_time += ds / c

            # Surface reflection
            if z_new <= self.z_surface and tau < 0:

                dr_hit = (self.z_surface - z) / tau

                r_new = r + dr_hit
                z_new = self.z_surface + eps
                tau_new = -tau

                ds = dr_hit / np.cos(theta)
                s += ds
                travel_time += ds / c

                theta_rel = theta
                R = self.reflection_coefficient(theta_rel, z, r, 'surface')

                R_hist.append('T')
                R_coeff_hist.append(R)

            # Bottom reflection
            elif self.bty_interp is not None:

                z_b = self.bty_interp(r_new)

                if z_new >= z_b and tau > 0:

                    dr_hit = (z_b - z) / tau

                    r_new = r + dr_hit
                    z_new = z_b - eps
                    slope = self.dbty_dr_interp(r_new)
                    phi = np.arctan(slope)
                    theta_new = 2 * phi - theta
                    tau_new = np.tan(theta_new)

                    ds = dr / np.cos(theta)
                    s += ds
                    travel_time += ds / c

                    theta_rel = theta - phi
                    R = self.reflection_coefficient(theta_rel, z, r, 'bottom')

                    R_hist.append('B')
                    R_coeff_hist.append(R)

            # update state
            r = r_new
            z = z_new
            tau = tau_new

            r_hist.append(r)
            z_hist.append(z)
            theta_hist.append(np.rad2deg(np.arctan(tau)))
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
    

    def eigenrays(self, dtheta=0.1, tol=1):

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

            # print(f"Shooting {theta0:.2f}°")

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

            # print(f"   depth at receiver = {z_at_rec:.3f}, error = {error:.3f}")

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
                        if miss_distance < tol:
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
                        tol=1e-10):

        r_rec = self.receiver_range

        if self.bty_ranges is not None:
            bty_max = np.max(self.bty_ranges)
            r_max = min(1.01*self.receiver_range, bty_max)
        else:
            r_max = 1.01*self.receiver_range

        # --- shoot lower ray ---
        r_hist, z_hist, *_ = self.propagate_ray(theta_low, r_max=r_max)
        if r_rec > r_hist[-1]:
            return None

        z_low = np.interp(r_rec, r_hist, z_hist)
        error_low = z_low - z_rec

        # --- shoot upper ray ---
        r_hist, z_hist, *_ = self.propagate_ray(theta_high, r_max=r_max)
        if r_rec > r_hist[-1]:
            return None

        z_high = np.interp(r_rec, r_hist, z_hist)
        error_high = z_high - z_rec

        # Ensure bracket
        if error_low * error_high > 0:
            return None

        for i in range(1000):

            theta_mid = 0.5 * (theta_low + theta_high)

            r_hist, z_hist, *_ = self.propagate_ray(theta_mid, r_max=r_max)
            if r_rec > r_hist[-1]:
                return None

            z_mid = np.interp(r_rec, r_hist, z_hist)
            error_mid = z_mid - z_rec

            # Determine which bracket contains receiver
            if error_low * error_mid < 0:
                theta_high = theta_mid
                error_high = error_mid
            else:
                theta_low = theta_mid
                error_low = error_mid

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
            plt.savefig(os.path.join(self.save_dir, 'ray_paths_1m.png'))

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

        # Distance Travelled
        s_closest = s[idx]
        arrival_t = time[idx]
        
        # Reflection Product
        R_total = 1.0
        for R in r_coeff_hist:
            R_total *= R

        # # Will worry about phase later
        # c_local = self.c0  # or interpolate SSP if range-dependent
        # k = 2*np.pi*self.freq / c_local

        p = np.real(R_total * (1 / s_closest) * np.exp(-self.water_atten * s_closest) * 10 ** (self.source_level / 20))
        # print(f"Pressure: {p} uPa, {20*np.log10(np.abs(p))} dB, Distance Travelled: {s_closest} m, Bounces: {r_hist}")

        return p, arrival_t
    

    def coherent_pressure(self, p_list, arrivals):

        p_list = np.asarray(p_list)
        arrivals = np.asarray(arrivals)

        # Align arrivals to earliest arrival
        t0 = np.min(arrivals)
        arrivals = arrivals - t0

        signal = self.signal
        t_signal = self.signal_time

        dt = t_signal[1] - t_signal[0]

        # Maximum delay in samples
        max_delay = int(np.ceil(np.max(arrivals) / dt))

        out_len = len(signal) + max_delay
        total_signal = np.zeros(out_len)

        for i in range(len(p_list)):

            amp = p_list[i]
            delay = int(np.round(arrivals[i] / dt))

            total_signal[delay:delay + len(signal)] += amp * signal

        # Sampling frequency
        fs = 1 / dt

        # Define window: 0.1s to 0.9s
        start_idx = int(0.1 * fs)
        end_idx = int(0.9 * fs)

        # Ensure bounds are valid
        start_idx = min(start_idx, len(total_signal))
        end_idx = min(end_idx, len(total_signal))

        window = total_signal[start_idx:end_idx]

        # RMS over window
        p_rms = np.sqrt(np.mean(window**2))
        rl = 20 * np.log10(p_rms)

        print(f"Pressure (RMS): {p_rms} uPa, Receiver Level: {rl} dB")

        return p_rms, total_signal


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

        p_list = []
        arrivals = []
        for theta0 in eigen_list:

            if theta0 is None:
                continue

            print(f"Propagating eigenray at {theta0:.6f}°")

            r, z, theta, time, s, r_hist, r_coeff_hist = self.propagate_ray(
                theta0,
                r_max=r_max
            )
            p, arrival_t = self.pressure_at_zr(theta0=theta0)

            p_list.append(p)
            arrivals.append(arrival_t)

            p_total += np.abs(p)

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

        # Coherent Pressure
        p_rms, total_signal = self.coherent_pressure(p_list, arrivals)

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

        x = np.arange(len(total_signal)) / 80000
        plt.figure()
        plt.plot(x, total_signal)
        plt.xlabel("Time (s)")
        plt.ylabel("Pressure (uPa)")
        plt.title("Received Signal at Depth = 55 meters")
        plt.savefig("signal.png", dpi=300)

        savemat(output_file, {
            "frequency": self.freq,
            "angles_deg": np.array(eigen_list),
            "source_depth": self.source_depth,
            "receiver_depth": self.receiver_depth,
            "receiver_range": self.receiver_range,
            "ray_data": ray_data,
            "total_signal": total_signal,
            "p_rms": p_rms
        })

        print(f"Eigenray data saved to {output_file}")


# if __name__ == "__main__":
#     data_dir = '/Users/justindiamond/Documents/Documents/UW-APL/Research/ARMS_DATA_MASTER/LIVEOCEAN_DAILY/2020-01-24'
#     ssp_file = os.path.join(data_dir, 'ssp.mat')
#     bty_file = os.path.join(data_dir, 'bty.mat')
#     ati_file = os.path.join(data_dir, 'ati.mat')
#     source_level = 195 # dB re 1 μPa @ 1 m
#     freq = 3500 # Hz
#     angle_min, angle_max = -35, 35
#     angle_precision = 1
#     source_depth = 35
#     # receiver_range = 5400 # Meters (m)
#     receiver_depth = 25
#     water_prop = (1026, 1e-9)  # density (kg/m^3), attenuation (dB/m kHz)
#     bottom_prop = (2000, 3000, 1e-6)  # density (kg/m^3), sound speed (m/s), attenuation (Np/m)
#     surface_prop = (200, 350, 1e-3) # density (kg/m^3), sound speed (m/s), attenuation (Np/m)
#     lon_start, lon_end = -122.8, -122.84
#     lat_start, lat_end = 47.78, 47.73
#     num_points = 1000

#     # Data for signal
#     sample_freq = 80000
#     desired_pulse_start_freq_khz = 3.450
#     desired_pulse_stop_freq_khz = 3.550
#     desired_pulse_duration_ms = 1000
#     percent_taper = 10

#     signal, time_s = gen_waveform(
#         sample_freq,
#         desired_pulse_start_freq_khz,
#         desired_pulse_stop_freq_khz,
#         desired_pulse_duration_ms,
#         percent_taper
#     )

#     # Run Ray Tracing
#     ray = Diamond_Ray_Code(ssp_file=ssp_file, 
#                            freq=freq,
#                            source_level=source_level,
#                            angle_min=angle_min, 
#                            angle_max=angle_max, 
#                            angle_precision=angle_precision,
#                            source_depth=source_depth, 
#                            receiver_range=None,
#                            receiver_depth=receiver_depth,
#                            bottom_prop=bottom_prop,
#                            surface_prop=surface_prop,
#                            water_prop=water_prop,
#                            lon_start=lon_start,
#                            lon_end=lon_end,
#                            lat_start=lat_start,
#                            lat_end=lat_end,
#                            num_points=num_points,
#                            bty_file=bty_file, 
#                            ati_file=ati_file,
#                            signal=signal,
#                            signal_time=time_s,
#                            save_dir=data_dir)
    
#     # bty_ranges = np.arange(0, 5400+0.01, 0.01)
#     # bty_depths = np.ones_like(bty_ranges)*200
#     # ray.bty_ranges = bty_ranges
#     # ray.bty_depths = bty_depths
#     # ray.bty_interp = interp1d(
#     #     ray.bty_ranges,
#     #     ray.bty_depths,
#     #     bounds_error=False,
#     #     fill_value="extrapolate"
#     # )
#     # dz_dr = np.gradient(ray.bty_depths, ray.bty_ranges)
#     # ray.dbty_dr_interp = interp1d(
#     #     ray.bty_ranges,
#     #     dz_dr,
#     #     bounds_error=False,
#     #     fill_value="extrapolate"
#     # )
#     # ray.angles = [1.79, 1.03, 4.47, 7.41, -1.64]

#     # ray.ray_fan()
#     eigen_list = ray.eigenrays()
#     ray.plot_eigenrays(eigen_list)

#     # ray.plot_ray_fan()

if __name__ == "__main__":
    data_dir = '/Users/justindiamond/Documents/Documents/UW-APL/Research/ARMS_DATA_MASTER/LIVEOCEAN_DAILY/'
    daily_files = ['2020-01-21', '2020-01-24', '2020-01-28', '2020-01-29']

    # Preliminary Information
    source_level = 195 # dB re 1 μPa @ 1 m
    freq = 3500 # Hz
    angle_min, angle_max = -35, 35
    angle_precision = 1
    source_depth = 35
    water_prop = (1026, 1e-9)  # density (kg/m^3), attenuation (dB/m kHz)
    bottom_prop = (2000, 3000, 1e-6)  # density (kg/m^3), sound speed (m/s), attenuation (Np/m)
    surface_prop = (200, 350, 1e-3) # density (kg/m^3), sound speed (m/s), attenuation (Np/m)
    lon_start, lon_end = -122.80275, -122.84
    lat_start, lat_end = 47.7729, 47.73
    num_points = 1000

    # Data for signal
    sample_freq = 80000
    desired_pulse_start_freq_khz = 3.450
    desired_pulse_stop_freq_khz = 3.550
    desired_pulse_duration_ms = 1000
    percent_taper = 10

    signal, time_s = gen_waveform(
        sample_freq,
        desired_pulse_start_freq_khz,
        desired_pulse_stop_freq_khz,
        desired_pulse_duration_ms,
        percent_taper
    )


    for i in range(len(daily_files)):
    # for i in range(1):

        print("----------------------")
        print(f"Working on: {daily_files[i]}")
        print("----------------------")
        
        # Data Files
        ssp_file = os.path.join(data_dir, daily_files[i], 'ssp.mat')
        bty_file = os.path.join(data_dir, daily_files[i], 'bty.mat')
        ati_file = os.path.join(data_dir, daily_files[i], 'ati.mat')

        if daily_files[i] == '2020-01-21':
            day = loadmat("/Users/justindiamond/Documents/Documents/UW-APL/Research/ARMS_DATA_MASTER/array_depths/depths_20200121.mat")
            unique_depths = day["unique_depths"]
        elif daily_files[i] == '2020-01-24':
            day = loadmat("/Users/justindiamond/Documents/Documents/UW-APL/Research/ARMS_DATA_MASTER/array_depths/depths_20200124.mat")
            unique_depths = day["unique_depths"]
        elif daily_files[i] == '2020-01-28':
            day = loadmat("/Users/justindiamond/Documents/Documents/UW-APL/Research/ARMS_DATA_MASTER/array_depths/depths_20200128.mat")
            unique_depths = day["unique_depths"]
        else:
            day = loadmat("/Users/justindiamond/Documents/Documents/UW-APL/Research/ARMS_DATA_MASTER/array_depths/depths_20200129.mat")
            unique_depths = day["unique_depths"]

        unique_depths = np.array(unique_depths).squeeze(axis=0)
        p = np.zeros_like(unique_depths)
        p_rms = np.zeros_like(unique_depths)


        for d in range(len(unique_depths)):
        # for d in range(1):
            # Run Ray Tracing
            ray = Diamond_Ray_Code(ssp_file=ssp_file, 
                                freq=freq,
                                source_level=source_level,
                                angle_min=angle_min, 
                                angle_max=angle_max, 
                                angle_precision=angle_precision,
                                source_depth=source_depth, 
                                receiver_range=None,
                                receiver_depth=int(unique_depths[d]),
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
                                signal=signal,
                                signal_time=time_s,
                                save_dir=data_dir)
            
            eigen_list = ray.eigenrays()

            p_list = []
            arrivals = []
            for theta in eigen_list:
                pressure, arrival_t = ray.pressure_at_zr(theta)
                p_list.append(np.abs(pressure))
                arrivals.append(arrival_t)

            p_coherent, _ = ray.coherent_pressure(p_list, arrivals)

            p[d] = sum(p_list)
            p_rms[d] = p_coherent

        save_path = os.path.join(data_dir, daily_files[i], 'diamond_ray.mat')

        print("----------------------")
        print(f"Saving MAT File for {daily_files[i]}")
        print("----------------------")
        savemat(save_path, {
            "depths": unique_depths,
            "pressure": p,
            "pressure_rms": p_rms
        })

        # --- Clear large variables ---
        del day, unique_depths, p, p_rms
        del ray, eigen_list, p_list, arrivals, pressure, arrival_t, p_coherent

        gc.collect()