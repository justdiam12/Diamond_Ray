from diamond_pyram import *
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from scipy.interpolate import interp1d
from scipy.signal import hilbert
from transmit_sig import *
import os

def convolve_time_domain_direct(source_sig, source_time, impulse_sig, impulse_time):
    """
    Direct (non-FFT) time-domain convolution.

    Returns:
        y      : convolved signal
        t_out  : corresponding time vector
    """

    # --- Step 1: ensure same sampling ---
    dt_imp = impulse_time[1] - impulse_time[0]
    dt_src = source_time[1] - source_time[0]

    if not np.isclose(dt_imp, dt_src):
        # Resample source onto impulse grid
        t_src_uniform = np.arange(source_time[0], source_time[-1], dt_imp)
        interp_func = interp1d(source_time, source_sig,
                               bounds_error=False, fill_value=0)
        source_sig = interp_func(t_src_uniform)
        source_time = t_src_uniform

    dt = dt_imp

    # --- Step 2: direct convolution ---
    y = np.convolve(source_sig, impulse_sig, mode='full')

    # --- Step 3: output time vector ---
    t_start = source_time[0] + impulse_time[0]
    t_out = t_start + np.arange(len(y)) * dt

    return y, t_out

if __name__ == "__main__":
    # --- Paths ---
    models_path = "/Users/justindiamond/Documents/Documents/UW-APL/Research/Models"
    os.environ["PATH"] = f"{models_path}:{os.environ['PATH']}"

    data_dir = '/Users/justindiamond/Documents/Documents/UW-APL/Research/Diamond_Ray/data_files'
    bty_dir = '/Users/justindiamond/Documents/Documents/UW-APL/Research/ARMS_DATA_MASTER/Dabob_Bathymetry'
    bty_file = os.path.join(bty_dir, 'bty.mat')
    ssp_file = os.path.join(data_dir, 'ssp.mat')
    save_file_dir = '/Users/justindiamond/Documents/Documents/UW-APL/Research/Diamond_Ray/PYRAM_Simulation'
    save_file = 'pyram_response.png'

    # Acoustic Properties
    ssp_depths = np.linspace(0, 200, 201)
    ssp = np.ones(len(ssp_depths)) * 1500
    num_points = 551
    bty_ranges = np.linspace(0, 5500, num_points)
    bty_depths = np.ones(len(bty_ranges)) * 200
    source_level = 195
    freqs = np.arange(3495, 3506, 1)   # 1 Hz spacing (important!)
    nf = len(freqs)
    source_depth = 43.5
    receiver_depth = 55
    water_prop = (1026, 0.0)
    bottom_prop = (2000, 2000, 0.0)
    lon_start, lon_end = -122.802983, -122.84
    lat_start, lat_end = 47.77295, 47.73

    z_idx = None
    r_idx = None

    P_f = []

    # PYRAM 
    pyram = Diamond_PyRAM()
    pyram.ssp_file = ssp_file
    pyram.ssp, pyram.ssp_depths, pyram.ssp_ranges = pyram.read_ssp()
    pyram.bty_file = bty_file
    pyram.source_level = source_level
    pyram.source_depth = source_depth
    pyram.receiver_depth = receiver_depth
    pyram.bottom_density = bottom_prop[0]
    pyram.bottom_ss = bottom_prop[1]
    pyram.bottom_atten = bottom_prop[2]
    pyram.lon_start = lon_start
    pyram.lon_end = lon_end
    pyram.lat_start = lat_start
    pyram.lat_end = lat_end
    pyram.num_points = num_points
    pyram.rbzb = pyram.read_bty()
    receiver_range = pyram.rbzb[-1, 0]

    psv = np.zeros_like(freqs, dtype=complex)

    # Loop through frequencies, run PYRAM, and extract complex pressure at receiver location
    for f_idx, f in enumerate(freqs):
        print(f"Running frequency {f} Hz")
        pyram.freq = f
        pyram_model = pyram.create_model(dr=0.5)
        results = pyram_model.run()
        z = results['Depths']
        r = results['Ranges']
        z_idx = np.argmin(np.abs(z - receiver_depth))
        r__idx = np.argmin(np.abs(r - receiver_range))  # receiver at range=max_range
        psv[f_idx] = results['CP Grid'][z_idx, r__idx] # Relative pressure (complex) at receiver location for this frequency
    
    source_sig, source_time = generate_arms_waveform(23)  # Generate source signal (ARMS waveform) for convolution
    dt_src = source_time[1] - source_time[0]

    # Time vector
    dt = 0.0001
    T = 1 / np.min(np.diff(freqs))  
    time = np.arange(0, 2 * T + dt, dt)
    sig = np.zeros_like(time, dtype=complex)
    df = freqs[1] - freqs[0]  # Frequency spacing

    # Multiply and sum over frequency axis
    for i in range(len(time)):
        # print(f"Calculating time response at t={time[i]:.4f} s")
        sig[i] = np.trapz(psv[:] * np.exp(-1j * 2 * np.pi * freqs[:] * time[i]), freqs)


    sig[i] *= df  # Scale by frequency spacing
    sig = sig[(time >= 0.75) & (time <= 1.25)]
    time = time[(time >= 0.75) & (time <= 1.25)] + 2.3

    y, t_out = convolve_time_domain_direct(
    source_sig, source_time,
    sig, time)

    y = np.real(y) * dt # Back to relative pressure in time domain


    # fig, axs = plt.subplots(2, 1, figsize=(14, 7))
    # axs[0].plot(time, np.abs(sig))
    # axs[0].set_ylabel("Amplitude")
    # axs[0].set_title("Windowed Impulse Response at Receiver Location")
    # axs[1].plot(source_time, source_sig)
    # axs[1].set_xlabel("Time (s)")
    # axs[1].set_ylabel("Amplitude")
    # axs[1].set_title("Source Signal (ARMS LFM 100 Hz)")
    # plt.savefig(os.path.join(save_file_dir, "time_source_response.png"), dpi=200,
    #          bbox_inches='tight',
    #          pad_inches=0)
    # plt.show()

    y_env = np.abs(hilbert(np.real(y)))

    plt.figure(figsize=(12,4))
    plt.plot(t_out, np.abs(y), linewidth=1, label='Convolved Signal')
    plt.plot(t_out, y_env, linewidth=2, label='Envelope')
    plt.title("Received Signal (Convolved with ARMS Source)")
    plt.xlabel("Time (s)")
    plt.ylabel("Relative Pressure")
    plt.legend()
    plt.savefig(os.path.join(save_file_dir, "time_response.png"), dpi=200,
             bbox_inches='tight',
             pad_inches=0)
    plt.show()


    # # -------------------------------------------
    # # Impulse Response
    # # -------------------------------------------
    # plt.figure(figsize=(12,4))
    # plt.plot(time, np.abs(sig))
    # plt.xlabel("Time (s)")
    # plt.ylabel("Relative Pressure")
    # plt.title("Received Signal (Convolved with ARMS Source)")
    # plt.savefig(os.path.join(save_file_dir, "pyram_time_response_10s.png"), dpi=200,
    #          bbox_inches='tight',
    #          pad_inches=0)  
    # plt.show()


    