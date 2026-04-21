from diamond_pyram import *
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from scipy.interpolate import interp1d
from scipy.signal import hilbert
from transmit_sig import *
import os

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
    freqs = np.arange(3450, 3551, 1)   # 1 Hz spacing
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
        print(f"Depth: {z[np.argmin(np.abs(z - receiver_depth))]}")
        print(f"Range: {r[np.argmin(np.abs(r - receiver_range))]}")
        z_idx = np.argmin(np.abs(z - receiver_depth))
        r_idx = np.argmin(np.abs(r - receiver_range))  # receiver at range=max_range
        psv[f_idx] = 10 ** ((195 - results['TL Grid'][z_idx, r_idx]) / 20) # µPa at receiver location for this frequency

    source_sig, source_time = generate_arms_waveform(23) # LFM 100 Hz source signal, centered at 3500 Hz

    plt.figure(figsize=(12,4))
    plt.plot(freqs, 20 * np.log10(np.abs(psv)))
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Receiver Pressure (dB re 1 µPa)")
    plt.title("Receiver Pressure Level (r = 5506 m, z = 55 m) vs. Frequency")
    plt.savefig(os.path.join(save_file_dir, "pl_vs_freq.png"))
    plt.show()