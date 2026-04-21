from diamond_pyram import *
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation
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
    freqs = np.arange(3450, 3551, 1)   # 1 Hz spacing (important!)
    nf = len(freqs)
    angle_min, angle_max = -10, 10
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
    # pyram.ssp = ssp
    # pyram.ssp_depths = ssp_depths
    # pyram.ssp_ranges = [0]
    pyram.ssp_file = ssp_file
    pyram.ssp, pyram.ssp_depths, pyram.ssp_ranges = pyram.read_ssp()
    pyram.bty_file = bty_file
    # pyram.freq = freq
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

    # pyram.rbzb = np.column_stack((bty_ranges, bty_depths))
    pyram.rbzb = pyram.read_bty()
    receiver_range = pyram.rbzb[-1, 0]

    psv = np.zeros_like(freqs, dtype=complex)
    for f in freqs:
        print(f"Running frequency {f} Hz")
        pyram.freq = f
        pyram_model = pyram.create_model(dr=0.5)
        results = pyram_model.run()  # returns complex pressure (nz, nr)
        z = results['Depths']
        r = results['Ranges']
        z_idx = np.argmin(np.abs(z - receiver_depth))
        r__idx = np.argmin(np.abs(r - receiver_range))  # receiver at range=max_range
        p = results['CP Grid'][z_idx, r__idx]  # complex pressure at receiver depth and range
        psv[freqs == f] = p  # store complex field for each frequency

    # Time vector
    t = np.arange(0, 2.0001, 0.0001)

    # Create exponential matrix: (nfreq, nt)
    exp_matrix = np.exp(-1j * 2 * np.pi * freqs[:, None] * t[None, :])

    # Multiply and sum over frequency axis
    sig = np.sum(psv[:, None] * exp_matrix, axis=0)

    plt.plot(t, np.abs(sig))
    plt.xlabel("Time (s)")
    plt.ylabel("Amplitude")
    plt.title("Time-Domain Response at Receiver")
    plt.savefig(os.path.join(save_file_dir, "pyram_time_response.png"), dpi=200,
             bbox_inches='tight',
             pad_inches=0)  
    plt.show()

    
    # # --------------------------
    # # Run over frequencies
    # # --------------------------
    # fields = []

    # for i, f in enumerate(freqs):
    #     print(f"Running freq {f} Hz ({i+1}/{nf})")

    #     pyram.freq = f
    #     pyram_model = pyram.create_model(dr=0.5)
    #     results = pyram_model.run()

    #     # Complex pressure field
    #     field = results['TL Grid']   # (nz, nr)

    #     # Convert to dB
    #     pyram_pl = pyram.source_level - field  # TL is subtracted from source level to get PL

    #     fields.append(pyram_pl)

    # # fields: list of 2D arrays (nz_i, nr), possibly differing nz
    # nf = len(fields)

    # # 1. Trim all fields to the minimum depth dimension
    # min_nz = min(f.shape[0] for f in fields)
    # min_nr = min(f.shape[1] for f in fields)

    # fields_trimmed = [f[:min_nz, :min_nr] for f in fields]  # now (min_nz, min_nr)
    # fields_array = np.stack(fields_trimmed, axis=0)    # shape = (nf, min_nz, nr)

    # # 3. Incoherent sum over frequency
    # fields_sum = np.sum(fields_array, axis=0)       # shape = (min_nz, nr)
    # incoherent_avg = fields_sum / nf

    # # Grids
    # pyram_ranges = results['Ranges'].flatten()
    # pyram_depths = results['Depths'].flatten()


    # # --------------------------
    # # Create animation
    # # --------------------------
    # fig, ax = plt.subplots(figsize=(12, 6))

    # im = ax.imshow(fields[0],
    #                extent=[pyram_ranges.min(), pyram_ranges.max(),
    #                        pyram_depths.max(), pyram_depths.min()],
    #                aspect='auto',
    #                origin='upper',
    #                cmap='viridis',
    #                vmin=100, vmax=150)

    # title = ax.set_title(f"Frequency: {freqs[0]} Hz")
    # ax.set_xlabel("Range (m)")
    # ax.set_ylabel("Depth (m)")

    # cbar = plt.colorbar(im, ax=ax, label="Pressure Level (dB re 1 µPa)")

    # # --------------------------
    # # Update function
    # # --------------------------
    # def update(frame):
    #     im.set_data(fields[frame])
    #     title.set_text(f"Frequency: {freqs[frame]} Hz")
    #     return [im, title]

    # # --------------------------
    # # Animate
    # # --------------------------
    # ani = animation.FuncAnimation(
    #     fig,
    #     update,
    #     frames=nf,
    #     interval=150,   # ms
    #     blit=False
    # )

    # plt.tight_layout()
    # ani.save(os.path.join(save_file_dir, "pyram_f_sweep.mp4"), writer="ffmpeg", fps=10)

    # fig = plt.figure(figsize=(12, 6))
    # plt.imshow(incoherent_avg,
    #            extent=[pyram_ranges.min(), pyram_ranges.max(),
    #                    pyram_depths.max(), pyram_depths.min()],
    #            aspect='auto',
    #            origin='upper',
    #            cmap='viridis',
    #            vmin=100, vmax=150)
    # plt.colorbar(label="Pressure Level (dB re 1 µPa)")
    # plt.xlabel("Range (m)")
    # plt.ylabel("Depth (m)")
    # plt.title("Incoherent Sum over Frequencies")
    # plt.tight_layout()
    # plt.savefig(os.path.join(save_file_dir, "pyram_incoherent_sum.png"), dpi=200)
    # plt.show()