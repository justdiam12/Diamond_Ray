from diamond_arl import *
from diamond_pyram import *
import numpy as np
from matplotlib import pyplot as plt
import os

if __name__ == "__main__":
    # --- Paths ---
    models_path = "/Users/justindiamond/Documents/Documents/UW-APL/Research/Models"
    os.environ["PATH"] = f"{models_path}:{os.environ['PATH']}"

    data_dir = '/Users/justindiamond/Documents/Documents/UW-APL/Research/Diamond_Ray/data_files'
    bty_file = os.path.join(data_dir, 'bty.mat')
    ssp_file = os.path.join(data_dir, 'ssp.mat')

    # Acoustic Properties
    ssp_depths = np.linspace(0, 200, 201)
    ssp = np.ones(len(ssp_depths)) * 1500
    bty_ranges = np.linspace(0, 10000, 1001)
    bty_depths = np.ones(len(bty_ranges)) * 200
    source_level = 195
    freq = 3500
    angle_min, angle_max = -30, 30
    source_depth = 25
    receiver_depth = 55
    water_prop = (1026, 0.1)
    bottom_prop = (2000, 3000, 0.2)
    lon_start, lon_end = -122.8, -122.84
    lat_start, lat_end = 47.78, 47.73
    num_points = 85

    # PYRAM 
    pyram = Diamond_PyRAM()
    pyram.ssp = ssp
    pyram.ssp_depths = ssp_depths
    pyram.ssp_ranges = [0]
    pyram.ssp_file = ssp_file
    pyram.ssp, pyram.ssp_depths, pyram.ssp_ranges = pyram.read_ssp()
    pyram.bty_file = bty_file
    pyram.freq = freq
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
    pyram.rbzb = np.column_stack((bty_ranges, bty_depths))
    pyram.rbzb = pyram.read_bty()
    pyram_model = pyram.create_model(dr=0.5)
    results = pyram_model.run()

    pyram_pl = pyram.source_level - results["TL Grid"]
    pyram_ranges = results["Ranges"] / 1000.0
    pyram_depths = results["Depths"]

    # BELLHOP
    arlpy = Diamond_ARL()
    arlpy.ssp = [[ssp_depths[i], ssp[i]] for i in range(len(ssp_depths))]
    arlpy.ssp_file = ssp_file
    arlpy.ssp = arlpy.read_ssp()
    arlpy.bty_file = bty_file
    arlpy.source_level = source_level
    arlpy.freq = freq
    arlpy.angle_min = angle_min
    arlpy.angle_max = angle_max
    arlpy.source_depth = source_depth
    arlpy.bottom_density = bottom_prop[0]
    arlpy.bottom_ss = bottom_prop[1]
    arlpy.bottom_atten = bottom_prop[2]
    arlpy.water_density = water_prop[0]
    arlpy.water_atten = water_prop[1]
    arlpy.lon_start = lon_start
    arlpy.lon_end = lon_end
    arlpy.lat_start = lat_start
    arlpy.lat_end = lat_end
    arlpy.range_points = len(pyram_ranges)
    arlpy.depth_points = len(pyram_depths)
    arlpy.num_points = num_points
    arlpy.bty, arlpy.max_range, arlpy.max_depth = [[bty_ranges[i], bty_depths[i]] for i in range(len(bty_ranges))], max(bty_ranges), max(bty_depths)
    arlpy.bty, arlpy.max_range, arlpy.max_depth = arlpy.read_bty()
    arlpy.angles = np.arange(angle_min, angle_max + 1, 0.1)
    arlpy.num_beams = len(arlpy.angles)

    env = arlpy.coherent_tl()
    # print(env)
    tloss = pm.compute_transmission_loss(env)

    arlpl = arlpy.source_level + 20*np.log10(np.abs(tloss.values) + 1e-12)
    arlpy_ranges = env['rx_range'] / 1000.0
    arlpy_depths = env['rx_depth']

    # FIGURE 
    fig = plt.figure(figsize=(18,8))
    gs = fig.add_gridspec(2, 2, width_ratios=[2, 1])

    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0], sharex=ax1)
    ax3 = fig.add_subplot(gs[:,1])  # right side (spans both rows)

    cmap = plt.get_cmap('viridis').copy()
    cmap.set_bad('k')

    # PYRAM FIELD 
    im1 = ax1.imshow(pyram_pl,
                     extent=[pyram_ranges.min(), pyram_ranges.max(),
                             pyram_depths.max(), pyram_depths.min()],
                     aspect='auto', origin='upper',
                     cmap=cmap, vmin=100, vmax=150)
    ax1.set_ylabel("Depth (m)")
    ax1.set_title("PYRAM Result")
    plt.colorbar(im1, ax=ax1, label="Pressure Level (dB)")

    # BELLHOP FIELD 
    im2 = ax2.imshow(arlpl,
                     extent=[arlpy_ranges.min(), arlpy_ranges.max(),
                             arlpy_depths.max(), arlpy_depths.min()],
                     aspect='auto', origin='upper',
                     cmap=cmap, vmin=100, vmax=150)
    ax2.set_xlabel("Range (km)")
    ax2.set_ylabel("Depth (m)")
    ax2.set_title("Bellhop Result")
    plt.colorbar(im2, ax=ax2, label="Pressure Level (dB)")

    # Depth at Final Range
    pyram_slice = pyram_pl[:, -1]
    arl_slice = arlpl[:, -1]

    ax3.plot(pyram_slice, pyram_depths, label="PYRAM", linewidth=2)
    ax3.plot(arl_slice, arlpy_depths, label="Bellhop", linewidth=2)
    ax3.invert_yaxis()
    ax3.set_xlabel("Pressure Level (dB)")
    ax3.set_ylabel("Depth (m)")
    ax3.set_title("Final Range Depth Slice")
    ax3.set_xlim(100,140)
    ax3.legend()

    plt.tight_layout()
    plt.savefig('pyram_vs_bellhop.png', dpi=300)
    plt.show()