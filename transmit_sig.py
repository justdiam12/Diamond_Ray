import os
import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import chirp, windows

def generate_arms_waveform(date):
    sample_freq = 80000

    if date == 21 or date == 22 or date == 25:
        signal, time = gen_cw(sample_freq)
    elif date == 23 or date == 24 or date == 26 or date == 29:
        signal, time = gen_lfm_100(sample_freq)
    elif date == 27 or date == 28:
        signal, time = gen_lfm_1000(sample_freq)
    else:
        signal, time = gen_lfm_3000(sample_freq)

    return signal, time


def plot_waveform(date, save_file_dir):

    signal, time = generate_arms_waveform(date)

    if date == 21 or date == 22 or date == 25:
        title = "Continuous Waveform Pulse"
        save_file = "cw.png"
    elif date == 23 or date == 24 or date == 26 or date == 29:
        title = "Linear Frequency Modulated (100 Hz)"
        save_file = "lfm100.png"
    elif date == 27 or date == 28:
        title = "Linear Frequency Modulated (1000 Hz)"
        save_file = "lfm1000.png"
    else:
        title = "Linear Frequency Modulated (3000 Hz)"
        save_file = "lfm3000.png"

    save = os.path.join(save_file_dir, save_file)
    plt.figure(figsize=(12,6))
    plt.plot(time, signal)
    plt.xlabel("Time")
    plt.ylabel("Signal")
    plt.title(title)
    plt.savefig(save)
    plt.show()



def gen_cw(sample_freq, center_freq=3500, duration=1.0, percent_taper=10):
    fs = sample_freq
    T = duration

    t = np.arange(0, T, 1/fs)

    pulse = np.cos(2 * np.pi * center_freq * t)

    alpha = percent_taper / 100
    window = windows.tukey(len(pulse), alpha)

    return pulse * window, t


def gen_lfm_100(sample_freq, center_freq=3500, duration=1.0, percent_taper=10):
    fs = sample_freq
    T = duration

    f0 = center_freq - 50
    f1 = center_freq + 50

    t = np.arange(0, T, 1/fs)

    pulse = chirp(t, f0=f0, f1=f1, t1=T, method='linear')

    alpha = percent_taper / 100
    window = windows.tukey(len(pulse), alpha)

    return pulse * window, t


def gen_lfm_1000(sample_freq, center_freq=3500, duration=1.0, percent_taper=10):
    fs = sample_freq
    T = duration

    f0 = center_freq - 500
    f1 = center_freq + 500

    t = np.arange(0, T, 1/fs)

    pulse = chirp(t, f0=f0, f1=f1, t1=T, method='linear')

    alpha = percent_taper / 100
    window = windows.tukey(len(pulse), alpha)

    return pulse * window, t


def gen_lfm_3000(sample_freq, center_freq=3500, duration=1.0, percent_taper=10):
    fs = sample_freq
    T = duration

    f0 = center_freq - 1500
    f1 = center_freq + 1500

    t = np.arange(0, T, 1/fs)

    pulse = chirp(t, f0=f0, f1=f1, t1=T, method='linear')

    alpha = percent_taper / 100
    window = windows.tukey(len(pulse), alpha)

    return pulse * window, t


# if __name__ == "__main__":
#     save_file_dir = '/Users/justindiamond/Documents/Documents/UW-APL/Research/Diamond_Ray/Signal Plots'
#     # date = 21
#     # date = 23
#     # date = 27
#     date = 30
    
#     plot_waveform(date, save_file_dir)
