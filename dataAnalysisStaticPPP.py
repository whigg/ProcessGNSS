# Created by Derek Pickell
# 11/27/21
# Use for: 
# Plotting CSV PPP Data with Error Bars

import numpy as np 
import os
import matplotlib.pyplot as plt
import argparse

print("CSV must be in format Date, Height (m), 95%, Uncertainty (m) ... ")
def dir_path(string):
    if os.path.isfile(string):
        return string
    else:
        raise NotADirectoryError(string)

parser = argparse.ArgumentParser()
parser.add_argument("file_1_path", type=dir_path, help="full file path to first data file")
args = parser.parse_args()

file_path1 = os.path.join(args.file_1_path)
file1_name = os.path.basename(file_path1)

print("Parsing CSV")
data = np.genfromtxt(file_path1, delimiter=",", skip_header=1).T
heights = data[1]
uncertainties = data[2] 
increment = np.arange(1, len(heights)+1)

######## PPP Static 24 Hours 
fig, ax = plt.subplots()
ax.set_title("Static Solution Time Series (24hr solution each)", fontsize=14, fontname='Baskerville', fontweight='bold')
ax.set_ylabel("ITRF Height (m) Referenced to WGS84", fontsize=12, fontname='Baskerville')
ax.scatter(increment, heights)
ax.axhline(np.median(heights), color='red', lw=1, ls='--', label=f'median elevation {np.mean(heights):.3f}m')
ax.ticklabel_format(useOffset=False)
plt.errorbar(increment, heights, yerr=uncertainties, fmt='o', color='xkcd:navy')
plt.legend()
fig.show()

####### GNSS IR 
IR_heights = [1.59,1.645,1.615,1.705,1.58,1.631,1.63,1.68,1.766,1.63,1.645,1.68,1.686,1.655,1.68,1.66,1.765,1.68,1.641,1.75,1.59,1.665,1.715,1.63,1.58,1.615,1.66,1.621,1.68,1.605,1.635,1.585,1.601,1.671,1.585,1.71,1.621,1.65,1.73,1.68,1.75,1.665,1.631,1.665,1.615,1.645,1.746,1.6,1.606,1.62,1.691,1.675,1.721,1.585,1.625,1.635,1.685,1.595,1.65,1.615,1.6,1.615,1.58,1.615,1.6,1.67,1.745,1.581,1.696,1.586,1.671,1.651,1.64,1.66,1.705,1.59,1.625,1.665,1.72,1.63,1.646,1.56,1.685,1.625,1.645,1.59,1.666,1.63,1.7,1.751,1.636,1.665,1.736,1.606,1.685,1.71,1.675,1.75,1.665,1.675,1.745,1.625,1.7,1.741,1.59,1.615,1.64,1.645,1.655,1.69,1.6,1.675,1.58,1.59,1.66,1.595,1.681,1.74,1.64,1.635,1.665,1.675,1.67,1.73,1.635,1.606,1.665,1.615,1.605,1.655,1.586,1.605,1.655,1.651,1.745,1.615,1.66,1.611,1.681,1.715,1.65,1.605,1.67,1.775,1.63,1.685,1.591,1.681]
mean_IR_height = np.mean(IR_heights)
fig2, ax2 = plt.subplots()
ax2.hist(np.asarray(IR_heights), bins=20, histtype='bar', color='xkcd:navy', orientation='horizontal') 
ax2.set_title("Histogram of Antenna Height Estimates", fontsize=14, fontname='Baskerville', fontweight='bold')
ax2.set_ylabel("Estimated Height (m)", fontsize=12, fontname='Baskerville')
ax2.set_xlabel("Counts", fontsize=12, fontname='Baskerville')
ax2.axhline(mean_IR_height, color='red', lw=3, ls='-', label=f'mean height {mean_IR_height:.1f}m')
ax2.axhline(mean_IR_height + np.std(IR_heights), color='red', lw=1, ls='--', label=f'1 \u03C3 {np.std(IR_heights):.2f}m')
ax2.axhline(mean_IR_height - np.std(IR_heights), color='red', lw=1, ls='--')

plt.legend()
fig2.show()

plt.show()