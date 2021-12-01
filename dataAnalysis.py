# DEREK PICKELL
# DATA ANALYSIS RTKLIB .POS FILES, ECEF
# INPUTS: /tempData requires at least 1 .pos file
# .POS Must be in ECEF, with time in YYYY/MM/DD HH:MM:SS.SSSS

import numpy as np
import statistics
import sys, os, glob
import argparse
import pymap3d as pm 
import datetime as dt

####### MATPLOTLIB PARAMETERS
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.dates as mdates
from matplotlib.patches import Circle
plt.rc('font', family='serif')
plt.rc('xtick', labelsize=6)
plt.rc('ytick', labelsize=6)
plot_lw = 0.5
title_size = 10
label_size = 8
lim_range = 10
########

########HELPER FUNCTIONS
def parse_file(file_name):
    """
    INPUT: .obs file with ECEF positions (x, y, z, t)
    """
    x = []
    y = []
    z = []
    t1 = []
    t2 = []
    with open(file_name) as f:
        for i in range(30): #skip first 8 lines of header
            next(f)
        for line in f:
            values = line.split()
            x.append(float(values[2]))
            y.append(float(values[3]))
            z.append(float(values[4]))
            t1.append(values[0])
            t2.append((values[1])) #t.append(float(values[1]))
    t3 = [a+":"+b for a, b in zip(t1, t2)]
    tfinal = [dt.datetime.strptime(d, '%Y/%m/%d:%H:%M:%S.%f') for d in t3]

    return x, y, z, tfinal
        
def det_sample(x0, y0, z0, t0, case=1): #time=86400
    """
    If sample rate is 1Hz, decide whether to decimate to 15s or maintain
    """
    if case==15: 
            x1 = x0[0::15]
            y1 = y0[0::15]
            z1 = z0[0::15]
            t1 = t0[0::15]
            print("Decimating from 1Hz to 15Hz")
    elif case==1:
            x1 = x0
            y1 = y0
            z1 = z0
            t1 = t0  
    
    # x1 = x1[0:time]
    # y1 = y1[0:time]
    # z1 = z1[0:time]
    #t1 = t1[0:time]
  
    return x1, y1, z1, t1

def ecef_lla(x, y, z):
    """
    """
    lattitude = []
    longitude = []
    altitude = []
    for i in range(len(x)):
        lat, lon, alt = pm.ecef2geodetic(x[i], y[i], z[i])
        lattitude.append(lat)
        longitude.append(lon)
        altitude.append(alt)

    return lattitude, longitude, altitude

def ecef_w_lines_plot(t, y, mean, mean_line, sd_line1=None, sd_line2=None, stdev=None, ax=None, **kwargs):
    if ax is None:
        ax.plt.gca()

    ax.plot_date(t, y, marker='None', ls='-', lw=plot_lw, **kwargs)
    ax.plot_date(t, mean_line, marker='None',color='red', lw=plot_lw, ls='--', label="average: "+'%.11g' % mean +"m")
    ax.legend(loc='lower left', fontsize = label_size)
    ax.set_ylim(float('%.11g' % mean) - lim_range, float('%.11g' % mean) + lim_range)
    ax.margins(x=0)

    if sd_line1 is not None:
        ax.plot_date(t, sd_line1, marker='None',color='green', lw=plot_lw, ls='--', label="3\u03C3: +/-"+ '%.5g' % stdev +"m")
        ax.plot_date(t, sd_line2, marker='None', color='green', lw=plot_lw, ls='--')

    return(ax)

def histogram_w_stdev_plot(data, mean, sd_line1, sd_line2, ax=None):
    if ax is None:
        ax.plt.gca()

    ax.hist(data, bins=300)
    ax.axvline(mean, color='red', lw=plot_lw, ls='--')
    ax.axvline(sd_line1[0], color='green', lw=plot_lw, ls='--')
    ax.axvline(sd_line2[0], color='green', lw=plot_lw, ls='--')
    ax.set_xlim(float('%.11g' % mean) - lim_range, float('%.11g' % mean) + lim_range)

    return(ax)

##########DATA ANALYSIS
def ecef_and_histogram(x, y, z, t):
    """
    """
    fig, ((ax1, ax4), (ax2, ax5), (ax3, ax6)) = plt.subplots(3, 2, figsize=(12,8), gridspec_kw={'width_ratios': [3, 1]}) #, sharex=True
    fig.suptitle(f"ECEF Positions from " +str(t[0]) + " to "+ str(t[-1]))

    l = [x, y, z]
    ax = [ax1, ax2, ax3]
    ax_hist = [ax4, ax5, ax6]
    for i in range(3):
        mean = statistics.mean(l[i])
        mean_line = [mean] * len(l[i])
        stdev = statistics.pstdev(l[i])
        stdev3 = 3*stdev
        std_line = [mean+stdev3] *len(l[i])
        std_line2 = [mean-stdev3] *len(l[i])
        #plot ECEF
        ecef_w_lines_plot(t, l[i], mean, mean_line, std_line, std_line2, stdev3, ax=ax[i])
        #plot Histogram
        histogram_w_stdev_plot(l[i], mean, std_line, std_line2, ax_hist[i])
        
    ax1.set_title("x position (m)", fontsize = title_size) 
    ax1.set_xticklabels([])
    ax2.set_title("y position (m)", fontsize = title_size) 
    ax2.set_xticklabels([])
    ax3.set_title("z position (m)", fontsize = title_size) 
    ax3.set_xticklabels([])
    ax3.set_xticklabels(ax3.get_xticks(), rotation = 45)
    ax3.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d %H:%M:%S"))
    ax3.set(xlabel="Time")

    plt.margins(x=0)
    plt.show()

def lla_plot(lat, lon):
    """
    """
    meanlat = statistics.mean(lat)
    meanlon = statistics.mean(lon)

    deltalat = [x-meanlat for x in lat]
    deltalon = [x-meanlon for x in lon]

    #standard dev
    sigmalat = statistics.stdev(deltalat)
    sigmalon = statistics.stdev(deltalon)

    #CEP 50% radius
    cep = 0.59*(sigmalat + sigmalon)
    print("CEP 50%: ", cep)

    #2DRMS (95%)
    rms2d = 2*((sigmalat**2 + sigmalon**2)**(1/2))
    print("2DRMS: ", rms2d)

    #Plot
    fig, ax = plt.subplots(1)
    fig.suptitle("Lat/Lon Distribution")

    maxilat = max(lat)
    minilat = min(lat)
    maxilon = max(lon)
    minilon = min(lon)
    
    ax.scatter(lon, lat, alpha=0.1, edgecolors='None', color='tab:blue')
    ax.scatter(meanlon, meanlat, color='tab:red', lw=2, label="mean lattitude: "+'%.11g' % meanlat +" mean longitude: "+'%.11g' % meanlon)
    ax.set_xlim([maxilon, minilon])
    ax.set_ylim([minilat, maxilat])
    circle = Circle((meanlon, meanlat), rms2d, facecolor = 'None', lw = plot_lw, ls = '--', edgecolor='r', zorder=10, label="2D RMS")
    ax.add_patch(circle)
    ax.legend(loc='lower left', fontsize = 8)
    plt.show()

def ecefPPK_compared_ecefOriginal(x, y, z, t, x1, y1, z1, t1):
    """
    x, y, z: ECEF of unprocessed file
    x1, y1, z1: PPK/RTK data (will be plotted on top of original)
    """

    fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(12,8)) #, sharex=True
    fig.suptitle("ECEF Positions")

    l = [x, y, z]
    ax = [ax1, ax2, ax3]
    ax_hist = [ax4, ax5, ax6]
    for i in range(3):
        mean = statistics.mean(l[i])
        mean_line = [mean] * len(l[i])
        stdev = statistics.pstdev(l[i])
        stdev3 = 3*stdev
        std_line = [mean+stdev3] *len(l[i])
        std_line2 = [mean-stdev3] *len(l[i])
        #plot ECEF
        ecef_w_lines_plot(t, l[i], mean, mean_line, std_line=None, std_line2=None, stdev3=None, ax=ax[i], label="Unprocessed Data")

    ax1.plot(x1, lw=2, color='green', label="PPK Data")
    ax2.plot(y1, lw=2, color='green', label="PPK Data")
    ax3.plot(y1, lw=2, color='green', label="PPK Data")
    ax3.set(xlabel="Time")

    plt.margins(x=0)
    plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--xyz', help='Plot xyz coordinates on individual scatter plots')
    parser.add_argument('--compare', help='Plots kinematic processed data vs original')
    args = parser.parse_args()

    print("PULLING .pos file[s] from /tempData folder...")
    script_dir = os.path.dirname(__file__)
    path = os.path.join(script_dir, "tempData", "*.pos") 
    file_name = glob.glob(path)
    print("FILE FOUND: ", file_name)

    # PARSE FILE(s)
    x, y, z, t = parse_file(file_name[0])
    sample_rate = (t[1]-t[0]).total_seconds()
    print("Sample rate: " + str(sample_rate) + "s")

    if args.xyz is not None:
        # ECEF PLOTS
        ecef_and_histogram(x, y, z, t)

    if args.compare is not None: 
        print("Make sure ")
        num = input("Which file is the receiver .pos? [0, 1, ...]? ")
        x, y, z, t = parse_file(file_name[num]) 
        num2 = input("Which file is the base .pos? ")
        x1, y1, z1, t1 = parse_file(file_name[num2])
        ecefPPK_compared_ecefOriginal(x, y, z, t, x1, y1, z1, t1)

    else:
        #default
        x2, y2, z2, t2 = det_sample(x, y, z, t, 1) 
        ecef_and_histogram(x2, y2, z2, t2)

    
    
     


