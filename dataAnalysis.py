# DEREK PICKELL
# DATA ANALYSIS RTKLIB .POS FILES, ECEF

#import scipy
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
########


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
    determines sample rate of data, if not 15s, then decimate
    also curtails array to fixed time NUM_OBSERVATIONS
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

def obs_time(obs_time, x, y, z, t):
    """
    shortens data to specific length of time
    """
    x2 = x1[:obs_time]
    y2 = y1[:obs_time]
    z2 = z1[:obs_time]
    t2 = t1[:obs_time]

def calc_ECEFstatistics(x, y, z, t):
    """
    What about # of satellites? 
    Ensure sample #/sample rate is same for Trimble Data...
    """

    meanx = statistics.mean(x) 
    meanx_line = [meanx] * len(x)
    stdx = statistics.pstdev(x)
    stdx3 = 3*stdx
    stdx_line = [meanx+stdx*3] *len(x)
    stdx_line2 = [meanx+stdx*-3] *len(x)

    meany = statistics.mean(y) 
    stdy = statistics.pstdev(y)
    meany_line = [meany] * len(y)
    stdy3 = 3*stdy
    stdy_line = [meany+stdy*3] *len(y)
    stdy_line2 = [meany+stdy*-3] *len(y)

    meanz = statistics.mean(z) 
    stdz = statistics.pstdev(z)
    meanz_line = [meanz] * len(z)
    stdz3 = 3*stdz
    stdz_line = [meanz+stdz*3] *len(z)
    stdz_line2 = [meanz+stdz*-3] *len(z)


    #ploting
    fig, ((ax1, ax4), (ax2, ax5), (ax3, ax6)) = plt.subplots(3, 2, figsize=(12,8), gridspec_kw={'width_ratios': [3, 1]}) #, sharex=True
    fig.suptitle(f"ECEF Positions from " +str(t[0]) + " to "+ str(t[-1]))
    lim_range = 10

    ax1.plot_date(t, x, marker='None', ls='-', lw=plot_lw)
    ax1.plot_date(t, meanx_line, marker='None',color='red', lw=plot_lw, ls='--', label="x average: "+'%.11g' % meanx +"m")
    ax1.plot_date(t, stdx_line, marker='None',color='green', lw=plot_lw, ls='--', label="3\u03C3: +/-"+ '%.5g' % stdx3 +"m")
    ax1.plot_date(t, stdx_line2, marker='None', color='green', lw=plot_lw, ls='--')
    ax1.set_title("x position (m)", fontsize = title_size) #, stdev: "'%.5g' % stdx
    ax1.legend(loc='lower left', fontsize = label_size)
    ax1.set_ylim(float('%.11g' % meanx) - lim_range, float('%.11g' % meanx) + lim_range)
    ax1.margins(x=0)
    ax1.set_xticklabels([])
    
    ax2.plot_date(t, y, marker='None', ls='-', lw=plot_lw)
    ax2.plot_date(t, meany_line, marker='None',color='red', lw=plot_lw, ls='--', label="y average: "+'%.11g' % meany +"m")
    ax2.plot_date(t, stdy_line, marker='None',color='green', lw=plot_lw, ls='--', label="3\u03C3: +/-"+ '%.5g' % stdy3 +"m")
    ax2.plot_date(t, stdy_line2, marker='None',color='green', lw=plot_lw, ls='--')
    ax2.set_title("y position (m)", fontsize = title_size)
    ax2.legend(loc='lower left', fontsize = label_size)
    ax2.set_ylim(float('%.11g' % meany) - lim_range, float('%.11g' % meany) + lim_range)
    ax2.margins(x=0)
    ax2.set_xticklabels([])

    ax3.plot_date(t, z, marker='None', ls='-', lw=plot_lw)
    ax3.plot_date(t, meanz_line, marker='None',color='red', lw=plot_lw, ls='--', label="z average: "+'%.11g' % meanz +"m")
    ax3.plot_date(t, stdz_line, marker='None',color='green', lw=plot_lw, ls='--', label="3\u03C3: +/-"+ '%.5g' % stdz3 +"m")
    ax3.plot_date(t, stdz_line2, marker='None', color='green', lw=plot_lw, ls='--')
    ax3.set_title("z position (m)", fontsize = title_size)
    ax3.legend(loc='lower left', fontsize = label_size)

    ax3.set_xticklabels(ax3.get_xticks(), rotation = 45)
    ax3.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d %H:%M:%S"))
    ax3.set(xlabel="Time")
    ax3.set_ylim(float('%.11g' % meanz) - lim_range, float('%.11g' % meanz) + lim_range)
    ax3.margins(x=0)

    #histograms
    ax4.hist(x, bins=300, align='left')
    ax4.axvline(meanx, color='red', lw=plot_lw, ls='--')
    ax4.axvline(stdx_line[0], color='green', lw=plot_lw, ls='--')
    ax4.axvline(stdx_line2[0], color='green', lw=plot_lw, ls='--')
    ax4.set_xlim(float('%.11g' % meanx) - lim_range, float('%.11g' % meanx) + lim_range)
    ax5.hist(y, bins=300)
    ax5.axvline(meany, color='red', lw=plot_lw, ls='--')
    ax5.axvline(stdy_line[0], color='green', lw=plot_lw, ls='--')
    ax5.axvline(stdy_line2[0], color='green', lw=plot_lw, ls='--')
    ax5.set_xlim(float('%.11g' % meany) - lim_range, float('%.11g' % meany) + lim_range)
    ax6.hist(z, bins=300)
    ax6.axvline(meanz, color='red', lw=plot_lw, ls='--')
    ax6.axvline(stdz_line[0], color='green', lw=plot_lw, ls='--')
    ax6.axvline(stdz_line2[0], color='green', lw=plot_lw, ls='--')
    ax6.set_xlim(float('%.11g' % meanz) - lim_range, float('%.11g' % meanz) + lim_range)

    plt.margins(x=0)
    plt.show()

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

def calc_ECEFstatisticsPPK(x, y, z, t, x1, y1, z1, t1):
    """
    x, y, z: ECEF of unprocessed file
    x1, y1, z1: PPK/RTK data (will be plotted on top of original)
    """
    meanx = statistics.mean(x) 
    meanx_line = [meanx] * len(x)
    stdx = statistics.pstdev(x)
    stdx3 = 3*stdx
    stdx_line = [meanx+stdx*3] *len(x)
    stdx_line2 = [meanx+stdx*-3] *len(x)

    meany = statistics.mean(y) 
    stdy = statistics.pstdev(y)
    meany_line = [meany] * len(y)
    stdy3 = 3*stdy
    stdy_line = [meany+stdy*3] *len(y)
    stdy_line2 = [meany+stdy*-3] *len(y)

    meanz = statistics.mean(z) 
    stdz = statistics.pstdev(z)
    meanz_line = [meanz] * len(z)
    stdz3 = 3*stdz
    stdz_line = [meanz+stdz*3] *len(z)
    stdz_line2 = [meanz+stdz*-3] *len(z)

    #ploting
    fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(12,8)) #, sharex=True
    fig.suptitle(f"ECEF Positions")
    lim_range = 10

    ax1.plot(x, lw=plot_lw, alpha=.5, label="Unprocessed Data")
    ax1.plot(x1, lw=2, color='green', label="PPK Data")
    ax1.plot(meanx_line, color='red', lw=plot_lw, ls='--', label="x average: "+'%.11g' % meanx +"m")
    ax1.set_title("x position (m)", fontsize = title_size) #, stdev: "'%.5g' % stdx
    ax1.legend(loc='lower left', fontsize = label_size)
    ax1.set_ylim(float('%.11g' % meanx) - lim_range, float('%.11g' % meanx) + lim_range)
    ax1.margins(x=0)
    
    ax2.plot(y, lw=plot_lw, alpha=.5, label="Unprocessed Data")
    ax2.plot(y1, lw=2, color='green', label="PPK Data")
    ax2.plot(meany_line, color='red', lw=plot_lw, ls='--', label="y average: "+'%.11g' % meany +"m")
    ax2.set_title("y position (m)", fontsize = title_size)
    ax2.legend(loc='lower left', fontsize = label_size)
    ax2.set_ylim(float('%.11g' % meany) - lim_range, float('%.11g' % meany) + lim_range)
    ax2.margins(x=0)

    ax3.plot(z, lw=plot_lw, alpha=.5, label="Unprocessed Data")
    ax3.plot(z1, lw=2, color='green', label="PPK Data")
    ax3.plot(meanz_line, color='red', lw=plot_lw, ls='--', label="z average: "+'%.11g' % meanz +"m")
    ax3.set_title("z position (m)", fontsize = title_size)
    ax3.legend(loc='lower left', fontsize = label_size)
    ax3.set(xlabel="Time (s)")
    ax3.set_ylim(float('%.11g' % meanz) - lim_range, float('%.11g' % meanz) + lim_range)
    ax3.margins(x=0)

    plt.margins(x=0)
    plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--xyz', help='Plot xyz coordinates on individual scatter plots')
    parser.add_argument('--lla', help='Plot lat vs lon')
    parser.add_argument('--histogram', help='Plot histogram of lattitude')
    parser.add_argument('--kinematic', type=int, help='Plots kinematic processed data')
    args = parser.parse_args()

    print("PULLING .pos file from /tempData folder...")
    script_dir = os.path.dirname(__file__)
    path = os.path.join(script_dir, "tempData", "*.pos") 
    file_name = glob.glob(path)
    print("FILE FOUND: ", file_name)

    # PARSE FILE(s)
    x, y, z, t = parse_file(file_name[0])

    if args.xyz is not None:
        # ECEF PLOTS
        x1, y1, z1, t1 = det_sample15(x, y, z, t) 
        calc_ECEFstatistics(x1, y1, z1, t1)

    if args.lla is not None:
        # LAT LON  DISTRIBUTION
        lat, lon, alt = ecef_lla(x1, y1, z1)
        lla_plot(lat, lon)

    if args.histogram is not None:
        # HISTOGRAM
        plt.hist(lat, bins=300)
        plt.show()

    if args.kinematic is not None: 
        print("PPK Solution: ")
        x1, y1, z1, t1 = parse_file(file_name[1]) #FIX THISS!!!!!!!!! 
        x1, y1, z1, t1 = det_sample(x1, y1, z1, t1)
        x, y, z, t = det_sample(x, y, z, t) 
        if args.kinematic == 1:
            calc_ECEFstatisticsPPK(x1, y1, z1, t1, x, y, z, t) # first ECEF's are unprocessed
        else:
            calc_ECEFstatisticsPPK(x, y, z, t, x1, y1, z1, t1)

    else:
        #default
        x2, y2, z2, t2 = det_sample(x, y, z, t, 1) 
        calc_ECEFstatistics(x2, y2, z2, t2)

    
    
     


