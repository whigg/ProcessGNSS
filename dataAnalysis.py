# DEREK PICKELL
# DATA ANALYSIS RTKLIB .POS FILES, ECEF

#import scipy
import numpy as np
import statistics
import sys, os, glob
import pymap3d as pm 

####### MATPLOTLIB PARAMETERS
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Circle
plt.rc('font', family='serif')
plt.rc('xtick', labelsize=6)
plt.rc('ytick', labelsize=6)
########


def parse_file(file_name):
    """
    INPUT: .obs file with ECEF positions (x, y, z, t)
    """
    x = []
    y = []
    z = []
    t = []
    with open(file_name) as f:
        for i in range(8): #skip first 8 lines of header
            next(f)
        for line in f:
            values = line.split()
            x.append(float(values[2]))
            y.append(float(values[3]))
            z.append(float(values[4]))
            t.append(float(values[1]))

    return x, y, z, t
        
def det_sample15(x0, y0, z0, t0):
    """
    determines sample rate of data, if not 15s, then decimate
    also curtails array to fixed time NUM_OBSERVATIONS
    """
    if ((t0[2]-t0[1])==1):
        x1 = x0[0::15]
        y1 = y0[0::15]
        z1 = z0[0::15]
        t1 = t0[0::15]
        print("SAMPLE RATE 1Hz: Decimating")
    elif ((t0[2]-t0[1])==15):
        x1 = x0
        y1 = y0
        z1 = z0
        print("SAMPLE RATE 1/15s")
    else:
        print("unable to decimate: unknown sample rate")

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
    meanx_av = [meanx] * len(x)
    stdx = statistics.stdev(x)

    meany = statistics.mean(y) 
    stdy = statistics.stdev(y)
    meany_av = [meany] * len(y)

    meanz = statistics.mean(z) 
    stdz = statistics.stdev(z)
    meanz_av = [meanz] * len(z)

    #ploting
    fig, ax = plt.subplots(3, sharex=True, figsize=(10,6))
    fig.suptitle(f"ECEF Positions")
    ylim_range = 10

    ax[0].plot(x, lw=.5)
    ax[0].plot(meanx_av, color='red', lw=.5, ls='--', label="x average: "+'%.11g' % meanx +"m")
    ax[0].set_title("x position (m), stdev: "'%.5g' % stdx, fontsize = 10)
    ax[0].legend(loc='lower left', fontsize = 8)
    ax[0].set_ylim(float('%.11g' % meanx) - ylim_range, float('%.11g' % meanx) + ylim_range)
    
    ax[1].plot(y, lw=.5)
    ax[1].plot(meany_av, color='red', lw=.5, ls='--', label="y average: "+'%.11g' % meany +"m")
    ax[1].set_title("y position (m), stdev: "'%.5g' % stdy, fontsize = 10)
    ax[1].legend(loc='lower left', fontsize = 8)
    ax[1].set_ylim(float('%.11g' % meany) - ylim_range, float('%.11g' % meany) + ylim_range)

    ax[2].plot(z, lw=.5)
    ax[2].plot(meanz_av, color='red', lw=.5, ls='--', label="z average: "+'%.11g' % meanz +"m")
    ax[2].set_title("z position (m), stdev: "'%.5g' % stdz, fontsize = 10)
    ax[2].legend(loc='lower left', fontsize = 8)
    ax[2].set(xlabel="Time (s)")
    ax[2].set_ylim(float('%.11g' % meanz) - ylim_range, float('%.11g' % meanz) + ylim_range)

    #for i in [0, 1, 2]:
        #ax[i].yaxis.set_major_locator(ticker.MultipleLocator(2))
        #ax[i].yaxis.set_minor_locator(ticker.MultipleLocator(.5))
    
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
    circle = Circle((meanlon, meanlat), rms2d, facecolor = 'None', lw = .5, ls = '--', edgecolor='r', zorder=10, label="2D RMS")
    ax.add_patch(circle)
    ax.legend(loc='lower left', fontsize = 8)
    plt.show()


if __name__ == '__main__':
    print("PULLING .pos file from /tempData folder...")
    script_dir = os.path.dirname(__file__)
    path = os.path.join(script_dir, "tempData", "*.pos") 
    file_name = glob.glob(path)
    print("FILE FOUND: ", file_name)

    # PARSE FILE
    x, y, z, t = parse_file(file_name[0])

    # ECEF PLOTS
    x1, y1, z1, t1 = det_sample15(x, y, z, t) 
    #calc_ECEFstatistics(x1, y1, z1, t1)

    # LAT LON  DISTRIBUTION
    lat, lon, alt = ecef_lla(x1, y1, z1)
    lla_plot(lat, lon)

    # HISTOGRAM
    plt.hist(lat, bins=300)
    plt.show()
    
     


