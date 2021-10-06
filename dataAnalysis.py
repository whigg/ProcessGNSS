# DEREK PICKELL
# DATA ANALYSIS RTKLIB .POS FILES, ECEF

#import scipy
#import numpy
import statistics
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import sys, os, glob
import pymap3d as pm #

NUM_OBSERVATIONS = 6700 #at 15s sample rate, approx 1 day

def parse_file(file_name):
    x = []
    y = []
    z = []
    t = []
    with open(file_name) as f:
        for i in range(15):
            next(f)
        for line in f:
            values = line.split()
            x.append(float(values[2]))
            y.append(float(values[3]))
            z.append(float(values[4]))
            t.append(float(values[1]))

    return x, y, z, t
        
def det_sample(t0, x0, y0, z0):
    """
    determines sample rate of data, if not 15s, then decimate
    also curtails array to fixed time
    """
    if ((t0[2]-t0[1])==1):
        x1 = x0[0::15]
        y1 = y0[0::15]
        z1 = z0[0::15]
        t1 = t0[0::15]
        print("SAMPLE RATE 1Hz: Decimating")
    else:
        x1 = x0
        y1 = y0
        z1 = z0
        print("SAMPLE RATE 1/15s")

    #curtail array length NUM_OBSERVATIONS
    x2 = x1[:NUM_OBSERVATIONS]
    y2 = y1[:NUM_OBSERVATIONS]
    z2 = z1[:NUM_OBSERVATIONS]
    t2 = t1[:NUM_OBSERVATIONS]

    return x2, y2, z2, t2

def calc_statistics(x, y, z, t):
    """
    What about # of satellites? 
    Ensure sample #/sample rate is same for Trimble Data...
    """
    meanx = statistics.mean(x) 
    stdx = statistics.stdev(x)
    meany = statistics.mean(y) 
    stdy = statistics.stdev(y)
    meanz = statistics.mean(z) 
    stdz = statistics.stdev(z)
    # print ("x mean (m): ", meanx)
    # print ("x stdev (m): ", stdx)
    # print ("x sample size: ", len(x))
    # print ("z mean (m): ", meanz)
    # print ("z stdev (m): ", stdz)
    # print ("z sample size: ", len(z))

    #ploting
    fig, ax = plt.subplots(3, sharex=True)
    fig.suptitle("ECEF positions")

    ax[0].plot(x)
    ax[0].set_title(f"x position (m), mean: {meanx}, stdev: {stdx}", fontsize = 8)

    ax[1].plot(y)
    ax[1].set_title(f"y position (m), mean: {meany}, stdev: {stdy}", fontsize = 8)

    ax[2].plot(z)
    ax[2].set_title(f"z position (m), mean: {meanz}, stdev: {stdz}", fontsize = 8)
    ax[2].set(xlabel="GPS date in seconds")
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

    #ploting
    # maxilat = max(lattitude)
    # minilat = min(lattitude)
    # maxilon = max(longitude)
    # minilon = min(longitude)

    # print("Number of Observations: ", len(lattitude))
    # plt.scatter(longitude, lattitude, alpha=0.1, edgecolors='none')
    # plt.xlim([maxilon, minilon])
    # plt.ylim([minilat, maxilat])
    # plt.xlabel("Longitude E/W")
    # plt.ylabel("Lattitude N/S")
    # plt.title("Lat/Lon Scatter")
    # plt.show()

    return lattitude, longitude

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
    fig.suptitle("LL")

    maxilat = max(lat)
    minilat = min(lat)
    maxilon = max(lon)
    minilon = min(lon)
    
    ax.scatter(lon, lat, alpha=0.1, edgecolors='None')
    ax.scatter(meanlon, meanlat, color='red')
    ax.set_xlim([maxilon, minilon])
    ax.set_ylim([minilat, maxilat])
    circle = Circle((meanlon, meanlat), rms2d, facecolor = 'None', edgecolor='r', zorder=10)
    ax.add_patch(circle)
    plt.show()


if __name__ == '__main__':
    print("PULLING .pos file from /tempData folder...")
    script_dir = os.path.dirname(__file__)
    path = os.path.join(script_dir, "tempData", "*.pos") 
    file_name = glob.glob(path)
    print("FILE FOUND: ", file_name)

    #Calculations:
    x, y, z, t = parse_file(file_name[0])
    x1, y1, z1, t1 = det_sample(t, x, y, z) 
    #calc_statistics(x1, y1, z1, t1)
    lat, lon = ecef_lla(x1, y1, z1)
    lla_plot(lat, lon)
    
     


