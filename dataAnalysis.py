# DEREK PICKELL
# DATA ANALYSIS RTKLIB .POS FILES, ECEF

#import scipy
#import numpy
import statistics
import matplotlib.pyplot as plt
import sys, os, glob

NUM_OBSERVATIONS = 5700 #at 15s sample rate, approx 1 day

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
        
def det_sample(t, x0, y0, z0):
    """
    determines sample rate of data, if not 15s, then decimate
    also curtails array to fixed time
    """
    if ((t[2]-t[1])==1):
        x1 = x0[0::15]
        y1 = y0[0::15]
        z1 = z0[0::15]
        print("Sample Rate 1Hz: Decimating")
    else:
        x1 = x0
        y1 = y0
        z1 = z0
        print("Sample Rate 15s")

    #curtail array length NUM_OBSERVATIONS
    x2 = x1[:NUM_OBSERVATIONS]
    y2 = y1[:NUM_OBSERVATIONS]
    z2 = z1[:NUM_OBSERVATIONS]

    return x2, y2, z2

def calc_statistics(x, y, z):
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
    print ("x mean (m): ", meanx)
    print ("x stdev (m): ", stdx)
    print ("x sample size: ", len(x))
    print ("z mean (m): ", meanz)
    print ("z stdev (m): ", stdz)
    print ("z sample size: ", len(z))

    #ploting
    fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True)
    fig.suptitle("ECEF positions")

    ax1.plot(x)
    ax1.set_title(f"x position (m), mean: {meanx}, stdev: {stdx}", fontsize = 8)
    #ax1.text(3, 8, "test", fontsize=15, color='red')

    ax2.plot(y)
    ax2.set_title(f"y position (m), mean: {meany}, stdev: {stdy}", fontsize = 8)

    ax3.plot(z)
    ax3.set_title(f"z position (m), mean: {meanz}, stdev: {stdz}", fontsize = 8)
    
    plt.show()


if __name__ == '__main__':
    print("Pullling .pos file from /tempData folder...")
    script_dir = os.path.dirname(__file__)
    path = os.path.join(script_dir, "tempData", "*.pos") 
    file_name = glob.glob(path)

    #Calculations:
    x, y, z, t = parse_file(file_name[0])
    x1, y1, z1 = det_sample(t, x, y, z) 
    calc_statistics(x1, y1, z1)
     


