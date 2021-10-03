# DEREK PICKELL
# DATA ANALYSIS RTKLIB .POS FILES, ECEF

#import scipy
#import numpy
import statistics
import matplotlib.pyplot as plt
import sys, os, glob

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

    return x, y, z
        
def plot_xyz(x, y, z):
    xplot = plt.figure(1)
    plt.plot(x)
    plt.title("x position (m)")

    yplot = plt.figure(2)
    plt.plot(y)
    plt.title("y position (m)")

    zplot = plt.figure(3)
    plt.plot(z)
    plt.title("z position (m)")
    
    plt.show()

def calc_statistics(x, y, z):
    """
    What about # of satellites? 
    Ensure sample #/sample rate is same for Trimble Data...
    """
    meanx = statistics.mean(x) 
    stdx = statistics.stdev(x)
    print ("x mean (m): ", meanx)
    print ("x stdev (m): ", stdx)
    print ("x sample size: ", len(x))



if __name__ == '__main__':
     script_dir = os.path.dirname(__file__)
     path = os.path.join(script_dir, "tempData", "*.pos") 
     file_name = glob.glob(path)
     x, y, z = parse_file(file_name[0])
     
     calc_statistics(x, y, z)
     #plot_xyz(x, y, z)
     


