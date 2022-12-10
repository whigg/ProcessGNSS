# Created by Derek Pickell
# 11/11/21
# Use for: 
# Comparison of two different GNSS instruments, Kinematic PPP (CSRS PROCESSED)
# Comparison of GNSS Kinematic PPP (CSRS) & ICESAT-2 (.txt, icepyx)
# Comparison of GNSS Static PPP (CSRS) & ICESAT-2    (.txt, icepyx)
# Comparison of PPK GNSS (RTKLIB .pos) and ICESAT-2  (.txt, icepyx)

import numpy as np 
import os
import matplotlib.pyplot as plt
import argparse
from geopy import distance
from sklearn.neighbors import BallTree

################## PARSE INPUTS ###########################
def dir_path(string):
    if os.path.isfile(string):
        return string
    else:
        raise NotADirectoryError(string)

parser = argparse.ArgumentParser()
parser.add_argument("file_1_path", type=dir_path, help="full file path to first data file")
parser.add_argument("type", choices = ['Kin_PPP', 'Stat_PPP', 'ICESat-2', 'PPK'], help="PPP is CSRS processed .csv, ICESat-2 is icepyx .txt, PPK is RTKLIB .pos")
parser.add_argument("file_2_path", type=dir_path, help="full file path to first data file")
parser.add_argument("type2", choices = ['Kin_PPP', 'Stat_PPP', 'ICESat-2', 'PPK'], help="PPP is CSRS processed .csv, ICESat-2 is icepyx .txt, PPK is RTKLIB .pos")

args = parser.parse_args()

file_path1 = os.path.join(args.file_1_path)
file_path2 = os.path.join(args.file_2_path)
file1_name = os.path.basename(file_path1)
file2_name = os.path.basename(file_path2)

if args.type == 'Kin_PPP':
    print("Parsing Kinematic PPP from CSRS")
    data = np.genfromtxt(file_path1, delimiter=',', skip_header=1, usecols=(0, 1, 2, 3)).T
    lattitudes1 = data[0]
    longitudes1 = data[1]
    ellipsoidal_heights1 = data[2]
    decimal_hour1 = data[3]

elif args.type == 'Stat_PPP':
    print("Static PPP")
elif args.type == 'ICESat-2':
    print("ICESat-2")
elif args.type =='PPK':
    """
    Parses Data from Emlid Studio RTKLIB, First 10 rows are metadata, lat/lon decimal degrees
    Filter for fixed solution only
    """
    print("Parsing PPK from RTKLIB")
    data = np.genfromtxt(file_path1, delimiter=None, skip_header=10, usecols=(1, 2, 3, 4, 5)).T
    fixed_flag1 = data[4]
    fixed = np.where(fixed_flag1 == 1, True, False)
    lattitudes1 = data[1][fixed]
    longitudes1 = data[2][fixed]
    ellipsoidal_heights1 = data[3][fixed]
    decimal_hour1 = data[0][fixed]  * 3600 - 162800

if args.type2 == 'Kin_PPP':
    print("Parsing Kinematic PPP from CSRS")
    data = np.genfromtxt(file_path2, delimiter=',', skip_header=1, usecols=(0, 1, 2, 3)).T
    lattitudes2 = data[0]
    longitudes2 = data[1]
    ellipsoidal_heights2 = data[2]
    decimal_hour2 = data[3]
elif args.type2 == 'Stat_PPP':
    print("Static PPP")
elif args.type2 == 'ICESat-2':
    print("ICESat-2")
elif args.type2 =='PPK':
    print("Parsing PPK from RTKLIB")
    data = np.genfromtxt(file_path2, delimiter=None, skip_header=10, usecols=(1, 2, 3, 4, 5)).T
    fixed_flag2 = data[4]
    fixed = np.where(fixed_flag2 == 1, True, False)
    lattitudes2 = data[1][fixed]
    longitudes2 = data[2][fixed]
    ellipsoidal_heights2 = data[3][fixed]
    decimal_hour2 = data[0][fixed] * 3600 - 162800

############################################################

def get_PSP_stats(lat, lon, h, plot=True):
    """
    finds pseudostatic points and computes standard deviation
    A pseudostatic point is defined as a cluser of points whose
    distance from each successive point is < "minimum_distance"
    Pseudostatic points are filtered such that they contain >10 points 
    but less than 50 points. This is to account for other idle moments..
    returns (1) mean of standard deviations
            (2) number of pseudostatic points
    """
    minimum_distance = 0.05 # in meters, minimum distance to be part of a stationary PSP
    pseudostatic_points = []
    location_PSP = []
    i = 0

    while i < len(lat) - 1:
        # get distance between first two points
        d = distance.distance((lat[i], lon[i]), (lat[i+1], lon[i+1])).meters # distance.distance calculates great circle distance of ellipsoid WGS84
        if d < minimum_distance: # potential cluster found
            temp_list = []
            temp_list.append(h[i]) # append first point
            j = i
            while distance.distance((lat[j], lon[j]), (lat[j+1], lon[j+1])).meters < minimum_distance and j < len(lat) -2:
                temp_list.append(h[j+1])
                j+=1

            if len(temp_list) > 10 and len(temp_list) < 50: # cluster found
                pseudostatic_points.append(temp_list)
                location_PSP.append([lat[j], lon[j]])
            i = j # jump outside of cluster

        i +=1

    number_PSP = len(pseudostatic_points)
    
    stdvs = []
    for i in range(0, number_PSP):
        stdvs.append(np.std(pseudostatic_points[i]))

    mean_1sigma = np.mean(stdvs)

    if plot:
        ## Verify Data ## 
        location_PSP_transformed = np.asarray(location_PSP).T
        plt.scatter(lon, lat, c='b')
        plt.scatter(location_PSP_transformed[1], location_PSP_transformed[0], c='r', s=20)
        plt.xlabel('Longitude', fontsize=12, fontweight='light', fontname='Baskerville')
        plt.ylabel('Latitude', fontsize=12, fontweight='light', fontname='Baskerville')
        plt.title(f"Location of PSPs (n={number_PSP})", fontsize=14, fontweight='bold', fontname='Baskerville')    
        plt.show()

    return number_PSP, mean_1sigma

def nearest_neighbor_compare(lat1, lon1, h1, lat2, lon2, h2, search_distance, bias1, bias2, plot=True):
    """
    - INPUTS: lat1/lon1/h1: single coordinate or lists; lat2/lon2/h2: list of coordinates; search_distance = 
        how close points need to be as "neighbors"
    - Find nearest neighbor between dataset 1 and dataset 2 points
    - Compute residuals & standard deviation of between these two points
    - RETURN: (1) residuals (2) standard deviation of residuals 
    *** NOTE: lat1/lon1/h1 is considered "truth"
    """
    # correct for antenna heights
    h1 = h1 - bias1
    h2 = h2 - bias2

    residuals = []
    residual_locations = []

    print("lat1 n#" , len(lat1))
    print("lat2 n#", len(lat2))

    array = np.asarray((lat1, lon1)).T
    array2 = np.asarray((lat2, lon2)).T

    tree = BallTree(array2, metric='euclidean')
    distances, indices = tree.query(array, k =1)

    lst = ["|","/","-","\\"]
    for i in range(0, len(lat1)):
        print(lst[i % 4], end="\r")
        
        if distance.distance((lat1[i], lon1[i]), (lat2[indices[i]], lon2[indices[i]])).meters < search_distance:
            residuals.append(h1[i] - h2[indices[i]])
            residual_locations.append([lat1[i], lon1[i]])

    number_residuals = len(residuals)
    res_sigma = np.std(residuals)
    res_mean = np.mean(residuals)
    res_median = np.median(residuals)
    new = [i for i in residuals if res_mean - 2* res_sigma<i<res_mean + 2*res_sigma] 
    new_sigma = np.std(new)
    print("filtered sigma: ", new_sigma)

    if plot:
        ## Verify Data ## 
        # plot where the residuals are, on top of original GPS points
        fig0, ax_a = plt.subplots()
        residual_locations_transformed = np.asarray(residual_locations).T
        ax_a.scatter(lon1, lat1, c='k', s=40, marker="o", label="dataset 1")
        ax_a.scatter(lon2, lat2, c='dimgray', s=20, marker="o", label="dataset 2")
        m = ax_a.scatter(residual_locations_transformed[1], residual_locations_transformed[0], marker="x", c=np.ndarray.flatten(np.asarray(residuals)), cmap='jet', vmin=res_mean-1.5*res_sigma, vmax=res_mean+1.5*res_sigma, s=5)
        ax_a.set_xlabel('Longitude', fontsize=14, fontweight='bold')
        ax_a.set_ylabel('Latitude', fontsize=14, fontweight='bold')
        ax_a.set_title("Found Pairs within %.2f m, n=%.2f" % (search_distance, len(residual_locations_transformed[0])))
        plt.legend()    
        fig0.colorbar(m, label='residuals (m)')
        fig0.show()

        # plot histogram of residuals
        fig1, ax_b = plt.subplots()
        ax_b.hist(np.asarray(residuals)*100, bins=50, histtype='bar', color='xkcd:navy') #range=(res_mean-3*res_sigma, res_mean+3*res_sigma)
        ax_b.minorticks_on()
        ax_b.tick_params(bottom=True, right=True, left=True, top=True, which='minor') 
        ax_b.tick_params(bottom=True, right=True, left=True, top=True, which='major') 
        ax_b.tick_params(labeltop=False, labelright=False, labelbottom=True, labelleft=True) 
        ax_b.set_title(f"Median Residual: {res_median*100:.1f}cm and 1\u03C3 SD: {res_sigma*100:.1f}cm", fontsize=14, fontname='Baskerville')
        ax_b.set_xlabel("Elevation Residual (cm)", fontsize=11, fontname='Baskerville', fontweight='light')
        ax_b.set_ylabel("Counts", fontsize=11, fontname='Baskerville', fontweight='light')
        fig1.show()

        # plot where the residuals are far from the mean residual
        fig2, ax_c = plt.subplots()
        outliers = []
        outliers_locations = []
        for i in range(0, len(residuals)):
            if residuals[i] > res_mean + 2*res_sigma or residuals[i] < res_mean - .5*res_sigma:
                outliers.append(residuals[i])
                outliers_locations.append(residual_locations[i])
        outliers_lcations_transformed = np.asarray(outliers_locations).T
        colors = np.ndarray.flatten(np.asarray(outliers))
        ax_c.scatter(lon1, lat1, c='y', s=1)
        n = ax_c.scatter(outliers_lcations_transformed[1], outliers_lcations_transformed[0], c=colors, cmap='jet', s=5)
        fig2.colorbar(n, label='residuals (m)')
        ax_c.set_title("Outlier Residual Locations")
        ax_c.set_xlabel('Longitude', fontsize=14, fontweight='bold')
        ax_c.set_ylabel('Latitude', fontsize=14, fontweight='bold')
        fig2.show()

        # plot elevation data through time both datasets
        fig3, ax_c = plt.subplots()
        ax_c.scatter(decimal_hour1, h1, s=1, label="dataset 1")
        ax_c.scatter(decimal_hour2, h2, s=1, label="dataset 2")
        plt.legend()
        fig3.show()

        plt.show()

    return number_residuals, res_median, res_sigma

def nearest_neighbor_compare_radius(lat1, lon1, h1, lat2, lon2, h2, radius):
    """
    - INPUTS: lat1/lon1/h1: single coordinate; lat2/lon2/h2: list of coordinates; search_distance = 
        how close points need to be as "neighbors"
    - Find points relative to lat/lon 1 within radius r 
    - RETURN: (1) residuals (2) standard deviation of residuals 
    *** NOTE: lat1/lon1/h1 is considered "truth"
    """

    lon_radius = []
    lat_radius = []
    height_radius = []

    for i in range(0, len(lat1)):
        d = distance.distance((lat1, lon1), (lat2[i], lon2[i])).meters
        if d < radius:
            lat_radius.append(lat2[i])
            lon_radius.append(lon2[i])
            height_radius.append(h2[i])

    if len(height_radius) == 0:
        print("no points in vicinity")
    else:
        print("number of points found: ", len(height_radius))
        plt.scatter(lon_radius, lat_radius, c=height_radius)
        plt.scatter(lon1, lat1, c='r')
        plt.rcParams.update({'font.size': 14})
        plt.xlabel('Longitude ($^o$)', fontsize=14, fontweight='bold')
        plt.ylabel('Latitude ($^o$)', fontsize=14, fontweight='bold') 
        plt.colorbar()
        plt.show()


# bias = distance from antenna base to compacted snow
bias2 =  0 # 1.797 + (.0175+.015)/2
bias1 = 0 # .245 + (.0175+.015)/2
############# PSEUDOSTATIC COMPARE
# num_PSP, mean = get_PSP_stats(lattitudes1, longitudes1, ellipsoidal_heights1, False)
# num_PSP2, mean2 = get_PSP_stats(lattitudes2, longitudes2, ellipsoidal_heights2, True)

# print("Data set 1 # PSPs: ", num_PSP)
# print("Data set 1 mean of 1sigma: ", mean)
# print("Data set 2 # PSPs: ", num_PSP2)
# print("Data set 2 mean of 1sigma: ", mean2)

number_residuals, res_median, res_sigma = nearest_neighbor_compare(lattitudes1, longitudes1, ellipsoidal_heights1, lattitudes2, longitudes2, ellipsoidal_heights2, 1, bias1, bias2, True)
print("res_sigma    ", res_sigma)
print("res_median   ", res_median)
print("num residuals", number_residuals)