# Created by Derek Pickell
# 11/11/21
# Use for: 
# Comparison of two different GNSS instruments, Kinematic PPP (CSRS PROCESSED)
# Comparison of two different GNSS instruments, Static PPP    (CSRS PROCESSED)
# Comparison of GNSS Kinematic PPP (CSRS) & ICESAT-2 (.txt, icepyx)
# Comparison of GNSS Static PPP (CSRS) & ICESAT-2    (.txt, icepyx)
# Comparison of PPK GNSS (RTKLIB .pos) and ICESAT-2  (.txt, icepyx)

import numpy as np 
import os
import matplotlib.pyplot as plt
import argparse
from geopy import distance
from sklearn.neighbors import BallTree

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

print(args.type)
print(args.type2)
print(args.file_1_path)

file_path1 = os.path.join(args.file_1_path)
file_path2 = os.path.join(args.file_2_path)

if args.type == 'Kin_PPP':
    print("Parsing Kinematic PPP from CSRS")
    data = np.genfromtxt(file_path1, delimiter=',', skip_header=1, usecols=(0, 1, 2)).T
    lattitudes1 = data[0]
    longitudes1 = data[1]
    ellipsoidal_heights1 = data[2]

elif args.type == 'Stat_PPP':
    print("Static PPP")
elif args.type == 'ICESat-2':
    print("ICESat-2")
elif args.type =='PPK':
    print("PPK")

if args.type2 == 'Kin_PPP':
    print("Parsing Kinematic PPP from CSRS")
    data2 = np.genfromtxt(file_path2, delimiter=',', skip_header=1, usecols=(0, 1, 2)).T
    lattitudes2 = data2[0]
    longitudes2 = data2[1]
    ellipsoidal_heights2 = data2[2]
    
elif args.type2 == 'Stat_PPP':
    print("Static PPP")
elif args.type2 == 'ICESat-2':
    print("ICESat-2")
elif args.type2 =='PPK':
    print("PPK")

def getPseudoPrecision(lat, lon, h, plot=True):
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
        d = distance.distance((lat[i], lon[i]), (lat[i+1], lon[i+1])).meters
        if d < minimum_distance: # potential cluster found
            temp_list = []
            temp_list.append(h[i]) # append first point
            j = i
            while distance.distance((lat[j], lon[j]), (lat[j+1], lon[j+1])).meters < minimum_distance and j < len(lat) -2:
                temp_list.append(h[j+1])
                j+=1

            if len(temp_list) > 10 and len(temp_list) < 50:
                # print("cluser found")
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
        print(location_PSP_transformed)
        plt.scatter(lon, lat, c='b')
        plt.scatter(location_PSP_transformed[1], location_PSP_transformed[0], c='r', s=20)
        plt.rcParams.update({'font.size': 14})
        plt.xlabel('Longitude', fontsize=14, fontweight='bold')
        plt.ylabel('Latitude', fontsize=14, fontweight='bold')
        plt.title("Location of PSPs")    
        plt.show()

    return number_PSP, mean_1sigma

def getInterPrecision(lat1, lon1, h1, lat2, lon2, h2, search_distance, plot=True):
    """
    Find nearest neighbor between dataset 1 and dataset 2 points
    Get residual between these two points
    Calculate standard deviation
    Return (1) residuals (2) standard deviation of residuals 
    *** NOTE: lat1/lon1/h1 is considered "truth"
    """
    minimum_distance = search_distance # in meters, minimum distance to be part of a stationary PSP
    residuals = []
    residual_locations = []
    i = 0
    print("lat1" , len(lat1))
    print("lat2", len(lat2))
    array = np.asarray((lat1, lon1)).T
    array2 = np.asarray((lat2, lon2)).T
    tree = BallTree(array2, metric='euclidean')
    distances, indices = tree.query(array, k =1)

    lst = ["|","/","-","\\"]
    distances = []
    for i in range(0, len(lat1)):
        print(lst[i % 4], end="\r")
        
        if distance.distance((lat1[i], lon1[i]), (lat2[indices[i]], lon2[indices[i]])).meters < search_distance:
            residuals.append(h1[i] - h2[indices[i]])
            residual_locations.append([lat1[i], lon1[i]])

    number_residuals = len(residuals)
    
    res_sigma = np.std(residuals)

    print(np.mean(residuals))

    if plot:
        ## Verify Data ## 
        residual_locations_transformed = np.asarray(residual_locations).T
        plt.scatter(lon1, lat1, c='r', s=20)
        plt.scatter(lon2, lat2, c='b', s=20)
        plt.scatter(residual_locations_transformed[1], residual_locations_transformed[0], c='g', s=5)
        plt.rcParams.update({'font.size': 14})
        plt.xlabel('Longitude', fontsize=14, fontweight='bold')
        plt.ylabel('Latitude', fontsize=14, fontweight='bold')
        plt.title("Found Pairs within %.3f m" % search_distance)    
        plt.show()

    return number_residuals, res_sigma



############# PSEUDOSTATIC COMPARE
# num_PSP, mean = getPseudoPrecision(lattitudes1, longitudes1, ellipsoidal_heights1)
# num_PSP2, mean2 = getPseudoPrecision(lattitudes2, longitudes2, ellipsoidal_heights2)

# print("Data set 1 # PSPs: ", num_PSP)
# print("Data set 1 mean of 1sigma: ", mean)
# print("Data set 2 # PSPs: ", num_PSP2)
# print("Data set 2 mean of 1sigma: ", mean2)

num_residuals, res_sigma = getInterPrecision(lattitudes1, longitudes1, ellipsoidal_heights1, lattitudes2, longitudes2, ellipsoidal_heights2, 1)
print("num residuals", num_residuals)
print("res_sigma", res_sigma)