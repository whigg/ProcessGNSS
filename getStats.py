# Created by Derek Pickell
# 11/11/21
# Use for: 
# Comparison of two different GNSS instruments, Kinematic PPP
# Comparison of two different GNSS instruments, Static PPP
# Comparison of GNSS & ICESAT-2 Kinematic PPP
# Comparison of GNSS & ICESAT-2 Static PPP
# Comparison of PPK GNSS and ICESAT-2

import numpy as np 
import os
import matplotlib.pyplot as plt
import argparse
from geopy import distance

def dir_path(string):
    if os.path.isfile(string):
        return string
    else:
        raise NotADirectoryError(string)



parser = argparse.ArgumentParser()
parser.add_argument("file_1_path", type=dir_path, help="full file path to first data file")
parser.add_argument("type", choices = ['Kin_PPP', 'Stat_PPP', 'ICESat-2', 'PPK'], help="\'Kin_PPP\' , \'Stat_PPP\', \'ICESat-2\', \'PPK\'")
parser.add_argument("file_2_path", type=dir_path, help="full file path to first data file")
parser.add_argument("type2", choices = ['Kin_PPP', 'Stat_PPP', 'ICESat-2', 'PPK'], help="\'Kin_PPP\' , \'Stat_PPP\', \'ICESat-2\', \'PPK\'")

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

def getPseudoPrecision(lat, lon, h):
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
            i = j # jump outside of cluster

        i +=1

    number_PSP = len(pseudostatic_points)
    
    stdvs = []
    for i in range(0, number_PSP):
        stdvs.append(np.std(pseudostatic_points[i]))

    mean_1sigma = np.mean(stdvs)

    return number_PSP, mean_1sigma

############# PSEUDOSTATIC COMPARE
num_PSP, mean = getPseudoPrecision(lattitudes1, longitudes1, ellipsoidal_heights1)
num_PSP2, mean2 = getPseudoPrecision(lattitudes2, longitudes2, ellipsoidal_heights2)

print("Data set 1 # PSPs: ", num_PSP)
print("Data set 1 mean of 1sigma: ", mean)
print("Data set 2 # PSPs: ", num_PSP2)
print("Data set 2 mean of 1sigma: ", mean2)