# Created by Derek Pickell
# 11/11/21
# Use for: 
# Comparison of two different GNSS instruments, Kinematic PPP (CSRS PROCESSED)
# Comparison of GNSS Kinematic PPP (CSRS) & ICESAT-2 (.txt, icepyx)
# Comparison of GNSS Static PPP (CSRS) & ICESAT-2    (.txt, icepyx)
# Comparison of PPK GNSS (RTKLIB .pos) and ICESAT-2  (.txt, icepyx)

import numpy as np 
import os
import argparse
import matplotlib.pyplot as plt
from getStatsHelperFuncs import *
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
    """
    Parses Data from text file, in format [lon, lat, h, h_err, g_err] spaced by blank space
    """
    print("Parsing ICESat-2 ATL06")
    data = np.genfromtxt(file_path1).T
    longitudes1 = data[0]
    lattitudes1 = data[1]
    ellipsoidal_heights1 = data[2]
    h_err1 = data[3]
    g_err1 = data[4]
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
    print("Parsing ICESat-2 ATL06")
    data = np.genfromtxt(file_path2).T
    longitudes2 = data[0]
    lattitudes2 = data[1]
    ellipsoidal_heights2 = data[2]
    h_err2 = data[3]
    g_err2 = data[4]
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

def get_PSP_stats(lat, lon, h, bias, plot=True):
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
    vertical_stats = []
    psp_full_locations = []
    horizontal_stats = []
    centroids = []
    i = 0

    # correct for antenna heights
    h = h - bias

    lst = ["|","/","-","\\"]
    while i < len(lat) - 1:
        print(lst[i % 4], end="\r")

        # get distance between first two points
        d = distance.distance((lat[i], lon[i]), (lat[i+1], lon[i+1])).meters # distance.distance calculates great circle distance of ellipsoid WGS84
        if d < minimum_distance: # potential cluster found
            temp_list = []
            temp_list_positions = []
            temp_list.append(h[i]) # append first point
            temp_list_positions.append([lat[i], lon[i]])
            j = i
            while distance.distance((lat[j], lon[j]), (lat[j+1], lon[j+1])).meters < minimum_distance and j < len(lat) -2:
                temp_list.append(h[j+1])
                temp_list_positions.append([lat[j+1], lon[j+1]])
                j+=1

            if len(temp_list) > 7: # and len(temp_list) < 100: # cluster found, 10-50 seconds
                vertical_stats.append(np.std(temp_list))
                psp_full_locations.extend(temp_list_positions)
                horizontal_sd, centroid = calculateHorizontals(temp_list_positions)
                horizontal_stats.append(horizontal_sd)
                centroids.append(centroid)
            i = j # jump outside of cluster

        i +=1

    # CALCULATE STATS IN Z 
    number_PSP = len(vertical_stats)
    mean_1_sigma_z = np.mean(vertical_stats)

    # CALCULATE STATS IN X, Y
    mean_1_sigma_xy = np.mean(horizontal_stats)

    # psp_full_locations_transformed = np.asarray(psp_full_locations).T
    if plot: plot_PSPs(lon, lat, centroids, psp_full_locations, number_PSP)
        
    return number_PSP, mean_1_sigma_z, mean_1_sigma_xy

def get_residuals(lat1, lon1, h1, lat2, lon2, h2, search_distance, bias1, bias2, plot=True):
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
    hi1 = []
    hi2 = []

    array = np.asarray((lat1, lon1)).T
    array2 = np.asarray((lat2, lon2)).T

    tree = BallTree(array2, metric='euclidean')
    distances, indices = tree.query(array, k =1)

    #get elements of second array that don't have neighbors from array
    prune_indices = []

    lst = ["|","/","-","\\"]
    for i in range(0, len(array)):
        print(lst[i % 4], end="\r")
        
        if distance.distance((lat1[i], lon1[i]), (lat2[indices[i]], lon2[indices[i]])).meters < search_distance:
            residuals.append(h1[i] - h2[indices[i]])
            residual_locations.append([lat2[indices[i]], lon2[indices[i]]])
            hi1.append(h1[i])
            hi2.append(h2[indices[i]])
            prune_indices.append(indices[i])

    if len(residuals) == 0:
        print("no residuals found")
        return

    #prune arrays for second iteration
    array2_pruned = np.delete(array2, prune_indices, 0)
    lat2_pruned = np.delete(lat2, prune_indices)
    lon2_pruned = np.delete(lon2, prune_indices)
    h2_pruned = np.delete(h2, prune_indices)

    tree2 = BallTree(array, metric='euclidean')
    distances2, indices2 = tree2.query(array2_pruned, k =1)

    lst = ["|","/","-","\\"]
    for i in range(0, len(array2_pruned)):
        print(lst[i % 4], end="\r")
        
        if distance.distance((lat1[indices2[i]], lon1[indices2[i]]), (lat2_pruned[i], lon2_pruned[i])).meters < search_distance:
            residuals.append(h1[indices2[i]] - h2_pruned[i])
            residual_locations.append([lat2_pruned[i], lon2_pruned[i]])
            # hi1.append(h1[i])
            # hi2.append(h2[indices[i]])

    number_residuals = len(residuals)
    res_sigma = np.std(residuals)
    res_mean = np.mean(residuals)
    res_median = np.median(residuals)
    new_res_mean = [i for i in residuals if res_mean - 3* res_sigma<i<res_mean + 3*res_sigma] 
    new_sigma = np.std(new_res_mean)
    # print("filtered sigma:             ", new_sigma)

    if plot:
        ## Verify Data ## 
        plot_residuals(residuals, residual_locations, search_distance, res_mean, new_res_mean, res_median, res_sigma, new_sigma, hi1, hi2 )
        
    return number_residuals, res_median, res_sigma

def get_residuals_radius(lat1, lon1, h1, lat2, lon2, h2, radius):
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

def get_residuals_at_PSPs(psp_array1, location_PSP1, psp_array2, location_PSP2):
    """
    computes the mean elevation at each PSP, 
    and finds the residual if a nearby PSP exists 
    from the second dataset
    """
    PSP_residuals = []
    lat_PSP_resid = []
    lon_PSP_resid = []
    array = np.asarray(location_PSP1).T
    array2 = np.asarray(location_PSP2).T
    array = np.concatenate(([np.deg2rad(array[0])], [np.deg2rad(array[1])])).T
    array2 = np.concatenate(([np.deg2rad(array2[0])], [np.deg2rad(array2[1])])).T
    tree = BallTree(array2, metric='haversine')
    distances2, indices = tree.query(array, k=1)
    search_distance = 6 # meters
    indices = np.ndarray.flatten(indices)

    lst = ["|","/","-","\\"]
    for i in range(0, len(location_PSP1)):
        print(lst[i % 4], end="\r")
        if distance.distance((location_PSP1[i][0], location_PSP1[i][1]), (location_PSP2[indices[i]][0], location_PSP2[indices[i]][1])).meters < search_distance:
            mean1 = np.mean(psp_array1[i])
            mean2 = np.mean(psp_array2[indices[i]])
            PSP_residuals.append(mean1 - mean2)
            lat_PSP_resid.append(location_PSP1[i][0])
            lon_PSP_resid.append(location_PSP1[i][1])


    print("median residual between PSPs", np.median(PSP_residuals))
    print("n=                          ", len(PSP_residuals))
    fig, ax = plt.subplots()
    fig0, ax_a = plt.subplots()
    res_mean = np.mean(PSP_residuals)
    res_sigma = np.std(PSP_residuals)
    m = ax_a.scatter(lon_PSP_resid, lat_PSP_resid, marker="x", c=np.ndarray.flatten(np.asarray(PSP_residuals)), cmap='jet', vmin=res_mean-1.5*res_sigma, vmax=res_mean+1.5*res_sigma, s=5)
    ax_a.set_xlabel('Longitude', fontsize=14, fontweight='bold')
    ax_a.set_ylabel('Latitude', fontsize=14, fontweight='bold')
    plt.legend()    
    fig0.colorbar(m, label='residuals (m)')
    fig0.show()

    new = [i for i in PSP_residuals if res_mean - 3* res_sigma<i<res_mean + 3*res_sigma] 
    ax.hist(np.asarray(new)*100, bins=40, histtype='bar', color='xkcd:navy') #range=(res_mean-3*res_sigma, res_mean+3*res_sigma)
    ax.minorticks_on()
    ax.tick_params(bottom=True, right=True, left=True, top=True, which='minor') 
    ax.tick_params(bottom=True, right=True, left=True, top=True, which='major') 
    ax.tick_params(labeltop=False, labelright=False, labelbottom=True, labelleft=True) 
    ax.set_title(f"Median Residual: {np.median(PSP_residuals)*100:.1f}cm and 1\u03C3 SD: {np.std(PSP_residuals)*100:.1f}cm", fontsize=14, fontname='Baskerville')
    ax.set_xlabel("Elevation Residual (cm)", fontsize=11, fontname='Baskerville', fontweight='light')
    ax.set_ylabel("Counts", fontsize=11, fontname='Baskerville', fontweight='light')
    fig.show()
    plt.show()

def prune_PSPs(lat, lon, h, bias):
    """
    """
    minimum_distance = 0.05 # in meters, minimum distance to be part of a stationary PSP
    h = np.asarray(h)
    lat = np.asarray(lat)
    lon = np.asarray(lon)
    i = 0

    # correct for antenna heights
    h = h - bias

    lst = ["|","/","-","\\"]
    while i < len(lat) - 1:
        print(lst[i % 4], end="\r")

        # get distance between first two points
        d = distance.distance((lat[i], lon[i]), (lat[i+1], lon[i+1])).meters # distance.distance calculates great circle distance of ellipsoid WGS84
        if d < minimum_distance: # potential cluster found
            indices = []
            indices.append(i) # append first point
            j = i
            while distance.distance((lat[j], lon[j]), (lat[j+1], lon[j+1])).meters < minimum_distance and j < len(lat) -2:
                indices.append(j+1)
                j+=1

            if len(indices) > 7: # and len(temp_list) < 100: # cluster found, 10-50 seconds
                h = np.delete(h, indices)
                lat = np.delete(lat, indices)
                lon = np.delete(lon, indices)
            i = j # jump outside of cluster

        i +=1
    return h, lat, lon

def get_residuals_no_PSPs(lat1, lon1, h1, lat2, lon2, h2, search_distance, bias1, bias2, plot=True):
    """
    Inputs: Lat/Lon/Elevation of Each dataset, plus bias above surface and search distance
    Output: Calculates all residuals within search_distance, not including PSPs
    """
    h1_pruned, lat1_pruned, lon1_pruned = prune_PSPs(lat1, lon1, h1, bias1)
    h2_pruned, lat2_pruned, lon2_pruned = prune_PSPs(lat2, lon2, h2, bias2)
    number_residuals, res_median, res_sigma = get_residuals(lat1_pruned, lon1_pruned, h1_pruned, lat2_pruned, lon2_pruned, h2_pruned, search_distance, 0, 0, plot=True)
    print("res_sigma_pruned_PSPs       ", res_sigma)
    print("res_median_pruned_PSPs      ", res_median)
    print("n: (pruned_PSPs)            ", number_residuals)
    

#################################################################
# bias = distance from antenna base to compacted snow [Polypod 1.797+snow depth; Sled: 0.245+snow depth]
bias1 = 0
bias2 =    .245 - .0825 #1.797 - 0.0825 - .046
search_distance = 15 # meters
#################################################################
############# PSEUDOSTATIC COMPARE #############
num_PSP, mean_1s_z, mean_1s_xy = get_PSP_stats(lattitudes1, longitudes1, ellipsoidal_heights1, bias1, True)
# num_PSP2, mean2_1s_z,  mean2_1s_xy = get_PSP_stats(lattitudes2, longitudes2, ellipsoidal_heights2, bias2, False)
print("")
print("Dataset 1 # PSPs:           ", num_PSP)
print("Dataset 1 mean of 1sigma z: ", mean_1s_z)
print("Dataset 1 mean of 1s    xy: ", mean_1s_xy)
# print("Dataset 2 # PSPs:           ", num_PSP2)
# print("Dataset 2 mean of 1sigma z: ", mean2_1s_z)
# print("Dataset 2 mean of 1s    xy: ", mean2_1s_xy)
print("")
############# RESIDUALS COMPARE #############
# get_residuals_at_PSPs(psp_array1, location_PSP1, psp_array2, location_PSP2)


# test = -1
# get_residuals_no_PSPs(lattitudes1[0:test], longitudes1[0:test], ellipsoidal_heights1[0:test], lattitudes2[0:test], longitudes2[0:test], ellipsoidal_heights2[0:test], search_distance, bias1, bias2, True)
# print("")
# number_residuals, res_median, res_sigma = get_residuals(lattitudes1, longitudes1, ellipsoidal_heights1, lattitudes2, longitudes2, ellipsoidal_heights2, search_distance, bias1, bias2, True)
# print("res_sigma                   ", res_sigma)
# print("res_median                  ", res_median)
# print("n=:                         ", number_residuals)

####################################################