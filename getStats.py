# Created by Derek Pickell
# 11/11/21
# Use for: 
# Comparison of two different GNSS instruments, Kinematic PPP (CSRS PROCESSED)
# Comparison of GNSS Kinematic PPP (CSRS) & ICESAT-2 (.txt, icepyx)
# Comparison of GNSS Static PPP (CSRS) & ICESAT-2    (.txt, icepyx)
# Comparison of PPK GNSS (RTKLIB .pos) and ICESAT-2  (.txt, icepyx)

import os
import argparse
from getStatsHelperFuncs import *

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
#################################################################
# bias = distance from antenna base to compacted snow [Polypod 1.797+snow depth; Sled: 0.245+snow depth]
bias1 = .245
bias2 =  1.797 - 0.06 #.245 - .0825 #1.797 - 0.0825 - .046
search_distance = 1 # meters
#################################################################
############# PSEUDOSTATIC COMPARE #############
num_PSP, mean_1s_z, mean_1s_xy, vertical_stats1, horizontal_stats1 = get_PSP_stats(lattitudes1, longitudes1, ellipsoidal_heights1, bias1, False)
num_PSP2, mean2_1s_z,  mean2_1s_xy, vertical_stats2, horizontal_stats2 = get_PSP_stats(lattitudes2, longitudes2, ellipsoidal_heights2, bias2, False)
# print("")
# print("Dataset 1 # PSPs:           ", num_PSP)
# print("Dataset 1 mean of 1sigma z: ", mean_1s_z)
# print("Dataset 1 mean of 1s    xy: ", mean_1s_xy)
# print("Dataset 2 # PSPs:           ", num_PSP2)
# print("Dataset 2 mean of 1sigma z: ", mean2_1s_z)
# print("Dataset 2 mean of 1s    xy: ", mean2_1s_xy)
# print("")
### pickle data
# pickleme("vertical_1", vertical_stats1, os.path.basename(file_path1))
# pickleme("vertical_2", vertical_stats2, os.path.basename(file_path2))
# pickleme("horizontal_1", horizontal_stats1, os.path.basename(file_path1))
# pickleme("horizontal_2", horizontal_stats2, os.path.basename(file_path2))
### plot pickle data
# plot_PSPs_all("vertical_1")
# plot_PSPs_all("vertical_2")
# plot_PSPs_all("horizontal_1")
# plot_PSPs_all("horizontal_1")

############# RESIDUALS COMPARE #############
# get_residuals, get_residuals_no_PSPs
# residuals, number_residuals, res_median, res_sigma = get_residuals_no_PSPs(lattitudes1, longitudes1, ellipsoidal_heights1, lattitudes2, longitudes2, ellipsoidal_heights2, search_distance, bias1, bias2, False)
# print("res_sigma                   ", res_sigma)
# print("res_median                  ", res_median)
# print("n=:                         ", number_residuals)
### pickle data
# pickleme("residuals", residuals, os.path.basename(file_path1))
### plot pickle data
plot_residuals_all("residuals")
####################################################