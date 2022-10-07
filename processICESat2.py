# DEREK PICKELL - Adapted FROM: https://github.com/glaciology/ICEPICS_examples/blob/master/Renland_data_access.ipynb
# BEAM OUTPUT FILES: does this account for 180 degree yaw change? 
# Data filtering: read_atl06() outputs 6 files filtered for quality flag and bounding box

import numpy as np 
import matplotlib.pyplot as plt
import h5py
import os
from geopy import distance
from glob import glob
from astropy.time import Time

def gps2dyr(time):
    """ Converte GPS time to decimal years. """
    return Time(time, format='gps').decimalyear

def track_type(time, lat, tmax=1):
    """
    Separate tracks into ascending and descending.
    
    Defines tracks as segments with time breaks > tmax,
    and tests whether lat increases or decreases w/time.
    """
    tracks = np.zeros(lat.shape)  # generate track segment
    tracks[0:np.argmax(np.abs(lat))] = 1  # set values for segment
    i_asc = np.zeros(tracks.shape, dtype=bool)  # output index array

    # Loop through individual secments
    for track in np.unique(tracks):
    
        i_track, = np.where(track == tracks)  # get all pts from seg
    
        if len(i_track) < 2: continue
    
        # Test if lat increases (asc) or decreases (des) w/time
        i_min = time[i_track].argmin()
        i_max = time[i_track].argmax()
        lat_diff = lat[i_track][i_max] - lat[i_track][i_min]
    
        # Determine track type
        if lat_diff > 0:  i_asc[i_track] = True
    
    return i_asc, np.invert(i_asc)  # index vectors

def read_atl06(fname, bbox):
    """ 
    Read 1 ATL06 file and output 6 reduced files. 
    
    Extract variables of interest and separate the ATL06 file 
    into each beam (ground track) and ascending/descending orbits.
    """

    # Each beam is a group
    group = ['/gt1l', '/gt1r', '/gt2l', '/gt2r', '/gt3l', '/gt3r']

    # Loop trough beams
    for k,g in enumerate(group):
    
        #-----------------------------------#
        # 1) Read in data for a single beam #
        #-----------------------------------#
    
        # Load variables into memory (more can be added!)
        with h5py.File(fname, 'r') as fi:
            lat = fi[g+'/land_ice_segments/latitude'][:]
            lon = fi[g+'/land_ice_segments/longitude'][:]
            h_li = fi[g+'/land_ice_segments/h_li'][:]
            s_li = fi[g+'/land_ice_segments/h_li_sigma'][:]
            t_dt = fi[g+'/land_ice_segments/delta_time'][:]
            q_flag = fi[g+'/land_ice_segments/atl06_quality_summary'][:]
            t_ref = fi['/ancillary_data/atlas_sdp_gps_epoch'][:]
            rgt = fi['/orbit_info/rgt'][:] * np.ones(len(lat))
            orb = np.full_like(h_li, k)

        #---------------------------------------------#
        # 2) Filter data according region and quality #
        #---------------------------------------------#
        
        # Select a region of interest
        if bbox:
            lonmin, lonmax, latmin, latmax = bbox
            bbox_mask = (lon >= lonmin) & (lon <= lonmax) & \
                        (lat >= latmin) & (lat <= latmax)
        else:
            bbox_mask = np.ones_like(lat, dtype=bool)  # get all
            
        # Only keep good data, and data inside bbox
        mask = (q_flag == 0) & (np.abs(h_li) < 10e3) & (bbox_mask == 1)
        
        # Update variables        
        lat, lon, h_li, s_li, t_dt, q_flag, \
             rgt, orb = \
                lat[mask], lon[mask], h_li[mask], s_li[mask], t_dt[mask], \
                q_flag[mask],  \
                rgt[mask], orb[mask]
        

        # Test for no data
        if len(h_li) == 0: continue

        #-------------------------------------#
        # 3) Convert time and separate tracks #
        #-------------------------------------#
        
        # Time in GPS seconds (secs sinde 1980...)
        t_gps = t_ref + t_dt

        # Time in decimal years
        t_year = gps2dyr(t_gps)

        # Determine orbit type
        i_asc, i_des = track_type(t_year, lat)
        
        #-----------------------#
        # 4) Save selected data #
        #-----------------------#
        
        # Define output file name
        ofile = fname.replace('.h5', '_'+g[1:]+'.h5')
                
        # Save variables
        with h5py.File(ofile, 'w') as f:
            f['orbit'] = orb
            f['lon'] = lon
            f['lat'] = lat
            f['h_elv'] = h_li
            f['t_year'] = t_year
            f['t_sec'] = t_gps
            f['s_elv'] = s_li
            f['q_flg'] = q_flag
            f['rgt'] = rgt
            f['trk_type'] = i_asc

            print('out ->', ofile)

def read_h5(fname, vnames=[]):
    """ Simple HDF5 reader. """
    with h5py.File(fname, 'r') as f:
        return [f[v][:] for v in vnames]

def get_nearby_points(stationLon, stationLat, icesat_filepath, radius):
    """
    returns a cropped list of ICESat elevations that are within defined radius of 
    the GNSS station
    """
    lon, lat, t, h = read_h5(file_path, ['lon', 'lat', 't_year', 'h_elv'])
    lon_radius = []
    lat_radius = []
    height_radius = []

    for i in range(0, lon.size):
        d = distance.distance((lat[i], lon[i]), (stationLat, stationLon)).meters
        if d < radius:
            print(d)
            lat_radius.append(lat[i])
            lon_radius.append(lon[i])
            height_radius.append(h[i])

    if len(height_radius) == 0:
        print("no points in vicinity")
    else:
        print("number of points found: ", len(height_radius))
        plt.scatter(lon_radius, lat_radius, c=height_radius)
        plt.scatter(stationLon, stationLat, c='r')
        plt.rcParams.update({'font.size': 14})
        plt.xlabel('Longitude ($^o$)', fontsize=14, fontweight='bold')
        plt.ylabel('Latitude ($^o$)', fontsize=14, fontweight='bold')
        plt.title(file_path)    
        plt.colorbar()
        plt.show()

        # plt.scatter(np.ones(len(height_radius)), height_radius, c='b')
        # plt.scatter(0, 3240, c='r')
        # plt.show()


    return height_radius

if __name__ == '__main__':

    station_coordinates = [[-38.5442387, 72.61694651], [-38.47053573, 72.57373086], [-38.01298382, 72.64268509], [-37.75339548, 72.64253483], [-38.29891839, 72.58576698], [-38.03950328, 72.58514], [-37.77985812, 72.5850767], [-38.53220871, 72.64324861], [-38.55849657, 72.58575389], [-38.81850178, 72.58508771], [-38.2726624, 72.64291666], [-38.79179363, 72.64336517]]
    coordinates_transformed = np.array(station_coordinates).T

    script_dir = os.path.dirname(__file__)
    tempdir = 'tempData'
    filename = 'ATL06_20220519155012_08791505_005_02.h5'
    filename2 = 'ATL06_20220519155012_08791505_005_02_gt1l.h5'
    file_path = os.path.join(script_dir, tempdir, filename2)

    # get_nearby_points(station_coordinates[1][0], station_coordinates[1][1], file_path, 200)

    # Get individual Files
    # bbox = -38.9,-37.7,72.5,72.7 # OGRENET BOUNDING BOX
    # read_atl06(os.path.join(script_dir, tempdir, filename), bbox)

    # lon, lat, t, h = read_h5(file_path, ['lon', 'lat', 't_year', 'h_elv'])
    # plt.scatter(lon, lat, c=h)
    # plt.scatter(coordinates_transformed[0], coordinates_transformed[1], c='r')
    # plt.rcParams.update({'font.size': 14})
    # # plt.Circle((station_coordinates[1][0], station_coordinates[1][1]), 2, color='b', fill=False)
    # plt.xlabel('Longitude ($^o$)', fontsize=14, fontweight='bold')
    # plt.ylabel('Latitude ($^o$)', fontsize=14, fontweight='bold')
    # plt.colorbar()
    # plt.show()

    # print("saving CSV")
    # np.savetxt(filename2 + '.txt',np.c_[lon, lat, h])


