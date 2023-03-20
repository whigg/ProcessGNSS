from geopy import distance
# from pyproj import transform
import numpy as np
import scipy
import matplotlib.pyplot as plt
import matplotlib as mpl
import pickle
from glob import glob
from sklearn.neighbors import BallTree

mpl.rcParams['axes.labelsize'] = 12
mpl.rcParams['axes.labelweight'] = 'light'
mpl.rcParams['axes.titlesize'] = 14
mpl.rcParams['axes.titleweight'] = 'bold'
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'Optima'
plt.rcParams["figure.figsize"] = (3.34,3.34)
plt.rcParams['figure.constrained_layout.use'] = 'True'


def calculateHorizontals(list_positions):
    """
    INPUTS: list of x, y coordinates, finds centroid and SD
    assumes distances are small enough to make euclidean approximation:
    https://gis.stackexchange.com/questions/58653/what-is-approximate-error-of-pythagorean-theorem-vs-haversine-formula-in-measur
    OUTPUTS: standard deviation of radial lengths between centroid and point; location of centroid
    """
    positions_T = np.asarray(list_positions).T
    x_avg = np.sum(positions_T[0])/len(positions_T[0])
    y_avg = np.sum(positions_T[1])/len(positions_T[1])

    centroid = [x_avg, y_avg]
    distances_from_centroid = []

    for i in range(0, len(positions_T[0])):
        d = distance.distance((x_avg, y_avg), (positions_T[0][i], positions_T[1][i]), ).meters # distance.distance calculates great circle distance of ellipsoid WGS84
        distances_from_centroid.append(d)

    #calculate stats
    horizontal_sd = np.std(distances_from_centroid)

    return horizontal_sd, centroid

def plot_PSPs(lon, lat, centroids, psp_full_locations, number_PSP):
    """plot PSPs on top of all data, plus centroid of PSPs and ID'd points"""
    centroids_T = np.asarray(centroids).T
    psps_T = np.asarray(psp_full_locations).T
    ## Verify Data ## 
    plt.scatter(lon, lat, c='b', s=5)
    plt.scatter(psps_T[1], psps_T[0], c='r', s=5)
    plt.scatter(centroids_T[1], centroids_T[0], c='g', s=20)
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title(f"Location of PSPs (n={number_PSP})")    
    plt.show()

def plot_residuals(residuals, residual_locations, search_distance, res_mean, new_res_mean, res_median, res_sigma, lon1, lon2, lat1, lat2, hi1, hi2 ):
    """plot residuals: location, histogram, outliers, etc."""
    # plot where the residuals are, on top of original GPS points
    fig0, ax_a = plt.subplots()
    residual_locations_transformed = np.asarray(residual_locations, dtype='object').T
    ax_a.scatter(lon1, lat1, c='k', s=25, marker=".", label="dataset 1")
    ax_a.scatter(lon2, lat2, c='dimgray', s=25, marker=".", label="dataset 2")
    m = ax_a.scatter(residual_locations_transformed[1], residual_locations_transformed[0], marker=".", c=np.ndarray.flatten(np.asarray(residuals)), cmap='jet', vmin=res_mean-3*res_sigma, vmax=res_mean+3*res_sigma, s=35)
    ax_a.set_xlabel('Longitude')
    ax_a.set_ylabel('Latitude')
    ax_a.set_title("Found Pairs within %.2f m, n=%.2f" % (search_distance, len(residual_locations_transformed[0])))
    plt.legend()    
    fig0.colorbar(m, label='residuals (m)')
    fig0.show()

    # plot histogram of residuals
    fig1, ax_b = plt.subplots()
    ax_b.hist(np.asarray(new_res_mean)*100, bins=25, histtype='bar', color='xkcd:navy') #range=(res_mean-3*res_sigma, res_mean+3*res_sigma)
    ax_b.minorticks_on()
    ax_b.tick_params(bottom=True, right=True, left=True, top=True, which='minor') 
    ax_b.tick_params(bottom=True, right=True, left=True, top=True, which='major') 
    ax_b.tick_params(labeltop=False, labelright=False, labelbottom=True, labelleft=True) 
    ax_b.set_title(f"Median Residual: {res_median*100:.1f}cm and 1\u03C3 SD: {res_sigma*100:.1f}cm \
        \n skew: {scipy.stats.skew(new_res_mean)}")
    ax_b.set_xlabel("Elevation Residual (cm)")
    ax_b.set_ylabel("Counts")
    fig1.show()

    # plot where the residuals are far from the mean residual
    fig2, ax_c = plt.subplots()
    outliers = []
    outliers_locations = []
    for i in range(0, len(residuals)):
        if  residuals[i] > res_mean + 3*res_sigma or residuals[i] < res_mean - 3*res_sigma:
            outliers.append(residuals[i])
            outliers_locations.append(residual_locations[i])
    outliers_locations_transformed = np.asarray(outliers_locations, dtype='object').T
    outliers_locations_transformed = outliers_locations_transformed.astype('float')
    if len(outliers)>0:
        print("outliers found              ", len(outliers))
        if np.shape(outliers_locations_transformed) == (1, 2, len(outliers)):
            outliers_locations_transformed = np.squeeze(outliers_locations_transformed)
        colors = np.ndarray.flatten(np.asarray(outliers))
        ax_c.scatter(lon1, lat1, c='y', s=1)
        ax_c.scatter(lon2, lat2, c='dimgray', s=1, marker="o")
        n = ax_c.scatter(outliers_locations_transformed[1], outliers_locations_transformed[0], c=colors, cmap='jet', s=5)
        fig2.colorbar(n, label='residuals (m)')
    ax_c.set_title("Outlier Residual Locations")
    ax_c.set_xlabel('Longitude')
    ax_c.set_ylabel('Latitude')
    fig2.show()

    # plot elevation data through time both datasets
    fig3, ax_c = plt.subplots()
    ax_c.scatter(np.arange(0, len(hi1), 1), hi1, s=1, label="dataset 1")
    ax_c.scatter(np.arange(0, len(hi2), 1), hi2, s=1, label="dataset 2")
    plt.legend()
    fig3.show()

    plt.show()

def plot_residuals_all(list_type):
    search = "**/" + list_type + '*.pickled'
    list_pickled = glob(search)
    print(list_pickled)
    aggregated_list = []
    for i in range(0, len(list_pickled)):
        with open(list_pickled[i], 'rb') as f:
            aggregated_list.extend(pickle.load(f))
        # plot histogram of residuals
    print("mean residual:              ", np.mean(aggregated_list))
    print("standard deviation:         ", np.std(aggregated_list))
    fig1, ax = plt.subplots()
    ax.hist(np.asarray(aggregated_list)*100, bins=30, histtype='bar', color='xkcd:navy') #range=(res_mean-3*res_sigma, res_mean+3*res_sigma)
    ax.minorticks_on()
    ax.tick_params(bottom=True, right=True, left=True, top=True, which='minor') 
    ax.tick_params(bottom=True, right=True, left=True, top=True, which='major') 
    ax.tick_params(labeltop=False, labelright=False, labelbottom=True, labelleft=True) 
    ax.set_title(f"All Residuals (n = {len(aggregated_list)}) \n skew: {scipy.stats.skew(aggregated_list)}")
    ax.set_xlabel("Residuals (cm)")
    ax.set_ylabel("Counts")
    fig1.show()
    plt.show()

def plot_PSPs_all(list_type1, list_type2):
    """Plots PSPs across all surveys, represented by pickled data that needs 
    to be pre-created for each survey"""
    list_pickled = glob(list_type1 + '*.pickled')
    aggregated_list = []
    for i in range(0, len(list_pickled)):
        with open(os.path.join("/tempData/", list_pickled[i]), 'rb') as f:
            aggregated_list.extend(pickle.load(f))

    list_pickled = glob(list_type2 + '*.pickled')
    aggregated_list2 = []       
    for i in range(0, len(list_pickled)):
        with open(os.path.join("/tempData/", list_pickled[i]), 'rb') as f:
            aggregated_list2.extend(pickle.load(f))   

    fig1, ax = plt.subplots()
    ax.hist(np.asarray(aggregated_list)*100, bins=20, histtype='bar', color='xkcd:blue', alpha=.5, label='OGRE') 
    ax.hist(np.asarray(aggregated_list2)*100, bins=20, histtype='bar', color='xkcd:yellow', alpha=.5, label='Trimble R7') 
    ax.minorticks_on()
    ax.tick_params(bottom=True, right=True, left=True, top=True, which='minor') 
    ax.tick_params(bottom=True, right=True, left=True, top=True, which='major') 
    ax.tick_params(labeltop=False, labelright=False, labelbottom=True, labelleft=True) 
    ax.set_title("All PSPs")
    ax.set_xlabel("1\u03C3 of PSPs (cm)")
    ax.set_ylabel("Counts")
    fig1.show()
    plt.legend() 
    plt.show()
    
def plot_IR():
    ####### GNSS IR 
    IR_heights = [1.59,1.645,1.615,1.705,1.58,1.631,1.63,1.68,1.766,1.63,1.645,1.68,1.686,1.655,1.68,1.66,1.765,1.68,1.641,1.75,1.59,1.665,1.715,1.63,1.58,1.615,1.66,1.621,1.68,1.605,1.635,1.585,1.601,1.671,1.585,1.71,1.621,1.65,1.73,1.68,1.75,1.665,1.631,1.665,1.615,1.645,1.746,1.6,1.606,1.62,1.691,1.675,1.721,1.585,1.625,1.635,1.685,1.595,1.65,1.615,1.6,1.615,1.58,1.615,1.6,1.67,1.745,1.581,1.696,1.586,1.671,1.651,1.64,1.66,1.705,1.59,1.625,1.665,1.72,1.63,1.646,1.56,1.685,1.625,1.645,1.59,1.666,1.63,1.7,1.751,1.636,1.665,1.736,1.606,1.685,1.71,1.675,1.75,1.665,1.675,1.745,1.625,1.7,1.741,1.59,1.615,1.64,1.645,1.655,1.69,1.6,1.675,1.58,1.59,1.66,1.595,1.681,1.74,1.64,1.635,1.665,1.675,1.67,1.73,1.635,1.606,1.665,1.615,1.605,1.655,1.586,1.605,1.655,1.651,1.745,1.615,1.66,1.611,1.681,1.715,1.65,1.605,1.67,1.775,1.63,1.685,1.591,1.681]
    mean_IR_height = np.mean(IR_heights)
    fig2, ax2 = plt.subplots()
    ax2.hist(np.asarray(IR_heights), bins=20, histtype='bar', color='xkcd:navy', orientation='horizontal') 
    ax2.set_title("Histogram of Antenna Height Estimates")
    ax2.set_ylabel("Estimated Height (m)")
    ax2.set_xlabel("Counts")
    ax2.axhline(mean_IR_height, color='red', lw=3, ls='-', label=f'mean height {mean_IR_height:.1f}m')
    ax2.axhline(mean_IR_height + np.std(IR_heights), color='red', lw=1, ls='--', label=f'1 \u03C3 {np.std(IR_heights):.2f}m')
    ax2.axhline(mean_IR_height - np.std(IR_heights), color='red', lw=1, ls='--')

    plt.legend()
    fig2.show()

    plt.show()

def plot_carrier():
    dataset1_phases = [1, 2, 3, 4, 5, 6]
    dataset2_phases = [1.1, 2.1, 3.1, 4.1, 5.1, 6.1]
    width = .35
    ind = np.arange(len(dataset1_phases))
    fig, ax = plt.subplots()
    ax.bar(ind, dataset1_phases, width, color='royalblue', label='OGRE')
    ax.bar(ind+width, dataset2_phases, width, color='seagreen', label='Trimble NetR9')
    ax.set_ylabel('Mean Carrier-phase Residuals')
    ax.set_xticks(ind + width / 2)
    ax.set_xticklabels( (f'5\u00B0', f'5\u00B0', f'5\u00B0', f'5\u00B0', f'5\u00B0', f'5\u00B0') )
    ax.set_xlabel("Elevation Angle")
    plt.legend(loc='best')
    fig.show()
    plt.show()

def calc_baseline_UTM(list1, list2):
    """
    list format: [[N, E, Alt]]
    """
    distances = []
    for i in range(0, len(list1)):
        # distances.append(math.dist(list1[i]-list2[i]))
        distances.append(np.linalg.norm(list1[i]- list2[i]))

    return distances

def calc_baseline_convert(list1, list2):
    """calculate the baseline length in ECEF (untested)"""
    distances = []
    transformer = pyproj.Transformer.from_crs(
    {"proj":'latlong', "ellps":'GRS80', "datum":'ITRF20'},
    {"proj":'geocent', "ellps":'GRS80', "datum":'ITRF20'})

    for i in range(0, len(list1)):
        x1 ,y1, z1 = transformer.transform(list1[i][0],list1[i][1],list1[i][2],radians = False)
        x2 ,y2, z2 = transformer.transform(list2[i][0],list2[i][1],list2[i][2],radians = False)
        point1 = np.array((x1 ,y1, z1))
        point2 = np.array((x2 ,y2, z2))
        distances.append(np.linalg.norm(point1 - point2))

    return distances

def pickleme(list_type, list, filename):
    filename_new = list_type + filename + ".pickled"
    path = os.path.join(os.getcwd(), "/tempData/", filename_new)
    with open(path, 'wb') as f:
        pickle.dump(list, f)

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
        
    return number_PSP, mean_1_sigma_z, mean_1_sigma_xy, vertical_stats, horizontal_stats

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
    # new_sigma = np.std(new_res_mean)
    # print("filtered sigma:             ", new_sigma)

    if plot:
        ## Verify Data ## 
        plot_residuals(residuals, residual_locations, search_distance, res_mean, new_res_mean, res_median, res_sigma, lon1, lon2, lat1, lat2, hi1, hi2 )
        
    return residuals, number_residuals, res_median, res_sigma

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

def prune_PSPs(lat, lon, h):
    """
    """
    minimum_distance = 0.05 # in meters, minimum distance to be part of a stationary PSP
    h = np.asarray(h)
    lat = np.asarray(lat)
    lon = np.asarray(lon)
    i = 0
    delete_indices = []

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
                delete_indices.extend(indices)
            i = j # jump outside of cluster

        i +=1
    h = np.delete(h, delete_indices)
    lat = np.delete(lat, delete_indices)
    lon = np.delete(lon, delete_indices)
    return h, lat, lon

def get_residuals_no_PSPs(lat1, lon1, h1, lat2, lon2, h2, search_distance, bias1, bias2, plot=True):
    """
    Inputs: Lat/Lon/Elevation of Each dataset, plus bias above surface and search distance
    Output: Calculates all residuals within search_distance, not including PSPs
    """
    h1_pruned, lat1_pruned, lon1_pruned = prune_PSPs(lat1, lon1, h1)
    h2_pruned, lat2_pruned, lon2_pruned = prune_PSPs(lat2, lon2, h2)
    return get_residuals(lat1_pruned, lon1_pruned, h1_pruned, lat2_pruned, lon2_pruned, h2_pruned, search_distance, bias1, bias2, plot)
    