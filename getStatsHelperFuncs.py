from geopy import distance
# from pyproj import transform
import numpy as np
import scipy
import matplotlib.pyplot as plt



def calculateHorizontals(list_positions):
    """
    takes list of x, y coordinates, finds centroid and SD
    assumes distances are small enough to make euclidean approximation:
    https://gis.stackexchange.com/questions/58653/what-is-approximate-error-of-pythagorean-theorem-vs-haversine-formula-in-measur
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
    plt.xlabel('Longitude', fontsize=12, fontweight='light', fontname='Baskerville')
    plt.ylabel('Latitude', fontsize=12, fontweight='light', fontname='Baskerville')
    plt.title(f"Location of PSPs (n={number_PSP})", fontsize=14, fontweight='bold', fontname='Baskerville')    
    plt.show()

def plot_residuals(residuals, residual_locations, search_distance, res_mean, new_res_mean, res_median, res_sigma, lon1, lon2, lat1, lat2, hi1, hi2 ):
    # plot where the residuals are, on top of original GPS points
    fig0, ax_a = plt.subplots()
    residual_locations_transformed = np.asarray(residual_locations, dtype='object').T
    ax_a.scatter(lon1, lat1, c='k', s=25, marker=".", label="dataset 1")
    ax_a.scatter(lon2, lat2, c='dimgray', s=25, marker=".", label="dataset 2")
    m = ax_a.scatter(residual_locations_transformed[1], residual_locations_transformed[0], marker=".", c=np.ndarray.flatten(np.asarray(residuals)), cmap='jet', vmin=res_mean-3*res_sigma, vmax=res_mean+3*res_sigma, s=35)
    ax_a.set_xlabel('Longitude', fontsize=14, fontweight='bold')
    ax_a.set_ylabel('Latitude', fontsize=14, fontweight='bold')
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
        \n skew: {scipy.stats.skew(new_res_mean)}", fontsize=14, fontname='Baskerville')
    ax_b.set_xlabel("Elevation Residual (cm)", fontsize=11, fontname='Baskerville', fontweight='light')
    ax_b.set_ylabel("Counts", fontsize=11, fontname='Baskerville', fontweight='light')
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
    ax_c.set_xlabel('Longitude', fontsize=14, fontweight='bold')
    ax_c.set_ylabel('Latitude', fontsize=14, fontweight='bold')
    fig2.show()

    # plot elevation data through time both datasets
    fig3, ax_c = plt.subplots()
    ax_c.scatter(np.arange(0, len(hi1), 1), hi1, s=1, label="dataset 1")
    ax_c.scatter(np.arange(0, len(hi2), 1), hi2, s=1, label="dataset 2")
    plt.legend()
    fig3.show()

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

def plot_IR():
    ####### GNSS IR 
    IR_heights = [1.59,1.645,1.615,1.705,1.58,1.631,1.63,1.68,1.766,1.63,1.645,1.68,1.686,1.655,1.68,1.66,1.765,1.68,1.641,1.75,1.59,1.665,1.715,1.63,1.58,1.615,1.66,1.621,1.68,1.605,1.635,1.585,1.601,1.671,1.585,1.71,1.621,1.65,1.73,1.68,1.75,1.665,1.631,1.665,1.615,1.645,1.746,1.6,1.606,1.62,1.691,1.675,1.721,1.585,1.625,1.635,1.685,1.595,1.65,1.615,1.6,1.615,1.58,1.615,1.6,1.67,1.745,1.581,1.696,1.586,1.671,1.651,1.64,1.66,1.705,1.59,1.625,1.665,1.72,1.63,1.646,1.56,1.685,1.625,1.645,1.59,1.666,1.63,1.7,1.751,1.636,1.665,1.736,1.606,1.685,1.71,1.675,1.75,1.665,1.675,1.745,1.625,1.7,1.741,1.59,1.615,1.64,1.645,1.655,1.69,1.6,1.675,1.58,1.59,1.66,1.595,1.681,1.74,1.64,1.635,1.665,1.675,1.67,1.73,1.635,1.606,1.665,1.615,1.605,1.655,1.586,1.605,1.655,1.651,1.745,1.615,1.66,1.611,1.681,1.715,1.65,1.605,1.67,1.775,1.63,1.685,1.591,1.681]
    mean_IR_height = np.mean(IR_heights)
    fig2, ax2 = plt.subplots()
    ax2.hist(np.asarray(IR_heights), bins=20, histtype='bar', color='xkcd:navy', orientation='horizontal') 
    ax2.set_title("Histogram of Antenna Height Estimates", fontsize=14, fontname='Baskerville', fontweight='bold')
    ax2.set_ylabel("Estimated Height (m)", fontsize=12, fontname='Baskerville')
    ax2.set_xlabel("Counts", fontsize=12, fontname='Baskerville')
    ax2.axhline(mean_IR_height, color='red', lw=3, ls='-', label=f'mean height {mean_IR_height:.1f}m')
    ax2.axhline(mean_IR_height + np.std(IR_heights), color='red', lw=1, ls='--', label=f'1 \u03C3 {np.std(IR_heights):.2f}m')
    ax2.axhline(mean_IR_height - np.std(IR_heights), color='red', lw=1, ls='--')

    plt.legend()
    fig2.show()

    plt.show()