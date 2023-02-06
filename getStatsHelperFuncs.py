from geopy import distance
import numpy as np
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

def plot_residuals(residuals, residual_locations, search_distance, res_mean, new_res_mean, res_median, res_sigma, new_sigma, hi1, hi2 ):
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
    ax_b.hist(np.asarray(new_res_mean)*100, bins=50, histtype='bar', color='xkcd:navy') #range=(res_mean-3*res_sigma, res_mean+3*res_sigma)
    ax_b.minorticks_on()
    ax_b.tick_params(bottom=True, right=True, left=True, top=True, which='minor') 
    ax_b.tick_params(bottom=True, right=True, left=True, top=True, which='major') 
    ax_b.tick_params(labeltop=False, labelright=False, labelbottom=True, labelleft=True) 
    ax_b.set_title(f"Median Residual: {res_median*100:.1f}cm and 1\u03C3 SD: {new_sigma*100:.1f}cm \
        \n skew: {scipy.stats.skew(new_res_mean)}", fontsize=14, fontname='Baskerville')
    ax_b.set_xlabel("Elevation Residual (cm)", fontsize=11, fontname='Baskerville', fontweight='light')
    ax_b.set_ylabel("Counts", fontsize=11, fontname='Baskerville', fontweight='light')
    fig1.show()

    # plot where the residuals are far from the mean residual
    fig2, ax_c = plt.subplots()
    outliers = []
    outliers_locations = []
    for i in range(0, len(residuals)):
        if residuals[i] < residuals[i] < res_mean + 3*res_sigma or residuals[i] < res_mean - 3*res_sigma:
            outliers.append(residuals[i])
            outliers_locations.append(residual_locations[i])
    outliers_locations_transformed = np.asarray(outliers_locations, dtype='object').T
    colors = np.ndarray.flatten(np.asarray(outliers))
    ax_c.scatter(lon1, lat1, c='y', s=1)
    ax_c.scatter(lon2, lat2, c='dimgray', s=1, marker="o")
    if len(outliers)>0:
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
