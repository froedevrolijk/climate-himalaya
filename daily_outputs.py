import os
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import numpy as np

def daily_outputs(dictionary, key, cs, nr_of_days):
    """
    Create map with daily output for precipitation and temperature
    :param dictionary: Dictionary containing information about file name, variable name, map title and scalebar text
    :param key: The key in the dictionary pointing to the file that should be used
    :param cs: The colour scheme to be used for the output maps, e.g. 'Blues' or 'Oranges'
    :param nr_of_days: Number of maps (one map for each day) that the function should return
    :return: One map for each day, showing daily precipitation or temperature values
    """

    # Read data from dictionary
    datas = read_data(dictionary, key)[0]
    scalebar = dictionary[key][2]

    # Define the date range that should be used for the output maps
    days = pd.date_range(start='1980-12-31', end='2010-12-31', freq='D', closed='right')
    days = pd.Series(days.format())

    # User defined data type
    type = raw_input('Which type of data should be used? (e.g. mean precipitation, maximum temperature)')

    for i in range(0, nr_of_days):
        # Define the base map projection, extent and appearance
        m = Basemap(projection='tmerc', resolution='l', height=2300000, width=3500000, lat_0=29, lon_0=82)
        m.drawmapboundary()
        m.drawparallels(np.arange(-80., 81., 10.), labels=[1, 0, 0, 0])
        m.drawmeridians(np.arange(-180., 181., 10.), labels=[0, 0, 0, 1])
        m.drawcountries()
        m.shadedrelief()

        # Obtain the latitudes and longitudes, point to data and read the minimum and maximum values for the scalebar
        xi, yi = m(tr_lon, tr_lat)
        netcdf = m.pcolormesh(xi, yi, datas[i], cmap=cs)
        m.colorbar(netcdf, location='bottom', pad='10%', label=scalebar)
        sb_min = read_data(dictionary, key)[1]
        sb_max = read_data(dictionary, key)[2]
        plt.clim(sb_min, sb_max)

        # Define location for map output and create filename and title for the map
        path = '/users/froedevrolijk/Geo/outputs_hiaware/entire_IGB/daily'
        filename = dictionary[key][1]+'{0}.png'.format(days[0:nr_of_days][i])
        plt.title(type + ' for {0}'.format(days[0:nr_of_days][i]))
        filename = os.path.join(path, filename)

        # If the output file already exists, remove the old file
        if os.path.isfile(filename):
            os.remove(filename)

        # Save the map, and close the plotting window
        plt.savefig(filename)
        plt.close()
