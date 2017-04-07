# Import modules
import mpl_toolkits.basemap.pyproj as pyproj
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# Read netCDF
file = 'prec_mean_1981-2010_yearly.nc'
file_handle = Dataset(file, mode='r')
lon = file_handle.variables['lon'][:]  # 320
lat = file_handle.variables['lat'][:]  # 190
pr = file_handle.variables['P'][:]  # (30, 190, 320)
file_handle.close()

# Source projection
UTM_45N = pyproj.Proj("+init=EPSG:32645")

# Target projection
WGS4 = pyproj.Proj("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# Transformation
lons, lats = np.meshgrid(lon, lat)
tr_lon, tr_lat = pyproj.transform(UTM_45N, WGS84, lons, lats)

# Basemap
plt.figure()
m = Basemap(projection='tmerc', resolution='l',  height=3000000, width=3000000, lat_0=27.5, lon_0=82.5)
m.drawcoastlines(linewidth=1.5)
xi, yi = m(tr_lon, tr_lat)
cs = m.pcolormesh(xi, yi, pr[0])
plt.show()