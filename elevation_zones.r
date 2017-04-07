library(raster)
library(sp)
library(maptools)
library(rgeos)

indir = '/users/froedevrolijk/Dropbox/Tailoring_climate_data/data/'
outdir = '/users/froedevrolijk/Dropbox/Tailoring_climate_data/data/'

extent <- raster('/Volumes/Naamloos/climate_data/reference/10km/output/cdd_17mm.nc')

shapefile_chitwan <- '/users/froedevrolijk/Dropbox/HIAWARE_tailoring/data/shapefiles/Chitwan_prj.shp'

basins <- raster('/users/froedevrolijk/Dropbox/HIAWARE_tailoring/data/basins10km.tif')

# Read elevation zones tif
elevation_zones <- raster('/users/froedevrolijk/Dropbox/HIAWARE_tailoring/data/dem10km.tif')

# Read shapefiles
shp <- readShapePoly(shapefile_chitwan)

sapply(lapply(list.files(pattern="*.shp"), readShapePoly), gArea)


raster_shapefile <- rasterize(shp ,extent, field="ID_1")

plot(raster_shapefile)

# Give zones a value for the elevation of < 500 (zone 1), 500 - 1800 (zone 2) and > 1800 (zone 3)
elevation_zones[elevation_zones <= 500] = 1
elevation_zones[elevation_zones > 500 & elevation_zones < 1800] = 2
elevation_zones[elevation_zones >= 1800] = 3
elevation_zones[elevation_zones < 0] = NA

writeRaster(elevation_zones, filename = paste(outdir, "elevation_zones.tif", sep=""), overwrite=T)

hwdi <- brick("/Volumes/Naamloos/climate_data/future/entire_IGB_10km/output/fut_hwdi_rcp45_5C_6days.nc",
              level=2, varname="heat_wave_duration_index_wrt_mean_of_reference_period")

hwdi
plot(hwdi)
summary(hwdi)
