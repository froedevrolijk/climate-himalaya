import os
import glob
import subprocess

extents_dict = {'chitwan': [166058.87, 3116189.07, 312721.71, 2998326.01],
                'faisalbad': [-915605.83, 3641383.32, -732202.73, 3441863.91],
                'gandaki': [-3146.22, 3347916.98, 480904.11, 2927078.34],
                'gilgit': [-851712.00, 4209257.48, -564539.79, 3986650.86],
                'pothwar': [-1035816.99, 3958256.70, -642195.94, 3578737.26]}

def clip_rasters(extents_dict):
    """
    Function to clip netCDF files to a new extent as used for the various districts.
    :param extents_dict: The coordinates of the districts in a dictionary, formatted as ulx uly lrx lry
    :return: Clipped files for each key in the dictionary
    """

    # Define input and output folder
    dfolder = "/Volumes/Naamloos/climate_data/final_files/"
    ofolder = "/Volumes/Naamloos/climate_data/final_files/ClipRaster/"

    # Remove xml files from input folder
    xml_files = glob.glob(dfolder + '*.xml')
    [os.remove(xml) for xml in xml_files]

    # Get the coordinates for each district
    for extent_key in extents_dict:
        coordinates = extents_dict.get(extent_key)

        # Clip each file in the dictionary by calling gdal translate via subprocess
        for key in dictionary:
            indexes = dictionary[key]
            ncfile = indexes[0]
            vars = indexes[1]

            # Call GDAL translate via subprocess to execute the command in shell
            subprocess.call('gdal_translate -of netCDF -projwin' + ' ' + str(coordinates[0]) + ' ' + str(coordinates[1]) + ' ' + str(
                            coordinates[2]) + ' ' + str(coordinates[3]) + ' NETCDF:"' + ncfile + ':' + vars + '" ' + ofolder + key +
                            '_' + extent_key + '_clip.tif', shell=True, cwd=dfolder)
        print 'All done!'
