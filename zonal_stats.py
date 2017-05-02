import os
import glob
import shutil
import processing
import pandas as pd

# Zonal statistics
def zonal_stats(extents_dict):
    """
    Calculate zonal statistics for a number of districts, using coordinates from a dictionary.
    :param extents_dict: Dictionary containing the coordinates of the districts
    :return: CSV output for each input file containing zonal statistical values for each time step
    """
    # Download folder for output
    dfolder = "/Volumes/Naamloos/climate_data/final_files/ZonalStats/"

    # Create folders per district for output, and remove old files if necessary
    try:
        if not os.path.exists(os.path.join(dfolder + str(extents_dict.keys))):
            [os.mkdir(os.path.join(dfolder + extents_key)) for extents_key in extents_dict]
    except:
        output_files = glob.glob('/Volumes/Naamloos/climate_data/final_files/ZonalStats/*')
        [shutil.rmtree(output_file, ignore_errors=True) for output_file in output_files]

    # List shapefiles to calculate statistics
    shp_path = '/Users/froedevrolijk/Dropbox/Tailoring_climate_data/data/all_proj_shp/'
    shps = glob.glob(shp_path + '*.shp')

    # List netCDF files
    local = os.path.expanduser('/Volumes/Naamloos/climate_data/final_files/')
    netcdf_path = glob.glob(local + '*.nc')

    # Find files in clipped directory
    clipped_files_path = os.path.join(local + 'ClipRaster/')
    clipped_files = glob.glob(clipped_files_path + '*.tif')

    # Remove xml files
    xml_files = glob.glob(dfolder + '*.xml')
    [os.remove(xml) for xml in xml_files]

    # Split files into district groups and create a list of the districts
    gandaki = [c for c in clipped_files if "gandaki" in c]
    faisalbad = [c for c in clipped_files if "faisalbad" in c]
    chitwan = [c for c in clipped_files if "chitwan" in c]
    pothwar = [c for c in clipped_files if "pothwar" in c]
    gilgit = [c for c in clipped_files if "gilgit" in c]

    # Read number of raster layers
    for clipped_file in clipped_files:
        nr_of_bands =

    ############################

    # in the Python Console, load the layer and ensure that it is valid:
    rasterLyr = QgsRasterLayer("/qgis_data/rasters/satimage.tif", "Sat Image")
    rasterLyr.isValid()

    # Count the number of bands for each file
    rasterLyr.bandCount()

    # for tif in gandaki:
    processing.runalg("saga:gridstatisticsforpolygons", tif, gandaki_shp, False, False, False, False,
                      False, True, False, False, 0, dfolder + 'cdd_cddi_per_time_period_chitwan_clip_b1.csv'[tif])

    ############################

    # Write csv output for each district
    for district_files in district_list:
        district = district_names[district_list.index(district_files)]
        # List csv files in output directory
        csv_files_path = os.path.join(local + 'ZonalStats/')

        # Loop through the tif files of each district to get csv files (linked to the number of years)
        for tif_file in district_files:
            csv_files = glob.glob(csv_files_path + '*' + tif_file.rstrip('.tif').split('/')[-1] + '*.csv')

            # Check if output is generated for the csv file
            if len(csv_files) > 0:
                summed_data = pd.DataFrame(columns=['Value'], index=range(1, len(csv_files) + 1))
                summed_data = summed_data.fillna(0.0)

                # Retrieve output value from csv files
                row = 0
                for csv in csv_files:
                    row += 1
                    csv = pd.read_csv(csv, header=0)
                    cleaned = csv.iloc[:, -1]

                    summed_data.loc[row, ] = round(float(cleaned), ndigits=1)

                print(summed_data)

                # Create folders for output
                try:
                    if not os.path.exists(os.path.join(csv_files_path + str(district))):
                        os.mkdir(os.path.join(csv_files_path + str(district)))
                except:
                    pass

                summed_data.to_csv(csv_files_path + str(district) + '/' + tif_file.rstrip('.tif').split('/')[-1] + '.csv')