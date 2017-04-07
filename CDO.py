# Import modules
from cdo import *
cdo = Cdo()
import numpy as np
from netCDF4 import Dataset
from os.path import expanduser

# Set home directory input files
ref_files = expanduser("/Volumes/Naamloos/climate_data/reference/10km/merged_files/")

# Set input files
ifile_prec = ref_files+'/prec.nc'
ifile_avg_temp = ref_files+'/avg_temp.nc'
ifile_max_temp = ref_files+'/max_temp.nc'

# Select input file
input_file = "/Volumes/Naamloos/climate_data/future/entire_IGB_10km/output/finished/ydaymean/prec/ensmean_prec_rcp45_fut1.nc"

# Obtain information about the netCDF file
file_handle = Dataset(input_file, mode='r')
lat = file_handle.variables['lat'][:]
lon = file_handle.variables['lon'][:]
time = file_handle.variables['time'][:]
var = file_handle.variables['Tavg'][:]
max_value_check = np.amax(var)
min_value_check = np.amin(var)

def consecutive_dry_days():
    """
    Count the number of consecutive days where the amount of precipitation is less than a certain threshold.
    Set the number of consective days below the threshold. Precipitation threshold, unit: mm;
    Minimum number of days exceeding the threshold.

    :return: Number of dry periods of more than [nr_of_days] with a threshold of [prec_threshold]
    """

    # Define output directory and ask user for input for CDO function
    ofile_cdd = '/Volumes/Naamloos/climate_data/reference/10km/output/cdd_'
    inputs = ('mm', 'days', '.nc')
    input_mm = raw_input('Enter the threshold in [mm] for calculating the periods with dry days:')

    # CDO function to calculate the number of consecutive dry days
    cdo.eca_cdd(input_mm, input=ifile_prec,
                output=ofile_cdd + input_mm + '{0}'.format(*inputs) + '{2}'.format(*inputs))

def heat_wave():
    """
    Calculate the number of periods in a time series. This is calculated by comparing a reference temperature against the maximum temperature.
    Counted are the number of periods with at least [n] consecutive days and a temperature offset [t]. Temperature offset default is 5 C, days is 6.

    :return: Number of periods with [n] consecutive days and [t] temperature offset from the mean.
    """

    # Define output directory and ask user for input for CDO function
    ofile_hwdi = '/Volumes/Naamloos/climate_data/reference/10km/output/hwdi_'
    inputs = ('C', 'days', '.nc')
    input_temp = raw_input('Enter the offset of temperature in [degrees C] for calculating the periods with a heat wave:')
    input_days = raw_input('Enter the minimum number of consecutive days for a heat wave:')

    # CDO function to calculate the number of heat wave
    cdo.eca_hwdi(input_days, input_temp, input=" ".join((ifile_max_temp, ifile_avg_temp)),
                 output=ofile_hwdi+input_temp + '{0}'.format(*inputs) + '_' + input_days + '{1}'.format(*inputs) + '{2}'.format(*inputs))
