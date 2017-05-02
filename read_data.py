import numpy as np
from netCDF4 import Dataset

def read_data(dictionary, key):
    """
    Reads data from a dictionary containing the file names and variable names of netCDF files
    :param dictionary: Dictionary containing the file names and variable names of each output file
    :param key: Key to be used to read data regarding variable information
    :return: Data of the (var), the minimum value (min_value) and maximum value (max_value) for the input file
    """

    file_name = dictionary[key][0]
    variable = dictionary[key][1]
    path = '/Volumes/Naamloos/climate_data/final_files/'
    data = path + file_name
    file_handle = Dataset(data, mode='r')
    var = file_handle.variables[variable][:]
    min_value = np.amin(var)
    max_value = np.amax(var)
    return var, min_value, max_value

read_scalebar(dictionary, 'prec_mean')

