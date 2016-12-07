'''
Created on Jul 14, 2015

auxiliary functions for SHERPA
@author: degraba
'''

from netCDF4 import Dataset
from numpy import zeros, power, sqrt


# create a dictionary with emission reductions per precursor and nuts
#--------------------------------------------------------------------

# file format:
# POLL    MS1    MS2    MS3    MS4    MS5    MS6    MS7    MS8    MS9    MS10
# NOx    0    0    0    0    0    100    0    0    0    0
# NMVOC    0    0    0    0    0    100    0    0    0    0
# NH3    0    0    0    0    0    100    0    0    0    0
# PM25    0    0    0    0    0    100    0    0    0    0
# SOx    0    0    0    0    0    100    0    0    0    0

def create_emission_reduction_dict(path_reduction_txt):
    # read emission reductions per precursor and macro sector
    f = open(path_reduction_txt, 'r')
    emission_reduction_dict = {}
    f.readline()
    while True:
        line = f.readline().rstrip()
        if len(line) == 0:
            break
        value_lst = line.split('\t')
        # sub dictionary per precursor
        precursor = value_lst[0]
        emission_reduction_dict[precursor] = {}
        for snap in range(1, 11):
            emission_reduction_dict[precursor][snap] = float(value_lst[snap]) / 100.0
    f.close()
    
    return emission_reduction_dict

# function making a list of the precursors that are reduced
def create_reduced_precursor_lst(emission_reduction_dict):
    reduced_precursor_lst = []
    for precursor in emission_reduction_dict.keys():
        sum_reductions = 0
        for snap in emission_reduction_dict[precursor].keys():
            sum_reductions += emission_reduction_dict[precursor][snap]
        if sum_reductions > 0:
            reduced_precursor_lst.append(precursor)
    return reduced_precursor_lst


# create a dictionary with emissions per precursor, macrosector and postion (lat, lon)
#-------------------------------------------------------------------------------------
def create_emission_dict(path_emission_cdf, precursor_lst):
    # open the emission netcdf
    rootgrp = Dataset(path_emission_cdf, 'r')
    
    emission_dict = {}
    for precursor in precursor_lst:
        emission_dict[precursor] = rootgrp.variables[precursor][:, :, :]
    
    # get snap, longitude and latitude arrays from emission file
    snap_array = range(1, 11)
    lon_array = rootgrp.variables['longitude'][:]
    lat_array = rootgrp.variables['latitude'][:]
    emission_dict['Nsnaps'] = snap_array
    emission_dict['lon_array'] = lon_array
    emission_dict['lat_array'] = lat_array
    
    # close the emission file
    rootgrp.close()
    
    return emission_dict


# make a window with cell distances to the central cell
# -----------------------------------------------------
def create_window(radius):
    # the window contains the distance between each cell and the centre of the window
    # the distance is expressed in cells
    n_lon_win = 2 * radius + 1
    n_lat_win = 2 * radius + 1
     
    window = zeros((n_lat_win, n_lon_win))
    i_centre = radius  
    j_centre = radius  
    for iw in range(n_lon_win):
        for jw in range(n_lon_win):
            cell_dist = sqrt((float(iw - i_centre)) ** 2 + (float(jw - j_centre)) ** 2) 
            window[iw, jw] = 1 / (1 + cell_dist) 
     
    return window

# convert to progress log file to a dictionary
def read_progress_log(progresslog):
    progress_dict = {}
    f_prog = open(progresslog, 'r')
    line = f_prog.readline().rstrip()
    [start, divisor] = line.split('\t')
    progress_dict['start'] = float(start)
    progress_dict['divisor'] = float(divisor) 
    progress_dict['netcdf_output'] = False 

    return progress_dict

# write progress log file
def write_progress_log(progress_log_filename, start, divisor):
    # write progress log file
    f_prog = open(progress_log_filename, 'w')
    f_prog.write('%f\t%f' % (start, divisor))
    f_prog.close()
    
if __name__ == '__main__':
    
    # check the window function
    radius = 200
    testwindow = create_window(radius)
    window_file = open('C:/temp/source_recptor_window.txt', 'w')
    for i in range(2 * radius + 1):
        for j in range(2 * radius + 1):
            window_file.write('%e\t' % testwindow[i,j])
        window_file.write('\n')
    window_file.close()
    pass


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
    