'''
Created on Jun 23, 2015

In this module the concentrations are calculated for a given emission reductions scenario
Inputs: - baseline emissions (per precursor and cell),
        - a netcdf defining the area where emission reduction will be applied
        - emission reductions (per precursor and macrosector)
        - model coefficients (per pollutant, precursor, cell)
        - path where the results will be written
        - eventually a progress log to be able to report on the progress of the whole calculation when module 1
        is called by another moduel
        
output: - netcdf with concentration changes per pollutant and cell
        - delta emission netcdf with emission changes per precursor and cell
       
The calculation is optimized using a flat weight over the whole domain. This allows to update only the scale
factor of this flat weight. The bell shape central weighting factors have to be recalculated for each cell.
        
@author: degraba
Enrico agrees with this nice explanation of module 1
'''

# imports
from netCDF4 import Dataset
from numpy import lib, zeros, sum, power, ones
from math import isnan
from sherpa_globals import path_result_cdf_test
# path_emission_cdf_test, path_area_cdf_test, path_reduction_txt_test, path_model_cdf_test,
from sherpa_auxiliaries import create_emission_reduction_dict, create_emission_dict, create_window, read_progress_log, write_progress_log
import sys
from time import time
from os import remove

# function that applies reductions per snap sector and precursor to the emission netcdf
def create_delta_emission(path_emission_cdf, precursor_lst, path_area_cdf, path_reduction_txt, path_result_cdf, write_netcdf_output):
    # create a dictionary with reductions per precursor and macro sector
    emission_reduction_dict = create_emission_reduction_dict(path_reduction_txt)
    
    # open the emission netcdf
    emission_dict = create_emission_dict(path_emission_cdf, precursor_lst)
    
    # open the area netcdf
    rootgrp = Dataset(path_area_cdf, 'r')
    reduction_area = rootgrp.variables['AREA'][:] / 100.0
    rootgrp.close()
    
    # calculate a dictionary with the emission reductions per pollutant, macrosector and position
    delta_emission_dict = {}
    for precursor in precursor_lst:
        delta_emission_dict[precursor] = zeros(emission_dict[precursor].shape)
        # calculate the emission reduction
        # reductions are positive!
        # make the sum over all snap sectors
        for snap in range(1, 11):
            delta_emission_dict[precursor][snap - 1, :, :] = emission_dict[precursor][snap - 1] * reduction_area * emission_reduction_dict[precursor][snap]
#             print(snap)
#             print(sum(delta_emission_dict[precursor][snap - 1, :, :]))
    

    # before summing over all snap sectors write the delta emissions per precursor and snap to a netcdf
    # create an output netcdf with delta emissions
    # --------------------------------------------
    if write_netcdf_output == True:
        filename_delta_emission_cdf = path_result_cdf + 'delta_emission.nc'
        rootgrp = Dataset(filename_delta_emission_cdf, 'w', format='NETCDF3_CLASSIC')
     
        # create dimensions in the netcdf file
        rootgrp.createDimension('latitude', len(emission_dict['lat_array']))
        rootgrp.createDimension('longitude', len(emission_dict['lon_array']))
        rootgrp.createDimension('Nsnaps', len(emission_dict['Nsnaps']))
        latitudes = rootgrp.createVariable('latitude', 'f4', ('latitude',))
        latitudes.units = "degrees_north"
        latitudes[:] = emission_dict['lat_array']
        longitudes = rootgrp.createVariable('longitude', 'f4', ('longitude',))
        longitudes.units = "degrees_east"
        longitudes[:] = emission_dict['lon_array']
        Nsnaps = rootgrp.createVariable('Nsnaps', 'f4', ('Nsnaps',))
        Nsnaps[:] = emission_dict['Nsnaps']
        
        # create delta emission data
        for precursor in precursor_lst:
            delta_emission_precursor = rootgrp.createVariable(precursor, 'f4', ('Nsnaps', 'latitude', 'longitude',))
            delta_emission_precursor.units = "Mg/km2"
            delta_emission_precursor[:] = delta_emission_dict[precursor]
         
        rootgrp.close()
        
    # sum over all snap sectors
    for precursor in precursor_lst:
        delta_emission_dict[precursor] = sum(delta_emission_dict[precursor], axis=0)
              
    return delta_emission_dict

# function definition of source receptor model
def module1(path_emission_cdf, path_area_cdf, path_reduction_txt, path_model_cdf, path_result_cdf, *progresslog):

    # check if a progess log file was passed as argument
    if progresslog:
        progress_dict = read_progress_log(progresslog[0])
        write_netcdf_output = False
    else:
        progress_dict = {'start': 0.0, 'divisor': 1.0}
        write_netcdf_output = True
    
    # read the model netcdf
    # ---------------------
    rootgrp = Dataset(path_model_cdf, 'r')
    longitude_array = rootgrp.variables['lon'][0, :]
    latitude_array = rootgrp.variables['lat'][:, 0]
    n_lon = len(longitude_array)  # len(rootgrp.dimensions['longitude'])
    n_lat = len(latitude_array)  # len(rootgrp.dimensions['latitude'])  
    inner_radius = int(getattr(rootgrp, 'Radius of influence'))
    precursor_lst = getattr(rootgrp, 'Order_Pollutant').split(', ')
    alpha = rootgrp.variables['alpha'][:, :, :]    
    omega = rootgrp.variables['omega'][:, :, :] 
    flatWeight = rootgrp.variables['flatWeight'][:, :, :]

    # put alpha and omega in a dictionary
    alpha_dict = {}
    omega_dict = {}
    flatWeight_dict = {}
    for i in range(len(precursor_lst)):
        alpha_dict[precursor_lst[i]] = alpha[i, :, :]
        omega_dict[precursor_lst[i]] = omega[i, :, :]
        flatWeight_dict[precursor_lst[i]] = flatWeight[i, :, :]
        
    # close model netcdf
    rootgrp.close()

    # calculate the delta emissions, dictionary per pollutant a matrix of dimension n_lat x n_lon
    delta_emission_dict = create_delta_emission(path_emission_cdf, precursor_lst, path_area_cdf, path_reduction_txt, path_result_cdf, write_netcdf_output)
    
    # make a window
    window = create_window(inner_radius)
    (n_lon_inner_win, n_lat_inner_win) = window.shape
    
    # create flat window and a inner window
    borderweight = window[inner_radius, 0]
#     print(flat_var_limit)
#     flatweight = (window * (window < flat_var_limit)).sum() / (window < flat_var_limit).sum()
#     print(flatweight)
    
#     # flat_window = ones((n_lat_win, n_lon_win)) * flatweight
#     n_inner = 2 * inner_radius + 1
#     inner_window = zeros((n_inner, n_inner))
#     for iw in range(n_inner):
#         for jw in range(n_inner):
#             cell_dist = 1 / (1 + sqrt((iw - inner_radius) ** 2 + (jw - inner_radius) ** 2))
#             if (cell_dist < flat_var_limit):
#                 inner_window[iw, jw] = flatweight
#             else:
#                 inner_window[iw, jw] = cell_dist
    
#     inner_window = window[(radius - inner_radius):(radius + inner_radius + 1), (radius - inner_radius):(radius + inner_radius + 1)] 
    
    window_ones = ones(window.shape)
    for i in range(n_lat_inner_win):
        for j in range(n_lon_inner_win):
            if window[i,j] < borderweight:
                window[i,j] = 0
                window_ones[i,j] = 0
       
    pad_delta_emission_dict = {}
    for precursor in precursor_lst:
        pad_delta_emission_dict[precursor] = lib.pad(delta_emission_dict[precursor], inner_radius, 'constant', constant_values=0)
    
    # apply source receptor relationships
    # -----------------------------------
    last_progress_print = time()
#     calculate weighted emissions for all precursors
#     norm_delta_conc = zeros((n_lat, n_lon))
    delta_conc = zeros((n_lat, n_lon))
    cell_counter = 0
    n_cell = n_lat * n_lon
    
    # dictionary with sum of emissions over full domain per precursor
    sum_emissions_flat = {}
    for precursor in precursor_lst:
        sum_emissions_flat[precursor] = delta_emission_dict[precursor].sum()   
    
    for ie in range(n_lat):
        if (time() - last_progress_print) > 1:
            if progress_dict['start'] >= 0:
                progress = progress_dict['start'] + float(cell_counter) / float(n_cell) * 100 / progress_dict['divisor']
                sys.stdout.write('\r')
                sys.stdout.flush()
                sys.stdout.write('progress:%f\r' % progress)
                sys.stdout.flush()
                last_progress_print = time()
            
        for je in range(n_lon):
            for precursor in precursor_lst:
                # apply averaging window
                alpha_ij = alpha_dict[precursor][ie, je]
                omega_ij = omega_dict[precursor][ie, je]
                flatWeight_ij = flatWeight_dict[precursor][ie, je]
                
                if not(isnan(alpha_ij)):
                    # update or recalculate flat weighted emissions for each precursor
#                     if sum_emissions_flat[precursor] == None:
#                         sum_emissions_flat[precursor] = pad_delta_emission_dict[precursor][ie:(ie + n_lon_outer_win), je:(je + n_lat_outer_win)].sum()
#                     else:
                        # update flat weight emissions
#                         sum_emissions_flat[precursor] -= pad_delta_emission_dict[precursor][ie:(ie + n_lon_outer_win), (je - 1)].sum()
#                         sum_emissions_flat[precursor] += pad_delta_emission_dict[precursor][ie:(ie + n_lon_outer_win), (je + n_lat_outer_win - 1)].sum()
                    
                    # apply the weight to the flat weighted emissions
                    weighted_emissions_flat = flatWeight_ij * sum_emissions_flat[precursor]  
                    
                    # calculate the inner variable weighted emissions
#                     ring = outer_radius - inner_radius
#                     emissions_centre = pad_delta_emission_dict[precursor][(ie + ring):(ie + n_lon_outer_win - ring), (je + ring):(je + n_lat_outer_win - ring)]
                    emissions_centre = pad_delta_emission_dict[precursor][ie:(ie + n_lon_inner_win), je:(je + n_lat_inner_win)]
                    
                    # weighted_emissions_centre = (power(weights_centre, omega_ij) * emissions_centre).sum()
                    weighted_emissions_centre = ((power(window, omega_ij) - window_ones * flatWeight_ij) * emissions_centre).sum()
                    # to avoid that the little triangles in the 4 corners of the centre area are not counted
                    # weighted_emissions_centre[weighted_emissions_centre < 0] = 0
                    # sum the contribution of the precursor
                    delta_conc[ie, je] = delta_conc[ie, je] + alpha_ij * (weighted_emissions_centre + weighted_emissions_flat)
            
#                 else:
                    # reset the sum after Nan alphas
                    # sum_emissions_flat[precursor] = None
            
            cell_counter += 1


    # create a result netcdf 
    # -----------------------
    if write_netcdf_output == True:
        filename_result_cdf = path_result_cdf + 'delta_concentration.nc'
        rootgrp = Dataset(filename_result_cdf, 'w', format='NETCDF3_CLASSIC')
        
        # create dimensions in the netcdf file
        rootgrp.createDimension('latitude', n_lat)
        rootgrp.createDimension('longitude', n_lon)
        latitudes = rootgrp.createVariable('latitude', 'f4', ('latitude',))
        latitudes.units = "degrees_north"
        longitudes = rootgrp.createVariable('longitude', 'f4', ('longitude',))
        longitudes.units = "degrees_east"
        latitudes[:] = latitude_array
        longitudes[:] = longitude_array
    
        # create delta concentration data
        delta_conc_pm25 = rootgrp.createVariable('delta_concentration', 'f4', ('latitude', 'longitude',))
        delta_conc_pm25.units = 'ug/m3'
        delta_conc_pm25[:] = delta_conc
        
        rootgrp.close()
        
    # create a results object
    mod1_res = {}
    mod1_res['delta_conc'] = delta_conc
    mod1_res['delta_emis_dict'] = delta_emission_dict
    mod1_res['n_lat'] = n_lat
    mod1_res['n_lon'] = n_lon
    mod1_res['latitude_array'] = latitude_array
    mod1_res['longitude_array'] = longitude_array
    
    
    return mod1_res

if __name__ == '__main__':
    
    # module 1 test inputs
    module = 1
    # if it doesn't exist strart=0 and dividsor=1
    progresslog = 'input/progress.log'
    
    # run module 1 without progress log
    start = time()
    start = time()
    emission_1611_test = 'input/20151116_SR_no2_pm10_pm25/BC_emi_PM25_Y.nc'
    area_1611_test = 'input/London_region.nc'
    model_1611_test = 'input/20151116_SR_no2_pm10_pm25/SR_PM25_Y.nc'
    # fullmodel = 'input/fullFunction/SR_PM25_Y_fullFunction.nc'
    output_test_1611 = 'output/'
 
    # run module 1 with progress log
    proglog_filename = path_result_cdf_test + 'proglog'
    write_progress_log(proglog_filename, 25, 2)
    start = time()
    # module1(emission_1611_test, area_1611_test, 'input/user_reduction_snap7.txt', model_1611_test, output_test_1611) #, proglog_filename)
    # debugging Denise's error
    france = '../cities/Paris/Paris_NUTS0.nc'
    module1('../bug_correction_20161205/one_emission_source.nc', france, 'input/user_reduction_snap7.txt', model_1611_test, '../bug_correction_20161205/') #, proglog_filename)
    
    stop = time()
    print('Module 1 run time: %s sec.' % (stop-start))
    remove(proglog_filename)
     
    pass



