'''
Created on Jun 23, 2015

Module 6 calculates for 1 cell the concentration change due to a 50 percent reductions in 
the snap sectors defined in the input file 'path_reduction_txt'. Emission are reduced in
each NUTS area in the input file 'path_area_cdf'
There are 2 outputs:
- a text file with all nuts codes and the DC/C/alpha (relative potential) as percent due to a reduction in that nuts area
- a map where each nuts has the value of the concentration change it causes in the target cell 

for compatibility the header is 'potency' in the output txt

@author: degraba
'''

# imports
from netCDF4 import Dataset
from numpy import lib, zeros, sum, power, ones
from math import isnan
# path_emission_cdf_test, path_area_cdf_test, path_reduction_txt_test, path_model_cdf_test,
from time import time
import sys
from sherpa_globals import alpha_potency
from sherpa_auxiliaries import create_emission_reduction_dict, create_emission_dict, create_window

# function that applies reductions per snap sector and precursor to the emission netcdf
def create_delta_emission(path_emission_cdf, precursor_lst, reduction_area_array, path_reduction_txt):
        
    # create a dictionary with reductions per precursor and macro sector
    emission_reduction_dict = create_emission_reduction_dict(path_reduction_txt)
    
    # open the emission netcdf
    emission_dict = create_emission_dict(path_emission_cdf, precursor_lst)
       
    # calculate a dictionary with the emission reductions per pollutant, macrosector and position
    delta_emission_dict = {}
    for precursor in precursor_lst:
        delta_emission_dict[precursor] = zeros(emission_dict[precursor].shape)
        # calculate the emission reduction
        # reductions are positive!
        # make the sum over all snap sectors
        for snap in range(1, 11):
            delta_emission_dict[precursor][snap - 1, :, :] = emission_dict[precursor][snap - 1] * reduction_area_array * emission_reduction_dict[precursor][snap]
        
    # sum over all snap sectors
    for precursor in precursor_lst:
        delta_emission_dict[precursor] = sum(delta_emission_dict[precursor], axis=0)
              
    return delta_emission_dict

# function definition of source receptor model
def module6(path_emission_cdf, path_area_cdf, target_cell_lat, target_cell_lon, path_reduction_txt, path_base_conc_cdf, path_model_cdf, path_result_cdf):
    
    # read the model netcdf
    # ---------------------
    rootgrp = Dataset(path_model_cdf, 'r')
    longitude_array = rootgrp.variables['lon'][0, :]
    latitude_array = rootgrp.variables['lat'][:, 0]
    n_lon = len(longitude_array)  
    n_lat = len(latitude_array)  
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
    
    # open the area netcdf and get lat lon indexes of target cell
    #-------------------------------------------------------------
    rootgrp_nuts = Dataset(path_area_cdf, 'r')
    n_nuts = len(rootgrp_nuts.dimensions['nuts_id'])
    nuts_codes_raw = rootgrp_nuts.variables['NUTS'][:]
    nuts_codes = []
    for i_code in range(len(nuts_codes_raw)):
        code = ''
        for letter in nuts_codes_raw[i_code]:
            code = code + letter
        nuts_codes.append(code)
    
    # convert latitude and longitude string in float
    target_cell_lat = float(target_cell_lat)
    target_cell_lon = float(target_cell_lon)
    
    # get row index of latitude and col index of longitude
    i_lat_target = 0
    lat_error = float('inf')
    for i in range(len(latitude_array)):
        lat_dist = abs(target_cell_lat - latitude_array[i])
        if lat_dist < lat_error:
            lat_error = lat_dist
            i_lat_target = i
    
    i_lon_target = 0
    lon_error = float('inf')
    for i in range(len(longitude_array)):
        lon_dist = abs(target_cell_lon - longitude_array[i])
        if lon_dist < lon_error:
            lon_error = lon_dist
            i_lon_target = i
    
    # read base concentrations and extract base case concentration in the target cell
    # -------------------------------------------------------------------------------
    rootgrp = Dataset(path_base_conc_cdf, 'r')
    target_conc_basecase = rootgrp.variables['conc'][i_lat_target, i_lon_target]
    # close model netcdf
    rootgrp.close()
    
    # make a window
    window = create_window(inner_radius)
    (n_lon_inner_win, n_lat_inner_win) = window.shape
    
    # create flat window and a inner window
    borderweight = window[inner_radius, 0]
    
    window_ones = ones(window.shape)
    for i in range(n_lat_inner_win):
        for j in range(n_lon_inner_win):
            if window[i,j] < borderweight:
                window[i,j] = 0
                window_ones[i,j] = 0
    
    delta_conc = {} 
    DC_target_arrray = zeros((n_lat, n_lon))
        
    # loop over all nuts in 
    for nuts_id in range(n_nuts):
        # initialize delta_conc
        nuts_code = nuts_codes[nuts_id]
        delta_conc[nuts_code] = 0
        # print the progress
        progress = float(nuts_id) / float(n_nuts) * 100
        sys.stdout.write('\r')
        sys.stdout.flush()
        sys.stdout.write('progress:%f\r' % progress)
        sys.stdout.flush()
    
        reduction_area_array = rootgrp_nuts.variables['AREA'][nuts_id,:,:] / 100.0
        
        # calculate the delta emissions, dictionary per pollutant a matrix of dimension n_lat x n_lon
        delta_emission_dict = create_delta_emission(path_emission_cdf, precursor_lst, reduction_area_array, path_reduction_txt)
       
        pad_delta_emission_dict = {}
        for precursor in precursor_lst:
            pad_delta_emission_dict[precursor] = lib.pad(delta_emission_dict[precursor], inner_radius, 'constant', constant_values=0)
        
        # apply source receptor relationships
        # -----------------------------------
        
        # dictionary with sum of emissions over full domain per precursor
        sum_emissions_flat = {}
        for precursor in precursor_lst:
            sum_emissions_flat[precursor] = delta_emission_dict[precursor].sum()   
                    
        for precursor in precursor_lst:
            # apply averaging window
            alpha_ij = alpha_dict[precursor][i_lat_target, i_lon_target]
            omega_ij = omega_dict[precursor][i_lat_target, i_lon_target]
            flatWeight_ij = flatWeight_dict[precursor][i_lat_target, i_lon_target]
            
            if not(isnan(alpha_ij)):
                
                # apply the weight to the flat weighted emissions
                weighted_emissions_flat = flatWeight_ij * sum_emissions_flat[precursor]  
                
                emissions_centre = pad_delta_emission_dict[precursor][i_lat_target:(i_lat_target + n_lon_inner_win), i_lon_target:(i_lon_target + n_lat_inner_win)]
                
                # weighted_emissions_centre = (power(weights_centre, omega_ij) * emissions_centre).sum()
                weighted_emissions_centre = ((power(window, omega_ij) - window_ones * flatWeight_ij) * emissions_centre).sum()
                delta_conc[nuts_code] = delta_conc[nuts_code] + alpha_ij * (weighted_emissions_centre + weighted_emissions_flat)
    
        # create an output map with in each nuts the DC in the target cell
        DC_target_arrray = DC_target_arrray + delta_conc[nuts_code] * reduction_area_array
        
    # close nuts cdf
    rootgrp_nuts.close()
    
    # sort nuts codes from delta_conc from high to low delta conc
    sorted_nuts_codes = sorted(delta_conc, key=lambda i: delta_conc[i], reverse=True) 
    
    # write the result to a netcdf file
    path_DC_target_cdf = path_result_cdf + 'radius_result.nc'
    rootgrp = Dataset(path_DC_target_cdf, 'w', format = 'NETCDF3_CLASSIC')
    rootgrp.createDimension('latitude', n_lat)
    rootgrp.createDimension('longitude', n_lon)
    latitudes = rootgrp.createVariable('latitude', 'f4', ('latitude',))
    latitudes.units = "degrees_north"
    latitudes[:] = latitude_array
    longitudes = rootgrp.createVariable('longitude', 'f4', ('longitude',))
    longitudes.units = "degrees_east"
    longitudes[:] = longitude_array
    area = rootgrp.createVariable('AREA', 'f4', ('latitude', 'longitude',))
    area[:] = DC_target_arrray
    rootgrp.close()
    
    # write a result file
    f_res = open(path_result_cdf + 'radius_result.txt', 'w')
    f_res.write('nuts_code\t%\n')
    for nuts_code in sorted_nuts_codes:
        f_res.write('%s\t%e\n' % (nuts_code, delta_conc[nuts_code] / target_conc_basecase / (alpha_potency / 100) * 100)) # rel potential in percentage
    f_res.close()

    # return delta_conc

if __name__ == '__main__':
    
    # module 1 test inputs
    module = 1
    # if it doesn't exist strart=0 and dividsor=1
    progresslog = 'input/progress.log'
    
    # run module 1 without progress log
    start = time()
    start = time()
    emission_1611_test = 'input/20151116_SR_no2_pm10_pm25/BC_emi_PM25_Y.nc'
    nuts2_netcdf = 'input/EMI_RED_ATLAS_NUTS2.nc'
    target_cell_lat = 51.51     # 51.51
    target_cell_lon = -0.13  # 9.19      #-0.13
    path_base_conc_cdf = 'input/20151116_SR_no2_pm10_pm25/BC_conc_PM25_Y.nc'
    model_1611_test = 'input/20151116_SR_no2_pm10_pm25/SR_PM25_Y.nc'
    path_result_cdf = 'output/'
     
    # run module 1 with progress log
    start = time()
    module6(emission_1611_test, nuts2_netcdf, target_cell_lat, target_cell_lon, 'input/user_reduction_snap7.txt', path_base_conc_cdf, model_1611_test, path_result_cdf)
    # print(DC)
    stop = time()
    print('Module 6 run time: %s sec.' % (stop-start))
     
    pass




