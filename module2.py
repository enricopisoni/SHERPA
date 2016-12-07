'''
Created on Jul 13, 2015

ATLAS mode
Input: - baseline emissions (per precursor and cell)
       - a netcdf defining the areas where emission reductions will be applied
       - txt file with reductions (per snap and precursor) to be applied in each nuts
       - model parameters
       - result path
Output: - netcdf with concentration reductions, one layer per NUTS

@author: degraba
'''

from netCDF4 import Dataset
from numpy import zeros, sum, lib, power
from sherpa_auxiliaries import create_emission_reduction_dict, create_emission_dict, \
    create_window, create_reduced_precursor_lst
from sherpa_globals import path_emission_cdf_test, path_nuts0_cdf_test, path_result_cdf_test, \
    path_nuts3_cdf_test
# path_reduction_txt_test, path_model_cdf_test, 
import time
from math import isnan
import sys

def module2(path_emission_cdf, path_nuts_cdf, path_reduction_txt, path_model_cdf, path_result_cdf): 
    
    # read the model netcdf
    #----------------------
    # open file
    rootgrp = Dataset(path_model_cdf, 'r')
    
    longitude_array = rootgrp.variables['lon'][0, :]
    latitude_array = rootgrp.variables['lat'][:, 0]
    n_lon = len(longitude_array)  # len(rootgrp.dimensions['longitude'])
    n_lat = len(latitude_array)  # len(rootgrp.dimensions['latitude'])   
    inner_radius = int(getattr(rootgrp, 'Radius of influence'))
    
    alpha = rootgrp.variables['alpha'][:, :, :]    
    omega = rootgrp.variables['omega'][:, :, :] 
    flatWeight = rootgrp.variables['flatWeight'][:, :, :]
    precursor_lst = getattr(rootgrp, 'Order_Pollutant').split(', ')
    
    alpha_dict = {}
    omega_dict = {}
    flatWeight_dict = {}
    
    for i in range(len(precursor_lst)):
        alpha_dict[precursor_lst[i]] = alpha[i, :, :]
        omega_dict[precursor_lst[i]] = omega[i, :, :]
        flatWeight_dict[precursor_lst[i]] = flatWeight[i, :, :]
    
    # close netcdf
    rootgrp.close()
    
    # create a dictionary with reductions per precursor and macro sector
    emission_reduction_dict = create_emission_reduction_dict(path_reduction_txt)
    
    # create list of reduced precursors
    reduced_precursor_lst = create_reduced_precursor_lst(emission_reduction_dict)

    # read nuts netcdf: values between 0 and 100 indicate which share of the cell has to be taken into account
    rootgrp_nuts = Dataset(path_nuts_cdf, 'r')
    n_nuts = len(rootgrp_nuts.dimensions['nuts_id'])
    
    # declare results array for delta concentration, dimension: n_nuts x lat x lon
    delta_conc_nuts = zeros((n_nuts, n_lat, n_lon))
    
    # create emission dictionary
    emission_dict = create_emission_dict(path_emission_cdf, precursor_lst)

    # make a window
    window = create_window(inner_radius)
    (n_lon_inner_win, n_lat_inner_win) = window.shape
    
    # create flat window and a inner window
    borderweight = window[inner_radius, 0]
    
    for i in range(n_lat_inner_win):
        for j in range(n_lon_inner_win):
            if window[i,j] < borderweight:
                window[i,j] = 0

    # result dictionary with emission reductions per nuts of the reduced precursor
    # only if one precursor is reduced this result has to be stored for potency calculations
    delta_emission_nuts_dict = {}

    # loop over all nuts
    nuts_counter = 0
    for nuts in range(n_nuts):
        
        # create delta emission dictionary for nuts area
        reduction_area = rootgrp_nuts.variables['AREA'][nuts, :, :] / 100
        pad_delta_emission_dict = {}
        # dictionary with sum of emissions over full domain per precursor
        sum_emissions_flat = {}
        
        for precursor in precursor_lst:
            # calculate emission reductions per nuts area and snap sector
            delta_emission_precursor_snap = zeros(emission_dict[precursor].shape)
            for snap in range(1, 11):
                delta_emission_precursor_snap[snap - 1, :, :] = emission_dict[precursor][snap - 1] * reduction_area * emission_reduction_dict[precursor][snap]
            # sum delta emission per snap over all snap sectors
            delta_emission_precursor = sum(delta_emission_precursor_snap, axis=0)
            
            # store delta emissions per nuts and precursor for potency calculation
            if len(reduced_precursor_lst) == 1 and precursor == reduced_precursor_lst[0]:
                delta_emission_nuts_dict[nuts] = delta_emission_precursor
            
            # pad the delta emission arrays with zeros
            pad_delta_emission_dict[precursor] = lib.pad(delta_emission_precursor, inner_radius, 'constant', constant_values=0)

            sum_emissions_flat[precursor] = delta_emission_precursor.sum()   

        # apply source receptor relationships
        # -----------------------------------
        delta_conc = zeros((n_lat, n_lon))
        for ie in range(n_lat):
            for je in range(n_lon):
                # test if the cell overlaps with the nuts sector
                if reduction_area[ie, je] > 0:
                    # apply the correlation between delta_emission and delta concentration
                    for precursor in precursor_lst:
                        alpha_ij = alpha_dict[precursor][ie, je]
                        omega_ij = omega_dict[precursor][ie, je]
                        flatWeight_ij = flatWeight_dict[precursor][ie, je]

                        if not(isnan(alpha_ij)):
                            # apply the weight to the flat weighted emissions
                            weighted_emissions_flat = flatWeight_ij * sum_emissions_flat[precursor]  
                            
                            emissions_centre = pad_delta_emission_dict[precursor][ie:(ie + n_lon_inner_win), je:(je + n_lat_inner_win)]
                            
                            # weighted_emissions_centre = (power(weights_centre, omega_ij) * emissions_centre).sum()
                            weighted_emissions_centre = ((power(window, omega_ij) - flatWeight_ij) * emissions_centre).sum()
                            delta_conc[ie, je] = delta_conc[ie, je] + alpha_ij * (weighted_emissions_centre + weighted_emissions_flat)

        
        # store the result
        delta_conc_nuts[nuts,:,:] = delta_conc
        # print('NUTS %d calculated' % (nuts))
        nuts_counter += 1
        progress = float(nuts_counter) / float(n_nuts) * 100
        sys.stdout.write('\r')
        sys.stdout.flush()
        sys.stdout.write('progress:%f' % progress)
        sys.stdout.flush()
            

    # create a result netcdf 
    # -----------------------
    filename_result_cdf = path_result_cdf + 'delta_concentration_nuts.nc'
    rootgrp = Dataset(filename_result_cdf, 'w', format='NETCDF3_CLASSIC')
    
    # create dimensions in the netcdf file
    rootgrp.createDimension('nuts_id', n_nuts)
    rootgrp.createDimension('latitude', n_lat)
    rootgrp.createDimension('longitude', n_lon)
    nuts_ids = rootgrp.createVariable('nuts_id', 'f4', ('nuts_id',))
    nuts_vector = range(1, n_nuts + 1)
    nuts_ids[:] = nuts_vector
    latitudes = rootgrp.createVariable('latitude', 'f4', ('latitude',))
    latitudes.units = "degrees_north"
    latitudes[:] = latitude_array
    longitudes = rootgrp.createVariable('longitude', 'f4', ('longitude',))
    longitudes.units = "degrees_east"
    longitudes[:] = longitude_array

    # create delta concentration data
    delta_conc_nuts_var = rootgrp.createVariable('delta_concentration_nuts', 'f4', ('nuts_id', 'latitude', 'longitude', ))
    delta_conc_nuts_var.units = "ug/m3"
    delta_conc_nuts_var[:] = delta_conc_nuts
    
    rootgrp.close()
    
    # create a results dictionary
    mod2_res = {}
    mod2_res['delta_conc_nuts'] = delta_conc_nuts
    mod2_res['delta_emis_precursor_nuts'] = delta_emission_nuts_dict
    mod2_res['n_nuts'] = n_nuts
    mod2_res['n_lat'] = n_lat
    mod2_res['n_lon'] = n_lon
    mod2_res['nuts_vector'] = nuts_vector
    mod2_res['latitude_array'] = latitude_array
    mod2_res['longitude_array'] = longitude_array
    
    # close rootgrp_nuts
    rootgrp_nuts.close()
    
    return mod2_res

if __name__ == '__main__':
       
    start = time.time()
    emission_1611_test = 'input/20151116_SR_no2_pm10_pm25/BC_emi_PM25_Y.nc'
    concentration_1611_test = 'input/20151116_SR_no2_pm10_pm25/BC_conc_PM25_Y.nc'
    model_1611_test = 'input/20151116_SR_no2_pm10_pm25/SR_PM25_Y.nc'

    module2(emission_1611_test, path_nuts0_cdf_test, 'input/user_reduction_snap10.txt', model_1611_test, path_result_cdf_test)
    stop = time.time()
    print('\nModule 2 calculation time = %f' % (stop - start))
    pass
