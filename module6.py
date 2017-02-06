'''
Created on Jun 23, 2015

!!!!!!!! this script is in the Accelerator Branch !!!!!!!!!!!!!

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
from numpy import lib, zeros, sum, power, ones, array, sqrt, savetxt
from math import isnan
# path_emission_cdf_test, path_area_cdf_test, path_reduction_txt_test, path_model_cdf_test,
from time import time
import sys
from sherpa_globals import alpha_potency
from sherpa_auxiliaries import create_emission_reduction_dict, create_emission_dict, create_window, deltaNOx_to_deltaNO2

# extra variables for the accelerator
n_lowres = 3   # width and height of the aggregation of emmissions, this number has to be odd!

class Emissions:
    def __init__(self, n_lowres, delta_emission_dict):
        self.dem = delta_emission_dict
        self.n_lat = delta_emission_dict['PPM'].shape[0]
        self.n_lon = delta_emission_dict['PPM'].shape[1]
        self.n_lowres = n_lowres
        self.agg_radius = (n_lowres - 1) / 2
        self.aggregated_emissions = {}  
        self.emissions_hires = {}    
    
    # generates or looks up hi and low resolution emissions for a precursor and coordinate
    def getHiLowEmis(self, precursor, i_lat_target, i_lon_target):
        
        i_lat_mod = i_lat_target % n_lowres
        i_lon_mod = i_lon_target % n_lowres
        i_tuple_mod = (i_lat_mod, i_lon_mod)
        
        # check if for the modulus of the row-column index tuple some aggregations are done
        if not(i_tuple_mod in self.aggregated_emissions.keys()):
            self.aggregated_emissions[i_tuple_mod] = {}
            self.emissions_hires[i_tuple_mod] = {}
            
            # define the breaks of the aggregation blocks along latitudes and longitudes
            if i_lon_mod == self.agg_radius:
                lon_breaks = [0]
            else:
                lon_breaks = [0, i_lon_mod - self.agg_radius]
            # as long as the last element is smaller than n_lon, add a break
            while lon_breaks[-1] < self.n_lon: 
                next_break = lon_breaks[-1] + self.n_lowres   # last element + width of aggregation window
                if next_break < self.n_lon:
                    lon_breaks.append(next_break)
                else:
                    lon_breaks.append(self.n_lon)
             
            # same for latitudes
            if i_lat_mod == self.agg_radius:
                lat_breaks = [0]
            else:
                lat_breaks = [0, i_lat_mod - self.agg_radius]
            while lat_breaks[-1] < self.n_lat: 
                next_break = lat_breaks[-1] + self.n_lowres   # last element + width of aggregation window
                if next_break < self.n_lat:
                    lat_breaks.append(next_break)
                else:
                    lat_breaks.append(self.n_lat)
            
            self.aggregated_emissions[i_tuple_mod]['n_lat_agg'] = len(lat_breaks) - 1
            self.aggregated_emissions[i_tuple_mod]['n_lon_agg'] = len(lon_breaks) - 1
            self.aggregated_emissions[i_tuple_mod]['lat_breaks'] = lat_breaks
            self.aggregated_emissions[i_tuple_mod]['lon_breaks'] = lon_breaks
            
        if not(precursor in (self.aggregated_emissions[i_tuple_mod]).keys()):
            n_lat_agg = self.aggregated_emissions[i_tuple_mod]['n_lat_agg']
            n_lon_agg = self.aggregated_emissions[i_tuple_mod]['n_lon_agg']
            lat_breaks = self.aggregated_emissions[i_tuple_mod]['lat_breaks']
            lon_breaks = self.aggregated_emissions[i_tuple_mod]['lon_breaks']
            
            self.aggregated_emissions[i_tuple_mod][precursor] = zeros((n_lat_agg, n_lon_agg))
            self.emissions_hires[i_tuple_mod][precursor] = zeros((self.n_lat, self.n_lon))
            for i_lat_agg in range(n_lat_agg):
                for i_lon_agg in range(n_lon_agg):
                    # select the block to be aggregated from delta_emissions[precursor]
                    emis2aggregate = self.dem[precursor][lat_breaks[i_lat_agg]:lat_breaks[i_lat_agg + 1], lon_breaks[i_lon_agg]:lon_breaks[i_lon_agg + 1]]
                    sum_emissions = emis2aggregate.sum()
                    n_agg = emis2aggregate.size
                    self.aggregated_emissions[i_tuple_mod][precursor][i_lat_agg, i_lon_agg] = sum_emissions
                    # high resolution delta emissions are delta emissions minus average aggregated emissions
                    self.emissions_hires[i_tuple_mod][precursor][lat_breaks[i_lat_agg]:lat_breaks[i_lat_agg + 1], lon_breaks[i_lon_agg]:lon_breaks[i_lon_agg + 1]] = emis2aggregate - sum_emissions / n_agg
#                     print(self.aggregated_emissions[i_tuple_mod][precursor])
#                     print(self.emissions_hires[i_tuple_mod][precursor])
        
        res_dict = {'emis_lowres': self.aggregated_emissions[i_tuple_mod][precursor], 'emis_hires': self.emissions_hires[i_tuple_mod][precursor]}
        
        return res_dict
                        
# Window class that returns aggregated weighting windows for a given omega
class OmegaWindows:
    def __init__(self, aggreg_size, hires_window_size, lowres_window_size):
        self.aggreg_size = aggreg_size            # number of rows/columns that are aggregated
        self.aggreg_radius = (aggreg_size - 1) / 2
        
        self.hires_window_size = hires_window_size
        self.hires_window_radius = (hires_window_size - 1) / 2
        self.hires_inverse_distance = zeros((self.hires_window_size, self.hires_window_size))
        for iw in range(self.hires_window_size):
            for jw in range(self.hires_window_size):
                cell_dist = sqrt((float(iw - self.hires_window_radius)) ** 2 + (float(jw - self.hires_window_radius)) ** 2) 
                self.hires_inverse_distance[iw, jw] = 1 / (1 + cell_dist)  
        
        self.lowres_window_size = lowres_window_size
        self.lowres_window_radius = (lowres_window_size - 1) / 2
        self.lowres_inverse_distance = zeros((self.lowres_window_size, self.lowres_window_size))
        for iw in range(self.lowres_window_size):
            for jw in range(self.lowres_window_size):
                cell_dist = sqrt((float(iw - self.lowres_window_radius)) ** 2 + (float(jw - self.lowres_window_radius)) ** 2) 
                self.lowres_inverse_distance[iw, jw] = 1 / (1 + cell_dist)  
        
        self.aggreg_lowres_window_size = self.lowres_window_size / self.aggreg_size
        
        self.hires_omega_windows = {} 
        self.lowres_omega_windows = {}

    def getOmegaWindow(self, omega):
        if omega in self.hires_omega_windows.keys():
            hires_omega_window = self.hires_omega_windows[omega]
            lowres_omega_window = self.lowres_omega_windows[omega]
        else:
            self.hires_omega_windows[omega] = power(self.hires_inverse_distance, omega)
            
            # aggregate the omega window in blocks of n_lowres x n_lowres cells
            lowres_omega_window = zeros((self.aggreg_lowres_window_size , self.aggreg_lowres_window_size))
            for iw in range(self.aggreg_lowres_window_size):
                for jw in range(self.aggreg_lowres_window_size):
                    lowres_omega_window[iw, jw] = (power(self.lowres_inverse_distance[(iw * self.aggreg_size):((iw + 1) * self.aggreg_size), (jw * self.aggreg_size):((jw + 1) * self.aggreg_size)], omega)).mean()
            self.lowres_omega_windows[omega] = lowres_omega_window
            
        res_dict = {'hires_omega_window': self.hires_omega_windows[omega], 'lowres_omega_window': lowres_omega_window} 
        
        return res_dict 
            


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
            
        # create the aggregated emissions at low resolution
        pad_delta_lowres_emission_dict = {}
        for precursor in precursor_lst:
            pad_delta_lowres_emission_dict[precursor] = lib.pad(delta_emission_dict[precursor], (n_lowres - 1) / 2, 'constant', constant_values=0)
            
       
        
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
            # flatWeight_ij = flatWeight_dict[precursor][i_lat_target, i_lon_target]
            
            if not(isnan(alpha_ij)):
                
                # apply the weight to the flat weighted emissions
                # weighted_emissions_flat = flatWeight_ij * sum_emissions_flat[precursor]  
                
                emissions_centre = pad_delta_emission_dict[precursor][i_lat_target:(i_lat_target + n_lon_inner_win), i_lon_target:(i_lon_target + n_lat_inner_win)]
                
                # weighted_emissions_centre = (power(weights_centre, omega_ij) * emissions_centre).sum()
                # weighted_emissions_centre = ((power(window, omega_ij)) * emissions_centre).sum()
                weighted_emissions_lowres = NaN
                weighted_emissions_hires = NaN
                delta_conc[nuts_code] = delta_conc[nuts_code] + alpha_ij * (weighted_emissions_lowres + weighted_emissions_hires)
                
        # In the case of NOx the NO2 concentrations have to be calculated with the NO2 fraction correlation
        if (path_model_cdf.find('NO2eq') > -1):
            rootgrp = Dataset(path_base_conc_cdf, 'r')
            base_conc_nox = array(rootgrp.variables['conc'][i_lat_target, i_lon_target])  
            base_conc_no2 = array(rootgrp.variables['NO2'][i_lat_target, i_lon_target])
            rootgrp.close() 
            delta_conc[nuts_code] = deltaNOx_to_deltaNO2(delta_conc[nuts_code], base_conc_nox, base_conc_no2)
           
    
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
    
    # test window function
    testwin = OmegaWindows(3, 3, 9)
    # print(testwin.hires_inverse_distance)
    savetxt("C:/temp/hires_inverse_distance.csv", testwin.hires_inverse_distance, delimiter="\t")
    savetxt("C:/temp/lowres_inverse_distance.csv", testwin.lowres_inverse_distance, delimiter="\t")
    # print(testwin.lowres_inverse_distance)
    omegawindows = testwin.getOmegaWindow(1)
    print(omegawindows['hires_omega_window'])
    savetxt("C:/temp/hires_omega_window.csv", omegawindows['hires_omega_window'], delimiter="\t")
    savetxt("C:/temp/lowres_omega_window.csv", omegawindows['lowres_omega_window'], delimiter="\t")
    print(omegawindows['lowres_omega_window'])
    
    
    # test the Emissions class
    delta_emission_dict = {}
    delta_emission_dict['PPM'] = ones((380,480))
    delta_emission_dict['PPM'][5,5] = 10
    delta_emission_dict['PPM'][0,1] = 10
    print(delta_emission_dict['PPM'])
    print(delta_emission_dict['PPM'].sum())
    
    test = Emissions(3, delta_emission_dict)
    start = time()
    HiLowEmis_dict = test.getHiLowEmis('PPM', 5, 5)
    HiEmis = HiLowEmis_dict['emis_hires']
    LowEmis = HiLowEmis_dict['emis_lowres']
    print(HiEmis.sum())
    print(LowEmis.sum())
    stop = time()
    print('First call: %f' % (stop - start))
    
    start = time()
    emis_low_res = test.getLowResEmis('PPM', 5, 5)
    stop = time()
    print('Second call: %f' % (stop - start))
    
    print(emis_low_res)
    emis_low_res = test.getLowResEmis('PPM', 1, 1)
    print(emis_low_res)
    
#     # module 1 test inputs
#     module = 1
#     # if it doesn't exist strart=0 and dividsor=1
#     progresslog = 'input/progress.log'
#     
#     # run module 1 without progress log
#     start = time()
#     # emissions = 'input/20151116_SR_no2_pm10_pm25/BC_emi_NO2_Y.nc'
#     emissions = 'input/20151116_SR_no2_pm10_pm25/BC_emi_PM25_Y.nc'
#     nuts2_netcdf = 'input/EMI_RED_ATLAS_NUTS2.nc'
#     target_cell_lat = 45.46         # Milan
#     target_cell_lon = 9.19          # Milan
#     path_reduction_txt = 'input/user_reduction_all50.txt'
#     # base_conc_cdf = 'input/20151116_SR_no2_pm10_pm25/BC_conc_NO2_NO2eq_Y_mgm3.nc'
#     base_conc_cdf = 'input/20151116_SR_no2_pm10_pm25/BC_conc_PM25_Y.nc'
#     # model_NO2eq = 'input/20151116_SR_no2_pm10_pm25/SR_NO2eq_Y.nc'
#     model_PM25old = 'input/20151116_SR_no2_pm10_pm25/SR_PM25_Y.nc'
#     model_PM25new = 'input/20151116_SR_no2_pm10_pm25/SR_PM25_Y_prctiles.nc'
#     output_path = 'output/NO2eq/Milan/'
#      
#     # run module 1 with progress log
#     start = time()
#     module6(emissions, nuts2_netcdf, target_cell_lat, target_cell_lon, path_reduction_txt, base_conc_cdf, model_PM25new, output_path)
#     # print(DC)
#     stop = time()
#     print('Module 6 run time: %s sec.' % (stop-start))
     
    pass




