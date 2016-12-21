# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 16:26:42 2016

@author: peduzem

INPUT:   - netcdf with concentration changes per pollutant and cell
         - total emissions CLE scenario (2010)
OUTPUT:  - netcdf with the MOLL for the CLE scenario
         - graphs with the PM2.5 concentraion and the MOLL for the CLE scenario

... (to be continued)

First trial for benefit assessment based on Enrico's matlab file"""

# imports

# for importing matlab files
import scipy.io as sio
# for using netCDF files
from netCDF4 import Dataset
# for scientific operators
import numpy as np
# for plotting
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
from sherpa_auxiliaries import create_emission_reduction_dict, create_emission_dict, create_window, read_progress_log, write_progress_log
# Definition of vehicle class
import measures as m
#from measures import av_dsl_pc_eu6
def write_reductions(path_reduction_txt, red):
    """ 
    Function to write the reduction file. 
    At the moment only for MS7 (**will be extended to the others)
    input:
        path_reduction_txt: path to the reduction file (update with sherpa globals)
        red: array with the reduction (>0) % per pollutant
    """ 
    text_file = open(path_reduction_txt, "w")
    text_file.write(" POLL	MS1	MS2	MS3	MS4	MS5	MS6	MS7	MS8	MS9	MS10\n")
    text_file.write(" NOx	0	0	0	0	0	0	%s	0	0	0\n" % red[0])
    text_file.write(" NMVOC	0	0	0	0	0	0	%s	0	0	0\n" % red[1])
    text_file.write(" NH3	0	0	0	0	0	0	%s	0	0	0\n" % red[2])
    text_file.write(" PPM	0	0	0	0	0	0	%s	0	0	0\n" % red[3])
    text_file.write(" SOx	0	0	0	0	0	0	%s	0	0	0\n" % red[4])  
    text_file.close()

def calc_impacts(deaths, pop, pop30plus, pyll, pm25_cle, area_area):
    # -------------------------------------------------------------------------
    # (Cost) Benefit Analysis 
    # -------------------------------------------------------------------------
    '''
    Methodology one (used by Enrico in previous work,  )
    Anenberg, S.C. et al., 2010. Environmental Health Perspectives, 118(9),
    pp.1189–1195. Available at: http://ehp.niehs.nih.gov/0901220
    [Accessed December 13, 2016].
    
    input:
        deaths: total deaths of population over 30
        pop: total population over 30
        pop30plus: distribution of population over 30
        pyll: potential years of life loss, PYLL per 100 000 -30+ average EU28
        pm25_cle: current pollutant concentration (units**) 
        ** will modify this to read directly module 1 output
        area_area: area of interest 
    Output: 
        deltayoll: delta yoll per grid cell between scenario and c
        , deltayll_reg, deltayoll_tot
        
    '''
    # open SCEnario (2030, MFR) and read PM2.5 values
    # **(This will not be necessary as it will be the result of module1)
    nc_file2 = 'netcdf/sherpaRiat-inerisEMI-pm25AQI-2030MFR.nc'
    fh_sce = Dataset(nc_file2, mode='r')
    pm25_sce = fh_sce.variables['PM25'][:]
    fh_sce.close()
    
    drate = deaths/pop
    
    # Incidence rate (as calculated by Enrico, not used here)
    ir= pyll / 100000 * pop / deaths

    # crf derived by RR in Anenberg, S.C. et al., 2010. 
    crf = 0.006

    # Mortality in baseline (CLE) and scenario. 
    mortcle = pm25_cle*crf*pop30plus*drate
    mortsce = pm25_sce*crf*pop30plus*drate
    
    # Years of life loss (yll) and months of life loss (mll) in baseline (CLE) and scenario. 
    yllcle = mortcle * pyll/100000 / drate
    yllsce = mortsce * pyll/100000 / drate
    mollcle = yllcle*12
    mollsce = yllsce*12
    
    # Years of life loss (yll) and months of life loss (mll) in baseline (CLE) and scenario. 
    deltamoll = mollsce-mollcle
    deltayoll = yllsce-yllcle
    
    # Calculate delta moll and yll in the selected region
    deltamoll_reg = np.sum(np.multiply(deltamoll, area_area)/100)
    deltayll_reg = deltamoll_reg / 12
    
    # Calculate delta moll and yll in total
    deltamoll_tot = np.sum(deltamoll)
    deltayoll_tot = deltamoll_tot / 12
    
    return deltayoll, deltayll_reg, deltayoll_tot

# main program starts here
if __name__ == '__main__':
    # -------------------------------------------------------------------------
    # Preparing the data (**will be translated to a def)
    # -------------------------------------------------------------------------

    # read the model netcdf (copied from Model1)
    rootgrp = Dataset('netcdf/20151116_SR_no2_pm10_pm25/SR_PM25_Y.nc', 'r')
    precursor_lst = getattr(rootgrp, 'Order_Pollutant').split(', ')
    # create a dictionary with emissions per precursor, macrosector and postion (lat, lon)
    emission_dict = create_emission_dict('netcdf/20151116_SR_no2_pm10_pm25/BC_emi_PM25_Y.nc', precursor_lst)

    # open 2010 CLE and read PM2.5 concentration values
    nc_file = 'netcdf/sherpaRiat-inerisEMI-pm25AQI-2010CLE.nc'
    fh = Dataset(nc_file, mode='r')
    lons = fh.variables['lon'][:]
    lats = fh.variables['lat'][:]
    pm25_cle = fh.variables['PM25'][:]
    pm25_units_cle = fh.variables['PM25'].units
    fh.close()      
       
    # Area of each cell in the domain    
    nc_file_area = 'netcdf/netcdf_with_area/JRC01.nc'
    fh_area_cells = Dataset(nc_file_area, mode='r')
    # lons = fh.variables['lon'][:]
    # lats = fh.variables['lat'][:]
    area_cell = fh_area_cells.variables['surface'][:]
    area_units = fh_area_cells.variables['surface'].units
    fh_area_cells.close() 
    
    
    # open area of interest
    nc_file_reg = 'netcdf/London_region.nc'
    fh_area = Dataset(nc_file_reg, mode='r')
    area_area = fh_area.variables['AREA'][:]
    fh_area.close() 
    
    # Load population file from Marco (treated in "regridPop" routine from Enrico).
    # The matlab file contains a structure called Anew,
    # I carried out the same operations as Enrico's matlab file.
    A = sio.loadmat('population.mat')
    Anew = A['Anew']
    Anew_T = Anew[np.newaxis]
    popall=np.fliplr(Anew_T)
    # total pop according to Marco (all age groups)
    sumpop=np.sum(popall)

    # Data for baseline  population
    # ICD codes: ICD-10: A00-B99,C00-D48,D50-D89,E00-E88,F01-F99,G00-G98,
    # H00-H59,H60-H93,I00-I99,J00-J98,K00-K92,L00-L98,M00-M99,N00-N98,
    # O00-O99,P00-P96,Q00-Q99,R00-R99
    # Age: '30 - 85 +'									
    # Sex: Both									
    # http://data.euro.who.int/dmdb/
    # [Accessed December 13, 2016].
    # Potential Years of life loss, PYLL per 100 000 -30+ average EU28
    pyll = 4694.26465
    # TOTAL POP above 30
    pop = 331923577
    # TOTAL DEATHS  above 30
    deaths = 4639244
    
    # Distribution of the population over 30
    # HP: it is the same distribution as all the population (as provided by Marco)
    pop30plus = (popall/sumpop) * pop

    # Death rate over 30
   

    # -------------------------------------------------------------------------
    # Definition of the base case
    # -------------------------------------------------------------------------

    # This is an example,so I am taking the same file as the total emissions 
    # above and I try to find the emissions that are due to GSL and MD 
    # vehicles. When they will be ready I will use Marco's files and sectors.
    # WARNING - the code between **~** has many assumptions ans will be substituted
    # once Marcos files are available.
    # **
    
    # open area for emissions activity calibration(corresponds to Austria)
    nc_file_at = 'netcdf/EMI_RED_ATLAS_NUTS0.nc'
    fh_area_at = Dataset(nc_file_at, mode='r')
    area_at = fh_area_at.variables['AREA'][0][:]
    fh_area_at.close()

    em_den_nh3_7 = emission_dict['NH3'][6]  # [Mg/km2]
    em_den_nmvoc_7 = emission_dict['NMVOC'][6]  # [Mg/km2]
    em_den_nox_7 = emission_dict['NOx'][6]  # [Mg/km2]
    em_den_ppm_7 = emission_dict['PPM'][6]  # [Mg/km2]
    em_den_sox_7 = emission_dict['SOx'][6]  # [Mg/km2]
    
    # Comparisong with the data from TREMOVE for Austria
    at_total_emnmvoc_7 = np.sum(np.multiply(em_den_nmvoc_7, area_cell,
                                           area_at)/100)  # [Mg]
    at_total_emsox_7 = np.sum(np.multiply(em_den_sox_7, area_cell,
                                         area_at)/100)  # [Mg]
    at_total_emnox_7 = np.sum(np.multiply(em_den_nox_7, area_cell,
                                         area_at)/100)  # [Mg]
    at_total_emppm_7 = np.sum(np.multiply(em_den_ppm_7, area_cell,
                                         area_at)/100)  # [Mg]

    # For Austria _ fuel activity from TREMOVE 225.302+14.534 PJ for M7 
    # M7 includes all transport means apart from ships. 
    # For Austria _activity from TREMOVE 169.773 PJ for M7 -PC 
    # M7 -PC includes only as MD and GSL passenger cars
    # The average for the EU28 are used to back calculate activities the 
    # activities as follows:
    # 0.63 is the ratio between the acticity of PC (GSL, MD) 
    # 261.37 is the emission factor for NOx [ton/PJ] 
    # 584.0085736 is the Mpkm/PJ for PC (GSL,MD) 
    
    act_den_7 = em_den_nox_7 * (1/261.37)  # [PJ/km2]
    act_den_7pc = em_den_nox_7 * (1/261.37) * 0.63  # [PJ/km2]
    ser_den_7 = act_den_7 * 584.0085736 # Mpkm/km2      
    ser_den_7pc = act_den_7pc * 584.0085736 # Mpkm/km2 
    
    at_total_act_7pc = np.sum(np.multiply(act_den_7pc, area_cell,
                                            (area_at))/100)  # [PJ]
    # results in 121.61398 PJ instead of 169.7 PJ obtained in TREMOVE, which
    # at the moment can be considered good enough
    # **

    # Low emission zone
    # emission density of pollutant from macrosector 
    # **(modify pollutant and macrosector as they should be variables)
    # em_den_nh3_7 = emission_dict['NH3'][7]  # [Mg/km2]
    # em_den_nmvoc_7 = emission_dict['NMVOC'][7]  # [Mg/km2]
    # em_den_nox_7 = emission_dict['NOx'][7]  # [Mg/km2]
    # em_den_ppm_7 = emission_dict['PPM'][7]  # [Mg/km2]
    # em_den_sox_7 = emission_dict['SOx'][7]  # [Mg/km2]   

    # total emission from  the area of interest of macrosector 7
    # **(do it better, in a loop)
    tot_em_nh3_7 = np.sum(np.multiply(em_den_nh3_7,
                               area_area, area_cell)/100)  # [Mg]
    tot_em_nmvoc_7 = np.sum(np.multiply(em_den_nmvoc_7,
                                 area_area, area_cell)/100)  # [Mg]
    tot_em_nox_7 = np.sum(np.multiply(em_den_nox_7,
                               area_area, area_cell)/100)  # [Mg]
    tot_em_ppm_7 = np.sum(np.multiply(em_den_ppm_7,
                               area_area, area_cell)/100)  # [Mg]
    tot_em_sox_7 = np.sum(np.multiply(em_den_sox_7,
                               area_area, area_cell)/100)  # [Mg]
    # total activities from  the area of interest of macrosector 7
    # **(do it better, in a loop)
    area_tot_act_7 = np.sum(np.multiply(act_den_7,
                               area_area, area_cell)/100)  # [PJ]
    area_tot_ser_7 = np.sum(np.multiply(ser_den_7,
                               area_area, area_cell)/100)  # [Mpkm]
    # -------------------------------------------------------------------------
    # Measures implementation (**will be translated to a def)
    # -------------------------------------------------------------------------
    
    # An example with a low emission zone:                                       
                                      
    # total service activity for PC (GSL,MD) in the area of interest   
    # occ = 1.65  # average accupancy for passenger cars.
    area_tot_ser_7pc= np.sum(np.multiply(ser_den_7pc, area_area, area_cell)/100)
    area_total_act_7pc = np.sum(np.multiply(act_den_7pc, area_cell,
                                            (area_at))/100)  # [PJ]
    # Case 1) Total activity is reduced by 20 % in the area of interest.
    red_ratio = 0.2
    tot_ser_7pc_red = area_tot_ser_7pc * (red_ratio)
    red_on_MS = tot_ser_7pc_red/area_tot_ser_7 * 100
    path_reduction_txt="reductions.txt"
    red = [red_on_MS, red_on_MS, red_on_MS, red_on_MS, red_on_MS]
    
    write_reductions(path_reduction_txt, red)
    # Calculate the reduction absolute and then percentage on the total
    # Write emission reduction file
    # Run model1
    # Substitute the SCE with the result of model 1

    
    # ** to be continued here ** 
    
    
    deltayoll, deltayll_reg, deltayoll_tot = calc_impacts(deaths, pop,
                                                          pop30plus, pyll,
                                                          pm25_cle, area_area)
    
    # -------------------------------------------------------------------------
    # Output of results (todo)
    # -------------------------------------------------------------------------
    
    # create new netcdf file for results
#    nc_file3 = 'netcdf/mollAQI-2010CLE.nc'
#    fh_sce_moll = Dataset(nc_file3, mode='w')
#    fh_sce_moll.createDimension('time', 1)
#    fh_sce_moll.createDimension('y', 448)
#    fh_sce_moll.createDimension('x', 384)
#    time = fh_sce_moll.createVariable('time', 'f8', ('time'))
#    latitude = fh_sce_moll.createVariable('lat', 'f4', ('y', 'x'))
#    longitude = fh_sce_moll.createVariable('lon', 'f4', ('y', 'x'))
#    moll = fh_sce_moll.createVariable('moll', 'f4', ('time', 'y', 'x'))
#    fh_sce_moll.variables['moll'].units = 'months'
#    fh_sce_moll.variables['moll'].long_name = 'Months of lost lives'
#    longitude[:] = lons_sce
#    latitude[:] = lats_sce
#    time[:] = 2.0091231958333332E7
#    moll[:] = mollcle
#    fh_sce_moll.close()
    # create new netcdf file for results
#    nc_file4 = 'netcdf/test.nc'
#    fh_pop30plus = Dataset(nc_file4, mode='w')
#    fh_pop30plus.createDimension('time', 1)
#    fh_pop30plus.createDimension('y', 448)
#    fh_pop30plus.createDimension('x', 384)
#    time = fh_pop30plus.createVariable('time', 'f8', ('time'))
#    latitude = fh_pop30plus.createVariable('lat', 'f4', ('y', 'x'))
#    longitude = fh_pop30plus.createVariable('lon', 'f4', ('y', 'x'))
#    pops = fh_pop30plus.createVariable('area_austria', 'f4', ('time', 'y', 'x'))
#    fh_pop30plus.variables['area_austria'].units = 'test'
#    fh_pop30plus.variables['area_austria'].long_name = 'test'
#    longitude[:] = lons
#    latitude[:] = lats
#    time[:] = 2.0091231958333332E7
#    pops[:] =   area_austria
#    fh_pop30plus.close()
    # Get some parameters for the Stereographic Projection
    #lon_0 = lons.mean()
    #lat_0 = lats.mean()
    
    # setup stereographic basemap.
    # lat_ts is latitude of true scale.
    # lon_0,lat_0 is central point.
    #m = Basemap(width=5000000, height=3500000,
    #            resolution='l', projection='stere',
    #            lat_ts=40, lat_0=lat_0, lon_0=lon_0)
    #m.drawcoastlines()
    #xi, yi = m(lons, lats)
    #
    #cs = m.pcolor(xi, yi, np.squeeze(pm25_cle))
    #cbar = m.colorbar(cs, location='bottom', pad="10%")
    #
    #plt.title("PM2.5")
    #plt.savefig('pm25.png')
    #plt.show()
    #
    #m2 = Basemap(width=5000000, height=3500000,
    #             resolution='l', projection='stere',
    #             lat_ts=40, lat_0=lat_0, lon_0=lon_0)
    #m2.drawcoastlines()
    #xi, yi = m2(lons, lats)
    #cs = m2.pcolor(xi, yi, np.squeeze(pop30plus), vmin=0, vmax=10000)
    #cbar = m2.colorbar(cs, location='bottom', pad="10%")
    #plt.title("population")
    #
    #
    #
    #plt.savefig('moll.png')
    #plt.show()
    pass
