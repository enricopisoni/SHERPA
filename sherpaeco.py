# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 16:26:42 2016

@author: peduzem

INPUT:
    netcdf/sherpaRiat-inerisEMI-pm25AQI-2010CLE.nc
    nc_file_area = 'netcdf/netcdf_with_area/JRC01.nc
    *todo* definition of measure implementation
    red_ratio: reduction of service
    or

OUTPUT:
    file of emission reductions by sector: input/sherpaeco_reduction.txt
    netcdf with the YOLL gained from the measures
    YOLL in the region of interest
    YOLL elsewhere
    as for module1:
    netcdf with concentration changes per pollutant and cell
    delta emission netcdf with emission changes per precursor and cell

... (to be continued)

First trial for benefit assessment based on Enrico's matlab file"""

# imports
# this is necessary because in py27 / is integer division
# whereas in py3.5 / is flat division 
from __future__ import division
                                    
# for importing matlab files
import scipy.io as sio
from PIL import Image
# gdal package to hanlde tiff, because of this package I had to switch to 
# python 2.7 
from osgeo import gdal, osr #conda install -c conda-forge gdal
# from pyGTiff import geotiff
from gdalconst import *
# math
# from numpy import lib, zeros, sum, power, ones
# for using netCDF files
from netCDF4 import Dataset
# for scientific operators
import numpy as np
# for plotting
from mpl_toolkits.basemap import Basemap  #conda install -c conda-forge basemap
import matplotlib.pyplot as plt
from sherpa_auxiliaries import create_emission_reduction_dict, create_emission_dict, create_window, read_progress_log, write_progress_log
from sherpa_globals import path_area_cdf_test, path_model_cdf_test, path_emission_cdf_test, path_result_cdf_test
# Definition of measures and technologies 
#import measures as m

import module1 as shrp
from time import time  # for module1
from os import remove  # for module1

import pandas as pd # conda install pandas
import csv


def write_reductions(path_reduction_txt, red):
    """ 
    Function to write the reduction file, that is the percentage reduction 
    per pollutant per macro sector. 
    At the moment only for MS7 (**will be extended to the other sectors)
    
    Input:
        - path_reduction_txt: path to the reduction file (from sherpa_globals)
        - red: array with the reduction (>0) % per pollutant
    Output:
        - path_retduction_txt.txt: with the % reduction of pollutants of MS7
    
    @author: peduzem
    """ 
    text_file = open(path_reduction_txt, "w")
    text_file.write("POLL	MS1	MS2	MS3	MS4	MS5	MS6	MS7	MS8	MS9	MS10\n")
    text_file.write("NOx	0	0	0	0	0	0	%s	0	0	0\n" % red[0])
    text_file.write("NMVOC	0	0	0	0	0	0	%s	0	0	0\n" % red[1])
    text_file.write("NH3	0	0	0	0	0	0	%s	0	0	0\n" % red[2])
    text_file.write("PPM	0	0	0	0	0	0	%s	0	0	0\n" % red[3])
    text_file.write("SOx	0	0	0	0	0	0	%s	0	0	0\n" % red[4])  
    text_file.close()

def calc_impacts(deaths, pop, pop30plus, pyll, pm25_cle, area_area):

    """
    Methodology (used by Enrico in previous work, slightly modified, check with him**)
    Anenberg, S.C. et al., 2010. Environmental Health Perspectives, 118(9),
    pp.1189–1195. Available at: http://ehp.niehs.nih.gov/0901220
    [Accessed December 13, 2016].
    
    Input:
        - deaths: total deaths of population over 30
        - pop: total population over 30
        - pop30plus: distribution of population over 30
        - pyll: potential years of life loss, PYLL per 100 000 -30+ average
          EU28
        - pm25_cle: current pollutant concentration (units**) 
        ** will modify this to read directly module 1 output
        - area_area: area of interest 
    Output: 
        - deltayoll: delta yoll per grid cell 
        - deltayll_reg, deltayoll_tot
    
    @author: peduzem
    """
    # open SCEnario (2030, MFR) and read PM2.5 values
    # **(This will not be necessary as it will be the result of module1)
    nc_file2 = 'netcdf/sherpaRiat-inerisEMI-pm25AQI-2030MFR.nc'
    fh_sce = Dataset(nc_file2, mode='r')
    pm25_sce = fh_sce.variables['PM25'][:]
    fh_sce.close()

    # Death rate over 30
    drate = deaths/pop

    # Incidence rate (as calculated by Enrico, not used here)
    # ir = pyll / 100000 * pop / deaths

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
    deltamoll_reg = np.sum(deltamoll * area_area / 100)
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

    # read the precursors list (as in model1)  
    rootgrp = Dataset(path_model_cdf_test, 'r')
    precursor_lst = getattr(rootgrp, 'Order_Pollutant').split(', ')
#    precursor_lst = ['NOx', 'NMVOC', 'NH3', 'PPM', 'SOx']
    # create a dictionary with emissions per precursor, macrosector and postion (lat, lon)
    emission_dict = create_emission_dict(path_emission_cdf_test, precursor_lst)

    # open 2010 CLE and read PM2.5 concentration values
    # ** this is not in the input, file with the current concentrations
    # ** may not be necessary in the future
    nc_file = 'netcdf/sherpaRiat-inerisEMI-pm25AQI-2010CLE.nc'
    fh = Dataset(nc_file, mode='r')
    lons = fh.variables['lon'][:]
    lats = fh.variables['lat'][:]
    pm25_cle = fh.variables['PM25'][:]
    pm25_units_cle = fh.variables['PM25'].units
    fh.close()      
       
    # Area of each cell in the domain    
    # ** this is not in the input folder, has to be added
    nc_file_area = 'netcdf/netcdf_with_area/JRC01.nc'
    # ** this is not in the input,
    fh_area_cells = Dataset(nc_file_area, mode='r')
    area_cell = fh_area_cells.variables['surface'][:]
    area_units = fh_area_cells.variables['surface'].units
    fh_area_cells.close() 
    
    
    # open area of interest
    # This is just to check data wiht austria, comment to consider the area of 
    # interest (see below)
    nc_file_at = 'input/EMI_RED_ATLAS_NUTS0.nc'
    fh_area_at = Dataset(nc_file_at, mode='r')
    area_area = fh_area_at.variables['AREA'][0][:] 
    fh_area_at.close()
    
    # ** Uncomment this to make it work with the area of interest. 
    #    nc_file_reg = path_area_cdf_test 
    #    fh_area = Dataset(nc_file_reg, mode='r')
    #    area_area = fh_area.variables['AREA'][:]
    #    fh_area.close() 
    
    ds   = gdal.Open('CO2_emiss/7km_eur_TRA_RD_LD4C_GSL_Mall.tif.tif')
    arr    = ds.ReadAsArray()
    [cols,rows] = arr.shape
    (Xarr, deltaX, rotation, Yarr, rotation, deltaY) = ds.GetGeoTransform()
    #CO2_TRA_RD_LD4C_GSL_Mall = np.array(ds.GetRasterBand(1).ReadAsArray())
    ds = None   

    # Define longitude and latitude for the data
    longl = []
    for i in range(0, rows): 
        longl.append(Xarr + i*deltaX + deltaX*0.5)
    lonsmtif=np.array(longl)
    latl = []
    for i in range(0, cols): 
        latl.append(Yarr + i*deltaY + deltaY*0.5)
    latsmtif=np.array(latl)
    X, Y = np.meshgrid(lonsmtif , latsmtif)

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

    # -------------------------------------------------------------------------
    # Rading data from Marco's inventory files
    # -------------------------------------------------------------------------
   
    # Input parameters
    # ** should be a function
    m_sector = 7 # macro sector, SNAP # 7: road transport 
    #pollutant = 'NH3'  
    pollutant_lst = ['NH3','NOx','VOC', 'SO2', 'PM10']
    # sector = 'TRA_RD_LD4C' #sub sector
    sector_lst = ['TRA_RD_LD4C','TRA_RD_HDB','TRA_RD_LD2','TRA_RD_LD4T','TRA_RD_HDT' ]
    # net = network (** will be all, urban, rural, motorway etc.)
    net = 'all' # network
    fuel_lst = ['GSL', 'MD', 'GAS', 'LPG', 'TYRE']

    # Take marco's .tif files and build corresponding arrays with the information on the 
    # emissions
    
    # initialiaze arrays
    emi = {} # emission inventory (for Sherpa)
    em_inv ={} # emission inventory (Marco's grid)
    act = {} # activity [PJ]
    ser = {} # service [Mvkm]
    ef_inv ={} # emission factor inventory (Marco)
    em_new ={} # emission after measure implementation
    red = {} # reduction at the macrosector level 
    
    # calculate the activity level starting from the CO2 emission inventory
    # Activity in PJ for each sector-fuel combination in Marcos inventory
    # calculating it from the CO2 because directly proportional to the fuel 
    # consumption (only ttw emissions), same emission factor for all countries
    # apart from france and the united kingdom (I don't know why)
    # ** the different emission factors for F and UK should be taken into account  
    # initialiaze arrays
    emi['CO2'] = {}
    for sector in sector_lst: 
        # initialiaze arrays
        emi['CO2'][sector] = {}
        act[sector] = {} 
        for fuel in fuel_lst:
            act[sector][fuel] = {}
            act[sector][fuel][net] = {}
            emi['CO2'][sector][fuel] = {}
            emi['CO2'][sector][fuel][net] = {} 
            # open CO2 emission inventory
            ds = None
            try:
                ds = gdal.Open('CO2_emiss/7km_eur_{}_{}_M{}.tif.tif'.format(sector, fuel, net))
                # write data in array                
                emco2  = np.array(ds.GetRasterBand(1).ReadAsArray())
                ds = None
                # re-arrange emission inventory file for Sherpa's grid 
                # (tried copying from Enrico's matlab):
                # initialize array
                Amineco2 = emco2
                for i in range(1,382):
                    ind1=2*(i-1)  # included
                    ind2=1+2*(i-1)+1 # excluded
                    Amineco2[:,i-1]=(np.sum(emco2[:,ind1:ind2],axis=1))
                Amineco2[:,382:384]=0
                # Cancelling the extra columns and extra rows 
                # (there has to be a better way to do this) 
                for deli in range (0,144):
                    Amineco2 = np.delete(Amineco2, (0), axis=0) # delete first 144 rows
                for delj in range (0,398): # delete last 398 columns
                    Amineco2 = np.delete(Amineco2, (383), axis=1)
                Amineco2_T = Amineco2[np.newaxis]
                Afinalco2=np.fliplr(Amineco2_T)
        
                # CO2 emissions (Sherpa's grid)
                emi['CO2'][sector][fuel][net] =  Afinalco2
                # Activity emissions (Sherpa's grid)
                act[sector][fuel][net] = np.sum(emi['CO2'][sector][fuel][net]*(1/66219.992) * area_area / 100)*1000000  # [PJ]
            except(RuntimeError, AttributeError):
                    pass
                
    emi['PM10'] = {}
    emftyre = {}

    emftyre['TRA_RD_LD4C']=0.0089238 # kton/Gvkm
    emftyre['TRA_RD_HDB']=0.0152918237 # kton/Gvkm
    emftyre['TRA_RD_LD2']=0.0038364 # kton/Gvkm
    emftyre['TRA_RD_LD4T']=0.0140946 # kton/Gvkm
    emftyre['TRA_RD_HDT']=0.0396 # kton/Gvkm
    emftyre['M4']=0.0038364 # kton/Gvkm

    for sector in sector_lst: 
        # initialiaze arrays
        emi['PM10'][sector] = {}
        ser[sector] = {} 
        tyre = 'TYRE'
        ser[sector] = {}
        ser[sector][net] = {}
        emi['PM10'][sector][tyre] = {}
        emi['PM10'][sector][tyre][net] = {} 
        # open PM10 emission inventory
        ds = None
        try:
            ds = gdal.Open('PM10_emiss/7km_eur_{}_{}{}.tif.tif'.format(sector, fuel, net))
#           # write data in array                
            emPPM  = np.array(ds.GetRasterBand(1).ReadAsArray())
            ds = None
            # re-arrange emission inventory file for Sherpa's grid 
            # (tried copying from Enrico's matlab):
            # initialize array
            Amineppm = emPPM
            for i in range(1,382):
                ind1=2*(i-1)  # included
                ind2=1+2*(i-1)+1 # excluded
                Amineppm[:,i-1]=(np.sum(emPPM[:,ind1:ind2],axis=1))
            Amineppm[:,382:384]=0
            # Cancelling the extra columns and extra rows 
            # (there has to be a better way to do this) 
            for deli in range (0,144):
                Amineppm = np.delete(Amineppm, (0), axis=0) # delete first 144 rows
            for delj in range (0,398): # delete last 398 columns
                Amineppm = np.delete(Amineppm, (383), axis=1)
            Amineppm_T = Amineppm[np.newaxis]
            Amineppm=np.fliplr(Amineppm_T)
#    
            # PPM emissions (Sherpa's grid)
            emi['PM10'][sector][tyre][net] =  Amineppm
            # Activity emissions (Sherpa's grid)
            ser[sector][net] = np.sum(emi['PM10'][sector][tyre][net]*(1/emftyre[sector]) * area_area / 100) *1000 # [Mvkm]
        except(RuntimeError, AttributeError):
                pass
            
    for pollutant in pollutant_lst: 
        emi[pollutant] = {}
        em_inv[pollutant]={}
        ef_inv[pollutant]={}        
    
        for sector in sector_lst: 
            # initialize arrays
            em_inv[pollutant][sector] = {} 
            ef_inv[pollutant][sector] = {} 
            emi[pollutant][sector] = {}
    
            for fuel in fuel_lst:
                # initialize arrays
                emi[pollutant][sector][fuel] = {}
                emi[pollutant][sector][fuel][net] = {} 
                em_inv[pollutant][sector][fuel] = {}
                em_inv[pollutant][sector][fuel][net] = {}
                ef_inv[pollutant][sector][fuel] = {}
                ef_inv[pollutant][sector][fuel][net] = {}
                # open emission inventory
                gdal.UseExceptions()
                ds = None
                try:
                    ds   = gdal.Open('{}_emiss/7km_eur_{}_{}{}.tif.tif'.format(pollutant, sector, fuel, net))
                    em  = np.array(ds.GetRasterBand(1).ReadAsArray())
                    ds = None
                    # re-arrange emission inventory file for Sherpa's grid 
                    # (tried copying from Enrico's matlab):
                    # initialize array
                    Amine = em # Amine : matrix
                    for i in range(1,382):
                        ind1=2*(i-1)  # included
                        ind2=1+2*(i-1)+1 # excluded
                        Amine[:,i-1]=(np.sum(em[:,ind1:ind2],axis=1))
                    Amine[:,382:384]=0
                    # Cancelling the extra columns and extra rows 
                    # (there has to be a better way to do this)            
                    for deli in range (0,144):
                        Amine = np.delete(Amine, (0), axis=0) # delete first 144 rows
                    for delj in range (0,398): # delete last 398 columns
                        Amine = np.delete(Amine, (383), axis=1)
                    Amine_T = Amine[np.newaxis]
                    Afinal=np.fliplr(Amine_T)
        
                    emi[pollutant][sector][fuel][net] = Afinal
                       
                    area_tot_em_dict = {}
                    for precursor in precursor_lst:
                        for snap in range(1, 11):
                            area_tot_em_dict[precursor, snap - 1] = np.sum(emission_dict[precursor][snap - 1]  * area_cell * area_area / 100)
                            # [Mg]         
                
                    # Emissions of the sector-fuel combination in Marcos inventory
                    try: 
                        em_inv[pollutant][sector][fuel][net] = np.sum(emi[pollutant][sector][fuel][net]* area_area / 100) * 1000  # ton    
                        if fuel is 'TYRE':
                            ef_inv[pollutant][sector][fuel][net] = em_inv[pollutant][sector][fuel][net]/ser[sector][net]
                        else:
                            ef_inv[pollutant][sector][fuel][net] = em_inv[pollutant][sector][fuel][net]/act[sector][fuel][net]
                    except RuntimeError, err:
                        pass
                except(RuntimeError, AttributeError):
                    pass
                # if the array is empty write zero (there's probably a better way to do this so that I can do the sum below)
                if not em_inv[pollutant][sector][fuel][net]:
                    em_inv[pollutant][sector][fuel][net]=0.0
    # write data in array 
    for pollutant in pollutant_lst: 
        em_inv[pollutant] = np.sum(np.sum(em_inv[pollutant][sector][fuel][net] for fuel in fuel_lst ) for sector in sector_lst)
#    # -------------------------------------------------------------------------
#    # Measures implementation (**will be translated to a def)
#    # -------------------------------------------------------------------------
#   
    # New emission factors after the implementation of measures:
    # first option - reduction of the emission factors for each sector/activity
    # read csv file with the emission factors
    df = pd.read_csv('input/ef_red_sherpaeco.csv',  index_col=[0,1], names=['POLL','ACT'].extend(sector_lst), skipinitialspace=True)
    # second option (to be implemented**)- emission factors for the best available technology
    
    em_bc ={} # emissions for the base case
    for precursor in precursor_lst:
        em_bc[precursor,m_sector-1] = np.sum(emission_dict[precursor][m_sector-1] * area_cell * area_area / 100)
    
    for pollutant in pollutant_lst: 
        em_new[pollutant]={}
        for sector in sector_lst:    
            em_new[pollutant][sector] = {}
            for fuel in fuel_lst:
                em_new[pollutant][sector][fuel] = {}
                em_new[pollutant][sector][fuel][net] = {}
                if ef_inv[pollutant][sector][fuel][net]:
                    ef = ef_inv[pollutant][sector][fuel][net] *(1 - df.loc[pollutant,fuel][sector])
                    if fuel is 'TYRE':
                        em_new[pollutant][sector][fuel][net] = np.sum(ser[sector][net] * area_area / 100) * ef
                    else:
                        em_new[pollutant][sector][fuel][net] = np.sum(act[sector][fuel][net] * area_area / 100) * ef
                     # if the array is empty write zero (there's probably a better way to do this so that I can do the sum below)
                if not em_new[pollutant][sector][fuel][net]:
                    em_new[pollutant][sector][fuel][net]=0.0

    for pollutant in pollutant_lst: 
        em_new[pollutant] = np.sum(np.sum(em_new[pollutant][sector][fuel][net] for fuel in fuel_lst) for sector in sector_lst)
#        
#        
    for precursor in precursor_lst: 
        red[precursor]={}
        for pollutant in pollutant_lst:
            if pollutant == precursor:
                red[precursor] = (em_inv[pollutant]-em_new[pollutant])/em_bc[precursor,m_sector-1]*100
        if precursor == 'NMVOC':
                red[precursor] = (em_inv['VOC']-em_new['VOC'])/em_bc[precursor,m_sector-1]*100
        if precursor == 'SOx': 
                red[precursor] = (em_inv['SO2']-em_new['SO2'])/em_bc[precursor,m_sector-1]*100
        if precursor == 'PPM':
                red[precursor] = (em_inv['PM10']-em_new['PM10'])/em_bc[precursor,m_sector-1]*100

    reductions = {}
    reductions[m_sector-1]={}
#    #                         NOx	     NMVOC      NH3        PPM        SOx
    reductions[m_sector-1] = [red['NOx'], red['NMVOC'], red['NH3'] , red['PPM'], red['SOx']]
    path_reduction_txt='input/sherpaeco_reduction.txt'
    write_reductions(path_reduction_txt, reductions[m_sector-1])

  
    # -------------------------------------------------------------------------
    # Running module1 
    # -------------------------------------------------------------------------

#    # if it doesn't exist start=0 and dividsor=1
#    progresslog = 'input/progress.log'
#    start = time()
#    output = 'output/'
#
#    # run module 1 with progress log
#    proglog_filename = path_result_cdf_test + 'proglog'
#    write_progress_log(proglog_filename, 25, 2)
#    start = time()
#    shrp.module1(path_emission_cdf_test, path_area_cdf_test,
#                 path_reduction_txt, path_model_cdf_test, output)
#    stop = time()
#    print('Module 1 run time: %s sec.' % (stop-start))
#    remove(proglog_filename)
#    
    # -------------------------------------------------------------------------
    # (Cost) Benefit Analysis 
    # -------------------------------------------------------------------------
#    deltayoll, deltayll_reg, deltayoll_tot = calc_impacts(deaths, pop,
#                                                          pop30plus, pyll,
#                                                          pm25_cle, area_area)
    
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
#    lon_0 = lons.mean()
#    lat_0 = lats.mean()
##    
##    # setup stereographic basemap.
##    # lat_ts is latitude of true scale.
##    # lon_0,lat_0 is central point.
#    check = Aco2new-Afinal
#    m = Basemap(width=5000000, height=4500000,
#                resolution='l', projection='stere',
#                lat_ts=40, lat_0=lat_0, lon_0=lon_0)
#    m.drawcoastlines()
#    xi, yi = m(lons, lats)
#    
#    cs = m.pcolor(xi, yi, np.squeeze(Aco2new), vmin=0, vmax=0.01)
#    cbar = m.colorbar(cs, location='bottom', pad="10%")
#    
#    plt.title("PM2.5")
#    plt.savefig('pm25.png')
#    plt.show()
#    
##    lon_0 = lonsmtif.mean()
##    lat_0 = latsmtif.mean()
#    m2 = Basemap(width=5000000, height=4500000,
#                 resolution='l', projection='stere',
#                 lat_ts=40, lat_0=lat_0, lon_0=lon_0)
#    m2.drawcoastlines()
##    xi, yi = m2(X, Y)
#    xi, yi = m(lons, lats)
#    cs = m2.pcolor(xi, yi, np.squeeze(check), vmin=0, vmax=0.00000000001)
#    cbar = m2.colorbar(cs, location='bottom', pad="10%")
#    plt.title("population")
#    
    #
    #
    #plt.savefig('moll.png')
#    plt.show()
#    pass
#           A = sio.loadmat('co2.mat')
#            Aco2 = A['Amine']
#            Aco2_T = Aco2[np.newaxis]
#            Aco2new=np.fliplr(Aco2_T)  
  # Comparisong with the data from TREMOVE for Austria
           
            # Austria, TREMOVE, M7 (all road transport) vs SHERPA vs Marco
            # Activity 221.679  vs ... vs 316 [PJ]
            # NMVOC 12210 [Mg] vs 14457 
            # SOX 104 [Mg] vs 144
            # NOX 61488 [Mg] vs 107133
            # PPM 2112+707 vs 5174
        
            # For Austria, GSL PC from Marco
            # activity 71.33405959 PJ for GLS PC (Marco) 
            # 0.225 is the ratio between the acticity of PC (GSL) and M7 (Marco) 316 PJ
            # 87.346743 is the emission factor for nmvoc [ton/PJ] GLS PC (Marco)  
            # 101.161882 is the emission factor for NOx [ton/PJ] GLS PC (Marco) 
            # 0.584308195 for PM2.5, 0.584308195 for PM10 [ton/PJ] GLS PC (Marco) 
            # 0.434782609 for SO2/SOx [ton/PJ] GLS PC (Marco) 
            # 24.55769089 for NH3 [ton/PJ]
            # 66219.992 for CO2 (ttw I think) [ton/PJ]
            

    #    nc_marco = '7km_eur_TRA_RD_LD4C_GSL_Mall.nc'
    #    fh_marco = Dataset(nc_marco, mode='r')
    #    CO2_TRA_RD_LD4C_GSL_Mall = fh_marco.variables['7km_eur_TRA_RD_LD4C_GSL_Mall.tif.tif'][:]
    #    lonsm = fh_marco.variables['lon'][:]
    #    latsm = fh_marco.variables['lat'][:]
    #    X, Y = np.meshgrid(lonsm , latsm)
    #    fh_marco.close()
        
    # Open Marco's tif files and write the information in an array   
#    ds   = gdal.Open('CO2_emiss/7km_eur_TRA_RD_LD4C_GSL_Mall.tif.tif')
#    arr    = ds.ReadAsArray()
#    [cols,rows] = arr.shape
#    (Xarr, deltaX, rotation, Yarr, rotation, deltaY) = ds.GetGeoTransform()
#    CO2_TRA_RD_LD4C_GSL_Mall = np.array(ds.GetRasterBand(1).ReadAsArray())
#    ds = None    
            
            
#    
#    # An example with a low emission zone:                                       
#                                      
##    # total service activity for PC (GSL,MD) in the area of interest   
##    area_tot_act_7pc = np.sum(act_den_7pc * area_cell *area_area / 100)  # [PJ]
##    area_tot_ser_7pc= np.sum(ser_den_7pc * area_area * area_cell / 100)
#    ser_tot_sce = {}
#
#    # Case 1) Total activity is reduced by 20 % in the area of interest.
#    red_ratio = 0.2  
#    ser_tot_sce[m_sector-1,sector] = ser_tot[m_sector-1,sector] * (1-red_ratio)
#    red_on_MS = (ser_tot[m_sector-1, sector]-ser_tot_sce[m_sector-1,sector])/ ser_tot[m_sector-1] * 100
#    
#    red={}
#    #                  NOx	     NMVOC      NH3        PPM        SOx
#    red[m_sector-1] = [red_on_MS, red_on_MS, 0  , red_on_MS, red_on_MS]
#    path_reduction_txt='input/sherpaeco_reduction.txt'
#    write_reductions(path_reduction_txt, red[m_sector-1])
#
#    # Case 2) Substitute mobility service and or increase p/v
#    occ = 1.65  # average accupancy for passenger cars.
#    # Data from TREMOVE show slighlty different occ
#    # Here the average value is assumed
#    newocc = 1.65  # p/v
#    m.av_gsl_pc_eu5.occ = newocc  # p/v
#    m.av_dsl_pc_eu6.occ = newocc  # p/v
#    dis_den[m_sector-1,sector] = ser_den[m_sector-1,sector] * newocc # Mpkm/km2 
#    dis_tot[m_sector-1,sector] = np.sum(dis_den[m_sector-1,sector]* area_cell * area_area /100)  # [Mvkm])
#    tot_sub_ratio = 1  # fraction of mobility substituted in the area
#    gsl_sub_ratio = tot_sub_ratio * 0.7  # fraction of tot_sub_ratio substituted with gsl
#    dsl_sub_ratio = tot_sub_ratio - gsl_sub_ratio 
#    gsl_ser = ser_tot[m_sector-1,sector] * tot_sub_ratio * gsl_sub_ratio  # Mpkm
#    dsl_ser = ser_tot[m_sector-1,sector] * tot_sub_ratio * dsl_sub_ratio  # Mpkm
#    gsl_dist = dis_tot[m_sector-1,sector] * tot_sub_ratio * gsl_sub_ratio  # Mvkm
#    dsl_dist = dis_tot[m_sector-1,sector] * tot_sub_ratio * dsl_sub_ratio  # Mvkm
#    gsl_act = gsl_ser / m.av_gsl_pc_eu5.service_eff()  # PJ
#    dsl_act = dsl_ser / m.av_dsl_pc_eu6.service_eff()  # PJ
#    
#    gls_em = {}
#    for k in m.av_gsl_pc_eu5.emf:
#        if k is 'PPM_nex':
#            gls_em[k] =  m.av_gsl_pc_eu5.emf[k]*gsl_dist
#        else: 
#            gls_em[k] = m.av_gsl_pc_eu5.emf[k]*gsl_act
#    
#    dsl_em = {}
#    for k in m.av_dsl_pc_eu6.emf:
#        if k is 'PPM_nex':
#            dsl_em[k] =  m.av_dsl_pc_eu6.emf[k]*gsl_dist
#        else: 
#            dsl_em[k] = m.av_dsl_pc_eu6.emf[k]*gsl_act
#    
#    new_em={}
#    for k in m.av_gsl_pc_eu5.emf:    
#        #for precursor in precursor_lst:
#        for precursor in precursor_lst:
#            if k in precursor_lst:
#                new_em[k]=gls_em[k]+dsl_em[k]
#            elif k is 'CO2_wtw':
#                new_em[k]=gls_em['CO2_wtw']+dsl_em['CO2_wtw']
#            elif k is 'PPM_nex':
#                new_em['PPM']=gls_em['PPM_nex']+gls_em['PPM_ex']+dsl_em['PPM_nex']+dsl_em['PPM_ex']
#            
#    finalpre_lst = ['NMVOC','NOx','PPM','SOx']
#
#    red = {}
#
#    ## start again from here, data needs to be more detailed to make sense. 
#    for finalprecursor in finalpre_lst:    
#        red[finalprecursor]= ((em_tot[finalprecursor,m_sector-1]*0.63- new_em[finalprecursor])/em_tot[finalprecursor,m_sector-1])*100
#    redPPM =   ((tot_em_ppm_7/area_tot_act_7)*(gsl_act+dsl_act)-(gsl_em_ppm+dsl_em_ppm)) / tot_em_ppm_7 * 100

#    # Case 3) Increase PC