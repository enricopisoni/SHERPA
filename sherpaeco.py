# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 16:26:42 2016
@author: peduzem

Module to calculate the impact of measures

INPUT:
     *to be completed/updated*
    Need to add an input directory with the following files:
        input/sherpaRiat-inerisEMI-pm25AQI-2010CLE.nc : using it to get the values of lat and lon, will not be necessary in the future
        input/population.mat: population file by Marco
        # input/20151116_SR_no2_pm10_pm25/BC_conc_PM25_Y.nc' : base case concentration 
        input/EMI_RED_ATLAS_NUTS0.nc: or NUTS1 to select the area of interest
        input/JRC01.nc : area of each cell 
        input/ef_reduction_sherpaeco : csv file with the reduction per sector/fuel (see example)
        

OUTPUT:
    *to be completed/updated*
    input/sherpaeco_reduction.txt: file of emission reductions per macro sector 
    Delta Years of Life Lost and Delta Mortality in the area of interest and elsewhere
    deltayll_reg, deltayll_tot, delta_mort_tot, delta_mort_reg
    as for module1:
    netcdf with concentration changes per pollutant and cell
    delta emission netcdf with emission changes per precursor and cell
... (to be continued)

"""

# imports
                              
# for importing matlab files
import scipy.io as sio
from PIL import Image
from osgeo import gdal, osr #conda install -c conda-forge gdal
from gdalconst import *
from netCDF4 import Dataset # for using netCDF files
import numpy as np # for scientific operators
# for plotting
from mpl_toolkits.basemap import Basemap  #conda install -c conda-forge basemap
import matplotlib.pyplot as plt
from sherpa_auxiliaries import create_emission_reduction_dict, create_emission_dict, create_window, read_progress_log, write_progress_log
from sherpa_globals import path_area_cdf_test, path_model_cdf_test, path_emission_cdf_test, path_result_cdf_test
# Definition of measures and technologies 
#import measures as m
import module1 as shrp
import module1ema as shrp1
from time import time  # for module1
from os import remove  # for module1
import pandas as pd # conda install pandas
import csv
# to save results direclty as python objects
import pickle
import collections

class data( object ):
    def __init__( self, pol, sec, aty, net, value ):
        self.pollutant= pol
        self.sector= sec
        self.aty= aty
        self.net= net
        self.value= value

def groupBy( facts, name ):
    total= collections.defaultdict( int )
    for f in facts:
        key= getattr( f, name )
        total[key] += f.count
# -----------------------------------------------------------------------------
def save_obj(obj, name ):
    with open('workdir/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
# -----------------------------------------------------------------------------
def load_obj(name ):
    with open('workdir/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)   
# -----------------------------------------------------------------------------
def write_reductions(path_reduction_txt, red):
    """ 
    Function to write the reduction file, that is the percentage reduction 
    per pollutant per macrosector. 
    At the moment only for MS7 (**will be extended to the other sectors)
    
    Input:
        - path_reduction_txt: path to the reduction file (** can be from sherpa_globals)
        - red: array with the reduction (>0) % per pollutant (** the order is important)
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
# -----------------------------------------------------------------------------
def calc_impacts(deltaconc, area_area):

    """
    Health impacts of PM2.5
    Methodology detailed in the REVIHAAP and HRAPIE studies led by WHO-Europe, as described in:
    Holland, M., 2014. Cost-benefit Analysis of Final Policy Scenarios for the EU Clean Air Package Version. Version 2
    
    Concentration-Response-Function taken from the software from AirQ+ (WHO)
    
    Years of life loss as a function of mortality are calculated according to: 
    Anenberg, S.C. et al., 2010. Environmental Health Perspectives, 118(9),
    pp.1189–1195. Available at: http://ehp.niehs.nih.gov/0901220
    [Accessed December 13, 2016].
    
    Years of life loss as calculated as a function of life expectancy according to  
    Holland, M., 2014. Cost-benefit Analysis of Final Policy Scenarios for the EU Clean Air Package Version. Version 2
    
    Data for baseline  population
    ICD codes: ICD-10: A00-B99,C00-D48,D50-D89,E00-E88,F01-F99,G00-G98,
    H00-H59,H60-H93,I00-I99,J00-J98,K00-K92,L00-L98,M00-M99,N00-N98,
    O00-O99,P00-P96,Q00-Q99,R00-R99
    Age: '30 - 85 +'									
    Sex: Both									
    http://data.euro.who.int/dmdb/
    [Accessed December 13, 2016].
     
    NB: Hereafter a positive delta means a reduction! 
    
    
    Input:
        - area_area: area of interest (nc file with 100 in the cells of the area of interest)
        - deltaconc: path to the file  with the delta concentrations (output of module1)
    Output: 
        - deltayll_reg, deltayll_tot, delta_mort_pp, delta_yll_pp
    
    @author: peduzem
    """  
    
    # Load population file from Marco (treated in "regridPop" routine from Enrico).
    # The matlab file contains a structure called Anew,
    # I carried out the same operations as in Enrico's matlab file.
    A = sio.loadmat('input/population.mat')
    Anew = A['Anew']
    Anew_T = Anew[np.newaxis]
    popall=np.fliplr(Anew_T)
    
    # total pop according to Marco (all age groups)
    sumpop=np.sum(popall)

    # TOTAL population above 30 data.euro.who.int/dmdb/
    pop = 331923577
    # TOTAL DEATHS  above 30
    deaths = 4639244
    # Potential Years of life loss, PYLL 30+ total in EU28
    ylltot = 14038453.71
    
    # Distribution of the population over 30
    # HP: assuming it has the same distribution as all the population (as provided by Marco)
    pop30plus = (popall/sumpop) * pop 

    # open delta concentration for PM25 - result of module1
    fh_deltapm25 = Dataset(deltaconc, mode='r')
    pm25_delta = fh_deltapm25.variables['delta_concentration'][:]
    fh_deltapm25.close()

    # Death rate over 30
    drate = deaths/pop

    # crf derived by RR in Anenberg, S.C. et al., 2010. 
    # crf = 0.006
 
    # crf Taken from AirQ+ (WHO)
#    cutoff = 10 # microg/m3 # Taken from AirQ+ (WHO)
    beta = 0.006015392281974714 # Taken from AirQ+ (WHO)
#    baseconfile = 'input/20151116_SR_no2_pm10_pm25/BC_conc_PM25_Y.nc'
#    fh_basecon = Dataset(baseconfile , mode='r')
#    baseconc = fh_basecon.variables['conc'][:]
#    fh_basecon.close()
#    

    # Delta mortality: uncomment the most suitable one
    
    # linear approximation
    # delta_mort = pm25_delta*crf*pop30plus*drate
#    # exponential approximation
#    delta_mort = (1-(np.exp(-beta*pm25_delta)))*pop30plus*drate  
    # formula with threshold, don't think it is necessary and I am not sure it is right**
    # create and array of 1s where the condition is met, 0s elsewhere. 
#    mask1 = (baseconc > cutoff).astype(int)
#    mask2 = (baseconc-pm25_delta > cutoff).astype(int)  
    delta_mort = (1-(np.exp(-beta*pm25_delta)))*pop30plus*drate
    delta_mort_tot = np.sum(delta_mort)
    delta_mort_reg = np.sum(delta_mort*area_area/100)
    
    # Delta Years of life loss (yll) according to Anenberg 2010
    # delta_yll = delta_mort * ylltot / deaths
    # ** I tried validating these values but I obtain very different results 
    # from the ones in Holland, M., 2014. 
    # 79.9 is the life expectancy, should be by country, this is the average from EUROSTAT 
    # http://ec.europa.eu/eurostat/statistics-explained/images/9/93/Life_expectancy_at_birth%2C_1980%E2%80%932014_%28years%29_YB16.png
    lyg = np.exp(8.161-(0.04478*79.9)) # life years gained per 100000 ppl for a unit concentration change
    delta_yll= pm25_delta * lyg / 100000* pop30plus
    # Calculate delta yll in the selected region
    deltayll_reg = np.sum(delta_yll * area_area / 100)
    
    # Calculate delta moll and yll in total
    deltayll_tot = np.sum(delta_yll)
    
    # Per person 
    delta_mort_pp=(1-(np.exp(-beta*pm25_delta)))*drate 
    delta_yll_pp = delta_mort_pp * ylltot / deaths  
    return deltayll_reg, deltayll_tot, delta_mort_tot, delta_mort_reg, delta_yll_pp
# -----------------------------------------------------------------------------  


def tiftosherpagrid(pollutant, variable, sector, aty, net, dct):   
    gdal.UseExceptions()
    ds = None
    try: 
        if pollutant is not 'CO2':  
            ds   = gdal.Open('{}_{}/7km_eur_{}_{}{}.tif.tif'.format(pollutant, variable, sector, aty, net)) 
        else:
            ds = gdal.Open('{}_{}/7km_eur_{}_{}_M{}.tif.tif'.format(pollutant, variable, sector, aty, net))
        em  = np.array(ds.GetRasterBand(1).ReadAsArray())
        ds = None
#                    # re-arrange emission inventory file for Sherpa's grid 
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
        dct=Afinal
        return dct
    except(RuntimeError, AttributeError):
        pass
    
def readinventories(m_sector, pollutant_lst, sector_lst, net, acttype_lst, nonfuels_lst, name_emi, name_gainsCO2, name_act_fossil,name_act_all,name_ser):
    """
    
    Calculate the energy use activity level starting from the CO2 emission inventory
    Activity in PJ for each sector-fuel combination in Marcos inventory
    I am calculating it from the CO2 because emissions are directly
    proportional to the fuel consumption (only ttw emissions), 
    this means that there is the same emission factor for all countries
    apart from france and the united kingdom (I don't know why)
    ** the different emission factors for F and UK should be taken into account  
         
    Calculate the service activity level starting from the PM10-TYRE emission inventory
    Activity in Gvkm for each sector in Marcos inventory
    I am calculating it from the PM10-TYRE because emission factors
    are the same for for all countries in the EU (though they are different for 
    other countries and I don't know why)
    
    Calculate the reference emission factors and emission ivnentories, 
    for each pollutant, sector, fuel combination from Marco's inventories    
    INPUT:      
        m_sector = 7 # macro sector, SNAP # 7: road transport 
        pollutant_lst = ['NH3','NOx','VOC', 'SO2', 'PM10', 'CO2'] # List of pollutants 
        sector_lst = ['TRA_RD_LD4C','TRA_RD_HDB','TRA_RD_LD2','TRA_RD_LD4T','TRA_RD_HDT','TRA_RD_M4' ]# List of subsectors of the m_sector
        net_lst = ['rur', 'urb', 'mot']  
        acttype_lst = ['GSL', 'MD', 'GAS', 'LPG','TYRE', 'ABRASION','BRAKE'] # activity types
        nonfuels_lst = ['TYRE', 'ABRASION','BRAKE'] # activities that are not fuels
        name_em_inv= 'em_inv' name of python binary file to save results 
        name_ser= 'ser' name of python binary file to save results
        name_act= 'act' name of python binary file to save results

    
    @author: peduzem
    """

    emfco2ipcc={}
    emfco2ipcc['GSL']= 69200 # ton/PJ  , IPCC2006 69200 
    emfco2ipcc['MD']= 74000 # ton/PJ, IPCC2006 69200 
    emfco2ipcc['LPG']= 63000 # ton/PJ, IPCC2006 69200 
    emfco2ipcc['GAS']= 56100 # ton/PJ, IPCC2006 69200 
    
    
    pollutant = 'CO2'
    variable = 'act_lev'
    act_fossil={}
    for sector in sector_lst: 
        act_fossil[sector]={}
        for aty in [aty for aty in acttype_lst if aty not in nonfuels_lst]:  
            act_fossil[sector][aty]={}
            for net in net_lst:
                value=tiftosherpagrid(pollutant, variable, sector, aty, net, act_fossil)
                if value is not None:
                    act_fossil[sector][aty][net]={}
                    act_fossil[sector][aty][net]=value
            if not act_fossil[sector][aty]: 
                act_fossil[sector].pop(aty, None)

    
    pollutant = 'PM25'
    variable = 'act_lev'
    act_all={}
    for sector in sector_lst: 
        act_all[sector]={}
        for aty in [aty for aty in acttype_lst if aty not in nonfuels_lst]: 
            act_all[sector][aty]={}
            for net in net_lst:
                value=tiftosherpagrid(pollutant, variable, sector, aty, net, act_all)
                if value is not None:
                    act_all[sector][aty][net]={}
                    act_all[sector][aty][net]=value
            if not act_all[sector][aty]: 
                act_all[sector].pop(aty, None)
                
    
    pollutant = 'PM25'
    variable = 'act_lev'
    ser={}
    for sector in sector_lst:
        ser[sector]={} 
        for aty in nonfuels_lst:
            ser[sector][aty]={}
            for net in net_lst:  
                value=tiftosherpagrid(pollutant, variable, sector, aty, net, act_all)
                if value is not None:
                    ser[sector][aty][net]={}
                    ser[sector][aty][net]=value
            if not ser[sector][aty]: 
                ser[sector].pop(aty, None)

    variable= 'emiss'
    emi={}  
    gainsCO2={}
    for pollutant in pollutant_lst: 
        emi[pollutant]={}
        for sector in sector_lst:
            emi[pollutant][sector]={}
            for aty in acttype_lst:
                emi[pollutant][sector][aty]={}
                for net in net_lst:
                    if pollutant is 'CO2' and sector in act_fossil and aty in act_fossil[sector]:
                        emi[pollutant][sector][aty][net]={};
                        emi[pollutant][sector][aty][net]=act_fossil[sector][aty][net]*emfco2ipcc[aty]
                        value=tiftosherpagrid(pollutant, variable, sector, aty, net, emi)
                        if value is not None:
                             gainsCO2[sector]={};gainsCO2[sector][aty]={};gainsCO2[sector][aty][net]={};
                             gainsCO2[sector][aty][net]=value
                    elif pollutant is not 'CO2': 
                        value=tiftosherpagrid(pollutant, variable, sector, aty, net, emi)
                        if value is not None:
                            emi[pollutant][sector][aty][net]={}
                            emi[pollutant][sector][aty][net]=value
                if not emi[pollutant][sector][aty]: 
                    emi[pollutant][sector].pop(aty, None)
                         
    save_obj((emi), name_emi) 
    save_obj((gainsCO2), name_gainsCO2)
    save_obj((act_fossil), name_act_fossil)
    save_obj((act_all),name_act_all)
    save_obj((ser), name_ser)
 
    return
# ------------------------------------------------------------------------- 
# main program starts here
# ------------------------------------------------------------------------- 

if __name__ == '__main__':
    
    # -------------------------------------------------------------------------
    # Preparing the input data
    # -------------------------------------------------------------------------

    # Open model file, I just need it to load lons and lats which are used later. 
    # ** there has to be a better way to do this
    rootgrp = Dataset('input/20151116_SR_no2_pm10_pm25/SR_PM25_Y.nc', 'r')
    lon_array = rootgrp.variables['lon'][0, :]
    lat_array = rootgrp.variables['lat'][:, 0]
    n_lon = len(lon_array)  # len(rootgrp.dimensions['longitude'])
    n_lat = len(lat_array)  # len(rootgrp.dimensions['latitude'])  
    rootgrp.close()
           
       
    # Area of each cell in the domain    
    nc_area = 'input/JRC01.nc'
    rootgrp = Dataset(nc_area, mode='r')
    area = rootgrp.variables['surface'][:]
    area_units = rootgrp.variables['surface'].units
    rootgrp.close() 
    
    
    # open area of interest selarea = selected area
    # This is just to check data with Austria (0), Italy (16) comment to consider the area of 
    # interest already defined(see below), NUTS1 - ile de france (41)
#    nc_file_at = 'input/EMI_RED_ATLAS_NUTS1.nc'
    nc_file_at = 'input/EMI_RED_ATLAS_NUTS0.nc'
    rootgrp = Dataset(nc_file_at, mode='r')
    area_area = rootgrp.variables['AREA'][0][:] 
    rootgrp.close()
    # ** Uncomment this to make it work with the area of interest. 
    #    nc_file_reg = path_area_cdf_test 
    #    fh_area = Dataset(nc_file_reg, mode='r')
    #    area_area = fh_area.variables['AREA'][:]
    #    fh_area.close() 
    
    # save the area of interest in a nc file so it can be used later by module 1
    # **this will not be necessary as this file is created by/provided to Sherpa 
    nc_selarea = 'workdir/selarea.nc'
    rootgrp = Dataset(nc_selarea, mode='w')
    rootgrp.createDimension('time', 1)
    rootgrp.createDimension('y', 448)
    rootgrp.createDimension('x', 384)
    lati = rootgrp.createVariable('lat', 'f4', ('y', 'x'))
    longi = rootgrp.createVariable('lon', 'f4', ('y', 'x'))
    selarea = rootgrp.createVariable('AREA', 'f4', ('y', 'x'))
    rootgrp.variables['AREA'].units = '%'
    rootgrp.variables['AREA'].long_name = '% cell area belonging to selected area'
    longi[0, :] = lon_array
    lati[:, 0] = lat_array
    selarea[:] = area_area
    rootgrp.close()

    # -------------------------------------------------------------------------
    # Rading data from Marco's inventory files
    # -------------------------------------------------------------------------
    # Take marco's .tif files and build corresponding arrays with the information on the 
    # pollutants - sector - activity - network
    
    # -------------------------------------------------------------------------
    # Input parameters
    m_sector = 7 # macro sector, SNAP # 7: road transport 
    # List of pollutants of interest  
    pollutant_lst = ['NH3','NOx','VOC', 'SO2', 'PM10', 'CO2']
    # List of subsectors of the m_sector
    sector_lst = ['TRA_RD_LD4C','TRA_RD_HDB','TRA_RD_LD2','TRA_RD_LD4T','TRA_RD_HDT','TRA_RD_M4' ]
    # network (** another option is 'all' but can be extended)
    net_lst = ['all'] # 'rur', 'urb', 'mot' or 'all'
    # list of activity types
    acttype_lst = ['GSL', 'MD', 'GAS','LPG', 'TYRE', 'ABRASION','BRAKE']#, 'TYRE', 'ABRASION','BRAKE'] # 
    nonfuels_lst = ['TYRE', 'ABRASION','BRAKE']

    name_emi= 'emi'
    name_gainsCO2 = 'gainsCO2'
    name_ser= 'ser'
    name_act_fossil='act_fossil'
    name_act_all='act_all'
    # -------------------------------------------------------------------------
    
    # uncomment to read inventories from Marco
    
    readinventories(m_sector, pollutant_lst, sector_lst, net, acttype_lst, nonfuels_lst, name_emi, name_gainsCO2, name_act_fossil,name_act_all,name_ser)
    # -------------------------------------------------------------------------
    # Measures implementation (**will be translated to a def)
    # -------------------------------------------------------------------------
    
#    # load results of readinventories
    emi=load_obj(name_emi)
    gainsCO2=load_obj(name_gainsCO2)
    ser=load_obj(name_ser)
    act_all=load_obj(name_act_all)
    act_fossil=load_obj(name_act_fossil)
    
#    # initialization
#    ser_tot = {} # total service [Mvkm] in the area of interest 
#    ef_inv ={} # emission factor inventory (Marco)
#    act_tot = {} # total activity [PJ] in the area of interest

    # calculate total activity in area of 
    act_all_area_t ={}
    act_fossil_area_t ={}
    for sector in act_all:
        for aty in act_all[sector]:
            for net in act_all[sector][aty]:    
                act_all_area_t[sector,aty,net] = np.sum(act_all[sector][aty][net] * area_area / 100)  # [PJ] 
                act_fossil_area_t[sector,aty,net] = np.sum(act_fossil[sector][aty][net] * area_area / 100)  # [PJ] 
                
    emi_area_t ={}
    for pollutant in emi:
        for sector in emi[pollutant]:
            for aty in emi[pollutant][sector]:
                for net in emi[pollutant][sector][aty]:
                    emi_area_t[pollutant,sector,aty,net] = np.sum(emi[pollutant][sector][aty][net] * area_area / 100)
            
    ser_area_t={}
    for sector in ser:
            for net in ser[sector]['TYRE']:
                ser_area_t[sector,net] = np.sum(ser[sector]['TYRE'][net] * area_area / 100) # [Mvkm]
    
    ef_area_t={}
    for pollutant in emi:
            for sector in emi[pollutant]:
                for aty in emi[pollutant][sector]:
                    for net in emi[pollutant][sector][aty]:
                        if aty in [aty for aty in acttype_lst if aty not in nonfuels_lst]:
                            if pollutant is not 'CO2':
                                ef_area_t[pollutant,sector,aty,net]=emi_area_t[pollutant,sector,aty,net]/act_all_area_t[sector,aty,net]
                            elif pollutant is 'CO2': 
                                ef_area_t[pollutant,sector,aty,net]=emi_area_t[pollutant,sector,aty,net]/act_fossil_area_t[sector,aty,net]
                        if aty in nonfuels_lst:
                            ef_area_t[pollutant,sector,aty,net]=emi_area_t[pollutant,sector,aty,net]/ser_area_t[sector,net]

##    # technical measures                      
##    # New emission factors after the implementation of measures:
##    # first option - reduction of the emission factors for each sector/activity
##    # read csv file with the emission factors
    df = pd.read_csv('input/ef_red_sherpaeco.csv',  index_col=[0,1], names=['POLL','ACT'].extend(sector_lst), skipinitialspace=True)
##    # second option (to be implemented**)- emission factors for the best available technology
##    
#    # read the precursors list (as in model1)  
    rootgrp = Dataset(path_model_cdf_test, 'r')
    precursor_lst = getattr(rootgrp, 'Order_Pollutant').split(', ')
    # create a dictionary with emissions per precursor, macrosector and postion (lat, lon)
    emission_dict = create_emission_dict(path_emission_cdf_test, precursor_lst)
#    
   # caluclate total emission in the area of interest for the base case (Sherpa's emission inventory)
    em_bc ={} # emissions for the base case
    for precursor in precursor_lst:
        em_bc[precursor,m_sector-1] = np.sum(emission_dict[precursor][m_sector-1] * area * area_area / 100)
    
    # calculate new emissions in the area of interest
    em_new ={} # emission after measure implementation    
    for pollutant in emi:
#        em_new[pollutant]={}
        for sector in emi[pollutant]:    
#            em_new[pollutant][sector] = {}
           for aty in emi[pollutant][sector]:
#                em_new[pollutant][sector][fuel] = {}
                for net in emi[pollutant][sector][aty]:
#                    em_new[pollutant][sector][fuel][net] = {}
                    ef = ef_area_t[pollutant,sector,aty,net] *(1 - df.loc[pollutant,aty][sector])
                    if aty not in nonfuels_lst:
                        if pollutant is not 'CO2':
                            em_new[pollutant,sector,aty,net] = (act_all[sector][aty][net] * ef * area_area / 100)
                        elif pollutant is 'CO2':
                            em_new[pollutant,sector,aty,net] = (act_fossil[sector][aty][net] * ef * area_area / 100)
#                        a = (np.sum(em_new[pollutant][sector][fuel][net]) - np.sum(em_tot[pollutant][sector][fuel][net]))/np.sum(em_tot[pollutant][sector][fuel][net])
#                        print(a, pollutant, sector, fuel, net)

                    else:
                        em_new[pollutant,sector,aty,net] = (ser[sector]['TYRE'][net]  * ef * area_area / 100)
#                        a = (np.sum(em_new[pollutant][sector][fuel][net]) - np.sum(em_tot[pollutant][sector][fuel][net]))/np.sum(em_tot[pollutant][sector][fuel][net])
#                        print(a, pollutant, sector, fuel, net)
#                        
##    # -------------------------------------------------------------------------
##    # Running module1 
##    # --------------------------------------------------------`-----------------
##                    
#    em_new_sum={}
#    em_tot_sum={}
#    em_new_percell={}
#    em_tot_percell={}
#    for pollutant in pollutant_lst: 
#        em_new_percell[pollutant] = np.sum(np.sum(np.sum(em_new[pollutant][sector][fuel][net] for net in em_new[pollutant][sector][fuel]) for fuel in em_new[pollutant][sector]) for sector in em_new[pollutant])
#        em_tot_percell[pollutant] = np.sum(np.sum(np.sum(em_tot[pollutant][sector][fuel][net] for net in em_tot[pollutant][sector][fuel]) for fuel in em_tot[pollutant][sector]) for sector in em_tot[pollutant])
#        em_tot_sum[pollutant] = np.sum(em_tot_percell[pollutant])
#        em_new_sum[pollutant] = np.sum(em_new_percell[pollutant])
#        act_tot_sum = np.sum(np.sum(np.sum(act_tot[sector][fuel][net] for net in act_tot[sector][fuel]) for fuel in act_tot[sector]) for sector in act_tot)
#    delta_emission_dict={}
#
#    for precursor in precursor_lst:     
#        for pollutant in pollutant_lst:
#            if pollutant == precursor:
#                delta_emission_dict[precursor] = (em_tot_percell[pollutant]-em_new_percell[pollutant])
##                checkdif[precursor]=(em_tot[pollutant]-em_bc[precursor,m_sector-1])/em_bc[precursor,m_sector-1]*100
#        if precursor == 'NMVOC':
#                delta_emission_dict[precursor] = (em_tot_percell['VOC']-em_new_percell['VOC'])
##                checkdif[precursor]=(em_tot['VOC']-em_bc[precursor,m_sector-1])/em_bc[precursor,m_sector-1]*100
#        if precursor == 'SOx': 
#                delta_emission_dict[precursor] = (em_tot_percell['SO2']-em_new_percell['SO2'])
##                checkdif[precursor]=(em_tot['SO2']-em_bc[precursor,m_sector-1])/em_bc[precursor,m_sector-1]*100
#        if precursor == 'PPM':
#                delta_emission_dict[precursor] = (em_tot_percell['PM10']-em_new_percell['PM10'])
##                checkdif[precursor]=(em_tot['PM10']-em_bc[precursor,m_sector-1])/em_bc[precursor,m_sector-1]*100
#    
#    for precursor in precursor_lst:
#        delta_emission_dict[precursor] = np.sum((delta_emission_dict[precursor]/area), axis=0)
#        
#    red = {}
##    checkdif={}   
##    for precursor in precursor_lst: 
##        red[precursor]={}
###        checkdif[precursor]={}
##        for pollutant in pollutant_lst:
##            if pollutant == precursor:
##                red[precursor] = (em_tot_sum[pollutant]-em_new_sum[pollutant])/em_bc[precursor,m_sector-1]*100
###                checkdif[precursor]=(em_tot[pollutant]-em_bc[precursor,m_sector-1])/em_bc[precursor,m_sector-1]*100
##        if precursor == 'NMVOC':
##                red[precursor] = (em_tot_sum['VOC']-em_new_sum['VOC'])/em_bc[precursor,m_sector-1]*100
###                checkdif[precursor]=(em_tot['VOC']-em_bc[precursor,m_sector-1])/em_bc[precursor,m_sector-1]*100
##        if precursor == 'SOx': 
##                red[precursor] = (em_tot_sum['SO2']-em_new_sum['SO2'])/em_bc[precursor,m_sector-1]*100
###                checkdif[precursor]=(em_tot['SO2']-em_bc[precursor,m_sector-1])/em_bc[precursor,m_sector-1]*100
##        if precursor == 'PPM':
##                red[precursor] = (em_tot_sum['PM10']-em_new_sum['PM10'])/em_bc[precursor,m_sector-1]*100
###                checkdif[precursor]=(em_tot['PM10']-em_bc[precursor,m_sector-1])/em_bc[precursor,m_sector-1]*100
#    for precursor in precursor_lst: 
#        red[precursor]={}
##        checkdif[precursor]={}
#        for pollutant in pollutant_lst:
#            if pollutant == precursor:
#                red[precursor] = (em_tot_sum[pollutant]-em_new_sum[pollutant])/em_tot_sum[pollutant]*100
##                checkdif[precursor]=(em_tot[pollutant]-em_bc[precursor,m_sector-1])/em_bc[precursor,m_sector-1]*100
#        if precursor == 'NMVOC':
#                red[precursor] = (em_tot_sum['VOC']-em_new_sum['VOC'])/em_tot_sum['VOC']*100
##                checkdif[precursor]=(em_tot['VOC']-em_bc[precursor,m_sector-1])/em_bc[precursor,m_sector-1]*100
#        if precursor == 'SOx': 
#                red[precursor] = (em_tot_sum['SO2']-em_new_sum['SO2'])/em_tot_sum['SO2']*100
##                checkdif[precursor]=(em_tot['SO2']-em_bc[precursor,m_sector-1])/em_bc[precursor,m_sector-1]*100
#        if precursor == 'PPM':
#                red[precursor] = (em_tot_sum['PM10']-em_new_sum['PM10'])/em_tot_sum['PM10']*100
##                checkdif[precursor]=(em_tot['PM10']-em_bc[precursor,m_sector-1])/em_bc[precursor,m_sector-1]*100
#                
#    area_tot_em_dict = {}
#    for precursor in precursor_lst:
#        for snap in range(1, 11):
#            area_tot_em_dict[precursor, snap - 1] = np.sum(emission_dict[precursor][snap - 1]  * area * area_area / 100)  # [Mg] 
#    
#
#    reductions = {}
#    reductions[m_sector-1]={}
##                                   NOx	        NMVOC       NH3           PPM        SOx
#    reductions[m_sector-1] = [red['NOx'], red['NMVOC'], red['NH3'] , red['PPM'], red['SOx']]
#    path_reduction_txt='input/sherpaeco_reduction.txt'
#    write_reductions(path_reduction_txt, reductions[m_sector-1])  
##
##    # if it doesn't exist start=0 and dividsor=1
##    progresslog = 'input/progress.log'
##    start = time()
##    output = 'output/'
##    output2 = 'output2/'
##    # run module 1 with progress log
##    proglog_filename = path_result_cdf_test + 'proglog'
##    write_progress_log(proglog_filename, 25, 2)
##    start = time()
##    shrp.module1(path_emission_cdf_test, nc_selarea,
##                    path_reduction_txt, path_model_cdf_test, output)
##    shrp1.module1(path_emission_cdf_test, nc_selarea,
##                    delta_emission_dict, path_model_cdf_test, output2)  
##    stop = time()
##    print('Module 1 run time: %s sec.' % (stop-start))
##    remove(proglog_filename)
##   
####     -------------------------------------------------------------------------
####     (Cost) Benefit Analysis 
####     -------------------------------------------------------------------------
###
##    deltaconc='output/delta_concentration.nc'
##    deltayll_reg, deltayll_tot, delta_mort_tot, delta_mort_reg, delta_yll_pp = calc_impacts(deltaconc, area_area)
##    
##    deltaconc2='output2/delta_concentration.nc'
##    deltayll_reg2, deltayll_tot2, delta_mort_tot2, delta_mort_reg2, delta_yll_pp2 = calc_impacts(deltaconc2, area_area)
##    
##    # -------------------------------------------------------------------------
##    # Output of results (todo)
##    # -------------------------------------------------------------------------
#  
##    for sec in sector_lst: dct[sec] = {}
#    [(x, y) for x in [1,2,3] for y in [3,1,4] if x != y]
##    dict([('sape', 4139), ('guido', 4127), ('jack', 4098)])
##    dct = defaultdict(dict)
#
#    dct1 = {}
## dict(pol,{}) 
#     #    dct[pollutant] = {} for pollutant in pollutant_lst
#    dct1 =[dict([((pol,sector,aty,net), {})]) 
#        for pol in pollutant_lst for sector in sector_lst 
#        for aty in acttype_lst for net in net_lst]
#        
#    dct2 = {}
##    for pol in pollutant_lst:
##        for sector in sector_lst: 
##            for aty in acttype_lst:  
##                for net in net_lst:
##                    dct2[pol,sector,aty,net] ={}
##setdefault(
##    d = dict()
##    d =dict(d.setdefault(pol, {}).setdefault(sector, {}).setdefault(aty, {}) for pol in pollutant_lst 
##        for sector in sector_lst 
##        for aty in acttype_lst)
#        
##     dct.update(dict(pol,{}) for pol in pollutant_lst)  
##    dct = [dict(pol, {}) for pol in pollutant_lst]   
##    dct['CO2']=3
##    em ={} # emission after measure implementation    
##    for pollutant in pollutant_lst:
##        em[pollutant]={}
#
##    target_dict = defaultdict(dict)
##    target_dict['PM10']['TRA_RD_HDT'] = 5
#
##    [aty for aty in acttype_lst if aty not in nonfuels_lst]

           