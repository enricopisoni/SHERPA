# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 16:46:00 2017

Allows to create the tables Marco needs to apply to his rasters the urban, motorway and non urban shares. 

Inputs:
    # define aggregation list (see excel file EMISSIONS_comparisonWithMarco_v2 for definitions)
    shares_lst=[1, 2, 3, 4, 5, 6, 7, 8]
    # pollutants of interest in TREMOVE: 
    pollutant_lst =  ['CO','NOx','VOC', 'SO2', 'PM', 'CO2', 'PM']
    # path to the folder where you have you excel file from TREMOVE: EMISSIONS_comparisonWithMarco_v2
    dirpath='D:/emissioninventory'
     
Output:
    all the the relevant tables in xlsx and csv formats
    
@author: peduzem
"""
import os
# for importing matlab files
import scipy.io as sio
import pandas as pd # conda install pandas
from openpyxl import load_workbook
import simpledbf as sdbf # pip install D:\pipinstall\simpledbf-0.2.6.tar.gz
import dbf #     pip install D:\pipinstall\aenum-1.4.7-py3-none-any.whl
           #     pip install D:\pipinstall\dbf-0.96.8-py3-none-any.whl

def readtre(df, share, pollutant, polltype, directory):
          
    df1=df.loc['Sum of {} {}'.format(pollutant, polltype).rstrip()]
    # change name of columns as in Marco's files
    df1 = df1.rename(columns={'urban road': 'URBAN', 'motorway': 'MOTORWAY', 'non urban road': 'NOT_URBAN'}) 
    # remove motorways from cyprus and malta
    df1.set_value('CY', 'MOTORWAY', 0) 
    df1.set_value('MT', 'MOTORWAY', 0)  
    # fill nan values with 0s
    df1.fillna(0, inplace=True)
    # calculate shares per road type
    df2=df1.apply(lambda x:100 * x / x.sum(), axis=1)	
    # change name of columns as in Marco's files
    df2 = df2.rename(columns={'URBAN': 'URBAN__', 'MOTORWAY': 'MOTORWAY__', 'NOT_URBAN': 'NOT_URBA_1'})
    # concatenate the two dataframes
    result = pd.concat([df1, df2], axis=1)
    # rearrenge column order as in Marcos files
    cols = result.columns.tolist()
    cols_ord = [cols[0]] + [cols[2]] + [cols[1]] + [cols[3]] + [cols[5]] + [cols[4]]
    result_ord=result[cols_ord]
    if polltype == 'non-exhaust':
        sfx='ne'
    else:
        sfx=''
    # Write to excel 
    result_ord.to_excel(directory + '/TRE{}{}.xlsx'.format(share, sfx))            
    # convert excel to csv
    data_xls = pd.read_excel(directory + '/TRE{}{}.xlsx'.format(share, sfx), 'Sheet1')
    data_xls.to_csv(directory + '/TRE{}{}.csv'.format(share, sfx), index=None) #.replace("-", "")
    return
               
if __name__ == '__main__':
    # define aggregation list (see excel file EMISSIONS_comparisonWithMarco_v2 for definitions)
    shares_lst=[1, 2, 3, 4, 5, 6, 7, 8]
    # pollutants of interest in TREMOVE: 
    pollutant_lst =  ['CO','NOx','VOC', 'SO2', 'PM', 'CO2', 'PM', 'vkm']
    # path to the folder where you have your excel file from TREMOVE
    dirpath='D:/emissioninventory'

    for share in shares_lst:
        # create pandas dataframe from excel file
        df=pd.read_excel(dirpath+'/EMISSIONS_comparisonWithMarco_v2.xlsm', 
                         sheetname='TRE{}'.format(share), skiprows=16,  header=[1],  
                         index_col=[0,1]).sortlevel() #,index_col='sig_name'
        for pollutant in pollutant_lst:     
            # create directories for the excel files, one for each pollutant
            directory=dirpath+'/{}_tre'.format(pollutant)
            try:
                os.stat(directory)
            except:
                os.mkdir(directory)    
            # for PM     
            if pollutant == 'PM':
                # 1 and 2 should only have exhaust part, 4 only non-exhaust
                if share != 1 and share != 2 and share != 4: 
                    polltype='exhaust'
                    readtre(df, share, pollutant, polltype, directory)     
                    polltype='non-exhaust'
                    readtre(df, share,pollutant, polltype, directory)
                # because the non-exhaust part is in 4
                elif share == 4: 
                    polltype='non-exhaust'
                    readtre(df, share, pollutant, polltype, directory)
                # all others should only have exhaust
                else:
                    polltype='exhaust'
                    readtre(df, share, pollutant, polltype, directory)  
            # all other pollutants only have exhaust    
            if pollutant == 'vkm':
                if share in [3, 4, 5, 6, 7, 8]:
                    polltype=''
                    readtre(df, share, pollutant, polltype, directory)        
            else:
                polltype='exhaust'
                readtre(df, share, pollutant, polltype, directory)
                