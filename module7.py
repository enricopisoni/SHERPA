'''
NAME
    Reads SHERPA ncdf file with Python
PURPOSE
    To read matrix data and put them in a dataframe of vectors
PROGRAMMER(S)
    Denise Pernigotti
REVISION HISTORY
    
REFERENCES
    
'''
from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/  
import numpy as np
import pandas as pd  
def read_list_nc(nc_file):
    #nc_file='input/20151116_SR_no2_pm10_pm25/BC_emi_PM25_Y.nc'
    #nc_file='input/20151116_SR_no2_pm10_pm25/SR_PM25_Y.nc'
    nc_data = Dataset(nc_file, 'r') 
    nc_dims = [dim for dim in nc_data.dimensions]
    nc_vars = [var for var in nc_data.variables]
    #sometimes the latitude is written just with lat as in model data
    latname=filter (lambda x: 'lat' in x, nc_vars)[0]
    lats = nc_data.variables[latname][:] 
    lonname=filter (lambda x: 'lon' in x, nc_vars)[0]
    lons = nc_data.variables[lonname][:] 
    #if there are three dimensional arrays
    if len(nc_dims)==3:
        ncz=str(list(set(nc_dims)-set(['latitude','longitude']))[0])
        nz=range(len(nc_data.dimensions[ncz]))
        if ncz=='pollutant':
            strpoll=nc_data.Order_Pollutant
            nznames=strpoll.split(', ')
        else:
            nznames=[ncz + s for s in map(str,range(1,len(nc_data.dimensions[ncz])+1))]         
    #create an index with lat and lon
    #latrep=map(str, np.repeat(lats,len(lons)))
    #lonrep=map(str, np.tile(lons,len(lats)))
    #trasform variables arrays in vectors
    #allvar={'lat_lon':map(lambda (x,y): x+'_'+y, zip(latrep, lonrep))}
    #create lat and lon info
    if len(lats.shape)==2 and len(lons.shape)==2:
        nrow=lats.shape[0] 
        ncol=lats.shape[1]
        lon=lons.ravel()
        lat=lats.ravel()
    else:
        nrow=len(lats)
        ncol=len(lons)
        lon=np.tile(lons,nrow)
        lat=np.repeat(lats,ncol)
         
    col=map(str, np.repeat(range(1, nrow+1),ncol))
    row=map(str, np.tile(range(1, ncol+1),nrow)) 
    index=map(lambda (x,y): x+'_'+y, zip(row, col))
    #index=range(0,ncol*nrow)
    allvar={'grid':index,'lat': lat,'lon':lon}
    nc_vars.remove(latname)
    nc_vars.remove(lonname)
    for var in nc_vars:
        if len(nc_dims)==3:
            for sn in nz:
                var_sn=nc_data.variables[var][sn]
                #print 'var=',var,' zval=',nznames[sn],' shape=',var_sn.shape,' len=',len(var_sn.ravel())
                if ncz=='pollutant':
                    allvar[var+'_'+nznames[sn]]=var_sn.ravel() 
                else:
                    allvar[nznames[sn] +'_'+ var]=var_sn.ravel() 
        else:
            vardata=nc_data.variables[var][:]
            #print 'var=',var,' shape=',vardata.shape,' len=',len(vardata.ravel())
            allvar[var]=vardata.ravel()
    df = pd.DataFrame(data=allvar, index=index)
    return df

'''
NAME
    Implementation of haversine formula (form lat lon to distances in km) for vectors
    Calculate the great-circle distance between two points on the Earth surface.
PURPOSE
    Implementation of haversine formula (form lat lon to distances in km) for vectors
    :input: one 2-tuples, and a vector of 2-tuples containing the latitude and longitude of a point
    in decimal degrees and a vector.
    Example: haversine((45.7597, 4.8422), (lat, lon))
    :output: Returns the distance in km bewteen the the point to all other points.PROGRAMMER(S)
    Denise Pernigotti from https://github.com/mapado/haversine/blob/master/haversine
REVISION HISTORY
    
REFERENCES
    
'''
import numpy as np
AVG_EARTH_RADIUS = 6371  # in km


def haversine_vec(lon1,lat1,lon_vec2,lat_vec2):

    # calculate haversine
    dlat = np.radians(lat_vec2) - np.radians(lat1)
    dlon = np.radians(lon_vec2) - np.radians(lon1)
    d = np.sin(dlat * 0.5) ** 2 + np.cos(np.radians(lat1)) * np.cos(np.radians(lat_vec2)) * np.sin(dlon * 0.5) ** 2
    h = 2 * AVG_EARTH_RADIUS * np.arcsin(np.sqrt(d))
    return h # in kilometers

'''
NAME
    Estimates distance in grid points
PURPOSE
    Estimates distance in grid points when dlat_res and dlon_res are constants
PROGRAMMER(S)
    Denise Pernigotti
REVISION HISTORY
    
REFERENCES
    
'''
import numpy as np
def cell_dist(lon1, lat1, lon2, lat2,dlon_res,dlat_res):
    dist=np.sqrt(((lon2 - lon1)/dlon_res)  ** 2 + ((lat2 - lat1)/dlat_res)  ** 2) 
    return dist   

    
'''
NAME
    Vectorized estimates distance in grid points using index
PURPOSE
    Estimates distance in grid points using indexes
PROGRAMMER(S)
    Denise Pernigotti
REVISION HISTORY
    
REFERENCES
    
'''
import numpy as np
def cell_dist_index_vec(x1_y1,x2_y2_vec):
    x2_y2=pd.DataFrame(x2_y2_vec.str.split('_',expand=True).tolist(), columns=['x2','y2'], index=x2_y2_vec)
    dx=pd.to_numeric(x2_y2['x2'])-float(x1_y1.split('_')[0])
    dy=pd.to_numeric(x2_y2['y2'])-float(x1_y1.split('_')[1])
    dist=np.sqrt(dx ** 2 + dy ** 2) 
    return dist   

    
'''
NAME
    Import info on grid points attribution to nuts or specific area type from ascii file
PURPOSE 
    Import info on grid points attribution to nuts or specific area type from ascii file/s. 
    If the file is single then it must contain the column 'Area [km2]' relative to % of the area in the finest nut, 
    this datum will be set to each nut but it will then aggregated for larger nuts when nutsarea will be calculated
    If the files are two, then each nut will have its own % area for each grid point, then the data will be merged here
PROGRAMMER(S)
    Denise Pernigotti
REVISION HISTORY
    
REFERENCES
    
'''    
import pandas as pd  
def read_nuts_area(filenuts):
    nuts_def= path_input + filenuts +'.txt'
    nuts_info = pd.read_csv(nuts_def,delimiter="\t")
    nuts_info=nuts_info.dropna(axis=1,how='all')
    nutsnames=list(nuts_info.columns[~nuts_info.columns.isin(['POP','COL','ROW','Area [km2]','LAT','LON'])])
    nutsnames.insert(0, 'ALL_'+nutsnames[0])
    nuts_info[nutsnames[0]]=nutsnames[0] 
    nuts_info['grid']=['_'.join(str(i) for i in z) for z in zip(nuts_info['COL'],nuts_info['ROW'])]      
    if 'Area [km2]' in nuts_info.columns:
        nuts_area=pd.concat(map(lambda p: nuts_info['Area [km2]'],nutsnames),axis=1)
        #nuts_area.index=nuts_info['grid']
        nuts_area.columns=nutsnames  
       #nuts_info=nuts_info[nutsnames]
    else:
        sys.exit("missing infos on grid cells area per nut")

    #aggregate data for each nut, create a dictionary
    nuts_info_all={}
    for nut in nutsnames:
        #create a multindex
        index = pd.MultiIndex.from_tuples(list(zip(nuts_info['grid'],nuts_info[nut])), names=['grid', 'nutname'])
        nut_info=pd.Series(list(nuts_area[nut]), index=index)
        nut_info=nut_info.to_frame(name='area')
        #aggregate data on these nuts if necessary
        nut_info_nut=nut_info.groupby(level=[0,1]).sum()    
        #find total area
        grid_area_tot=nut_info_nut.groupby(level=['grid']).sum()
        nut_info_nut['parea']=100.*nut_info_nut/grid_area_tot
        nuts_info_all[nut]=nut_info_nut
    return nuts_info_all
    
'''
NAME
    Core of sherpa model calculation
PURPOSE
    Calculate the delta conc for one precursor on all snaps
PROGRAMMER(S)
    Denise Pernigotti
REVISION HISTORY
    
REFERENCES
    
'''
import numpy as np
def sherpa_model(prec,model_point,inner_radius,emissions):
    #prec='NOx'
    #model_point=model.loc[target_idx]
    alpha_point=model_point['alpha_'+ prec]
    omega_point=model_point['omega_'+ prec]
    flatWeight_point=model_point['flatWeight_'+ prec]
    alldists_point=emissions['dist_cells']
    #emissions['DeltaE'] act as a mask on grid points: 1 on selected area and 0 elsewhere
    #print "SR model for precursor ",prec,"on snaps ", filter (lambda x: '_'+ prec in x, list(emissions))
    emissiondata=emissions[filter (lambda x: '_'+ prec in x, list(emissions))].multiply(emissions['deltaE'],axis=0)
    window_all=(1.+alldists_point)**(-omega_point)
    #apply model to all snaps 
    dc_cells=alpha_point*emissiondata.multiply(window_all,axis=0)
    #apply flatweight approximation for coherence with Sour Receptor relations, even if it is not needed here 
    #identify which points require the approximation
    outer_cells=emissions[(emissions['dist_cells']>inner_radius) & (emissions['deltaE']>0)].index
    inner_cells=emissions[emissions['dist_cells']<=inner_radius].index
    #calculate approximate total value if there are points with non null emissions out of inner radius
    if len(outer_cells)>0:
        print "fo prec ",prec,"applying flatweight approximation on ",len(outer_cells)," points beyond ",inner_radius,"cells"
        weighted_emissions_flat = alpha_point*(flatWeight_point * emissiondata.loc[outer_cells].sum())
    #disaggregate the value on all involved points
        flat_value=weighted_emissions_flat/len(outer_cells)
        flat_matrix=pd.DataFrame(index=outer_cells, columns=flat_value.index)
        flat_matrix=flat_matrix.fillna(flat_value)
        dc_cells.loc[outer_cells,:]=flat_matrix
    return dc_cells       

# main 
if __name__ == '__main__':
    
    import timeit
    import sys
    import os
    import pandas as pd  
    import matplotlib.pyplot as plt
    import openpyxl
    path_work='C:\\Users\\pernide\\AppData\\Local\\Sherpa\\app\\data\\'
    path_input=path_work + 'input\\'
    path_output=path_work + 'output\\'
    #from geopy.distance import vincenty
    #emission_nc = 'input/20151116_SR_no2_pm10_pm25/BC_emi_PM25_Y.nc'
    #concentration_nc = 'input/20151116_SR_no2_pm10_pm25/BC_conc_PM25_Y.nc'
    #model_nc = 'input/20151116_SR_no2_pm10_pm25/SR_PM25_Y.nc'
    emission_nc = path_input + 'base_emissions\\BC_emi_PM25_Y.nc'
    concentration_nc = path_input + 'base_concentrations\\BC_conc_PM25_Y.nc'
    model_nc = path_input + 'source_receptors\\SR_PM25_Y.nc'
 
    #grab nuts/are info from txt
    nuts_info=read_nuts_area('selection\\grid_intersect')
    nuts_info.update(read_nuts_area('nuts_selection\\grid_int_gcities'))
    nuts_info.update(read_nuts_area('nuts_selection\\grid_int_fua'))

    
    #read netcdf files, put the data in two dataframes and check consistency 
    emissions = read_list_nc(emission_nc) 
    concentration = read_list_nc(concentration_nc) 
    model= read_list_nc(model_nc)
    inner_radius = int(getattr(Dataset(model_nc,'r'), 'Radius of influence'))
    model=concentration.merge(model,left_on='grid', right_on='grid', how='outer')
    model.index=model['grid']
    model = model.rename(columns={'lon_x': 'lon', 'lat_x': 'lat'})
    #remove duplicated columns for latitude and longitude in model
    if model['lat'].equals(model['lat_y']) and model['lon'].equals(model['lon_y']): 
        model=model.drop(['lat_y','lon_y'],axis=1)      
    else:
        sys.exit("latitude is different in loaded matrices")

    #check consistency of model and emissions
    if model['lat'].equals(emissions['lat']) and model['lon'].equals(emissions['lon']): 
        print 'OK latitude and longitude in matrices emissions, model and conc are the same'
    else:
        sys.exit("latitude is different in loaded matrices")    
        
    #from modelled points remove unused rows (grid cells where alpha is zero, model is not defined)
    alpha_ind=filter (lambda x: 'alpha' in x, list(model))
    print 'Remove ',len(np.where(model[alpha_ind].sum(axis=1)==0)[0]), 'rows(cells) with all alphas to zero'
    model.drop(model[model[alpha_ind].sum(axis=1)==0].index,inplace=True)
    
    #identify precursors not modelled
    alpha_p=model[alpha_ind].sum(axis=0)
    zero_prec=[alpha_p.index[i][0].split('_')[1] for i in np.where(alpha_p==0)]
    #remove model and emicsions data  for unused precursors
    if len(zero_prec)>0:
        print 'Remove columns associated with a non modelled precursor in model ',[s for s in model.columns if any(xs in s for xs in zero_prec)]
        model.drop([s for s in model.columns if any(xs in s for xs in zero_prec)],inplace=True,axis=1)
        print 'Remove columns associated with a non modelled precursor in emissions ',[s for s in emissions.columns if any(xs in s for xs in zero_prec)]
        emissions.drop([s for s in emissions.columns if any(xs in s for xs in zero_prec)],inplace=True,axis=1)
    
    #precursors names from model
    precursors_alpha= np.unique([i.split('_')[1] for i in filter(lambda x: 'alpha' in x, list(model))]) 
    #snaps names and precursors names from emissions
    snaps=np.unique([i.split('_')[0] for i in filter(lambda x: 'Nsnaps' in x, list(emissions))]) 
    precursors=np.unique([i.split('_')[1] for i in filter(lambda x: 'Nsnaps' in x, list(emissions))]) 
    if np.array_equal(precursors,precursors_alpha):
        print 'OK model and emissions contain the same precursors'
    else:
       sys.exit("model and emissions to not contain the same precursors ") 

    #for emissions data remove points with all emissions set to zero
    #emissions_ind=filter (lambda x: 'Nsnaps' in x, list(emissions))
    #print 'Remove ',len(np.where(emissions[emissions_ind].sum(axis=1)==0)[0]), 'rows(cells) with all all emissions to zero'
    #emissions.drop(emissions[emissions[emissions_ind].sum(axis=1)==0].index,inplace=True)

     #set the point to be tested
    #test_point={'lon':8.8251,'lat':45.8206}#varese
    #test_point={'lon':8.8125,'lat':45.84375}
    #test_point={'lon':9.191383,'lat':45.464211} #Milan
    #test_point={'lon':9.1875,'lat':45.46875} #closer grid point of Milan
    test_point={'lon':2.349014,'lat':48.864716} #Paris 
    #find grid resolution in degrees 
    dlon_res=min(v for v in [abs(j-i) for i, j in zip(emissions['lon'][:-1], emissions['lon'][1:])] if not v in (None,0))
    dlat_res=min(v for v in [abs(j-i) for i, j in zip(emissions['lat'][:-1], emissions['lat'][1:])] if not v in (None,0))
    #find distance from all other points in emissions data in km
    #emissions['dist_greatc']=map(lambda x:great_circle((test_point['lat'],test_point['lon']),(emissions.loc[x,'lat'],emissions.loc[x,'lon'])).km,emissions.index)
    emissions['dist_haver']=haversine_vec(test_point['lon'],test_point['lat'],emissions['lon'],emissions['lat'])
    #find distance from all other points in cell number
    target_idx=emissions['dist_haver'].idxmin()
    emissions['dist_cells']=cell_dist(emissions.loc[target_idx,'lon'],emissions.loc[target_idx,'lat'],emissions['lon'],emissions['lat'],dlon_res,dlat_res)
    #emissions['dist_cells_idx']=cell_dist_index_vec(target_idx,emissions.index)
    print "closer grid point to ",test_point," is ",zip(emissions.loc[target_idx,['lon','lat']])
         
    #identify point within predefined radii
    #test_nut=[20,80,1000,1000000]
    #identify points with NETCDF area specification
    #test_nut=nuts_names
    all_target_nut={}
    alldc_snaps={}
    nutnames=nuts_info.keys()
    nutnames.sort()
    for rads in nuts_info.keys():
        #identify the e nut of the target, eventually the one with the bigger area in
        target_nut=list(nuts_info[rads].loc[target_idx,].sort_values('parea',ascending=False).index)[0]
        all_target_nut[rads]=target_nut
        print "for nut ",rads," the tested point is for the majority in ",target_nut
        #select grid point in the target_nut
        narea_data=nuts_info[rads].loc[(slice(None),target_nut),'parea']
        narea=pd.Series(list(narea_data),index= narea_data.index.get_level_values('grid'))
        rad1=narea.index
        #rad1=nuts_parea.index[np.where(nuts_parea[rads]>0)]
        print 'there are ',len(rad1),' points with nut ',rads,'=',target_nut
         #create a fake emission to zero
        emissions['deltaE']=0.
         #identify which points to apply a deltaE=E
        emissions.loc[rad1,'deltaE']=1.
        #add the calulations on area percentage in nuts
        #careful rad1 is not a unique index in nuts_parea and nuts_info, a point may have multiple nuts in %
        #emissions.loc[rad1,'deltaE']= emissions.loc[rad1,'deltaE']*nuts_parea[rads][nuts_info[rads]==test_nut[rads]]/100.
        emissions.loc[rad1,'deltaE']= emissions.loc[rad1,'deltaE']*narea[rad1]/100.
        #apply model, create a dataframe, with a column for each snap and each precursor
        dc=pd.concat(map(lambda p: sherpa_model(p, model.loc[target_idx],inner_radius,emissions),precursors),axis=1)
        #aggregate data summarizing dc
        alldc=dc.sum()
        #aggregate on snaps 
        alldc_snaps[rads]=map(lambda s: alldc[filter (lambda x: s+'_' in x, list(dc))].sum(),snaps)
     #alldc_snaps.columns=all_target_nut 
    
    alldc_snaps=pd.DataFrame.from_dict(alldc_snaps)    
    print model.loc[target_idx,'conc']
    print alldc_snaps
    alldc_snaps_perc=alldc_snaps*100/model.loc[target_idx,'conc']
    alldc_snaps_perc.T.plot.barh(stacked=True,xlim=(0,100),legend=False);
    print "max perc ",alldc_snaps_perc['ALL_NUTS_Lv0'].sum()
    #load an existing xls file
    #xfile = openpyxl.load_workbook('test.xlsx')
    #sheet = xfile.get_sheet_by_name('Sheet1')
    #sheet['A1'] = 'hello world'
    #xfile.save('text2.xlsx')

    
    
    
