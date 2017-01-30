'''
Created on Jun 23, 2015

simulates the main tool developped by teraria and call the source-receptor model

The sherpa tool can be ran in different modes

1) Module 1: scen_nuts
----------------------
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

2) Module 2: Scen_Atlas
-----------------------
This module calculates the concentration reduction in each cell of a nuts sector due to emission reductions in that same nuts sector.
Input: - baseline emissions (per precursor and cell)
       - a netcdf defining the areas where emission reductions will be applied
       - txt file with reductions (per snap and precursor) to be applied in each nuts
       - model parameters
       - result path
Output: - netcdf with concentration reductions, one layer per NUTS

3) Module 3: source apportionment
---------------------------------
Module 3 executes module 4 for for:
- all macrosectros individually and all together for 1 or more precursor => module 3a
- all precursors individually and all together for 1 or more macrosectors => module 3b

The calculated potencies are:
- DC / C
- DC / C / alfa

4) Module 4: potencies
-----------------------
Calculate potencies for an emission reduction in a given area. 

5) Module 5: potency atlas
--------------------------

6) Module 6: radius calculation
-------------------------------

A target cell is selected. For this cell the concentration changes is calculated due to a defined emission reduction
applied in each NUTS individually
Input: -

Output: - text file with NUTS code and corresponding concentration change in the target cell.
        - netcdf with each NUTS with the vallue of its concentration change caused in the target cell.


@author: degraba
'''

from module1 import module1
from module2 import module2
from module3 import module3a, module3b
from module4 import module4
from module5 import module5
from module6 import module6
from sherpa_globals import path_emission_cdf_test, path_area_cdf_test, path_reduction_txt_test, \
    path_model_cdf_test, path_result_cdf_test, path_nuts0_cdf_test, path_nuts2_cdf_test, path_base_conc_cdf_test
from sherpa_auxiliaries import is_number
from sys import argv
import os.path
from time import time



if __name__ == '__main__':
    
    if len(argv) == 1:
        # no command arguments are provided, sherpa is ran in test mode with fixed input arguments
        
        # run module 1 with test inputs
        start = time()
        module1(path_emission_cdf_test, path_area_cdf_test, path_reduction_txt_test, path_model_cdf_test, path_result_cdf_test)
        stop = time()
        print('Module 1 run time: %s sec.' % (stop-start))
        
        
        # run module 2 with test inputs
        start = time()
        module2(path_emission_cdf_test, path_nuts0_cdf_test, path_reduction_txt_test, path_model_cdf_test, path_result_cdf_test)
        stop = time()
        print('Module 2 run time: %s sec.' % (stop-start))
        
        # run module 3a test inputs with test inputs
        start = time()
        module3a(path_emission_cdf_test, path_area_cdf_test, path_reduction_txt_test, path_base_conc_cdf_test, path_model_cdf_test, path_result_cdf_test)
        stop = time()
        print('Module 3a calculation time = %f' % (stop - start))

        # run module 3b with test inputs
        start = time()
        module3b(path_emission_cdf_test, path_area_cdf_test, path_reduction_txt_test, path_base_conc_cdf_test, path_model_cdf_test, path_result_cdf_test)
        stop = time()
        print('Module 3b calculation time = %f' % (stop - start))

        # run module 4  with test inputs
        start = time()
        module4(path_emission_cdf_test, path_area_cdf_test, path_reduction_txt_test, path_base_conc_cdf_test, path_model_cdf_test, path_result_cdf_test)
        stop = time()
        print('Module 4 run time: %s sec.' % (stop-start))
        
        # run module 5 test inputs
        start = time()
        module5(path_emission_cdf_test, path_nuts0_cdf_test, path_reduction_txt_test, path_base_conc_cdf_test, path_model_cdf_test, path_result_cdf_test)
        stop = time()
        print('Module 5 calculation time = %f' % (stop - start))

        # run module 6 test inputs
        start = time()
        # paris
        target_cell_lat = 48.85     # 51.51
        target_cell_lon = 2.35  # 9.19      #-0.13
        module6(path_emission_cdf_test, path_nuts2_cdf_test, target_cell_lat, target_cell_lon, 'input/user_reduction_snap7.txt', path_base_conc_cdf_test, path_model_cdf_test, path_result_cdf_test)
        stop = time()
        print('Module 6 calculation time = %f' % (stop - start))
        
    else:
        # check which module has to be ran
        module = int(argv[1])

        # check if all files (emissions, emission reduction, model) exist
        for i_input in range(2, len(argv)):
            if not(os.path.exists(argv[i_input])) and not(is_number(argv[i_input])):
                print('Error: %s does not exist!' % argv[i_input])
        
        # ---------#
        # module 1 #
        # ---------#
        if module == 1:
            path_emission_cdf = argv[2]     
            path_area_cdf = argv[3]
            path_reduction_txt = argv[4]
            path_base_conc_cdf = argv[5]
            path_model_cdf = argv[6]
            path_result_cdf = argv[7]
            
            # run module 1
            print('running module1')
            module1(path_emission_cdf, path_area_cdf, path_reduction_txt, path_base_conc_cdf, path_model_cdf, path_result_cdf)
            
        # ---------#
        # module 2 #
        # ---------#
        elif module == 2:
            path_emission_cdf = argv[2]     
            path_nuts_cdf = argv[3]
            path_reduction_txt = argv[4]
            path_model_cdf = argv[5]
            path_result_cdf = argv[6]
            
            module2(path_emission_cdf, path_nuts_cdf, path_reduction_txt, path_model_cdf, path_result_cdf)

        # ---------#
        # module 3a #
        # ---------#
        elif module == 31:
            path_emission_cdf = argv[2]     
            path_area_cdf = argv[3]
            path_reduction_txt = argv[4]
            path_base_conc_cdf = argv[5]
            path_model_cdf = argv[6]
            path_result_cdf = argv[7]
            
            module3a(path_emission_cdf, path_area_cdf, path_reduction_txt, path_base_conc_cdf, path_model_cdf, path_result_cdf)

        # ----------#
        # module 3b #
        # ----------#
        elif module == 32:
            path_emission_cdf = argv[2]     
            path_area_cdf = argv[3]
            path_reduction_txt = argv[4]
            path_base_conc_cdf = argv[5]
            path_model_cdf = argv[6]
            path_result_cdf = argv[7]
            
            module3b(path_emission_cdf, path_area_cdf, path_reduction_txt, path_base_conc_cdf, path_model_cdf, path_result_cdf)
            
        # ---------#
        # module 4 #
        # ---------#
        elif module == 4:
            path_emission_cdf = argv[2]     
            path_area_cdf = argv[3]
            path_reduction_txt = argv[4]
            path_base_conc_cdf = argv[5]
            path_model_cdf = argv[6]
            path_result_cdf = argv[7]
            
            module4(path_emission_cdf, path_area_cdf, path_reduction_txt, path_base_conc_cdf, path_model_cdf, path_result_cdf)
            
        # ---------#
        # module 5 #
        # ---------#
        elif module == 5:
            path_emission_cdf = argv[2]     
            path_nuts_cdf = argv[3]
            path_reduction_txt = argv[4]
            path_base_conc_cdf = argv[5]
            path_model_cdf = argv[6]
            path_result_cdf = argv[7]
            
            module5(path_emission_cdf, path_nuts_cdf, path_reduction_txt, path_base_conc_cdf, path_model_cdf, path_result_cdf)

        # ---------#
        # module 6 #
        # ---------#
        elif module == 6:
            path_emission_cdf = argv[2]     
            path_nuts_cdf = argv[3]
            cell_lat = argv[4]
            cell_lon = argv[5]
            path_reduction_txt = argv[6]
            path_base_conc_cdf = argv[7]
            path_model_cdf = argv[8]
            path_result_cdf = argv[9]
            
            module6(path_emission_cdf, path_nuts_cdf, cell_lat, cell_lon, path_reduction_txt, path_base_conc_cdf, path_model_cdf, path_result_cdf)
    
        else:
            print('unknown module %d' % module)
        

    pass



