# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 17:00:37 2016

Measures and technologies database, 
I have to find a better way to do this. 

@author: peduzem
"""
# create abstract methods and objects
from abc import ABCMeta, abstractmethod
import numpy as np

class Vehicle(object):
    """Vehicle metaclass
        Attributes:
            sector: the type of the vehicle is a string
            emf_voc
            emf_nox
            emf_ppm
            emf_ppm_nex
            emf_sox
            emf_co2_wtw
    """

    __metaclass__ = ABCMeta

    base_sale_price = 0
    wheels = 0

    def __init__(self, sector, vtype, activity, network, emf_voc, emf_nox, emf_ppm,
                 emf_ppm_nex, emf_sox, emf_co2_wtw, eff, occ):
        self.sector = sector
        self.vtype = vtype
        self.activity = activity
        self.network = network
        self.emf_voc = emf_voc  # ton/PJ
        self.emf_nox = emf_nox
        self.emf_ppm = emf_ppm
        self.emf_ppm_nex = emf_ppm_nex
        self.emf_sox = emf_sox
        self.emf_co2_wtw = emf_co2_wtw
        self.eff = eff  # Mvkm/PJ
        self.occ = occ  # p/v

    def service_eff(self):
        """Return the service efficiency."""
        return self.eff * self.occ  # Mpkm/PJ       
   
    @abstractmethod
    def vehicle_type(self):
        """"Return a string representing the type of vehicle this is."""
        pass


class PC(Vehicle):
    """A car for sale by Jeffco Car Dealership."""

#    base_sale_price = 8000
#    wheels = 4
#    max_occ = 5
#    av_occ = 1.65

    def vehicle_type(self):
        """"Return a string representing the type of vehicle this is."""
        return 'passenger car'


class PT(Vehicle):
    """Public Transport"""
#
#    base_sale_price = 10000
#    wheels = 4

    def vehicle_type(self):
        """"Return a string representing the type of vehicle this is."""
        return 'public transport'

      

# -------------------------------------------------------------------------
# Data specific to technologies/measures (**have to find a better way)
# -------------------------------------------------------------------------

# Vehicle technologies

# car - GSL - Euro 5 - allnet data from TREMOVE (2010)
# p/v	1.669679702
#            ton/PJ	ton/Mpkm	ton/Mvkm
# nh3	      n.a. 	n.a. 	n.a.
# em NMVOC	10.321432	0.018041	0.030122
# em Nox	10.994458	0.019217	0.032086
# em PPM ex	0.345389	0.000604	0.001008
# em PPM n-ex	2.398561	0.004192	0.007000
# em Sox	0.454410	0.000794	0.001326
# em CO2 wtw	82239	144	240
#            Mvkm/PJ	Mpkm/PJ
# eff	     342.6523843	572.1197309

sector = 'TRA_RD_LD4C'
vtype = 'EU5'
activity = 'GSL'
network = 'all'
emf_voc = 10.321432
emf_nox = 10.994458
emf_ppm = 0.345389
emf_ppm_nex = 0.007000  # ton/Mvkm
emf_sox = 0.454410
emf_co2_wtw = 82239
eff = 342.6523843  # Mvkm/PJ
occ = 1.669679702  # p/c

av_gsl_pc_eu5 = PC(sector, vtype, activity, network, emf_voc, emf_nox, emf_ppm,
                   emf_ppm_nex, emf_sox, emf_co2_wtw, eff, occ)
# car - diesel - Euro 6 - allnet data from TREMOVE (2010)
# p/v	1.665896994
#            ton/PJ	ton/Mpkm	ton/Mvkm
# nh3	      n.a. 	n.a. 	n.a.
# em NMVOC	10.198840	0.015669	0.026103
# em Nox	47.344232	0.072737	0.121172
# em PPM ex	1.635748	0.002513	0.004187
# em PPM n-ex	2.735029	0.004202	0.007000
# em Sox	0.472869	0.000726	0.001210
# em CO2 wtw	83847  	129	      215
#            Mvkm/PJ	Mpkm/PJ
# eff	     390.7185221	650.8968115

sector = 'TRA_RD_LD4C'
vtype = 'EU6'
activity = 'MD'
network = 'all'
emf_voc = 10.198840
emf_nox = 47.344232
emf_ppm = 1.635748
emf_ppm_nex = 0.007000  # ton/Mvkm
emf_sox = 0.472869
emf_co2_wtw = 83847
eff = 390.7185221  # Mvkm/PJ
occ = 1.665896994  # p/c

av_dsl_pc_eu6 = PC(sector, vtype, activity, network, emf_voc, emf_nox, emf_ppm,
                   emf_ppm_nex, emf_sox, emf_co2_wtw, eff, occ)

# car - diesel - av - allnet data from TREMOVE (2010)
#    p/v	1.676010185
#            ton/PJ	ton/Mpkm	ton/Mvkm
# nh3	      n.a.   	n.a. 	      n.a.
# em NMVOC	11.215239	0.017207	0.028839
# em Nox	206.725224	0.317164	0.531569
# em PPM ex	11.532097	0.017693	0.029653
# em PPM n-ex	2.722270	0.004177	0.007000
# em Sox	0.472864	0.000725	0.001216
# em CO2 wtw	83847	      129	     216
#
#       Mvkm/PJ	Mpkm/PJ
# eff	388.8959626	651.7935942

sector = 'TRA_RD_LD4C'
vtype = 'av'
activity = 'MD'
network = 'all'
emf_voc = 11.215239
emf_nox = 206.725224
emf_ppm = 11.532097
emf_ppm_nex = 0.007000  # ton/Mvkm
emf_sox = 0.472864
emf_co2_wtw = 83847
eff = 388.8959626  # Mvkm/PJ
occ = 1.676010185  # p/c

av_dsl_pc = PC(sector, vtype, activity, network, emf_voc, emf_nox, emf_ppm,
               emf_ppm_nex, emf_sox, emf_co2_wtw, eff, occ)

# car - petrol - av - allnet data from TREMOVE (2010)
# p/v	1.650005954
#            ton/PJ	ton/Mpkm	ton/Mvkm
# nh3	n.a. 	n.a. 	n.a.
# em NMVOC	89.436955	0.159527	0.263220
# em Nox	99.336201	0.177184	0.292354
# em PPM ex	0.493155	0.000880	0.001451
# em PPM n-ex	2.378448	0.004242	0.007000
# em Sox	0.454409	0.000811	0.001337
# em CO2 wtw	82239	147	242
#
#       Mvkm/PJ	Mpkm/PJ
# eff	339.780127	560.6392326

sector = 'TRA_RD_LD4C'
vtype = 'av'
activity = 'GSL'
network = 'all'
emf_voc = 89.436955
emf_nox = 99.336201
emf_ppm = 0.493155
emf_ppm_nex = 0.007000  # ton/Mvkm
emf_sox = 0.454409
emf_co2_wtw = 82239
eff = 339.780127  # Mvkm/PJ
occ = 1.650005954  # p/c

av_gsl_pc = PC(sector, vtype, activity, network, emf_voc, emf_nox, emf_ppm,
               emf_ppm_nex, emf_sox, emf_co2_wtw, eff, occ)

#av_gsl_pc=pc()


# PC CNG
    # define the CNG car with the class PC

# 



 
# Example for low emission zone
# substituting current diesel cars with low emission
# increasing number of passengers per vehicles
