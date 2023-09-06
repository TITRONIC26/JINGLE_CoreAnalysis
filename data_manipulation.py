"""
Python script aimed at modifying the data directly from the CSV file
"""

#import libraries here
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy as sci
import math as mt

#import other scripts here
import gather_data as GD
import base_plotter as BPLT
import core_analysis as CA
import formatting_functions as FF
import ratio_plots as RPLT
import constants as C

#declare global variables here

#define functions here
def JINGLE_main(src):
    #replace with null values
    src.loc[src['Z_PP04_O3N2'] < 0, 'Z_PP04_O3N2'] = np.nan
    src.loc[src['LOGMHALO_TEMPEL'] < 0, 'LOGMHALO_TEMPEL'] = np.nan
    src.loc[src['LOGMH2_RYAN'] < 0.1, 'LOGMH2_RYAN'] = np.nan
    src.loc[src['H1_FLAG'] < 0, ['LOGMH1_MATT','LOGMH1_MATT_ERR']] = np.nan

    src = CA.generate_std_error(src, src['LOGMH2_RYAN'])
    
    src['LOGMDUST_DELOOZE_ERR'] = [(abs(g)+abs(h))/2 for g,h in zip(src['dex-'], src['dex+'])]
    
    #generate the total gas masses and metals
    src = CA.find_Mmetals(src)

    return src

def VERTICO_main(src):
    #convert the non log base dustpedia masses to log base for comparing to Brown.
    src['LOGMSTAR_DP_ERR'] = (1/np.log(10)) * (src['LOGMSTAR_DP_ERR'] / src['LOGMSTAR_DP'])
    src['LOGMSTAR_DP'] = np.log10(src['LOGMSTAR_DP'])

    src['LOGSFR_DP_ERR'] = (1/np.log(10)) * (src['LOGSFR_DP_ERR'] / src['LOGSFR_DP'])
    src['LOGSFR_DP'] = np.log10(src['LOGSFR_DP'])

    src['LOGMDUST_DP_ERR'] = (1/np.log(10)) * (src['LOGMDUST_DP_ERR'] / src['LOGMDUST_DP'])
    src['LOGMDUST_DP'] = np.log10(src['LOGMDUST_DP'])

    return src