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
    #fixing upper limit duplicates
    src.loc[src['LOGMH2'] != 0, 'LOGMH2_PRED'] = 0

    #generate average errors for DeLooze Dust Masses
    src['LOGMDUST_DELOOZE_ERR'] = [(abs(g)+abs(h))/2 for g,h in zip(src['-dex'], src['+dex'])]
    src['LOGMH2_ALL'] = src['LOGMH2'] + src['LOGMH2_PRED']

    src = CA.generate_std_error(src, src['LOGMH1'])
    src = CA.generate_std_error(src, src['LOGMH2_ALL'])

    src.loc[src['LOGMH1'] <= 0.5, 'LOGMH1'] = np.nan
    src.loc[src['LOGMH2_ALL'] <= 0.5, 'LOGMH2_ALL'] = np.nan

    #generate the total gas masses and metals
    src = CA.find_Mmetals(src)

    return src

def XCOLDGASS_main(src):
    #generate average errors for stellar masses
    src = CA.generate_std_error(src, src['LOGMSTAR'])

    src['LOGMH2_ALL'] = src['LOGMH2'] + src['LOGMH2_PRED']

    src = CA.generate_std_error(src, src['LOGMH2_ALL'])

    src.loc[src['LOGMH2'] <= 0.5, 'LOGMH2'] = np.nan
    src.loc[src['LOGMH2_PRED'] <= 0.5, 'LOGMH2_PRED'] = np.nan
    src.loc[src['LOGMH2_ERR'] > 0.01, 'LOGMH2_ALL_ERR'] = src['LOGMH2_ERR']

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