"""
This script is intended for functions that will modify the original values and perform calculations to be used in later sections.
"""

#import libraries here
import pandas as pd
import numpy as np
import matplotlib as mpl
import scipy as sci
import math as mt

#import scripts here
import constants as C

#global variables here


#function definitions here
def Reference_Scaling(src, z_type):
    z = src[z_type].copy()

    src['LOG_GDR'] = [C.REF_ALPHA + (C.REF_METALICITY - x) for x in z]

    return src['LOG_GDR']

def Calc_Gas_Content_From_Dust(src, dust_type, z_type):
    d = dust_type.copy()

    gdr = Reference_Scaling(src, z_type)

    src['LOGMGAS'] = d + gdr

    return src['LOGMGAS']

def Calc_Gas_Content_Total(src):
    h1 = src['LOGMH1'].copy()
    h2 = src['LOGMH2'].copy()
    h2_upperLimit = src['LOGMH2_PRED'].copy()

    gasTotal = np.log10(np.float_power(10,h1) + np.float_power(10,h2))
    gasTotal_upperLimit = np.log10(np.float_power(10,h1) + np.float_power(10,h2_upperLimit))

    return (gasTotal, gasTotal_upperLimit)

def Group_By_Density(src):
    src['GALACTIC_ENV'] = pd.DataFrame(np.zeros(len(src['IDNUM'])))
    
    #"""
    for x in src['IDNUM']:
        idx = src.index[src['IDNUM'] == x][0]

        if (src['LOGMHALO'][idx] >= C.HALOMASS_HIGH) and (src['NGAL'][idx] >= C.NGAL_HIGH):
            src['GALACTIC_ENV'][idx] = 'Cluster'
        elif (src['LOGMHALO'][idx] >= C.HALOMASS_LOW) and (src['NGAL'][idx] >= C.NGAL_LOW):
            src['GALACTIC_ENV'][idx] = 'Group'
        elif (src['LOGMHALO'][idx] >= C.HALOMASS_LOW):
            src['GALACTIC_ENV'][idx] = 'Small Group'
        elif (src['NGAL'][idx] == 2):
            src['GALACTIC_ENV'][idx] = 'Binary Pair'
        elif (src['NGAL'][idx] == 1):
            src['GALACTIC_ENV'][idx] = 'Isolated'
        else:
            src['GALACTIC_ENV'][idx] = 'Poor'
    #"""

    """
    src.loc[(src['LOGMHALO'] >= C.HALOMASS_HIGH) and (src['NGAL'] >= C.NGAL_HIGH), 'GALACTIC_ENV'] = 'Cluster'
    src.loc[(src['LOGMHALO'] >= C.HALOMASS_LOW) and (src['NGAL'] >= C.NGAL_LOW) and (src['GALACTIC_ENV'] == 0), 'GALACTIC_ENV'] = 'Group'
    src.loc[(src['LOGMHALO'] >= C.HALOMASS_LOW) and (src['GALACTIC_ENV'] == 0), 'GALACTIC_ENV'] = 'Small Group'
    src.loc[(src['NGAL'] == 2) and (src['GALACTIC_ENV'] == 0), 'GALACTIC_ENV'] = 'Binary Pair'
    src.loc[(src['NGAL'] == 1) and (src['GALACTIC_ENV'] == 0), 'GALACTIC_ENV'] = 'Isolated'
    src.loc[(src['GALACTIC_ENV'] == 0), 'GALACTIC_ENV'] = 'Poor'
    """

    return src