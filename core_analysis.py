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