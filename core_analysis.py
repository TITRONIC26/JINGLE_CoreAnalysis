"""
This script is intended for functions that will modify the original values and perform calculations to be used in later sections.
"""

#import libraries here
import pandas as pd
import numpy as np
import math as mt

#import scripts here
import constants as C
import math_funcs as MF

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

def Group_By_Env(src):
    src['GALACTIC_ENV'] = pd.DataFrame(np.zeros(len(src['IDNUM'])))
    
    #"""
    for x in src['IDNUM']:
        idx = src.index[src['IDNUM'] == x][0]

        if (src['LOGMHALO_TEMPEL'][idx] >= C.HALOMASS_HIGH) and (src['NGAL'][idx] >= C.NGAL_HIGH):
            src['GALACTIC_ENV'][idx] = 'Cluster'
        elif (src['LOGMHALO_TEMPEL'][idx] >= C.HALOMASS_LOW) and (src['NGAL'][idx] >= C.NGAL_LOW):
            src['GALACTIC_ENV'][idx] = 'Group'
        elif (src['LOGMHALO_TEMPEL'][idx] >= C.HALOMASS_LOW):
            src['GALACTIC_ENV'][idx] = 'Small Group'
        elif (src['NGAL'][idx] == 2):
            src['GALACTIC_ENV'][idx] = 'Binary Pair'
        elif (src['NGAL'][idx] == 1):
            src['GALACTIC_ENV'][idx] = 'Isolated'
        else:
            src['GALACTIC_ENV'][idx] = 'Poor'
    #"""

    return src

def Group_By_Dens(src):
    src['GALACTIC_DENS'] = pd.DataFrame(np.zeros(len(src['IDNUM'])))
    
    #"""
    for x in src['IDNUM']:
        idx = src.index[src['IDNUM'] == x][0]

        if (src['LOGMHALO_TEMPEL'][idx] >= C.HALOMASS_LOW) and (src['NGAL'][idx] >= C.NGAL_LOW):
            src['GALACTIC_DENS'][idx] = 'High-Density'
        #elif (src['NGAL'][idx] <= 2):
        #    src['GALACTIC_DENS'][idx] = 'Low-Density'
        else:
            src['GALACTIC_DENS'][idx] = 'Medium-Density'
    #"""
    
    return src

def generate_std_error(src, column):
    df = column.to_frame()
    name = str(df.columns.values[0])

    stdev = src.describe().loc['std', name]
    sterr = stdev/mt.sqrt(len(src[name]))

    src[name+str('_ERR')] = pd.DataFrame(np.zeros(len(src[name])))
    src[name+str('_ERR')] = src[name+str('_ERR')].mask(src[name+str('_ERR')] == 0, sterr)

    return src

def find_Mgas(src):
    vals = MF.error_add(np.power(10,src['LOGMH1_MATT']), np.log(10)*np.power(10,src['LOGMH1_MATT'])*(src['LOGMH1_MATT_ERR']), np.power(10,src['LOGMH2_RYAN']), np.log(10)*np.power(10,src['LOGMH2_RYAN'])*(src['LOGMH2_RYAN_ERR']))
    src['LOGMGAS'] = np.log10(vals[0])
    src['LOGMGAS_ERR'] = (1/np.log(10)) * ((vals[1])/(vals[0]))
    src['MGAS_FLAG'] = 0

    src.loc[(pd.notnull(src['LOGMH1_MATT'])) & (pd.isnull(src['LOGMH2_RYAN'])), 'MGAS_FLAG'] = 1
    src.loc[(pd.isnull(src['LOGMH1_MATT'])) & (pd.notnull(src['LOGMH2_RYAN'])), 'MGAS_FLAG'] = 2
    src.loc[(pd.notnull(src['LOGMH1_MATT'])) & (pd.notnull(src['LOGMH2_RYAN'])), 'MGAS_FLAG'] = 3

    return src

def find_Mmetals(src):
    fz = 27.36*np.float_power(10, (src['Z_PP04_O3N2'] - 12))

    src = find_Mgas(src)

    src['LOGMZ'] = fz * src['LOGMGAS'] + src['LOGMDUST_DELOOZE']
    src['LOGMZ_ERR'] = MF.error_add(fz*src['LOGMGAS'], fz*src['LOGMGAS_ERR'], src['LOGMDUST_DELOOZE'], src['LOGMDUST_DELOOZE_ERR'])[1]

    return src
