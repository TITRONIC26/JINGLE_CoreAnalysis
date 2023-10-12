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
    src['LOGMH1'] = src['LOGMH1_MATT']
    src['LOGMH2'] = src['LOGMH2_RYAN']

    src.loc[src['LOGMH1_TING'] == 0, 'LOGMH1_TING'] = np.nan
    src.loc[src['LOGMH2_TING'] == 0, 'LOGMH2'] = np.nan

    src.loc[pd.isnull(src['LOGMH1']), 'LOGMH1'] = src['LOGMH1_TING']
    src.loc[pd.isnull(src['LOGMH2']), 'LOGMH2'] = src['LOGMH2_TING']

    vals = MF.error_add(np.power(10,src['LOGMH1']), np.log(10)*np.power(10,src['LOGMH1'])*(src['LOGMH1_MATT_ERR']), np.power(10,src['LOGMH2']), np.log(10)*np.power(10,src['LOGMH2'])*(src['LOGMH2_RYAN_ERR']))
    src['LOGMGAS'] = np.log10(vals[0])
    src['LOGMGAS_ERR'] = (1/np.log(10)) * ((vals[1])/(vals[0]))
    src['MGAS_FLAG'] = 0

    src.loc[(pd.notnull(src['LOGMH1'])) & (pd.notnull(src['LOGMH2'])), 'MGAS_FLAG'] = 3

    flag = src.index[(src['LOGMGAS'].notnull())].tolist()
    flag_0 = src.index[(src['MGAS_FLAG'] != 3)].tolist()

    G2DR = (src['LOGMGAS'][flag] - src['LOGMDUST_DELOOZE'][flag]).mean()
    G2DR_ERR = (np.sqrt(np.power(src['LOGMGAS_ERR'][flag], 2) + np.power(src['LOGMDUST_DELOOZE_ERR'][flag], 2))).mean()

    y_val = G2DR + src['LOGMDUST_DELOOZE']
    y_val_err = y_val * np.sqrt(mt.pow(G2DR_ERR/G2DR,2) + np.power(src['LOGMDUST_DELOOZE_ERR']/src['LOGMDUST_DELOOZE'],2))

    src.loc[src['MGAS_FLAG'] != 3, 'MGAS_FLAG'] = 4
    src.loc[src['MGAS_FLAG'] == 4, 'LOGMGAS'] = y_val
    src.loc[src['MGAS_FLAG'] == 4, 'LOGMGAS_ERR'] = y_val_err

    src.loc[pd.isnull(src['LOGMGAS']), 'MGAS_FLAG'] = 5
    src.loc[(pd.notnull(src['LOGMH1'])) & (pd.isnull(src['LOGMH2'])), 'MGAS_FLAG'] = 1
    src.loc[(pd.isnull(src['LOGMH1'])) & (pd.notnull(src['LOGMH2'])), 'MGAS_FLAG'] = 2
    src.loc[(pd.isnull(src['LOGMH1_MATT'])) & (pd.isnull(src['LOGMH2_RYAN'])) & (src['MGAS_FLAG'] == 3), 'MGAS_FLAG'] = 6

    h2s = np.log10(np.power(10,src['LOGMGAS']) - np.power(10,src['LOGMH1']))
    h1s = np.log10(np.power(10,src['LOGMGAS']) - np.power(10,src['LOGMH2']))

    src.loc[pd.isnull(src['LOGMH1']), 'LOGMH1'] = src['LOGMH1'] + h1s
    src.loc[pd.isnull(src['LOGMH2']), 'LOGMH2'] = src['LOGMH2'] + h2s

    #MGAS_FLAGS
    #1 == H1 detected
    #2 == H2 detected
    #3 == Both H1 and H2 detected by Ryan and Matt
    #4 == Dust Derived
    #5 == Non detection
    #6 == Detected by Ting

    return src

def find_Mmetals(src):
    fz = 27.36*np.float_power(10, (src['Z_PP04_O3N2'] - 12))

    src = find_Mgas(src)

    src['LOGMZ'] = fz * src['LOGMGAS'] + src['LOGMDUST_DELOOZE']
    src['LOGMZ_ERR'] = MF.error_add(fz*src['LOGMGAS'], fz*src['LOGMGAS_ERR'], src['LOGMDUST_DELOOZE'], src['LOGMDUST_DELOOZE_ERR'])[1]

    return src
