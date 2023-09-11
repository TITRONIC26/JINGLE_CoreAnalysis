"""
New plotting script for plotting the data gathered from the 4 new datafiles
"""

#libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#old plotting scripts
import base_plotter as BPLT
import ratio_plots as RPLT

#essential scripts
import gather_data as GD
import callable_functions as CF
import constants as C

#data analysis and manipulation
import core_analysis as CA
import data_manipulation as DM
import linmixer as LM

import math_funcs as MF

#formatting
import formatting_functions as FF

#define functions here
def gas_plots(src, vertico=False):
    #'LOGMSTAR_MAGPHYS','LOGMSTAR_MAGPHYS_ERR','LOGMH2_RYAN','LOGMH1_MATT','LOGMH1_MATT_ERR','H1_FLAG','LOGMH2_RYAN_ERR','LOGMGAS','LOGMGAS_ERR','MGAS_FLAG'

    flag_0 = src.index[src['MGAS_FLAG'] == 0].tolist()
    
    flag_1 = src.index[(src['MGAS_FLAG'] == 1) & (src['H1_FLAG'] == 1)].tolist()
    flag_1_up = src.index[(src['MGAS_FLAG'] == 1) & (src['H1_FLAG'] == 0)].tolist()

    flag_2 = src.index[src['MGAS_FLAG'] == 2].tolist()

    flag_3 = src.index[(src['MGAS_FLAG'] == 3) & (src['H1_FLAG'] == 1)].tolist()
    flag_3_up = src.index[(src['MGAS_FLAG'] == 3) & (src['H1_FLAG'] == 0)].tolist()

    fig,ax = plt.subplots()

    #Both H1 and H2 are known (flag 3)
    ax.errorbar(src['LOGMSTAR_MAGPHYS'][flag_3], src['LOGMGAS'][flag_3], src['LOGMGAS_ERR'][flag_3], src['LOGMSTAR_MAGPHYS_ERR'][flag_3], fmt='.', c='Blue', ecolor='Blue', label='H1 + H2 ('+str(len(flag_3)+len(flag_3_up))+')')
    ax.errorbar(src['LOGMSTAR_MAGPHYS'][flag_3_up], src['LOGMGAS'][flag_3_up], src['LOGMGAS_ERR'][flag_3_up], src['LOGMSTAR_MAGPHYS_ERR'][flag_3_up], fmt='.', c='Blue', ecolor='Blue', uplims=True)

    #H1 is known (flag 1)
    ax.errorbar(src['LOGMSTAR_MAGPHYS'][flag_1], src['LOGMH1_MATT'][flag_1], src['LOGMH1_MATT_ERR'][flag_1], src['LOGMSTAR_MAGPHYS_ERR'][flag_1], fmt='.', c='Green', ecolor='Green', label='H1 only ('+str(len(flag_1)+len(flag_1_up))+')')
    ax.errorbar(src['LOGMSTAR_MAGPHYS'][flag_1_up], src['LOGMH1_MATT'][flag_1_up], src['LOGMH1_MATT_ERR'][flag_1_up], src['LOGMSTAR_MAGPHYS_ERR'][flag_1_up], fmt='.', c='Green', ecolor='Green', uplims=True)

    #H2 is known (flag 2)
    ax.errorbar(src['LOGMSTAR_MAGPHYS'][flag_2], src['LOGMH2_RYAN'][flag_2], src['LOGMH2_RYAN_ERR'][flag_2], src['LOGMSTAR_MAGPHYS_ERR'][flag_2], fmt='.', c='Red', ecolor='Red', label='H2 only ('+str(len(flag_2))+')')

    #No gas content known (flag 0)
    ax.errorbar(src['LOGMSTAR_MAGPHYS'][flag_0], src['LOGMGAS'][flag_0], src['LOGMGAS_ERR'][flag_0], src['LOGMSTAR_MAGPHYS_ERR'][flag_0], fmt='.', c='Grey', ecolor='Grey', lolims=True, label='No Data ('+str(len(flag_0))+')')

    ax.set_ylabel('$M_{Gas}$ [Log ($M_{\odot}$)]')
    ax.set_xlabel('$M_{Star}$ [Log ($M_{\odot}$)]')

    #formatting plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.show()

    fig,ax = plt.subplots()

    drop = 0.5

    if vertico == True:
        ver = CF.src4.copy()

        ver = ver.loc[(ver['LOGMH2'].notnull()) & (ver['LOGMH1_DP'].notnull())]

        flag_ver_H2 = ver.index[ver['L_LOGMH2'].notnull()].tolist()
        flag_ver = ver.index[ver['L_LOGMH2'].isnull()].tolist()

        ax.errorbar(ver['LOGMH2'][flag_ver], ver['LOGMH1_DP'][flag_ver], ver['LOGMH1_DP_ERR'][flag_ver], ver['LOGMH2_ERR'][flag_ver], fmt='s', c='Black', ecolor='Black', alpha=0.2, label='('+str(len(flag_ver)+len(flag_ver_H2))+')')
        ax.errorbar(ver['LOGMH2'][flag_ver_H2], ver['LOGMH1_DP'][flag_ver_H2], ver['LOGMH1_DP_ERR'][flag_ver_H2], ver['LOGMH2_ERR'][flag_ver_H2], fmt='<', c='Black', ecolor='Black', alpha=0.2, xuplims=True)
        LM.curvefitting(ver['LOGMH2'][flag_ver], ver['LOGMH1_DP'][flag_ver], axs=ax, color='Grey')

        drop = 1.5
        
    #Both
    ax.errorbar(src['LOGMH2_RYAN'][flag_3], src['LOGMH1_MATT'][flag_3], src['LOGMH1_MATT_ERR'][flag_3], src['LOGMH2_RYAN_ERR'][flag_3], fmt='.', c='Blue', ecolor='Blue', label='('+str(len(flag_3)+len(flag_3_up))+')')
    ax.errorbar(src['LOGMH2_RYAN'][flag_3_up], src['LOGMH1_MATT'][flag_3_up], src['LOGMH1_MATT_ERR'][flag_3_up], src['LOGMH2_RYAN_ERR'][flag_3_up], fmt='.', c='Blue', ecolor='Blue', uplims=True)
    LM.curvefitting(src['LOGMH2_RYAN'][flag_3], src['LOGMH1_MATT'][flag_3], axs=ax, color='Black')

    #H1
    src['H2'] = src['LOGMH2_RYAN'][flag_3].min() - drop
    ax.errorbar(src['H2'][flag_1], src['LOGMH1_MATT'][flag_1], src['LOGMH1_MATT_ERR'][flag_1], fmt='>', c='Green', ecolor='Green', label='('+str(len(flag_1)+len(flag_1_up))+')')
    ax.errorbar(src['H2'][flag_1_up], src['LOGMH1_MATT'][flag_1_up], src['LOGMH1_MATT_ERR'][flag_1_up], fmt='>', c='Green', ecolor='Green', uplims=True)
                
    #H2
    src['H1'] = src['LOGMH1_MATT'][flag_3].min() - drop
    ax.errorbar(src['LOGMH2_RYAN'][flag_2], src['H1'][flag_2], xerr=src['LOGMH2_RYAN_ERR'][flag_2], fmt='^', c='Red', ecolor='Red', label='('+str(len(flag_2))+')')

    #Neither
    ax.errorbar(src['H2'][flag_0], src['H1'][flag_0], fmt='x', c='Grey', ecolor='Grey', label='('+str(len(flag_0))+')')

    ax.set_ylabel('$M_{HI}$ [Log ($M_{\odot}$)]')
    ax.set_xlabel('$M_{HII}$ [Log ($M_{\odot}$)]')

    #formatting plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.show()

    return

def GalacticDensity_gas_plots(src, sSFR=False):
    #'ID','RA','DE','Vel','S/N','L_LOGMH2','LOGMH2','LOGMH2_ERR','LOGMSTAR','LOGMSTAR_ERR','LOGMH1','LOGMH1_ERR'
    vert = CF.src5.copy()
    verDP = CF.src4[['ID','LOGSFR_DP','LOGSFR_DP_ERR']].copy()

    ver = pd.merge(vert, verDP, on='ID')

    FF.print_full(ver)

    #'LOGMSTAR_MAGPHYS','LOGMSTAR_MAGPHYS_ERR','LOGSFR_MAGPHYS','LOGSFR_MAGPHYS_ERR','LOGMH2_RYAN','LOGMH1_MATT','LOGMH1_MATT_ERR','H1_FLAG','LOGMH2_RYAN_ERR','LOGMGAS','LOGMGAS_ERR','MGAS_FLAG'

    #plot of sfr vs stellar mass
    fig, ax = plt.subplots()

    index_up = src.index[(src['H1_FLAG'] == 0)].tolist()
    index = src.index[(src['H1_FLAG'] == 1)].tolist()

    lenth = str(len(index))

    if sSFR == True:
        vals = MF.error_subtract(src['LOGMH1_MATT'], src['LOGMH1_MATT_ERR'], src['LOGMSTAR_MAGPHYS'], src['LOGMSTAR_MAGPHYS_ERR'])
        values = vals[0]
        value_err = vals[1]
    else:
        values = src['LOGMH1_MATT']
        value_err = src['LOGMH1_MATT_ERR']

    ax.errorbar(src['LOGMSTAR_MAGPHYS'][index], values[index], value_err[index], src['LOGMSTAR_MAGPHYS_ERR'][index], fmt='o', c='Blue', ecolor='Blue', label='JINGLE '+lenth, alpha=0.75)
    ax.errorbar(src['LOGMSTAR_MAGPHYS'][index_up], values[index_up], value_err[index_up], src['LOGMSTAR_MAGPHYS_ERR'][index_up], fmt='o', mfc='white', c='Blue', ecolor='Blue')
    LM.curvefitting(src['LOGMSTAR_MAGPHYS'][index], values[index], axs=ax, color='Blue')

    if sSFR == True:
        ver_vals = MF.error_subtract(ver['LOGMH1'], ver['LOGMH1_ERR'], ver['LOGMSTAR'], ver['LOGMSTAR_ERR'])
        ver_values = ver_vals[0]
        ver_value_err = ver_vals[1]
    else:
        ver_values = ver['LOGMH1']
        ver_value_err = ver['LOGMH1_ERR']

    ax.errorbar(ver['LOGMSTAR'], ver_values, ver_value_err, ver['LOGMSTAR_ERR'], fmt='s', c='Grey', ecolor='Grey', label='VERTICO', alpha=0.75)
    flag = ver.index[(ver['LOGMSTAR'].notnull()) & (ver['LOGMH1'].notnull())].tolist()
    LM.curvefitting(ver['LOGMSTAR'][flag], ver_values[flag], axs=ax, color='Black')

    ax.set_ylabel('$M_{HI}$/$M_{star}$ [Log ()]')
    ax.set_xlabel('$M_{star}}$ [Log ($M_{\odot}$)]')

    ax.set_xlim(8.75, 11.5)
    ax.set_ylim(-3.5, 1)

    #formatting plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.show()

    