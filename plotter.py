"""
New plotting script for plotting the data gathered from the 4 new datafiles
"""

#libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math as mt

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

def GalacticDensity_field_plots(src):
    xcg = CF.src6[['LOGMSTAR','LOGMH2','LOGMH2_LIM','LOGMH1','LOGSFR']].copy()
    her = CF.src7[['LOGMSTAR','LOGMSTAR_ERR','LOGSFR','LOGSFR_ERR','LOGMH2','LOGMH2_ERR','LOGMH1','LOGMH1_ERR']].copy()

    xcg['LOGMGAS'] = np.log10(np.power(10,xcg['LOGMH1'])+np.power(10,xcg['LOGMH2']))
    xcg['LOGMGAS_LIM'] = np.log10(np.power(10,xcg['LOGMH1'])+np.power(10,xcg['LOGMH2_LIM']))


    her_ind = her.index[(her['MGAS/MSTAR'].notnull()) & her['LOGMSTAR'].notnull()].tolist()

    flag_1 = src.index[(src['MGAS_FLAG'] == 1) & (src['H1_FLAG'] == 1)].tolist()
    flag_1_up = src.index[(src['MGAS_FLAG'] == 1) & (src['H1_FLAG'] == 0)].tolist()

    flag_2 = src.index[src['LOGMH2_RYAN'].notnull()].tolist()

    flag_3 = src.index[(src['MGAS_FLAG'] == 3) & (src['H1_FLAG'] == 1)].tolist()
    flag_3_up = src.index[(src['MGAS_FLAG'] == 3) & (src['H1_FLAG'] == 0)].tolist()

    xcg_ind = xcg.index[(xcg['LOGMH2'].notnull()) & (xcg['LOGMSTAR'].notnull()) & (xcg['LOGMH1'].notnull())].tolist()
    xcg_ind_up = xcg.index[(xcg['LOGMH2_LIM'].notnull()) & (xcg['LOGMSTAR'].notnull()) & (xcg['LOGMH1'].notnull())].tolist()

    FF.print_full(src)
    FF.print_full(xcg.loc)
    FF.print_full(her)

    fig, ax = plt.subplots()

    y_vals = xcg['LOGMGAS'] - xcg['LOGMSTAR']
    ax.errorbar(xcg['LOGMSTAR'][xcg_ind], y_vals[xcg_ind], fmt='s', c='Grey', label='xColdGass '+str(len(xcg_ind)+len(xcg_ind_up)), alpha=0.75)
    y_vals = xcg['LOGMGAS_LIM'] - xcg['LOGMSTAR']
    ax.errorbar(xcg['LOGMSTAR'][xcg_ind_up], y_vals[xcg_ind_up], fmt='s', mfc='White', c='Grey', alpha=0.75)
    #LM.curvefitting(xcg['LOGMSTAR'][xcg_ind], y_vals[xcg_ind], axs=ax, color='Black')

    ax.errorbar(her['LOGMSTAR'][her_ind], her['MGAS/MSTAR'][her_ind], fmt='D', c='Red', label='HERACLES '+str(len(her_ind)), alpha=0.75)
    LM.curvefitting(her['LOGMSTAR'][her_ind], her['MGAS/MSTAR'][her_ind], axs=ax, color='darkred')

    vals = MF.error_subtract(src['LOGMGAS'], src['LOGMGAS_ERR'], src['LOGMSTAR_MAGPHYS'], src['LOGMSTAR_MAGPHYS_ERR'])
    y_vals = vals[0]
    y_err = vals[1]
    ax.errorbar(src['LOGMSTAR_MAGPHYS'][flag_3], y_vals[flag_3], y_err[flag_3], src['LOGMSTAR_MAGPHYS_ERR'][flag_3], fmt='o', c='Blue', ecolor='Blue', label='JINGLE '+str(len(flag_3)+len(flag_3_up)), alpha=0.75)
    ax.errorbar(src['LOGMSTAR_MAGPHYS'][flag_3_up], y_vals[flag_3_up], y_err[flag_3_up], src['LOGMSTAR_MAGPHYS_ERR'][flag_3_up], fmt='o', mfc='White', c='Blue', ecolor='Blue', alpha=0.75)
    LM.curvefitting(src['LOGMSTAR_MAGPHYS'][flag_3], y_vals[flag_3], axs=ax, color='Blue')

    ax.set_ylabel('$M_{Gas}$/$M_{star}$ [Log ()]')
    ax.set_xlabel('$M_{star}}$ [Log ($M_{\odot}$)]')

    ax.set_xlim(8.75, 11.5)
    ax.set_ylim(-1.5, 0.75)

    #formatting plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.show()

def JINGLE_total_gas(vertico=True, xcoldgass=True, heracles=True):
    jngl = CF.src1[['LOGMSTAR_MAGPHYS','LOGMSTAR_MAGPHYS_ERR','LOGMH2_RYAN','LOGMH1_MATT','LOGMH1_MATT_ERR','H1_FLAG','LOGMH2_RYAN_ERR','LOGMGAS','LOGMGAS_ERR','MGAS_FLAG','LOGMDUST_DELOOZE','LOGMDUST_DELOOZE_ERR']].copy()

    flag_3 = jngl.index[(jngl['MGAS_FLAG'] == 3) & (jngl['H1_FLAG'] == 1)].tolist()
    flag_3_up = jngl.index[(jngl['MGAS_FLAG'] == 3) & (jngl['H1_FLAG'] == 0)].tolist()

    flag = jngl.index[(jngl['LOGMGAS'].notnull())].tolist()

    flag_0 = jngl.index[(jngl['MGAS_FLAG'] != 3)].tolist()

    fig,ax = plt.subplots()

    ax.axline((0,0), slope=1, c='Black', alpha=0.1)

    ax.errorbar(jngl['LOGMDUST_DELOOZE'][flag_3], jngl['LOGMGAS'][flag_3], jngl['LOGMGAS_ERR'][flag_3], jngl['LOGMDUST_DELOOZE_ERR'][flag_3], fmt='o', c='Blue', ecolor='Blue', label='JINGLE '+str(len(flag_3)+len(flag_3_up)), alpha=0.75)
    ax.errorbar(jngl['LOGMDUST_DELOOZE'][flag_3_up], jngl['LOGMGAS'][flag_3_up], jngl['LOGMGAS_ERR'][flag_3_up], jngl['LOGMDUST_DELOOZE_ERR'][flag_3_up], fmt='o', mfc='White', c='Blue', ecolor='Blue', alpha=0.75)
    LM.linmixing(jngl['LOGMDUST_DELOOZE'][flag], jngl['LOGMGAS'][flag], jngl['LOGMDUST_DELOOZE_ERR'][flag], jngl['LOGMGAS_ERR'][flag], axs=ax, color='Blue')
    LM.curvefitting(jngl['LOGMDUST_DELOOZE'][flag], jngl['LOGMGAS'][flag], axs=ax)
    LM.pearson(jngl['LOGMDUST_DELOOZE'][flag], jngl['LOGMGAS'][flag], axs=ax)

    ax.set_ylabel('$M_{Gas}$ [Log ($M_{\odot}$)]')      
    ax.set_xlabel('$M_{Dust}}$ [Log ($M_{\odot}$)]')

    ax.set_xlim(6, 9)
    ax.set_ylim(8.5, 11)

    #formatting plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.show()
    
    G2DR = (jngl['LOGMGAS'][flag] - jngl['LOGMDUST_DELOOZE'][flag]).mean()
    G2DR_ERR = (np.sqrt(np.power(jngl['LOGMGAS_ERR'][flag], 2) + np.power(jngl['LOGMDUST_DELOOZE_ERR'][flag], 2))).mean()

    fig,ax = plt.subplots()

    ax.axline((0,0), slope=1, c='Black', alpha=0.1)

    y_val = G2DR + jngl['LOGMDUST_DELOOZE'][flag_0]
    y_val_err = y_val * np.sqrt(mt.pow(G2DR_ERR/G2DR,2) + np.power(jngl['LOGMDUST_DELOOZE_ERR'][flag_0]/jngl['LOGMDUST_DELOOZE'][flag_0],2))
    
    ax.errorbar(jngl['LOGMSTAR_MAGPHYS'][flag_0], y_val, y_val_err, jngl['LOGMSTAR_MAGPHYS_ERR'][flag_0], fmt='D', mfc='White', c='blue', ecolor='blue', label='JINGLE '+str(len(flag_0)), alpha=0.5)

    ax.errorbar(jngl['LOGMSTAR_MAGPHYS'][flag_3], jngl['LOGMGAS'][flag_3], jngl['LOGMGAS_ERR'][flag_3], jngl['LOGMSTAR_MAGPHYS_ERR'][flag_3], fmt='o', c='Blue', ecolor='Blue', label='JINGLE '+str(len(flag_3)+len(flag_3_up)), alpha=0.75)
    ax.errorbar(jngl['LOGMSTAR_MAGPHYS'][flag_3_up], jngl['LOGMGAS'][flag_3_up], jngl['LOGMGAS_ERR'][flag_3_up], jngl['LOGMSTAR_MAGPHYS_ERR'][flag_3_up], fmt='o', mfc='White', c='Blue', ecolor='Blue', alpha=0.75)
    
    if vertico == True:
        vrtc = CF.src5[['LOGMH1','LOGMSTAR','LOGMSTAR_ERR','LOGMH2','LOGMH2_ERR']].copy()
        vrtc = DM.VERTICO_main(vrtc)

        vals = MF.error_add_logs(vrtc['LOGMH2'], vrtc['LOGMH2_ERR'], vrtc['LOGMH1'], vrtc['LOGMH1_ERR'])

        y_val = vals[0]
        y_val_err = vals[1]

        ax.errorbar(vrtc['LOGMSTAR'], y_val, y_val_err, vrtc['LOGMSTAR_ERR'], fmt='o', c='Red', label='VERTICO '+str(len(y_val)), alpha=0.75)
    
    if xcoldgass == True:
        xcld = CF.src6[['LOGMSTAR','LOGMH2','LOGMH2_ERR','LOGMH2_LIM','LOGMH1','LOGMH1_ERR']].copy()

        xcld_flag = xcld.index[(xcld['LOGMH2'].notna()) & (xcld['LOGMSTAR'].notna()) & (xcld['LOGMH1'].notna())].tolist()
        xcld_flag_up = xcld.index[(xcld['LOGMH2_LIM'].notna()) & (xcld['LOGMSTAR'].notna()) & (xcld['LOGMH1'].notna())].tolist()

        xcld.loc[xcld['LOGMH2'].isna(), 'LOGMH2'] = xcld['LOGMH2_LIM']

        vals = MF.error_add_logs(xcld['LOGMH2'], xcld['LOGMH2_ERR'], xcld['LOGMH1'], xcld['LOGMH1_ERR'])

        y_val = vals[0]
        y_val_err = vals[1]

        ax.errorbar(xcld['LOGMSTAR'][xcld_flag], y_val[xcld_flag], y_val_err[xcld_flag], fmt='s', c='Grey', label='XCOLDGASS '+str(len(xcld_flag)+len(xcld_flag_up)), alpha=0.75)
        ax.errorbar(xcld['LOGMSTAR'][xcld_flag_up], y_val[xcld_flag_up], y_val_err[xcld_flag_up], fmt='s', mfc='White', c='Grey', alpha=0.75)

    if heracles == True:
        hrcl = CF.src7[['LOGMSTAR','LOGMH1','LOGMH2']].copy()

        y_val = np.log10((np.power(10,hrcl['LOGMH1']) + np.power(10,hrcl['LOGMH2'])))

        ax.errorbar(hrcl['LOGMSTAR'], y_val, fmt='h', c='Green', label='HERACLES '+str(len(y_val)), alpha=0.75)

    ax.set_ylabel('$M_{Gas}$ [Log ($M_{\odot}$)]')      
    ax.set_xlabel('$M_{Star}}$ [Log ($M_{\odot}$)]')

    ax.set_xlim(8, 11.5)
    ax.set_ylim(7, 11)

    #formatting plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.show()


#call on the main function when the script is executed
if __name__ == '__main__':
    JINGLE_total_gas()