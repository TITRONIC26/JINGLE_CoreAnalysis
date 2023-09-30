"""
Script meant for showing special case plots. 
"""
#import libraries here
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math as mt

#plotting scripts
import base_plotter as BPLT
import ratio_plots as RPLT

#essential scripts
import gather_data as GD
import constants as C
import callable_functions as CF

#data analysis and manipulation
import core_analysis as CA
import data_manipulation as DM
import linmixer as LM

#formatting
import formatting_functions as FF

#get data frames
#df_ryan = GD.get_data(GD.RYAN_ISM, keep=True)
#df_delooze = GD.get_data(GD.DELOOZE, keep=True)

def ryan_delooze_dustMasses_comparison():
    src1 = df_ryan.copy()
    src2 = df_delooze.copy()

    src = pd.merge(src1, src2, on='IDNUM')

    fig,ax = plt.subplots()

    ax.scatter(src['LOGMSTAR_MAGPHYS'], src['LOGMDUST_MAGPHYS'], s=5, color='Red', label='Ryan')
    LM.curvefitting(src['LOGMSTAR_MAGPHYS'], src['LOGMDUST_MAGPHYS'], axs=ax, color='Black')
    ax.scatter(src['LOGMSTAR_MAGPHYS'], src['LOGMDUST'], s=5, color='Grey', label='DeLooze')
    LM.curvefitting(src['LOGMSTAR_MAGPHYS'], src['LOGMDUST'], axs=ax, color='Blue')

    ax.set_title('Comparison of Ryan and DeLooze Dust Masses')
    ax.set_xlabel('$M_{Star}$ [Log ($M_{\odot}$)]')
    ax.set_ylabel('$M_{Dust}$ [Log ($M_{\odot}$)]')

    ax.legend()
    plt.show()

    separation = src['LOGMDUST_MAGPHYS'] - src['LOGMDUST']

    fig,ax = plt.subplots()

    ax.bar(src['IDNUM'], separation)

    ax.hlines(separation.mean(), xmin=0, xmax=193, linestyles='--', colors='Black', label=str("{0:.3g}".format(separation.mean())))

    ax.set_title('$M_{Dust}$ Separation between Ryan and DeLooze')
    ax.set_xlabel('Galaxy ID Number')
    ax.set_ylabel('$M_{Ryan} / M_{DeLooze}$ [Log ()]')

    ax.legend()
    plt.show()

    return

def slope_fixed(x):
    #y = a + (x0 - x)
    a = 2.21
    x0 = 8.69

    return a + (x0 - x)

def slope_free(x):
    #y = a + A(x0 - x)
    a = 2.21
    A = 2.02
    A_err = 0.28
    x0 = 8.69

    return a + A*(x0 - x)

def broken_power_law(x):
    #y = a + A(x0 - x) for x > xt
    #y = b + B(x0 - x) for x <= xt
    a = 2.21
    A = 1.00
    b = 0.96
    B = 3.10
    B_err = 1.33
    xt = 8.10
    xt_err = 0.43
    x0 = 8.69

    ys = []

    for xs in x:
        if xs <= xt:
            ys.append(b + B*(x0-xs))
        else:
            ys.append(a + A*(x0-xs))

    return ys

def dustMass_comp():
    src = CF.src1[['Z_PP04_O3N2','LOGMSTAR_MAGPHYS','LOGMSTAR_MAGPHYS_ERR','LOGMH2_RYAN','LOGMH2_RYAN_ERR','LOGMH1_MATT','LOGMH1_MATT_ERR','H1_FLAG','LOGMGAS','LOGMGAS_ERR','MGAS_FLAG','LOGMDUST_DELOOZE','LOGMDUST_DELOOZE_ERR']]

    GD1 = slope_fixed(src['Z_PP04_O3N2'])
    GD2 = slope_free(src['Z_PP04_O3N2'])
    GD3 = broken_power_law(src['Z_PP04_O3N2'])

    fixed = GD1 + src['LOGMDUST_DELOOZE']
    free = GD2 + src['LOGMDUST_DELOOZE']
    broken = GD3 + src['LOGMDUST_DELOOZE']

    flag_3 = src.index[(src['MGAS_FLAG'] == 3) & (src['H1_FLAG'] == 1)].tolist()
    flag_3_up = src.index[(src['MGAS_FLAG'] == 3) & (src['H1_FLAG'] == 0)].tolist()

    flag = src.index[(src['LOGMGAS'].notnull())].tolist()
    flag_0 = src.index[(src['MGAS_FLAG'] != 3)].tolist()

    G2DR = (src['LOGMGAS'][flag] - src['LOGMDUST_DELOOZE'][flag]).mean()
    G2DR_ERR = (np.sqrt(np.power(src['LOGMGAS_ERR'][flag], 2) + np.power(src['LOGMDUST_DELOOZE_ERR'][flag], 2))).mean()

    y_val = G2DR + src['LOGMDUST_DELOOZE']
    y_val_err = y_val * np.sqrt(mt.pow(G2DR_ERR/G2DR,2) + np.power(src['LOGMDUST_DELOOZE_ERR']/src['LOGMDUST_DELOOZE'],2))

    fig, ax = plt.subplots()

    ax.scatter(src['LOGMSTAR_MAGPHYS'], fixed, marker='+', s=15, color='Red', label='Reference Scaling', alpha = 0.75, zorder=10)
    ax.scatter(src['LOGMSTAR_MAGPHYS'], free, marker='H', s=15, color='Green', label='Slope Free Power Law', alpha = 0.75, zorder=8)
    ax.scatter(src['LOGMSTAR_MAGPHYS'], broken, marker='^', s=15, color='Blue', label='Broken Power Law', alpha = 0.75, zorder=6)

    ax.errorbar(src['LOGMSTAR_MAGPHYS'][flag_0], y_val[flag_0], y_val_err[flag_0], src['LOGMSTAR_MAGPHYS_ERR'][flag_0], markersize=5, fmt='D', mfc='White', c='Grey', ecolor='Grey', label='From Known G/D '+str(len(flag_0)), alpha=0.5, zorder=1)
    ax.errorbar(src['LOGMSTAR_MAGPHYS'][flag_3], src['LOGMGAS'][flag_3], src['LOGMGAS_ERR'][flag_3], src['LOGMSTAR_MAGPHYS_ERR'][flag_3], markersize=5, fmt='o', c='Grey', ecolor='Grey', label='JINGLE Gas '+str(len(flag_3)+len(flag_3_up)), alpha=0.75, zorder=3)
    ax.errorbar(src['LOGMSTAR_MAGPHYS'][flag_3_up], src['LOGMGAS'][flag_3_up], src['LOGMGAS_ERR'][flag_3_up], src['LOGMSTAR_MAGPHYS_ERR'][flag_3_up], markersize=5, fmt='o', mfc='White', c='Grey', ecolor='Grey', alpha=0.75, zorder=2)

    ax.set_ylabel('$M_{Gas}$ [Log ($M_{\odot}$)]')      
    ax.set_xlabel('$M_{Star}}$ [Log ($M_{\odot}$)]')

    ax.set_xlim(8.5, 11.5)
    ax.set_ylim(8.25, 11)

    #formatting plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.show()

    src.loc[src['LOGMGAS'].isna(), 'LOGMGAS'] = y_val

    fig,ax = plt.subplots(3, sharex=True)

    lst = [fixed, free, broken]
    names = ['fixed','free','broken']
    c = 0

    for n in lst:
        avg = (((src['LOGMGAS'][flag_3]-n[flag_3])).mean() + ((src['LOGMGAS'][flag_3_up]-n[flag_3_up])).mean())/2

        ax[c].scatter(src['LOGMSTAR_MAGPHYS'][flag_3], (src['LOGMGAS'] - n)[flag_3], s=15, alpha = 0.75, marker='o', color='Blue')
        ax[c].errorbar(src['LOGMSTAR_MAGPHYS'][flag_3_up], (src['LOGMGAS'] - n)[flag_3_up], markersize=5, alpha = 0.75, fmt='o', mfc='White', c='Blue')
        ax[c].hlines(avg, xmin = 8.5, xmax = 11.5, linestyle='--', color='Black', label=str("{0:.5g}".format(avg)))
        ax[c].set_ylabel('Log(JINGLE) - Log('+str(names[c])+')')
        ax[c].set_xlim(8.75, 11.35)
        ax[c].legend()
        c+=1

    plt.xlabel('$M_{Star}}$ [Log ($M_{\odot}$)]')
    plt.show()

    

    return

#call on the main function when the script is executed
if __name__ == '__main__':
    #ryan_delooze_dustMasses_comparison()
    dustMass_comp()
    #print()