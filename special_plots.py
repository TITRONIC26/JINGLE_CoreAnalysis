"""
Script meant for showing special case plots. 
"""
#import libraries here
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

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
df_ryan = GD.get_data(GD.RYAN_ISM, keep=True)
df_delooze = GD.get_data(GD.DELOOZE, keep=True)

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

#call on the main function when the script is executed
if __name__ == '__main__':
    ryan_delooze_dustMasses_comparison()
