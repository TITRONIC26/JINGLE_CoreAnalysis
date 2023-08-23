"""
Script reserved exclusively for plotting data from other sections.
Includes the plots from every section, so be mindful of the length of this script.
"""

#import libraries here
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy as sci
import statistics as stat
import math as mt

#import other scripts and constants here
import gather_data as GD
import constants as C

#global variables here

#begin functions for plotting here
def Mstar_vs_Mdust(src):
    fig,ax = plt.subplots()

    #Mstar vs DeLooze
    ax.scatter(src['LOGMSTAR_MAGPHYS'], src['LOGMDUST_DELOOZE'], marker='*', s=C.SIZE, alpha=C.ALPHA, color='Purple', label='DeLooze, 2020')
    #Mstar vs MAGPHYS
    ax.scatter(src['LOGMSTAR_MAGPHYS'], src['LOGMDUST_MAGPHYS'], marker='o', s=C.SIZE, alpha=C.ALPHA, color='Blue', label='Ryan, 2018')
    #Mstar vs SMBB
    ax.scatter(src['LOGMSTAR_MAGPHYS'], src['logMc_SMBB'], marker='s', s=C.SIZE, alpha=C.ALPHA, color='Green', label='SMBB')
    #Mstar vs BMBB
    ax.scatter(src['LOGMSTAR_MAGPHYS'], src['logMc_BMBB'], marker='p', s=C.SIZE, alpha=C.ALPHA, color='Orange', label='BMBB')
    #Mstar vs TMBB
    ax.scatter(src['LOGMSTAR_MAGPHYS'], src['logMc_TMBB'], marker='h', s=C.SIZE, alpha=C.ALPHA, color='Red', label='TMBB')

    #plot labels
    ax.set_title('Evaluated Dust Masses')
    ax.set_xlabel('$M_{star}$ [Log($M_{\odot}$)]')
    ax.set_ylabel('$M_{dust}$ [Log($M_{\odot}$)]')

    #formatting plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    plt.show()

    return

