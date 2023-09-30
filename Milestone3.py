"""
Script meant to plot and analyze the points outlined for Milestone 3 of IndStudy.
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

"""
Plots to be completed.
----------------------
1) Compare the gas contents of JINGLE in multiple different density environments and compare to Vertico, Heracles, and COLD GASS. 
2) Mimic plots represented within the Catinella paper
    a) Gas ratios relative to M* and SFR / sSFR
    b) Repeat above plots with binning (200 galaxies of JINGLE, can bin into numerous different types depending on M*, SFR)
        i) save the trends for each bin into a plot for future reference
    c) look at total gas content along with binning
"""

def MH1():
    jngl = CF.src1.copy()
    #['LOGMSTAR_MAGPHYS','LOGMSTAR_MAGPHYS_ERR','LOGSFR_MAGPHYS','LOGSFR_MAGPHYS_ERR','LOGMDUST_DELOOZE','LOGMH2_RYAN','LOGMH1_MATT',,'H1_FLAG','LOGMGAS','LOGMGAS_ERR','MGAS_FLAG']
    #FF.print_full(jngl)

    xcg = CF.src6.copy()
    #['LOGMSTAR','LOGSFR','LOGSFR_ERR','LOGMH1','H1_FLAG','LOGMH2','LOGMH2_ERR','LOGMH2_LIM','GROUPID','ENV_CODE','NGAL','LOGMH']
    FF.print_full(xcg)

    Ms = jngl['LOGMSTAR_MAGPHYS']
    Ms_e = jngl['LOGMSTAR_MAGPHYS_ERR']

    sfr = jngl['LOGSFR_MAGPHYS']
    sfr_e = jngl['LOGSFR_MAGPHYS_ERR']

    Mh1 = jngl['LOGMH1_MATT']
    Mh1_e = jngl['LOGMH1_MATT_ERR']

    fig, ax = plt.subplots()

    size = len(Mh1.index[jngl['LOGMH1_MATT'].notnull()].tolist())
    ax.errorbar(Ms, Mh1-Ms, markersize=5, fmt='o', c='Grey', ecolor='Grey', alpha=0.5, zorder=10, label=size)

    size1 = len(xcg.index[xcg['LOGMH1'].notnull()].tolist())
    ax.errorbar(xcg['LOGMSTAR'], xcg['LOGMH1']-xcg['LOGMSTAR'], markersize=5, fmt='d', c='blue', ecolor='Grey', alpha=0.5, zorder=15, label=size1)

    ax.set_ylabel('log $M_{HI}/M_{*}$')      
    ax.set_xlabel('log $M_{*}$')

    #ax.set_xlim(8.5, 11.5)
    #ax.set_ylim(8.25, 11)

    #formatting plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.show()

    #------------------------

    fig, ax = plt.subplots()

    size = len(Mh1.index[jngl['LOGMH1_MATT'].notnull()].tolist())
    ax.errorbar(sfr-Ms, Mh1-Ms, markersize=5, fmt='o', c='Grey', ecolor='Grey', alpha=0.5, zorder=10, label=size)

    size1 = len(xcg.index[xcg['LOGMH1'].notnull()].tolist())
    ax.errorbar(xcg['LOGSFR']-xcg['LOGMSTAR'], xcg['LOGMH1']-xcg['LOGMSTAR'], markersize=5, fmt='d', c='blue', ecolor='Grey', alpha=0.5, zorder=15, label=size1)

    ax.set_ylabel('log $M_{HI}/M_{*}$')      
    ax.set_xlabel('log $sSFR$')

    #ax.set_xlim(8.5, 11.5)
    #ax.set_ylim(8.25, 11)

    #formatting plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.show()




if __name__ == '__main__':
    MH1()