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

def all():
    jngl = CF.src1.copy()
    xcg = CF.src6.copy()
    vrt = CF.src5.copy()
    hrc = CF.src7.copy()

    Ms_list = [jngl['LOGMSTAR_MAGPHYS'],xcg['LOGMSTAR'],vrt['LOGMSTAR'],hrc['LOGMSTAR']]
    sfr_list = [jngl['LOGSFR_MAGPHYS'],xcg['LOGSFR'],vrt['LOGSFR'],hrc['LOGSFR']]
    Mh1_list = [jngl['LOGMH1_MATT'],xcg['LOGMH1'],vrt['LOGMH1'],hrc['LOGMH1']]

    xcg.loc[xcg['LOGMH2'].isnull(), 'LOGMH2'] = 0
    xcg.loc[xcg['LOGMH2_LIM'].isnull(), 'LOGMH2_LIM'] = 0
    xcg['LOGMH2_ALL'] = xcg['LOGMH2'] + xcg['LOGMH2_LIM']

    Mh2_list = [jngl['LOGMH2_RYAN'],xcg['LOGMH2_ALL'],vrt['LOGMH2'],hrc['LOGMH2']]

    Ms = pd.concat(Ms_list, ignore_index=True)
    sfr = pd.concat(sfr_list, ignore_index=True)
    Mh1 = pd.concat(Mh1_list, ignore_index=True)
    Mh2 = pd.concat(Mh2_list, ignore_index=True)

    all = pd.concat([Ms,sfr,Mh1,Mh2], axis=1)
    all = all.rename(columns={0:'LOGMSTAR',1:'LOGSFR',2:'LOGMH1',3:'LOGMH2'})
    #FF.print_full(all)

    return all

def binning():
    all_data = all()

    #Mstar binning
    range_all = 12-7
    steps = 9
    step = range_all/steps

    bins = np.arange(7,12,step)
    labels = [0,1,2,3,4,5,6,7]

    bin = (pd.cut(all_data['LOGMSTAR'], bins, labels=labels))

    all_data['MSTAR_BIN'] = bin

    #sSFR binning
    range_all = (-8)-(-13)
    steps = 6
    step = range_all/steps

    bins = np.arange(-13,-8,step)
    labels = [0,1,2,3,4]

    bin = (pd.cut(all_data['LOGSFR']-all_data['LOGMSTAR'], bins, labels=labels))
    print(bin)

    all_data['SFR_BIN'] = bin

    return all_data

def weighting(col,y,x):
    all_data = binning()

    groups = all_data.groupby(col)

    ave = []
    ave_e = []
    med = []
    xs = []

    for key, group in groups:
        group = group.dropna()
        print(group)
        ave.append(np.average(group[y]-group[x]))#, weights=group[x]))
        ave_e.append(np.std(group[y]-group[x]) / np.sqrt(len(group[y])))
        med.append(np.median(group[y]-group[x]))
        xs.append(np.average(group[x]))

    return (xs, ave, ave_e, med)

def weighting2(col,y,x1,x2):
    all_data = binning()

    groups = all_data.groupby(col)

    ave = []
    ave_e = []
    med = []
    xs = []

    for key, group in groups:
        group = group.dropna()
        ave.append(np.average(group[y]-group[x1]))#, weights=group[x2]-group[x1]))
        ave_e.append(np.std(group[y]-group[x1]) / np.sqrt(len(group[y])))
        med.append(np.median(group[y]-group[x1]))
        xs.append(np.average(group[x2]-group[x1]))

    return (xs, ave, ave_e, med)

def MH1():
    jngl = CF.src1.copy()
    #['LOGMSTAR_MAGPHYS','LOGMSTAR_MAGPHYS_ERR','LOGSFR_MAGPHYS','LOGSFR_MAGPHYS_ERR','LOGMDUST_DELOOZE','LOGMH2_RYAN','LOGMH1_MATT',,'H1_FLAG','LOGMGAS','LOGMGAS_ERR','MGAS_FLAG']
    #FF.print_full(jngl)
    j = jngl.index[jngl['H1_FLAG']==1].tolist()
    jup = jngl.index[jngl['H1_FLAG']==0].tolist()

    xcg = CF.src6.copy()
    #['LOGMSTAR','LOGSFR','LOGSFR_ERR','LOGMH1','H1_FLAG','LOGMH2','LOGMH2_ERR','LOGMH2_LIM','GROUPID','ENV_CODE','NGAL','LOGMH']
    #FF.print_full(xcg)
    x = xcg.index[xcg['H1_FLAG']!=99].tolist()
    xup = xcg.index[xcg['H1_FLAG']==99].tolist()

    vrt = CF.src5.copy()
    #['ID','RA','DEC','z','LOGMH1','LOGMSTAR','LOGMSTAR_ERR','LOGSFR','LOGSFR_ERR','LOGMH2','LOGMH2_ERR','LOGMH1_ERR]
    #FF.print_full(vrt)

    hrc = CF.src7.copy()
    #['ID','RA','DEC','Vel','z','LOGMSTAR','LOGMSTAR_ERR','LOGSFR','LOGSFR_ERR','LOGMH2','LOGMH2_ERR','LOGMH1','LOGMH1_ERR']
    #FF.print_full(hrc)

    Ms = jngl['LOGMSTAR_MAGPHYS']
    Ms_e = jngl['LOGMSTAR_MAGPHYS_ERR']

    sfr = jngl['LOGSFR_MAGPHYS']
    sfr_e = jngl['LOGSFR_MAGPHYS_ERR']

    Mh1 = jngl['LOGMH1_MATT']

    #-----------------------------

    fig, ax = plt.subplots()

    size = len(Mh1.index[jngl['LOGMH1_MATT'].notnull()].tolist())
    ax.errorbar(Ms[j], Mh1[j]-Ms[j], markersize=5, fmt='o', c='Grey', ecolor='Grey', alpha=0.5, zorder=15, label=str(size)+' - JINGLE')
    ax.errorbar(Ms[jup], Mh1[jup]-Ms[jup], markersize=5, fmt='v', c='Grey', ecolor='Grey', zorder=15)

    size1 = len(xcg.index[xcg['LOGMH1'].notnull()].tolist())
    ax.errorbar(xcg['LOGMSTAR'][x], xcg['LOGMH1'][x]-xcg['LOGMSTAR'][x], markersize=5, fmt='o', c='grey', ecolor='Grey', alpha=0.5, zorder=10, label=str(size1)+' - xCOLDGASS')
    ax.errorbar(xcg['LOGMSTAR'][xup], xcg['LOGMH1'][xup]-xcg['LOGMSTAR'][xup], markersize=5, fmt='v', c='grey', ecolor='Grey', zorder=10)

    size2 = len(vrt.index[vrt['LOGMH1'].notnull()].tolist())
    ax.errorbar(vrt['LOGMSTAR'], vrt['LOGMH1']-vrt['LOGMSTAR'], markersize=5, fmt='o', c='grey', ecolor='Grey', alpha=0.5, zorder=20, label=str(size2)+' - VERTICO')
    
    size3 = len(hrc.index[hrc['LOGMH1'].notnull()].tolist())
    ax.errorbar(hrc['LOGMSTAR'], hrc['LOGMH1']-hrc['LOGMSTAR'], markersize=5, fmt='o', c='grey', ecolor='Grey', alpha=0.5, zorder=25, label=str(size3)+' - HERACLES')

    ax.set_ylabel('log $M_{HI}/M_{*}$')      
    ax.set_xlabel('log $M_{*}$ $[M_{\odot}]$')

    #binned stuff
    values = weighting('MSTAR_BIN','LOGMH1','LOGMSTAR')
    ax.errorbar(values[0], values[1], values[2], markersize=10, fmt='o', c='black', mfc='red', zorder=100)
    ax.errorbar(values[0], values[3], markersize=8, fmt='D', c='black', mfc='dodgerblue', zorder=100)

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
    ax.errorbar(sfr[j]-Ms[j], Mh1[j]-Ms[j], markersize=5, fmt='o', c='Grey', ecolor='Grey', alpha=0.5, zorder=15, label=str(size)+' - JINGLE')
    ax.errorbar(sfr[jup]-Ms[jup], Mh1[jup]-Ms[jup], markersize=5, fmt='v', c='Grey', ecolor='Grey', zorder=15)

    size1 = len(xcg.index[xcg['LOGMH1'].notnull()].tolist())
    ax.errorbar(xcg['LOGSFR'][x]-xcg['LOGMSTAR'][x], xcg['LOGMH1'][x]-xcg['LOGMSTAR'][x], markersize=5, fmt='o', c='grey', ecolor='Grey', alpha=0.5, zorder=10, label=str(size1)+' - xCOLDGASS')
    ax.errorbar(xcg['LOGSFR'][xup]-xcg['LOGMSTAR'][xup], xcg['LOGMH1'][xup]-xcg['LOGMSTAR'][xup], markersize=5, fmt='v', c='grey', ecolor='Grey', zorder=10)

    size2 = len(vrt.index[vrt['LOGMH1'].notnull()].tolist())
    ax.errorbar(vrt['LOGSFR']-vrt['LOGMSTAR'], vrt['LOGMH1']-vrt['LOGMSTAR'], markersize=5, fmt='o', c='grey', ecolor='Grey', alpha=0.5, zorder=20, label=str(size2)+' - VERTICO')
    
    size3 = len(hrc.index[hrc['LOGMH1'].notnull()].tolist())
    ax.errorbar(hrc['LOGSFR']-hrc['LOGMSTAR'], hrc['LOGMH1']-hrc['LOGMSTAR'], markersize=5, fmt='o', c='grey', ecolor='Grey', alpha=0.5, zorder=25, label=str(size3)+' - HERACLES')

    ax.set_ylabel('log $M_{HI}/M_{*}$')      
    ax.set_xlabel('log $sSFR$  $[yr^{-1}]$')

    #binned stuff
    values = weighting2('SFR_BIN','LOGMH1','LOGMSTAR','LOGSFR')
    ax.errorbar(values[0], values[1], values[2], markersize=10, fmt='o', c='black', mfc='red', zorder=100)
    ax.errorbar(values[0], values[3], markersize=8, fmt='D', c='black', mfc='dodgerblue', zorder=100)

    #ax.set_xlim(8.5, 11.5)
    #ax.set_ylim(8.25, 11)

    #formatting plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.show()

def MH2():
    jngl = CF.src1.copy()
    #['LOGMSTAR_MAGPHYS','LOGMSTAR_MAGPHYS_ERR','LOGSFR_MAGPHYS','LOGSFR_MAGPHYS_ERR','LOGMDUST_DELOOZE','LOGMH2_RYAN','LOGMH1_MATT',,'H1_FLAG','LOGMGAS','LOGMGAS_ERR','MGAS_FLAG']
    #FF.print_full(jngl)

    xcg = CF.src6.copy()
    #['LOGMSTAR','LOGSFR','LOGSFR_ERR','LOGMH1','H1_FLAG','LOGMH2','LOGMH2_ERR','LOGMH2_LIM','GROUPID','ENV_CODE','NGAL','LOGMH']
    #FF.print_full(xcg)
    x = xcg.index[xcg['LOGMH2'].notnull()].tolist()
    xup = xcg.index[xcg['LOGMH2_LIM'].notnull()].tolist()

    vrt = CF.src5.copy()
    #['ID','RA','DEC','z','LOGMH1','LOGMSTAR','LOGMSTAR_ERR','LOGSFR','LOGSFR_ERR','LOGMH2','LOGMH2_ERR','LOGMH1_ERR]
    #FF.print_full(vrt)

    hrc = CF.src7.copy()
    #['ID','RA','DEC','Vel','z','LOGMSTAR','LOGMSTAR_ERR','LOGSFR','LOGSFR_ERR','LOGMH2','LOGMH2_ERR','LOGMH1','LOGMH1_ERR']
    #FF.print_full(hrc)

    Ms = jngl['LOGMSTAR_MAGPHYS']
    Ms_e = jngl['LOGMSTAR_MAGPHYS_ERR']

    sfr = jngl['LOGSFR_MAGPHYS']
    sfr_e = jngl['LOGSFR_MAGPHYS_ERR']

    Mh2 = jngl['LOGMH2_RYAN']

    #-----------------------------

    fig, ax = plt.subplots()

    size = len(Mh2.index[jngl['LOGMH2_RYAN'].notnull()].tolist())
    ax.errorbar(Ms, Mh2-Ms, markersize=5, fmt='o', c='Grey', ecolor='Grey', alpha=0.5, zorder=15, label=str(size)+' - JINGLE')
    
    size1 = len(xcg.index[xcg['LOGMH2'].notnull()].tolist()) + len(xcg.index[xcg['LOGMH2_LIM'].notnull()].tolist())
    ax.errorbar(xcg['LOGMSTAR'][x], xcg['LOGMH2'][x]-xcg['LOGMSTAR'][x], markersize=5, fmt='o', c='Grey', ecolor='Grey', alpha=0.5, zorder=10, label=str(size1)+' - xCOLDGASS')
    ax.errorbar(xcg['LOGMSTAR'][xup], xcg['LOGMH2_LIM'][xup]-xcg['LOGMSTAR'][xup], markersize=5, fmt='v', c='Grey', ecolor='Grey', zorder=10)

    size2 = len(vrt.index[vrt['LOGMH2'].notnull()].tolist())
    ax.errorbar(vrt['LOGMSTAR'], vrt['LOGMH2']-vrt['LOGMSTAR'], markersize=5, fmt='o', c='Grey', ecolor='Grey', alpha=0.5, zorder=20, label=str(size2)+' - VERTICO')
    
    size3 = len(hrc.index[hrc['LOGMH2'].notnull()].tolist())
    ax.errorbar(hrc['LOGMSTAR'], hrc['LOGMH2']-hrc['LOGMSTAR'], markersize=5, fmt='o', c='Grey', ecolor='Grey', alpha=0.5, zorder=25, label=str(size3)+' - HERACLES')

    ax.set_ylabel('log $M_{H2}/M_{*}$')      
    ax.set_xlabel('log $M_{*}$ $[M_{\odot}]$')

    #binned stuff
    values = weighting('MSTAR_BIN','LOGMH2','LOGMSTAR')
    ax.errorbar(values[0], values[1], values[2], markersize=10, fmt='o', c='black', mfc='red', zorder=100)
    ax.errorbar(values[0], values[3], markersize=8, fmt='D', c='black', mfc='dodgerblue', zorder=100)

    #ax.set_xlim(8.5, 11.5)
    #ax.set_ylim(8.25, 11)

    #formatting plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.show()

    #------------------------

    fig, ax = plt.subplots()

    size = len(Mh2.index[jngl['LOGMH2_RYAN'].notnull()].tolist())
    ax.errorbar(sfr-Ms, Mh2-Ms, markersize=5, fmt='o', c='Grey', ecolor='Grey', alpha=0.5, zorder=15, label=str(size)+' - JINGLE')
    
    size1 = len(xcg.index[xcg['LOGMH2'].notnull()].tolist()) + len(xcg.index[xcg['LOGMH2_LIM'].notnull()].tolist())
    ax.errorbar(xcg['LOGSFR'][x]-xcg['LOGMSTAR'][x], xcg['LOGMH2'][x]-xcg['LOGMSTAR'][x], markersize=5, fmt='o', c='Grey', ecolor='Grey', alpha=0.5, zorder=10, label=str(size1)+' - xCOLDGASS')
    ax.errorbar(xcg['LOGSFR'][xup]-xcg['LOGMSTAR'][xup], xcg['LOGMH2_LIM'][xup]-xcg['LOGMSTAR'][xup], markersize=5, fmt='v', c='Grey', ecolor='Grey', zorder=10)

    size2 = len(vrt.index[vrt['LOGMH2'].notnull()].tolist())
    ax.errorbar(vrt['LOGSFR']-vrt['LOGMSTAR'], vrt['LOGMH2']-vrt['LOGMSTAR'], markersize=5, fmt='o', c='Grey', ecolor='Grey', alpha=0.5, zorder=20, label=str(size2)+' - VERTICO')
    
    size3 = len(hrc.index[hrc['LOGMH2'].notnull()].tolist())
    ax.errorbar(hrc['LOGSFR']-hrc['LOGMSTAR'], hrc['LOGMH2']-hrc['LOGMSTAR'], markersize=5, fmt='o', c='Grey', ecolor='Grey', alpha=0.5, zorder=25, label=str(size3)+' - HERACLES')

    ax.set_ylabel('log $M_{H2}/M_{*}$')      
    ax.set_xlabel('log $sSFR$  $[yr^{-1}]$')

    #binned stuff
    values = weighting2('SFR_BIN','LOGMH2','LOGMSTAR','LOGSFR')
    ax.errorbar(values[0], values[1], values[2], markersize=10, fmt='o', c='black', mfc='red', zorder=100)
    ax.errorbar(values[0], values[3], markersize=8, fmt='D', c='black', mfc='dodgerblue', zorder=100)

    #ax.set_xlim(8.5, 11.5)
    #ax.set_ylim(8.25, 11)

    #formatting plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.show()




if __name__ == '__main__':
    MH2()