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

cc = ['Blue','Grey','Red','Green','Purple','Orange','Yellow']

def all():
    jngl = CF.src1.copy()
    xcg = CF.src6.copy()
    vrt = CF.src5.copy()
    hrc = CF.src7.copy()

    Ms_list = [jngl['LOGMSTAR_MAGPHYS'],xcg['LOGMSTAR'],vrt['LOGMSTAR'],hrc['LOGMSTAR']]
    sfr_list = [jngl['LOGSFR_MAGPHYS'],xcg['LOGSFR'],vrt['LOGSFR'],hrc['LOGSFR']]
    Mh1_list = [jngl['LOGMH1'],xcg['LOGMH1'],vrt['LOGMH1'],hrc['LOGMH1']]

    xcg.loc[xcg['LOGMH2'].isnull(), 'LOGMH2'] = 0
    xcg.loc[xcg['LOGMH2_LIM'].isnull(), 'LOGMH2_LIM'] = 0
    xcg['LOGMH2_ALL'] = xcg['LOGMH2'] + xcg['LOGMH2_LIM']

    Mh2_list = [jngl['LOGMH2'],xcg['LOGMH2_ALL'],vrt['LOGMH2'],hrc['LOGMH2']]

    Ms = pd.concat(Ms_list, ignore_index=True)
    sfr = pd.concat(sfr_list, ignore_index=True)
    Mh1 = pd.concat(Mh1_list, ignore_index=True)
    Mh2 = pd.concat(Mh2_list, ignore_index=True)

    MhR = Mh2 - Mh1

    all = pd.concat([Ms,sfr,Mh1,Mh2,MhR], axis=1)
    all = all.rename(columns={0:'LOGMSTAR',1:'LOGSFR',2:'LOGMH1',3:'LOGMH2',4:'RATIO'})
    FF.print_full(all)

    return all

def binning():
    all_data = all()

    #Mstar binning
    range_all = 11.5-8.5
    steps = 9
    step = range_all/steps

    bins = np.arange(8.5,11.5,step)
    labels = [0,1,2,3,4,5,6,7]

    bin = (pd.cut(all_data['LOGMSTAR'], bins, labels=labels))

    all_data['MSTAR_BIN'] = bin

    #sSFR binning
    range_all = (-8.5)-(-12.5)
    steps = 6
    step = range_all/steps

    bins = np.arange(-12.5,-8.5,step)
    labels = [0,1,2,3,4]

    bin = (pd.cut(all_data['LOGSFR']-all_data['LOGMSTAR'], bins, labels=labels))
    #print(bin)

    all_data['SFR_BIN'] = bin
    """
    #MH1/MSTAR binning
    range_all = (1)-(-3.5)
    steps = 9
    step = range_all/steps

    bins = np.arange(-3.5,1,step)
    labels = [0,1,2,3,4,5,6,7]

    bin = (pd.cut(all_data['LOGMH1']-all_data['LOGMSTAR'], bins, labels=labels))
    #print(bin)

    all_data['MH1_BIN'] = bin
    """
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
        ave.append(np.average(group[y]-group[x], weights=group[x]))
        ave_e.append(np.std(group[y]) / np.sqrt(len(group[y])-group[x]) / np.sqrt(len(group[y])))
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
        ave.append(np.average(group[y]))#-group[x1]))#, weights=group[x2]-group[x1]))
        ave_e.append(np.std(group[y]) / np.sqrt(len(group[y])))#-group[x1]) / np.sqrt(len(group[y])))
        med.append(np.median(group[y]))#-group[x1]))
        xs.append(np.average(group[x2]-group[x1]))

    return (xs, ave, ave_e, med)

def weighting3(col,y1,y2,x1,x2):
    all_data = binning()

    groups = all_data.groupby(col)

    ave = []
    ave_e = []
    med = []
    xs = []

    for key, group in groups:
        group = group.dropna()
        ave.append(np.average(group[y2]-group[y1]))#, weights=group[x2]-group[x1]))
        ave_e.append(np.std(group[y2]-group[y1]) / np.sqrt(len(group[y2])))#-group[x1]) / np.sqrt(len(group[y])))
        med.append(np.median(group[y2]-group[y1]))#-group[x1]))
        xs.append(np.average(group[x2]-group[x1]))

    return (xs, ave, ave_e, med)

def MH1(jngl, show=False):
    #['LOGMSTAR_MAGPHYS','LOGMSTAR_MAGPHYS_ERR','LOGSFR_MAGPHYS','LOGSFR_MAGPHYS_ERR','LOGMDUST_DELOOZE','LOGMH2_RYAN','LOGMH1_MATT','H1_FLAG','LOGMGAS','LOGMGAS_ERR','MGAS_FLAG']
    #FF.print_full(jngl)
    j = jngl.index[jngl['H1_FLAG']==1].tolist()
    jup = jngl.index[(jngl['H1_FLAG']==0)].tolist()

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

    Mh1 = jngl['LOGMH1']

    #-----------------------------

    fig, ax = plt.subplots()

    size = len(Mh1.index[jngl['LOGMH1'].notnull()].tolist())
    ax.errorbar(Ms[j], Mh1[j]-Ms[j], markersize=5, fmt='o', c=cc[0], ecolor='Grey', alpha=0.5, zorder=15, label=str(size)+' - JINGLE')
    ax.errorbar(Ms[jup], Mh1[jup]-Ms[jup], markersize=5, fmt='v', c=cc[0], ecolor='Grey', zorder=15)

    if show == True:
        size1 = len(xcg.index[xcg['LOGMH1'].notnull()].tolist())
        ax.errorbar(xcg['LOGMSTAR'][x], xcg['LOGMH1'][x]-xcg['LOGMSTAR'][x], markersize=5, fmt='o', c=cc[1], ecolor='Grey', alpha=0.5, zorder=10, label=str(size1)+' - xCOLDGASS')
        ax.errorbar(xcg['LOGMSTAR'][xup], xcg['LOGMH1'][xup]-xcg['LOGMSTAR'][xup], markersize=5, fmt='v', c=cc[1], ecolor='Grey', zorder=10)

        size2 = len(vrt.index[vrt['LOGMH1'].notnull()].tolist())
        ax.errorbar(vrt['LOGMSTAR'], vrt['LOGMH1']-vrt['LOGMSTAR'], markersize=5, fmt='o', c=cc[2], ecolor='Grey', alpha=0.5, zorder=20, label=str(size2)+' - VERTICO')
        
        size3 = len(hrc.index[hrc['LOGMH1'].notnull()].tolist())
        ax.errorbar(hrc['LOGMSTAR'], hrc['LOGMH1']-hrc['LOGMSTAR'], markersize=5, fmt='o', c=cc[3], ecolor='Grey', alpha=0.5, zorder=25, label=str(size3)+' - HERACLES')

    ax.set_ylabel('log $M_{HI}/M_{*}$')      
    ax.set_xlabel('log $M_{*}$ $[M_{\odot}]$')

    #ax.set_xlim(8.5, 11.5)
    #ax.set_ylim(8.25, 11)

    #formatting plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.show()

    #------------------------

    fig, ax = plt.subplots()

    size = len(Mh1.index[jngl['LOGMH1'].notnull()].tolist())
    ax.errorbar(sfr[j]-Ms[j], Mh1[j]-Ms[j], markersize=5, fmt='o', c=cc[0], ecolor='Grey', alpha=0.5, zorder=15, label=str(size)+' - JINGLE')
    ax.errorbar(sfr[jup]-Ms[jup], Mh1[jup]-Ms[jup], markersize=5, fmt='v', c=cc[0], ecolor='Grey', zorder=15)

    if show == True:
        size1 = len(xcg.index[xcg['LOGMH1'].notnull()].tolist())
        ax.errorbar(xcg['LOGSFR'][x]-xcg['LOGMSTAR'][x], xcg['LOGMH1'][x]-xcg['LOGMSTAR'][x], markersize=5, fmt='o', c=cc[1], ecolor='Grey', alpha=0.5, zorder=10, label=str(size1)+' - xCOLDGASS')
        ax.errorbar(xcg['LOGSFR'][xup]-xcg['LOGMSTAR'][xup], xcg['LOGMH1'][xup]-xcg['LOGMSTAR'][xup], markersize=5, fmt='v', c=cc[1], ecolor='Grey', zorder=10)

        size2 = len(vrt.index[vrt['LOGMH1'].notnull()].tolist())
        ax.errorbar(vrt['LOGSFR']-vrt['LOGMSTAR'], vrt['LOGMH1']-vrt['LOGMSTAR'], markersize=5, fmt='o', c=cc[2], ecolor='Grey', alpha=0.5, zorder=20, label=str(size2)+' - VERTICO')
        
        size3 = len(hrc.index[hrc['LOGMH1'].notnull()].tolist())
        ax.errorbar(hrc['LOGSFR']-hrc['LOGMSTAR'], hrc['LOGMH1']-hrc['LOGMSTAR'], markersize=5, fmt='o', c=cc[3], ecolor='Grey', alpha=0.5, zorder=25, label=str(size3)+' - HERACLES')

    ax.set_ylabel('log $M_{HI}/M_{*}$')      
    ax.set_xlabel('log $sSFR$  $[yr^{-1}]$')

    #ax.set_xlim(8.5, 11.5)
    #ax.set_ylim(8.25, 11)

    #formatting plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.show()

def MH2(jngl, show=False):
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

    Mh2 = jngl['LOGMH2']

    #-----------------------------

    fig, ax = plt.subplots()

    size = len(Mh2.index[jngl['LOGMH2'].notnull()].tolist())
    ax.errorbar(Ms, Mh2-Ms, markersize=5, fmt='o', c=cc[0], ecolor='Grey', alpha=0.5, zorder=15, label=str(size)+' - JINGLE')
    
    if show == True:
        size1 = len(xcg.index[xcg['LOGMH2'].notnull()].tolist()) + len(xcg.index[xcg['LOGMH2_LIM'].notnull()].tolist())
        ax.errorbar(xcg['LOGMSTAR'][x], xcg['LOGMH2'][x]-xcg['LOGMSTAR'][x], markersize=5, fmt='o', c=cc[1], ecolor='Grey', alpha=0.5, zorder=10, label=str(size1)+' - xCOLDGASS')
        ax.errorbar(xcg['LOGMSTAR'][xup], xcg['LOGMH2_LIM'][xup]-xcg['LOGMSTAR'][xup], markersize=5, fmt='v', c=cc[1], ecolor='Grey', zorder=10)

        size2 = len(vrt.index[vrt['LOGMH2'].notnull()].tolist())
        ax.errorbar(vrt['LOGMSTAR'], vrt['LOGMH2']-vrt['LOGMSTAR'], markersize=5, fmt='o', c=cc[2], ecolor='Grey', alpha=0.5, zorder=20, label=str(size2)+' - VERTICO')
        
        size3 = len(hrc.index[hrc['LOGMH2'].notnull()].tolist())
        ax.errorbar(hrc['LOGMSTAR'], hrc['LOGMH2']-hrc['LOGMSTAR'], markersize=5, fmt='o', c=cc[3], ecolor='Grey', alpha=0.5, zorder=25, label=str(size3)+' - HERACLES')

    ax.set_ylabel('log $M_{H2}/M_{*}$')      
    ax.set_xlabel('log $M_{*}$ $[M_{\odot}]$')

    #ax.set_xlim(8.5, 11.5)
    #ax.set_ylim(8.25, 11)

    #formatting plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.show()

    #------------------------

    fig, ax = plt.subplots()

    size = len(Mh2.index[jngl['LOGMH2'].notnull()].tolist())
    ax.errorbar(sfr-Ms, Mh2-Ms, markersize=5, fmt='o', c=cc[0], ecolor='Grey', alpha=0.5, zorder=15, label=str(size)+' - JINGLE')

    if show == True:
        size1 = len(xcg.index[xcg['LOGMH2'].notnull()].tolist()) + len(xcg.index[xcg['LOGMH2_LIM'].notnull()].tolist())
        ax.errorbar(xcg['LOGSFR'][x]-xcg['LOGMSTAR'][x], xcg['LOGMH2'][x]-xcg['LOGMSTAR'][x], markersize=5, fmt='o', c=cc[1], ecolor='Grey', alpha=0.5, zorder=10, label=str(size1)+' - xCOLDGASS')
        ax.errorbar(xcg['LOGSFR'][xup]-xcg['LOGMSTAR'][xup], xcg['LOGMH2_LIM'][xup]-xcg['LOGMSTAR'][xup], markersize=5, fmt='v', c=cc[1], ecolor='Grey', zorder=10)

        size2 = len(vrt.index[vrt['LOGMH2'].notnull()].tolist())
        ax.errorbar(vrt['LOGSFR']-vrt['LOGMSTAR'], vrt['LOGMH2']-vrt['LOGMSTAR'], markersize=5, fmt='o', c=cc[2], ecolor='Grey', alpha=0.5, zorder=20, label=str(size2)+' - VERTICO')
        
        size3 = len(hrc.index[hrc['LOGMH2'].notnull()].tolist())
        ax.errorbar(hrc['LOGSFR']-hrc['LOGMSTAR'], hrc['LOGMH2']-hrc['LOGMSTAR'], markersize=5, fmt='o', c=cc[3], ecolor='Grey', alpha=0.5, zorder=25, label=str(size3)+' - HERACLES')

    ax.set_ylabel('log $M_{H2}/M_{*}$')      
    ax.set_xlabel('log $sSFR$  $[yr^{-1}]$')

    #ax.set_xlim(8.5, 11.5)
    #ax.set_ylim(8.25, 11)

    #formatting plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.show()

def MH1MH2(jngl, show=False):
    #['LOGMSTAR_MAGPHYS','LOGMSTAR_MAGPHYS_ERR','LOGSFR_MAGPHYS','LOGSFR_MAGPHYS_ERR','LOGMDUST_DELOOZE','LOGMH2_RYAN','LOGMH1_MATT',,'H1_FLAG','LOGMGAS','LOGMGAS_ERR','MGAS_FLAG']
    #FF.print_full(jngl)
    j = jngl.index[jngl['H1_FLAG']==1].tolist()
    jup = jngl.index[jngl['H1_FLAG']==0].tolist()

    xcg = CF.src6.copy()
    #['LOGMSTAR','LOGSFR','LOGSFR_ERR','LOGMH1','H1_FLAG','LOGMH2','LOGMH2_ERR','LOGMH2_LIM','GROUPID','ENV_CODE','NGAL','LOGMH']
    #FF.print_full(xcg)
    x = xcg.index[(xcg['H1_FLAG']!=99) & (xcg['LOGMH2'].notnull())].tolist()
    xup = xcg.index[(xcg['H1_FLAG']==99) & (xcg['LOGMH2_LIM'].notnull())].tolist()
    
    vrt = CF.src5.copy()
    #['ID','RA','DEC','z','LOGMH1','LOGMSTAR','LOGMSTAR_ERR','LOGSFR','LOGSFR_ERR','LOGMH2','LOGMH2_ERR','LOGMH1_ERR]
    #FF.print_full(vrt)

    hrc = CF.src7.copy()
    #['ID','RA','DEC','Vel','z','LOGMSTAR','LOGMSTAR_ERR','LOGSFR','LOGSFR_ERR','LOGMH2','LOGMH2_ERR','LOGMH1','LOGMH1_ERR']
    #FF.print_full(hrc)

    Ms = jngl['LOGMSTAR_MAGPHYS']
    sfr = jngl['LOGSFR_MAGPHYS']
    Mh1 = jngl['LOGMH1']
    Mh2 = jngl['LOGMH2']
    MhR = Mh2 - Mh1

    #print(MhR)

    #-----------------------------

    fig, ax = plt.subplots()

    s1 = len(MhR.index[MhR.notnull()].tolist())
    ax.errorbar(Ms[j], MhR[j], markersize=5, fmt='o', c=cc[0], ecolor='Grey', alpha=0.5, zorder=15, label=str(s1)+' - JINGLE')
    ax.errorbar(Ms[jup], MhR[jup], markersize=5, fmt='v', c=cc[0], ecolor='Grey', zorder=15)

    if show == True:
        MhR = xcg['LOGMH2'] - xcg['LOGMH1']
        MhRup = xcg['LOGMH2_LIM'] - xcg['LOGMH1']

        size1 = len(xcg.index[MhR.notnull()].tolist()) + len(xcg.index[MhRup.notnull()].tolist())
        ax.errorbar(xcg['LOGMSTAR'][x], MhR[x], markersize=5, fmt='o', c=cc[1], ecolor='Grey', alpha=0.5, zorder=10, label=str(size1)+' - xCOLDGASS')
        ax.errorbar(xcg['LOGMSTAR'][xup], MhRup[xup], markersize=5, fmt='v', c=cc[1], ecolor='Grey', zorder=10)

        MhR = vrt['LOGMH2'] - vrt['LOGMH1']

        size2 = len(vrt.index[MhR.notnull()].tolist())
        ax.errorbar(vrt['LOGMSTAR'], MhR, markersize=5, fmt='o', c=cc[2], ecolor='Grey', alpha=0.5, zorder=20, label=str(size2)+' - VERTICO')
        
        MhR = hrc['LOGMH2'] - hrc['LOGMH1']

        size3 = len(hrc.index[MhR.notnull()].tolist())
        ax.errorbar(hrc['LOGMSTAR'], MhR, markersize=5, fmt='o', c=cc[3], ecolor='Grey', alpha=0.5, zorder=25, label=str(size3)+' - HERACLES')

    ax.set_ylabel('log $M_{H2}/M_{HI}$')      
    ax.set_xlabel('log $M_{*}$ $[M_{\odot}]$')

    #ax.set_xlim(8.5, 11.5)
    #ax.set_ylim(8.25, 11)

    #formatting plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.show()

    #-----------------------------

    Ms = jngl['LOGMSTAR_MAGPHYS']
    sfr = jngl['LOGSFR_MAGPHYS']
    Mh1 = jngl['LOGMH1']
    Mh2 = jngl['LOGMH2']
    MhR = Mh2 - Mh1

    fig, ax = plt.subplots()

    s1 = len(MhR.index[MhR.notnull()].tolist())
    ax.errorbar((sfr-Ms)[j], MhR[j], markersize=5, fmt='o', c=cc[0], ecolor='Grey', alpha=0.5, zorder=15, label=str(s1)+' - JINGLE')
    ax.errorbar((sfr-Ms)[jup], MhR[jup], markersize=5, fmt='v', c=cc[0], ecolor='Grey', zorder=15)

    if show == True:
        MhR = xcg['LOGMH2'] - xcg['LOGMH1']
        MhRup = xcg['LOGMH2_LIM'] - xcg['LOGMH1']

        size1 = len(xcg.index[MhR.notnull()].tolist()) + len(xcg.index[MhRup.notnull()].tolist())
        ax.errorbar((xcg['LOGSFR'] - xcg['LOGMSTAR'])[x], MhR[x], markersize=5, fmt='o', c=cc[1], ecolor='Grey', alpha=0.5, zorder=10, label=str(size1)+' - xCOLDGASS')
        ax.errorbar((xcg['LOGSFR'] - xcg['LOGMSTAR'])[xup], MhRup[xup], markersize=5, fmt='v', c=cc[1], ecolor='Grey', zorder=10)

        MhR = vrt['LOGMH2'] - vrt['LOGMH1']

        size2 = len(vrt.index[MhR.notnull()].tolist())
        ax.errorbar((vrt['LOGSFR'] - vrt['LOGMSTAR']), MhR, markersize=5, fmt='o', c=cc[2], ecolor='Grey', alpha=0.5, zorder=20, label=str(size2)+' - VERTICO')
        
        MhR = hrc['LOGMH2'] - hrc['LOGMH1']

        size3 = len(hrc.index[MhR.notnull()].tolist())
        ax.errorbar((hrc['LOGSFR'] - hrc['LOGMSTAR']), MhR, markersize=5, fmt='o', c=cc[3], ecolor='Grey', alpha=0.5, zorder=25, label=str(size3)+' - HERACLES')

    ax.set_ylabel('log $M_{H2}/M_{HI}$')      
    ax.set_xlabel('log $sSFR$  $[yr^{-1}]$')

    #ax.set_xlim(8.5, 11.5)
    #ax.set_ylim(8.25, 11)

    #formatting plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.show()

def MH1byMH2(jngl, show=False):
    #['LOGMSTAR_MAGPHYS','LOGMSTAR_MAGPHYS_ERR','LOGSFR_MAGPHYS','LOGSFR_MAGPHYS_ERR','LOGMDUST_DELOOZE','LOGMH2_RYAN','LOGMH1_MATT',,'H1_FLAG','LOGMGAS','LOGMGAS_ERR','MGAS_FLAG']
    #FF.print_full(jngl)
    j = jngl.index[jngl['H1_FLAG']==1].tolist()
    jup = jngl.index[jngl['H1_FLAG']==0].tolist()

    xcg = CF.src6.copy()
    #['LOGMSTAR','LOGSFR','LOGSFR_ERR','LOGMH1','H1_FLAG','LOGMH2','LOGMH2_ERR','LOGMH2_LIM','GROUPID','ENV_CODE','NGAL','LOGMH']
    #FF.print_full(xcg)
    x = xcg.index[(xcg['H1_FLAG']!=99) & (xcg['LOGMH2'].notnull())].tolist()
    xup = xcg.index[(xcg['H1_FLAG']==99) & (xcg['LOGMH2_LIM'].notnull())].tolist()
    
    vrt = CF.src5.copy()
    #['ID','RA','DEC','z','LOGMH1','LOGMSTAR','LOGMSTAR_ERR','LOGSFR','LOGSFR_ERR','LOGMH2','LOGMH2_ERR','LOGMH1_ERR]
    #FF.print_full(vrt)

    hrc = CF.src7.copy()
    #['ID','RA','DEC','Vel','z','LOGMSTAR','LOGMSTAR_ERR','LOGSFR','LOGSFR_ERR','LOGMH2','LOGMH2_ERR','LOGMH1','LOGMH1_ERR']
    #FF.print_full(hrc)

    Ms = jngl['LOGMSTAR_MAGPHYS']
    sfr = jngl['LOGSFR_MAGPHYS']
    Mh1 = jngl['LOGMH1']
    Mh2 = jngl['LOGMH2']

    #-----------------------------

    fig, ax = plt.subplots()

    ax.axline((0,0), slope=1, c='Black', linestyle='--')

    s1 = len(jngl.index[(Mh1.notnull()) & (Mh2.notnull())].tolist())
    ax.errorbar((Mh1-Ms)[j], (Mh2-Ms)[j], markersize=5, fmt='o', c=cc[0], ecolor='Grey', alpha=0.5, zorder=15, label=str(s1)+' - JINGLE')
    ax.errorbar((Mh1-Ms)[jup], (Mh2-Ms)[jup], markersize=5, fmt='v', c=cc[0], ecolor='Grey', zorder=15)

    if show == True:
        Mh1 = xcg['LOGMH1'] - xcg['LOGMSTAR']
        Mh2 = xcg['LOGMH2'] - xcg['LOGMSTAR']
        Mh2up = xcg['LOGMH2_LIM'] - xcg['LOGMSTAR']

        size1 = len(xcg.index[(Mh1.notnull()) & (Mh2.notnull())].tolist()) + len(xcg.index[(Mh1.notnull()) & (Mh2up.notnull())].tolist())
        ax.errorbar(Mh1[x], Mh2[x], markersize=5, fmt='o', c=cc[1], ecolor='Grey', alpha=0.5, zorder=10, label=str(size1)+' - xCOLDGASS')
        ax.errorbar(Mh1[xup], Mh2up[xup], markersize=5, fmt='v', c=cc[1], ecolor='Grey', zorder=10)

        Mh1 = vrt['LOGMH1'] - vrt['LOGMSTAR']
        Mh2 = vrt['LOGMH2'] - vrt['LOGMSTAR']

        size2 = len(vrt.index[(Mh1.notnull()) & (Mh2.notnull())].tolist())
        ax.errorbar(Mh1, Mh2, markersize=5, fmt='o', c=cc[2], ecolor='Grey', alpha=0.5, zorder=20, label=str(size2)+' - VERTICO')
        
        Mh1 = hrc['LOGMH1'] - hrc['LOGMSTAR']
        Mh2 = hrc['LOGMH2'] - hrc['LOGMSTAR']

        size3 = len(hrc.index[(Mh1.notnull()) & (Mh2.notnull())].tolist())
        ax.errorbar(Mh1, Mh2, markersize=5, fmt='o', c=cc[3], ecolor='Grey', alpha=0.5, zorder=25, label=str(size3)+' - HERACLES')

    ax.set_ylabel('log $M_{H2}/M_{*}$')      
    ax.set_xlabel('log $M_{HI}/M_{*}$')

    #ax.set_xlim(8.5, 11.5)
    #ax.set_ylim(8.25, 11)

    #formatting plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.show()

def MGAS(jngl, show=False):
    #['LOGMSTAR_MAGPHYS','LOGMSTAR_MAGPHYS_ERR','LOGSFR_MAGPHYS','LOGSFR_MAGPHYS_ERR','LOGMDUST_DELOOZE','LOGMH2_RYAN','LOGMH1_MATT',,'H1_FLAG','LOGMGAS','LOGMGAS_ERR','MGAS_FLAG']
    #FF.print_full(jngl)
    j = jngl.index[(jngl['MGAS_FLAG']==3) | (jngl['MGAS_FLAG']==6)].tolist()
    jh1 = jngl.index[(jngl['MH_FLAG']==1)].tolist()
    jh2 = jngl.index[(jngl['MH_FLAG']==2)].tolist()
    jd = jngl.index[(jngl['MGAS_FLAG']==4)].tolist()
    jnon = jngl.index[(jngl['MGAS_FLAG']==5) | (jngl['MH_FLAG']==0)].tolist()

    xcg = CF.src6.copy()
    #['LOGMSTAR','LOGSFR','LOGSFR_ERR','LOGMH1','H1_FLAG','LOGMH2','LOGMH2_ERR','LOGMH2_LIM','GROUPID','ENV_CODE','NGAL','LOGMH']
    #FF.print_full(xcg)
    x = xcg.index[(xcg['H1_FLAG']!=99) & (xcg['LOGMH2'].notnull())].tolist()
    xup = xcg.index[(xcg['H1_FLAG']==99) & (xcg['LOGMH2_LIM'].notnull())].tolist()
    
    vrt = CF.src5.copy()
    #['ID','RA','DEC','z','LOGMH1','LOGMSTAR','LOGMSTAR_ERR','LOGSFR','LOGSFR_ERR','LOGMH2','LOGMH2_ERR','LOGMH1_ERR]
    #FF.print_full(vrt)

    hrc = CF.src7.copy()
    #['ID','RA','DEC','Vel','z','LOGMSTAR','LOGMSTAR_ERR','LOGSFR','LOGSFR_ERR','LOGMH2','LOGMH2_ERR','LOGMH1','LOGMH1_ERR']
    #FF.print_full(hrc)

    Ms = jngl['LOGMSTAR_MAGPHYS']
    sfr = jngl['LOGSFR_MAGPHYS']
    Mg = jngl['LOGMGAS'] - Ms

    #print(MhR)

    #-----------------------------

    fig, ax = plt.subplots()

    s1 = len(Mg.index[Mg.notnull()].tolist())
    ax.errorbar(Ms[j], Mg[j], markersize=5, fmt='o', c=cc[0], ecolor='Grey', alpha=0.5, zorder=15, label=str(len(j)+len(jd))+' - JINGLE')
    ax.errorbar(Ms[jd], Mg[jd], markersize=5, fmt='v', c=cc[0], ecolor='Grey', zorder=15)#, label=str(len(jd))+' - JINGLE')
    #ax.errorbar(Ms[jh1], Mg[jh1], markersize=5, fmt='s', c=cc[2], ecolor='Grey', alpha=0.5, zorder=15, label=str(len(jh1))+' - JINGLE')
    #ax.errorbar(Ms[jh2], Mg[jh2], markersize=5, fmt='d', c=cc[3], ecolor='Grey', alpha=0.5, zorder=15, label=str(len(jh2))+' - JINGLE')
    #ax.errorbar(Ms[jnon], Mg[jnon], markersize=5, fmt='h', c=cc[4], ecolor='Grey', alpha=0.5, zorder=15, label=str(len(jnon))+' - JINGLE')

    if show == True:
        #"""
        MhR = xcg['LOGMH2'] - xcg['LOGMH1']
        MhRup = xcg['LOGMH2_LIM'] - xcg['LOGMH1']

        size1 = len(xcg.index[MhR.notnull()].tolist()) + len(xcg.index[MhRup.notnull()].tolist())
        ax.errorbar(xcg['LOGMSTAR'][x], MhR[x], markersize=5, fmt='o', c=cc[1], ecolor='Grey', alpha=0.5, zorder=10, label=str(size1)+' - xCOLDGASS')
        ax.errorbar(xcg['LOGMSTAR'][xup], MhRup[xup], markersize=5, fmt='v', c=cc[1], ecolor='Grey', zorder=10)

        MhR = vrt['LOGMH2'] - vrt['LOGMH1']

        size2 = len(vrt.index[MhR.notnull()].tolist())
        ax.errorbar(vrt['LOGMSTAR'], MhR, markersize=5, fmt='o', c=cc[2], ecolor='Grey', alpha=0.5, zorder=20, label=str(size2)+' - VERTICO')
        
        MhR = hrc['LOGMH2'] - hrc['LOGMH1']

        size3 = len(hrc.index[MhR.notnull()].tolist())
        ax.errorbar(hrc['LOGMSTAR'], MhR, markersize=5, fmt='o', c=cc[3], ecolor='Grey', alpha=0.5, zorder=25, label=str(size3)+' - HERACLES')

        #"""

    ax.set_ylabel('log $M_{Gas}/M_{*}$')      
    ax.set_xlabel('log $M_{*}$ $[M_{\odot}]$')

    #ax.set_xlim(8.5, 11.5)
    #ax.set_ylim(8.25, 11)

    #formatting plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.show()

    #-----------------------------

    Ms = jngl['LOGMSTAR_MAGPHYS']
    sfr = jngl['LOGSFR_MAGPHYS']
    Mg = jngl['LOGMGAS'] - Ms

    fig, ax = plt.subplots()

    s1 = len(Mg.index[Mg.notnull()].tolist())
    ax.errorbar((sfr-Ms)[j], Mg[j], markersize=5, fmt='o', c=cc[0], ecolor='Grey', alpha=0.5, zorder=15, label=str(len(j)+len(jd))+' - JINGLE')
    ax.errorbar((sfr-Ms)[jd], Mg[jd], markersize=5, fmt='v', c=cc[0], ecolor='Grey', zorder=15)#, label=str(len(jd))+' - JINGLE')
    #ax.errorbar((sfr-Ms)[jh1], Mg[jh1], markersize=5, fmt='s', c=cc[2], ecolor='Grey', alpha=0.5, zorder=15, label=str(len(jh1))+' - JINGLE')
    #ax.errorbar((sfr-Ms)[jh2], Mg[jh2], markersize=5, fmt='d', c=cc[3], ecolor='Grey', alpha=0.5, zorder=15, label=str(len(jh2))+' - JINGLE')
    #ax.errorbar((sfr-Ms)[jnon], Mg[jnon], markersize=5, fmt='h', c=cc[4], ecolor='Grey', alpha=0.5, zorder=15, label=str(len(jnon))+' - JINGLE')

    if show == True:
        #"""

        MhR = xcg['LOGMH2'] - xcg['LOGMH1']
        MhRup = xcg['LOGMH2_LIM'] - xcg['LOGMH1']

        size1 = len(xcg.index[MhR.notnull()].tolist()) + len(xcg.index[MhRup.notnull()].tolist())
        ax.errorbar((xcg['LOGSFR'] - xcg['LOGMSTAR'])[x], MhR[x], markersize=5, fmt='o', c=cc[1], ecolor='Grey', alpha=0.5, zorder=10, label=str(size1)+' - xCOLDGASS')
        ax.errorbar((xcg['LOGSFR'] - xcg['LOGMSTAR'])[xup], MhRup[xup], markersize=5, fmt='v', c=cc[1], ecolor='Grey', zorder=10)

        MhR = vrt['LOGMH2'] - vrt['LOGMH1']

        size2 = len(vrt.index[MhR.notnull()].tolist())
        ax.errorbar((vrt['LOGSFR'] - vrt['LOGMSTAR']), MhR, markersize=5, fmt='o', c=cc[2], ecolor='Grey', alpha=0.5, zorder=20, label=str(size2)+' - VERTICO')
        
        MhR = hrc['LOGMH2'] - hrc['LOGMH1']

        size3 = len(hrc.index[MhR.notnull()].tolist())
        ax.errorbar((hrc['LOGSFR'] - hrc['LOGMSTAR']), MhR, markersize=5, fmt='o', c=cc[3], ecolor='Grey', alpha=0.5, zorder=25, label=str(size3)+' - HERACLES')

        #"""

    ax.set_ylabel('log $M_{Gas}/M_{*}$')      
    ax.set_xlabel('log $sSFR$  $[yr^{-1}]$')

    #ax.set_xlim(8.5, 11.5)
    #ax.set_ylim(8.25, 11)

    #formatting plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.show()

def main(src):
    MH1(src)
    MH2(src)
    MH1MH2(src)
    MH1byMH2(src)
    MGAS(src)

if __name__ == '__main__':
    main()