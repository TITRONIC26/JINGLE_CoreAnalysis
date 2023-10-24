"""
Script meant to plot and analyze the points outlined for Milestone 3 of IndStudy.
"""
#import libraries here
import pandas as pd
import numpy as np
import matplotlib
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
-Plot SFR by Mstar with JINGLE galaxies color coordinated by MH2/M* ratio
-Plot G/D over M* and SFR for 80 detections
-Plot SFR by M* color coded  by G/D
-Plot the jingle and xcg gas/m* plots again
"""

def dist(src):
    points = [(src['LOGMSTAR_MAGPHYS'][m], (src['LOGSFR_MAGPHYS']-src['LOGMSTAR_MAGPHYS'])[m]) for m in range(len(src['LOGSFR_MAGPHYS']))]
    P1 = np.array([9,-9.822])
    P2 = np.array([12,-10.854])

    distance = [-(np.cross(P2 - P1, P1 - P3)/ np.linalg.norm(P2 - P1)) for P3 in points]

    return distance

def f(x):
    return -0.344*(x - 9) - 9.822

def binning(src):
    src['DIST'] = dist(src)

    bins = np.arange(-2,2,0.5)
    print(len(bins))
    labels = [0,1,2,3,4,5,6]
    bin = (pd.cut(src['DIST'], bins, labels=labels))

    src['DIST_BINS'] = bin

    return src

def G2DR(src):
    j = src.index[(src['MGAS_FLAG']==3) | (src['MGAS_FLAG']==6)].tolist()
    bins = [0,1,2,3,4,5,6]

    gds = []

    for x in bins:
        ratio = (src['LOGMGAS'] - src['LOGMDUST_DELOOZE'])[j].loc[src['DIST_BINS']==x]

        GD_ratio = ratio.mean()

        gds.append(GD_ratio)

    m = (2.525 - 2.298) / (6-2)
    b = 2.525 - m*6

    for i in range(len(gds)):
        if mt.isnan(gds[i]):
            gds[i] = m*i + b

    src['LOGMGAS_NEW'] = src['LOGMGAS'].copy()

    for x in bins:
        src.loc[(src['DIST_BINS'] == x) & (src['MGAS_FLAG']==4), 'LOGMGAS_NEW'] = gds[x] + src['LOGMDUST_DELOOZE']
    
    #FF.print_full(src)

    return src

def plot1(src):
    ratio = src['LOGMH2'] - src['LOGMSTAR_MAGPHYS']

    color_map = matplotlib.cm.get_cmap('rainbow')
    plt.set_cmap(color_map)

    plt.scatter(src['LOGMSTAR_MAGPHYS'], src['LOGSFR_MAGPHYS'], c=ratio, edgecolors= "black")
    
    plt.colorbar(label='$M_{H2}/M_{*}$')

    plt.xlabel('log $M_{*}$')      
    plt.ylabel('log SFR $[M_{\odot} yr^{-1}]$')

    plt.show()

def plot2(src):
    j = src.index[(src['MGAS_FLAG']==3) | (src['MGAS_FLAG']==6)].tolist()
    jd = src.index[(src['MGAS_FLAG']==4)].tolist()

    ratio = src['LOGMGAS'] - src['LOGMDUST_DELOOZE']

    fig, ax = plt.subplots()

    ax.errorbar(src['LOGMSTAR_MAGPHYS'][j], ratio[j], markersize=5, fmt='o', c='blue', ecolor='Grey', alpha=0.5, label=str(len(j))+' - JINGLE')

    ax.set_xlabel('log $M_{*}$')
    ax.set_ylabel('$G/D$')

    #formatting plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.show()

    fig, ax = plt.subplots()

    ax.errorbar(src['LOGSFR_MAGPHYS'][j], ratio[j], markersize=5, fmt='o', c='blue', ecolor='Grey', alpha=0.5, label=str(len(j))+' - JINGLE')

    ax.set_xlabel('log SFR $[M_{\odot} yr^{-1}]$')
    ax.set_ylabel('$G/D$')

    #formatting plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.show()

    color_map = matplotlib.cm.get_cmap('rainbow')
    plt.set_cmap(color_map)

    plt.scatter(src['LOGMSTAR_MAGPHYS'][j], (src['LOGSFR_MAGPHYS']-src['LOGMSTAR_MAGPHYS'])[j], c=ratio[j], edgecolors= "black")

    plt.colorbar(label='G/D')

    plt.xlabel('log $M_{*}$')      
    plt.ylabel('log sSFR $[yr^{-1}]$')

    plt.xlim(8.75, 11.5)

    plt.show()

def plot3(src):
    j = src.index[(src['MGAS_FLAG']==3) | (src['MGAS_FLAG']==6)].tolist()
    jd = src.index[(src['MGAS_FLAG']==4)].tolist()

    ratio = src['LOGMGAS'] - src['LOGMDUST_DELOOZE']

    color_map = matplotlib.cm.get_cmap('rainbow')
    plt.set_cmap(color_map)

    plt.scatter(src['LOGMSTAR_MAGPHYS'][j], src['Z_PP04_O3N2'][j], c=ratio[j], edgecolors= "black")
    
    plt.colorbar(label='G/D')

    plt.xlabel('log $M_{*}$')      
    plt.ylabel('Z_PP04_O3N2')

    plt.show()

def plot4(src):
    j = src.index[(src['MGAS_FLAG']==3) | (src['MGAS_FLAG']==6)].tolist()
    jd = src.index[(src['MGAS_FLAG']==4)].tolist()

    ratio = src['LOGMGAS'] - src['LOGMDUST_DELOOZE']

    color_map = matplotlib.cm.get_cmap('rainbow')
    plt.set_cmap(color_map)

    plt.scatter(src['LOGMSTAR_MAGPHYS'][j], (src['LOGSFR_MAGPHYS']-src['LOGMSTAR_MAGPHYS'])[j], c=ratio[j], edgecolors= "black")
    
    xs = np.linspace(8.5,12, 100)

    plt.plot(xs, f(xs), linestyle='--', color='Grey')

    plt.colorbar(label='G/D')

    plt.xlabel('log $M_{*}$')      
    plt.ylabel('log sSFR $[yr^{-1}]$')

    plt.xlim(8.75, 11.5)

    plt.show()

def plot5(src):
    j = src.index[(src['MGAS_FLAG']==3) | (src['MGAS_FLAG']==6)].tolist()
    jd = src.index[(src['MGAS_FLAG']==4)].tolist()

    ratio = src['DIST_BINS']

    color_map = matplotlib.cm.get_cmap('rainbow')
    plt.set_cmap(color_map)

    plt.scatter(src['LOGMSTAR_MAGPHYS'], (src['LOGSFR_MAGPHYS']-src['LOGMSTAR_MAGPHYS']), c=ratio, edgecolors= "black")
    
    xs = np.linspace(8.5,12, 100)

    plt.plot(xs, f(xs), linestyle='--', color='Grey')

    plt.colorbar(label='Distance')

    plt.xlabel('log $M_{*}$')      
    plt.ylabel('log sSFR $[yr^{-1}]$')

    plt.xlim(8.75, 11.5)

    plt.show()

def plot6(src):
    j = src.index[(src['MGAS_FLAG']==3) | (src['MGAS_FLAG']==6)].tolist()
    jd = src.index[(src['MGAS_FLAG']==4)].tolist()

    src = G2DR(src)

    Ms = src['LOGMSTAR_MAGPHYS']
    sfr = src['LOGSFR_MAGPHYS']
    Mg = src['LOGMGAS_NEW'] - Ms

    fig, ax = plt.subplots()

    s1 = len(Mg.index[Mg.notnull()].tolist())
    ax.errorbar(Ms[j], Mg[j], markersize=5, fmt='o', c='Blue', ecolor='Grey', alpha=0.5, zorder=15, label=str(len(j))+' - JINGLE')
    ax.errorbar(Ms[jd], Mg[jd], markersize=5, fmt='v', c='Grey', ecolor='Grey', zorder=15, label=str(len(jd))+' - JINGLE')

    ax.set_ylabel('log $M_{Gas}/M_{*}$')      
    ax.set_xlabel('log $M_{*}$ $[M_{\odot}]$')

    #ax.set_xlim(8.5, 11.5)
    #ax.set_ylim(8.25, 11)

    #formatting plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.show()

    fig, ax = plt.subplots()

    ax.errorbar((sfr-Ms)[j], Mg[j], markersize=5, fmt='o', c='Blue', ecolor='Grey', alpha=0.5, zorder=15, label=str(len(j))+' - JINGLE')
    ax.errorbar((sfr-Ms)[jd], Mg[jd], markersize=5, fmt='v', c='Grey', ecolor='Grey', zorder=15, label=str(len(jd))+' - JINGLE')

    ax.set_ylabel('log $M_{Gas}/M_{*}$')      
    ax.set_xlabel('log $M_{*}$ $[M_{\odot}]$')

    #ax.set_xlim(8.5, 11.5)
    #ax.set_ylim(8.25, 11)

    #formatting plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    plt.show()


def main():
    df = CF.src1.copy()
    df = binning(df)
    #FF.print_full(df)

    #plot1(df)
    #plot2(df)
    #plot3(df)
    #plot4(df)
    #plot5(df)
    plot6(df)

    
if __name__ == '__main__':
    main()