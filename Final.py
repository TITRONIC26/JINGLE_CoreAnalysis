#import libraries here
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
import numpy as np

from scipy.stats import kstest
from scipy.stats import ks_2samp

#plotting scripts
import base_plotter as BPLT
import ratio_plots as RPLT
import plotter as PLT

#essential scripts
import gather_data as GD
import constants as C
import Milestone3 as M3
import callable_functions as CF
import main as MAIN

#data analysis and manipulation
import core_analysis as CA
import data_manipulation as DM
import linmixer as LM

#formatting
import formatting_functions as FF

def f(x):
    return -0.344*(x - 9) - 9.822

def dist(src):
    points = [(src['LOGMSTAR_MAGPHYS'][m], (src['LOGSFR_MAGPHYS']-src['LOGMSTAR_MAGPHYS'])[m]) for m in range(len(src['LOGSFR_MAGPHYS']))]
    P1 = np.array([9,-9.822])
    P2 = np.array([12,-10.854])

    distance = [-(np.cross(P2 - P1, P1 - P3)/ np.linalg.norm(P2 - P1)) for P3 in points]
    
    src['CC'] = distance

    return src

def median(x, y):
    nots = y.index[y.notnull()].tolist()
    
    x = x[nots]
    y = y[nots]

    total_bins = int(len(y)/15)
    print(total_bins)

    bins = np.linspace(x.min(), x.max(), total_bins)
    delta = bins[1]-bins[0]
    idx = np.digitize(x, bins)
    running_median = [np.median(y[idx==k]) for k in range(total_bins)]

    plt.plot(bins-delta/2, running_median, color='White', linewidth=5)
    plt.plot(bins-delta/2, running_median, color='Black', linewidth=3)

def seven_odds():
    src = MAIN.jingle.copy()

    seven = src.index[(src['LOGMGAS']-src['LOGMSTAR_MAGPHYS'] > 0.5)
                      & (src['LOGSFR_MAGPHYS']-src['LOGMSTAR_MAGPHYS'] >= -10.5)
                      & (src['LOGSFR_MAGPHYS']-src['LOGMSTAR_MAGPHYS'] <= -9.4)].tolist()
    
    print(src['JINGLEID'][seven])

    return seven
    
def sSFR_MS(src, xcg = True, seven = True):
    src = dist(src)

    color_map = matplotlib.cm.get_cmap('rainbow_r')
    plt.set_cmap(color_map)

    if xcg == True:
        xcoldgass = CF.src6.copy()

        xn = xcoldgass.index[(xcoldgass['H1_FLAG']!=99) & (xcoldgass['LOGMH2'].notnull())].tolist()
        xup = xcoldgass.index[(xcoldgass['H1_FLAG']==99) & (xcoldgass['LOGMH2_LIM'].notnull())].tolist()

        plt.scatter(xcoldgass['LOGMSTAR'][xn], (xcoldgass['LOGSFR']-xcoldgass['LOGMSTAR'])[xn], c='Black', alpha=0.25, marker='o')
        plt.scatter(xcoldgass['LOGMSTAR'][xup], (xcoldgass['LOGSFR']-xcoldgass['LOGMSTAR'])[xup], c='Black', alpha=0.25, marker=r'$\downarrow$')
        

    plt.scatter(src['LOGMSTAR_MAGPHYS'], (src['LOGSFR_MAGPHYS']-src['LOGMSTAR_MAGPHYS']), c=src['CC'], edgecolors= "black")

    xs = np.linspace(8.5,12, 100)

    plt.plot(xs, f(xs), linestyle='--', color='Grey')

    plt.colorbar(label='MS dist')

    plt.xlabel('log $M_{*}$ [$M_{\odot}$]')      
    plt.ylabel('log sSFR $[yr^{-1}]$')

    plt.xlim(8.75, 11.5)

    if seven == True:
        plt.scatter(src['LOGMSTAR_MAGPHYS'][seven_odds()], (src['LOGSFR_MAGPHYS']-src['LOGMSTAR_MAGPHYS'])[seven_odds()], c='Black', marker='x')

    plt.show()

def Catinella_1(src, xcg = True, seven = True):
    src = dist(src)
    
    color_map = matplotlib.cm.get_cmap('rainbow_r')
    plt.set_cmap(color_map)

    if xcg == True:
        xcoldgass = CF.src6.copy()

        xn = xcoldgass.index[(xcoldgass['H1_FLAG']!=99) & (xcoldgass['LOGMH2'].notnull())].tolist()
        xup = xcoldgass.index[(xcoldgass['H1_FLAG']==99) & (xcoldgass['LOGMH2_LIM'].notnull())].tolist()

        plt.scatter(xcoldgass['LOGMSTAR'][xn], (xcoldgass['LOGMH2']-xcoldgass['LOGMH1'])[xn], c='Black', alpha=0.25, marker='o')
        plt.scatter(xcoldgass['LOGMSTAR'][xup], (xcoldgass['LOGMH2_LIM']-xcoldgass['LOGMH1'])[xup], c='Black', alpha=0.25, marker=r'$\downarrow$')

    x = src['LOGMSTAR_MAGPHYS']
    y = src['LOGMH2'] - src['LOGMH1']

    norm = src.index[(src['H1_F']==1) & (src['H2_F']==1)].tolist()
    up = src.index[(src['H1_F']==2)].tolist()
    down = src.index[(src['H2_F']==2)].tolist()

    plt.scatter(x[up], y[up], c=src['CC'][up], marker=r'$\uparrow$')
    plt.scatter(x[down], y[down], c=src['CC'][down], marker=r'$\downarrow$')
    plt.scatter(x[norm], y[norm], c=src['CC'][norm], edgecolors= "black")

    median(x, y)

    plt.colorbar(label='MS dist')

    plt.ylabel('log $M_{H2}$/$M_{HI}$')      
    plt.xlabel('log $M_{*}$ [$M_{\odot}$]')

    if seven == True:
        plt.scatter(x[seven_odds()], y[seven_odds()], c='Black', marker='x') 

    plt.show()

def Catinella_2(src, xcg = True, seven = True):
    src = dist(src)

    color_map = matplotlib.cm.get_cmap('rainbow_r')
    plt.set_cmap(color_map)

    if xcg == True:
        xcoldgass = CF.src6.copy()

        xn = xcoldgass.index[(xcoldgass['H1_FLAG']!=99) & (xcoldgass['LOGMH2'].notnull())].tolist()
        xleft = xcoldgass.index[(xcoldgass['H1_FLAG']==99) & (xcoldgass['LOGMH2'].notnull())].tolist()
        xup = xcoldgass.index[(xcoldgass['LOGMH2_LIM'].notnull())].tolist()
        xs = xcoldgass.index[(xcoldgass['H1_FLAG']==99) & (xcoldgass['LOGMH2_LIM'].notnull())].tolist()

        plt.scatter((xcoldgass['LOGMH2']-xcoldgass['LOGMSTAR'])[xn], (xcoldgass['LOGMH1']-xcoldgass['LOGMSTAR'])[xn], c='Black', alpha=0.25, marker='o')
        plt.scatter((xcoldgass['LOGMH2']-xcoldgass['LOGMSTAR'])[xleft], (xcoldgass['LOGMH1']-xcoldgass['LOGMSTAR'])[xleft], c='Black', alpha=0.25, marker=r'$\leftarrow$')
        plt.scatter((xcoldgass['LOGMH2_LIM']-xcoldgass['LOGMSTAR'])[xup], (xcoldgass['LOGMH1']-xcoldgass['LOGMSTAR'])[xup], c='Black', alpha=0.25, marker=r'$\downarrow$')
        plt.scatter((xcoldgass['LOGMH2_LIM']-xcoldgass['LOGMSTAR'])[xs], (xcoldgass['LOGMH1']-xcoldgass['LOGMSTAR'])[xs], c='Black', alpha=0.25, marker='x')

    x = src['LOGMH2'] - src['LOGMSTAR_MAGPHYS']
    y = src['LOGMH1'] - src['LOGMSTAR_MAGPHYS']

    norm = src.index[(src['H1_F']==1) & (src['H2_F']==1)].tolist()
    left = src.index[(src['H1_F']==2)].tolist()
    down = src.index[(src['H2_F']==2)].tolist()

    plt.scatter(x[left], y[left], c=src['CC'][left], marker=r'$\leftarrow$')
    plt.scatter(x[down], y[down], c=src['CC'][down], marker=r'$\downarrow$')
    plt.scatter(x[norm], y[norm], c=src['CC'][norm], edgecolors= "black")

    plt.colorbar(label='MS dist')

    plt.ylabel('log $M_{H2}$/$M_{*}$')      
    plt.xlabel('log $M_{HI}$/$M_{*}$')

    if seven == True:
        plt.scatter(x[seven_odds()], y[seven_odds()], c='Black', marker='x') 

    plt.show()

def Catinella_3(src, xcg = True, seven = True):
    src = dist(src)

    color_map = matplotlib.cm.get_cmap('rainbow_r')
    plt.set_cmap(color_map)

    if xcg == True:
        xcoldgass = CF.src6.copy()

        xn = xcoldgass.index[(xcoldgass['H1_FLAG']!=99) & (xcoldgass['LOGMH2'].notnull())].tolist()
        xup = xcoldgass.index[(xcoldgass['H1_FLAG']==99) & (xcoldgass['LOGMH2_LIM'].notnull())].tolist()

        Mg = np.log10(1.3*(np.power(10,xcoldgass['LOGMH2'])+np.power(10,xcoldgass['LOGMH1'])))
        Mgup = np.log10(1.3*(np.power(10,xcoldgass['LOGMH2_LIM'])+np.power(10,xcoldgass['LOGMH1'])))

        plt.scatter((xcoldgass['LOGMSTAR'])[xn], (Mg-xcoldgass['LOGMSTAR'])[xn], c='Black', alpha=0.25, marker='o')
        plt.scatter((xcoldgass['LOGMSTAR'])[xup], (Mgup-xcoldgass['LOGMSTAR'])[xup], c='Black', alpha=0.25, marker=r'$\downarrow$')

    y = src['LOGMGAS'] - src['LOGMSTAR_MAGPHYS']
    x = src['LOGMSTAR_MAGPHYS']

    norm = src.index[(src['H1_F']==1) & (src['H2_F']==1)].tolist()
    down = src.index[(src['H1_F']==2) | (src['H2_F']==2)].tolist()

    plt.scatter(x[down], y[down], c=src['CC'][down], marker=r'$\downarrow$')
    plt.scatter(x[norm], y[norm], c=src['CC'][norm], edgecolors= "black")

    median(x, y)

    plt.colorbar(label='MS dist')
    
    plt.ylabel('log $M_{gas}$/$M_{*}$')     
    plt.xlabel('log $M_{*}$ [$M_{\odot}$]')

    if seven == True:
        plt.scatter(x[seven_odds()], y[seven_odds()], c='Black', marker='x') 

    plt.show()

    if xcg == True:
        xcoldgass = CF.src6.copy()

        xn = xcoldgass.index[(xcoldgass['H1_FLAG']!=99) & (xcoldgass['LOGMH2'].notnull())].tolist()
        xup = xcoldgass.index[(xcoldgass['H1_FLAG']==99) & (xcoldgass['LOGMH2_LIM'].notnull())].tolist()

        Mg = np.log10(1.3*(np.power(10,xcoldgass['LOGMH2'])+np.power(10,xcoldgass['LOGMH1'])))
        Mgup = np.log10(1.3*(np.power(10,xcoldgass['LOGMH2_LIM'])+np.power(10,xcoldgass['LOGMH1'])))

        plt.scatter((xcoldgass['LOGSFR'] - xcoldgass['LOGMSTAR'])[xn], (Mg-xcoldgass['LOGMSTAR'])[xn], c='Black', alpha=0.25, marker='o')
        plt.scatter((xcoldgass['LOGSFR'] - xcoldgass['LOGMSTAR'])[xup], (Mgup-xcoldgass['LOGMSTAR'])[xup], c='Black', alpha=0.25, marker=r'$\downarrow$')

    y = src['LOGMGAS'] - src['LOGMSTAR_MAGPHYS']
    x = src['LOGSFR_MAGPHYS'] - src['LOGMSTAR_MAGPHYS']

    norm = src.index[(src['H1_F']==1) & (src['H2_F']==1)].tolist()
    down = src.index[(src['H1_F']==2) | (src['H2_F']==2)].tolist()

    plt.scatter(x[down], y[down], c=src['CC'][down], marker=r'$\downarrow$')
    plt.scatter(x[norm], y[norm], c=src['CC'][norm], edgecolors= "black")

    median(x, y)

    plt.colorbar(label='MS dist')
    
    plt.ylabel('log $M_{gas}$/$M_{*}$')      
    plt.xlabel('log $sSFR$ [$yr^{-1}$]')

    if seven == True:
        plt.scatter(x[seven_odds()], y[seven_odds()], c='Black', marker='x') 

    plt.show()

def Catinella_4(src, xcg = True, seven = True):
    src['CC'] = src['LOGMH2'] - src['LOGMH1']

    color_map = matplotlib.cm.get_cmap('rainbow')
    plt.set_cmap(color_map)

    if xcg == True:
        xcoldgass = CF.src6.copy()

        xn = xcoldgass.index[(xcoldgass['H1_FLAG']!=99) & (xcoldgass['LOGMH2'].notnull())].tolist()
        xup = xcoldgass.index[(xcoldgass['H1_FLAG']==99) & (xcoldgass['LOGMH2_LIM'].notnull())].tolist()

        Mg = np.log10(1.3*(np.power(10,xcoldgass['LOGMH2'])+np.power(10,xcoldgass['LOGMH1'])))
        Mgup = np.log10(1.3*(np.power(10,xcoldgass['LOGMH2_LIM'])+np.power(10,xcoldgass['LOGMH1'])))

        plt.scatter((xcoldgass['LOGMSTAR'])[xn], (Mg-xcoldgass['LOGMSTAR'])[xn], c='Black', alpha=0.25, marker='o')
        plt.scatter((xcoldgass['LOGMSTAR'])[xup], (Mgup-xcoldgass['LOGMSTAR'])[xup], c='Black', alpha=0.25, marker=r'$\downarrow$')

    y = src['LOGMGAS'] - src['LOGMSTAR_MAGPHYS']
    x = src['LOGMSTAR_MAGPHYS']

    norm = src.index[(src['H1_F']==1) & (src['H2_F']==1)].tolist()
    down = src.index[(src['H1_F']==2) | (src['H2_F']==2)].tolist()

    plt.scatter(x[down], y[down], c=src['CC'][down], marker=r'$\downarrow$')
    plt.scatter(x[norm], y[norm], c=src['CC'][norm], edgecolors= "black")

    median(x, y)

    plt.colorbar(label='log Rmol')
    
    plt.ylabel('log $M_{gas}$/$M_{*}$')     
    plt.xlabel('log $M_{*}$ [$M_{\odot}$]')

    if seven == True:
        plt.scatter(x[seven_odds()], y[seven_odds()], c='Black', marker='x') 

    plt.show()

    if xcg == True:
        xcoldgass = CF.src6.copy()

        xn = xcoldgass.index[(xcoldgass['H1_FLAG']!=99) & (xcoldgass['LOGMH2'].notnull())].tolist()
        xup = xcoldgass.index[(xcoldgass['H1_FLAG']==99) & (xcoldgass['LOGMH2_LIM'].notnull())].tolist()

        Mg = np.log10(1.3*(np.power(10,xcoldgass['LOGMH2'])+np.power(10,xcoldgass['LOGMH1'])))
        Mgup = np.log10(1.3*(np.power(10,xcoldgass['LOGMH2_LIM'])+np.power(10,xcoldgass['LOGMH1'])))

        plt.scatter((xcoldgass['LOGSFR'] - xcoldgass['LOGMSTAR'])[xn], (Mg-xcoldgass['LOGMSTAR'])[xn], c='Black', alpha=0.25, marker='o')
        plt.scatter((xcoldgass['LOGSFR'] - xcoldgass['LOGMSTAR'])[xup], (Mgup-xcoldgass['LOGMSTAR'])[xup], c='Black', alpha=0.25, marker=r'$\downarrow$')

    y = src['LOGMGAS'] - src['LOGMSTAR_MAGPHYS']
    x = src['LOGSFR_MAGPHYS'] - src['LOGMSTAR_MAGPHYS']

    norm = src.index[(src['H1_F']==1) & (src['H2_F']==1)].tolist()
    down = src.index[(src['H1_F']==2) | (src['H2_F']==2)].tolist()

    plt.scatter(x[down], y[down], c=src['CC'][down], marker=r'$\downarrow$')
    plt.scatter(x[norm], y[norm], c=src['CC'][norm], edgecolors= "black")

    median(x, y)

    plt.colorbar(label='log Rmol')
    
    plt.ylabel('log $M_{gas}$/$M_{*}$')      
    plt.xlabel('log $sSFR$ [$yr^{-1}$]') 

    if seven == True:
        plt.scatter(x[seven_odds()], y[seven_odds()], c='Black', marker='x')

    plt.show()

def AGN_classes():
    JINGLE = GD.get_data(r"C:\Users\jsmes\OneDrive\Documents\Temp\Thesis\Data\JINGLE\JINGLE_FULL.csv", keep=True)

    print(JINGLE['AGNCLASS'][seven_odds()])

def Catinella_5(src, seven = True):
    #src = dist(src)
    src['CC'] = src['LOGMH2'] - src['LOGMH1']

    color_map = matplotlib.cm.get_cmap('rainbow')
    plt.set_cmap(color_map)

    y = src['LOGMDUST'] - src['LOGMSTAR_MAGPHYS']
    x = src['LOGMSTAR_MAGPHYS']

    plt.scatter(x, y, c=src['CC'], edgecolors= "black")

    median(x, y)

    plt.colorbar(label='log Rmol')
    
    plt.ylabel('log $M_{dust}$/$M_{*}$')     
    plt.xlabel('log $M_{*}$ [$M_{\odot}$]')

    if seven == True:
        plt.scatter(x[seven_odds()], y[seven_odds()], c='Black', marker='x') 

    plt.show()

    y = src['LOGMDUST'] - src['LOGMSTAR_MAGPHYS']
    x = src['LOGSFR_MAGPHYS'] - src['LOGMSTAR_MAGPHYS']

    plt.scatter(x, y, c=src['CC'], edgecolors= "black")

    median(x, y)

    plt.colorbar(label='log Rmol')
    
    plt.ylabel('log $M_{dust}$/$M_{*}$')      
    plt.xlabel('log $sSFR$ [$yr^{-1}$]')

    if seven == True:
        plt.scatter(x[seven_odds()], y[seven_odds()], c='Black', marker='x') 

    plt.show()

def ks_testing():
    src=MAIN.jingle

    data = src['LOGMGAS'] - src['LOGMSTAR_MAGPHYS']
    data = data.dropna()
    print(data)

    xcoldgass=MAIN.xcoldgass

    xn = xcoldgass.index[(xcoldgass['H1_FLAG']!=99) & (xcoldgass['LOGMH2'].notnull())].tolist()

    Mg = np.log10(1.3*(np.power(10,xcoldgass['LOGMH2'])+np.power(10,xcoldgass['LOGMH1'])))

    data2 = Mg - xcoldgass['LOGMSTAR'][xn]
    data2 = data2.dropna()
    print(data2)
    
    return ks_2samp(data, data2)