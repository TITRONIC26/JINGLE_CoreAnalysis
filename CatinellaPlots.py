#import libraries here
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
import numpy as np
import collections as coll

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

def dist2(src):
    points = [(src['LOGMSTAR'][m], (src['LOGSFR']-src['LOGMSTAR'])[m]) for m in range(len(src['LOGSFR']))]
    P1 = np.array([9,-9.822])
    P2 = np.array([12,-10.854])

    distance = [-(np.cross(P2 - P1, P1 - P3)/ np.linalg.norm(P2 - P1)) for P3 in points]
    
    src['CC'] = distance

    return src

def median_dot(x, y, num):
    nots = y.index[y.notnull()].tolist()
    
    x = x[nots]
    y = y[nots]

    total_bins = num+1
    print(total_bins)

    bins = np.linspace(x.min(), x.max(), total_bins)
    delta = bins[1]-bins[0]
    idx = np.digitize(x, bins)
    running_median = [np.median(y[idx==k]) for k in range(total_bins)]

    plt.scatter(bins-delta/2, running_median, color='Blue', edgecolors='Black', marker='D', s=100)

def mean_dot(x, y, num):
    nots = y.index[y.notnull()].tolist()
    
    x = x[nots]
    y = y[nots]

    total_bins = num+1
    print(total_bins)

    bins = np.linspace(x.min(), x.max(), total_bins)
    delta = bins[1]-bins[0]
    idx = np.digitize(x, bins)
    running_mean = [np.mean(y[idx==k]) for k in range(total_bins)]

    plt.scatter(bins-delta/2, running_mean, color='Red', edgecolors='Black', marker='o', s=125)

def median(x, y, num):
    nots = y.index[y.notnull()].tolist()
    
    x = x[nots]
    y = y[nots]

    total_bins = int(len(y)/num)
    print(total_bins)

    bins = np.linspace(x.min(), x.max(), total_bins)
    delta = bins[1]-bins[0]
    idx = np.digitize(x, bins)
    running_median = [np.median(y[idx==k]) for k in range(total_bins)]

    plt.plot(bins-delta/2, running_median, color='White', linewidth=5)
    plt.plot(bins-delta/2, running_median, color='Black', linewidth=3)

def Figure1(src, xcg=True):
    if xcg == True:
        xcoldgass = CF.src6.copy()

        plt.scatter(0.0667*xcoldgass['RA'], xcoldgass['DEC'], c='Blue', alpha=0.75, marker='o', edgecolors= "blue", label='xCOLDGASS')

    plt.scatter(0.0667*src['RA'], src['DEC'], c='DarkOrange', alpha=0.75, marker='o', edgecolors= "DarkOrange", label='JINGLE')

    plt.xlabel('RA [hr]')
    plt.ylabel('DEC [deg]')
    plt.title('JINGLE & xCOLDGASS sky')

    plt.legend()

    plt.show()

def Figure5(src, xcg=True):
    if xcg == True:
        xcoldgass = CF.src6.copy()

        x = xcoldgass['LOGMSTAR']
        y = xcoldgass['LOGMH1']-xcoldgass['LOGMSTAR']

        norm = xcoldgass.index[(xcoldgass['H1_FLAG']!=99)].tolist()
        xup = xcoldgass.index[(xcoldgass['H1_FLAG']==99)].tolist()

        plt.scatter(x[norm], y[norm], c='Blue', alpha=0.75, marker='o', edgecolors= "blue")
        plt.scatter(x[xup], y[xup], c='SkyBlue', alpha=0.75, marker=r'$\downarrow$')
    
    x = src['LOGMSTAR_MAGPHYS']
    y = src['LOGMH1']-src['LOGMSTAR_MAGPHYS']

    norm = src.index[(src['H1_F']==1)].tolist()
    xup = src.index[(src['H1_F']==3)].tolist()
    xder = src.index[(src['H1_F']==2)].tolist()

    plt.scatter(x[norm], y[norm], c='DarkOrange', alpha=0.75, marker='o', edgecolors= "DarkOrange")
    plt.scatter(x[xup], y[xup], c='DarkOrange', alpha=0.75, marker=r'$\downarrow$')
    plt.scatter(x[xder], y[xder], c='DarkOrange', alpha=0.75, marker='x')

    plt.xlabel('log $M_{*}$ [$M_{\odot}$]')
    plt.ylabel('log $M_{HI}$/$M_{*}$')

    plt.show()

    #---------------------

    if xcg == True:
        xcoldgass = CF.src6.copy()

        x = xcoldgass['LOGSFR']-xcoldgass['LOGMSTAR']
        y = xcoldgass['LOGMH1']-xcoldgass['LOGMSTAR']

        norm = xcoldgass.index[(xcoldgass['H1_FLAG']!=99)].tolist()
        xup = xcoldgass.index[(xcoldgass['H1_FLAG']==99)].tolist()

        plt.scatter(x[norm], y[norm], c='Blue', alpha=0.75, marker='o', edgecolors= "blue")
        plt.scatter(x[xup], y[xup], c='SkyBlue', alpha=0.75, marker=r'$\downarrow$')
    
    x = src['LOGSFR_MAGPHYS']-src['LOGMSTAR_MAGPHYS']
    y = src['LOGMH1']-src['LOGMSTAR_MAGPHYS']

    norm = src.index[(src['H1_F']==1)].tolist()
    xup = src.index[(src['H1_F']==3)].tolist()
    xder = src.index[(src['H1_F']==2)].tolist()

    plt.scatter(x[norm], y[norm], c='DarkOrange', alpha=0.75, marker='o', edgecolors= "DarkOrange")
    plt.scatter(x[xup], y[xup], c='DarkOrange', alpha=0.75, marker=r'$\downarrow$')
    plt.scatter(x[xder], y[xder], c='DarkOrange', alpha=0.75, marker='x')

    plt.xlabel('log sSFR [$yr^{-1}$]')
    plt.ylabel('log $M_{HI}$/$M_{*}$')

    plt.show()

def Figure6(src, xcg=True):
    if xcg == True:
        xcoldgass = CF.src6.copy()

        x = xcoldgass['LOGMSTAR']
        y = xcoldgass['LOGMH1']-xcoldgass['LOGMSTAR']

        norm = xcoldgass.index[(xcoldgass['H1_FLAG']!=99)].tolist()
        xup = xcoldgass.index[(xcoldgass['H1_FLAG']==99)].tolist()

        plt.scatter(x[norm], y[norm], c='Grey', alpha=0.25, marker='o', edgecolors= "Grey")
        plt.scatter(x[xup], y[xup], c='Grey', alpha=0.25, marker=r'$\downarrow$')

        mean_dot(x,y,8)
        median_dot(x,y,8)
    else:
        x = src['LOGMSTAR_MAGPHYS']
        y = src['LOGMH1']-src['LOGMSTAR_MAGPHYS']

        norm = src.index[(src['H1_F']==1)].tolist()
        xup = src.index[(src['H1_F']==3)].tolist()
        xder = src.index[(src['H1_F']==2)].tolist()

        plt.scatter(x[norm], y[norm], c='Grey', alpha=0.25, marker='o', edgecolors= "Grey")
        plt.scatter(x[xup], y[xup], c='Grey', alpha=0.25, marker=r'$\downarrow$')
        plt.scatter(x[xder], y[xder], c='Grey', alpha=0.25, marker='x')

        mean_dot(x,y,4)
        median_dot(x,y,4)

    plt.xlabel('log $M_{*}$ [$M_{\odot}$]')
    plt.ylabel('log $M_{HI}$/$M_{*}$')

    plt.show()

    #---------------------

    if xcg == True:
        xcoldgass = CF.src6.copy()

        x = xcoldgass['LOGSFR']-xcoldgass['LOGMSTAR']
        y = xcoldgass['LOGMH1']-xcoldgass['LOGMSTAR']

        norm = xcoldgass.index[(xcoldgass['H1_FLAG']!=99)].tolist()
        xup = xcoldgass.index[(xcoldgass['H1_FLAG']==99)].tolist()

        plt.scatter(x[norm], y[norm], c='Grey', alpha=0.25, marker='o', edgecolors= "Grey")
        plt.scatter(x[xup], y[xup], c='Grey', alpha=0.25, marker=r'$\downarrow$')

        nots = x.index[x >= -13].tolist()

        mean_dot(x[nots],y[nots],5)
        median_dot(x[nots],y[nots],5)
    else:    
        x = src['LOGSFR_MAGPHYS']-src['LOGMSTAR_MAGPHYS']
        y = src['LOGMH1']-src['LOGMSTAR_MAGPHYS']

        norm = src.index[(src['H1_F']==1)].tolist()
        xup = src.index[(src['H1_F']==3)].tolist()
        xder = src.index[(src['H1_F']==2)].tolist()

        plt.scatter(x[norm], y[norm], c='Grey', alpha=0.25, marker='o', edgecolors= "Grey")
        plt.scatter(x[xup], y[xup], c='Grey', alpha=0.25, marker=r'$\downarrow$')
        plt.scatter(x[xder], y[xder], c='Grey', alpha=0.25, marker='x')

        mean_dot(x,y,3)
        median_dot(x,y,3)

    plt.xlabel('log sSFR [$yr^{-1}$]')
    plt.ylabel('log $M_{HI}$/$M_{*}$')

    plt.show()

def Figure7(src, xcg=True):
    if xcg == True:
        xcoldgass = CF.src6.copy()

        x = xcoldgass['LOGMSTAR']
        y = xcoldgass['LOGSFR']-xcoldgass['LOGMSTAR']

        h1_det = xcoldgass.index[(xcoldgass['H1_FLAG']!=99)&(xcoldgass['LOGMH2'].isnull())&(xcoldgass['LOGMH2_LIM'].isnull())].tolist()
        h1_non = xcoldgass.index[(xcoldgass['H1_FLAG']==99)].tolist()
        h2_det = xcoldgass.index[(xcoldgass['H1_FLAG']!=99)&(xcoldgass['LOGMH2'].notnull())].tolist()
        h2_non = xcoldgass.index[(xcoldgass['H1_FLAG']!=99)&(xcoldgass['LOGMH2_LIM'].notnull())].tolist()

        plt.scatter(x[h1_det], y[h1_det], c='SkyBlue', alpha=0.75, marker='o', edgecolors= "SkyBlue")
        plt.scatter(x[h2_det], y[h2_det], c='Blue', alpha=0.75, marker='o', edgecolors= "Blue")
        #plt.scatter(x[h1_non], y[h1_non], c='White', alpha=0.75, marker='o', edgecolors= "Blue")
        plt.scatter(x[h2_non], y[h2_non], c='Magenta', alpha=0.75, marker='+')
    else:
        y = src['LOGSFR_MAGPHYS']-src['LOGMSTAR_MAGPHYS']
        x = src['LOGMSTAR_MAGPHYS']

        h1_det = src.index[(src['H1_F']==1)&(src['H2_F']==0)].tolist()
        h1_non = src.index[(src['H1_F']==3)&(src['H2_F']!=1)].tolist()
        h1_der = src.index[(src['H1_F']==2)&(src['H2_F']!=1)].tolist()

        h2_det = src.index[(src['H1_F']==1)&(src['H2_F']==1)].tolist()
        h2_non = src.index[(src['H1_F']==3)&(src['H2_F']==1)].tolist()
        h2_der = src.index[(src['H1_F']==2)&(src['H2_F']==1)].tolist()

        pinks = src.index[(src['H1_F']==1)&(src['H2_F']==2)].tolist()

        plt.scatter(x[h1_det], y[h1_det], c='SkyBlue', alpha=0.75, marker='o', edgecolors= "SkyBlue")
        plt.scatter(x[h1_non], y[h1_non], c='SkyBlue', alpha=0.75, marker='v', edgecolors= "SkyBlue")
        plt.scatter(x[h1_der], y[h1_der], c='White', alpha=0.75, marker='o', edgecolors= "SkyBlue")

        plt.scatter(x[h2_det], y[h2_det], c='Blue', alpha=0.75, marker='o', edgecolors= "Blue")
        plt.scatter(x[h2_der], y[h2_der], c='White', alpha=0.75, marker='o', edgecolors= "Blue")
        plt.scatter(x[h2_non], y[h2_non], c='Blue', alpha=0.75, marker='v', edgecolors='Blue')
    
        plt.scatter(x[pinks], y[pinks], c='Magenta', alpha=0.75, marker='+')

    plt.ylabel('log sSFR [$yr^{-1}$]')
    plt.xlabel('log $M_{*}$ [$M_{\odot}$]')

    plt.show()

def Figure8(src, xcg=True):
    if xcg == True:
        xcoldgass = CF.src6.copy()

        xcoldgass['LOGMGAS'] = np.log10(1.3*(np.power(10,xcoldgass['LOGMH2'])+np.power(10,xcoldgass['LOGMH1'])))
        xcoldgass['LOGMGAS_UP'] = np.log10(1.3*(np.power(10,xcoldgass['LOGMH2_LIM'])+np.power(10,xcoldgass['LOGMH1'])))

        x = xcoldgass['LOGMSTAR']
        y = xcoldgass['LOGMGAS']-xcoldgass['LOGMSTAR']
        yup = xcoldgass['LOGMGAS_UP']-xcoldgass['LOGMSTAR']

        norm = xcoldgass.index[(xcoldgass['H1_FLAG']!=99)].tolist()
        xup = xcoldgass.index[(xcoldgass['H1_FLAG']==99)].tolist()

        plt.scatter(x[norm], y[norm], c='Grey', alpha=0.25, marker='o', edgecolors= "Grey")
        plt.scatter(x[norm], yup[norm], c='LightGrey', alpha=0.75, marker='o', edgecolors= "Black")
        plt.scatter(x[xup], y[xup], c='Grey', alpha=0.25, marker=r'$\downarrow$')
        plt.scatter(x[xup], yup[xup], c='Black', alpha=0.25, marker=r'$\downarrow$')

        mean_dot(x,y,8)
        median_dot(x,y,8)
    else:
        x = src['LOGMSTAR_MAGPHYS']
        y = src['LOGMGAS']-src['LOGMSTAR_MAGPHYS']

        h1_det = src.index[(src['H1_F']==1)&(src['H2_F']==0)].tolist()
        h1_non = src.index[(src['H1_F']==3)&(src['H2_F']!=1)].tolist()
        h1_only_der = src.index[(src['H1_F']==2)&(src['H2_F']!=1)].tolist()

        both = src.index[(src['H1_F']==1)&(src['H2_F']==1)].tolist()
        both_non = src.index[(src['H1_F']==3)&(src['H2_F']==1)].tolist()

        h1_der = src.index[(src['H1_F']==2)&(src['H2_F']==1)].tolist()
        h2_der = src.index[(src['H1_F']==1)&(src['H2_F']==2)].tolist()

        plt.scatter(x[both], y[both], c='Grey', alpha=0.25, marker='o', edgecolors= "Grey")
        plt.scatter(x[both_non], y[both_non], c='Grey', alpha=0.75, marker=r'$\downarrow$', edgecolors= "Grey")

        plt.scatter(x[h1_det], y[h1_det], c='LightGrey', alpha=0.75, marker='o', edgecolors= "Black")

        plt.scatter(x[h1_der], y[h1_der], c='Grey', alpha=0.25, marker='^', edgecolors= "Grey")
        plt.scatter(x[h2_der], y[h2_der], c='Grey', alpha=0.25, marker='^', edgecolors= "Grey")

        plt.scatter(x[h1_non], y[h1_non], c='LightGrey', alpha=0.75, marker=r'$\downarrow$', edgecolors= "Black")
        plt.scatter(x[h1_only_der], y[h1_only_der], c='LightGrey', alpha=0.75, marker='^', edgecolors= "Black")

        mean_dot(x,y,8)
        median_dot(x,y,8)

    plt.xlabel('log $M_{*}$ [$M_{\odot}$]')
    plt.ylabel('log $M_{gas}$/$M_{*}$')

    plt.show()

    #---------------------

    if xcg == True:
        xcoldgass = CF.src6.copy()

        xcoldgass['LOGMGAS'] = np.log10(1.3*(np.power(10,xcoldgass['LOGMH2'])+np.power(10,xcoldgass['LOGMH1'])))
        xcoldgass['LOGMGAS_UP'] = np.log10(1.3*(np.power(10,xcoldgass['LOGMH2_LIM'])+np.power(10,xcoldgass['LOGMH1'])))


        x = xcoldgass['LOGSFR']-xcoldgass['LOGMSTAR']
        y = xcoldgass['LOGMGAS']-xcoldgass['LOGMSTAR']
        yup = xcoldgass['LOGMGAS_UP']-xcoldgass['LOGMSTAR']

        norm = xcoldgass.index[(xcoldgass['H1_FLAG']!=99)].tolist()
        xup = xcoldgass.index[(xcoldgass['H1_FLAG']==99)].tolist()

        plt.scatter(x[norm], y[norm], c='Grey', alpha=0.25, marker='o', edgecolors= "Grey")
        plt.scatter(x[norm], yup[norm], c='LightGrey', alpha=0.75, marker='o', edgecolors= "Black")
        plt.scatter(x[xup], y[xup], c='Grey', alpha=0.25, marker=r'$\downarrow$')
        plt.scatter(x[xup], yup[xup], c='Black', alpha=0.25, marker=r'$\downarrow$')

        nots = x.index[x >= -13].tolist()

        mean_dot(x[nots],y[nots],5)
        median_dot(x[nots],y[nots],5)
    else:    
        x = src['LOGSFR_MAGPHYS']-src['LOGMSTAR_MAGPHYS']
        y = src['LOGMGAS']-src['LOGMSTAR_MAGPHYS']

        h1_det = src.index[(src['H1_F']==1)&(src['H2_F']==0)].tolist()
        h1_non = src.index[(src['H1_F']==3)&(src['H2_F']!=1)].tolist()
        h1_only_der = src.index[(src['H1_F']==2)&(src['H2_F']!=1)].tolist()

        both = src.index[(src['H1_F']==1)&(src['H2_F']==1)].tolist()
        both_non = src.index[(src['H1_F']==3)&(src['H2_F']==1)].tolist()

        h1_der = src.index[(src['H1_F']==2)&(src['H2_F']==1)].tolist()
        h2_der = src.index[(src['H1_F']==1)&(src['H2_F']==2)].tolist()

        plt.scatter(x[both], y[both], c='Grey', alpha=0.25, marker='o', edgecolors= "Grey")
        plt.scatter(x[both_non], y[both_non], c='Grey', alpha=0.75, marker=r'$\downarrow$', edgecolors= "Grey")
        
        plt.scatter(x[h1_det], y[h1_det], c='LightGrey', alpha=0.75, marker='o', edgecolors= "Black")

        plt.scatter(x[h1_der], y[h1_der], c='Grey', alpha=0.25, marker='^', edgecolors= "Grey")
        plt.scatter(x[h2_der], y[h2_der], c='Grey', alpha=0.25, marker='^', edgecolors= "Grey")

        plt.scatter(x[h1_non], y[h1_non], c='LightGrey', alpha=0.75, marker=r'$\downarrow$', edgecolors= "Black")
        plt.scatter(x[h1_only_der], y[h1_only_der], c='LightGrey', alpha=0.75, marker='^', edgecolors= "Black")

        mean_dot(x,y,5)
        median_dot(x,y,5)

    plt.xlabel('log sSFR [$yr^{-1}$]')
    plt.ylabel('log $M_{gas}$/$M_{*}$')

    plt.show()

def Figure9(src, xcg=True):
    label = 'MS dist'
    cc = 'jet_r'

    color_map = matplotlib.cm.get_cmap(cc)
    plt.set_cmap(color_map)

    if xcg == True:
        xcoldgass = CF.src6.copy()

        xcoldgass = dist2(xcoldgass)

        y = xcoldgass['LOGMH2'] - xcoldgass['LOGMH1']
        yup = xcoldgass['LOGMH2_LIM'] - xcoldgass['LOGMH1']
        x = xcoldgass['LOGMSTAR']

        dot = xcoldgass.index[(xcoldgass['H1_FLAG']!=99) & (xcoldgass['LOGMH2'].notnull())].tolist()
        down = xcoldgass.index[(xcoldgass['H1_FLAG']!=99) & (xcoldgass['LOGMH2_LIM'].notnull())].tolist()
        up = xcoldgass.index[(xcoldgass['H1_FLAG']==99) & (xcoldgass['LOGMH2'].notnull())].tolist()
        hollow = xcoldgass.index[(xcoldgass['H1_FLAG']==99) & (xcoldgass['LOGMH2_LIM'].notnull())].tolist()

        plt.scatter(x[dot], y[dot], c=xcoldgass['CC'][dot], alpha=0.75, marker='o', edgecolors='Black')
        plt.scatter(x[hollow], yup[hollow], c='White', alpha=0.75, marker='o', edgecolors='Black')
        plt.scatter(x[up], y[up], c=xcoldgass['CC'][up], alpha=0.75, marker=r'$\uparrow$')
        plt.scatter(x[down], yup[down], c=xcoldgass['CC'][down], alpha=0.75, marker=r'$\downarrow$')

        median(x, y, 35)

    else:
        src = dist(src)

        y = src['LOGMH2'] - src['LOGMH1']
        x = src['LOGMSTAR_MAGPHYS']

        dot = src.index[(src['H1_F']==1) & (src['H2_F']==1)].tolist()
        up = src.index[(src['H1_F']==3) & (src['H2_F']==1)].tolist()
        
        hollow_up = src.index[(src['H1_F']==3) & (src['H2_F']==2)].tolist()

        tup = src.index[(src['H1_F']==1) & (src['H2_F']==2)].tolist()
        tdown = src.index[(src['H1_F']==2) & (src['H2_F']==1)].tolist()

        plt.scatter(x[dot], y[dot], c=src['CC'][dot], alpha=0.75, marker='o', edgecolors='Black')
        plt.scatter(x[up], y[up], c=src['CC'][up], alpha=0.75, marker=r'$\uparrow$')
        
        plt.scatter(x[hollow_up], y[hollow_up], c='White', alpha=0.75, marker='o', edgecolors='Black')
        plt.scatter(x[tup], y[tup], c=src['CC'][tup], alpha=0.75, marker='^', edgecolors='Black')
        plt.scatter(x[tdown], y[tdown], c=src['CC'][tdown], alpha=0.75, marker='v', edgecolors='Black')
        
        median(x, y, 15)
    
    plt.colorbar(label=label)

    plt.ylabel('log $M_{H2}$/$M_{HI}$')      
    plt.xlabel('log $M_{*}$ [$M_{\odot}$]')

    plt.show()

def Figure10(src, xcg=True):
    label = 'MS dist'
    cc = 'jet_r'

    color_map = matplotlib.cm.get_cmap(cc)
    plt.set_cmap(color_map)

    if xcg == True:
        xcoldgass = CF.src6.copy()

        xcoldgass = dist2(xcoldgass)

        y = xcoldgass['LOGMH2'] - xcoldgass['LOGMSTAR']
        yup = xcoldgass['LOGMH2_LIM'] - xcoldgass['LOGMSTAR']
        x = xcoldgass['LOGMH1'] - xcoldgass['LOGMSTAR']

        dot = xcoldgass.index[(xcoldgass['H1_FLAG']!=99) & (xcoldgass['LOGMH2'].notnull())].tolist()
        H2up = xcoldgass.index[(xcoldgass['H1_FLAG']!=99) & (xcoldgass['LOGMH2_LIM'].notnull())].tolist()
        H1up = xcoldgass.index[(xcoldgass['H1_FLAG']==99) & (xcoldgass['LOGMH2'].notnull())].tolist()
        both = xcoldgass.index[(xcoldgass['H1_FLAG']==99) & (xcoldgass['LOGMH2_LIM'].notnull())].tolist()

        plt.scatter(x[dot], y[dot], c=xcoldgass['CC'][dot], alpha=0.75, marker='o', edgecolors='Black')
        plt.scatter(x[both], yup[both], c='White', alpha=0.75, marker='o', edgecolors='Black')
        plt.scatter(x[H1up], y[H1up], c=xcoldgass['CC'][H1up], alpha=0.75, marker=r'$\leftarrow$')
        plt.scatter(x[H2up], yup[H2up], c=xcoldgass['CC'][H2up], alpha=0.75, marker=r'$\downarrow$')

    else:
        src = dist(src)

        y = src['LOGMH2'] - src['LOGMSTAR_MAGPHYS']
        x = src['LOGMH1'] - src['LOGMSTAR_MAGPHYS']

        dot = src.index[(src['H1_F']==1) & (src['H2_F']==1)].tolist()
        H1up = src.index[(src['H1_F']==3) & (src['H2_F']==1)].tolist()
        
        hollow_up = src.index[(src['H1_F']==3) & (src['H2_F']==2)].tolist()
        
        tup = src.index[(src['H1_F']==1) & (src['H2_F']==2)].tolist()
        tdown = src.index[(src['H1_F']==2) & (src['H2_F']==1)].tolist()

        plt.scatter(x[dot], y[dot], c=src['CC'][dot], alpha=0.75, marker='o', edgecolors='Black')
        plt.scatter(x[H1up], y[H1up], c=src['CC'][H1up], alpha=0.75, marker=r'$\leftarrow$')
        plt.scatter(x[hollow_up], y[hollow_up], c='White', alpha=0.75, marker='o', edgecolors='Black')

        plt.scatter(x[tup], y[tup], c=src['CC'][tup], alpha=0.75, marker='^', edgecolors='Black')
        plt.scatter(x[tdown], y[tdown], c=src['CC'][tdown], alpha=0.75, marker='>', edgecolors='Black')
    
    plt.colorbar(label=label)

    plt.ylabel('log $M_{H2}$/$M_{*}$')      
    plt.xlabel('log $M_{HI}$/$M_{*}$')

    plt.show()

def Figure11(src, xcg=True):
    label = 'log Rmol'
    cc = 'jet'

    color_map = matplotlib.cm.get_cmap(cc)
    plt.set_cmap(color_map)

    if xcg == True:
        xcoldgass = CF.src6.copy()

        xcoldgass['CC'] = xcoldgass['LOGMH2'] - xcoldgass['LOGMH1']
        xcoldgass.loc[xcoldgass['CC'].isnull(), 'CC'] = xcoldgass['LOGMH2_LIM'] - xcoldgass['LOGMH1']

        xcoldgass['LOGMGAS'] = np.log10(1.3*(np.power(10,xcoldgass['LOGMH2'])+np.power(10,xcoldgass['LOGMH1'])))
        xcoldgass['LOGMGAS_UP'] = np.log10(1.3*(np.power(10,xcoldgass['LOGMH2_LIM'])+np.power(10,xcoldgass['LOGMH1'])))

        y = xcoldgass['LOGMGAS'] - xcoldgass['LOGMSTAR']
        yup = xcoldgass['LOGMGAS_UP'] - xcoldgass['LOGMSTAR']
        x = xcoldgass['LOGMSTAR']

        dot = xcoldgass.index[(xcoldgass['H1_FLAG']!=99) & (xcoldgass['LOGMH2'].notnull())].tolist()
        H2up = xcoldgass.index[(xcoldgass['H1_FLAG']!=99) & (xcoldgass['LOGMH2_LIM'].notnull())].tolist()
        H1up = xcoldgass.index[(xcoldgass['H1_FLAG']==99) & (xcoldgass['LOGMH2'].notnull())].tolist()
        both = xcoldgass.index[(xcoldgass['H1_FLAG']==99) & (xcoldgass['LOGMH2_LIM'].notnull())].tolist()

        plt.scatter(x[dot], y[dot], c=xcoldgass['CC'][dot], alpha=0.75, marker='o', edgecolors='Black')
        plt.scatter(x[both], yup[both], c='Grey', alpha=0.75, marker=r'$\downarrow$')
        plt.scatter(x[H1up], y[H1up], c=xcoldgass['CC'][H1up], alpha=0.75, marker=r'$\downarrow$')
        plt.scatter(x[H2up], yup[H2up], c=xcoldgass['CC'][H2up], alpha=0.75, marker='v', edgecolors='Black')

        xcoldgass.loc[xcoldgass['LOGMGAS'].isnull(), 'LOGMGAS'] = xcoldgass['LOGMGAS_UP']
        y = xcoldgass['LOGMGAS'] - xcoldgass['LOGMSTAR']

        median(x,y,55)

    else:
        src['CC'] = src['LOGMH2'] - src['LOGMH1']

        y = src['LOGMGAS'] - src['LOGMSTAR_MAGPHYS']
        x = src['LOGMSTAR_MAGPHYS']

        dot = src.index[(src['H1_F']==1) & (src['H2_F']==1)].tolist()
        H1up = src.index[(src['H1_F']==3) & (src['H2_F']==1)].tolist()
        hollow_up = src.index[(src['H1_F']==3) & (src['H2_F']==2)].tolist()
        tup = src.index[(src['H1_F']==1) & (src['H2_F']==2)].tolist()
        tdown = src.index[(src['H1_F']==2) & (src['H2_F']==1)].tolist()

        plt.scatter(x[dot], y[dot], c=src['CC'][dot], alpha=0.75, marker='o', edgecolors='Black')
        plt.scatter(x[H1up], y[H1up], c=src['CC'][H1up], alpha=0.75, marker=r'$\downarrow$')

        plt.scatter(x[hollow_up], y[hollow_up], c='White', alpha=0.75, marker='o', edgecolors='Black')

        plt.scatter(x[tup], y[tup], c=src['CC'][tup], alpha=0.75, marker='^', edgecolors='Red')
        plt.scatter(x[tdown], y[tdown], c=src['CC'][tdown], alpha=0.75, marker='^', edgecolors='Blue')

        median(x,y,20)
    
    plt.colorbar(label=label)

    plt.ylabel('log $M_{gas}$/$M_{*}$')      
    plt.xlabel('log $M_{*}$ [$M_{\odot}$]')

    plt.show()

    #------------------------------------------

    if xcg == True:
        xcoldgass = CF.src6.copy()

        xcoldgass['CC'] = xcoldgass['LOGMH2'] - xcoldgass['LOGMH1']
        xcoldgass.loc[xcoldgass['CC'].isnull(), 'CC'] = xcoldgass['LOGMH2_LIM'] - xcoldgass['LOGMH1']

        xcoldgass['LOGMGAS'] = np.log10(1.3*(np.power(10,xcoldgass['LOGMH2'])+np.power(10,xcoldgass['LOGMH1'])))
        xcoldgass['LOGMGAS_UP'] = np.log10(1.3*(np.power(10,xcoldgass['LOGMH2_LIM'])+np.power(10,xcoldgass['LOGMH1'])))

        y = xcoldgass['LOGMGAS'] - xcoldgass['LOGMSTAR']
        yup = xcoldgass['LOGMGAS_UP'] - xcoldgass['LOGMSTAR']
        x = xcoldgass['LOGSFR'] - xcoldgass['LOGMSTAR']

        dot = xcoldgass.index[(xcoldgass['H1_FLAG']!=99) & (xcoldgass['LOGMH2'].notnull())].tolist()
        H2up = xcoldgass.index[(xcoldgass['H1_FLAG']!=99) & (xcoldgass['LOGMH2_LIM'].notnull())].tolist()
        H1up = xcoldgass.index[(xcoldgass['H1_FLAG']==99) & (xcoldgass['LOGMH2'].notnull())].tolist()
        both = xcoldgass.index[(xcoldgass['H1_FLAG']==99) & (xcoldgass['LOGMH2_LIM'].notnull())].tolist()

        plt.scatter(x[dot], y[dot], c=xcoldgass['CC'][dot], alpha=0.75, marker='o', edgecolors='Black')
        plt.scatter(x[both], yup[both], c='Grey', alpha=0.75, marker=r'$\downarrow$')
        plt.scatter(x[H1up], y[H1up], c=xcoldgass['CC'][H1up], alpha=0.75, marker=r'$\downarrow$')
        plt.scatter(x[H2up], yup[H2up], c=xcoldgass['CC'][H2up], alpha=0.75, marker='v', edgecolors='Black')

        xcoldgass.loc[xcoldgass['LOGMGAS'].isnull(), 'LOGMGAS'] = xcoldgass['LOGMGAS_UP']
        y = xcoldgass['LOGMGAS'] - xcoldgass['LOGMSTAR']

        nots = x.index[x >= -13].tolist()

        median(x[nots],y[nots],55)

    else:
        src['CC'] = src['LOGMH2'] - src['LOGMH1']

        y = src['LOGMGAS'] - src['LOGMSTAR_MAGPHYS']
        x = src['LOGSFR_MAGPHYS'] - src['LOGMSTAR_MAGPHYS']

        dot = src.index[(src['H1_F']==1) & (src['H2_F']==1)].tolist()
        H1up = src.index[(src['H1_F']==3) & (src['H2_F']==1)].tolist()
        hollow_up = src.index[(src['H1_F']==3) & (src['H2_F']==2)].tolist()
        tup = src.index[(src['H1_F']==1) & (src['H2_F']==2)].tolist()
        tdown = src.index[(src['H1_F']==2) & (src['H2_F']==1)].tolist()

        plt.scatter(x[dot], y[dot], c=src['CC'][dot], alpha=0.75, marker='o', edgecolors='Black')
        plt.scatter(x[H1up], y[H1up], c=src['CC'][H1up], alpha=0.75, marker=r'$\downarrow$')
        plt.scatter(x[hollow_up], y[hollow_up], c='White', alpha=0.75, marker='o', edgecolors='Black')

        plt.scatter(x[tup], y[tup], c=src['CC'][tup], alpha=0.75, marker='^', edgecolors='Red')
        plt.scatter(x[tdown], y[tdown], c=src['CC'][tdown], alpha=0.75, marker='^', edgecolors='Blue')

        median(x,y,20)
    
    plt.colorbar(label=label)

    plt.ylabel('log $M_{gas}$/$M_{*}$')      
    plt.xlabel('log sSFR [$y^{-1}$]')

    plt.show()

def New_Figure1(src, xcg=True, selection=True):
    if xcg == True:
        xcoldgass = CF.src6.copy()

        if selection == True:
            xcoldgass['LOGMGAS'] = np.log10(1.3*(np.power(10,xcoldgass['LOGMH2'])+np.power(10,xcoldgass['LOGMH1'])))
            xcoldgass['LOGMGAS_UP'] = np.log10(1.3*(np.power(10,xcoldgass['LOGMH2_LIM'])+np.power(10,xcoldgass['LOGMH1'])))

            ind = xcoldgass.index[(xcoldgass['LOGSFR'] - xcoldgass['LOGMSTAR'] >= -11) & ((xcoldgass['LOGMGAS'] - xcoldgass['LOGMSTAR'] >= -0.5) | (xcoldgass['LOGMGAS_UP'] - xcoldgass['LOGMSTAR'] >= -0.5))].tolist()
            #ind = xcoldgass.index[xcoldgass['LOGSFR'] >= -1].tolist()
            print(len(ind))
    
            x = xcoldgass['LOGMSTAR'][ind]
            y = xcoldgass['LOGSFR'][ind]
        else:
            x = xcoldgass['LOGMSTAR']
            y = xcoldgass['LOGSFR']    

        plt.scatter(x, y, c='Blue', alpha=0.75, marker='o', edgecolors= "Blue", label='xCOLD GASS')
    
    if selection == True:
        ind = src.index[(src['LOGSFR_MAGPHYS'] - src['LOGMSTAR_MAGPHYS'] >= -11) & (src['LOGMGAS'] - src['LOGMSTAR_MAGPHYS'] >= -0.5)].tolist()
        #ind = src.index[src['LOGSFR_MAGPHYS'] >= -1].tolist()
        print(len(ind))

        x = src['LOGMSTAR_MAGPHYS'][ind]
        y = src['LOGSFR_MAGPHYS'][ind]

    else:
        x = src['LOGMSTAR_MAGPHYS']
        y = src['LOGSFR_MAGPHYS']

    plt.scatter(x, y, c='DarkOrange', alpha=0.75, marker='o', edgecolors= "DarkOrange", label='JINGLE')

    plt.xlabel('log $M_{*}$ [$M_{\odot}$]')
    plt.ylabel('log $SFR$ [$M_{\odot} y^{-1}$]')

    plt.legend()

    plt.show()

def New_Figure2(src, xcg=True):
    if xcg == True:
        xcoldgass = CF.src6.copy()

        xcoldgass['LOGMGAS'] = np.log10(1.3*(np.power(10,xcoldgass['LOGMH2'])+np.power(10,xcoldgass['LOGMH1'])))
        xcoldgass['LOGMGAS_UP'] = np.log10(1.3*(np.power(10,xcoldgass['LOGMH2_LIM'])+np.power(10,xcoldgass['LOGMH1'])))

        y = xcoldgass['LOGMGAS'] - xcoldgass['LOGMSTAR']
        yup = xcoldgass['LOGMGAS_UP'] - xcoldgass['LOGMSTAR']
        x = xcoldgass['LOGSFR'] - xcoldgass['LOGMSTAR']

        dot = xcoldgass.index[(xcoldgass['H1_FLAG']!=99) & (xcoldgass['LOGMH2'].notnull())].tolist()
        H2up = xcoldgass.index[(xcoldgass['H1_FLAG']!=99) & (xcoldgass['LOGMH2_LIM'].notnull())].tolist()
        H1up = xcoldgass.index[(xcoldgass['H1_FLAG']==99) & (xcoldgass['LOGMH2'].notnull())].tolist()
        both = xcoldgass.index[(xcoldgass['H1_FLAG']==99) & (xcoldgass['LOGMH2_LIM'].notnull())].tolist()

        plt.scatter(x[dot], y[dot], c='Blue', alpha=0.75, marker='o', edgecolors='Blue')
        plt.scatter(x[both], yup[both], c='Grey', alpha=0.75, marker=r'$\downarrow$')
        plt.scatter(x[H1up], y[H1up], c='Blue', alpha=0.75, marker=r'$\downarrow$')
        plt.scatter(x[H2up], yup[H2up], c='Blue', alpha=0.75, marker='v', edgecolors='Blue')

    y = src['LOGMGAS'] - src['LOGMSTAR_MAGPHYS']
    x = src['LOGSFR_MAGPHYS'] - src['LOGMSTAR_MAGPHYS']

    dot = src.index[(src['H1_F']==1) & (src['H2_F']==1)].tolist()
    H1up = src.index[(src['H1_F']==3) & (src['H2_F']==1)].tolist()
    hollow_up = src.index[(src['H1_F']==3) & (src['H2_F']==2)].tolist()
    tup = src.index[(src['H1_F']==1) & (src['H2_F']==2)].tolist()
    tdown = src.index[(src['H1_F']==2) & (src['H2_F']==1)].tolist()

    plt.scatter(x[dot], y[dot], c='DarkOrange', alpha=0.75, marker='o', edgecolors='DarkOrange')
    plt.scatter(x[H1up], y[H1up], c='DarkOrange', alpha=0.75, marker=r'$\downarrow$')
    plt.scatter(x[hollow_up], y[hollow_up], c='White', alpha=0.75, marker='o', edgecolors='DarkOrange')

    plt.scatter(x[tup], y[tup], c='DarkOrange', alpha=0.75, marker='^', edgecolors='Red')
    plt.scatter(x[tdown], y[tdown], c='DarkOrange', alpha=0.75, marker='^', edgecolors='Blue')

    plt.ylabel('log $M_{gas}$/$M_{*}$')      
    plt.xlabel('log sSFR [$y^{-1}$]')

    plt.show()

def New_Figure3(src, xcg=CF.src6.copy(), showNon=False):
    xcg['LOGMGAS'] = np.log10(1.3*(np.power(10,xcg['LOGMH2'])+np.power(10,xcg['LOGMH1'])))
    xcg['LOGMGAS_UP'] = np.log10(1.3*(np.power(10,xcg['LOGMH2_LIM'])+np.power(10,xcg['LOGMH1'])))

    indXCG = xcg.index[(xcg['LOGSFR'] - xcg['LOGMSTAR'] >= -11) & ((xcg['LOGMGAS'] - xcg['LOGMSTAR'] >= -0.5) | (xcg['LOGMGAS_UP'] - xcg['LOGMSTAR'] >= -0.5))].tolist()
    indJGL = src.index[(src['LOGSFR_MAGPHYS'] - src['LOGMSTAR_MAGPHYS'] >= -11) & (src['LOGMGAS'] - src['LOGMSTAR_MAGPHYS'] >= -0.5)].tolist()

    notXCG = xcg.index[(xcg['LOGSFR'] - xcg['LOGMSTAR'] < -11) | ((xcg['LOGMGAS'] - xcg['LOGMSTAR'] < -0.5) | (xcg['LOGMGAS_UP'] - xcg['LOGMSTAR'] < -0.5))].tolist()
    notJGL = src.index[(src['LOGSFR_MAGPHYS'] - src['LOGMSTAR_MAGPHYS'] < -11) | (src['LOGMGAS'] - src['LOGMSTAR_MAGPHYS'] < -0.5)].tolist()

    #SFR vs Mstar
    if showNon == True:
        x = xcg['LOGMSTAR'][notXCG]
        y = xcg['LOGSFR'][notXCG]

        plt.scatter(x, y, c='Grey', alpha=0.75, marker='x')

        x = src['LOGMSTAR_MAGPHYS'][notJGL]
        y = src['LOGSFR_MAGPHYS'][notJGL]

        plt.scatter(x, y, c='Grey', alpha=0.75, marker='o')

    x = xcg['LOGMSTAR'][indXCG]
    y = xcg['LOGSFR'][indXCG]

    plt.scatter(x, y, c='Blue', alpha=0.75, marker='o', edgecolors= "Blue", label='xCOLD GASS')

    x = src['LOGMSTAR_MAGPHYS'][indJGL]
    y = src['LOGSFR_MAGPHYS'][indJGL]

    plt.scatter(x, y, c='DarkOrange', alpha=0.75, marker='o', edgecolors= "DarkOrange", label='JINGLE')

    plt.xlabel('log $M_{*}$ [$M_{\odot}$]')
    plt.ylabel('log $SFR$ [$M_{\odot} y^{-1}$]')

    plt.legend()
    plt.show()

    #SFR vs Mgas
    if showNon == True:
        x = xcg['LOGMGAS'][notXCG]
        xup = xcg['LOGMGAS_UP'][notXCG]
        y = xcg['LOGSFR'][notXCG]

        plt.scatter(x, y, c='Grey', alpha=0.75, marker='x')
        plt.scatter(xup, y, c='Grey', alpha=0.75, marker='x')

        x = src['LOGMGAS'][notJGL]
        y = src['LOGSFR_MAGPHYS'][notJGL]

        plt.scatter(x, y, c='Grey', alpha=0.75, marker='o')

    x = xcg['LOGMGAS'][indXCG]
    xup = xcg['LOGMGAS_UP'][indXCG]
    y = xcg['LOGSFR'][indXCG]

    dot = list(coll.Counter(xcg.index[(xcg['H1_FLAG']!=99) & (xcg['LOGMH2'].notnull())].tolist()) & coll.Counter(indXCG))
    H2up = list(coll.Counter(xcg.index[(xcg['H1_FLAG']!=99) & (xcg['LOGMH2_LIM'].notnull())].tolist()) & coll.Counter(indXCG))
    H1up = list(coll.Counter(xcg.index[(xcg['H1_FLAG']==99) & (xcg['LOGMH2'].notnull())].tolist()) & coll.Counter(indXCG))
    both = list(coll.Counter(xcg.index[(xcg['H1_FLAG']==99) & (xcg['LOGMH2_LIM'].notnull())].tolist()) & coll.Counter(indXCG))

    plt.scatter(x[dot], y[dot], c='Blue', alpha=0.75, marker='o', edgecolors= "Blue", label='xCOLD GASS')
    plt.scatter(xup[both], y[both], c='White', alpha=0.75, marker='o', edgecolors= "Blue")
    plt.scatter(x[H1up], y[H1up], c='Blue', alpha=0.75, marker='<', edgecolors='Red')
    plt.scatter(xup[H2up], y[H2up], c='Blue', alpha=0.75, marker='<', edgecolors='Black')

    x = src['LOGMGAS'][indJGL]
    y = src['LOGSFR_MAGPHYS'][indJGL]

    dot = list(coll.Counter(src.index[(src['H1_F']==1) & (src['H2_F']==1)].tolist()) & coll.Counter(indJGL))
    H1up = list(coll.Counter(src.index[(src['H1_F']==3) & (src['H2_F']==1)].tolist()) & coll.Counter(indJGL))
    hollow_up = list(coll.Counter(src.index[(src['H1_F']==3) & (src['H2_F']==2)].tolist()) & coll.Counter(indJGL))
    H1dust = list(coll.Counter(src.index[(src['H1_F']==2) & (src['H2_F']==1)].tolist()) & coll.Counter(indJGL))
    H2dust = list(coll.Counter(src.index[(src['H1_F']==1) & (src['H2_F']==2)].tolist()) & coll.Counter(indJGL))
    
    plt.scatter(x[dot], y[dot], c='DarkOrange', alpha=0.75, marker='o', edgecolors='DarkOrange', label='JINGLE')
    plt.scatter(x[H1up], y[H1up], c='DarkOrange', alpha=0.75, marker='<')
    plt.scatter(x[hollow_up], y[hollow_up], c='White', alpha=0.75, marker='o', edgecolors='DarkOrange')
    plt.scatter(x[H1dust], y[H1dust], c='DarkOrange', alpha=0.75, marker='>', edgecolors='Red')
    plt.scatter(x[H2dust], y[H2dust], c='DarkOrange', alpha=0.75, marker='>', edgecolors='Black')

    plt.xlabel('log $M_{gas}$ [$M_{\odot}$]')
    plt.ylabel('log $SFR$ [$M_{\odot} y^{-1}$]')

    plt.legend()
    plt.show()

    #Fgas vs Mstar
    if showNon == True:
        y = xcg['LOGMGAS'][notXCG] - xcg['LOGMSTAR'][notXCG]
        yup = xcg['LOGMGAS_UP'][notXCG] - xcg['LOGMSTAR'][notXCG]
        x = xcg['LOGMSTAR'][notXCG]

        plt.scatter(x, y, c='Grey', alpha=0.75, marker='x')
        plt.scatter(x, yup, c='Grey', alpha=0.75, marker='x')

        y = src['LOGMGAS'][notJGL] - src['LOGMSTAR_MAGPHYS'][notJGL]
        x = src['LOGMSTAR_MAGPHYS'][notJGL]

        plt.scatter(x, y, c='Grey', alpha=0.75, marker='o')

    y = xcg['LOGMGAS'][indXCG] - xcg['LOGMSTAR'][indXCG]
    yup = xcg['LOGMGAS_UP'][indXCG] - xcg['LOGMSTAR'][indXCG]
    x = xcg['LOGMSTAR'][indXCG]

    dot = list(coll.Counter(xcg.index[(xcg['H1_FLAG']!=99) & (xcg['LOGMH2'].notnull())].tolist()) & coll.Counter(indXCG))
    H2up = list(coll.Counter(xcg.index[(xcg['H1_FLAG']!=99) & (xcg['LOGMH2_LIM'].notnull())].tolist()) & coll.Counter(indXCG))
    H1up = list(coll.Counter(xcg.index[(xcg['H1_FLAG']==99) & (xcg['LOGMH2'].notnull())].tolist()) & coll.Counter(indXCG))
    both = list(coll.Counter(xcg.index[(xcg['H1_FLAG']==99) & (xcg['LOGMH2_LIM'].notnull())].tolist()) & coll.Counter(indXCG))

    plt.scatter(x[dot], y[dot], c='Blue', alpha=0.75, marker='o', edgecolors= "Blue", label='xCOLD GASS')
    plt.scatter(x[both], yup[both], c='White', alpha=0.75, marker='o', edgecolors= "Blue")
    plt.scatter(x[H1up], y[H1up], c='Blue', alpha=0.75, marker='v', edgecolors='Red')
    plt.scatter(x[H2up], yup[H2up], c='Blue', alpha=0.75, marker='v', edgecolors='Black')

    y = src['LOGMGAS'][indJGL] - src['LOGMSTAR_MAGPHYS'][indJGL]
    x = src['LOGMSTAR_MAGPHYS'][indJGL]

    dot = list(coll.Counter(src.index[(src['H1_F']==1) & (src['H2_F']==1)].tolist()) & coll.Counter(indJGL))
    H1up = list(coll.Counter(src.index[(src['H1_F']==3) & (src['H2_F']==1)].tolist()) & coll.Counter(indJGL))
    hollow_up = list(coll.Counter(src.index[(src['H1_F']==3) & (src['H2_F']==2)].tolist()) & coll.Counter(indJGL))
    H1dust = list(coll.Counter(src.index[(src['H1_F']==2) & (src['H2_F']==1)].tolist()) & coll.Counter(indJGL))
    H2dust = list(coll.Counter(src.index[(src['H1_F']==1) & (src['H2_F']==2)].tolist()) & coll.Counter(indJGL))
    
    plt.scatter(x[dot], y[dot], c='DarkOrange', alpha=0.75, marker='o', edgecolors='DarkOrange', label='JINGLE')
    plt.scatter(x[H1up], y[H1up], c='DarkOrange', alpha=0.75, marker='v')
    plt.scatter(x[hollow_up], y[hollow_up], c='White', alpha=0.75, marker='o', edgecolors='DarkOrange')
    plt.scatter(x[H1dust], y[H1dust], c='DarkOrange', alpha=0.75, marker='^', edgecolors='Red')
    plt.scatter(x[H2dust], y[H2dust], c='DarkOrange', alpha=0.75, marker='^', edgecolors='Black')

    plt.xlabel('log $M_{star}$ [$M_{\odot}$]')
    plt.ylabel('log $M_{gas}$/$M_{*}$')

    plt.legend()
    plt.show()

    #Fgas vs sSFR
    if showNon == True:
        y = xcg['LOGMGAS'][notXCG] - xcg['LOGMSTAR'][notXCG]
        yup = xcg['LOGMGAS_UP'][notXCG] - xcg['LOGMSTAR'][notXCG]
        x = xcg['LOGSFR'][notXCG] - xcg['LOGMSTAR'][notXCG]

        plt.scatter(x, y, c='Grey', alpha=0.75, marker='x')
        plt.scatter(x, yup, c='Grey', alpha=0.75, marker='x')

        y = src['LOGMGAS'][notJGL] - src['LOGMSTAR_MAGPHYS'][notJGL]
        x = src['LOGSFR_MAGPHYS'][notJGL] - src['LOGMSTAR_MAGPHYS'][notJGL]

        plt.scatter(x, y, c='Grey', alpha=0.75, marker='o')

    y = xcg['LOGMGAS'][indXCG] - xcg['LOGMSTAR'][indXCG]
    yup = xcg['LOGMGAS_UP'][indXCG] - xcg['LOGMSTAR'][indXCG]
    x = xcg['LOGSFR'][indXCG] - xcg['LOGMSTAR'][indXCG]

    dot = list(coll.Counter(xcg.index[(xcg['H1_FLAG']!=99) & (xcg['LOGMH2'].notnull())].tolist()) & coll.Counter(indXCG))
    H2up = list(coll.Counter(xcg.index[(xcg['H1_FLAG']!=99) & (xcg['LOGMH2_LIM'].notnull())].tolist()) & coll.Counter(indXCG))
    H1up = list(coll.Counter(xcg.index[(xcg['H1_FLAG']==99) & (xcg['LOGMH2'].notnull())].tolist()) & coll.Counter(indXCG))
    both = list(coll.Counter(xcg.index[(xcg['H1_FLAG']==99) & (xcg['LOGMH2_LIM'].notnull())].tolist()) & coll.Counter(indXCG))

    plt.scatter(x[dot], y[dot], c='Blue', alpha=0.75, marker='o', edgecolors= "Blue", label='xCOLD GASS')
    plt.scatter(x[both], yup[both], c='White', alpha=0.75, marker='o', edgecolors= "Blue")
    plt.scatter(x[H1up], y[H1up], c='Blue', alpha=0.75, marker='v', edgecolors='Red')
    plt.scatter(x[H2up], yup[H2up], c='Blue', alpha=0.75, marker='v', edgecolors='Black')

    y = src['LOGMGAS'][indJGL] - src['LOGMSTAR_MAGPHYS'][indJGL]
    x = src['LOGSFR_MAGPHYS'][indJGL] - src['LOGMSTAR_MAGPHYS'][indJGL]

    dot = list(coll.Counter(src.index[(src['H1_F']==1) & (src['H2_F']==1)].tolist()) & coll.Counter(indJGL))
    H1up = list(coll.Counter(src.index[(src['H1_F']==3) & (src['H2_F']==1)].tolist()) & coll.Counter(indJGL))
    hollow_up = list(coll.Counter(src.index[(src['H1_F']==3) & (src['H2_F']==2)].tolist()) & coll.Counter(indJGL))
    H1dust = list(coll.Counter(src.index[(src['H1_F']==2) & (src['H2_F']==1)].tolist()) & coll.Counter(indJGL))
    H2dust = list(coll.Counter(src.index[(src['H1_F']==1) & (src['H2_F']==2)].tolist()) & coll.Counter(indJGL))
    
    plt.scatter(x[dot], y[dot], c='DarkOrange', alpha=0.75, marker='o', edgecolors='DarkOrange', label='JINGLE')
    plt.scatter(x[H1up], y[H1up], c='DarkOrange', alpha=0.75, marker='v')
    plt.scatter(x[hollow_up], y[hollow_up], c='White', alpha=0.75, marker='o', edgecolors='DarkOrange')
    plt.scatter(x[H1dust], y[H1dust], c='DarkOrange', alpha=0.75, marker='^', edgecolors='Red')
    plt.scatter(x[H2dust], y[H2dust], c='DarkOrange', alpha=0.75, marker='^', edgecolors='Black')

    plt.xlabel('log sSFR [$y^{-1}$]')
    plt.ylabel('log $M_{gas}$/$M_{*}$')

    plt.legend()
    plt.show()

