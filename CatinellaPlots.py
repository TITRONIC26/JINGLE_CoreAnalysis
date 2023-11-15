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