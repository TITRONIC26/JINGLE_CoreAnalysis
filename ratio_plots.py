"""
Plots generated using the ratios of different parameters within the JINGLE data set
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
import core_analysis as CA
import linmixer as LM
import formatting_functions as FF
import math_funcs as MF

#declare global parameters here

#define functions here
def remove_nans(input):
    df = input.to_frame()
    index_nan = np.where(df.notnull())[0]

    return index_nan

def four_plot(x, x_err, y, y_err, names, x_label):
    fig,axs = plt.subplots(2,2, sharex=True)

    fig.subplots_adjust(hspace=0.35, wspace=0.05)

    axes = [axs[0,0], axs[0,1], axs[1,0], axs[1,1]]
    counter = 0

    for ax in axes:
        index = remove_nans(y[counter])
        
        xs = x[index]
        xs_err = x_err[index]
        ys = y[counter][index]
        ys_err = y_err[counter][index]

        LM.pearson(xs, ys, ax)
        LM.curvefitting(xs, ys, ax)
        LM.linmixing(xs, ys, xs_err, ys_err, ax)
        ax.errorbar(xs, ys, ys_err, xs_err, fmt='o', color='Purple', ecolor='Black', markersize=C.SIZE/3, alpha=C.ALPHA)
    
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, borderaxespad=0., mode='expand')
        
        ax.set_ylabel(names[counter])

        counter+=1

    for ax in [axs[1,0], axs[1,1]]:
        ax.set_xlabel(x_label)
    
    plt.show()

def x_Mstar(src):
    star = src['LOGMSTAR_MAGPHYS']
    star_err = src['LOGMSTAR_MAGPHYS_ERR']

    dust = src['LOGMDUST_DELOOZE']
    dust_err = src['LOGMDUST_DELOOZE_ERR']

    h_alpha = src['LOGMH1']
    h_alpha_err = src['LOGMH1_ERR']

    h_mol = src['LOGMH2_ALL']
    h_mol_err = src['LOGMH2_ALL_ERR']

    metal = src['LOGMMETAL']
    metal_err = src['LOGMMETAL_ERR']

    #Mstar plots
    values = [dust, h_alpha, h_mol, metal]
    values_err = [dust_err, h_alpha_err, h_mol_err, metal_err]
    names = ['$M_{dust}$ [Log($M_{\odot}$)]', '$M_{H1}$ [Log($M_{\odot}$)]', '$M_{H2}$ [Log($M_{\odot}$)]', '$M_{metal}$ [Log($M_{\odot}$)]']
    four_plot(star, star_err, values, values_err, names, '$M_{star}$ [Log($M_{\odot}$)]')