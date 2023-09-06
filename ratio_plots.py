"""
OUTDATED
Plots generated using the ratios of different parameters within the JINGLE data set
"""

#import libraries here
import numpy as np
import matplotlib.pyplot as plt

#import other scripts and constants here
import constants as C
import linmixer as LM
import math_funcs as MF

#declare global parameters here

#define functions here
def remove_nans(input):
    df = input.to_frame()
    index_nan = np.where(df.notnull())[0]

    return index_nan

def four_plot(x, x_err, y, y_err, names, x_label):
    fig,axs = plt.subplots(2,2, sharex=True)

    fig.subplots_adjust(hspace=0.35, wspace=0.1)

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

def x_Mstar(src, x, x_err, xtitle):
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

    gas = src['LOGMGAS']
    gas_err = src['LOGMGAS_ERR']

    #Mstar plots
    values = [dust, h_alpha, h_mol, metal]
    values_err = [dust_err, h_alpha_err, h_mol_err, metal_err]
    names = ['$M_{dust}$ [Log($M_{\odot}$)]', '$M_{H1}$ [Log($M_{\odot}$)]', '$M_{H2}$ [Log($M_{\odot}$)]', '$M_{metal}$ [Log($M_{\odot}$)]']
    four_plot(x, x_err, values, values_err, names, xtitle)

    #Mstar plots normalized by Mstar
    values = [dust-star, h_alpha-star, h_mol-star, metal-star]
    values_err = [MF.addErrors(dust_err, star_err), MF.addErrors(h_alpha_err, star_err), MF.addErrors(h_mol_err, star_err), MF.addErrors(metal_err, star_err)]
    names = ['$M_{dust}$/$M_{star}$ [Log]', '$M_{H1}$/$M_{star}$ [Log]', '$M_{H2}$/$M_{star}$ [Log]', '$M_{metal}$/$M_{star}$ [Log]']
    four_plot(x, x_err, values, values_err, names, xtitle)

    #Mstar plots normalized by Mdust
    values = [gas-dust, h_alpha-dust, h_mol-dust, metal-dust]
    values_err = [MF.addErrors(gas_err, dust_err), MF.addErrors(h_alpha_err, dust_err), MF.addErrors(h_mol_err, dust_err), MF.addErrors(metal_err, dust_err)]
    names = ['$M_{gas}$/$M_{dust}$ [Log]', '$M_{H1}$/$M_{dust}$ [Log]', '$M_{H2}$/$M_{dust}$ [Log]', '$M_{metal}$/$M_{dust}$ [Log]']
    four_plot(x, x_err, values, values_err, names, xtitle)

    #Mstar plots normalized by Mgas
    values = [dust-gas, h_alpha-gas, h_mol-gas, metal-gas]
    values_err = [MF.addErrors(dust_err, gas_err), MF.addErrors(h_alpha_err, gas_err), MF.addErrors(h_mol_err, gas_err), MF.addErrors(metal_err, gas_err)]
    names = ['$M_{dust}$/$M_{gas}$ [Log]', '$M_{H1}$/$M_{gas}$ [Log]', '$M_{H2}$/$M_{gas}$ [Log]', '$M_{metal}$/$M_{gas}$ [Log]']
    four_plot(x, x_err, values, values_err, names, xtitle)

    #Mstar plots normalized by Mmetal
    values = [dust-metal, h_alpha-metal, h_mol-metal, gas-metal]
    values_err = [MF.addErrors(dust_err, metal_err), MF.addErrors(h_alpha_err, metal_err), MF.addErrors(h_mol_err, metal_err), MF.addErrors(gas_err, metal_err)]
    names = ['$M_{dust}$/$M_{metal}$ [Log]', '$M_{H1}$/$M_{metal}$ [Log]', '$M_{H2}$/$M_{metal}$ [Log]', '$M_{gas}$/$M_{metal}$ [Log]']
    four_plot(x, x_err, values, values_err, names, xtitle)