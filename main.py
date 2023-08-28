"""
Main function. Run this code to access the data from the local JINGLE csv files.
This code is meant to give a "quick" look at the data within the JINGLE files, showing a brief overlook at the main information.
"""

#import libraries here
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy as sci
import math as mt
import warnings

from pandas.errors import SettingWithCopyWarning

#import other scripts here
import gather_data as GD
import base_plotter as BPLT
import core_analysis as CA
import formatting_functions as FF
import ratio_plots as RPLT
import constants as C

#set-up parameters and globals here
src1 = GD.get_data(GD.JINGLE_MASTER)
src2 = GD.get_data(GD.SED_FITTINGS)
src3 = GD.get_data(GD.JINGLE_TEMPEL)

#ignore pandas warning
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)

#fixing upper limit duplicates
src1.loc[src1['LOGMH2'] != 0, 'LOGMH2_PRED'] = 0

#generate average errors for DeLooze Dust Masses
src1['LOGMDUST_DELOOZE_ERR'] = [(abs(g)+abs(h))/2 for g,h in zip(src1['-dex'], src1['+dex'])]
src1['LOGMH2_ALL'] = src1['LOGMH2'] + src1['LOGMH2_PRED']

src1 = CA.generate_std_error(src1, src1['LOGMH1'])
src1 = CA.generate_std_error(src1, src1['LOGMH2_ALL'])

src1.loc[src1['LOGMH1'] <= 0.5, 'LOGMH1'] = np.nan
src1.loc[src1['LOGMH2_ALL'] <= 0.5, 'LOGMH2_ALL'] = np.nan

#generate the total gas masses and metals
src1 = CA.find_Mmetals(src1)

#main function for analyzing data
def main():
    print("Hello and welcome to the plotting software!")
    print("-------------------------------------------")
    return

def jingle_galaxy_base_parameters(Dust = True, Gas = True, SFR = True):
    if Dust == True:
        df1 = src1[['JINGLEID','IDNUM','LOGMSTAR_MAGPHYS','LOGMSTAR_MAGPHYS_ERR','LOGMDUST_MAGPHYS','LOGMDUST_DELOOZE','-dex','+dex']].copy()
        df2 = src2[['IDNUM','logMc_SMBB','e_logMc_SMBB','logMc_BMBB','e_logMc_BMBB','logMc_TMBB','e_logMc_TMBB']].copy()
        df = pd.merge(df1, df2, on='IDNUM')

        BPLT.Mstar_vs_Mdust(src=df)
    
    if Gas == True:
        df = src1[['JINGLEID','IDNUM','LOGMSTAR_MAGPHYS','LOGMSTAR_MAGPHYS_ERR', 'LOGMH1','LOGMH1_FRAC','LOGMH2','LOGMH2_PRED']].copy()

        BPLT.Mstar_vs_Gas(src=df)
    
    if SFR == True:
        df = src1[['JINGLEID','IDNUM','LOGMSTAR_MAGPHYS','LOGMSTAR_MAGPHYS_ERR', 'LOGSFR_MAGPHYS', 'LOGSFR_MAGPHYS_ERR']].copy()

        BPLT.Mstar_vs_SFR(src=df)

def gas_content_comparisons():
    df = src1[['JINGLEID','IDNUM','Z_PP04_N2','Z_PP04_O3N2','Z_MZR','LOGMSTAR_MAGPHYS','LOGMSTAR_MAGPHYS_ERR','LOGMDUST_MAGPHYS','LOGMH1','LOGMH2','LOGMH2_PRED']].copy()

    BPLT.Gas_Content(src=df)

    BPLT.H1_vs_H2(src=df)

def specific_SFR(Mstar = True, Mdust = True, Mh1 = True, Mh2 = True, Mgas = True):
    df = src1[['JINGLEID','IDNUM','LOGSFR_MAGPHYS','LOGSFR_MAGPHYS_ERR','LOGMSTAR_MAGPHYS','LOGMSTAR_MAGPHYS_ERR','LOGMDUST_MAGPHYS','LOGMH1','LOGMH2','LOGMH2_PRED']].copy()
    
    if Mstar == True:
        BPLT.SSFR(src=df, s=df['LOGMSTAR_MAGPHYS'])
    if Mdust == True:
        BPLT.SSFR(src=df, s=df['LOGMDUST_MAGPHYS'])
    if Mh1 == True:
        BPLT.SSFR(src=df, s=df['LOGMH1'])
    if Mh2 == True:
        BPLT.SSFR(src=df, s=df['LOGMH2'])
    if Mgas == True:
        gas = CA.Calc_Gas_Content_Total(src=df)
        BPLT.SSFR(src=df, s=gas)

def grouped_by(Env = False, Den = False):
    df1 = src1.copy()
    df3 = src3[['IDNUM','IDGAL','NGAL']].copy()

    df = pd.merge(df1, df3, on='IDNUM')

    df = CA.Group_By_Env(df)
    df = CA.Group_By_Dens(df)

    if Env == True and Den == False:
        grouper = df.groupby('GALACTIC_ENV')
    elif Env == False and Den == True:
        grouper = df.groupby('GALACTIC_DENS')
    else:
        grouper = df.groupby('GALACTIC_DENS')

    #FF.print_full(df)

    #for key, group in grouper:
        #fig,ax = plt.subplots()
        #BPLT.linmix_plots(key, group, group['LOGMSTAR_MAGPHYS'], group['LOGSFR_MAGPHYS'], group['LOGMSTAR_MAGPHYS_ERR'], group['LOGSFR_MAGPHYS_ERR'], ax=ax)
        #BPLT.linmix_plots(key, group, group['LOGMSTAR_MAGPHYS'], group['LOGMDUST_DELOOZE'], group['LOGMSTAR_MAGPHYS_ERR'], group['LOGMDUST_DELOOZE_ERR'], ax=ax)
        #BPLT.linmix_plots(key, group, group['LOGMSTAR_MAGPHYS'], group['LOGMH1'], group['LOGMSTAR_MAGPHYS_ERR'], group['LOGMH1_ERR'], ax=ax)
        #BPLT.linmix_plots(key, group, group['LOGMSTAR_MAGPHYS'], group['LOGMH2_ALL'], group['LOGMSTAR_MAGPHYS_ERR'], group['LOGMH2_ALL_ERR'], ax=ax)
        #plt.show()

    counter = 0
    fig,ax = plt.subplots()

    for key, group in grouper:
        BPLT.linmix_plots_multi(key, group, group['LOGMSTAR_MAGPHYS'], group['LOGMDUST_DELOOZE'], group['LOGMSTAR_MAGPHYS_ERR'], group['LOGMDUST_DELOOZE_ERR'], ax=ax, color=C.COLORS[counter], LMcolor=C.COLORS[counter])
        counter+=1
    
    plt.show()
    
    
def galaxy_ratios():
    df1 = src1.copy()

    #RPLT.x_Mstar(df1, df1['LOGMSTAR_MAGPHYS'], df1['LOGMSTAR_MAGPHYS_ERR'], '$M_{star}$ [Log($M_{\odot}$)]')
    RPLT.x_Mstar(df1, df1['LOGSFR_MAGPHYS'], df1['LOGSFR_MAGPHYS_ERR'], '$SFR$ [Log($M_{\odot}$/yr)]')

#call on the main function when the script is executed
if __name__ == '__main__':
    #main()
    #jingle_galaxy_base_parameters()
    #gas_content_comparisons()
    #specific_SFR()
    grouped_by(Env=False)
    #galaxy_ratios()