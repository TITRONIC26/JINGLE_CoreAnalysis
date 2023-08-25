"""
Main function. Run this code to access the data from the local JINGLE csv files.
This code is meant to give a "quick" look at the data within the JINGLE files, showing brief a brief overlook at the main information.
"""

#import libraries here
import pandas as pd
import numpy as np
import matplotlib as mpl
import scipy as sci
import math as mt
import warnings

from pandas.errors import SettingWithCopyWarning

#import other scripts here
import gather_data as GD
import base_plotter as BPLT
import core_analysis as CA
import formatting_functions as FF

#set-up parameters and globals here
src1 = GD.get_data(GD.JINGLE_MASTER)
src2 = GD.get_data(GD.SED_FITTINGS)
src3 = GD.get_data(GD.JINGLE_TEMPEL)

#ignore pandas warning
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)

#fixing upper limit duplicates
src1.loc[src1['LOGMH2'] != 0, 'LOGMH2_PRED'] = 0

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

def linmix_datasets():
    df = src1[['JINGLEID','IDNUM','LOGSFR_MAGPHYS','LOGSFR_MAGPHYS_ERR','LOGMSTAR_MAGPHYS','LOGMSTAR_MAGPHYS_ERR']].copy()

    BPLT.linmix_plots(df, df['LOGMSTAR_MAGPHYS'], df['LOGSFR_MAGPHYS'], df['LOGMSTAR_MAGPHYS_ERR'], df['LOGSFR_MAGPHYS_ERR'])

def grouped_by_density():
    df1 = src1[['JINGLEID','IDNUM','LOGMHALO']].copy()
    df3 = src3[['IDNUM','IDGAL','NGAL']].copy()

    df = pd.merge(df1, df3, on='IDNUM')

    df = CA.Group_By_Density(df)

    FF.print_full(df)



#call on the main function when the script is executed
if __name__ == '__main__':
    #main()
    #jingle_galaxy_base_parameters()
    #gas_content_comparisons()
    #specific_SFR()
    #linmix_datasets()
    grouped_by_density()
