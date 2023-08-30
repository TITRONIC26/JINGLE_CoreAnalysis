"""
This script is intended to be populated with functions that will be called directly from the main directory.
This script is not intended to perform explicit plotting, calculations, or other types of analysis.
"""

#import libraries here
import pandas as pd
import matplotlib.pyplot as plt

#import other scripts here
import gather_data as GD
import base_plotter as BPLT
import core_analysis as CA
import ratio_plots as RPLT
import constants as C
import data_manipulation as DM

#set-up parameters and globals here
src1 = GD.get_data(GD.JINGLE_MASTER)
src2 = GD.get_data(GD.SED_FITTINGS)
src3 = GD.get_data(GD.JINGLE_TEMPEL)
src4 = GD.get_data(GD.XCOLDGASS)
src5 = GD.get_data(GD.VERTICO)

#fix src1 data
src1 = DM.JINGLE_main(src1)
src4 = DM.XCOLDGASS_main(src4)
src5 = DM.VERTICO_main(src5)

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
        BPLT.linmix_plots_multi(key, group, group['LOGMSTAR_MAGPHYS'], group['LOGMMETAL'], group['LOGMSTAR_MAGPHYS_ERR'], group['LOGMMETAL_ERR'], ax=ax, color=C.COLORS[counter], LMcolor=C.COLORS[counter])
        counter+=1
    
    plt.show()
    
def galaxy_ratios():
    df1 = src1.copy()

    RPLT.x_Mstar(df1, df1['LOGMSTAR_MAGPHYS'], df1['LOGMSTAR_MAGPHYS_ERR'], '$M_{star}$ [Log($M_{\odot}$)]')
    RPLT.x_Mstar(df1, df1['LOGSFR_MAGPHYS'], df1['LOGSFR_MAGPHYS_ERR'], '$SFR$ [Log($M_{\odot}$/yr)]')

def compare():
    df1 = src1.copy()
    df2 = src4.copy()

    fig,ax = plt.subplots()

    #compare Mstar vs SFR
    BPLT.linmix_plots_multi('xCOLDGASS', df2, df2['LOGMSTAR'], df2['LOGMH2_ALL'], df2['LOGMSTAR_ERR'], df2['LOGMH2_ALL_ERR'], ax=ax, color=C.COLORS[1], LMcolor=C.COLORS[1])
    BPLT.linmix_plots_multi('JINGLE', df1, df1['LOGMSTAR_MAGPHYS'], df1['LOGMH2_ALL'], df1['LOGMSTAR_MAGPHYS_ERR'], df1['LOGMH2_ALL_ERR'], ax=ax, color=C.COLORS[0], LMcolor=C.COLORS[0])

    plt.show()