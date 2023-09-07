"""
This script is intended to be populated with functions that will be called directly from the main directory.
This script is not intended to perform explicit plotting, calculations, or other types of analysis.
"""

#import libraries here
import pandas as pd
import matplotlib.pyplot as plt

#plotting scripts
import base_plotter as BPLT
import ratio_plots as RPLT
import plotter as PLT

#essential scripts
import gather_data as GD
import constants as C

#data analysis and manipulation
import core_analysis as CA
import data_manipulation as DM
import linmixer as LM

#formatting
import formatting_functions as FF

#set-up parameters and globals here
src1 = GD.get_data(GD.JINGLE_GALAXY)
src2 = GD.get_data(GD.JINGLE_FLUX)
src3 = GD.get_data(GD.JINGLE_TEMPEL)
src4 = GD.get_data(GD.VERTICO_DP)

#fix src data
src1 = DM.JINGLE_main(src1)
src4 = DM.VERTICO_main(src4)

def main():
    #FF.print_full(src1)
    total_gas()

def total_gas():
    df1 = src1[['LOGMSTAR_MAGPHYS','LOGMSTAR_MAGPHYS_ERR','LOGMH2_RYAN','LOGMH1_MATT','LOGMH1_MATT_ERR','H1_FLAG','LOGMH2_RYAN_ERR','LOGMGAS','LOGMGAS_ERR','MGAS_FLAG','JINGLEID']].copy()
    df2 = src3[['JINGLEID','NGAL','LOGMHALO_TEMPEL']].copy()
    df3 = src2[['JINGLEID','IDNUM']].copy()

    src = pd.merge(df1, df2, on='JINGLEID')
    src = pd.merge(src, df3, on='JINGLEID')

    src = CA.Group_By_Dens(src)
    grouper = src.groupby('GALACTIC_DENS')

    for key, group in grouper:
        print(key)
        PLT.gas_plots(group, vertico=True)


#call on the main function when the script is executed
if __name__ == '__main__':
    main()

"""
Old functions, for archival purposes...
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
    df3 = src5.copy()

    fig,ax = plt.subplots()

    #compare Mstar vs SFR
    BPLT.linmix_plots_multi('VERTICO', df3, df3['LOGMSTAR_DP'], df3['LOGSFR_DP'], df3['LOGMSTAR_DP_ERR'], df3['LOGSFR_DP_ERR'], ax=ax, color=C.COLORS[2], LMcolor=C.COLORS[2])
    BPLT.linmix_plots_multi('xCOLDGASS', df2, df2['LOGMSTAR'], df2['LOGSFR_BEST'], df2['LOGMSTAR_ERR'], df2['LOGSFR_BEST_ERR'], ax=ax, color=C.COLORS[1], LMcolor=C.COLORS[1])
    BPLT.linmix_plots_multi('JINGLE', df1, df1['LOGMSTAR_MAGPHYS'], df1['LOGSFR_MAGPHYS'], df1['LOGMSTAR_MAGPHYS_ERR'], df1['LOGSFR_MAGPHYS_ERR'], ax=ax, color=C.COLORS[0], LMcolor=C.COLORS[0])

    #compare Mstar vs Mdust
    BPLT.linmix_plots_multi('VERTICO', df3, df3['LOGMSTAR_DP'], df3['LOGMDUST_DP'], df3['LOGMSTAR_DP_ERR'], df3['LOGMDUST_DP_ERR'], ax=ax, color=C.COLORS[2], LMcolor=C.COLORS[2])
    BPLT.linmix_plots_multi('JINGLE', df1, df1['LOGMSTAR_MAGPHYS'], df1['LOGMDUST_DELOOZE'], df1['LOGMSTAR_MAGPHYS_ERR'], df1['LOGMDUST_DELOOZE_ERR'], ax=ax, color=C.COLORS[0], LMcolor=C.COLORS[0])

    #compare Mstar vs MH1
    BPLT.linmix_plots_multi('VERTICO', df3, df3['LOGMSTAR_DP'], df3['LOGMH1_DP'], df3['LOGMSTAR_DP_ERR'], df3['LOGMH1_DP_ERR'], ax=ax, color=C.COLORS[2], LMcolor=C.COLORS[2])
    BPLT.linmix_plots_multi('JINGLE', df1, df1['LOGMSTAR_MAGPHYS'], df1['LOGMH1'], df1['LOGMSTAR_MAGPHYS_ERR'], df1['LOGMH1_ERR'], ax=ax, color=C.COLORS[0], LMcolor=C.COLORS[0])

    #compare Mstar vs MH2
    BPLT.linmix_plots_multi('VERTICO', df3, df3['LOGMSTAR_DP'], df3['LOGMH2_BROWN'], df3['LOGMSTAR_DP_ERR'], df3['LOGMH2_BROWN_ERR'], ax=ax, color=C.COLORS[2], LMcolor=C.COLORS[2])
    BPLT.linmix_plots_multi('xCOLDGASS', df2, df2['LOGMSTAR'], df2['LOGMH2_ALL'], df2['LOGMSTAR_ERR'], df2['LOGMH2_ALL_ERR'], ax=ax, color=C.COLORS[1], LMcolor=C.COLORS[1])
    BPLT.linmix_plots_multi('JINGLE', df1, df1['LOGMSTAR_MAGPHYS'], df1['LOGMH2_ALL'], df1['LOGMSTAR_MAGPHYS_ERR'], df1['LOGMH2_ALL_ERR'], ax=ax, color=C.COLORS[0], LMcolor=C.COLORS[0])

    #compare Mstar vs Mmetal
    BPLT.linmix_plots_multi('VERTICO', df3, df3['LOGMSTAR_DP'], df3['LOGMMETAL_DP'], df3['LOGMSTAR_DP_ERR'], df3['LOGMMETAL_DP_ERR'], ax=ax, color=C.COLORS[2], LMcolor=C.COLORS[2])
    BPLT.linmix_plots_multi('JINGLE', df1, df1['LOGMSTAR_MAGPHYS'], df1['LOGMMETAL'], df1['LOGMSTAR_MAGPHYS_ERR'], df1['LOGMMETAL_ERR'], ax=ax, color=C.COLORS[0], LMcolor=C.COLORS[0])

    plt.show()
"""