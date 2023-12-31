"""
OUTDATED
Script reserved exclusively for plotting data from other sections.
Includes the plots from every section, so be mindful of the length of this script.
"""

#import libraries here
import pandas as pd
import matplotlib.pyplot as plt
#import other scripts and constants here
import constants as C
import core_analysis as CA
import linmixer as LM

#global variables here

#begin functions for plotting here
def Mstar_vs_Mdust(src):
    fig,ax = plt.subplots()

    #Mstar vs DeLooze
    ax.scatter(src['LOGMSTAR_MAGPHYS'], src['LOGMDUST_DELOOZE'], marker='*', s=C.SIZE, alpha=C.ALPHA, color='Purple', label='DeLooze, 2020')
    #Mstar vs MAGPHYS
    ax.scatter(src['LOGMSTAR_MAGPHYS'], src['LOGMDUST_MAGPHYS'], marker='o', s=C.SIZE, alpha=C.ALPHA, color='Blue', label='Ryan, 2018')
    #Mstar vs SMBB
    ax.scatter(src['LOGMSTAR_MAGPHYS'], src['logMc_SMBB'], marker='s', s=C.SIZE, alpha=C.ALPHA, color='Green', label='SMBB')
    #Mstar vs BMBB
    ax.scatter(src['LOGMSTAR_MAGPHYS'], src['logMc_BMBB'], marker='p', s=C.SIZE, alpha=C.ALPHA, color='Orange', label='BMBB')
    #Mstar vs TMBB
    ax.scatter(src['LOGMSTAR_MAGPHYS'], src['logMc_TMBB'], marker='h', s=C.SIZE, alpha=C.ALPHA, color='Red', label='TMBB')

    #plot labels
    ax.set_title('Evaluated Dust Masses')
    ax.set_xlabel('$M_{star}$ [Log($M_{\odot}$)]')
    ax.set_ylabel('$M_{dust}$ [Log($M_{\odot}$)]')

    #formatting plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    plt.show()

    return

def Mstar_vs_Gas(src):
    fig,ax = plt.subplots()

    #Mstar vs MH1
    ax.scatter(src['LOGMSTAR_MAGPHYS'], src['LOGMH1'], marker='*', s=C.SIZE, alpha=C.ALPHA, color='Purple', label='Ryan, 2020')
    #plot labels
    ax.set_title(r'$H_{\alpha}$ Mass')
    ax.set_xlabel('$M_{star}$ [Log($M_{\odot}$)]')
    ax.set_ylabel(r'$M_{H_{\alpha}}$ [Log($M_{\odot}$)]')
    #formatting plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_ylim(8.5, 10.5)
    plt.show()

    fig,ax = plt.subplots()

    #Mstar vs MH2
    ax.scatter(src['LOGMSTAR_MAGPHYS'], src['LOGMH2'], marker='*', s=C.SIZE*2, alpha=C.ALPHA, color='Purple', label='Ryan, 2020')
    ax.scatter(src['LOGMSTAR_MAGPHYS'], src['LOGMH2_PRED'], marker='v', s=C.SIZE, alpha=C.ALPHA, color='Grey', label='Upper Limits')
    #plot labels
    ax.set_title('$H_{2}$ Mass')
    ax.set_xlabel('$M_{star}$ [Log($M_{\odot}$)]')
    ax.set_ylabel('$M_{H2}}$ [Log($M_{\odot}$)]')
    #formatting plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_ylim(7.75, 10.75)
    plt.show()

def Mstar_vs_SFR(src):
    fig,ax = plt.subplots()

    #Mstar vs SFR
    ax.scatter(src['LOGMSTAR_MAGPHYS'], src['LOGSFR_MAGPHYS'], marker='*', s=C.SIZE*2, alpha=C.ALPHA, color='Purple', label='Ryan, 2020')
    #plot labels
    ax.set_title('Star Formation Rate')
    ax.set_xlabel('$M_{star}$ [Log($M_{\odot}$)]')
    ax.set_ylabel('SFR [Log($M_{\odot}$/yr)]')
    #formatting plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.show()

def Gas_Content(src):
    dust = src['LOGMDUST_MAGPHYS'].copy()
    metalicity = src[['Z_PP04_N2','Z_PP04_O3N2','Z_MZR']].copy()
    gas = src[['LOGMH1','LOGMH2','LOGMH2_PRED']].copy()

    refScalingGas = [CA.Calc_Gas_Content_From_Dust(src, dust, x) for x in metalicity]
    totalGas = CA.Calc_Gas_Content_Total(gas)

    fig,ax = plt.subplots()

    #Mstar vs totalGas
    ax.scatter(src['LOGMSTAR_MAGPHYS'], totalGas[0], marker='*', s=C.SIZE*2, alpha=C.ALPHA, color='Purple', label='Total Gas Content (Ryan, 2020)')
    ax.scatter(src['LOGMSTAR_MAGPHYS'], totalGas[1], marker='v', s=C.SIZE*2, alpha=C.ALPHA, color='Purple', label='Upper Limits (Ryan, 2020)')
    #Mstar vs refScalingGas
    ax.scatter(src['LOGMSTAR_MAGPHYS'], refScalingGas[0], marker='o', s=C.SIZE, alpha=C.ALPHA, color='Blue', label='Z_PP04_O3N2 (Remy-Ruyer, 2014)')
    #ax.scatter(src['LOGMSTAR_MAGPHYS'], refScalingGas[1], marker='s', s=C.SIZE, alpha=C.ALPHA, color='Green', label='Z_PP04_N2 (Remy-Ruyer, 2014)')
    #ax.scatter(src['LOGMSTAR_MAGPHYS'], refScalingGas[2], marker='p', s=C.SIZE, alpha=C.ALPHA, color='Orange', label='Z_MZR (Remy-Ruyer, 2014)')
    #plot labels
    ax.set_title('Total Gas Content Comparissons')
    ax.set_xlabel('$M_{star}$ [Log($M_{\odot}$)]')
    ax.set_ylabel('$M_{gas}$ [Log($M_{\odot}$)]')
    #formatting plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_ylim(7.75, 10.75)
    plt.show()

def H1_vs_H2(src):
    fig,ax = plt.subplots()

    #MH1 vs MH2
    ax.scatter(src['LOGMH1'], src['LOGMH2'], marker='*', s=C.SIZE*2, alpha=C.ALPHA, color='Purple', label='Ryan, 2020')
    ax.scatter(src['LOGMH1'], src['LOGMH2_PRED'], marker='v', s=C.SIZE*2, alpha=C.ALPHA, color='Purple', label='Upper Limits')
    #plot labels
    ax.set_title(r'$H_{\alpha}$ vs $H_{2}$')
    ax.set_xlabel(r'$M_{H_{\alpha}}$ [Log($M_{\odot}$)]')
    ax.set_ylabel('$M_{H_{2}}$ [Log($M_{\odot}$)]')
    #formatting plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    #set limits
    ax.set_xlim(8.5,10.5)
    ax.set_ylim(8.25, 10.25)
    plt.show()
    
def SSFR(src, s):
    if len(s) == 2:
        s1 = s[1]
        s = s[0]
        column_name = str(['LOGMGAS'])
    else:
        dfs=s.to_frame()
        column_name = str(dfs.columns.values)
    
    fig,ax = plt.subplots()

    y = src['LOGSFR_MAGPHYS']-s
    x_vals = src['LOGMSTAR_MAGPHYS'][y < -2]
    y_vals = y[y < -2]

    #SFR/s vs Mstar
    ax.scatter(x_vals, y_vals, marker='*', s=C.SIZE*2, alpha=C.ALPHA, color='Purple', label=column_name)

    if column_name == str(['LOGMH2']):
        y1 = src['LOGSFR_MAGPHYS']-src['LOGMH2_PRED']
        x1_vals = src['LOGMSTAR_MAGPHYS'][y1 < -2]
        y1_vals = y1[y1 < -2]
        
        ax.scatter(x1_vals, y1_vals, marker='v', s=C.SIZE*2, alpha=C.ALPHA, color='Purple', label='Upper Limits')
    
    if column_name == str(['LOGMGAS']):
        y1 = src['LOGSFR_MAGPHYS']-s1
        x1_vals = src['LOGMSTAR_MAGPHYS'][y1 < -2]
        y1_vals = y1[y1 < -2]

        ax.scatter(x1_vals, y1_vals, marker='v', s=C.SIZE*2, alpha=C.ALPHA, color='Purple', label='Upper Limits')   
    
    #plot labels
    ax.set_xlabel('$M_{star}$ [Log($M_{\odot}$)]')
    ax.set_ylabel('${SFR}/'+column_name+'$ [Log(1/yr)]')
    #formatting plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    #set limits
    #ax.set_xlim(8.5,10.5)
    #ax.set_ylim(8.25, 10.25)
    plt.show()
    
def clean_set(x, y, x_err, y_err):
    dfx = x.to_frame()
    dfy = y.to_frame()
    dfx_err = x_err.to_frame()
    dfy_err = y_err.to_frame()

    name1 = str(dfx.columns.values)
    name2 = str(dfy.columns.values)
    name3 = str(dfx_err.columns.values)
    name4 = str(dfy_err.columns.values)

    df = pd.DataFrame({name1: x, name2: y, name3: x_err, name4: y_err})
    df = df.dropna()

    return (df[name1], df[name2], df[name3], df[name4])

def linmix_plots(key,src,x,y,x_err,y_err, ax=plt):
    dfx=x.to_frame()
    x_name = str(dfx.columns.values)
    dfy=y.to_frame()
    y_name = str(dfy.columns.values)

    (x,y,x_err,y_err) = clean_set(x,y,x_err,y_err)

    LM.pearson(x,y, ax)
    LM.curvefitting(x,y, ax)
    LM.linmixing(x,y,x_err,y_err, ax)

    ax.errorbar(x, y, y_err, x_err, fmt='o', color='Purple', ecolor='Black', markersize=C.SIZE/3, alpha=C.ALPHA)

    #plot labels
    ax.set_xlabel(x_name)
    ax.set_ylabel(y_name)
    ax.set_title(key)
    #formatting plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

def linmix_plots_multi(key,src,x,y,x_err,y_err,ax=plt,color='Purple',LMcolor='Red'):
    dfx=x.to_frame()
    x_name = str(dfx.columns.values)
    dfy=y.to_frame()
    y_name = str(dfy.columns.values)

    (x,y,x_err,y_err) = clean_set(x,y,x_err,y_err)

    LM.pearson(x,y,ax)
    LM.curvefitting(x,y,ax,color=LMcolor)
    LM.linmixing(x,y,x_err,y_err,ax,color=LMcolor)

    ax.errorbar(x, y, y_err, x_err, fmt='o', color=color, ecolor='Black', markersize=C.SIZE/3, alpha=C.ALPHA, label=key)

    #plot labels
    ax.set_xlabel(x_name)
    ax.set_ylabel(y_name)
    #formatting plot
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
