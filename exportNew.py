#import libraries here
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

#plotting scripts
import base_plotter as BPLT
import ratio_plots as RPLT
import plotter as PLT

#essential scripts
import gather_data as GD
import constants as C
import Milestone3 as M3

#data analysis and manipulation
import core_analysis as CA
import data_manipulation as DM
import linmixer as LM

#formatting
import formatting_functions as FF

def JINGLE():
    src = GD.get_data(GD.JINGLE_GALAXY)
    #src2 = GD.get_data(GD.JINGLE_FLUX)
    #src3 = GD.get_data(GD.JINGLE_TEMPEL)

    #src = pd.merge(src1, src2, on='JINGLEID')
    #jngl = pd.merge(src, src3, on='JINGLEID')

    src.loc[src['Z_PP04_O3N2'] < 0, 'Z_PP04_O3N2'] = np.nan
    src['LOGMDUST_DELOOZE_ERR'] = [(abs(g)+abs(h))/2 for g,h in zip(src['dex-'], src['dex+'])]

    src.loc[src['H1_FLAG'] < 0, ['LOGMH1_MATT','LOGMH1_MATT_ERR']] = np.nan
    src.loc[src['LOGMH1_TING'] < 0.01, 'LOGMH1_TING'] = np.nan

    src.loc[src['LOGMH2_RYAN'] < 0.1, 'LOGMH2_RYAN'] = np.nan
    src.loc[src['LOGMH2_TING'] < 0.01, 'LOGMH2_TING'] = np.nan

    src['LOGMGAS'] = np.log10(1.3*(np.power(10,src['LOGMH1_MATT'])+np.power(10,src['LOGMH2_TING'])))
    #1.3 comes from Chown and Catinella

    src['GAS_FLAG'] = 0
    src.loc[src['LOGMGAS'].notnull(), 'GAS_FLAG'] = 1

    G2DR = (src['LOGMGAS'] - src['LOGMDUST_DELOOZE']).mean()

    src.loc[src['GAS_FLAG'] == 0, 'LOGMGAS'] = G2DR + src['LOGMDUST_DELOOZE']
    src.loc[(src['LOGMGAS'].notnull()) & (src['GAS_FLAG'] != 1), 'GAS_FLAG'] = 2

    src['LOGMH1'] = src['LOGMH1_MATT']
    src['LOGMH2'] = src['LOGMH2_TING']

    src['H1_F'] = 0
    src['H2_F'] = 0

    src.loc[src['LOGMH1'].notnull(), 'H1_F'] = 1
    src.loc[src['LOGMH2'].notnull(), 'H2_F'] = 1
    
    h1_reverse = np.log10(np.power(10, src['LOGMGAS']) / 1.3 - np.power(10, src['LOGMH2']))
    h2_reverse = np.log10(np.power(10, src['LOGMGAS']) / 1.3 - np.power(10, src['LOGMH1']))

    src.loc[src['H1_F'] == 0, 'LOGMH1'] = h1_reverse
    src.loc[src['H2_F'] == 0, 'LOGMH2'] = h2_reverse

    src.loc[(src['LOGMH1'].notnull()) & (src['H1_F'] != 1), 'H1_F'] = 2
    src.loc[(src['LOGMH2'].notnull()) & (src['H2_F'] != 1), 'H2_F'] = 2
    
    fz = 27.36*np.float_power(10, (src['Z_PP04_O3N2'] - 12))
    src['LOGMZ'] = fz * src['LOGMGAS'] + src['LOGMDUST_DELOOZE']
    
    FF.print_full(src)

    print(len(src.index[src['H1_F']!=0].tolist()))
    print(len(src.index[src['H2_F']!=0].tolist()))

    src.to_csv('JINGLE_DATA.csv')

    #Keep Matt for H1
    #Keep Ting for H2