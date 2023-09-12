"""
Python script aimed at modifying the data directly from the CSV file
"""

#import libraries here
import numpy as np

#import other scripts here
import core_analysis as CA

#declare global variables here

#define functions here
def JINGLE_main(src):
    #replace with null values
    src.loc[src['Z_PP04_O3N2'] < 0, 'Z_PP04_O3N2'] = np.nan
    src.loc[src['LOGMHALO_TEMPEL'] < 0, 'LOGMHALO_TEMPEL'] = np.nan
    src.loc[src['LOGMH2_RYAN'] < 0.1, 'LOGMH2_RYAN'] = np.nan
    src.loc[src['H1_FLAG'] < 0, ['LOGMH1_MATT','LOGMH1_MATT_ERR']] = np.nan

    src = CA.generate_std_error(src, src['LOGMH2_RYAN'])
    
    src['LOGMDUST_DELOOZE_ERR'] = [(abs(g)+abs(h))/2 for g,h in zip(src['dex-'], src['dex+'])]
    
    #generate the total gas masses and metals
    src = CA.find_Mmetals(src)

    return src

def VERTICO_main(src):
    #convert the non log base dustpedia masses to log base for comparing to Brown.
    src['LOGMSTAR_DP_ERR'] = (1/np.log(10)) * (src['LOGMSTAR_DP_ERR'] / src['LOGMSTAR_DP'])
    src['LOGMSTAR_DP'] = np.log10(src['LOGMSTAR_DP'])

    src['LOGSFR_DP_ERR'] = (1/np.log(10)) * (src['LOGSFR_DP_ERR'] / src['LOGSFR_DP'])
    src['LOGSFR_DP'] = np.log10(src['LOGSFR_DP'])

    src['LOGMDUST_DP_ERR'] = (1/np.log(10)) * (src['LOGMDUST_DP_ERR'] / src['LOGMDUST_DP'])
    src['LOGMDUST_DP'] = np.log10(src['LOGMDUST_DP'])

    return src

def HERACLES_main(src):
    #convert normal to log base 10 values.
    src['MH2/MH1'] = np.log10(src['MH2/MH1'])
    src['MH2/MSTAR'] = np.log10(src['MH2/MSTAR'])
    src['MGAS/MSTAR'] = np.log10(src['MGAS/MSTAR'])
    src['MG/MS_5L'] = np.log10(src['MG/MS_5L'])
    src['LOGSFR'] = np.log10(src['LOGSFR'])

    return src

def XCG_main(src):
    #make 0s NaNs
    src.loc[src['LOGMH2'] < 0.1, 'LOGMH2'] = np.nan
    src.loc[src['LOGMH2_LIM'] < 0.1, 'LOGMH2_LIM'] = np.nan

    return src