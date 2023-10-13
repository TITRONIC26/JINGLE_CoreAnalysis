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
    src.loc[src['LOGMH2_RYAN'] < 0.1, 'LOGMH2_RYAN'] = np.nan
    src.loc[src['H1_FLAG'] < 0, ['LOGMH1_MATT','LOGMH1_MATT_ERR']] = np.nan

    src = CA.generate_std_error(src, src['LOGMH2_RYAN'])
    
    src['LOGMDUST_DELOOZE_ERR'] = [(abs(g)+abs(h))/2 for g,h in zip(src['dex-'], src['dex+'])]
    
    #generate the total gas masses and metals
    src = CA.find_Mmetals(src)

    return src

def VERTICO_main(src):
    #convert the non log base dustpedia masses to log base for comparing to Brown.
    src = CA.generate_std_error(src, src['LOGMH1'])

    return src

def HERACLES_main(src):
    #nothing yet
    

    return src

def XCG_main(src):
    src['IDNUM'] = src.index
    #make 0s NaNs
    src.loc[src['LOGMH2'] < 0.1, 'LOGMH2'] = np.nan
    src.loc[src['LOGMH2_LIM'] < 0.1, 'LOGMH2_LIM'] = np.nan

    src = CA.generate_std_error(src, src['LOGMH1'])

    return src