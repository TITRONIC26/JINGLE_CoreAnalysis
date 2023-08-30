"""
Code used explicitly for collecting and processing the data stored locally
"""

#import libraries here
import pandas as pd

#file locations go here
JINGLE_MASTER = r"C:\Users\jsmes\OneDrive\Documents\Temp\Thesis\Data\CombinedDatasets\JINGLE_MASTER.csv"
JINGLE_TEMPEL = r"C:\Users\jsmes\OneDrive\Documents\Temp\Thesis\Data\CombinedDatasets\JINGLE_TEMPEL_FULL.csv"
SED_FITTINGS = r"C:\Users\jsmes\OneDrive\Documents\Temp\Thesis\Data\CombinedDatasets\SED_FITTINGS.csv"
XCOLDGASS = r"C:\Users\jsmes\OneDrive\Documents\Temp\Thesis\Data\xCOLDGASS\xCOLDGASS_Near.csv"
VERTICO = r"C:\Users\jsmes\OneDrive\Documents\Temp\Thesis\Data\CombinedDatasets\VERTICO_DP.csv"

#header lists go here
JM_header = ['JINGLEID','IDNUM','RA','DEC','z','z_ERR','Z_PP04_N2','Z_PP04_O3N2','Z_MZR','LOGMSTAR_MAGPHYS','LOGMSTAR_MAGPHYS_ERR','LOGSFR_MAGPHYS','LOGSFR_MAGPHYS_ERR','LOGMDUST_MAGPHYS','LOGMDUST_DELOOZE','-dex','+dex','LOGMH1','LOGMH1_FRAC','LOGMH2','LOGMH2_PRED','LOGMHALO','GROUP_TYPE']
JT_header = ['IDNUM','GROUP_TYPE','IDGAL','NGAL','DGAL','MNFW','DEN1','FAMID','FAM_SIZE']
SF_header = ['IDNUM','logMc_SMBB','e_logMc_SMBB','Tc_SMBB','e_Tc_SMBB','logMc_BMBB','e_logMc_BMBB','Tc_BMBB','e_Tc_BMBB','logMc_TMBB','e_logMc_TMBB','Tc_TMBB','e_Tc_TMBB','logMw_TMBB','e_logMw_TMBB','Tw_TMBB','e_Tw_TMBB']
XCG_header = ['IDNUM','RA','DEC','z','LOGMSTAR','LOGSFR_BEST','LOGSFR_BEST_ERR','Z_PP04_N2','Z_PP04_O3N2','Z_MZR','LOGMH2','LOGMH2_ERR','LOGMH2_PRED']
V_header = ['NAME','RA','DEC','z','LOGMH2','LOGMH2_ERR','LOGMH1_DP','LOGMH1_DP_ERR','MH2_DP','MH2_DP_ERR','LOGMMETAL_DP', 'LOGMMETAL_DP_ERR', 'Z_DP_O3N2']

#functions go here
def get_header(file):
    """
    Selects the appropriate header for the selected file.
    """
    if file == JINGLE_MASTER:
        header = JM_header
    elif file == JINGLE_TEMPEL:
        header = JT_header
    elif file == SED_FITTINGS:
        header = SF_header
    elif file == XCOLDGASS:
        header = XCG_header
    elif file == VERTICO:
        header = V_header
    else:
        print("Invalid selection")
        header = []
    
    return header

def get_data(file):
    """
    Extract the data from the selected file location.
    """
    try:
        #grab the appropriate header file
        name = get_header(file)
        #extract the data from the csv
        data = pd.read_csv(file, sep=',', header=0, names=name)
        #print a summary of the key information within the dataset to confirm the selection
        #print(data.describe())
    
    except Exception as e:
        #return the error causing the issue
        print(e)
    
    return data