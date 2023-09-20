"""
Code used explicitly for collecting and processing the data stored locally
"""

#import libraries here
import pandas as pd

#file locations go here
JINGLE_FLUX = r"C:\Users\jsmes\OneDrive\Documents\Temp\Thesis\Data\Combined Datasets\JINGLE_FluxMeasurements.csv"
JINGLE_GALAXY = r"C:\Users\jsmes\OneDrive\Documents\Temp\Thesis\Data\Combined Datasets\JINGLE_GalaxyProperties.csv"
JINGLE_TEMPEL = r"C:\Users\jsmes\OneDrive\Documents\Temp\Thesis\Data\Combined Datasets\JINGLE_TEMPEL.csv"
VERTICO_DP = r"C:\Users\jsmes\OneDrive\Documents\Temp\Thesis\Data\Combined Datasets\VERTICO_DP.csv"
XCOLD_GASS = r"C:\Users\jsmes\OneDrive\Documents\Temp\Thesis\Data\Combined Datasets\xCOLD_GASS.csv"
VERTICO = r"C:\Users\jsmes\OneDrive\Documents\Temp\Thesis\Data\Combined Datasets\VERTICO.csv"
HERACLES = r"C:\Users\jsmes\OneDrive\Documents\Temp\Thesis\Data\Combined Datasets\HERACLES.csv"

#minor file locations go here
RYAN_ISM = r"c:\Users\jsmes\OneDrive\Documents\Temp\Thesis\Data\Special_Datasets\Ryan_ISM_Masses.csv"
DELOOZE = r"C:\Users\jsmes\OneDrive\Documents\Temp\Thesis\Data\Special_Datasets\DeLooze_DustMasses.csv"


#header lists go here
JF_header = ['JINGLEID','IDNUM','FLAG','SCUBA2_850_PeakSNR','SCUBA2_850_Flux','SCUBA2_850_Err',
             'DETFLAG_850','SN850','F850_CORR','F850ERR_CORR','ICO21','SNCO21','ALPHA_CO']

JG_header = ['JINGLEID','RA','DEC','z','z_ERR','Z_PP04_O3N2','LOGMSTAR_MAGPHYS',
             'LOGMSTAR_MAGPHYS_ERR','LOGSFR_MAGPHYS','LOGSFR_MAGPHYS_ERR','LOGMHALO_TEMPEL',
             'GROUPRANK_TEMPEL','GROUPCLASS','LOGMDUST_DELOOZE','dex-','dex+','LOGMH2_RYAN',
             'LOGMH1_MATT','LOGMH1_MATT_ERR','H1_FLAG']

JT_header = ['JINGLEID','GROUPID','NGAL','Rank','Dist.g','pE','pS0','pSa','pSc','RAgroup',
             'DEgroup','zgroup','sig.v','Rvir','Rmax','LOGMHALO_TEMPEL','LOGMHALO_TEMPEL_ERR','Den1','Den2','Den4',
             'Den8']

V_header = ['ID','RA','DE','Vel','S/N','L_LOGMH2','LOGMH2','LOGMH2_ERR','LOGSFR_DP',
            'LOGSFR_DP_ERR','LOGMSTAR_DP','LOGMSTAR_DP_ERR','LOGMDUST_DP','LOGMDUST_DP_ERR',
            'LOGMH1_DP','LOGMH1_DP_ERR','H1_Flag','Z_O3N2','Z_O3N2_ERRDWN','Z_O3N2_ERRUP',
            'LOGMZ_DP','LOGMZ_DP_ERRDWN','LOGMZ_DP_ERRUP']

XCG_header = ['ID','RA','DEC','z','LOGMSTAR','LOGSFR','LOGSFR_ERR','Z_PP04_O3N2','LOGMH2','LOGMH2_ERR','LOGMH2_LIM','LOGMH1','LOGMH1_Q']

VERT_header = ['ID','RA','DEC','z','LOGMH1','LOGMSTAR','LOGMSTAR_ERR','LOGSFR','LOGSFR_ERR','LOGMH2','LOGMH2_ERR']

HER_header = ['ID','RA','DEC','Vel','z','LOGMSTAR','LOGMSTAR_ERR','LOGSFR','LOGSFR_ERR','LOGMH2','LOGMH2_ERR','LOGMH1','LOGMH1_ERR']

#functions go here
def get_header(file):
    """
    Selects the appropriate header for the selected file.
    """
    if file == JINGLE_FLUX:
        header = JF_header
    elif file == JINGLE_TEMPEL:
        header = JT_header
    elif file == JINGLE_GALAXY:
        header = JG_header
    elif file == VERTICO_DP:
        header = V_header
    elif file == XCOLD_GASS:
        header = XCG_header
    elif file == VERTICO:
        header = VERT_header
    elif file == HERACLES:
        header = HER_header
    else:
        print("Invalid selection")
        header = []
    
    return header

def get_data(file, keep=False):
    """
    Extract the data from the selected file location.
    """
    try:
        if keep == True:
            #extract the data from the csv
            data = pd.read_csv(file, sep=',')
        else:
            #grab the appropriate header file
            name = get_header(file)
            #extract the data from the csv
            data = pd.read_csv(file, sep=',', header=0, names=name)
    
    except Exception as e:
        #return the error causing the issue
        print(e)
    
    return data