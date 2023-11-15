"""
Main function. Run this code to access the data from the local JINGLE csv files.
This code is meant to give a "quick" look at the data within the JINGLE files, showing a brief overlook at the main information.
"""

#import libraries here
import warnings
from pandas.errors import SettingWithCopyWarning

#plotting scripts
import base_plotter as BPLT
import ratio_plots as RPLT

#essential scripts
import gather_data as GD
import callable_functions as CF
import constants as C
import Final as F
import CatinellaPlots as CP

#data analysis and manipulation
import core_analysis as CA
import data_manipulation as DM

#formatting
import formatting_functions as FF

#ignore pandas warning
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)

#Variables and Data
JINGLE = r"C:\Users\jsmes\OneDrive\Documents\Temp\Thesis\Data\FINAL\JINGLE_DATA.csv"
jingle = GD.get_data(JINGLE, keep=True)

xcoldgass = CF.src6.copy()

#main function for analyzing data
def main():
    print("Hello and welcome to the plotting software!")
    print("-------------------------------------------")
    
    #preamble
    UL = CF.main()
    jingle['H1_F'][UL] = 3

    #CP.Figure1(jingle)
    #CP.Figure5(jingle)
    #CP.Figure6(jingle, xcg=True)
    #CP.Figure7(jingle, xcg=False)
    #CP.Figure8(jingle, xcg=False)
    #CP.Figure9(jingle, xcg=False)
    #CP.Figure10(jingle, xcg=False)
    CP.Figure11(jingle, xcg=False)


    #FF.print_full(jingle)

    #FF.print_full(jingle['LOGMH1']-jingle['LOGMSTAR_MAGPHYS'])
    #print(len(jingle['LOGMH1']-jingle['LOGMSTAR_MAGPHYS']))
    #FF.print_full(xcoldgass)

    #F.sSFR_MS(jingle, xcg = True, seven = False, RMOL = False, DUST = True, Mean = True, Median = True)
    #F.Catinella_1(jingle)
    #F.Catinella_2(jingle, xcg = False, seven = False, RMOL=True, DUST=False, Mean=False, Median=False)

    #F.Catinella_3(jingle)
    #F.Catinella_4(jingle)
    #F.Catinella_5(jingle, seven=False)

    return
    



#call on the main function when the script is executed
if __name__ == '__main__':
    main()
    