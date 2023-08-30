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

#data analysis and manipulation
import core_analysis as CA
import data_manipulation as DM

#formatting
import formatting_functions as FF

#ignore pandas warning
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)

#main function for analyzing data
def main():
    print("Hello and welcome to the plotting software!")
    print("-------------------------------------------")
    return

#call on the main function when the script is executed
if __name__ == '__main__':
    main()
    #CF.jingle_galaxy_base_parameters()
    #CF.gas_content_comparisons()
    #CF.specific_SFR()
    #CF.grouped_by(Env=False)
    CF.galaxy_ratios()
    #CF.compare()

    #FF.print_full(CF.src5)