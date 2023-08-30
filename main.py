"""
Main function. Run this code to access the data from the local JINGLE csv files.
This code is meant to give a "quick" look at the data within the JINGLE files, showing a brief overlook at the main information.
"""

#import libraries here
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy as sci
import math as mt
import warnings

from pandas.errors import SettingWithCopyWarning

#import other scripts here
import gather_data as GD
import base_plotter as BPLT
import core_analysis as CA
import formatting_functions as FF
import ratio_plots as RPLT
import constants as C
import data_manipulation as DM
import callable_functions as CF

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
    #CF.galaxy_ratios()
    CF.compare()