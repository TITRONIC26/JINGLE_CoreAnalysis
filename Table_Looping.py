#import libraries here
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib
import numpy as np

from scipy.stats import kstest
from scipy.stats import ks_2samp

#plotting scripts
import base_plotter as BPLT
import ratio_plots as RPLT
import plotter as PLT

#essential scripts
import gather_data as GD
import constants as C
import Milestone3 as M3
import callable_functions as CF
import main as MAIN

#data analysis and manipulation
import core_analysis as CA
import data_manipulation as DM
import linmixer as LM

#formatting
import formatting_functions as FF

def print_table(src):
    FF.print_full(src)

def H1s(src):
    detected = src.index[src['H1_F'] == 1].tolist()
    upperlimit = src.index[src['H1_F'] == 3].tolist()
    derived = src.index[src['H1_F'] == 2].tolist()
    null = src.index[src['H1_F'] == 0].tolist()

    FF.print_full(src['LOGMH1'][detected])
    FF.print_full(src['LOGMH1'][upperlimit])
    FF.print_full(src['LOGMH1'][derived])
    FF.print_full(src['LOGMH1'][null])

    print(len(detected))
    print(len(upperlimit))
    print(len(derived))
    print(len(null))

def H2s(src):
    detected = src.index[src['H2_F'] == 1].tolist()
    derived = src.index[src['H2_F'] == 2].tolist()
    null = src.index[src['H2_F'] == 0].tolist()

    FF.print_full(src['LOGMH2'][detected])
    FF.print_full(src['LOGMH2'][derived])
    FF.print_full(src['LOGMH2'][null])

    print(len(detected))
    print(len(derived))
    print(len(null))

def Gas(src):
    detected = src.index[(src['GAS_FLAG'] == 1)].tolist()
    derived = src.index[src['GAS_FLAG'] == 2].tolist()
    null = src.index[src['GAS_FLAG'] == 0].tolist()

    FF.print_full(src['LOGMGAS'][detected])
    FF.print_full(src['LOGMGAS'][derived])
    FF.print_full(src['LOGMGAS'][null])

    print(len(detected))
    print(len(derived))
    print(len(null))

    return (detected, derived)

def Dust(src):
    detected = src.index[(src['LOGMDUST'].notnull())].tolist()
    null = src.index[src['LOGMDUST'].isnull()].tolist()

    FF.print_full(src['LOGMGAS'][detected])
    FF.print_full(src['LOGMGAS'][null])

    print(len(detected))
    print(len(null))


def GDR(src):
    vals = Gas(src)
    detected = (src['LOGMGAS'] - src['LOGMDUST'])[vals[0]].mean()
    derived = (src['LOGMGAS'] - src['LOGMDUST'])[vals[1]].mean()

    print(detected)
    print(derived)

def find_dust(src):
    gas = np.log10(1.3*(np.power(10,src['LOGMH1'])+np.power(10,src['LOGMH2'])))

    gas_diff = src['LOGMGAS'] - gas
    #FF.print_full(gas_diff)

    G2DR = src['LOGMGAS'] - src['LOGMDUST']

    dusts = gas-G2DR

    FF.print_full(src['LOGMDUST'] - dusts)





if __name__ == '__main__':
    src = MAIN.jingle.copy()
    UL = CF.main()
    src['H1_F'][UL] = 3

    print_table(src)
    #find_dust(src)
    Dust(src)