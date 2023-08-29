"""
Script to provide the LINMIX function operation on a specified dataset with errors
"""

#import libraries here
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy as sci
import math as mt
import linmix

from scipy.optimize import curve_fit as cf
from scipy import stats

#import scripts here
import constants as C

#global variable here

#define functions here
def linmixing(x, y, x_err, y_err, axs=plt, color='Red'):
    lm = linmix.LinMix(x, y, x_err, y_err, K=2, seed=2)
    lm.run_mcmc(silent=True)

    xs = np.linspace(x.min()-1, x.max()+1, 10000)
    y_high = []
    y_low = []
    
    for x in xs:
        ys = lm.chain['alpha'] + x * lm.chain['beta']
        y_high.append(ys.max())
        y_low.append(ys.min())
    
    axs.fill_between(xs, y_high, y_low, color=color, alpha=0.2)

def func(x,m,b):
    return m*x+b

def curvefitting(x,y, axs=plt, color='Black'):
    popt,pcov = cf(func, x, y, maxfev=10000)
    xdata = np.linspace(x.min()-1, x.max()+1, 100)
    val1 = str("{0:.3g}".format(popt[0]))
    val2 = str("{0:.1g}".format(pcov[0,0]))
    eqn = 'm = '+val1+r' $\pm$ '+val2

    axs.plot(xdata, func(xdata, *popt), '--', color=color, label=eqn)
    axs.set_xlim(x.min()-0.25, x.max()+0.25)
    axs.set_ylim(y.min()-0.25, y.max()+0.25)

def pearson(x,y, axs=plt):
    results = stats.pearsonr(x,y)
    p_val = str("{0:.1g}".format(results[1]))

    axs.scatter(x.mean(), y.mean(), color='White', s=1, label='p-value = '+p_val)