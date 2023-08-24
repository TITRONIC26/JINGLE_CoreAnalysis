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
def linmixing(x, y, x_err, y_err):
    lm = linmix.LinMix(x, y, x_err, y_err, K=2, seed=2)
    lm.run_mcmc(silent=True)

    xs = np.linspace(x.min()-0.25, x.max()+0.25, 10000)
    y_high = []
    y_low = []
    
    for x in xs:
        ys = lm.chain['alpha'] + x * lm.chain['beta']
        y_high.append(ys.max())
        y_low.append(ys.min())
    
    plt.fill_between(xs, y_high, y_low, color='Red', alpha=0.2)

def func(x,m,b):
    return m*x+b

def curvefitting(x,y):
    popt,pcov = cf(func, x, y, maxfev=10000)
    xdata = np.linspace(x.min()-0.25, x.max()+0.25, 100)
    val1 = str("{0:.3g}".format(popt[0]))
    val2 = str("{0:.1g}".format(pcov[0,0]))
    eqn = 'm = '+val1+r' $\pm$ '+val2

    plt.plot(xdata, func(xdata, *popt), '--', color='Black', label=eqn)

def pearson(x,y):
    results = stats.pearsonr(x,y)
    p_val = str("{0:.1g}".format(results[1]))

    plt.scatter(x.mean(), y.mean(), color='White', label='p-value = '+p_val)

