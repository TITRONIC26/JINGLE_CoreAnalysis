"""
Functions for doing common mathematics
"""

#import libraries here
import pandas as pd
import numpy as np
import matplotlib as mpl
import scipy as sci
import math as mt

#import scripts here
import constants as C
import formatting_functions as FF

#define functions here
def addErrors(a,b):
    z = np.sqrt(np.float_power(a,2) + np.float_power(b,2))

    return z

def error_add(A, a, B, b):
    Z = A + B
    z = addErrors(a,b)

    return (Z,z)

def error_subtract(A, a, B, b):
    Z = A - B
    z = addErrors(a,b)

    return (Z,z)
