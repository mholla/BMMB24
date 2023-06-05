import csv
import numpy as np
import pandas as pd
import statistics as st
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import ttest_ind


def interpolateStretch(stretchFit, wbar, ave):
    """ interpolates to find pre-stretch that corresponds to incision opening """
    for j in stretchFit.index:
        if stretchFit[j] < stretchFit.iloc[-1]:
            if ave >= wbar[j] and ave <= wbar[j+1]:
                stretch = stretchFit[j]+(stretchFit[j+1]-stretchFit[j])*((ave-wbar[j])/(wbar[j+1]-wbar[j]))
    return stretch


def polyfit(x, a, b, c, d, e, f):
    """ polynomial regression fit """
    polyfit = a*(x**5)+b*(x**4)+c*(x**3)+d*(x**2)+e*(x)+f
    return polyfit