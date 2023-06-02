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


def extrapolate(df_sim):
    # fit polynomial to simulation data
    lambda_orig = df_sim.lambda_p
    wbar_orig = df_sim.wbar
    [[a, b, c, d, e, f], covariance] = curve_fit(f=polyfit, xdata=lambda_orig, ydata=wbar_orig, bounds=(-np.inf, np.inf))
    
    # extrapolate remaining points
    lambda_extrapolated = [1.135, 1.140, 1.145, 1.150, 1.155, 1.16]
    wbar_extrapolated = np.zeros(len(lambda_extrapolated))
    for i in range(len(lambda_extrapolated)):
        wbar_extrapolated[i] = polyfit(lambda_extrapolated[i], a, b, c, d, e, f)
    
    # add to simulation data
    df_ext = pd.DataFrame({'lambda_p': lambda_extrapolated, 'wbar': wbar_extrapolated})
    df_sim = pd.concat([df_sim, df_ext], ignore_index=True)

    return df_sim
