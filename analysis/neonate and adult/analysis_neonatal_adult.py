""" Last updated: 2023/05/31

Script for plotting the change in incision opening (w/l) as stretch increases,
for the transverse and longitudinal cuts in the neonatal and adult murine 
cranial dura mater simulations.

"""
import csv
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def experimentalStretch(stretchFit, wOverl, ave):
    global stretch
    for j in stretchFit.index:
        if stretchFit[j] < stretchFit.iloc[-1]:
            if ave >= wOverl[j] and ave <= wOverl[j+1]:
                stretch = stretchFit[j]+(stretchFit[j+1]-stretchFit[j])*((ave-wOverl[j])/(wOverl[j+1]-wOverl[j]))

def polyRegFit(x, a, b, c, d, e, f):
    polyReg = a*(x**5)+b*(x**4)+c*(x**3)+d*(x**2)+e*(x)+f
    return polyReg

def stretchVsIncisionOpening(expWoverL, odb, title, xticks, yticks):
    # figure settings 
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = "serif"
    plt.rcParams['axes.autolimit_mode'] = 'round_numbers'
    plt.rcParams['axes.xmargin'] = 0
    plt.rcParams['axes.ymargin'] = 0
    plt.rcParams['font.size'] = 13


    # extract stretch and w/l data
    dataframe = pd.read_csv(odb, encoding='ISO-8859-1')

    # curve fit of adult transverse data for standard deviation higher bound of stretch
    if title == 'adult_transverse':
        orginalDFLam = dataframe.lambda_p
        orginalDFRatio = dataframe.wOverlHalf
        parameters, covariance = curve_fit(f=polyRegFit, xdata=orginalDFLam, ydata=orginalDFRatio, bounds=(-np.inf, np.inf))
        a = parameters[0]
        b = parameters[1]
        c = parameters[2]
        d = parameters[3]
        e = parameters[4]
        f = parameters[5]
        stretches = [1.135, 1.140, 1.145, 1.150, 1.155, 1.16]
        wlRatio = []
        for lam in stretches:
            ratio = a*(lam**5)+b*(lam**4)+c*(lam**3)+d*(lam**2)+e*(lam)+f
            wlRatio.append(ratio)
        df_add = pd.DataFrame({'lambda_p': stretches, 'wOverlHalf': wlRatio})
        dataframe = pd.concat([dataframe, df_add], ignore_index=True)

    # estimate pre-stretch
    global expStAve 
    global expStStd
    experimentalStretch(dataframe.lambda_p, dataframe.wOverlHalf, expWoverL[0])
    expStAve = stretch
    experimentalStretch(dataframe.lambda_p, dataframe.wOverlHalf, (expWoverL[0] - expWoverL[1]))
    expStStdNeg = np.abs(expStAve - stretch)
    experimentalStretch(dataframe.lambda_p, dataframe.wOverlHalf, (expWoverL[0] + expWoverL[1]))
    expStStdPos = np.abs(expStAve - stretch)
    expStStd = np.mean([expStStdNeg,expStStdPos])

    # plot
    plt.figure()
    if title == 'adult_transverse':
        plt.plot(orginalDFLam, orginalDFRatio, linestyle='-',color='k', linewidth=2)
    else:
        plt.plot(dataframe.lambda_p, dataframe.wOverlHalf, linestyle='-',color='k', linewidth=2)
    plt.plot([1.00, expStAve], [expWoverL[0], expWoverL[0]], linestyle='--', color='k', linewidth=1)
    plt.plot([expStAve,expStAve], [0, expWoverL[0]], linestyle='--', color='k', linewidth=1)
    plt.axhspan((expWoverL[0]-expWoverL[1]), (expWoverL[0]+expWoverL[1]), facecolor='#DCDCDC')
    plt.fill_between(dataframe.lambda_p, dataframe.wOverlHalf, color='#7F7F7F', alpha=1)
    plt.fill_between([1.00, (expStAve-expStStdNeg)], [(expWoverL[0]-expWoverL[1]), (expWoverL[0]-expWoverL[1])], color='white', alpha=1)
    plt.fill_between([(expStAve+expStStdPos), 1.16], [(expWoverL[0]+expWoverL[1]), 0.35], color='white', alpha=1)
    plt.xlabel(r'$\lambda_p$')
    plt.ylabel(r'$\overline{w}_S$')
    plt.savefig("figure-" + title, dpi=400)

def plot(sim):
    xticks = np.arange(1.000, 1.180, step=0.020)
    yticks = np.arange(0.000, 0.400, step=0.050)
    stretchVsIncisionOpening(sim.expWoverL, sim.odb, sim.title, xticks, yticks)
    print('Simulation:', sim.title)
    print('Estimated stretch:', expStAve,'+/-',expStStd)

def run(sims):
    # run simulations
    experiments = []
    expWoverLAve = []
    expWoverLStd = []
    expStretchAve = []
    expStretchStd = []
    for sim in sims:
        plot(sim)
        experiments.append(sim.title)
        expWoverLAve.append(sim.expWoverL[0])
        expWoverLStd.append(sim.expWoverL[1])
        expStretchAve.append(expStAve)
        expStretchStd.append(expStStd)
    # write experimental stretch information to file
    outFile = open('results-estimated_stretches.csv', 'w+')
    header = ['Experiment','w/l (exp, ave)','w/l (exp, std)','stretch (ave)','stretch (std)']
    output = zip(experiments, expWoverLAve, expWoverLStd, expStretchAve, expStretchStd)
    writer = csv.writer(outFile)
    writer.writerow(header)
    writer.writerows(output)
    outFile.close()

class neonateTran():
    def __init__(self):
        self.expWoverL = [0.171, 0.073]
        self.odb = '../simulations/neonate and adult .csv files/neonateICRMouseDuraMaterTranCutResults.csv'
        self.title = 'neonate_transverse'
class neonateLong():
    def __init__(self):
        self.expWoverL = [0.106, 0.045]
        self.odb = '../simulations/neonate and adult .csv files/neonateICRMouseDuraMaterLongCutResults.csv'
        self.title = 'neonate_longitudinal'
class adultTran():
    def __init__(self):
        self.expWoverL = [0.244, 0.072]
        self.odb = '../simulations/neonate and adult .csv files/adultICRMouseDuraMaterTranCutResults.csv'
        self.title = 'adult_transverse'
class adultLong():
    def __init__(self):
        self.expWoverL = [0.094, 0.043]
        self.odb = '../simulations/neonate and adult .csv files/adultICRMouseDuraMaterLongCutResults.csv'
        self.title = 'adult_longitudinal'

if __name__ == '__main__':
    sims = [neonateTran(), neonateLong(), adultTran(), adultLong()]
    run(sims)
