""" Last updated: 2023/05/31

Script for plotting the effect of curvature and cut length on incision 
opening (w/l). These are plotted for the change in applied stretch.

"""
import csv
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def polyRegFit(x, a, b, c, d, e, f):
    polyReg = a*(x**5)+b*(x**4)+c*(x**3)+d*(x**2)+e*(x)+f
    return polyReg

def plot(odb1, odb2, odb3, odb4, odb5, xticks, yticks, legend, name, colors):
    cmap = matplotlib.colormaps["Blues"]
    
    # figure settings   
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = "serif"
    plt.rcParams['axes.autolimit_mode'] = 'round_numbers'
    plt.rcParams['axes.xmargin'] = 0
    plt.rcParams['axes.ymargin'] = 0
    plt.rcParams['font.size'] = 13
    
    # extract data
    dataframe1 = pd.read_csv(odb1, encoding='ISO-8859-1')
    dataframe2 = pd.read_csv(odb2, encoding='ISO-8859-1')
    dataframe3 = pd.read_csv(odb3, encoding='ISO-8859-1')
    dataframe4 = pd.read_csv(odb4, encoding='ISO-8859-1')
    dataframe5 = pd.read_csv(odb5, encoding='ISO-8859-1')
    dataframeHenderson = pd.read_csv('../simulations/curvature and cut length .csv files/Henderson_2005_Flat_case.csv', encoding='ISO-8859-1')
    
    # curve fit data so that graph properly ends at 1.14 stretch
    # note: the reason the graph does not start from zero is because the nodes are not directly aligned prior to deformation
    dataframes = [dataframe1, dataframe2, dataframe3, dataframe4, dataframe5]
    endLam = 1.140
    for i in range(len(dataframes)):
        parameters, covariance = curve_fit(f=polyRegFit, xdata=dataframes[i].lambda_p, ydata=dataframes[i].wOverlHalf, bounds=(-np.inf, np.inf))
        a = parameters[0]
        b = parameters[1]
        c = parameters[2]
        d = parameters[3]
        e = parameters[4]
        f = parameters[5]
        ratio= a*(endLam**5)+b*(endLam**4)+c*(endLam**3)+d*(endLam**2)+e*(endLam)+f
        df = pd.DataFrame({'lambda_p': [endLam], 'wOverlHalf': [ratio]})
        if i == 0:
            for stre, j in zip(dataframes[i].lambda_p, dataframes[i].index):
                if stre > 1.140:
                    dataframe1 = dataframe1.drop(j)
            dataframe1 = pd.concat([dataframe1, df], ignore_index=True)
        elif i == 1:
            for stre, j in zip(dataframes[i].lambda_p, dataframes[i].index):
                if stre > 1.140:
                    dataframe2 = dataframe2.drop(j)
            dataframe2 = pd.concat([dataframe2, df], ignore_index=True)
        elif i == 2:
            for stre, j in zip(dataframes[i].lambda_p, dataframes[i].index):
                if stre > 1.140:
                    dataframe3 = dataframe3.drop(j)
            dataframe3 = pd.concat([dataframe3, df], ignore_index=True)
        elif i == 3:
            for stre, j in zip(dataframes[i].lambda_p, dataframes[i].index):
                if stre > 1.140:
                    dataframe4 = dataframe4.drop(j)
            dataframe4 = pd.concat([dataframe4, df], ignore_index=True)
        elif i == 4:
            for stre, j in zip(dataframes[i].lambda_p, dataframes[i].index):
                if stre > 1.140:
                    dataframe5 = dataframe5.drop(j)
            dataframe5 = pd.concat([dataframe5, df], ignore_index=True)

    # create plot
    plt.figure()
    # note that we have larger values (darker color) to smaller values (lighter color)
    plt.plot(dataframe1.lambda_p, dataframe1.wOverlHalf, linestyle='-', color=cmap(colors[0]), label=legend[0], linewidth=2)
    plt.plot(dataframe2.lambda_p, dataframe2.wOverlHalf, linestyle='-', color=cmap(colors[1]), label=legend[1], linewidth=2)
    plt.plot(dataframe3.lambda_p, dataframe3.wOverlHalf, linestyle='-', color=cmap(colors[2]), label=legend[2], linewidth=2)
    plt.plot(dataframe4.lambda_p, dataframe4.wOverlHalf, linestyle='-', color=cmap(colors[3]), label=legend[3], linewidth=2)
    if len(legend) == 6:
        plt.plot(dataframe5.lambda_p, dataframe5.wOverlHalf, linestyle='-',color=cmap(colors[4]), label=legend[4], linewidth=2)
        plt.plot(dataframeHenderson.lambda_p, dataframeHenderson.wOverlHalf, linestyle='-',color=cmap(colors[5]), label=legend[5], linewidth=2)

    plt.xlabel(r'$\lambda_p$')
    plt.ylabel(r'$\overline{w}_S$')
    plt.legend()
    
    plt.savefig(name, dpi=400)

def transverseCurvatureVsCutOpening():
    odba2Tran  = '../simulations/curvature and cut length .csv files/a=2MouseDuraMaterTranCutResults.csv'
    odba4Tran  = '../simulations/curvature and cut length .csv files/a=4MouseDuraMaterTranCutResults.csv'
    odba6Tran  = '../simulations/curvature and cut length .csv files/a=6MouseDuraMaterTranCutResults.csv'
    odba8Tran  = '../simulations/curvature and cut length .csv files/a=8MouseDuraMaterTranCutResults.csv'
    odba10Tran = '../simulations/curvature and cut length .csv files/a=10MouseDuraMaterTranCutResults.csv'
    legend = [r'$a=2$', r'$a=4$', r'$a=6$', r'$a=8$', r'$a=10$', 'Henderson et al.']
    name = 'figure-curvature_transverse.png'
    xticks = np.arange(1.000, 1.160, step=0.020)
    yticks = np.arange(0.000, 0.500, step=0.050)
    colors = [0.25, 0.40, 0.55, 0.70, 0.85, 1.00]
    plot(odba2Tran, odba4Tran, odba6Tran, odba8Tran, odba10Tran, xticks, yticks, legend, name, colors)

def longitudinalCurvatureVsCutOpening():
    odba2Long  = '../simulations/curvature and cut length .csv files/a=2MouseDuraMaterLongCutResults.csv'
    odba4Long  = '../simulations/curvature and cut length .csv files/a=4MouseDuraMaterLongCutResults.csv'
    odba6Long  = '../simulations/curvature and cut length .csv files/a=6MouseDuraMaterLongCutResults.csv'
    odba8Long  = '../simulations/curvature and cut length .csv files/a=8MouseDuraMaterLongCutResults.csv'
    odba10Long = '../simulations/curvature and cut length .csv files/a=10MouseDuraMaterLongCutResults.csv'
    legend = [r'$a=2$', r'$a=4$', r'$a=6$', r'$a=8$', r'$a=10$', 'Henderson et al.']
    name = 'figure-curvature_longitudinal.png'
    xticks = np.arange(1.000, 1.160, step=0.020)
    yticks = np.arange(0.000, 0.500, step=0.050)
    colors = [0.25, 0.40, 0.55, 0.70, 0.85, 1.00]
    plot(odba2Long, odba4Long, odba6Long, odba8Long, odba10Long, xticks, yticks, legend, name, colors)

def transverseCutLengthVsCutOpening():
    odbHalfTran = '../simulations/curvature and cut length .csv files/adultICRMouseDuraMaterTranCutHalfAxisResults.csv'
    odb4thTran  = '../simulations/curvature and cut length .csv files/adultICRMouseDuraMaterTranCut4thAxisResults.csv'
    odb8thTran  = '../simulations/curvature and cut length .csv files/adultICRMouseDuraMaterTranCut8thAxisResults.csv'
    odb16thTran = '../simulations/curvature and cut length .csv files/adultICRMouseDuraMaterTranCut16thAxisResults.csv'
    odb32ndTran = '../simulations/curvature and cut length .csv files/adultICRMouseDuraMaterTranCut32ndAxisResults.csv'
    legend = [r'$l_0=b/2$', r'$l_0=b/4$', r'$l_0=b/8$', r'$l_0=b/16$', r'$l_0=b/32$']
    name = 'figure-cutLength_transverse.png'
    xticks = np.arange(1.000, 1.160, step=0.020)
    yticks = np.arange(0.000, 0.500, step=0.050)
    colors = [1.00, 0.80, 0.60, 0.40]
    plot(odbHalfTran, odb4thTran, odb8thTran, odb16thTran, odb32ndTran, xticks, yticks, legend, name, colors)

def longitudinalCutLengthVsCutOpening():
    odbHalfLong = '../simulations/curvature and cut length .csv files/adultICRMouseDuraMaterLongCutHalfAxisResults.csv'
    odb4thLong  = '../simulations/curvature and cut length .csv files/adultICRMouseDuraMaterLongCut4thAxisResults.csv'
    odb8thLong  = '../simulations/curvature and cut length .csv files/adultICRMouseDuraMaterLongCut8thAxisResults.csv'
    odb16thLong = '../simulations/curvature and cut length .csv files/adultICRMouseDuraMaterLongCut16thAxisResults.csv'
    odb32ndLong = '../simulations/curvature and cut length .csv files/adultICRMouseDuraMaterLongCut32ndAxisResults.csv'
    legend = [r'$l_0=a/2$', r'$l_0=a/4$', r'$l_0=a/8$', r'$l_0=a/16$', r'$l_0=a/32$']
    name = 'figure-cutLength_longitudinal.png'
    xticks = np.arange(1.000, 1.160, step=0.020)
    yticks = np.arange(0.000, 0.500, step=0.050)
    colors = [1.00, 0.80, 0.60, 0.40]
    plot(odbHalfLong, odb4thLong, odb8thLong, odb16thLong, odb32ndLong, xticks, yticks, legend, name, colors)

if __name__ == '__main__':
    
    transverseCurvatureVsCutOpening()
    longitudinalCurvatureVsCutOpening()
    transverseCutLengthVsCutOpening()
    longitudinalCutLengthVsCutOpening()
