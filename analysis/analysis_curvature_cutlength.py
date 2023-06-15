""" Last updated: 2023/06/15

Script for plotting the effect of curvature and cut length on incision 
opening (w/l). These are plotted for the change in applied stretch.

"""
import csv
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def plot(path, csv1, csv2, csv3, csv4, csv5, ylimit, legend, name, colors):
        
    # extract data
    dataframe1 = pd.read_csv(path+csv1, encoding='ISO-8859-1')
    dataframe2 = pd.read_csv(path+csv2, encoding='ISO-8859-1')
    dataframe3 = pd.read_csv(path+csv3, encoding='ISO-8859-1')
    dataframe4 = pd.read_csv(path+csv4, encoding='ISO-8859-1')
    dataframe5 = pd.read_csv(path+csv5, encoding='ISO-8859-1')
    dataframeHenderson = pd.read_csv(path+'Henderson_2005_data.csv', encoding='ISO-8859-1')

    # create plot
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = "serif"
    plt.rcParams['axes.autolimit_mode'] = 'round_numbers'
    plt.rcParams['axes.xmargin'] = 0
    plt.rcParams['axes.ymargin'] = 0
    plt.rcParams['font.size'] = 13
    cmap = plt.cm.get_cmap("Blues")

    plt.figure()
    # note that we have larger values (darker color) to smaller values (lighter color)
    plt.plot(dataframe1.lambda_p, dataframe1.wbar, linestyle='-', color=cmap(colors[0]), label=legend[0], linewidth=2)
    plt.plot(dataframe2.lambda_p, dataframe2.wbar, linestyle='-', color=cmap(colors[1]), label=legend[1], linewidth=2)
    plt.plot(dataframe3.lambda_p, dataframe3.wbar, linestyle='-', color=cmap(colors[2]), label=legend[2], linewidth=2)
    plt.plot(dataframe4.lambda_p, dataframe4.wbar, linestyle='-', color=cmap(colors[3]), label=legend[3], linewidth=2)
    if len(legend) == 6:
        plt.plot(dataframe5.lambda_p, dataframe5.wbar, linestyle='-',color=cmap(colors[4]), label=legend[4], linewidth=2)
        plt.plot(dataframeHenderson.lambda_p, dataframeHenderson.wbar, linestyle='-',color=cmap(colors[5]), label=legend[5], linewidth=2)

    plt.xlim([1, 1.14])
    plt.ylim(ylimit)
    plt.xlabel(r'$\lambda_p$')
    plt.ylabel(r'$\overline{w}_S$')

    # reverse order of the legend to match the plot
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend(reversed(handles), reversed(labels))
    
    plt.savefig(name, dpi=400)

def transverseCurvatureVsCutOpening():
    path = '../simulations/results/'
    csva2Tran  = 'curvature_2_tran.csv'
    csva4Tran  = 'curvature_4_tran.csv'
    csva6Tran  = 'curvature_6_tran.csv'
    csva8Tran  = 'curvature_8_tran.csv'
    csva10Tran = 'curvature_10_tran.csv'
    legend = [r'$a=2$', r'$a=4$', r'$a=6$', r'$a=8$', r'$a=10$', 'Henderson et al.']
    name = 'figure-curvature_transverse.png'
    ylimit = [0.00, 0.45]
    colors = [0.25, 0.40, 0.55, 0.70, 0.85, 1.00]
    plot(path, csva2Tran, csva4Tran, csva6Tran, csva8Tran, csva10Tran, ylimit, legend, name, colors)

def longitudinalCurvatureVsCutOpening():
    path = '../simulations/results/'
    csva2Long  = 'curvature_2_long.csv'
    csva4Long  = 'curvature_4_long.csv'
    csva6Long  = 'curvature_6_long.csv'
    csva8Long  = 'curvature_8_long.csv'
    csva10Long = 'curvature_10_long.csv'
    legend = [r'$a=2$', r'$a=4$', r'$a=6$', r'$a=8$', r'$a=10$', 'Henderson et al.']
    name = 'figure-curvature_longitudinal.png'
    ylimit = [0.00, 0.45]
    colors = [0.25, 0.40, 0.55, 0.70, 0.85, 1.00]
    plot(path, csva2Long, csva4Long, csva6Long, csva8Long, csva10Long, ylimit, legend, name, colors)

def transverseCutLengthVsCutOpening():
    path = '../simulations/results/'
    csvHalfTran = 'cut-length_2_tran.csv'
    csv4thTran  = 'cut-length_4_tran.csv'
    csv8thTran  = 'cut-length_8_tran.csv'
    csv16thTran = 'cut-length_16_tran.csv'
    csv32ndTran = 'cut-length_32_tran.csv'
    legend = [r'$l_0=b/2$', r'$l_0=b/4$', r'$l_0=b/8$', r'$l_0=b/16$', r'$l_0=b/32$']
    name = 'figure-cutLength_transverse.png'
    ylimit = [0.00, 0.35]
    colors = [1.00, 0.80, 0.60, 0.40]
    plot(path, csvHalfTran, csv4thTran, csv8thTran, csv16thTran, csv32ndTran, ylimit, legend, name, colors)

def longitudinalCutLengthVsCutOpening():
    path = '../simulations/results/'
    csvHalfLong = 'cut-length_2_long.csv'
    csv4thLong  = 'cut-length_4_long.csv'
    csv8thLong  = 'cut-length_8_long.csv'
    csv16thLong = 'cut-length_16_long.csv'
    csv32ndLong = 'cut-length_32_long.csv'
    legend = [r'$l_0=a/2$', r'$l_0=a/4$', r'$l_0=a/8$', r'$l_0=a/16$', r'$l_0=a/32$']
    name = 'figure-cutLength_longitudinal.png'
    ylimit = [0.00, 0.35]
    colors = [1.00, 0.80, 0.60, 0.40]
    plot(path, csvHalfLong, csv4thLong, csv8thLong, csv16thLong, csv32ndLong, ylimit, legend, name, colors)

if __name__ == '__main__':
    
    transverseCurvatureVsCutOpening()
    longitudinalCurvatureVsCutOpening()
    transverseCutLengthVsCutOpening()
    longitudinalCutLengthVsCutOpening()
