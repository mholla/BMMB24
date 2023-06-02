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

import analysis_functions

def stretchVsIncisionOpening(wbar_exp, path, odb, title, xticks, yticks):

    # extract stretch and w/l data
    df_sim = pd.read_csv(path+odb, encoding='ISO-8859-1')

    # curve fit of adult transverse data for standard deviation higher bound of stretch
    if title == 'adult_transverse':
        df_sim = analysis_functions.extrapolate(df_sim)

    # estimate pre-stretch
    lambda_avg = analysis_functions.interpolateStretch(df_sim.lambda_p, df_sim.wbar, wbar_exp[0])
    lambda_low = analysis_functions.interpolateStretch(df_sim.lambda_p, df_sim.wbar, (wbar_exp[0] - wbar_exp[1]))
    lambda_high = analysis_functions.interpolateStretch(df_sim.lambda_p, df_sim.wbar, (wbar_exp[0] + wbar_exp[1]))
    
    lambda_std_low = np.abs(lambda_avg - lambda_low)
    lambda_std_high = np.abs(lambda_avg - lambda_high)
    lambda_std = np.mean([lambda_std_low,lambda_std_high])

    # plot
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = "serif"
    plt.rcParams['axes.autolimit_mode'] = 'round_numbers'
    plt.rcParams['axes.xmargin'] = 0
    plt.rcParams['axes.ymargin'] = 0
    plt.rcParams['font.size'] = 13

    plt.figure()
    plt.axhspan((wbar_exp[0]-wbar_exp[1]), (wbar_exp[0]+wbar_exp[1]), facecolor='#DCDCDC', zorder=1)
    plt.fill_between(df_sim.lambda_p, df_sim.wbar, color='#7F7F7F', alpha=1, zorder=2)
    plt.fill_between([1.0005, (lambda_avg-lambda_std_low)], [(wbar_exp[0]-wbar_exp[1]), (wbar_exp[0]-wbar_exp[1])], color='white', alpha=1, zorder=3)
    plt.fill_between([(lambda_avg+lambda_std_high), 1.16], [(wbar_exp[0]+wbar_exp[1]), 0.35], color='white', alpha=1, zorder=4)
    
    plt.plot(df_sim.lambda_p, df_sim.wbar, linestyle='-',color='k', linewidth=2, zorder=5)
    plt.plot([1.00, lambda_avg], [wbar_exp[0], wbar_exp[0]], linestyle='--', color='k', linewidth=1, zorder=5)
    plt.plot([lambda_avg,lambda_avg], [0, wbar_exp[0]], linestyle='--', color='k', linewidth=1, zorder=5)
    
    plt.xlabel(r'$\lambda_p$')
    plt.ylabel(r'$\overline{w}_S$')
    plt.savefig("figure-" + title, dpi=400)

    return lambda_avg, lambda_std

def run(sims):
    
    # aggregate results for each group
    experiments = []
    wbar_avgs = []
    wbar_stds = []
    lambda_avgs = []
    lambda_stds = []
    for sim in sims:
        xticks = np.arange(1.000, 1.180, step=0.020)
        yticks = np.arange(0.000, 0.400, step=0.050)
        [lambda_avg, lambda_std] = stretchVsIncisionOpening(sim.wbar_exp, sim.path, sim.odb, sim.title, xticks, yticks)
        print('Simulation:', sim.title)
        print('Estimated stretch:', lambda_avg,'+/-',lambda_std)
        experiments.append(sim.title)
        wbar_avgs.append(sim.wbar_exp[0])
        wbar_stds.append(sim.wbar_exp[1])
        lambda_avgs.append(lambda_avg)
        lambda_stds.append(lambda_std)
    
    # write experimental stretch information to file
    outFile = open('results-estimated_stretches.csv', 'w+')
    header = ['Experiment','w/l (exp, ave)','w/l (exp, std)','stretch (ave)','stretch (std)']
    output = zip(experiments, wbar_avgs, wbar_stds, lambda_avgs, lambda_stds)
    writer = csv.writer(outFile)
    writer.writerow(header)
    writer.writerows(output)
    outFile.close()


class ExperimentalGroup(object):
    def __init__(self, wbar_exp, path, odb, title):
        self.wbar_exp = wbar_exp
        self.path = path
        self.odb = odb
        self.title = title

neonateTran = ExperimentalGroup(
        [0.171, 0.073],
        '../simulations/results/',
        'neonate_tran_sim.csv',
        'neonate_transverse'
    )

neonateLong = ExperimentalGroup(
        [0.106, 0.045],
        '../simulations/results/',
        'neonate_long_sim.csv',
        'neonate_longitudinal'
    )

adultTran = ExperimentalGroup(
        [0.244, 0.072],
        '../simulations/results/',
        'adult_tran_sim.csv',
        'adult_transverse'
    )

adultLong = ExperimentalGroup(
        [0.094, 0.043],
        '../simulations/results/',
        'adult_long_sim.csv',
        'adult_longitudinal'
    )

# class neonateTran():
#     def __init__(self):
#         self.wbar_exp = [0.171, 0.073]
#         self.path = '../simulations/results/'
#         self.odb = 'neonate_tran_sim.csv'
#         self.title = 'neonate_transverse'
# class neonateLong():
#     def __init__(self):
#         self.wbar_exp = [0.106, 0.045]
#         self.path = '../simulations/results/'
#         self.odb = 'neonate_long_sim.csv'
#         self.title = 'neonate_longitudinal'
# class adultTran():
#     def __init__(self):
#         self.wbar_exp = [0.244, 0.072]
#         self.path = '../simulations/results/'
#         self.odb = 'adult_tran_sim.csv'
#         self.title = 'adult_transverse'
# class adultLong():
#     def __init__(self):
#         self.wbar_exp = [0.094, 0.043]
#         self.path = '../simulations/results/'
#         self.odb = 'adult_long_sim.csv'
#         self.title = 'adult_longitudinal'

if __name__ == '__main__':
    sims = [neonateTran, neonateLong, adultTran, adultLong]
    run(sims)
