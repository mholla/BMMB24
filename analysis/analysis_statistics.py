""" Last updated: 2023/05/31

Statistics for comparing measured cut openings and estimated prestretch 
between age groups and cut orientations 

"""

import csv
import numpy as np
import pandas as pd
import statistics as st
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import ttest_ind

import analysis_functions

def estimatedStretches(path, csv):
    """ estimate stretch and strain corresponding to each experimental wbar value """
    # extract stretch and w/l data 
    df_sim = pd.read_csv(path+csv, encoding='ISO-8859-1')
    df_exp = pd.read_csv('../experimental data/experimental_data.csv', encoding='ISO-8859-1')

    # curve fit of adult transverse data for standard deviation higher bound of stretch
    if csv == 'adult_tran_sim.csv':
        df_sim = analysis_functions.extrapolate(df_sim)
    
    group_name = []
    wbar_exp = []
    lambda_est = []
    strain_est = []

    for i in range(len(df_exp.wbar)):
        age_group = csv.split("_")[0]
        cut_orientation = csv.split("_")[1]

        if df_exp.age_group[i] == age_group and df_exp.cut_orientation[i] == cut_orientation:
                # extract experimental data for relevant group
                group_name.append(df_exp.age_group[i][0] + df_exp.cut_orientation[i][0])
                wbar_exp.append(df_exp.wbar_5min[i])
                
                # estimate corresponding stretch
                stretch = analysis_functions.interpolateStretch(df_sim.lambda_p, df_sim.wbar, df_exp.wbar_5min[i])
                
                # append stretch and strain values
                lambda_est.append(stretch)
                strain_est.append(stretch-1)

    return group_name, wbar_exp, lambda_est, strain_est

def statisticalAnalysis(path, csv1, csv2):
    csvs = [csv1, csv2]
    
    # gather all data for both comparison groups
    wbar_sets = []
    lambda_sets  = []
    strain_sets  = []
    group_names = []
    
    for csv in csvs:
        [group_name, wbar_exp, lambda_est, strain_est] = estimatedStretches(path, csv)
        wbar_sets.append(wbar_exp)
        lambda_sets.append(lambda_est)
        strain_sets.append(strain_est)
        group_names.append(group_name)
    
    # perform statistics on both quantities (measured wbar and estimated stretch)
    p_values  = []
    effect_size = []
    comparison = []
    
    for test in [wbar_sets, lambda_sets]:
        comparison.append(group_names[0][0]+group_names[1][0])
        
        # statistics w/ ttest_ind
        [t, p] = ttest_ind(np.array(test[0]), np.array(test[1]), equal_var=False)
        p_values.append(p)

        # effect size
        mean = [st.mean(test[0]), st.mean(test[1])]
        stdv = [np.std(test[0], ddof=1), np.std(test[1], ddof=1)]
        s = np.sqrt((stdv[0]**2 + stdv[1]**2)/2)
        effect_size.append((mean[0]-mean[1])/s)
    
    # print results
    print('w/l:                p = %g, d = %g' % (p_values[0], effect_size[0]))
    print('stretch:            p = %g, d = %g' % (p_values[1], effect_size[1]))

    # write to results file file
    output = zip(comparison, ['w/l', 'stretch'], p_values, effect_size)
    writer.writerows(output)

def plot(path, csvs):

    wbar_sets = []
    lambda_sets  = []
    strain_sets  = []
    group_names = []

    wbar_mean = []
    lambda_mean = []
    wbar_stdv = []
    lambda_stdv = []

    for csv in csvs:
        [group_name, wbar_exp, lambda_est, strain_est] = estimatedStretches(path, csv)
        wbar_sets.append(wbar_exp)
        lambda_sets.append(lambda_est)
        strain_sets.append(strain_est)
        group_names.append(group_name[0])

        wbar_mean.append(st.mean(wbar_exp))
        lambda_mean.append(st.mean(lambda_est))
        wbar_stdv.append(np.std(wbar_exp, ddof=1))
        lambda_stdv.append(np.std(lambda_est, ddof=1))

    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = "serif"
    plt.rcParams['axes.autolimit_mode'] = 'round_numbers'
    plt.rcParams['axes.xmargin'] = 0
    plt.rcParams['axes.ymargin'] = 0
    plt.rcParams['font.size'] = 13
    plt.rcParams['hatch.linewidth'] = 2

    cmap = plt.cm.get_cmap("Blues")
    neonate = cmap(0.2)
    adult = cmap(1.0)

    colors = [neonate, neonate, adult, adult]
    hatches = ['--', '||', '--', '||']

    plt.figure()
    bars = plt.bar(group_names, wbar_mean, yerr=wbar_stdv, capsize=4, color=colors)
    for i in range(len(group_names)):
        bars[i].set(hatch = hatches[i], edgecolor='white')
    plt.ylabel(r'$\overline{w}_E$')
    plt.savefig("figure-bar_wbar", dpi=400)

    plt.figure()
    bars = plt.bar(group_names, lambda_mean, yerr=lambda_stdv, capsize=4, color=colors)
    for i in range(len(group_names)):
        bars[i].set(hatch = hatches[i], edgecolor='white')
    plt.ylim([1, 1.2])
    plt.ylabel(r'$\lambda_p$')
    plt.savefig("figure-bar_lambda", dpi=400)



if __name__ == '__main__':

    # input .csv files
    path = '../simulations/results/'
    csv1 = 'neonate_tran_sim.csv'
    csv2 = 'neonate_long_sim.csv'
    csv3 = 'adult_tran_sim.csv'
    csv4 = 'adult_long_sim.csv'

    plot(path, [csv1, csv2, csv3, csv4])

    # analysis
    outFile = open('results-statistics.csv', 'w+')
    header = ['comparison', 'metric', 'test statistic', 'p-value', 'effect size']
    writer = csv.writer(outFile)
    writer.writerow(header)
    print('------------------------------------------------------------------')
    print('neonate transverse vs neonate longitudinal')
    statisticalAnalysis(path, csv1, csv2)
    print('------------------------------------------------------------------')
    print('adult transverse vs adult longitudinal')
    statisticalAnalysis(path, csv3, csv4)
    print('------------------------------------------------------------------')
    print('neonate transverse vs adult transverse')
    statisticalAnalysis(path, csv3, csv1)
    print('------------------------------------------------------------------')
    print('neonate longitudinal vs adult longitudinal')
    statisticalAnalysis(path, csv2, csv4)
    print('------------------------------------------------------------------')
    outFile.close()

    plt.show()