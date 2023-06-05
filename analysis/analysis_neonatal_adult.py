""" Last updated: 2023/06/05

Script for plotting the change in incision opening (w/l) as stretch increases,
for the transverse and longitudinal cuts in the neonatal and adult murine 
cranial dura mater simulations.

"""
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

def exp_and_est_data(path, csvs):
    group_names = []
    wbar_sets = []
    lambda_sets  = []
    strain_sets  = []
    
    for csv in csvs:
        [group_name, wbar_exp, lambda_est, strain_est] = estimate_individual_stretches(path, csv)
        group_names.append(group_name[0])
        wbar_sets.append(wbar_exp)
        lambda_sets.append(lambda_est)
        strain_sets.append(strain_est)

    return group_names, wbar_sets, lambda_sets, strain_sets

def estimate_individual_stretches(path, csv):
    """ estimate stretch and strain corresponding to each experimental wbar value """
    # extract stretch and w/l data 
    df_sim = pd.read_csv(path+csv, encoding='ISO-8859-1')
    df_exp = pd.read_csv('../experiments/experimental_data.csv', encoding='ISO-8859-1')
    
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
                wbar_exp.append(df_exp.wbar[i])
                
                # estimate corresponding stretch
                stretch = interpolateStretch(df_sim.lambda_p, df_sim.wbar, df_exp.wbar[i])
                
                # append stretch and strain values
                lambda_est.append(stretch)
                strain_est.append(stretch-1)

    return group_name, wbar_exp, lambda_est, strain_est

def exp_and_est_means(csvs, wbar_sets, lambda_sets):
    wbar_mean = []
    lambda_mean = []
    wbar_stdv = []
    lambda_stdv = []

    for i in range(len(csvs)):
        wbar_mean.append(st.mean(wbar_sets[i]))
        wbar_stdv.append(np.std(wbar_sets[i], ddof=1))
        lambda_mean.append(st.mean(lambda_sets[i]))
        lambda_stdv.append(np.std(lambda_sets[i], ddof=1))

    return wbar_mean, wbar_stdv, lambda_mean, lambda_stdv

def plot(groups):
    
    [group_names, wbar_sets, lambda_sets, strain_sets] = exp_and_est_data(path, csvs)
    [wbar_mean, wbar_stdv, lambda_mean, lambda_stdv] = exp_and_est_means(csvs, wbar_sets, lambda_sets)
    
    for i in range(len(groups)):
        print('Simulation:', groups[i].title)
        print('Experimental wbar:', wbar_mean[i],'+/-',wbar_stdv[i])
        print('Estimated stretch:', lambda_mean[i],'+/-',lambda_stdv[i])

        # extract stretch and w/l data from simulations
        df_sim = pd.read_csv(path+csvs[i], encoding='ISO-8859-1')

        # plot
        plt.rcParams['text.usetex'] = True
        plt.rcParams['font.family'] = "serif"
        plt.rcParams['axes.autolimit_mode'] = 'round_numbers'
        plt.rcParams['axes.xmargin'] = 0
        plt.rcParams['axes.ymargin'] = 0
        plt.rcParams['font.size'] = 13

        plt.figure()
        # horizontal
        plt.axhspan((wbar_mean[i]-wbar_stdv[i]), (wbar_mean[i]+wbar_stdv[i]), facecolor='#DCDCDC', zorder=1)
        # fill in under curve
        plt.fill_between(df_sim.lambda_p, df_sim.wbar, color='#7F7F7F', alpha=1, zorder=2)
        # white out left
        plt.fill_between([1.000, (lambda_mean[i]-lambda_stdv[i])], [(wbar_mean[i]-wbar_stdv[i]), (wbar_mean[i]-wbar_stdv[i])], color='white', alpha=1, zorder=3)
        # white out right side
        plt.fill_between([(lambda_mean[i]+lambda_stdv[i]), 1.16], [(wbar_mean[i]+wbar_stdv[i]), 0.35], color='white', alpha=1, zorder=4)
    
        plt.plot(df_sim.lambda_p, df_sim.wbar, linestyle='-',color='k', linewidth=2, zorder=5)
        plt.plot([1.00, lambda_mean[i]], [wbar_mean[i], wbar_mean[i]], linestyle='--', color='k', linewidth=1, zorder=5)
        plt.plot([lambda_mean[i],lambda_mean[i]], [0, wbar_mean[i]], linestyle='--', color='k', linewidth=1, zorder=5)
    
        plt.xlabel(r'$\lambda_p$')
        plt.ylabel(r'$\overline{w}_S$')
        ax1 = plt.twiny()
        ax1.axes.get_xaxis().set_visible(False)
        ax1.spines['bottom'].set_linewidth(1)
        plt.savefig("figure-" + groups[i].title, dpi=400)
    
    # write experimental stretch information to file
    outFile = open('results-estimated_stretches.csv', 'w+')
    header = ['Experiment','w/l (exp, ave)','w/l (exp, std)','stretch (ave)','stretch (std)']
    output = zip(group_names, wbar_mean, wbar_stdv, lambda_mean, lambda_stdv)
    writer = csv.writer(outFile)
    writer.writerow(header)
    writer.writerows(output)
    outFile.close()

def statistical_analysis(path, csv1, csv2):
    csvs = [csv1, csv2]
    
    # gather all data for both comparison groups
    [group_names, wbar_sets, lambda_sets, strain_sets] = exp_and_est_data(path, csvs)
    
    # perform statistics on both quantities (measured wbar and estimated stretch)
    p_values  = []
    effect_size = []
    comparison = []
    
    for test in [wbar_sets, lambda_sets]:
        comparison.append(group_names[0]+group_names[1])
        
        # statistics w/ ttest_ind
        [t, p] = ttest_ind(np.array(test[0]), np.array(test[1]), equal_var=False)
        p_values.append(p)

        # effect size
        mean = [st.mean(test[0]), st.mean(test[1])]
        stdv = [np.std(test[0], ddof=1), np.std(test[1], ddof=1)]
        s = np.sqrt((stdv[0]**2 + stdv[1]**2)/2)
        effect_size.append((mean[0]-mean[1])/s)
    
    # print results
    print('------------------------------------------------------------------')
    print(comparison)
    print('w/l:                p = %g, d = %g' % (p_values[0], effect_size[0]))
    print('stretch:            p = %g, d = %g' % (p_values[1], effect_size[1]))

    # write to results file file
    output = zip(comparison, ['w/l', 'stretch'], p_values, effect_size)
    writer.writerows(output)

def bar_plot(path, csvs):

    [group_names, wbar_sets, lambda_sets, strain_sets] = exp_and_est_data(path, csvs)

    [wbar_mean, wbar_stdv, lambda_mean, lambda_stdv] = exp_and_est_means(csvs, wbar_sets, lambda_sets)

    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = "serif"
    plt.rcParams['axes.autolimit_mode'] = 'round_numbers'
    plt.rcParams['axes.xmargin'] = 0
    plt.rcParams['axes.ymargin'] = 0
    plt.rcParams['font.size'] = 13
    plt.rcParams['hatch.linewidth'] = 2
    plt.rcParams["figure.autolayout"] = True

    cmap = plt.cm.get_cmap("Blues")
    neonate = cmap(0.2)
    adult = cmap(1.0)

    group_names = ['nT', 'nL', 'aT', 'aL']
    colors = [neonate, neonate, adult, adult]
    hatches = ['--', '||', '--', '||']

    plt.figure(figsize=[4,5])
    bars = plt.bar(group_names, wbar_mean, yerr=wbar_stdv, capsize=4, color=colors)
    for i in range(len(group_names)):
        bars[i].set(hatch = hatches[i], edgecolor='white')
    plt.ylabel(r'$\overline{w}_E$')
    plt.savefig("figure-bar_wbar", dpi=400)

    plt.figure(figsize=[4,5])
    bars = plt.bar(group_names, lambda_mean, yerr=lambda_stdv, capsize=4, color=colors)
    for i in range(len(group_names)):
        bars[i].set(hatch = hatches[i], edgecolor='white')
    plt.ylim([1, 1.2])
    plt.ylabel(r'$\lambda_p$')
    plt.savefig("figure-bar_lambda", dpi=400)

class ExperimentalGroup(object):
    def __init__(self, path, csv, title):
        self.path = path
        self.csv = csv
        self.title = title

neonateTran = ExperimentalGroup(
        
        '../simulations/results/',
        'neonate_tran_sim.csv',
        'neonate_transverse'
    )

neonateLong = ExperimentalGroup(
        
        '../simulations/results/',
        'neonate_long_sim.csv',
        'neonate_longitudinal'
    )

adultTran = ExperimentalGroup(
        
        '../simulations/results/',
        'adult_tran_sim.csv',
        'adult_transverse'
    )

adultLong = ExperimentalGroup(
        
        '../simulations/results/',
        'adult_long_sim.csv',
        'adult_longitudinal'
    )

if __name__ == '__main__':
    
    path = '../simulations/results/'

    groups = [neonateTran, neonateLong, adultTran, adultLong]
    csvs = []
    for group in groups: 
        csvs.append(group.csv)

    plot(groups)

    bar_plot(path, csvs)

    # statistical analysis
    outFile = open('results-statistics.csv', 'w+')
    header = ['comparison', 'metric', 'p-value', 'effect size']
    writer = csv.writer(outFile)
    writer.writerow(header)
    statistical_analysis(path, csvs[0], csvs[1])
    statistical_analysis(path, csvs[2], csvs[3])
    statistical_analysis(path, csvs[0], csvs[2])
    statistical_analysis(path, csvs[1], csvs[3])
    outFile.close()

    plt.show()
