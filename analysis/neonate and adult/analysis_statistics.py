""" Last updated: 2023/05/31

Script for statistically analyzing the change in incision opening 
(w/l) as it relaxes across five minutes. Here we determine the p-value 
and effect size of the experimental measurements.

p-value (statistical significance)------------------------------------------------

In statistics, the p-value is the probability of obtaining results at least 
as extreme as the observed results of a statistical hypothesis test, assuming 
that the null hypothesis is correct. The p-value serves as an alternative to 
rejection points to provide the smallest level of significance at which the 
null hypothesis would be rejected. A smaller p-value means that there is stronger 
evidence in favor of the alternative hypothesis. Bascially, a high p value says 
that your results are mostlikely due to chance.

Sources: 
1) What is a p-value: https://www.investopedia.com/terms/p/p-value.asp 
2) How to calculate student t statistic for unequal variance (by hand): 
https://www.youtube.com/watch?v=N2w6fb6O_Lg&t=75s
3) Python script for unequal variance: 
https://stackoverflow.com/questions/22611446/perform-2-sample-t-test

effect size (practical significance)----------------------------------------------

Effect size tells you how meaningful the relationship between variables or the 
difference between groups is. It indicates the practical significance of a 
research outcome. A large effect size means that a research finding has practical 
significance, while a small effect size indicates limited practical applications.

Using Cohen's d, for comparing two groups, we determine effect size from: 
d = (xbar_1 - xbar_2)/s, where xbar_1 and x_bar_2 are the means of the two groups,
respectively, and s is the standard deviation. s can be a pooled distribution based
on both groups, the control st. dev. if your study has a control, and st. dev. of 
pre-test data (if you have that). In our case, s is the pooled, and is the square
root of the average squared standard deviations 
s = sqrt((s_1^2 + s_2^2)/2).

Sources:
1) What is effect size: https://www.scribbr.com/statistics/effect-size/
2) Determining s for Cohen's d formula: https://www.uv.es/~friasnav/EffectSizeBecker.pdf
3) Online calculator: https://www.socscistatistics.com/effectsize/default3.aspx


"""
import csv
import numpy as np
import pandas as pd
import statistics as st
from scipy.optimize import curve_fit
from scipy.stats import ttest_ind, ttest_ind_from_stats
from scipy.special import stdtr

def experimentalStretch(stretchFit, wOverl, ave):
    global stretch
    for j in stretchFit.index:
        if stretchFit[j] < stretchFit.iloc[-1]:
            if ave >= wOverl[j] and ave <= wOverl[j+1]:
                stretch = stretchFit[j]+(stretchFit[j+1]-stretchFit[j])*((ave-wOverl[j])/(wOverl[j+1]-wOverl[j]))

def polyRegFit(x, a, b, c, d, e, f):
    polyReg = a*(x**5)+b*(x**4)+c*(x**3)+d*(x**2)+e*(x)+f
    return polyReg

def estimatedStretches(path, odb):
    # extract stretch and w/l data
    df = pd.read_csv(path+odb, encoding='ISO-8859-1')
    dfExp = pd.read_csv(path+'statistics_experimental_data.csv', encoding='ISO-8859-1')

    # curve fit of adult transverse data for standard deviation higher bound of stretch
    if odb == 'adultICRMouseDuraMaterTranCutResults.csv':
        orginalDFLam = df.lambda_p
        orginalDFRatio = df.wOverlHalf
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
        df = pd.concat([df, df_add], ignore_index=True)

    # estimate pre-stretch
    global estSt
    global estStrain
    global expWoverL
    global comp
    estSt = []
    estStrain = []
    expWoverL = []
    comp = []
    for i in range(len(dfExp.wOverl_t5min)):
        if dfExp.age_group[i][0] == odb[0]:
            if dfExp.cut_orientation[i][0] == odb[24] or dfExp.cut_orientation[i][0] == odb[22]:
                comp.append(dfExp.age_group[i][0] + dfExp.cut_orientation[i][0])
                experimentalStretch(df.lambda_p, df.wOverlHalf, dfExp.wOverl_t5min[i])
                estSt.append(stretch)
                estStrain.append(stretch-1)
                expWoverL.append(dfExp.wOverl_t5min[i])

def statisticalAnalysis(path, odb1, odb2):
    odbs = [odb1, odb2]
    samplewOverl = []
    sampleestSt  = []
    sampleestStrain  = []
    compList = []
    for odb in odbs:
        estimatedStretches(path, odb)
        samplewOverl.append(expWoverL)
        sampleestSt.append(estSt)
        sampleestStrain.append(estStrain)
        compList.append(comp)
    tData  = []
    tStats = []
    tF     = []
    eS     = []
    compare = []
    # p-value
    for test in [samplewOverl, sampleestSt, sampleestStrain]:
        compare.append(compList[0][0]+compList[1][0])
        # statistics w/ ttest_ind
        t, p = ttest_ind(np.array(test[0]), np.array(test[1]), equal_var=False)
        tData.append([t,p])
        # statistics w/ formula
        mean = []
        var  = []
        n    = []
        dofs = []
        stdv = []
        for metric in test:
            x = np.array(metric)
            xbar = x.mean()
            xvar = x.var(ddof=1)
            xn = x.size
            mean.append(xbar)
            var.append(xvar)
            n.append(xn)
            dofs.append(xn-1)
        # note: I wanted to calculate stdev.S just to confirm that the other tests using variance are correct
            nums = 0
            for i in range(xn):
                num = (x[i] - xbar)**2
                nums = nums + num
            std = np.sqrt(nums/(xn-1))
            stdv.append(std)
        t   = (mean[0] - mean[1]) / np.sqrt(var[0]/n[0] + var[1]/n[1])
        dof = (var[0]/n[0]+ var[1]/n[1])**2 / (var[0]**2/(n[0]**2*dofs[0]) + var[1]**2/(n[1]**2*dofs[1]))
        p = 2*stdtr(dof, -np.abs(t))
        tF.append([t,p])
        # statistics w/ ttest_ind_from_stats
        t, p = ttest_ind_from_stats(mean[0], stdv[0], n[0], mean[1], stdv[1], n[1], equal_var=False)
        tStats.append([t,p])
    # effect size
        s = np.sqrt((stdv[0]**2 + stdv[1]**2)/2)
        d = (mean[0]-mean[1])/s
        eS.append(d)
    # compare results
    print('w/l (ttest_ind):                t = %g p = %g' % (tData[0][0], tData[0][1]))
    print('w/l (ttest_ind_from_stats):     t = %g p = %g' % (tStats[0][0], tStats[0][1]))
    print('w/l (formula):                  t = %g p = %g' % (tF[0][0], tF[0][1]))
    print('w/l, effect size:               d = %g' % eS[0])
    print('stretch (ttest_ind):            t = %g p = %g' % (tData[1][0], tData[1][1]))
    print('stretch (ttest_ind_from_stats): t = %g p = %g' % (tStats[1][0], tStats[1][1]))
    print('stretch (formula):              t = %g p = %g' % (tF[1][0], tF[1][1]))
    print('stretch, effect size:           d = %g' % eS[1])
    print('strain (ttest_ind):             t = %g p = %g' % (tData[2][0], tData[2][1]))
    print('strain (ttest_ind_from_stats):  t = %g p = %g' % (tStats[2][0], tStats[2][1]))
    print('strain (formula):               t = %g p = %g' % (tF[2][0], tF[2][1]))
    print('strain, effect size:            d = %g' % eS[2])
    # write to .csv file
    output = zip(compare, ['w/l', 'stretch', 'strain'], [tF[0][0], tF[1][0], tF[2][0]], [tF[0][1], tF[1][1], tF[2][1]], eS)
    writer.writerows(output)

# input .csv files
path = '../../simulations/neonate and adult .csv files/'
odb1 = 'neonateICRMouseDuraMaterTranCutResults.csv'
odb2 = 'neonateICRMouseDuraMaterLongCutResults.csv'
odb3 = 'adultICRMouseDuraMaterTranCutResults.csv'
odb4 = 'adultICRMouseDuraMaterLongCutResults.csv'

# analysis
outFile = open('results-p-value_and_effect_size.csv', 'w+')
header = ['comparison', 'metric', 'test statistic', 'p-value', 'effect size']
writer = csv.writer(outFile)
writer.writerow(header)
print('------------------------------------------------------------------')
print('neonate transverse vs neonate longitudinal')
statisticalAnalysis(path, odb1, odb2)
print('------------------------------------------------------------------')
print('adult transverse vs adult longitudinal')
statisticalAnalysis(path, odb3, odb4)
print('------------------------------------------------------------------')
print('neonate transverse vs adult transverse')
statisticalAnalysis(path, odb3, odb1)
print('------------------------------------------------------------------')
print('neonate longitudinal vs adult longitudinal')
statisticalAnalysis(path, odb2, odb4)
print('------------------------------------------------------------------')
outFile.close()
