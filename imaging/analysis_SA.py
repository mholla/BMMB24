import numpy as np
import matplotlib.pyplot as plt


if __name__ == '__main__':

    # read data file
    data_file = 'data_SA.csv'

    with open('data_SA.csv','r') as data_file:
        data_lines = data_file.readlines()        
         
    weeks = [4,8,12,16,20]
    subjects = ["F1", "F2", "F3", "M1", "M2", "M3"]

    data_SA = np.zeros((len(subjects), len(weeks)))

    for i in range(len(subjects)): 
        data = data_lines[i+1].split(',')
        for j in range(len(weeks)): 
            data_SA[i][j] = float(data[j+1])

    # zeroing out F1 Week 8 for procedural errors
    data_SA[0][1] = None
    data_SA_masked = np.ma.masked_array(data_SA, np.isnan(data_SA))

    # focus on 4-12 weeks
    max_week = 3
    woi = weeks[0:max_week]
    
    # average data by week
    data_SA_transposed = np.transpose(data_SA_masked)

    weekly_average = np.zeros(max_week)
    weekly_stdev = np.zeros(max_week)

    for i in range(max_week): 
        weekly_average[i] = np.mean(data_SA_transposed[i])
        weekly_stdev[i] = np.std(data_SA_transposed[i])

    # plot individual data along with average and st. dev.
    
    plt.rcParams['text.usetex'] = True
    plt.rcParams['font.family'] = "serif"
    plt.rcParams['font.size'] = 13

    plt.figure()
    for i in range(len(subjects)):
        plt.plot(woi,data_SA[i][0:max_week], color="lightgray", linestyle=":")
    plt.errorbar(woi, weekly_average, yerr=weekly_stdev, fmt='k-', ecolor='k', capsize=5)

    plt.xticks(woi)
    plt.xlabel("Weeks")

    plt.gca().set_ylim(400,500)
    plt.ylabel(r"Surface area [$mm^2$]")

    plt.savefig("SA.png")
    plt.show()
