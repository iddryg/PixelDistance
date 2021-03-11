# This script analyzes the csv files output by main / pixel_distance.py
# uses some methods from pixel_distance.py

# pixel_distance.py actually performs the measurement of minimum distance
# between tumor and lyve-1 pixels, and outputs the results for each image.
# PixDistStats.py performs stats and makes plots on ALL the data separated by sample group. However,
# this is insufficient because it isn't split up into biological replicates, or normalized.
# PixDistStats2.py separates the data into biological replicates instead of aggregating
# all data for each sample group, and experiments with plots.
# PixDistStats3.py takes data from PixDistStats2, normalizes it to total pixels for each animal,
# does statistical comparisons and makes plots.

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import pixel_distance as pxd
import pandas as pd
from scipy.stats import stats
from statsmodels.stats.multicomp import pairwise_tukeyhsd, MultiComparison



def analyze(data, save_dir, labels):
    save_name = 'Results'

    # Scale the data
    # 0.754 microns per pixel
    print('data before scaling: ')
    print(data.head(15))
    data = 0.754*data
    print('data after scaling: ')
    print(data.head(15))
    data.replace(0, np.nan, inplace=True)
    print('adat after removing zeros: ')
    print(data.head(15))

    # Run One-way ANOVA
    # F is the F-statistic (an array, one for each group)
    # p is the p-value (an array, one for each group)
    # f_stat, p = stats.f_oneway(data[labels[0]][~np.isnan(data[labels[0]])],
    #                 data[labels[1]][~np.isnan(data[labels[1]])],
    #                 data[labels[2]][~np.isnan(data[labels[2]])])

    f_stat, p = stats.f_oneway(data[labels[0]].dropna(),
                               data[labels[1]].dropna(),
                               data[labels[2]].dropna())

    # Make new dataframe to stack the data for tukey comparisons...
    # print(data.head())
    data_stacked = data.stack().reset_index()
    data_stacked = data_stacked.dropna().rename(columns = {'level_0': 'id', 'level_1': 'group', 0: 'distance'})
    print(data_stacked.head(20))

    # Multiple comparisons... Tukey Test:
    mc = MultiComparison(data_stacked['distance'], data_stacked['group'])
    tukey = mc.tukeyhsd(alpha=0.05)
    print('Tukey results: ')
    print(tukey)
    print('Unique groups: {}'.format(mc.groupsunique))

    # Save ANOVA & Tukey results in a text file
    file0 = open(save_dir + 'ANOVA_Results.txt', 'a+')
    file0.write('ANOVA Results: \n')
    file0.write('F Statistic: ' + str(f_stat) + '\n')
    file0.write('p-value: ' + str(p) + '\n')
    file0.write('Tukey results: ' + '\n')
    file0.write(str(tukey) + '\n')
    file0.write('Unique groups: {}'.format(mc.groupsunique))

    # print('F and its shape and type: ')
    # print(f_stat)
    # print(np.shape(f_stat))
    # print(type(f_stat))
    # print('p and its shape and type: ')
    # print(p)
    # print(np.shape(p))
    # print(type(p))

    sns.set_theme(style="whitegrid")
    # ax = sns.boxplot(y=min_dists)
    # ax = sns.swarmplot(y=min_dists, color=".25")

    # violin & strip plot
    ax = sns.violinplot(data=data.dropna(), inner='quartile')
    # sns.stripplot(y=min_dists)

    # Update Feb 2021: Haley doesn't want text on the plots. ---------------------------------------
    # # Prep x and y position for text drawing (for p-values and sample sizes)
    # ypos = data.max()
    # xpos_list = [0, 1, 2]  # 0 1 2
    # print('y: ')
    # print(np.shape(ypos))
    # print(ypos)
    # print('x: ')
    # print(np.shape(xpos_list))
    # print(xpos_list)
    #
    # # Draw sample sizes and p-values onto the violin plot
    # for i in range(len(xpos_list)):
    #     ax.text(xpos_list[i], ypos[1], 'Tumor Pixels: \n' + str(data[labels[i]].dropna().count()) + '\n p = ' + str(p))
    # -----------------------------------------------------------------------------------------------

    # Save violin plot
    fig = ax.get_figure()
    fig.savefig(save_dir + save_name + 'violin.png')
    fig.clf()
    # plt.savefig(save_dir + save_name + 'violin.png')
    # plt.clf()

    # For loop for Histograms and txt file summaries for each group
    counter = 0
    for label in labels:
        # histogram plot
        sns.histplot(data=data[label].dropna())
        # ax.set(ylabel='Pixels')
        plt.title('Pixel Distance between tumor and lymphatic vessels for ' + label)
        plt.xlabel('Distance (pixels)')
        # plt.show()
        plt.savefig(save_dir + label + 'histogram.png')
        plt.clf()

        # Save some basic stats in a text file
        file1 = open(save_dir + label + '.txt', 'a+')
        file1.write('Stats for the sample: ' + label + '\n')
        file1.write('Number of Pixels used: ' + str(len(data[label].dropna())) + '\n')
        file1.write('Mean: ' + str(np.mean(data[label].dropna())) + '\n')
        file1.write('Standard Deviation: ' + str(np.std(data[label].dropna())) + '\n')
        file1.write('Percentiles: \n')
        file1.write('10: ' + str(np.percentile(data[label].dropna(), 10)) + '\n')
        file1.write('25: ' + str(np.percentile(data[label].dropna(), 25)) + '\n')
        file1.write('50: ' + str(np.percentile(data[label].dropna(), 50)) + '\n')
        file1.write('75: ' + str(np.percentile(data[label].dropna(), 75)) + '\n')
        file1.write('90: ' + str(np.percentile(data[label].dropna(), 90)) + '\n')
        counter += 1


    print('Violin plot, histogram plot, stats.txt and data file saved.')
    print('Part 6 - Plot Histogram of Distances: Complete')
    return

def process_csv_files(maindir):
    # initialize arrays...
    naive_array = np.zeros(0)
    naive_size_array = np.zeros(45)
    naive_size_count = 0
    ndln_array = np.zeros(0)
    ndln_size_array = np.zeros(45)
    ndln_size_count = 0
    tdln_array = np.zeros(0)
    tdln_size_array = np.zeros(45)
    tdln_size_count = 0

    print('start1 ------------------------------------')
    for root, dirs, files in os.walk(maindir):
        for name in files:
            if 'csv' in name and 'naive' in name:
                # flattened append
                naive_array = np.append(naive_array, np.genfromtxt(os.path.join(root, name), delimiter=','))
                print('naive array shape:')
                print(np.shape(naive_array))
                naive_size_array[naive_size_count] = len(np.genfromtxt(os.path.join(root, name), delimiter=','))
                naive_size_count += 1
                print('naive count: ' + str(naive_size_count))
            if 'csv' in name and 'ndLN' in name:
                ndln_array = np.append(ndln_array, np.genfromtxt(os.path.join(root, name), delimiter=','))
                print('ndLN array shape:')
                print(np.shape(ndln_array))
                ndln_size_array[ndln_size_count] = len(np.genfromtxt(os.path.join(root, name), delimiter=','))
                ndln_size_count += 1
                print('ndLN count: ' + str(ndln_size_count))
            if 'csv' in name and 'tdLN' in name:
                tdln_array = np.append(tdln_array, np.genfromtxt(os.path.join(root, name), delimiter=','))
                print('tdLN array shape:')
                print(np.shape(tdln_array))
                tdln_size_array[tdln_size_count] = len(np.genfromtxt(os.path.join(root, name), delimiter=','))
                tdln_size_count += 1
                print('tdLN count: ' + str(tdln_size_count))
    print(tdln_size_array)

    print('data collection complete.')
    print('beginning analysis.')

    # Create a data dictionary to load arrays into a pandas dataframe
    data = {
        'naive' : naive_array,
        'ndln' : ndln_array,
        'tdln' : tdln_array
    }
    # Want all arrays to be same size for pandas
    maxsize = max([a.size for a in data.values()])
    # pad the shorter arrays with NaNs to make them the same length as the longest array
    data_pad = {k: np.pad(v, pad_width=(0, maxsize - v.size,), mode='constant', constant_values=np.nan) for k, v in
                data.items()}
    df = pd.DataFrame(data_pad)
    analyze(df, maindir, ['naive', 'ndln', 'tdln'])

    # analyze(naive_array, naive_size_array, maindir, 'naive')
    # analyze(ndln_array, ndln_size_array, maindir, 'ndln')
    # analyze(tdln_array, tdln_size_array, maindir, 'tdln')


    return


dirname = pxd.file_import(prompt='Choose the directory containing tiff folders: ')
print(dirname)

process_csv_files(dirname)

