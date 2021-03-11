# This script analyzes the csv files output by PixDistStats2.py
# Updated Feb 2021.
# PixDistStats2 separates the data into biological replicates instead of aggregating all data for each sample group.
# This script takes those data and does stats and makes plots.

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
import joypy as jpy


def load_datas(dir):
    distbypercentiles = pd.read_csv(dir + 'dist_by_percentiles.csv', index_col='percentiles')
    numpixbydistbins = pd.read_csv(dir + 'numpix_by_dist_bins.csv', index_col='distance bins')
    normnumpixbydistbins = pd.read_csv(dir + 'norm_numpix_by_dist_bins.csv', index_col='distance bins')
    print('dist by percentiles: ')
    print(distbypercentiles.head(10))
    print('numpix by dist bins: ')
    print(numpixbydistbins.head(11))
    print('normalized numpix by dist bins: ')
    print(normnumpixbydistbins.head(11))

    # return datas as a list
    return [distbypercentiles, numpixbydistbins, normnumpixbydistbins]


def run_anova(data, savedir, labels):
    # ANOVA
    f_stat, p = stats.f_oneway(data[labels[0:4]],
                               data[labels[5:9]],
                               data[labels[10:14]])
    # Multiple comparisons... Tukey Test:
    # need a stacked dataframe with consistent labels...
    # measurement is the data, group is naive, tdLN, disLN
    # mc = MultiComparison(data_stacked['measurement'], data_stacked['group'])

    # stack data
    # data_stacked = data.stack().reset_index()
    # data is in Series form... so it's already stacked. Just reset_index()
    data_stacked = data.to_frame()
    # data_stacked = data_stacked.rename(columns={'level_0': 'id', 'level_1': 'group', 0: 'distance'})
    print(data_stacked.head(20))
    # make new column with supergroups (naive, disLN, tdLN)
    # data_stacked['supergroup'] = data_stacked['group'].map(lambda x: x.rstrip('12345'))
    data_stacked['supergroup'] = data_stacked.index.map(lambda x: x.rstrip('12345'))
    print(data_stacked.head(20))
    # mc = MultiComparison(data_stacked['distance'], data_stacked['supergroup'])
    mc = MultiComparison(data_stacked[data.name], data_stacked['supergroup'])
    tukey = mc.tukeyhsd(alpha=0.05)

    print(data_stacked[data.name])
    # Save ANOVA & Tukey results in a text file
    file0 = open(savedir + data.name + '_ANOVA.txt', 'a+')
    file0.write('Stats: \n')
    file0.write('Mean: ' + str(data_stacked.groupby(['supergroup']).mean()) + '\n')
    file0.write('Standard Deviation: ' + str(data_stacked.groupby(['supergroup']).std()) + '\n')
    file0.write('ANOVA Results: \n')
    file0.write('F Statistic: ' + str(f_stat) + '\n')
    file0.write('p-value: ' + str(p) + '\n')
    file0.write('Tukey results: ' + '\n')
    file0.write(str(tukey) + '\n')
    file0.write('Unique groups: {}'.format(mc.groupsunique))

    return


def transpose_data(data):
    # remove name of indexes
    data.index.names = [None]
    transposed_data = data.transpose()
    print('after transposing: ')
    print(transposed_data)
    # drop the number from the end of the indexes (make supergroups)
    transposed_data = transposed_data.rename(index=lambda x: x.rstrip('12345'))
    print('after renaming indexes: ')
    print(transposed_data)
    # stack based on supergroup
    # transposed_stacked_data = transposed_data.stack()
    # print('after stacking: ')
    # print(transposed_stacked_data)
    return transposed_data


def make_plots(dist_percentile, norm_numpix, savedir, labels):
    sns.set_theme(style="whitegrid")
    # keep x-axis as distance consistent across plots.
    # For dist_by_percentiles_transposed, try plotting a bar graph which shows the distance tumor cells invaded to
    # at each percentile
    # 10th %ile |---|               10% of cells invaded less than this distance
    # 20th %ile |-------|           20% of cells invaded less than this distance
    # 30th %ile |------------|      30% of cells invaded less than this distance
    # For norm_numpix_by_dist_bins_transposed, try plotting a histogram...
    # Proportion (normalized #) of pixels at each distance.
    # Can overlay all three histograms in different colors, slightly opaque
    # ------------------------------------------------------------------------------------------
    # bar plots for dist_percentile:
    print('initial assessment: ')
    dist_percentile.index.names = ['Group']
    print(dist_percentile)
    print(dist_percentile.index)
    # convert indexes to a column so we can melt it
    dist_percentile.reset_index(level=dist_percentile.index.names, inplace=True)
    print('after reset index: ')
    print(dist_percentile)
    melt_dist_percentile = pd.melt(dist_percentile, id_vars='Group', var_name='Percentile',
                               value_name='Distance (microns)')
    ax2 = sns.barplot(x='Distance (microns)', y='Percentile', hue='Group', data=melt_dist_percentile)
    fig2 = ax2.get_figure()
    fig2.set_size_inches(11, 8.5)  # increase figure size
    plt.gca().legend().set_title(None)  # remove legend title
    plt.gca().set_title('Distance from Lymphatics by Percentile')  # set plot title
    # Add annotations for statistical significance based on earlier anova & tukey comparisons (see txt files)
    # which comparisons were significant? by tukey:
    # 30th: disLN & tdLN. p-adj = 0.0401
    # 40th: disLN & tdLN. p-adj = 0.0191
    # 50th: disLN & tdLN. p-adj = 0.0126, naive & tdLN. p-adj = 0.0369
    # 60th: disLN & tdLN. p-adj = 0.012, naive & tdLN. p-adj = 0.0177
    # 70th: disLN & tdLN. p-adj = 0.0153, naive & tdLN. p-adj = 0.0122
    # 80th: disLN & tdLN. p-adj = 0.0221, naive & tdLN. p-adj = 0.011
    fig2.savefig(savedir + 'dist_by_percentiles.png')
    fig2.clf()
    # -----------------------------------------------------------------------------------------------------
    # histograms for norm_numpix:
    # this isn't actually a histogram... since I already have the x-labels as bins and
    # the counts (proportions) for each sample. What I really want to do is create a bunch of bar plots.
    # fig, ax = plt.subplots()
    # for a in [x, y]:
    #     sns.distplot(a, bins=range(1, 110, 10), ax=ax, kde=False)
    # ax.set_xlim([0, 100])
    # Try melting...
    print('before index rename attempt: ')
    print(norm_numpix.index)
    norm_numpix.index.names = ['Group']
    print('after index rename attempt: ')
    print(norm_numpix)
    print(norm_numpix.index)
    # convert indexes to a column so we can melt it
    norm_numpix.reset_index(level=norm_numpix.index.names, inplace=True)
    print('after reset index: ')
    print(norm_numpix)
    melt_norm_numpix = pd.melt(norm_numpix, id_vars='Group', var_name='Distance (microns)',
                               value_name='% of total pixels within group')
    print('after melting: ')
    print(melt_norm_numpix.head())
    # # Stack Data
    # norm_numpix = norm_numpix.stack()
    # print('after stacking: ')
    # print(norm_numpix)
    # print('indexes: ')
    # print(norm_numpix.index)
    # # samples = ['tdLN', 'disLN', 'naive']
    # # dist_bins = ['0-10um', '10-20um', '20-30um', '30-40um', '40-50um',
    # #              '50-60um', '60-70um', '70-80um', '80-90um', '90-100um', '100um+']
    # # norm_numpix.index = pd.MultiIndex.from_product([samples, dist_bins], names=['sample', 'dist_bin'])
    # # norm_numpix.rename_axis(index=['sample', 'dist_bin'])
    # norm_numpix.index.names = ['sample', 'dist_bin']
    # print('after rename attempt: ')
    # print(norm_numpix)
    # print(norm_numpix.index)
    # # g = sns.FacetGrid(norm_numpix, hue='sample', palette='coolwarm')
    ax = sns.barplot(x='Distance (microns)', y='% of total pixels within group', hue='Group', data=melt_norm_numpix)
    fig = ax.get_figure()
    fig.set_size_inches(11, 8.5)  # increase figure size
    plt.gca().legend().set_title(None)  # remove legend title
    plt.gca().set_title('% of Tumor+ Pixels vs. Distance from Lymphatics')  # set plot title
    # Add annotations for statistical significance based on earlier anova & tukey comparisons (see txt files)
    # which comparisons were significant? by tukey:
    # in general... 0-20um: tdLN sig lower. 30-50um: tdLN sig higher.
    # 0-10um: disLN & tdLN. p-adj = 0.0472
    # 10-20um: naive & tdLN. p-adj = 0.0306
    # 30-40um: naive & tdLN. p-adj = 0.0014
    # 40-50um: disLN & tdLN. p-adj = 0.0019. naive & tdLN. p-adj = 0.001

    fig.savefig(savedir + 'numpix_by_dist_bins.png')
    fig.clf()
    return


# -------------------------------------------------------------------------------------
# MAIN --------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------

# choose directory
dirname = pxd.file_import(prompt='Choose the directory containing tiff folders: ')
print(dirname)
save_dir1 = dirname + '/anova_outputs_feb2020/'
print(save_dir1)

# load datas
[dist_by_percentiles, numpix_by_dist_bins, norm_numpix_by_dist_bins] = load_datas(dirname)
data_labels = ['naive1', 'naive2', 'naive3', 'naive4', 'naive5',
                   'disLN1', 'disLN2', 'disLN3', 'disLN4', 'disLN5',
                   'tdLN1', 'tdLN2', 'tdLN3', 'tdLN4', 'tdLN5']

# STATS
# We'll try ANOVA. However, looking at the data distributions, it looks like the tumor data has a higher variance.
# One of the assumptions of ANOVA is that groups have the same variance. If this is a problem, we could try
# using Welch's ANOVA through Pengouin. https://pingouin-stats.org/generated/pingouin.welch_anova.html
# Also use this for loop to create a new dataframe where each measure is stacked (essentially a transpose)
dist_by_percentiles_transposed = pd.DataFrame()
numpix_by_dist_bins_transposed = pd.DataFrame()
for row in range(11):
    print('current row: ' + str(row))
    if row < 10:
        print('dist_by_percentiles row name: ')
        print(dist_by_percentiles.iloc[row].name) # print the current index (which is a name because it's a series)
        print('row: ')
        print(dist_by_percentiles.iloc[row])
        run_anova(dist_by_percentiles.iloc[row], save_dir1, data_labels) # only has 10 rows
    print('norm_numpix_by_dist_bins row name: ')
    print(norm_numpix_by_dist_bins.iloc[row].name) # print the current index (which is a name because it's a series)
    print('row: ')
    print(norm_numpix_by_dist_bins.iloc[row])
    run_anova(norm_numpix_by_dist_bins.iloc[row], save_dir1, data_labels) # has 11 rows

# transpose for plotting
dist_by_percentiles_transposed = transpose_data(dist_by_percentiles)
norm_numpix_by_dist_bins_transposed = transpose_data(norm_numpix_by_dist_bins)
# make plots
data_labels2 = ['tdLN', 'disLN', 'naive']
make_plots(dist_by_percentiles_transposed, norm_numpix_by_dist_bins_transposed, save_dir1, data_labels2)

