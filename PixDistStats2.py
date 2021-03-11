# This script analyzes the csv files output by main / pixel_distance.py
# uses some methods from pixel_distance.py
# Updated Feb 2021.
# This version separates the data into biological replicates instead of aggregating all data for each sample group.

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



def analyze(data, save_dir, labels, bin_size):
    save_name = 'Results'
    print('Bin size is: ')
    print(bin_size)

    # Scale the data
    # 0.754 microns per pixel
    print('data before scaling: ')
    print(data.head(15))
    data = 0.754*data
    print('data after scaling: ')
    print(data.head(15))
    data.replace(0, np.nan, inplace=True)
    print('data after removing zeros: ')
    print(data.head(15))

    stats1 = 0
    if stats1 == 1:
        # Run One-way ANOVA
        # F is the F-statistic (an array, one for each group)
        # p is the p-value (an array, one for each group)
        # f_stat, p = stats.f_oneway(data[labels[0]][~np.isnan(data[labels[0]])],
        #                 data[labels[1]][~np.isnan(data[labels[1]])],
        #                 data[labels[2]][~np.isnan(data[labels[2]])])

        # On ALL the data, unbinned but *by animal*
        f_stat, p = stats.f_oneway(data[labels[0]].dropna(),
                                   data[labels[1]].dropna(),
                                   data[labels[2]].dropna(),
                                   data[labels[3]].dropna(),
                                   data[labels[4]].dropna(),
                                   data[labels[5]].dropna(),
                                   data[labels[6]].dropna(),
                                   data[labels[7]].dropna(),
                                   data[labels[8]].dropna(),
                                   data[labels[9]].dropna(),
                                   data[labels[10]].dropna(),
                                   data[labels[11]].dropna(),
                                   data[labels[12]].dropna(),
                                   data[labels[13]].dropna(),
                                   data[labels[14]].dropna())

        # Make new dataframe to stack the data for tukey comparisons...
        # print(data.head())
        print('number of NaNs in data before stacking: ' + str(data.isna().sum()))
        data_stacked = data.stack().reset_index()
        print('number of NaNs in data_stacked initially: ' + str(data_stacked.isna().sum()))
        data_stacked = data_stacked.dropna().rename(columns={'level_0': 'id', 'level_1': 'group', 0: 'distance'})
        print('number of NaNs in data_stacked after dropna(): ' + str(data_stacked.isna().sum()))
        # make new column with supergroups (naive, ndln, tdln)
        data_stacked['supergroup'] = data_stacked['group'].map(lambda x: x.rstrip('12345'))
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

    plots = 0
    if plots == 1:
        # Violin plot ---------------------------------------------------------------------------------
        sns.set_theme(style="whitegrid")
        ax = sns.violinplot(data=data, inner='quartile', scale='width',
                            cut=0, bw='scott', width=0.9)
        # Save violin plot
        ax.set(xlabel=None)  # remove the axis label
        fig = ax.get_figure()
        fig.set_size_inches(14, 6)
        fig.tight_layout()
        fig.savefig(save_dir + save_name + 'violin.png')
        fig.clf()
        # plt.savefig(save_dir + save_name + 'violin.png')
        # plt.clf()

        # Violin with stacked data -------------------------------------------------------
        ax = sns.violinplot(data=data_stacked, x="group", y="distance", hue="supergroup", inner='quartile', scale='area',
                            cut=0, bw='scott', width=1, scale_hue=False)

        ax.set(xlabel=None)  # remove the axis labelfig = ax.get_figure()
        plt.gca().legend().set_title(None)
        fig.set_size_inches(20, 5)
        # fig.tight_layout()
        fig.savefig(save_dir + save_name + 'violin_supergroup.png')
        fig.clf()

        # Violin with stacked data 2 -------------------------------------------------------
        ax = sns.violinplot(data=data_stacked, x="supergroup", y="distance", hue="supergroup", inner='quartile',
                            scale='count', cut=0, bw='scott', width=0.8)

        ax.set(xlabel=None)  # remove the axis label
        fig = ax.get_figure()
        plt.gca().legend().set_title(None)
        fig.set_size_inches(10, 5)
        fig.tight_layout()
        fig.savefig(save_dir + save_name + 'violin_supergroup2.png')
        fig.clf()

        # Swarm Plot ---------------------------------------------------------------------
        # Swarmplot takes too long... try strip plot. If that doesn't work try scatter plot
        # ax = sns.swarmplot(data=data_stacked, x="group", y="distance", hue="supergroup", size=0.5)
        ax = sns.stripplot(data=data_stacked, x="group", y="distance", hue="supergroup", size=2)
        ax.set(xlabel=None)  # remove the axis label
        fig = ax.get_figure()
        fig.set_size_inches(13, 6)
        plt.gca().legend().set_title(None)
        # fig.savefig(save_dir + save_name + 'swarm_supergroup.png')
        fig.savefig(save_dir + save_name + 'strip_supergroup.png')
        fig.clf()

        # Swarm Plot with stacked data ---------------------------------------------------------------------
        # Swarmplot takes too long... try strip plot. If that doesn't work try scatter plot
        # ax = sns.swarmplot(data=data_stacked, x="supergroup", y="distance", hue="supergroup", size=0.5)
        ax = sns.stripplot(data=data_stacked, x="supergroup", y="distance", hue="supergroup", size=2.5)
        ax.set(xlabel=None)  # remove the axis label
        fig = ax.get_figure()
        fig.set_size_inches(10, 6)
        # fig.savefig(save_dir + save_name + 'swarm_all.png')
        plt.gca().legend().set_title(None)
        fig.savefig(save_dir + save_name + 'strip_all.png')
        fig.clf()

        # RIDGE PLOT ----------------------------------------------------------------------
        # Uses STACKED data
        # Initialize the FacetGrid object
        pal = sns.cubehelix_palette(15, rot=-.25, light=.7)
        g = sns.FacetGrid(data_stacked, row="group", hue="group", aspect=15, height=.5, palette=pal, sharey=True)

        # Draw the densities in a few steps
        g.map(sns.kdeplot, "distance", bw_adjust=.5, clip_on=False, fill=True, alpha=1, linewidth=1.5)
        g.map(sns.kdeplot, "distance", clip_on=False, color="w", lw=2, bw_adjust=.5)
        g.map(plt.axhline, y=0, lw=2, clip_on=False)

        # Set x-limits
        g.set(xlim=(-5, 80))
        # Define and use a simple function to label the plot in axes coordinates
        def label_axes(x, color, label):
            ax = plt.gca()
            ax.text(0, .7, label, fontweight="bold", color=color,
                    ha="left", va="center", transform=ax.transAxes)

        g.map(label_axes, "distance")

        # Set the subplots to overlap
        g.fig.subplots_adjust(hspace=0)

        # Remove axes details that don't play well with overlap
        g.set_titles("")
        g.set(yticks=[])
        g.despine(bottom=True, left=True)

        # title
        plt.suptitle('Tumor Invasion Distance by Animal')
        # fig.set_size_inches(9, 15)
        # uncomment the following line if there's a tight layout warning
        # g.fig.tight_layout()
        # show plot
        # plt.show()
        # Save file
        plt.savefig(save_dir + 'ridgeplot.png')
        plt.clf()

        # RIDGE PLOT  with supergroups ----------------------------------------------------------------------
        # Uses STACKED data
        # Initialize the FacetGrid object
        pal = sns.cubehelix_palette(4, rot=-.4, light=.7)
        g = sns.FacetGrid(data_stacked, row="supergroup", hue="supergroup", aspect=5, height=2, palette=pal, sharey=True)

        # Draw the densities in a few steps
        g.map(sns.kdeplot, "distance", bw_adjust=.5, clip_on=False, fill=True, alpha=1, linewidth=1.5)
        # g.map(sns.kdeplot, "distance", clip_on=False, color="w", lw=2, bw_adjust=.5)
        g.map(plt.axhline, y=0, lw=2, clip_on=False)

        # Set x-limits
        g.set(xlim=(-5, 80))

        # Define and use a simple function to label the plot in axes coordinates
        def label_axes(x, color, label):
            ax = plt.gca()
            ax.text(0, .5, label, fontweight="bold", color=color,
                    ha="left", va="center", transform=ax.transAxes)

        g.map(label_axes, "distance", label="supergroup")

        # Set the subplots to overlap
        g.fig.subplots_adjust(hspace=-0.25)

        # Remove axes details that don't play well with overlap
        g.set_titles("")
        g.set(yticks=[])
        g.despine(bottom=True, left=True)

        # title
        plt.suptitle('Tumor Invasion Distance by Group')
        # fig.set_size_inches(7, 6)
        # uncomment the following line if there's a tight layout warning
        # g.fig.tight_layout()
        # show plot
        # plt.show()
        # Save file
        plt.savefig(save_dir + 'ridgeplot2.png')
        plt.clf()

        # JOYPLOT ------------------------------------------------------------------------------
        # fig, axes = jpy.joyplot(data_stacked, by='group', ylim='own', overlap=1, range_style='own', x_range=[-5, 100])
        fig, axes = jpy.joyplot(data_stacked, by='group', overlap=1, range_style='own', x_range=[-5, 100])
        fig.savefig(save_dir + 'joyplot.png')
        plt.clf()

        # JOYPLOT all ------------------------------------------------------------------------------
        # fig, axes = jpy.joyplot(data_stacked, by='supergroup', ylim='own', overlap=1, range_style='own', x_range=[-5, 100])
        fig, axes = jpy.joyplot(data_stacked, by='supergroup', overlap=1, range_style='own', x_range=[-5, 100])
        fig.savefig(save_dir + 'joyplot_all.png')
        plt.clf()

    # INDIVIDUAL PLOTS ----------------------------------------------------------------------
    # For loop for Histograms and txt file summaries for each group
    # Update Feb 2021: Adding in new analyses using the percentiles
    # dist_by_percentiles:
    # 1. calculate the mean distance for each 10th percentile for each animal
    # 2. then, using those means as biological replicate data points, do statistical comparisons for each group
    #    at each percentile.
    percentile_index = ['10th', '20th', '30th', '40th', '50th', '60th', '70th', '80th', '90th', '100th']
    # dist_by_percentiles = data.iloc[0:0, :].copy()
    # dist_by_percentiles = pd.DataFrame(index=percentile_index, columns=data.columns).fillna(0)
    dist_by_percentiles = pd.DataFrame()
    # percent_by_dist_bins
    # 1. count the number of observations within each distance bin (for example, 0-20 microns, 20-40 microns, etc)
    # 2. normalize those counts to the total counts for each animal. (essentially a % or proportion of cells invading each distance)
    # 3. then do statistical comparisons for each distance comparing those normalized counts.
    dist_bins_index = ['0-10um', '10-20um', '20-30um', '30-40um', '40-50um', '50-60um', '60-70um', '70-80um', '80-90um', '90-100um', '100um+']
    dist_bins = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 500]
    # numpix_by_dist_bins = data.iloc[0:0, :].copy()
    # numpix_by_dist_bins = pd.DataFrame(index=dist_bins_index, columns=data.columns).fillna(0)
    numpix_by_dist_bins = pd.DataFrame()
    # initialize dataframe for normalized numpix_by_dist_bins
    norm_numpix_by_dist_bins = pd.DataFrame()

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

        print('run number: ' + str(counter + 1))
        # Add percentile data into dist_by_percentiles
        percentiles = pd.Series([0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        i = 0
        for p in list(range(10, 110, 10)):
            percentiles.iloc[i] = np.nanpercentile(data[label], p)
            # print(str(percentiles.iloc[i]))
            i += 1

        # this isn't working correctly, try just populating a fresh dataframe! - Feb 3rd 2021
        # dist_by_percentiles.insert(loc=0, column=label, value=percentiles)
        dist_by_percentiles.insert(loc=0, column=label, value=percentiles)
        # dist_by_percentiles[label] = pd.Series([np.nanpercentile(data[label], p) for p in list(range(10, 110, 10))])
        # print(dist_by_percentiles.dtypes)

        # Add distance bin data into numpix_by_dist_bins
        pix_count = pd.Series([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        pix_count = pd.cut(data[label], bins=dist_bins).value_counts()
        numpix_by_dist_bins.insert(loc=0, column=label, value=pix_count)
        #pd.cut(df['ext price'], bins=4).value_counts()
        # numpix_by_dist_bins[label] = pd.cut(data[label], bins=dist_bins).value_counts()
        # print('numpix by dist bins: ')
        # print(numpix_by_dist_bins.head(11))

        # normalized numpix_by_dist_bins
        # normalize each sample by the total number of pixels in that sample
        norm_numpix_by_dist_bins.insert(loc=0, column=label, value=100*pix_count/pix_count.sum())

        counter += 1

    # Now that dataframes are populated, set indexes
    # dist_by_percentiles = dist_by_percentiles.reindex(labels=percentile_index) # doesn't work, all values are NaNs
    # Insert a new column with the index names
    dist_by_percentiles.insert(loc=0, column='percentiles', value=percentile_index)
    # Use the new column as the index names
    dist_by_percentiles = dist_by_percentiles.set_index('percentiles')
    print('dist by percentiles: ')
    print(dist_by_percentiles.head(10))
    print('dist_by_percentiles row names: ')
    print(list(dist_by_percentiles.index))

    # numpix_by_dist_bins.reindex(labels=dist_bins_index)
    print('numpix by dist bins: ')
    print(numpix_by_dist_bins.head(11))
    print('numpix_by_dist_bins row names: ')
    print(list(numpix_by_dist_bins.index))
    numpix_by_dist_bins = numpix_by_dist_bins.sort_index()
    print('numpix by dist bins after sorting by index: ')
    print(numpix_by_dist_bins.head(11))
    # Insert a new column with the index names
    numpix_by_dist_bins.insert(loc=0, column='distance bins', value=dist_bins_index)
    # Use the new column as the index names
    numpix_by_dist_bins = numpix_by_dist_bins.set_index('distance bins')
    # numpix_by_dist_bins = numpix_by_dist_bins.reindex(labels=dist_bins_index) # doesn't work, all values are NaNs
    print('numpix by dist bins after re-indexing: ')
    print(numpix_by_dist_bins.head(11))

    # norm_numpix_by_dist_bins.reindex(labels=dist_bins_index)
    print('norm numpix by dist bins: ')
    print(numpix_by_dist_bins.head(11))
    print('norm_numpix_by_dist_bins row names: ')
    print(list(numpix_by_dist_bins.index))
    norm_numpix_by_dist_bins = norm_numpix_by_dist_bins.sort_index()
    print('norm numpix by dist bins after sorting by index: ')
    print(norm_numpix_by_dist_bins.head(11))
    # Insert a new column with the index names
    norm_numpix_by_dist_bins.insert(loc=0, column='distance bins', value=dist_bins_index)
    # Use the new column as the index names
    norm_numpix_by_dist_bins = norm_numpix_by_dist_bins.set_index('distance bins')
    # numpix_by_dist_bins = numpix_by_dist_bins.reindex(labels=dist_bins_index) # doesn't work, all values are NaNs
    print('norm numpix by dist bins after re-indexing: ')
    print(norm_numpix_by_dist_bins.head(11))

    print('Violin plot, histogram plot, stats.txt and data file saved.')
    print('Part 6 - Plot Histogram of Distances: Complete')

    # Save data as csv file for later analysis
    # np.savetxt(save_dir + 'dist_by_percentiles.csv', dist_by_percentiles, delimiter=", ", fmt='% s')
    # np.savetxt(save_dir + 'numpix_by_dist_bins.csv', numpix_by_dist_bins, delimiter=", ", fmt='% s')
    dist_by_percentiles.to_csv(save_dir + 'dist_by_percentiles.csv')
    numpix_by_dist_bins.to_csv(save_dir + 'numpix_by_dist_bins.csv')
    norm_numpix_by_dist_bins.to_csv(save_dir + 'norm_numpix_by_dist_bins.csv')

    return

def process_csv_files(maindir):
    # There are 3 groups: naive, ndln, tdln
    # Within each group, there are 5 animals
    # Within each animal, there are 3 fields imaged.
    # Want to aggregate the 3 fields imaged and keep animals separate.
    # Naming conventions:
    # 20201103 naive LN 1 20x 1_CH4.tif - naive, animal 1, field 1
    # 20201103 naive LN 1 20x 2_CH4.tif - naive, animal 1, field 2
    # 20201103 ndLN 1 20x 1_CH4.tif - ndln, animal 1, field 1
    # 20201103 tdLN 1 20x 1_CH4.tif - tdln, animal 1, field 1



    # initialize arrays...
    naive_array1 = np.zeros(0)
    naive_size_array1 = np.zeros(45)
    naive_size_count1 = 0
    naive_array2 = np.zeros(0)
    naive_size_array2 = np.zeros(45)
    naive_size_count2 = 0
    naive_array3 = np.zeros(0)
    naive_size_array3 = np.zeros(45)
    naive_size_count3 = 0
    naive_array4 = np.zeros(0)
    naive_size_array4 = np.zeros(45)
    naive_size_count4 = 0
    naive_array5 = np.zeros(0)
    naive_size_array5 = np.zeros(45)
    naive_size_count5 = 0

    ndln_array1 = np.zeros(0)
    ndln_size_array1 = np.zeros(45)
    ndln_size_count1 = 0
    ndln_array2 = np.zeros(0)
    ndln_size_array2 = np.zeros(45)
    ndln_size_count2 = 0
    ndln_array3 = np.zeros(0)
    ndln_size_array3 = np.zeros(45)
    ndln_size_count3 = 0
    ndln_array4 = np.zeros(0)
    ndln_size_array4 = np.zeros(45)
    ndln_size_count4 = 0
    ndln_array5 = np.zeros(0)
    ndln_size_array5 = np.zeros(45)
    ndln_size_count5 = 0

    tdln_array1 = np.zeros(0)
    tdln_size_array1 = np.zeros(45)
    tdln_size_count1 = 0
    tdln_array2 = np.zeros(0)
    tdln_size_array2 = np.zeros(45)
    tdln_size_count2 = 0
    tdln_array3 = np.zeros(0)
    tdln_size_array3 = np.zeros(45)
    tdln_size_count3 = 0
    tdln_array4 = np.zeros(0)
    tdln_size_array4 = np.zeros(45)
    tdln_size_count4 = 0
    tdln_array5 = np.zeros(0)
    tdln_size_array5 = np.zeros(45)
    tdln_size_count5 = 0

    print('start1 ------------------------------------')
    for root, dirs, files in os.walk(maindir):
        for name in files:
            if 'csv' in name and 'naive LN 1' in name:
                # Naive animal 1
                # flattened append
                naive_array1 = np.append(naive_array1, np.genfromtxt(os.path.join(root, name), delimiter=','))
                print('naive array shape 1:')
                print(np.shape(naive_array1))
                naive_size_array1[naive_size_count1] = len(np.genfromtxt(os.path.join(root, name), delimiter=','))
                naive_size_count1 += 1
                print('naive count 1: ' + str(naive_size_count1))
            if 'csv' in name and 'naive LN 2' in name:
                # Naive animal 2
                # flattened append
                naive_array2 = np.append(naive_array2, np.genfromtxt(os.path.join(root, name), delimiter=','))
                print('naive array shape 2:')
                print(np.shape(naive_array2))
                naive_size_array2[naive_size_count2] = len(np.genfromtxt(os.path.join(root, name), delimiter=','))
                naive_size_count2 += 1
                print('naive count 2: ' + str(naive_size_count2))
            if 'csv' in name and 'naive LN 3' in name:
                # Naive animal 3
                # flattened append
                naive_array3 = np.append(naive_array3, np.genfromtxt(os.path.join(root, name), delimiter=','))
                print('naive array shape 3:')
                print(np.shape(naive_array3))
                naive_size_array3[naive_size_count3] = len(np.genfromtxt(os.path.join(root, name), delimiter=','))
                naive_size_count3 += 1
                print('naive count 3: ' + str(naive_size_count3))
            if 'csv' in name and 'naive LN 4' in name:
                # Naive animal 4
                # flattened append
                naive_array4 = np.append(naive_array4, np.genfromtxt(os.path.join(root, name), delimiter=','))
                print('naive array shape 4:')
                print(np.shape(naive_array4))
                naive_size_array4[naive_size_count4] = len(np.genfromtxt(os.path.join(root, name), delimiter=','))
                naive_size_count4 += 1
                print('naive count 4: ' + str(naive_size_count4))
            if 'csv' in name and 'naive LN 5' in name:
                # Naive animal 5
                # flattened append
                naive_array5 = np.append(naive_array5, np.genfromtxt(os.path.join(root, name), delimiter=','))
                print('naive array shape 5:')
                print(np.shape(naive_array5))
                naive_size_array5[naive_size_count5] = len(np.genfromtxt(os.path.join(root, name), delimiter=','))
                naive_size_count5 += 1
                print('naive count 5: ' + str(naive_size_count5))
            if 'csv' in name and 'ndLN 1' in name:
                # ndLN animal 1
                ndln_array1 = np.append(ndln_array1, np.genfromtxt(os.path.join(root, name), delimiter=','))
                print('ndLN array shape 1:')
                print(np.shape(ndln_array1))
                ndln_size_array1[ndln_size_count1] = len(np.genfromtxt(os.path.join(root, name), delimiter=','))
                ndln_size_count1 += 1
                print('ndLN count 1: ' + str(ndln_size_count1))
            if 'csv' in name and 'ndLN 2' in name:
                # ndLN animal 2
                ndln_array2 = np.append(ndln_array2, np.genfromtxt(os.path.join(root, name), delimiter=','))
                print('ndLN array shape 2:')
                print(np.shape(ndln_array2))
                ndln_size_array2[ndln_size_count2] = len(np.genfromtxt(os.path.join(root, name), delimiter=','))
                ndln_size_count2 += 1
                print('ndLN count 2: ' + str(ndln_size_count2))
            if 'csv' in name and 'ndLN 3' in name:
                # ndLN animal 3
                ndln_array3 = np.append(ndln_array3, np.genfromtxt(os.path.join(root, name), delimiter=','))
                print('ndLN array shape 3:')
                print(np.shape(ndln_array3))
                ndln_size_array3[ndln_size_count3] = len(np.genfromtxt(os.path.join(root, name), delimiter=','))
                ndln_size_count3 += 1
                print('ndLN count 3: ' + str(ndln_size_count3))
            if 'csv' in name and 'ndLN 4' in name:
                # ndLN animal 4
                ndln_array4 = np.append(ndln_array4, np.genfromtxt(os.path.join(root, name), delimiter=','))
                print('ndLN array shape 4:')
                print(np.shape(ndln_array4))
                ndln_size_array4[ndln_size_count4] = len(np.genfromtxt(os.path.join(root, name), delimiter=','))
                ndln_size_count4 += 1
                print('ndLN count 4: ' + str(ndln_size_count4))
            if 'csv' in name and 'ndLN 5' in name:
                # ndLN animal 5
                ndln_array5 = np.append(ndln_array5, np.genfromtxt(os.path.join(root, name), delimiter=','))
                print('ndLN array shape 5:')
                print(np.shape(ndln_array5))
                ndln_size_array5[ndln_size_count5] = len(np.genfromtxt(os.path.join(root, name), delimiter=','))
                ndln_size_count5 += 1
                print('ndLN count 5: ' + str(ndln_size_count5))
            if 'csv' in name and 'tdLN 1' in name:
                # tdLN animal 1
                tdln_array1 = np.append(tdln_array1, np.genfromtxt(os.path.join(root, name), delimiter=','))
                print('tdLN array shape 1:')
                print(np.shape(tdln_array1))
                tdln_size_array1[tdln_size_count1] = len(np.genfromtxt(os.path.join(root, name), delimiter=','))
                tdln_size_count1 += 1
                print('tdLN count 1: ' + str(tdln_size_count1))
            if 'csv' in name and 'tdLN 2' in name:
                # tdLN animal 2
                tdln_array2 = np.append(tdln_array2, np.genfromtxt(os.path.join(root, name), delimiter=','))
                print('tdLN array shape 2:')
                print(np.shape(tdln_array2))
                tdln_size_array2[tdln_size_count2] = len(np.genfromtxt(os.path.join(root, name), delimiter=','))
                tdln_size_count2 += 1
                print('tdLN count 2: ' + str(tdln_size_count2))
            if 'csv' in name and 'tdLN 3' in name:
                # tdLN animal 3
                tdln_array3 = np.append(tdln_array3, np.genfromtxt(os.path.join(root, name), delimiter=','))
                print('tdLN array shape 3:')
                print(np.shape(tdln_array3))
                tdln_size_array3[tdln_size_count3] = len(np.genfromtxt(os.path.join(root, name), delimiter=','))
                tdln_size_count3 += 1
                print('tdLN count 3: ' + str(tdln_size_count3))
            if 'csv' in name and 'tdLN 4' in name:
                # tdLN animal 4
                tdln_array4 = np.append(tdln_array4, np.genfromtxt(os.path.join(root, name), delimiter=','))
                print('tdLN array shape 4:')
                print(np.shape(tdln_array4))
                tdln_size_array4[tdln_size_count4] = len(np.genfromtxt(os.path.join(root, name), delimiter=','))
                tdln_size_count4 += 1
                print('tdLN count 4: ' + str(tdln_size_count4))
            if 'csv' in name and 'tdLN 5' in name:
                # tdLN animal 5
                tdln_array5 = np.append(tdln_array5, np.genfromtxt(os.path.join(root, name), delimiter=','))
                print('tdLN array shape 5:')
                print(np.shape(tdln_array5))
                tdln_size_array5[tdln_size_count5] = len(np.genfromtxt(os.path.join(root, name), delimiter=','))
                tdln_size_count5 += 1
                print('tdLN count 5: ' + str(tdln_size_count5))
    print(tdln_size_array5)

    print('data collection complete.')
    print('beginning analysis.')

    # Create a data dictionary to load arrays into a pandas dataframe
    data = {
        'naive1': naive_array1,
        'naive2': naive_array2,
        'naive3': naive_array3,
        'naive4': naive_array4,
        'naive5': naive_array5,
        'disLN1': ndln_array1,
        'disLN2': ndln_array2,
        'disLN3': ndln_array3,
        'disLN4': ndln_array4,
        'disLN5': ndln_array5,
        'tdLN1': tdln_array1,
        'tdLN2': tdln_array2,
        'tdLN3': tdln_array3,
        'tdLN4': tdln_array4,
        'tdLN5': tdln_array5
    }
    # Want all arrays to be same size for pandas
    maxsize = max([a.size for a in data.values()])
    # pad the shorter arrays with NaNs to make them the same length as the longest array
    data_pad = {k: np.pad(v, pad_width=(0, maxsize - v.size,), mode='constant', constant_values=np.nan) for k, v in
                data.items()}
    df = pd.DataFrame(data_pad)
    data_labels = ['naive1', 'naive2', 'naive3', 'naive4', 'naive5',
                   'disLN1', 'disLN2', 'disLN3', 'disLN4', 'disLN5',
                   'tdLN1', 'tdLN2', 'tdLN3', 'tdLN4', 'tdLN5']
    analyze(df, maindir, data_labels, bin_size=10)

    # analyze(naive_array, naive_size_array, maindir, 'naive')
    # analyze(ndln_array, ndln_size_array, maindir, 'ndln')
    # analyze(tdln_array, tdln_size_array, maindir, 'tdln')

    return


dirname = pxd.file_import(prompt='Choose the directory containing tiff folders: ')
print(dirname)

process_csv_files(dirname)

