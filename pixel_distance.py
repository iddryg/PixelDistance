# PixelDistance3
#
# updated 12/3/2020
#
# This script (pixel_distance.py) actually performs the measurement of minimum distance
# between tumor and lyve-1 pixels, and outputs the results for each image.
# PixDistStats.py performs stats and makes plots on ALL the data separated by sample group. However,
# this is insufficient because it isn't split up into biological replicates, or normalized.
# PixDistStats2.py separates the data into biological replicates instead of aggregating
# all data for each sample group, and experiments with plots.
# PixDistStats3.py takes data from PixDistStats2, normalizes it to total pixels for each animal,
# does statistical comparisons and makes plots.

from aicsimageio import AICSImage, imread
import numpy as np
# import napari
import seaborn as sns
import matplotlib.pyplot as plt
from skimage import data, exposure, filters, viewer
from scipy import spatial
import tkinter as tk
from tkinter import filedialog
from PIL import Image
import os
import shutil


# Part 1
# Get user input for files to use (now: directory)
def file_import(prompt):
    print(prompt)
    root = tk.Tk()
    root.lift()
    root.wm_attributes('-topmost', True)
    root.withdraw()

    # file_path = tk.filedialog.askopenfilename()
    root.update()
    # old: asks for file
    # file_name = tk.filedialog.askopenfilename(parent=root, title=prompt)
    # new: asks for directory (containing all folders with pairs of tiffs inside)
    dir_name = tk.filedialog.askdirectory(parent=root, title=prompt)

    root.destroy()
    # old:
    # return file_name
    # new:
    return dir_name


# image_testing(tumor_name, lyve1_name, dir_path, save_dir, save_name):
# Part 2
# Image handling... uses img_subtract_background and img_process
# Note: this function is called by the main function for every tumor-lyve1 image pair
def image_handling(tumor_name, lyve1_name, dir_path, save_dir, save_name):
    if os.path.exists(save_dir):
        shutil.rmtree(save_dir)
    print('Beginning program...')
    # With the Keyence, images will be saved in RGB format. Need to convert to grayscale with .convert('L')
    tumor = np.array(Image.open(tumor_name).convert('L'))
    lyve1 = np.array(Image.open(lyve1_name).convert('L'))
    print('shape of tumor array: ' + str(np.shape(tumor)))
    print('shape of tumor array: ' + str(np.shape(lyve1)))
    thresh_tumor = img_process(tumor)
    thresh_lyve1 = img_process(lyve1)

    # Check Image Processing Quality (Optional)
    # Plot original and binarized images
    threshplot = 0
    if threshplot == 1:
        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12, 8))
        axes[0, 0].imshow(tumor, aspect="auto")
        axes[1, 0].imshow(thresh_tumor, aspect="auto")
        axes[0, 1].imshow(lyve1, aspect="auto")
        axes[1, 1].imshow(thresh_lyve1, aspect="auto")
        plt.tight_layout()
        plt.show()

    print('Part 2 - Process Images: Complete')
    # optional... if there is channel bleedover. Seems to be better now with new antibodies (Dec 2020)
    rmv_duplicates = 0
    if rmv_duplicates == 1:
        thresh_tumor = remove_duplicates(thresh_tumor, thresh_lyve1)
        print('Part 3 - Remove Duplicates from Channel Bleed-over: Complete')
    else:
        print('Part 3 - Remove Duplicates from Channel Bleed-over: SKIPPED!')

    thresh_tumor_pos = np.nonzero(thresh_tumor)
    thresh_lyve1_pos = np.nonzero(thresh_lyve1)
    print('Part 4 - Determine Indices of Nonzero Elements: Complete')

    # Change data shape for use in distances method
    # print('Attempting to stack Index arrays...')
    arr1_xy = np.column_stack([thresh_tumor_pos[1], thresh_tumor_pos[0]])
    # print('Stacked array1')
    arr2_xy = np.column_stack([thresh_lyve1_pos[1], thresh_lyve1_pos[0]])
    # print('Stacked array2')

    # Calculate distances between all combinations of points
    print('Calculating distances...')
    # print(np.shape(arr1_xy))
    # print(np.shape(arr2_xy))
    min_dists = distances_cKDTree(arr1_xy, arr2_xy)
    print('Part 5 - Calculate Minimum Distances Between Tumor & Lymphatic Pixels: Complete')

    # plot and save
    plot_save(min_dists, save_dir, save_name)

    return


# Part 2
# Background Subtraction Function
def img_subtract_background(image, radius=50, light_bg=False):
    from skimage.morphology import white_tophat, black_tophat, disk, ball
    # you can also use 'ball' here to get a slightly smoother result at the cost of increased computing time
    # str_el = disk(radius)
    str_el = disk(radius)
    if light_bg:
        return black_tophat(image, str_el)
    else:
        return white_tophat(image, str_el)


# Part 2
# gaussian blur and binarize the image with Otsu method
# uses img_subtract_background
def img_process(image):
    import skimage
    # subtract background
    sub_image = img_subtract_background(image)
    print("Background subtraction: done")
    # get kernel size from command line
    # sigma = float(sys.argv[2])
    # blur image
    # blur_image = skimage.filters.gaussian(
    #     sub_image, sigma=(sigma, sigma), truncate=3.5, multichannel=False)
    blur_image = skimage.filters.gaussian(
        sub_image, truncate=3.5, multichannel=False)
    print("Gaussian filter blurring: done")
    # determine value at which to threshold image using Otsu method
    threshold = skimage.filters.threshold_otsu(blur_image)
    thresh_image = blur_image > threshold
    print("Thresholding: done")
    return thresh_image


# Part 3
# If both channels (ch3 and ch4) are positive in the same pixels,
# sets channel 3 to zero (avoids bleedover from ch4)
# Assume: image1 and image2 are the same size and shape
def remove_duplicates(image1, image2):
    for row in range(np.shape(image1)[1]):  # y-size
        for column in range(np.shape(image1)[0]):  # x-size
            if image1[row][column] != 0 and image2[row][column] != 0:
                image1[row][column] = 0

    return image1


# Part 5
# Determine distances between all combinations of x,y coordinates
# Inputs should be ndarrays where each row is a different point coordinate [x, y]
# Inputs should only include coordinates of points that are nonzero
def distances_cKDTree(xy1, xy2):
    # This solution is optimal when xy2 is very large
    tree = spatial.cKDTree(xy2)
    mindist, minid = tree.query(xy1)
    return mindist


# plot and save
def plot_save(min_dists, save_dir, save_name):
    # make save directory, replace existing directory if it already exists.
    if os.path.exists(save_dir):
        shutil.rmtree(save_dir)
    os.mkdir(save_dir)
    sns.set_theme(style="whitegrid")
    # ax = sns.boxplot(y=min_dists)
    # ax = sns.swarmplot(y=min_dists, color=".25")

    # violin & strip plot
    sns.violinplot(y=min_dists, inner='stick', color=".8")
    # sns.stripplot(y=min_dists)
    plt.savefig(save_dir + save_name + 'violin.png')
    plt.clf()

    # histogram plot
    sns.histplot(x=min_dists)
    # ax.set(ylabel='Pixels')
    plt.title('Pixel Distance between tumor and lymphatic vessels');
    plt.xlabel('Distance (pixels)')
    # plt.show()
    plt.savefig(save_dir + save_name + 'histogram.png')

    # Save data as csv file for later analysis
    np.savetxt(save_dir + save_name + '.csv', min_dists, delimiter=", ", fmt='% s')

    # Save some basic stats in a text file
    file1 = open(save_dir + save_name + '.txt', 'a+')
    file1.write('Stats for the sample: ' + save_name + '\n')
    file1.write('Mean: ' + str(len(min_dists)) + '\n')
    file1.write('Mean: ' + str(np.mean(min_dists)) + '\n')
    file1.write('Standard Deviation: ' + str(np.std(min_dists))+ '\n')
    file1.write('Percentiles: \n')
    file1.write('10: ' + str(np.percentile(min_dists, 10)) + '\n')
    file1.write('25: ' + str(np.percentile(min_dists, 25)) + '\n')
    file1.write('50: ' + str(np.percentile(min_dists, 50)) + '\n')
    file1.write('75: ' + str(np.percentile(min_dists, 75)) + '\n')
    file1.write('90: ' + str(np.percentile(min_dists, 90)) + '\n')

    print('Violin plot, histogram plot, stats.txt and data file saved.')
    print('Part 6 - Plot Histogram of Distances: Complete')
    return



