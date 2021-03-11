# This script uses pixel_distance.py

import pixel_distance as pxd
import os

# User selects a directory which contains all the directories with tiffs.
# Each folder contains tumor and lyve1 tiff files that go together...
# Tumor will be CH4
# Lyve1 will be CH2
# The program will process all of these pairs one at a time.

def image_testing(tumor_name, lyve1_name, dir_path, save_dir, save_name):
    print('begin image testing method....')
    print('dir_path: ' + dir_path)
    print('save_dir: ' + save_dir)
    print('save_name: ' + save_name)
    os.mkdir(save_dir)
    file1 = open(save_dir + save_name + '.txt', 'a+')
    file1.write('Testing run.....\n')
    file1.write('tumor file: ' + tumor_name + '\n')
    file1.write('lyve1 file: ' + lyve1_name + '\n')
    file1.write('main directiory: ' + dir_path)


# Process all files within sub-directories of the main directory
def process_files(maindir):
    tumor_list = []
    lyve1_list = []

    # Iterate through all files in the directory to create a list
    # If the file contains C2 or C4, add to the list. those are the lyve1 and tumor tiffs.
    # Find the pairs that go together, analyze those and remove them from the list.
    # Once the list is empty, you've analyzed all the files!
    print('start1 ------------------------------------')
    for root, dirs, files in os.walk(maindir):
        for name in files:
            if 'CH2' in name:
                lyve1_list.append(os.path.join(root, name))
            if 'CH4' in name:
                tumor_list.append(os.path.join(root, name))

    # sort lists
    tumor_list.sort()
    lyve1_list.sort()

    print('tumor list and length:')
    print(tumor_list)
    print(len(tumor_list))
    print('lyve1 list and length:')
    print(lyve1_list)
    print(len(lyve1_list))
    print('start2 -----------------------------------')
    # Now go through the file list and process matches.
    counter = 0
    while len(tumor_list) != 0 and len(lyve1_list) != 0:
        for tumor_tif in tumor_list:
            for lyve1_tif in lyve1_list:
                if lyve1_tif.replace('CH2', 'CH4') == tumor_tif:
                    print('tumor file: ' + tumor_tif)
                    print('lyve1 file: ' + lyve1_tif)
                    counter += 1
                    print('PROCESSING SAMPLE.... ' + str(counter) + '-----------------------------------')
                    # Change this call to go between initial testing and actually running
                    pxd.image_handling(tumor_tif,
                                  lyve1_tif,
                                  root + os.sep,
                                  tumor_tif.replace('.tif', 'output' + os.sep),
                                  os.path.split(tumor_tif)[1].replace('C4.tif', ''))
                    # image_testing(tumor_tif,
                    #               lyve1_tif,
                    #               root + os.sep,
                    #               tumor_tif.replace('.tif', 'output' + os.sep),
                    #               os.path.split(tumor_tif)[1].replace('CH4', ''))
                    print('Program completed for sample: ' + lyve1_tif + ' with ' + tumor_tif)
                    lyve1_list.remove(lyve1_tif)
                    tumor_list.remove(tumor_tif)

        print('after processing: ' + str(counter) + 'times...')
        # re-sort lists
        tumor_list.sort()
        lyve1_list.sort()
        print('tumor list and length:')
        print(tumor_list)
        print(len(tumor_list))
        print('lyve1 list and length:')
        print(lyve1_list)
        print(len(lyve1_list))

    # return r
    return

# prompt user for the main directory
# The main directory should contain a subdirectory for all samples.
# Each subdirectory should contain two tiff files:
# one tiff file for the tumor cell marker, and one tiff file for the lyve1 channel.
dirname = pxd.file_import(prompt='Choose the directory containing tiff folders: ')
print(dirname)

# Run program for all files:
process_files(dirname)


