# PixelDistance
Using multi-channel microscopy images, calculates the minimum distance between two markers. In this case, it was used to calculate how invasive tumor cells were in lymph nodes by measuring the minimum distance between tumor-positive pixels and pixels positive for a lymphatic vessel marker (LYVE-1).

main.py runs pixel_distance.py. 

pixel_distance.py actually performs the measurement of minimum distance between tumor and lyve-1 pixels, and outputs the results for each image.

PixDistStats.py performs stats and makes plots on ALL the data separated by sample group. However, this is insufficient because it isn't split up into biological replicates, or normalized. Use PixDistStats2.py + PixDistStats3.py instead. 

PixDistStats2.py separates the data into biological replicates instead of aggregating all data for each sample group, and experiments with plots.

PixDistStats3.py takes data from PixDistStats2, normalizes it to total pixels for each animal, does statistical comparisons and makes plots.

These are un-polished scripts written by Ian Dryg and used for data analysis in the Lund Lab (thelundlab.com) in ~Feb2021. 
