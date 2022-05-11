from os.path import expanduser, join, abspath
import os
import numpy as np
import matplotlib.pyplot as plt

from backend import Image


# segmentation paramters
num_sizes = 10    # number sizes parameter
min_rad_nucleus = 5    # approximate radius of the smallest cell in pixels
max_rad_nucleus = 25    # approximate radius of the largest cell in pixels
dapi_threshold = 500    # filter out dim nuclei (for 16-bit image max is 65536, so 500 should be low enough)
shape_channel_num = 1    # put 0 to use dapi for morphology metrics calculations

# output parameters
dpi = 600    # output image quality

# specify input path for image-set and output directory
dir_path = "~/Research/TDA/R/TDA-Microscopy-Pipeline/Example/Images/"

# get absolute paths for input and output directories
dir_path = abspath(expanduser(dir_path))


# doxycycline concentrations
for dox in [0, 25]:
    # iterate over wells/locations
    for well in [1, 2, 3]:
        # get back to directory
        image_path = dir_path + f"/dox_{dox}_well_{well}/"

        # get paths to image files
        dapi_path = image_path + "DAPI.tif"
        nanog_path = image_path + "NANOG.tif"
        gata6_path = image_path + "GATA6.tif"
        ch_paths = [nanog_path, gata6_path]

        # create Image object, read DAPI and other channels
        image = Image(dapi_path, 16, ch_paths, image_path)
        image.read()

        # perform image segmentation and calculate cell metrics
        image.segment_pearson(min_rad_nucleus, max_rad_nucleus, num_sizes, channel=shape_channel_num)
        image.getcellprops()

        # calculate channel expression values then save to CSV
        image.overlay_cells()
        image.save_csv("intensities.csv")
        
        # optional method to adjust parameters
        #   - green is false positive cells
        #   - light blue is correct segmentation
        #   - dark blue is leftover dapi staining (false negative)
        # image.check_segmentation()

