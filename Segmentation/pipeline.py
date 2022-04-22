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
input_path = "~/Research/TDA/R/TDA-Microscopy-Pipeline/Example"
output_path = "~/Research/TDA/R/TDA-Microscopy-Pipeline/Example"

# get absolute paths for input and output directories
input_path = os.path.expanduser(input_path)
output_path = os.path.expanduser(output_path)

# set path for CSV output
csv_output_path = os.path.join(output_path, "intensities_bgc.csv")

# DAPI nuclear stain for identifying cells
dapi_file = "DAPI.tif"
ch_paths = ["NANOG.tif", "GATA6.tif"]

# make image paths full paths
dapi_path = os.path.join(input_path, dapi_file)
for i in range(len(ch_paths)):
    ch_paths[i] = os.path.join(input_path, ch_paths[i])

# create Image object, read DAPI and other channels
image = Image(dapi_path, 16, ch_paths)
image.read()

# perform image segmentation and calculate cell metrics
image.segment_pearson(min_rad_nucleus, max_rad_nucleus, num_sizes, channel=shape_channel_num)
image.getcellprops(output_path, dpi)

# calculate channel expression values then save to CSV
image.overlay_cells()
image.df.to_csv(csv_output_path, index=False)

#print("overlaid mask and image: blue is leftover dapi staining (false negative?),")
#print("green is false positive cells, light blue is correct segmentation")

#mask = np.zeros((image.h ,image.w , 3))
#mask[:, :, 1] = np.array(image.mask, dtype="uint16")
#mask[:, :, 2] = image.dapi / np.mean(image.dapi)

#fig = plt.figure(figsize=(10,10))
#plt.imshow(mask)
#plt.show()

