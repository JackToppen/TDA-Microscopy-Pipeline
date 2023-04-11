import os
import numpy as np
import tifffile as tiff


# specify input/output directories and number of slices (e.g., n=4 leads to 16 images)
input_dir = "../Nanog_Gata6"
output_dir = "out"
match = ".tif"   # search criteria for images (e.g. ".tif")
n = 4

# expand paths
input_dir = os.path.abspath(os.path.expanduser(input_dir))
output_dir = os.path.abspath(os.path.expanduser(output_dir))

# recursively search through directory
images = []
for root, dirs, files in os.walk(input_dir):
    for file in files:
        # check if search criteria is a subset of file name
        if match in file:
            # get relative path and add to holder
            relative = os.path.relpath(root, input_dir)
            images.append(os.path.join(relative, file))

# iterate through images
for image_file in images:
    # read the image (pip install tifffile)
    image = tiff.imread(os.path.join(input_dir, image_file))

    # get number of rows and columns
    n_rows = image.shape[1]
    n_cols = image.shape[2]

    # get slice sizes
    row_slice = int(n_rows / n)
    col_slice = int(n_cols / n)

    # split the image up to get the relative path and name without extension
    image_rel, image_name = os.path.split(image_file)
    no_ext = os.path.splitext(image_name)[0]

    # create path to new directory for each image
    out_path = os.path.join(output_dir, os.path.join(image_rel, no_ext))
    os.makedirs(out_path, exist_ok=True)

    # loop over rows
    for i in range(n):
        # calculate bounds for slice
        row_low = i * row_slice
        row_up = (i+1) * row_slice

        # don't go over/under for last slice
        if i == n-1:
            row_up = n_rows-1

        # loop over columns
        for j in range(n):
            # calculate bounds for slice
            col_low = j * col_slice
            col_up = (j+1) * col_slice

            # don't go over/under for last slice
            if j == n-1:
                col_up = n_cols-1

            # slice the image and write to new file
            im_slice = image[:, row_low:row_up, col_low:col_up]
            tiff.imwrite(os.path.join(out_path, f"{no_ext}_{i}{j}{match}"), im_slice)

