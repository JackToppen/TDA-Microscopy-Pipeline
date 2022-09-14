import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import ndimage
from matplotlib import colors
from PIL import Image
from tifffile import tifffile
from skimage.measure import regionprops
from skimage.filters import gaussian, threshold_local
from skimage.feature import peak_local_max
from skimage.segmentation import watershed
from skimage.morphology import label


# location of images to segment (directory will be searched recursively)
input_path = "/Users/jacktoppen/Research/TDA/Microscopy_images/Test/"
match = "merged.tif"   # search criteria for images (e.g. ".tif")

# output location for segmentation CSVs
output_path = "/Users/jacktoppen/Research/TDA/Microscopy_images/Test/"

# nuclear stain position
nuc_name = "dapi"
nuc_index = 2

# indicate channels (beyond cell nuclei stain) and their channel positions
channels = {"gfp": 1, "rfp": 0}


bit_depth = 16    # change to 8 if your image intensity maximum is 255


def get_images(path, match):
    """ Returns a list of file names after searching
        a directory recursively with search criteria.
    """
    # holder for images
    images = []

    # recursively search through directory
    for root, dirs, files in os.walk(path):
        for file in files:
            # check if search criteria is a subset of file name
            if match in file:
                # get relative path and add to holder
                relative = os.path.relpath(root, path)
                images.append(os.path.join(relative, file))

    return images


def segmentation(image):
    """ Performs the image segmentation computations.
    """
    # segmentation parameters
    blur_radius = 0.5
    segmentation_tile_size = 101
    substract_background = 0

    # generate thresholds for masks
    dapi_thresh = threshold_local(image, 3501)
    seg_thresh = threshold_local(image, segmentation_tile_size, offset=substract_background)

    # use thresholds to generate masks
    dapi_mask = image>dapi_thresh
    seg_mask = image>seg_thresh

    # take element-wise AND of the two masks
    mask = dapi_mask&seg_mask

    # get coordinates of local max signal
    distance = gaussian(ndimage.distance_transform_edt(mask), sigma=blur_radius)
    local_max_coords = peak_local_max(distance, exclude_border=False, labels=mask)

    # create mask for local max signals
    local_max_mask = np.zeros(distance.shape, dtype=bool)
    local_max_mask[tuple(local_max_coords.T)] = True

    # find individual cells
    markers = label(local_max_mask)
    nuclei_labels = watershed(-distance, markers, mask=mask)
    cytoprops = regionprops(nuclei_labels)

    return cytoprops


def identify_cells(cellprops):
    """ Determines which detected objects are cells.
    """
    # get number of objects detected
    n = len(cellprops)

    # create arrays to hold cell values
    x_pos = np.zeros(n)
    y_pos = np.zeros(n)
    area = np.zeros(n)
    cells = np.array(cellprops)

    # get the per-cell properties (x, y, area)
    for i in range(n):
        x_pos[i] = cellprops[i]["centroid"][0]
        y_pos[i] = cellprops[i]["centroid"][1]
        area[i] = cellprops[i]["area"]

    # get standard deviation and mean for cell areas
    std, mean = np.std(area), np.mean(area)

    # go through
    mask = np.zeros(n, dtype=bool)
    for i in range(n):
        if area[i] < (mean + 3 * std) and area[i] > max(5, mean - 3 * std):
            mask[i] = 1

    # apply mask to positions and detected objects
    x_pos, y_pos = x_pos[mask], y_pos[mask]
    cells = cells[mask]

    return x_pos, y_pos, cells


def overlay_cells(image, cells):
    """ Based on segmentation results, quantify individual
        cell signals.
    """
    # create image for outputting intensities
    out_image = np.zeros(image.shape)

    # count the number of cells and create array to store per cell intensities
    n = len(cells)
    intensities = np.zeros(n)

    # iterate through all cells
    for i in range(n):
        # add signal to a per cell holder for each pixel of the cell
        for coord in cells[i]["coords"]:
            intensities[i] += image[coord[0], coord[1]]

        # divide the intensities by the cell area
        intensities[i] /= cells[i]["area"]

        # set the cell pixels to the average expression
        for coord in cells[i]["coords"]:
            out_image[coord[0], coord[1]] = intensities[i]

    return intensities, out_image


def plot(image, name, bits=0):
    """ Outputs an image with overlayed segmentation/
        intensities to gauge segmentation accuracy.
    """
    # output image quality
    dpi = 600

    # based bounds on image bit-depth
    if bits == 0:
        m = np.max(image)
    else:
        m = 2**bits

    # create a custom colormap
    bounds=[0, 1e-9, m/12, m/6, m/4, m/3, m/2.4, m/2, 0.58*m, m/1.5, 0.75*m, m/1.2, 0.91*m, m]
    cmap = colors.ListedColormap(["white","darkblue","blue","cornflowerblue",
                                  "cyan","aquamarine","lime","greenyellow",
                                  "yellow","gold","orange","red","brown"])
    norm = colors.BoundaryNorm(bounds, cmap.N)

    # display the image and save it
    image = plt.imshow(image, cmap=cmap, norm=norm)
    plt.colorbar(image, cmap='jet')
    plt.savefig(name,dpi=dpi,bbox_inches = "tight")
    plt.close()


# get all images in input directory
image_paths = get_images(input_path, match)

# run the segmentation calculation on each image
for path in image_paths:
    # read the image
    image = tifffile.imread(input_path + path)

    # get the nuclear stain channel and segmen the image based on this
    nuc_channel = image[nuc_index]
    nuc_seg = segmentation(nuc_channel)

    # create path to CSV file
    csv_path = output_path + os.path.splitext(path)[0] + ".csv"

    # make sure directories exist for CSV file
    dir_path = os.path.split(csv_path)[0]
    os.makedirs(dir_path, exist_ok=True)

    # return the detected objects that are identified as cells
    x_pos, y_pos, cells = identify_cells(nuc_seg)

    # get the DAPI intensities and generated image
    nuc_intensities, out_image = overlay_cells(nuc_channel, cells)
    out_path = os.path.join(dir_path, nuc_name + ".png")
    # plot(out_image, out_path, bit_depth)    # save the image

    # create dictionary for storing cell specific values
    data = {"X": x_pos, "Y": y_pos, nuc_name: nuc_intensities}

    # iterate through the channels, quantifying signal in each
    for name in channels.keys():
        # get expression channel from image
        index = channels[name]
        channel = image[index]

        # calculate intensities
        intensities, out_image = overlay_cells(channel, cells)
        data[name] = intensities
        out_path = os.path.join(dir_path, name + ".png")
        # plot(out_image, out_path, bit_depth)    # save the image

    # save cell locations and signal intensities to CSV
    df = pd.DataFrame(data)
    df.to_csv(csv_path, index=False)

