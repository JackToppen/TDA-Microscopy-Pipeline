""":
Adapated from Nitika et al 2022.
"""
import os
import math
import cv2

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from scipy import ndimage
from matplotlib import colors
from matplotlib import cm
from skimage.measure import regionprops
from skimage.filters import gaussian, threshold_local, threshold_minimum
from skimage.feature import peak_local_max
from skimage.segmentation import watershed
from skimage.morphology import label, area_closing


def read_image(path):
    """ Read microscopy images such that the origin
        for the cell locations is in the bottom-left.
    """
    # read the entire image and rotate to such that origin is bottom-left
    image = cv2.imread(path, -1)
    return cv2.rotate(image, cv2.ROTATE_90_CLOCKWISE)


class Channel:
    """ Object for storing each expression channel.
    """
    def __init__(self, image, name):
        self.image = image
        self.name = name
        self.plot = None
        self.intensities = None


class Image:
    """ Object for processing each microscopy image.
    """
    def __init__(self, dapi_path, bit_depth, channels_paths, output_path): 
        self.dapi_path = dapi_path
        self.bit_depth = bit_depth
        self.channels_paths = channels_paths
        self.output_path = output_path

        self.dpi = 600
        self.channels = []
        self.cytoprops = []
        self.cell_count = 0
        self.dapi = None
        self.X = []
        self.Y = []
        self.df = pd.DataFrame()

    def read(self):
        """ Loads all expression channels into separate
            channel objects.
        """
        # read image and get shape
        self.dapi = read_image(self.dapi_path)
        self.h, self.w = self.dapi.shape

        # store DAPI as the first channel and normalize the image
        self.channels = [Channel(self.dapi, "dapi")]
        cv2.normalize(self.dapi, self.dapi, 0, 2 ** self.bit_depth, cv2.NORM_MINMAX)

        # read the other images to separate Channel objects
        for ch_path in self.channels_paths:
            self.channels.append(Channel(read_image(ch_path), os.path.splitext(os.path.basename(ch_path))[0]))

    def getcellprops(self):
        """ Updates data frame with cells that have area larger than 5.
        """
        # create temporary variables for updating data frame
        cells, area, perimeter, eccentricity = [],[],[],[]

        # iterate through cells
        for i in range(self.cell_count):
            # if the cell area is larger than 5, continue adding to lists
            if self.cytoprops[i]["area"] > 5:
                # append the cell metric data to each list
                cells.append(self.cytoprops[i])
                self.X.append(self.cytoprops[i]["centroid"][0])
                self.Y.append(self.cytoprops[i]["centroid"][1])
                area.append(self.cytoprops[i]["area"])
                perimeter.append(self.cytoprops[i]["perimeter"])
                eccentricity.append(self.cytoprops[i]["eccentricity"])
        
        # update the cell count and list of cells
        self.cell_count = len(cells)
        self.cytoprops = cells

        # update cell metrics in data frame
        self.df = pd.DataFrame({"X":self.X, "Y":self.Y, "area":area, "perimeter":perimeter, "eccentricity":eccentricity})

    def save_csv(self, filename):
        """ Saves a CSV with cell metrics data.
        """
        self.df.to_csv(self.output_path + filename)

    def check_segmentation(self):
        """ Optional method to check segmentation of confocal image.
        """
        # merge segmentation mask and DAPI confocal image
        mask = np.zeros((self.h, self.w, 3))
        mask[:, :, 1] = np.array(self.mask, dtype="uint16")
        mask[:, :, 2] = self.dapi / np.mean(self.dapi)
        mask = np.clip(mask, 0, 1)

        # create and save plot
        fig = plt.figure(figsize=(10, 10))
        plt.imshow(mask)
        plt.savefig(self.output_path + "segmentation.png", dpi=self.dpi, bbox_inches="tight")
        plt.close("all")

    def overlay_cells(self):
        """ Measure the channel intensities per cell.
        """
        # go through each Channel object
        for channel in self.channels:
            # empty arrays for plotting and storing intensities 
            channel.plot = np.zeros((self.h, self.w))
            channel.intensities = np.zeros(self.cell_count)

            # go through each cell
            for i in range(self.cell_count):
                # add intensity from each pixel the cell covers
                for coord in self.cytoprops[i]["coords"]:
                    channel.intensities[i] += channel.image[coord[0], coord[1]]
                
                # normalize the cell's total intensity based on its area
                channel.intensities[i] = channel.intensities[i] / self.cytoprops[i]["area"]

                # mark the intensity on the plot
                for coord in self.cytoprops[i]["coords"]:
                    channel.plot[coord[0], coord[1]] = channel.intensities[i]

            # save a plot of the channel intensities and update the data frame
            self.plot(channel.plot, os.path.join(self.output_path, channel.name + ".png"))
            self.df[channel.name] = channel.intensities

    def plot(self, image, name, max_value=0):
        """ Saves a plot of the segmentated image with respective
            channel intensities.
        """
        # set maximum value
        if max_value == 0:
            m = np.max(image)
        else:
            m = max_value

        # create custom color map
        cmap = colors.ListedColormap(['white', 'darkblue', 'blue', 'cornflowerblue', 'cyan', 'aquamarine', 'lime',
                                      'greenyellow','yellow','gold','orange','red','brown'])
        bounds = [0, 1e-9, m/12, m/6, m/4, m/3, m/2.4, m/2, 0.58*m, m/1.5, 0.75*m, m/1.2, 0.91*m, m]
        norm = colors.BoundaryNorm(bounds, cmap.N)
        
        # transform image to display correctly
        image = cv2.flip(image, 0)
        image = cv2.rotate(image, cv2.ROTATE_90_CLOCKWISE)

        # read the image to a plot and save the figure
        plt.imshow(image, cmap=cmap, norm=norm, origin="lower")
        plt.savefig(name, dpi=self.dpi, bbox_inches="tight")
        plt.close('all')

    def adaptive_mask(self, img, high, num_sizes):
        """ Perfroms adaptive mask on confocal image.
        """
        # apply Gaussian blur to the image then normalize the values
        confocal_blur = cv2.GaussianBlur(img, (3, 3), 0)
        confocal_blur = 255 * confocal_blur.astype(np.float64) / np.max(confocal_blur.astype(np.float64))

        # mask for performing adaptive threshold at multiple window sizes
        total_mask = np.zeros((self.h,self.w))

        # choose equally spaced kernel sizes between max cell diameter and minimum image side
        for sz in np.linspace(high * 2, np.min((self.h, self.w)), num_sizes):
            # round floats and make sure they are odd
            sz = int(round(sz))
            if sz % 2 == 0:
                sz += 1

            # create confocal mask
            confocal_mask = cv2.adaptiveThreshold(confocal_blur.astype('uint8'), 255, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY, sz, 0)

            # bitwise "or" to combine this mask with the total mask
            total_mask = np.bitwise_or(total_mask.astype(np.uint8), confocal_mask.astype(np.uint8))

        return total_mask

    def segment_pearson(self, low, high, num_sizes, channel=0):
        """ Performs the image segmentation.
        """
        # step size for edge detection
        step = 1

        # get the median and std deviation the confocal image DAPI intensities
        med = np.median(self.dapi)
        std = np.std(self.dapi)

        # set thresholds for outliers (to eliminate), make sure lower threshold is at least 0
        up_thr = int(round(med + .95 * std))
        low_thr = int(round(med - .95 * std))
        if low_thr < 0:
            low_thr = 0

        # limit values to upper threshold
        confocal_2 = cv2.threshold(self.dapi, up_thr, 255, cv2.THRESH_TRUNC)[1] # only used for entire colony mask

        # generate approximate mask of entire colony to later remove noise
        kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (high * 4 + 1, high * 4 + 1))
        closed = cv2.morphologyEx(confocal_2, cv2.MORPH_CLOSE, kernel, iterations=1)

        # normalize and apply Otsu's Thresholding 
        closed = 255 * closed.astype(np.float64) / np.max(closed)
        m = cv2.threshold(closed.astype(np.uint8), 0, 255, cv2.THRESH_OTSU + cv2.THRESH_BINARY)[1]

        # Gaussian blurring in preparation for adaptive thresh
        total_mask = self.adaptive_mask(self.dapi, high, num_sizes)

        # find range of sigmas using range of radii -> std deviation sigma of gaussian = radius/sqrt(2)
        sigma = np.arange(low, high + step, step) / math.sqrt(2)
        agg = np.zeros((confocal_2.shape[0], confocal_2.shape[1], len(sigma)))

        # Laplacian of Gaussian (LoG) for edge detection
        for i, s in enumerate(sigma):
            # apply Gaussian then Lapacian
            im = gaussian(self.dapi, s)
            im = cv2.Laplacian(im, cv2.CV_64F, ksize=1)

            # normalize by dividing by squared sigma
            n_im_norm = im / s**2

            # invert so blob centers are local maxima rather than minima
            n_im_inv = np.max(n_im_norm) - n_im_norm

            # add to array storing all LoG results
            agg[:, :, i] = n_im_inv

        # find maximum pixel value across all LoG results at different sigma values
        max_agg = np.amax(agg, axis=2)

        # find local minima, minimum distance between maxima is minimum radius of nucleus
        loc_inds = peak_local_max(max_agg, min_distance=low)
        locs = np.zeros(max_agg.shape)
        locs[loc_inds[:, 0], loc_inds[:, 1]] = 1

        # use approximate mask of colony, remove false seeds generated by noise
        locs[m==0] = 0

        # label seeds randomly
        labels = label(locs)

        # perform watershed to get cell metrics
        if channel == 0:
            nuclei_labels = watershed(total_mask, labels, mask=total_mask)
        else:
            temp_image = cv2.GaussianBlur(-self.channels[channel - 1].image, (3,3), 0)
            temp_mask = self.adaptive_mask(self.channels[channel - 1].image, 1.5 * high, num_sizes)
            nuclei_labels = watershed(temp_image, labels, mask=temp_mask, watershed_line=True)
    
        # combine total mask and approximate mask
        self.mask = total_mask & m

        # update cell metrics
        self.cytoprops = regionprops(nuclei_labels)
        self.cell_count = len(self.cytoprops)

