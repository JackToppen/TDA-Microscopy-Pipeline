# Pipeline for Microscopy Image Analysis

This pipeline was introduced in the *Hartsock et al 2023* manuscript. It consists of three modules: (1) segmentation; (2) discretization; and (3) topological data anaysis (TDA). The pipeline's input is a collection of microscopy images and it main outputs are topological descriptors of multicellular patterns (persistence landscapes and average persistence landscapes).

## Segmentation (Python)
This section has been adapted from *Nikitina et al 2020* to interface with
the subsequent sections of the pipeline. Once downloaded locally, the
locations of the microscopy images should be updated, and some of the
configuration parameters may need to be changed.

Please install the necessary Python dependencies for this section.

```
pip install -r requirements.txt
```

Using Anaconda, the segmentation section can be run using the command below.

```
jupyter notebook
```

## Discretization (R)
It is recommended that *pipeline.Rmd* is run in RStudio. The *cell_colors.csv*
file is used to specify how cells are labeled in the discretized images. Binary
strings encode high (1) or low (0) signal of a marker from the microscopy images.
The mapping from binary string to cell color should follow the order of the
measured signals from the segmentation section CSVs.

## TDA (R)
There are two R notebook files within this section. *general_pipeline.Rmd* should
be used for the example image set and for most applications of the pipeline. The
second notebook was configured specifically to the images used within the Hartsock
et al 2023 manuscript.

### Other
Example microscopy images can be found on [Figshare](https://figshare.com/projects/TDA_Microscopy_Data/148855). 
Additionally, *slicer.py* can divide larger microscopy images, as with
this work, into smaller slices.
