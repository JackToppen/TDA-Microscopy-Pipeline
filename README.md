# Pipeline for Microscopy Image Analysis using TDA

The pipeline has been introduced in *Hartsock
et al 2023* manuscript. It consists of three sequential modules: (1) Segmentation; (2) Cell type identification; and (3) Topological data anaysis (TDA). The pipeline's input is a collection of microscopy images and it main outputs are topological descriptors of multicellular patterns (persistence landscapes and average persistence landscapes).

### Microscopy data
Two datasets and an example from *Hartsock et al 2023* manuscript can be found on [Figshare](https://figshare.com/projects/TDA_Microscopy_Data/148855). The "Dox Treatment Groups" dataset consists of multiple microscopy images of cells labeled with DAPI, NANOG and either pan-GATA6 or HA markers across four Doxycycline (Dox) treatment groups from 0 ng/mL to 25 ng/mL, while the "Additional Co-Stained" dataset includes images where the cells are stained with all four markers. The "Example" images only have two of the four Dox groups (low: 0, high: 25) each containing 16 small images sampled from the "Dox Treatment Groups" dataset.

Provided script `slicer.py` can divide larger microscopy images, as was done in *Hartsock et al 2023* manuscript, into smaller patches. This is how the "Example" dateset was generated. For our pipeline, we suggest using images/patches of size at most 1000x1000 pixels.

## 1. Segmentation - Python
This section has been adapted from *Nikitina et al 2020* to interface with
the subsequent sections of the pipeline. Once downloaded locally, the
locations of the microscopy images should be updated, and some of the
configuration parameters may need to be changed.

Please install the necessary Python dependencies for this section.

```
pip install -r requirements.txt
```

The segmentation section should be run with a [Jupyter Notebook](https://jupyter.org/). After installation, open Jupyter with the command line below.

```
jupyter notebook
```

Note: if using alternative segmentation methods, the cell type identification section requires
CSVs formatted from left to right as follows.

```
| X coord | Y coord | nuclear marker | marker #1 | marker #2 | ... |
```
An example CSV following proper segmentation of the Example dataset, is provided in the module.


## 2. (Cell type) Identification - R
We recommend running all R Markdown files in [RStudio](https://posit.co/download/rstudio-desktop/). Similar to the Segmentation module, please update local paths to the segmentation outputs

The `cell_colors.csv` file is necessary to specify how cell types are shown in the psuedo-microscopy images. Binary strings encode low (0) or high (1) discretized signals of one or more markers from the segmented microscopy images. The mapping from binary string to cell type color should follow the order of the measured signals from left to right in the segmentation CSVs as shown in the example header above.

An example CSV and image following proper cell type identification of the Example dataset, is provided in the module.


## 3. Topological Data Analysis (TDA) - R
This module contains two R folders: `general_pipeline` and `hipsc_pipeline`.

1) `general_pipeline`
   
   This folder should be used for the example dataset and for the most applications of the pipeline. The file `general_pipeline.Rmd` is designed to analyse two discretized microscopy data groups, one at a time, and in     the end compare their tda outputs. Different stages of the analysis are separated into distinct code blocks.

   - **If using our cell type indentification method:** After specifying location of CSV files of a discretized microscopy data group, choose which cell types to use for computations.
   - **If using different cell type identification method:** Ensure that (x, y) coordinates of cells of different cell types are saved in different CSV files. The TDA analysis can then be applied to the entire CSV file      of each specific cell type.

   For each image, the pipeline:

   - Generates and plots the persistence diagram.
   - Plots representative cycles of detected holes in multicellular patterns.
   - Computes the persistence landscape based on the persistence diagram, and plots it (there is an option to save persistence landscapes in CSV files, to facilitate further analysis and reuse).
     
   Then, for each data group, the pipeline produces and plots the average persistence landscape.

   After the average persistence landscapes are computed for both data groups, the difference of two average persistence landscapes is plotted.

   Lastly, the permutation test is performed on the persistence landscapes of the two data groups, and the p-value is calculated.

   Sample outputs for the example dataset can be viewed in `general_pipeline.html`.    

   **Note:** All plots have an option to be saved automatically, enabling to retain visual results for each analysis. 
   
3) `hipsc_pipeline`

   This folder was configured specifically to the complete dataset used within the *Hartsock
   et al 2023* manuscript. It requires installing [tda-tools](https://github.com/jjbouza/tda-tools) by Jose Bouza. Besides TDA analysis, `hipsc_pipeline.Rmd` also contains machine learning and statistical hypothesis       testing computations.

   **Note:** `general_pipeline` can be used for the discretized complete dataset, but first, images need to be sliced into smaller patches with `slicer.py`. 
