# Pipeline for Microscopy Image Analysis using TDA

The pipeline has been introduced in *Hartsock
et al 2023* manuscript. It consists of three sequential modules: (1) segmentation; (2) cell type identification; and (3) topological data anaysis (TDA). The pipeline's input is a collection of microscopy images and it main outputs are topological descriptors of multicellular patterns (persistence landscapes and average persistence landscapes).

### Microscopy data
The `example dataset` and the `complete dataset` used in *Hartsock
et al 2023* manuscript can be found on [Figshare](https://figshare.com/projects/TDA_Microscopy_Data/148855). The example dataset consists of two different Dox (doxycycline) treatment groups each containing 16 small images. The complete dataset consists of two image groups based on biological markers pan-GATA6 and HA. Each group contains 4 various Dox treatment subgroups and every one of them consists of 15 large images.

`slicer.py` can divide larger microscopy images, as was done in *Hartsock
et al 2023* manuscript, into smaller patches. For our pipeline, we suggest using images/patches of size at most $1000 \times 1000$ pixels.

## 1. Segmentation (Python)
This section has been adapted from *Nikitina et al 2020* to interface with
the subsequent sections of the pipeline. Once downloaded locally, the
locations of the microscopy images should be updated, and some of the
configuration parameters may need to be changed.

Please install the necessary Python dependencies for this section.

```
pip install -r requirements.txt
```

The segmentation section can be run in the Jupyter Notebook. After instalation, open Jupyter Notebook with the command line below.

```
jupyter notebook
```

Note: if using alternative segmentation methods, the cell type identification section requires
CSVs formatted from left to right as follows.

```
| X coord | Y coord | nuclear marker | marker #1 | marker #2 | ... |
```

## 2. Cell type identification (R)
We recommend running all Rmd files in RStudio.
The `cell_colors.csv` file is used to specify how cells are labeled in the discretized images. Binary strings encode high (1) or low (0) signal of a marker from the microscopy images. The mapping from binary string to cell color should follow the order of the measured signals from the segmentation section CSVs. For instance, by running `pipeline.Rmd` on a segmented example dataset, each cell is assigned to one of the four cell types: 0, 1, 10, or 11.  

## 3. TDA (R)
There are two R folders within this section. `general_pipeline` folder should
be used for the example dataset and for the most applications of the pipeline. The file `general_pipeline.Rmd` is designed to analyse two discretized microscopy data groups one at a time and in the end compare their tda outputs. Note that different stages of the analysis in the file are splitted into separate blocks. After specifying location of CSV files of a discretized microscopy data group, choose which cell types to use for computations. Next, for each image, the persistence diagram is generated and then plotted. The representative cycles of detected holes in multicellular patterns are also plotted. Later, using persistence diagram computations, the persistence landscapes are generated and then ploted. Furthermore, for each data group, the average persistence landscape is produced and plotted. After persistence diagrams, persistence landscapes, and average persistence landscapes are computed for both data groups, the difference of two average persistence landscapes is plotted. Also, the permutation test is performed on the persistence landscapes of the two data groups, and the p-value is calculated. Lastly, in `general_pipeline.html` some of the outputs for the example dataset are shown.    

The second folder `hipsc_pipeline` was configured specifically to the complete dataset used within the *Hartsock
et al 2023* manuscript. It requires installing [tda-tools](https://github.com/jjbouza/tda-tools) by Jose Bouza. Besides TDA analysis, this notebook also contains machine learning and statistical hypothesis testing computations. Note that you can also use `general_pipeline` for the discretized complete dataset but you need to slice images into smaller patches first with `slicer.py`. 
