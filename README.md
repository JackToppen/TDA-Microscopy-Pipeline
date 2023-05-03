# Pipeline for quantifying cell patterning with TDA

The pipeline consists of three modules: (1) segmentation; (2) discretization; and (3) topological data anaysis (TDA). The pipeline's input is a collection of microscopy images and it main outputs are topological descriptors of multicellular patterns (persistence landscapes and average persistence landscapes).

### Microscopy data
The `example image dataset` and the `complete microscopy dataset` used in *Hartsock
et al 2023* manuscript can be found on [Figshare](https://figshare.com/projects/TDA_Microscopy_Data/148855). The example  image dataset consists of two different Dox (doxycycline) treatment groups each containing 16 small images. The complete microscopy dataset consists of two image groups based on biological markers pan-GATA6 and HA. Each group contains 4 various Dox treatment subgroups and every one of them consists of 15 large images.

`slicer.py` can divide larger microscopy images, as was done in *Hartsock
et al 2023* manuscript, into smaller patches. For our pipeline, we suggest using images/patches of size at most $1000 \times 1000$ pixels.

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
It is recommended that `pipeline.Rmd` is run in RStudio. The `cell_colors.csv`
file is used to specify how cells are labeled in the discretized images. Binary
strings encode high (1) or low (0) signal of a marker from the microscopy images.
The mapping from binary string to cell color should follow the order of the
measured signals from the segmentation section CSVs.

## TDA (R)
There are two R notebook files within this section. `general_pipeline.Rmd` should
be used for the example image set and for the most applications of the pipeline. After loading CSV files of discretized microscopy data, choose which cell types (which are based on signal codes) to use for TDA computations. Note that each block of code of this notebook correspond to different stage of the analysis. By running the notebook on the example image dataset, for each image, the persistence landscape is generated and plotted. Also, for each Dox treatment group, the average persistence landscape is produced and plotted. Additionally, in the final block of code, the permutation test is performed on the persistence landscapes of the two Dox treatment groups, and the p-value is calculated. 

The second notebook `hipsc_pipeline.Rmd` was configured specifically to the images used within the *Hartsock
et al 2023* manuscript. It requires installing [tda-tools](https://github.com/jjbouza/tda-tools) by Jose Bouza. Besides TDA analysis, this notebook also contains machine learning and statistical hypothesis testing computations.


