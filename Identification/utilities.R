library(ggplot2)
library(ggforce)
library(cowplot)


normalize <- function(path) {
  # get list of segmentation CSVs to normalize
  files <- list.files(path, pattern="\\.csv$", include.dirs=TRUE, recursive=TRUE)
  
  # iterate through all files
  for (file in files) {
    # read the segmentation CSV to a data frame
    csv_file <- concat_path(path, file)
    df <- read.csv(csv_file)
    
    # iterate through other channels
    for (col_index in 4:ncol(df)) {
      # get column name
      channel_name <- colnames(df)[col_index]
      
      # check that column hasn't been normalized already
      if (!grepl("_normalized", channel_name, fixed=TRUE)) {
        # get new column name for normalized channel signal
        norm_name <- paste0(channel_name, "_normalized")
        
        # create new column with normalized signal (dividing by nuclear signal)
        df[, norm_name] <- df[, col_index] / df[, 3]
      }
    }
    
    # save the CSV file with normalized intensity columns
    write.csv(df, csv_file, row.names=FALSE)
  }
}


get_thresholds <- function(percentile, paths) {
  # get number of channels
  channels <- length(paths)
  
  # create vector for storing generated threshold for each channel
  thresholds <- numeric(channels)
  
  # iterate through channels
  for (i in 1:channels) {
    # read baseline CSVs used to generate thresholds
    files <- list.files(paths[i], pattern="\\.csv$", include.dirs=TRUE, recursive=TRUE)
    
    # temporary variables for weighting channel threshold by number of cells
    temp_sum <- 0
    temp_mult <- 0
    
    # iterate through all files
    for (file in files) {
      # read the segmentation CSV to a data frame
      csv_file <- concat_path(paths[i], file)
      df <- read.csv(csv_file)
      
      # get the channel intensities
      col_index <- ncol(df) - channels + i
      intensities <- df[, col_index]
      
      # generate temp threshold to weight based on number of cells
      temp_thresh <- quantile(intensities, probs=percentile/100)
      num_cells <- nrow(df)
      temp_sum <- temp_sum + num_cells
      temp_mult <- temp_mult + num_cells * temp_thresh
    }
    
    # get the weighted channel threshold
    thresholds[i] <- temp_mult / temp_sum
  }
  
  # return vector with each index corresponding to channel threshold
  return(thresholds)
}


discretize <- function(path, disc_path, thresholds) {
  # get number of channels
  channels <- length(paths)
  
  # make output directory for discretized CSVs
  make_dir(disc_path)
  
  # get list of normalized segmentation CSVs to discretize
  files <- list.files(path, pattern="\\.csv$", include.dirs=TRUE, recursive=TRUE)
  
  # iterate through all files
  for (file in files) {
    # read the normalized segmentation CSV to a data frame
    csv_file <- concat_path(path, file)
    df <- read.csv(csv_file)
    
    # get total number of columns
    num_cols <- ncol(df)
    
    # create a column for storing the cell type as a Boolean string
    df[,"type"] <- ""
      
    # iterate through channels
    for (i in 1:channels) {
      # get index of channel column
      col_index <- num_cols - channels + i
      
      # determine discretized signal and append to Boolean string
      temp <- ifelse(df[, col_index] > thresholds[i], 1, 0)
      df[,"type"] <- paste0(df[,"type"], temp)
    }
    
    # get path to new CSV and make directory for it
    disc_file <- concat_path(disc_path, file)
    make_dir(dirname(disc_file))
    
    # save the CSV file with Boolean string cell type
    write.csv(df, disc_file, row.names=FALSE)
  }
}


generate_images <- function(disc_path, image_path, colors_csv, size, outline, resolution) {
  # make output directory for images
  make_dir(image_path)
  
  # read the cell colors CSV to a data frame
  cell_colors <- read.csv(colors_csv)

  # get list of discretized CSVs to image
  files <- list.files(disc_path, pattern="\\.csv$", include.dirs=TRUE, recursive=TRUE)
  
  # iterate through all files
  for (file in files) {
    # read the discretized CSV to a data frame
    disc_file <- concat_path(disc_path, file)
    df <- read.csv(disc_file)
    
    # create a column in the data frame for cell colors
    df["color"] <- NA
    
    # map cell type to color
    num_types <- nrow(cell_colors)
    for (i in 1:num_types) {
      df["color"][df["type"] == cell_colors[i,1]] <- cell_colors[i,2]
    }

    # create a scatter plot to display the cells
    x_max <- max(df["X"])
    y_max <- max(df["Y"])
    ggplot(df, aes(X, Y, fill=I(color))) + 
      geom_point(size=(size/5.9), stroke=(outline/11.8), shape=21) +
      coord_cartesian(xlim=c(0, x_max), ylim=c(0, y_max), expand=FALSE) + 
      labs(x=NULL, y=NULL) +
      theme_nothing() +
      theme(panel.background = element_rect(fill='black', color='black'))
    
    # convert filename from CSV to PNG
    no_ext <- sub('\\.csv$', '', file)
    image_local <- paste0(no_ext, ".png")
    
    # get path to new PNG and make directory for it
    image_file <- concat_path(image_path, image_local)
    make_dir(dirname(image_file))
    
    # convert width to inches and get height
    width <- resolution / 300    # 300 pixels per inch default
    height <- width * (x_max/y_max)
    
    # write image file
    ggsave(image_file, width=width, height=height)
  }
}


make_dir <- function(path) {
  # check if directory exists and if not, make it
  if (!dir.exists(path)) {
    dir.create(path, recursive=TRUE)
  }
}


concat_path <- function(path, filename) {
  # make sure path has file separator at the end
  n <- nchar(path)
  if (substr(path, n, n) != .Platform$file.sep) {
    path <- paste0(path, .Platform$file.sep)
  }
  
  # return correct joined path to file
  return(paste0(path, filename))
}

