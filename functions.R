library(ggplot2)
library(ggforce)
library(cowplot)


save_to_csv <- function(df, file_name) {
  # write current data frame to a CSV file
  write.csv(df, file_name, row.names=FALSE)
}


normalize <- function(df) {
  # get the values for DAPI nuclear stain
  dapi <- df$dapi
  
  # normalize the GFP and RFP values
  df$gfp_normalized <- df$gfp_bgc / dapi
  df$rfp_normalized <- df$rfp_bgc / dapi
  
  return(df)
}


gen_threshold <- function(percentile, baseline) {
  # temporary variables for weighting threshold by cell numbers
  temp_sum <- 0
  temp_mult <- 0
  
  # iterate through baseline expression values
  for (values in baseline) {
    temp_thresh <- quantile(values, probs=percentile / 100)
    num_cells <- length(values)
    temp_sum <- temp_sum + num_cells
    temp_mult <- temp_mult + num_cells * temp_thresh
  }
  
  # get the weighted threshold
  return(temp_mult / temp_sum)
}


discretize <- function(df, gfp_thresh, rfp_thresh) {
  # normalize the GFP and RFP values
  df$gfp_discrete = ifelse(df$gfp_normalized > gfp_thresh, 1, 0)
  df$rfp_discrete = ifelse(df$rfp_normalized > rfp_thresh, 1, 0)
  
  return(df)
}


# used by mapply() in assign_color()
colorscheme <- function(gfp, rfp) {
  # assign colors based on the discrete GFP and RFP values
  if (rfp == 1) {
    if (gfp == 1) {
      # if both high, color white
      return("white")
      
    } else {
      # if RFP high and GFP low, color red
      return("red")
    }
  } else {
    # otherwise green
    return("green")
  }
}


# used by mapply() in assign_color()
alpha_values <- function(gfp, rfp) {
  # if RFP and GFP are both low
  if (gfp == 0 && rfp == 0) {
    # set alpha to 0
    return(0)
  } else {
    # set alpha to max
    return(1)
  }
}


assign_color <- function(df, show_double_low=TRUE) {
  # apply color scheme based on GFP and RFP values
  df$color <- mapply(colorscheme, df$gfp_discrete, df$rfp_discrete)
  
  # if imaging double-low cells
  if (show_double_low) {
    # set alpha value for transparency to max
    df$show <- 1
  } else {
    # set transparency values for double-low to zero
    df$show <- mapply(alpha_values, df$gfp_discrete, df$rfp_discrete)
  }
  
  return(df)
}


save_image <- function(df, image_file) {
  # create a scatter plot
  ggplot(df, aes(X, Y, color=I(color), alpha=I(show))) + 
    geom_point(size=0.1, stroke=0.5, shape=16) +
    coord_cartesian(xlim=range(df$X), ylim=range(df$Y), expand=FALSE) + 
    labs(x=NULL, y=NULL) +
    theme_nothing() +
    theme(panel.background = element_rect(fill='black', color='black'))
  
  # save the plot
  ggsave(image_file, width=5, height=5)
}


gen_table <- function(df, gfp_thresh, rfp_thresh) {
  # count the number of cells and get the GFP/RFP data frame columns
  total_cells = nrow(df)
  GFP <- df[, 6]
  RFP <- df[, 7]
  
  # get the distribution of cell types
  M00 = sum(((RFP < rfp_thresh) + (GFP < gfp_thresh)) == 2) / total_cells
  M11 = sum(((RFP > rfp_thresh) + (GFP > gfp_thresh)) == 2) / total_cells
  M01 = sum(((RFP < rfp_thresh) + (GFP > gfp_thresh)) == 2) / total_cells
  M10 = sum(((RFP > rfp_thresh) + (GFP < gfp_thresh)) == 2) / total_cells
  
  # store values in a 2x2 matrix
  return(matrix(c(M10, M00, M11, M01), ncol=2))
}
