library(ggplot2)
library(ggforce)
library(cowplot)

# set directory
setwd("~/Research/TDA/R")

for (percent in list(75, 80, 85, 90, 95)) {
  for (rfp in list("Gata6", "HA")) {
    for (conc in list(0, 5, 15, 25)) {
      # get path to current directory
      dir_path <- paste("~/Research/TDA/R/discretized_data/", percent, "/Nanog_", rfp, "/concentration_", conc, "/", sep="")
      image_path <- paste("~/Research/TDA/R/images/", percent, "/Nanog_", rfp, "/concentration_", conc, "/", sep="")
      
      # make directory if it doesn't exist
      dir.create(image_path, recursive=TRUE, showWarnings=FALSE)
      
      for (well in list(1, 2, 3)) {
        for (loc in list(1, 2, 3, 4, 5)) {
          for (both_low in list("present", "absent")) {
            # get file names
            name <- paste("Nanog_", rfp, "_", conc, "_", well, "_", loc, sep="")
            csv_name <- paste(dir_path, name, ".csv", sep="")
            image_name <- paste(image_path, name, "_", both_low, ".png", sep="")
            
            # open CSV file
            data <- read.csv(csv_name, colClasses=rep("numeric", 9))
            
            # find min/max for x and y data
            x_max <- max(data[, 1])
            y_max <- max(data[, 2])
            x_min <- min(data[, 1])
            y_min <- min(data[, 2])
            
            # get masks
            red_mask <- data[,9] > data[,8]
            white_mask <- data[,9] == 1 & data[,8] == 1
            if (both_low == "present") {
              green_mask <- data[,9] < data[,8] | (data[,9] == 0 & data[,8] == 0)
            } else {
              green_mask <- data[,9] < data[,8]
            }
            
            # create data frame to plot
            red_cells <- data.frame(
              x_r <- data[, 1][red_mask],
              y_r <- data[, 2][red_mask]
            )
            white_cells <- data.frame(
              x_w <- data[, 1][white_mask],
              y_w <- data[, 2][white_mask]
            )
            green_cells <- data.frame(
              x_g <- data[, 1][green_mask],
              y_g <- data[, 2][green_mask]
            )
            
            # create a scatter plot
            ggplot() + 
              geom_point(aes(x=x_r, y=y_r, color=I(rgb(1,0,0))), red_cells, size=0.1, stroke=0.5, shape=16) +
              geom_point(aes(x=x_w, y=y_w, color=I(rgb(1,1,1))), white_cells, size=0.1, stroke=0.5, shape=16) +
              geom_point(aes(x=x_g, y=y_g, color=I(rgb(0,1,0))), green_cells, size=0.1, stroke=0.5, shape=16) +
              coord_cartesian(xlim=c(x_min, x_max), ylim=c(y_min, y_max), expand=FALSE) + 
              labs(x=NULL, y=NULL) +
              theme_nothing() +
              theme(panel.background = element_rect(fill = 'black', color = 'black'))
            
            # save the plot
            ggsave(image_name, width=5, height=5)
          }
        }
      }
      print(conc)
    }
    print(rfp)
  }
  print(percent)
}
