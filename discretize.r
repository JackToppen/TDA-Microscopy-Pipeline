# set directory
setwd("~/Research/TDA/R")

for (percent in list(75, 80, 85, 90, 95)) {
  for (rfp in list("Gata6", "HA")) {
    
    gfp_thresh_path <- paste("~/Research/TDA/R/thresholds/", toupper(rfp), "_GFP_", percent, ".csv", sep="")
    rfp_thresh_path <- paste("~/Research/TDA/R/thresholds/", toupper(rfp), "_RFP_", percent, ".csv", sep="")
    
    # open CSV files
    gfp_thresh_data <- read.csv(gfp_thresh_path, header=FALSE) 
    rfp_thresh_data <- read.csv(rfp_thresh_path, header=FALSE) 
    
    for (conc in list(0, 5, 15, 25)) {
      # get path to current directory
      old_path <- paste("~/Research/TDA/R/data/Nanog_", rfp, "/concentration_", conc, "/", sep="")
      new_path <- paste("~/Research/TDA/R/discretized_data/", percent, "/Nanog_", rfp, "/concentration_", conc, "/", sep="")
      
      # make directory if it doesn't exist
      dir.create(new_path, recursive=TRUE, showWarnings=FALSE)
      
      for (well in list(1, 2, 3)) {
        for (loc in list(1, 2, 3, 4, 5)) {
          # get thresholds based on well
          gfp_thresh <- gfp_thresh_data[loc, 2]
          rfp_thresh <- rfp_thresh_data[loc, 2]
          
          # get file name
          name <- paste("Nanog_", rfp, "_", conc, "_", well, "_", loc, ".csv", sep="")
          old_name <- paste(old_path, name, sep="")
          new_name <- paste(new_path, name, sep="")
          
          # open CSV file
          data <- read.csv(old_name, colClasses=rep("numeric", 7))
          
          # convert the GFP/RFP floats to discrete values 
          data$gfp_discrete <- ifelse(data[, 6] < gfp_thresh, 0, 1)
          data$rfp_discrete <- ifelse(data[, 7] < rfp_thresh, 0, 1)
          
          # write the CSV
          write.csv(data, new_name, row.names=FALSE)
        }
      }
      print(conc)
    }
    print(rfp)
  }
  print(percent)
}
