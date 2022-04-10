
------------------------- Normalize GFP/RFP Values --------------------------
setwd("~/Documents/Research/TDA/R/")

# iterate over experimental conditions
for (signal in c("Nanog_Gata6")) {
  for (conc in c(0, 5, 15, 25)) {
    for (well in c(1, 2, 3)) {
      for (loc in c(1, 2, 3, 4, 5)) {
        # get path to CSV file
        old_path <- paste0("data/", signal, "/concentration_", conc, "/")

        # make new directory for storing normalized GFP/RFP values
        new_path <- paste0("normalized/", signal, "/concentration_", conc, "/")
        dir.create(new_path, recursive=TRUE, showWarnings=FALSE)

        # get CSV filename
        file_name <- paste0(signal, "_", conc, "_", well, "_", loc, ".csv")

        # read the CSV, normalize the GFP/RFP, and save to new CSV
        df <- read.csv(paste0(old_path, file_name))
        df <- normalize(df)
        save_to_csv(df, paste0(new_path, file_name))
      }
    }
  }
}


# --------------------- Reconstruct Microscopy Images -------------------------
setwd("~/Documents/Research/TDA/R/TDA-Microscopy-Pipeline/")
source("functions.R")
setwd("~/Documents/Research/TDA/R/")

# dox concentrations
concentrations <- c(0, 5, 15, 25)
percentiles <- c(75, 80, 85, 90, 95)

# matrices to store distributions of cell types between dox concentrations
M <- list()
for (perc_index in c(1, 2, 3, 4, 5)) {
  M[[perc_index]] <- list()
  for (dox_index in c(1, 2, 3, 4)) {
    M[[perc_index]][[dox_index]] <- matrix(rep(0,4), ncol=2)
  }
}

# iterate over RFP signal types (GATA6/HA)
for (sig in c("Nanog_Gata6")) {
  # iterate over well locations (1-5)
  for (loc in c(1, 2, 3, 4, 5)) {
    # store location specific normalized expression values
    gfp_values <- list()
    rfp_values <- list()
    
    # iterate through wells
    for (well in c(1, 2, 3)) {
      # directory path for baseline expression values
      dir_gfp <- paste0("normalized/", sig, "/concentration_25/")
      dir_rfp <- paste0("normalized/", sig, "/concentration_0/")
      
      # get CSV file names for baseline
      gfp_file <- paste0("Nanog_Gata6_25_", well, "_", loc, ".csv")
      rfp_file <- paste0("Nanog_Gata6_0_", well, "_", loc, ".csv")
      
      # read the CSV to data frames
      gfp_df <- read.csv(paste0(dir_gfp, gfp_file))
      rfp_df <- read.csv(paste0(dir_rfp, rfp_file))
      
      # hold location specific expression values by well for GFP/RFP
      gfp_values[[well]] <- gfp_df[, 6]
      rfp_values[[well]] <- rfp_df[, 7]
    }
    
    # iterate over threshold percentiles
    for (perc_index in c(1)) {
      # get dox concentration value
      perc <- percentiles[[perc_index]]
      
      # generate thresholds based percentile
      gfp_threshold <- gen_threshold(perc, gfp_values)
      rfp_threshold <- gen_threshold(perc, rfp_values)
      
      # iterate over doxycycline treatment concentrations
      for (dox_index in c(1, 2, 3, 4)) {
        # get dox concentration value
        conc <- concentrations[[dox_index]]
          
        # iterate through wells again after generating the thresholds
        for (well in c(1, 2, 3)) {
          # get path to CSV file and file name
          csv_path <- paste0("normalized/", sig, "/concentration_", conc, "/")
          csv_name <- paste0(sig, "_", conc, "_", well, "_", loc, ".csv")
          
          # read the CSV to a data frame
          df <- read.csv(paste0(csv_path, csv_name))
          
          # update the data frame with discretized GFP/RFP values
          df <- discretize(df, gfp_threshold, rfp_threshold)
          
          # assign colors based on discrete GFP/RFP values
          df <- assign_color(df, show_double_low=TRUE)
          
          # get/make directory for storing normalized GFP/RFP values
          image_path <- paste0("images/", perc, "/", sig, "/concentration_", conc, "/")
          dir.create(image_path, recursive=TRUE, showWarnings=FALSE)
          
          # create image and save
          image_name <- paste0(sig, "_", conc, "_", well, "_", loc, ".png")
          save_image(df, paste0(image_path, image_name))
          
          # generate distribution tables
          M[[perc_index]][[dox_index]] <- M[[perc_index]][[dox_index]] + gen_table(df, gfp_threshold, rfp_threshold)
        }
      }
    }
  }
}

# output distribution tables
print(M[[1]][[4]] / sum(M[[1]][[4]]))

