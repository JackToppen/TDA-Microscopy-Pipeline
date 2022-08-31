#Iryna Hartsock
#load the files

data_file_location <- "./data/"
save_computation_file_location <- "./patches_persistence_65p/"
save_computation_filename <- "16patches_vec.RData"
saved_computation_filepath <- paste0(save_computation_file_location,save_computation_filename)
#load_file <- saved_computation_filepath

percentile <- 65
concentrations <- c(0,5,15,25)
cohorts <- c("Gata6","HA")
nconcentrations <- length(concentrations)
ncohorts <- length(cohorts)
nimages <- 15 #number of images per concentration
npatches <- 16 #number of patches per image
nsamples <- nimages*npatches #number of samples per concentration 


files <- list()
for (c in cohorts){
  for (i in concentrations){ # gata6
    for (j in 1:3){
      for (k in 1:5){
        files <- append(files, paste0(percentile,"/Nanog_", c,"/concentration_", i, "/Nanog_", c,"_",i, "_",j,"_", k))
      }
    }
  }
}

data_filenames <- list()
for (i in 1:length(files)){
  data_filenames[i] <- paste0(files[i],".csv")
}



