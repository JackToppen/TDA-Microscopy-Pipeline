#Iryna Hartsock
#load the files

data_file_location <- "./data/"
save_computation_file_location <- "./new_patches_75p/"
save_computation_filename <- "norm_knn_SVR.RData"
saved_computation_filepath <- paste0(save_computation_file_location,save_computation_filename)
#load_file <- saved_computation_filepath

percentile <- 75
concentrations <- c(0,5,15,25)
cohorts <- c("Gata6","HA")
nconcentrations <- length(concentrations)
ncohorts <- length(cohorts)
nimages <- 15 #number of images per concentration
npatches <- 16 #number of patches per image
nsamples <- nimages*npatches #number of samples per concentration 


files <- list()
for (j in cohorts){
  for (c in concentrations){ # gata6
    for (w in 1:3){
      for (p in 1:5){
        files <- append(files, paste0(percentile,"/Nanog_", j,"/", p, "/", c,"/Nanog_",j, "_",c,"_", w, "_", p))
      }
    }
  }
}

data_filenames <- list()
for (i in 1:length(files)){
  data_filenames[i] <- paste0(files[i],"/merged.csv")
}



