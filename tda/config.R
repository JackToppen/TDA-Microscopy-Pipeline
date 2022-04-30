#Iryna Hartsock
# February 2022
#load the files

data_file_location <- "./data/"
#save_computation_file_location <- "./patches_persistence_75p/"
#save_computation_filename <- "vec_red_green.RData"
#saved_computation_filepath <- paste0(save_computation_file_location,save_computation_filename)
#load_file <- saved_computation_filepath

percentile <- 90
concentrations <- c(0,5,15,25)
cohorts <- c("Gata6","HA")
files <- list()
for (c in c("Gata6","HA")){
  for (i in concentrations){ # gata6
    for (j in 1:3){
      for (k in 1:5){
        files <- append(files, paste0(percentile,"/Nanog-", c,"/concentration-", i, "/Nanog-", c,"-",i, "-",j,"-", k))
      }
    }
  }
}

data_filenames <- list()
for (i in 1:length(files)){
  data_filenames[i] <- paste0(files[i],".csv")
}


nconcentrations <- length(concentrations)
ncohorts <- length(cohorts)
npatches <- 4
nsamples <- length(data_filenames)/(length(cohorts)*length(concentrations)) #number of samples per concentration 

