#Iryna Hartsock
#compute persistence diagrams

library(tdatools)
library(TDA)
library(tictoc)

tic()
source("config.R")
source("source/intro_tda.R")
source("source/landscape_utilities.R")
save_computations <- TRUE

# moved these into the config file
save_file_location <- save_computation_file_location
save_filename <- save_computation_filename

PH_red <- list()
PH_green <- list()
PH_pos <- list()
PH_neg <- list()
for (i in 1:length(data_filenames)){
  tic()
  print(sprintf("Processing file %s", data_filenames[[i]]))
  PH_red[[i]] <- list()
  PH_green[[i]] <- list()
  PH_pos[[i]] <- list()
  PH_neg[[i]] <- list()
  for (k in 1:sqrt(npatches)){ #y-axis
    for (j in 1:sqrt(npatches)){ #x-axis
      all_cells <- read.csv(file=paste0(data_file_location,data_filenames[[i]]), header=TRUE)
      #every image is 3000x3000 
      all_cells <- all_cells[which( (3000/sqrt(npatches))*(j-1) <= all_cells[,1] & all_cells[,1] < (3000/sqrt(npatches))*j & (3000/sqrt(npatches))*(k-1) <= all_cells[,2] & all_cells[,2] < (3000/sqrt(npatches))*k),]
      red_cells <- all_cells[ which(all_cells[,8]==0 & all_cells[,9]==1), ]      #red cells           
      green_cells <- all_cells[ which(all_cells[,8]==1 & all_cells[,9]==0), ]    #green cells
      pos_cells <- all_cells[ which(all_cells[,8]==1 & all_cells[,9]==1), ]      #double positive cells
      neg_cells <- all_cells[ which(all_cells[,8]==0 & all_cells[,9]==0), ]      #double negative cells
      #compute persistence homology using Alpha complex which is also known as Delaunay complex
      PH_red[[i]][[j+(k-1)*sqrt(npatches)]] <-  alphaComplexDiag(red_cells[,1:2], maxdimension = 1, library = c("GUDHI", "Dionysus"), 
                                            location = TRUE)
      PH_green[[i]][[j+(k-1)*sqrt(npatches)]] <-  alphaComplexDiag(green_cells[,1:2], maxdimension = 1, library = c("GUDHI", "Dionysus"), 
                                              location = TRUE)
      PH_pos[[i]][[j+(k-1)*sqrt(npatches)]] <-  alphaComplexDiag(pos_cells[,1:2], maxdimension = 1, library = c("GUDHI", "Dionysus"), 
                                            location = TRUE)
      PH_neg[[i]][[j+(k-1)*sqrt(npatches)]] <-  alphaComplexDiag(neg_cells[,1:2], maxdimension = 1, library = c("GUDHI", "Dionysus"), 
                                            location = TRUE)
      print(sprintf("    computations time for %s:", data_filenames[i]))
    }
  } 
  toc()
}

PDs_red <-list()
PDs_green <-list()
PDs_pos <- list()
PDs_neg <- list()
for (i in 1:length(data_filenames)){
  PDs_red[[i]] <-list()
  PDs_green[[i]] <-list()
  PDs_pos[[i]] <-list()
  PDs_neg[[i]] <- list()
  for (k in 1:npatches){
    #convert from gudhi to tdatools diagram data structure
    PDs_red[[i]][[k]] <- gudhi2tdatools(PH_red[[i]][[k]]$diagram) 
    PDs_green[[i]][[k]] <- gudhi2tdatools(PH_green[[i]][[k]]$diagram)
    PDs_pos[[i]][[k]] <- gudhi2tdatools(PH_pos[[i]][[k]]$diagram) 
    PDs_neg[[i]][[k]] <- gudhi2tdatools(PH_neg[[i]][[k]]$diagram) 
    #birth and death values have been squared, so take the square root
    #homology in degree 0
    PDs_red[[i]][[k]]$pairs[[1]] <- sqrt(PDs_red[[i]][[k]]$pairs[[1]]) 
    PDs_green[[i]][[k]]$pairs[[1]] <- sqrt(PDs_green[[i]][[k]]$pairs[[1]]) 
    PDs_pos[[i]][[k]]$pairs[[1]] <- sqrt(PDs_pos[[i]][[k]]$pairs[[1]])
    PDs_neg[[i]][[k]]$pairs[[1]] <- sqrt(PDs_neg[[i]][[k]]$pairs[[1]]) 
    #homology in degree 1
    PDs_green[[i]][[k]]$pairs[[2]] <- sqrt(PDs_green[[i]][[k]]$pairs[[2]]) 
    PDs_red[[i]][[k]]$pairs[[2]] <- sqrt(PDs_red[[i]][[k]]$pairs[[2]]) 
    PDs_pos[[i]][[k]]$pairs[[2]] <- sqrt(PDs_pos[[i]][[k]]$pairs[[2]]) 
    PDs_neg[[i]][[k]]$pairs[[2]] <- sqrt(PDs_neg[[i]][[k]]$pairs[[2]]) 
  }
}

 toc()

#save data
if (save_computations){
  README <- "computed persistence diagrams for red, green, double positive, double negative cells"
  save(vec_red,vec_gr, vec_pos, vec_neg, 
    PL_length_red, PL_length_gr, PL_length_pos, PL_length_neg,
    #PDs_red,PDs_green, PDs_pos, PDs_neg,
    README, 
    file=paste0(save_file_location,save_filename))
  print(sprintf("finished computations saved to %s%s", save_file_location, save_filename))
}

