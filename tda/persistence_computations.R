#Iryna Hartsock
#March 2021


library(tdatools)
library(TDA)
library(tictoc)

tic()
source("config.R")
save_computations <- TRUE

# moved these into the config file
#save_file_location <- save_computation_file_location
#save_filename <- save_computation_filename

PH_red <- list()
PH_green <- list()
PH_pos <- list()
PH_neg <- list()
red_cells <- list()
green_cells <- list()
pos_cells <- list()
neg_cells <- list()
for (i in 1:length(data_filenames)){
  tic()
  print(sprintf("Processing file %s", data_filenames[[i]]))
  PH_red[[i]] <- list()
  PH_green[[i]] <- list()
  PH_pos[[i]] <- list()
  PH_neg[[i]] <- list()
  red_cells[[i]] <- list()
  green_cells[[i]] <- list()
  pos_cells[[i]] <- list()
  neg_cells[[i]] <- list()
  for (k in 1:npatches){
    all_cells <- read.csv(file=paste0(data_file_location,data_filenames[[i]]), header=TRUE)
    #select a patch
    if (k==1){
      all_cells <- all_cells[which( 500 <= all_cells[,1] & all_cells[,1] < 1500 & 500 <= all_cells[,2] & all_cells[,2] < 1500),]
    } else if (k==2){
      all_cells <- all_cells[which( 1500 <= all_cells[,1] & all_cells[,1] < 2500 & 500 <= all_cells[,2] & all_cells[,2] < 1500),]
    } else if (k==3){
      all_cells <- all_cells[which( 500 <= all_cells[,1] & all_cells[,1] < 1500 & 1500 <= all_cells[,2] & all_cells[,2] < 2500),]
    } else {
      all_cells <- all_cells[which( 1500 <= all_cells[,1] & all_cells[,1] < 2500 & 1500 <= all_cells[,2] & all_cells[,2] < 2500),]
    }
    red_cells[[i]][[k]] <- all_cells[ which(all_cells[,8]==0 & all_cells[,9]==1), ]      #red cells           
    green_cells[[i]][[k]] <- all_cells[ which(all_cells[,8]==1 & all_cells[,9]==0), ]    #green cells
    pos_cells[[i]][[k]] <- all_cells[ which(all_cells[,8]==1 & all_cells[,9]==1), ]      #double positive cells
    neg_cells[[i]][[k]] <- all_cells[ which(all_cells[,8]==0 & all_cells[,9]==0), ]      #double negative cells
    #compute persistence homology using Alpha complex which is also known as Delaunay complex
    PH_red[[i]][[k]] <-  alphaComplexDiag(red_cells[[i]][[k]][,1:2], maxdimension = 1, library = c("GUDHI", "Dionysus"), 
                                         location = TRUE)
    PH_green[[i]][[k]] <-  alphaComplexDiag(green_cells[[i]][[k]][,1:2], maxdimension = 1, library = c("GUDHI", "Dionysus"), 
                                            location = TRUE)
    PH_pos[[i]][[k]] <-  alphaComplexDiag(pos_cells[[i]][[k]][,1:2], maxdimension = 1, library = c("GUDHI", "Dionysus"), 
                                            location = TRUE)
    PH_neg[[i]][[k]] <-  alphaComplexDiag(neg_cells[[i]][[k]][,1:2], maxdimension = 1, library = c("GUDHI", "Dionysus"), 
                                          location = TRUE)
    print(sprintf("    computations time for %s:", data_filenames[i]))
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


# save data
#if (save_computations){
#  README <- "computed persistence diagrams for red, green, reg+green+double pos (4 patches)"
#  save(PDs_red, PDs_green, PDs_pos, PDs_neg,
#       data_file_location, data_filenames, 
#       README, 
#       file=paste0(save_file_location,save_filename))
#  print(sprintf("finished computations saved to %s%s", save_file_location, save_filename))
#}

print(sprintf("total computation time for %i files with 4 patches each:", length(data_filenames)))
toc()


























