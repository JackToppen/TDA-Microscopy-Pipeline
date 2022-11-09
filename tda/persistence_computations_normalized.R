#Iryna Hartsock
#compute compute persistence homology using Alpha complex
#which is also known as Delaunay complex 
#and persistence diagrams

library(tdatools)
library(TDA)
library(tictoc)
library(FNN)

tic()
source("config.R")
source("utilities.R")
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
  print(sprintf("%1.0f Processing file %s", i, data_filenames[[i]]))
  PH_red[[i]] <- list()
  PH_green[[i]] <- list()
  PH_pos[[i]] <- list()
  PH_neg[[i]] <- list()
  for (y in 1:sqrt(npatches)){ #y-axis
    for (x in 1:sqrt(npatches)){ #x-axis
      all_cells <- read.csv(file=paste0(data_file_location,data_filenames[[i]]), header=TRUE)
      #every image is approximately 2850x2850 
      all_cells <- all_cells[which( (2850/sqrt(npatches))*(x-1) <= all_cells[,1] & all_cells[,1] < (2850/sqrt(npatches))*x & (2850/sqrt(npatches))*(y-1) <= all_cells[,2] & all_cells[,2] < (2850/sqrt(npatches))*y),]
      #red cells 
      red_cells <- all_cells[ which(all_cells[,8]== 1), ]  
      if (nrow(red_cells) >10){
        k_dist <- knn.dist(red_cells[,1:2], k=10, algorithm = "kd_tree")
        avg_dist <- mean(k_dist[,10])
      } else {avg_dist <- (2850/sqrt(npatches))*sqrt(2)}
      #M <- as.matrix(pdist(red_cells[,1:2]))
      #M <- M[upper.tri(M)]
      #avg_dist <- mean(M)
      PH_red[[i]][[x+(y-1)*sqrt(npatches)]] <-  alphaComplexDiag(red_cells[,1:2]*(1/avg_dist), maxdimension = 1, library = c("GUDHI", "Dionysus"), 
                                                                 location = TRUE)
      #green cells
      green_cells <- all_cells[ which(all_cells[,8]==10), ]   
      if (nrow(green_cells) >10){
        k_dist <- knn.dist(green_cells[,1:2], k=10, algorithm = "kd_tree")
        avg_dist <- mean(k_dist[,10])
      } else {avg_dist <- (2850/sqrt(npatches))*sqrt(2)}
      #M <- as.matrix(pdist(green_cells[,1:2]))
      #M <- M[upper.tri(M)]
      #avg_dist <- mean(M)
      PH_green[[i]][[x+(y-1)*sqrt(npatches)]] <-  alphaComplexDiag(green_cells[,1:2]*(1/avg_dist), maxdimension = 1, library = c("GUDHI", "Dionysus"), 
                                                                   location = TRUE)
      #double positive cells
      pos_cells <- all_cells[ which(all_cells[,8]==11), ] 
      if (nrow(pos_cells) >10){
        k_dist <- knn.dist(pos_cells[,1:2], k=10, algorithm = "kd_tree")
        avg_dist <- mean(k_dist[,10])
      } else {avg_dist <- (2850/sqrt(npatches))*sqrt(2)}
      #M <- as.matrix(pdist(pos_cells[,1:2]))
      #M <- M[upper.tri(M)]
      #avg_dist <- mean(M)
      PH_pos[[i]][[x+(y-1)*sqrt(npatches)]] <-  alphaComplexDiag(pos_cells[,1:2]*(1/avg_dist), maxdimension = 1, library = c("GUDHI", "Dionysus"), 
                                                                 location = TRUE)
      #double negative cells
      neg_cells <- all_cells[ which(all_cells[,8]==0), ]
      if (nrow(neg_cells) >10){
        k_dist <- knn.dist(neg_cells[,1:2], k=10, algorithm = "kd_tree")
        avg_dist <- mean(k_dist[,10])
      } else {avg_dist <- (2850/sqrt(npatches))*sqrt(2)}
      #M <- as.matrix(pdist(neg_cells[,1:2]))
      #M <- M[upper.tri(M)]
      #avg_dist <- mean(M)
      PH_neg[[i]][[x+(y-1)*sqrt(npatches)]] <-  alphaComplexDiag(neg_cells[,1:2]*(1/avg_dist), maxdimension = 1, library = c("GUDHI", "Dionysus"), 
                                                                 location = TRUE)
      #print(sprintf("    computations time for %s:", data_filenames[i]))
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
  save(#prediction_H,
       #SVR_H_med,
       #SVR_H_mean,
       #SVR_H_mean6_11,
       #SVR_H_med_sc0_30,
       #SVR_H_mean_sc0_30,
       #SVR_H_mean_6_11_sc0_30,
       #SVR_G_mean,
       #SVR_G_med,
       #SVR_H_mean, 
       #SVR_H_med,
       #vec_red, PL_length_red, vec_gr, PL_length_gr, vec_pos, PL_length_pos, vec_neg, PL_length_neg,
       #vec_gr, vec_pos, vec_neg, PL_length_gr, PL_length_pos, PL_length_neg,
    PDs_red,PDs_green, PDs_pos, PDs_neg,
    README, 
    file=paste0(save_file_location,save_filename))
  print(sprintf("finished computations saved to %s%s", save_file_location, save_filename))
}

