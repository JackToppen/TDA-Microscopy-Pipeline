#Iryna Hartsock
#Pick the best patches based on cell number count
#March 2022

library(tictoc)

tic()
source("config.R")
#save_computations <- TRUE

# moved these into the config file
#save_file_location <- save_computation_file_location
#save_filename <- save_computation_filename

#create a list of all patches, for each file there are 4 patches
cell_number_r <- list()
red_patch <- list()
red_patches <- list()
cell_number_gr <- list()
green_patch <- list()
green_patches <- list()
cell_number_r_pos <- list()
r_pos_patch <- list()
r_pos_patches <- list()
cell_number_gr_neg <- list()
gr_neg_patch <- list()
gr_neg_patches <- list()
for (j in 1:ncohorts){                    #j=1 Gata6, j=2 HA
  cell_number_r[[j]] <- list()
  red_patch[[j]] <- list()
  red_patches[[j]] <- list()
  cell_number_gr[[j]] <- list()
  green_patch[[j]] <- list()
  green_patches[[j]] <- list()
  cell_number_r_pos[[j]] <- list()
  r_pos_patch[[j]] <- list()
  r_pos_patches[[j]] <- list()
  cell_number_gr_neg[[j]] <- list()
  gr_neg_patch[[j]] <- list()
  gr_neg_patches[[j]] <- list()
  for (c in 1:nconcentrations){                  #c corresponds to 4 concentrations
    red_patch[[j]][[c]] <- list()
    cell_number_r[[j]][[c]]<- list()
    red_patches[[j]][[c]] <- list()
    green_patch[[j]][[c]] <- list()
    cell_number_gr[[j]][[c]]<- list()
    green_patches[[j]][[c]] <- list()
    r_pos_patch[[j]][[c]] <- list()
    cell_number_r_pos[[j]][[c]]<- list()
    r_pos_patches[[j]][[c]] <- list()
    gr_neg_patch[[j]][[c]] <- list()
    cell_number_gr_neg[[j]][[c]]<- list()
    gr_neg_patches[[j]][[c]] <- list()
    for (i in 1:nsamples){               #k corresponds to 4 patches
      red_patch[[j]][[c]][[i]] <- list()
      cell_number_r[[j]][[c]][[i]] <- list()
      green_patch[[j]][[c]][[i]] <- list()
      cell_number_gr[[j]][[c]][[i]] <- list()
      r_pos_patch[[j]][[c]][[i]] <- list()
      cell_number_r_pos[[j]][[c]][[i]] <- list()
      gr_neg_patch[[j]][[c]][[i]] <- list()
      cell_number_gr_neg[[j]][[c]][[i]] <- list()
      for (k in 1:1){
        all_cells <- read.csv(file=paste0(data_file_location,data_filenames[[i+(c-1)*15+(j-1)*60]]), header=TRUE)
        all_cells <- all_cells[which( 500 <= all_cells[,1] & all_cells[,1] < 2500 & 500 <= all_cells[,2] & all_cells[,2] < 2500),]
        #if (k==1){
        #  all_cells <- all_cells[which( 500 <= all_cells[,1] & all_cells[,1] < 1500 & 500 <= all_cells[,2] & all_cells[,2] < 1500),]
       # } else if (k==2){
       #   all_cells <- all_cells[which( 1500 <= all_cells[,1] & all_cells[,1] < 2500 & 500 <= all_cells[,2] & all_cells[,2] < 1500),]
       # } else if (k==3){
       #   all_cells <- all_cells[which( 500 <= all_cells[,1] & all_cells[,1] < 1500 & 1500 <= all_cells[,2] & all_cells[,2] < 2500),]
       # } else {
       #   all_cells <- all_cells[which( 1500 <= all_cells[,1] & all_cells[,1] < 2500 & 1500 <= all_cells[,2] & all_cells[,2] < 2500),]
       # }
        #red
        red_patch[[j]][[c]][[i]][[k]] <- all_cells[ which((all_cells[,8]==0 & all_cells[,9]==1)), ] 
        cell_number_r[[j]][[c]][[i]][[k]] <- nrow(red_patch[[j]][[c]][[i]][[k]])
        #red_patches[[j]][[c]][[i+(k-1)*15]] <- red_patch[[j]][[c]][[k]][[i]] #list of all patches for a fixed concentration and cohort
        #green
        green_patch[[j]][[c]][[i]][[k]] <- all_cells[ which((all_cells[,8]==1 & all_cells[,9]==0)), ] 
        cell_number_gr[[j]][[c]][[i]][[k]] <- nrow(green_patch[[j]][[c]][[i]][[k]])
        #green_patches[[j]][[c]][[i+(k-1)*15]] <- green_patch[[j]][[c]][[k]][[i]]
        #double positive
        r_pos_patch[[j]][[c]][[i]][[k]] <- all_cells[ which((all_cells[,8]==1 & all_cells[,9]==1)), ] 
        cell_number_r_pos[[j]][[c]][[i]][[k]] <- nrow(r_pos_patch[[j]][[c]][[i]][[k]])
        #r_pos_patches[[j]][[c]][[i+(k-1)*15]] <- r_pos_patch[[j]][[c]][[k]][[i]] 
        #double negative
        gr_neg_patch[[j]][[c]][[i]][[k]] <- all_cells[ which((all_cells[,8]==0 & all_cells[,9]==0)), ] 
        cell_number_gr_neg[[j]][[c]][[i]][[k]] <- nrow(gr_neg_patch[[j]][[c]][[i]][[k]])
        #gr_neg_patches[[j]][[c]][[i+(k-1)*15]] <- gr_neg_patch[[j]][[c]][[k]][[i]]  
      }
    }
  }
}


counts <- list()
for (j in 1:ncohorts){ 
  counts[[j]] <- list()
  for (c in 1:nconcentrations){
    counts[[j]][[c]] <- vector()
    for (i in 1:nsamples){
      counts[[j]][[c]] <- rbind(counts[[j]][[c]],c(cell_number_r[[j]][[c]][[i]][[1]], cell_number_gr[[j]][[c]][[i]][[1]], cell_number_r_pos[[j]][[c]][[i]][[k]], cell_number_gr_neg[[j]][[c]][[i]][[k]]))
    }
  }
}


